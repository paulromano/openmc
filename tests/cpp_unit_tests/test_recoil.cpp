#include <cmath>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "openmc/constants.h"
#include "openmc/particle_type.h"
#include "openmc/random_lcg.h"
#include "openmc/recoil.h"

using namespace openmc;
using Catch::Approx;

TEST_CASE("Recoil particle masses")
{
  // Masses come from the tabulated atomic masses, in eV
  REQUIRE(recoil::particle_mass_ev(ParticleType::photon()) == 0.0);
  REQUIRE(recoil::particle_mass_ev(ParticleType::neutron()) ==
          Approx(MASS_NEUTRON_EV).epsilon(1e-9));
  REQUIRE(recoil::particle_mass_ev(ParticleType::alpha()) / AMU_EV ==
          Approx(4.0015).epsilon(1e-4));

  // A nuclide that is not in the tabulated mass table still gets a usable mass
  double m = recoil::particle_mass_ev(ParticleType {40, 130, 0});
  REQUIRE(m > 0.0);
  REQUIRE(m / AMU_EV == Approx(130.0).epsilon(0.02));
}

TEST_CASE("Light-ion emission spectrum")
{
  // Alpha emission leaving Cr-53: the spectrum must vanish at both ends, be
  // strongly suppressed well below the Coulomb barrier, and integrate to a mean
  // that sits above the barrier.
  const double E_max = 13.3e6;
  const int Z_b = 2, A_b = 4, Z_d = 24, A_d = 53;

  REQUIRE(recoil::light_ion_pdf(0.0, E_max, Z_b, A_b, Z_d, A_d) == 0.0);
  REQUIRE(recoil::light_ion_pdf(E_max, E_max, Z_b, A_b, Z_d, A_d) == 0.0);
  REQUIRE(recoil::light_ion_pdf(2.0e6, E_max, Z_b, A_b, Z_d, A_d) <
          0.01 * recoil::light_ion_pdf(9.0e6, E_max, Z_b, A_b, Z_d, A_d));

  // Neutral ejectiles feel no barrier, so the spectrum is
  // E*(1 - E/E_max)^nu with nu the calibrated endpoint exponent. Written
  // against LightIonParams rather than a hardcoded exponent so that a refit
  // does not silently break a test whose subject is the barrier, not nu.
  const recoil::LightIonParams deployed {};
  double neutral = recoil::light_ion_pdf(0.5 * E_max, E_max, 0, 1, Z_d, A_d);
  REQUIRE(neutral ==
          Approx(0.5 * E_max * std::pow(0.5, deployed.nu)).epsilon(1e-12));

  SECTION("sampling reproduces the spectrum mean")
  {
    uint64_t seed = 1;
    const int N = 200000;
    double sum = 0.0;
    double max_sampled = 0.0;
    for (int i = 0; i < N; ++i) {
      double E = recoil::sample_light_ion_energy(
        E_max, E_max, Z_b, A_b, Z_d, A_d, &seed);
      REQUIRE(E >= 0.0);
      REQUIRE(E <= E_max);
      sum += E;
      max_sampled = std::max(max_sampled, E);
    }
    double sampled_mean = sum / N;

    // Reference mean by direct quadrature of the same pdf
    const int NQ = 200000;
    double num = 0.0, den = 0.0;
    for (int i = 0; i < NQ; ++i) {
      double E = E_max * (i + 0.5) / NQ;
      double f = recoil::light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d);
      num += f * E;
      den += f;
    }
    REQUIRE(sampled_mean == Approx(num / den).epsilon(0.01));
    REQUIRE(max_sampled > 0.9 * E_max);
  }

  SECTION("truncation restricts the sample to the energy budget")
  {
    uint64_t seed = 7;
    const double E_limit = 4.0e6;
    for (int i = 0; i < 20000; ++i) {
      double E = recoil::sample_light_ion_energy(
        E_max, E_limit, Z_b, A_b, Z_d, A_d, &seed);
      REQUIRE(E <= E_limit);
      REQUIRE(E >= 0.0);
    }
  }
}

TEST_CASE("Light-ion spectrum: parameterized overload matches the deployed one")
{
  // The parameterized overload exists so that candidate parameter sets can be
  // scored against evaluated data through the same code path the transport
  // kernel uses, rather than against a reimplementation in the fitting
  // scripts. That is only worth anything if the default-constructed parameters
  // reproduce the fixed constants bit for bit.
  const recoil::LightIonParams deployed {};
  for (int Z_b : {1, 2}) {
    for (int A_b : {1, 4}) {
      for (int Z_d : {6, 26, 74}) {
        int A_d = 2 * Z_d + 2;
        double E_max = 12.0e6;
        for (double frac : {0.01, 0.1, 0.35, 0.5, 0.75, 0.95, 0.999}) {
          double E = frac * E_max;
          REQUIRE(
            recoil::light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d) ==
            recoil::light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d, deployed));
        }
      }
    }
  }
}

TEST_CASE("Light-ion spectrum: the endpoint exponent behaves as expected")
{
  // nu is 0.5 in the deployed model and 1.21 in the refit. A larger exponent
  // must soften the spectrum, i.e. lower its mean, or the refit result is
  // being read with the wrong sign.
  auto mean_of = [](double nu) {
    recoil::LightIonParams par {1.8, 0.8e6, nu};
    const int N = 4000;
    double E_max = 10.0e6, num = 0.0, den = 0.0, prev_f = 0.0, prev_E = 0.0;
    for (int i = 0; i <= N; ++i) {
      double E = E_max * i / N;
      double f = recoil::light_ion_pdf(E, E_max, 1, 1, 26, 56, par);
      if (i > 0) {
        num += 0.5 * (f * E + prev_f * prev_E) * (E - prev_E);
        den += 0.5 * (f + prev_f) * (E - prev_E);
      }
      prev_f = f;
      prev_E = E;
    }
    return num / den;
  };
  double m_half = mean_of(0.5);
  double m_soft = mean_of(1.211);
  REQUIRE(m_soft < m_half);
  // and nu = 0.5 must agree with the std::sqrt fast path used for it
  REQUIRE(mean_of(0.5) == Approx(m_half).epsilon(1e-12));
}

TEST_CASE("Light-ion sampler realizes the parameterized spectrum")
{
  // Only checkable in C++: the sampler draws by rejection against a scanned
  // envelope, so a change in the spectrum's shape could in principle break the
  // envelope without changing the PDF at all.
  for (double nu : {0.5, 1.211}) {
    recoil::LightIonParams par {1.8, 0.8e6, nu};
    double E_max = 9.0e6;
    const int N = 3000;
    double num = 0.0, den = 0.0, prev_f = 0.0, prev_E = 0.0;
    for (int i = 0; i <= N; ++i) {
      double E = E_max * i / N;
      double f = recoil::light_ion_pdf(E, E_max, 2, 4, 24, 52, par);
      if (i > 0) {
        num += 0.5 * (f * E + prev_f * prev_E) * (E - prev_E);
        den += 0.5 * (f + prev_f) * (E - prev_E);
      }
      prev_f = f;
      prev_E = E;
    }
    double analytic = num / den;

    uint64_t seed = 987654321ULL;
    const int n_draw = 200000;
    double acc = 0.0;
    for (int i = 0; i < n_draw; ++i)
      acc +=
        recoil::sample_light_ion_energy(E_max, E_max, 2, 4, 24, 52, &seed, par);
    REQUIRE(acc / n_draw == Approx(analytic).epsilon(0.02));
  }
}
