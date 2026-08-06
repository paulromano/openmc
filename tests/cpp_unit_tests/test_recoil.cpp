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

  // Neutral ejectiles feel no barrier, so the spectrum is E*sqrt(1 - E/E_max)
  double neutral = recoil::light_ion_pdf(0.5 * E_max, E_max, 0, 1, Z_d, A_d);
  REQUIRE(neutral == Approx(0.5 * E_max * std::sqrt(0.5)).epsilon(1e-12));

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
