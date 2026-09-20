#include "openmc/recoil.h"

#include "openmc/atomic_mass.h"
#include "openmc/constants.h"
#include "openmc/distribution_multi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/reaction.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"

#include <algorithm> // for clamp, lower_bound, max, min
#include <cassert>
#include <cmath>  // for sqrt, exp, cbrt, log
#include <limits> // for numeric_limits

#include <fmt/core.h>

namespace openmc {
namespace recoil {

namespace {

//==============================================================================
// Tunable constants of the light-ion emission model
//
// The authoritative calibration is stored in the defaults of LightIonParams
// and AngularParams. Parameterized overloads allow validation to exercise the
// same formulas as transport.
//==============================================================================

//! Conversion from angstroms to femtometers
constexpr double ANGSTROM_TO_FM = 1.0e5;

//! Coulomb constant e^2/(4 pi epsilon_0) = alpha hbar c in [eV fm], derived
//! from OpenMC's CODATA 2018 constants
//! (https://physics.nist.gov/cuu/Constants/archive2018.html)
constexpr double COULOMB_EV_FM =
  PLANCK_C * ANGSTROM_TO_FM / (2.0 * PI * FINE_STRUCTURE);

//! Points in the tabulated cumulative the light-ion sampler inverts
//!
//! Sixty-four intervals put the sampled mean within 0.1% of the analytic one
//! over the validation cases while keeping the work per sample fixed.
constexpr int N_TABLE = 64;

//! Negative infinity, for a log density that is identically zero
constexpr double INFTY = std::numeric_limits<double>::infinity();

//! Nuclear numbers of a neutron
constexpr NuclearNumbers NEUTRON_NUMBERS {0, 1};

//! Maximum attempts to resample a modeled product inside the energy budget
//!
//! Evaluated inclusive neutron spectra can put appreciable probability in the
//! part of the marginal that is infeasible after another neutron has already
//! been drawn. A high bound makes exhaustion negligible without changing the
//! common-case cost, since feasible candidates are still accepted immediately.
constexpr int MAX_BUDGET_TRIES = 1000;

//==============================================================================
// Momentum helpers. Masses and energies are in eV and momenta in eV with c = 1,
// so that p = sqrt(2 m E) for a massive particle and p = E for a photon.
//==============================================================================

Direction momentum_from_energy(double mass, double E, Direction u)
{
  if (mass <= 0.0 || E <= 0.0)
    return {};
  return std::sqrt(2.0 * mass * E) * u;
}

Direction neutron_momentum(double E, Direction u)
{
  return momentum_from_energy(MASS_NEUTRON_EV, E, u);
}

//==============================================================================
// Exit-channel bookkeeping
//==============================================================================

ParticleType ion_type(NuclearNumbers za)
{
  if (za.Z == 0 && za.A == 1)
    return ParticleType::neutron();
  if (za.Z == 1 && za.A == 1)
    return ParticleType::proton();
  return ParticleType {za.Z, za.A, 0};
}

//! Push the charged ions of an exit channel onto a fixed-capacity list
struct ChargedProducts {
  static constexpr int MAX = 8;
  NuclearNumbers za[MAX];
  int size {0};

  explicit ChargedProducts(const ExitChannel& e)
  {
    add({1, 1}, e.proton);
    add({1, 2}, e.deuteron);
    add({1, 3}, e.triton);
    add({2, 3}, e.he3);
    add({2, 4}, e.alpha);
  }

  void add(NuclearNumbers value, int count)
  {
    assert(count >= 0 && count <= MAX - size);
    for (int i = 0; i < count; ++i)
      za[size++] = value;
  }

  span<NuclearNumbers> to_span()
  {
    return {za, static_cast<std::size_t>(size)};
  }
};

//! Fill \p out with the emitted particles of a channel, light ions and all
//!
//! The complement of ChargedProducts, which keeps only the charged ones. A mass
//! budget needs every particle that leaves.
int all_products(const ExitChannel& e, NuclearNumbers* out, int capacity)
{
  const NuclearNumbers kinds[6] = {
    NEUTRON_NUMBERS, {1, 1}, {1, 2}, {1, 3}, {2, 3}, {2, 4}};
  const int counts[6] = {
    e.neutron, e.proton, e.deuteron, e.triton, e.he3, e.alpha};
  int n = 0;
  for (int k = 0; k < 6; ++k) {
    if (counts[k] < 0 || counts[k] > capacity - n)
      return -1;
    for (int i = 0; i < counts[k]; ++i) {
      out[n++] = kinds[k];
    }
  }
  return n;
}

//==============================================================================
// Kalbach-Mann angular distribution
//
// The emitted ion direction uses the same functional form as evaluated MF=6
// LANG=2 data,
//
//     f(mu) = a/(2 sinh a) [cosh(a mu) + r sinh(a mu)],
//
// with the slope a from Kalbach's systematics and the pre-equilibrium fraction
// r taken as the fraction of the available energy carried by the ion. That
// reproduces the qualitative behavior of evaluated r values, which rise from
// nearly zero at low outgoing energy to 0.5-0.9 near the kinematic maximum.
//==============================================================================

//! Kalbach's semi-empirical separation energy in [MeV]
//!
//! The liquid-drop coefficients and emitted-particle binding correction are
//! from C. Kalbach, Phys. Rev. C 37, 2350 (1988),
//! https://doi.org/10.1103/PhysRevC.37.2350.
double kalbach_separation_energy(
  NuclearNumbers compound, NuclearNumbers daughter, NuclearNumbers emitted)
{
  double A_c = compound.A;
  double Z_c = compound.Z;
  double N_c = compound.A - compound.Z;
  double A_a = daughter.A;
  double Z_a = daughter.Z;
  double N_a = daughter.A - daughter.Z;

  // Binding energy of the emitted particle from the CODATA bare-particle
  // masses in atomic_mass.h, converted to MeV.
  double emitted_mass = nuclear_mass(emitted.Z, emitted.A);
  double binding_energy =
    emitted_mass > 0.0
      ? (emitted.Z * MASS_PROTON + (emitted.A - emitted.Z) * MASS_NEUTRON -
          emitted_mass) *
          AMU_EV / 1.0e6
      : 0.0;

  return 15.68 * (A_c - A_a) -
         28.07 *
           ((N_c - Z_c) * (N_c - Z_c) / A_c - (N_a - Z_a) * (N_a - Z_a) / A_a) -
         18.56 * (std::pow(A_c, 2.0 / 3.0) - std::pow(A_a, 2.0 / 3.0)) +
         33.22 * ((N_c - Z_c) * (N_c - Z_c) / std::pow(A_c, 4.0 / 3.0) -
                   (N_a - Z_a) * (N_a - Z_a) / std::pow(A_a, 4.0 / 3.0)) -
         0.717 * (Z_c * Z_c / std::cbrt(A_c) - Z_a * Z_a / std::cbrt(A_a)) +
         1.211 * (Z_c * Z_c / A_c - Z_a * Z_a / A_a) - binding_energy;
}

//! Kalbach angular-distribution slope parameter

//==============================================================================
// Emission bookkeeping for one collision
//==============================================================================

//! Running state of a reaction's product accounting
//!
//! Tracks the laboratory momentum, mass, and remaining internal (excitation
//! plus kinetic) energy of the system that has not yet decayed. Emitting a
//! particle removes its momentum and its share of the internal energy, so the
//! sum of the modeled kinetic energies can never exceed the budget.
struct EmissionState {
  Direction momentum {};         //!< lab momentum of undecayed system [eV]
  double mass {0.0};             //!< mass of undecayed system [eV]
  double energy_available {0.0}; //!< rest-frame kinetic energy [eV]
  double energy_emitted {0.0};   //!< lab kinetic energy given to products [eV]
  NuclearNumbers za {};          //!< identity of undecayed system
};

//! Modeled light ion kept until the event can be completed and banked
struct SampledIon {
  ParticleType type;
  Direction direction;
  double energy;
};

//! Update the internal energy and reject a materially impossible state
bool update_internal_energy(EmissionState& state, double budget)
{
  if (state.mass <= 0.0 || !std::isfinite(state.mass))
    return false;
  double translation = state.momentum.dot(state.momentum) / (2.0 * state.mass);
  if (!std::isfinite(translation) || translation < 0.0)
    return false;
  double scale = std::max(
    {1.0, std::abs(budget), std::abs(state.energy_emitted), translation});
  double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  double available = remaining_internal_energy(
    budget, state.energy_emitted, state.momentum, state.mass);
  if (!std::isfinite(available) || available < -tolerance)
    return false;
  state.energy_available = std::max(0.0, available);
  return true;
}

//! How far an evaluated level Q may sit above the ground-state mass budget
//!
//! A level Q above the ground-state budget means a negative excitation, which
//! is unphysical. This tolerance accommodates differences between evaluated
//! and AME2020 masses while remaining well below a physical level spacing that
//! would indicate an inconsistent evaluation.
constexpr double Q_LEVEL_TOLERANCE = 0.25e6;

} // namespace

//==============================================================================
// Kinematic contract
//
// Shared expressions, documented in recoil.h, mirrored by recoil.kinematics in
// the analysis repository and checked against a shared fixture.
//==============================================================================

double nuclear_mass_ev(NuclearNumbers za)
{
  return nuclear_mass(za.Z, za.A) * AMU_EV;
}

double mass_excess_ev(NuclearNumbers za)
{
  double m = nuclear_mass_ev(za);
  return m <= 0.0 ? 0.0 : m - za.A * AMU_EV;
}

std::optional<double> mass_difference_q(
  NuclearNumbers target, span<const NuclearNumbers> emitted)
{
  NuclearNumbers daughter {target.Z, target.A + 1};
  // Mass excesses rather than masses. A W-184 channel differences four numbers
  // of order 1.7e11 eV to reach one of order 1e6, and a double carries about
  // sixteen digits, so the direct subtraction returns a Q good to only ten.
  // The mass numbers cancel identically because A_T + 1 = A_D + sum A_j, so
  // working in excesses of order 1e7 eV gives the same Q to fifteen.
  double exit_excess = 0.0;
  for (auto product : emitted) {
    daughter -= product;
    if (nuclear_mass_ev(product) <= 0.0)
      return std::nullopt;
    exit_excess += mass_excess_ev(product);
  }
  if (!target.is_valid() || !daughter.is_valid())
    return std::nullopt;
  if (nuclear_mass_ev(target) <= 0.0 || nuclear_mass_ev(daughter) <= 0.0)
    return std::nullopt;

  return mass_excess_ev(target) + mass_excess_ev(NEUTRON_NUMBERS) -
         mass_excess_ev(daughter) - exit_excess;
}

void initialize_reaction(const Nuclide& nuc, Reaction& rx)
{
  rx.recoil_ = {};
  auto emitted = reaction_exit_channel(rx.mt_);
  if (!emitted)
    return;

  NuclearNumbers products[16];
  int n_products = all_products(*emitted, products, 16);
  if (n_products < 0)
    return;

  NuclearNumbers target {nuc.Z_, nuc.A_};
  NuclearNumbers residual {target.Z, target.A + 1};
  for (int i = 0; i < n_products; ++i)
    residual -= products[i];
  if (!residual.is_valid())
    return;

  auto q_mass =
    mass_difference_q(target, span<const NuclearNumbers> {products,
                                static_cast<std::size_t>(n_products)});
  double q = q_mass.value_or(rx.q_value_);
  if (is_discrete_level(rx.mt_)) {
    if (!q_mass || rx.q_value_ <= *q_mass) {
      q = rx.q_value_;
    } else if (rx.q_value_ <= *q_mass + Q_LEVEL_TOLERANCE) {
      q = *q_mass;
    } else {
      if (mpi::master) {
        warning(fmt::format(
          "Recoil production is disabled for {} MT={} because its evaluated "
          "Q value ({:.6g} MeV) exceeds the ground-state mass-difference Q "
          "value ({:.6g} MeV) by more than {:.3g} MeV.",
          nuc.name_, rx.mt_, rx.q_value_ / 1.0e6, *q_mass / 1.0e6,
          Q_LEVEL_TOLERANCE / 1.0e6));
      }
      return;
    }
  }

  rx.recoil_.supported = true;
  rx.recoil_.emitted = *emitted;
  rx.recoil_.residual = ion_type(residual);
  rx.recoil_.q_value = q;
}

double final_state_internal_energy(double E_in, double m_final, double q)
{
  if (m_final <= 0.0)
    return E_in + q;
  return E_in + q - E_in * MASS_NEUTRON_EV / m_final;
}

double two_body_endpoint(double u, double m_b, double m_d)
{
  if (u <= 0.0 || m_b <= 0.0 || m_d <= 0.0)
    return 0.0;
  return u * m_d / (m_b + m_d);
}

double residual_excitation(double u, double e_cm, double m_b, double m_d)
{
  if (m_d <= 0.0)
    return u;
  return u - e_cm * (1.0 + m_b / m_d);
}

double remaining_internal_energy(
  double budget, double energy_emitted, Direction momentum, double mass)
{
  if (mass <= 0.0 || !std::isfinite(mass))
    return -INFTY;
  return budget - energy_emitted - momentum.dot(momentum) / (2.0 * mass);
}

double shape_endpoint(double E_in, NuclearNumbers target, NuclearNumbers ion)
{
  auto q = mass_difference_q(target, span<const NuclearNumbers> {&ion, 1});
  if (!q)
    return 0.0;
  NuclearNumbers daughter {target.Z - ion.Z, target.A + 1 - ion.A};
  double m_b = nuclear_mass_ev(ion);
  double m_d = nuclear_mass_ev(daughter);
  double u = final_state_internal_energy(E_in, m_b + m_d, *q);
  return two_body_endpoint(u, m_b, m_d);
}

//==============================================================================
// Public helpers
//==============================================================================

double kalbach_slope(
  double E_in, double E_cm, NuclearNumbers emitted, int Z_t, int A_t)
{
  NuclearNumbers target {Z_t, A_t};
  NuclearNumbers compound {target.Z, target.A + 1};
  NuclearNumbers recoil = compound - emitted;
  if (!recoil.is_valid())
    return 0.0;

  double epsilon_a = E_in * target.A / (target.A + 1.0) / 1.0e6;
  double epsilon_b = E_cm * (recoil.A + emitted.A) / (recoil.A * 1.0e6);
  double e_a =
    epsilon_a + kalbach_separation_energy(compound, target, NEUTRON_NUMBERS);
  double e_b = epsilon_b + kalbach_separation_energy(compound, recoil, emitted);
  if (e_a <= 0.0 || e_b <= 0.0 || !std::isfinite(e_a) || !std::isfinite(e_b))
    return 0.0;

  double r_1 = std::min(e_a, 130.0);
  double r_3 = std::min(e_a, 41.0);
  double x_1 = r_1 * e_b / e_a;
  double x_3 = r_3 * e_b / e_a;
  double m = (emitted.Z == 2 && emitted.A == 4) ? 2.0 : 1.0;
  return std::max(0.0,
    0.04 * x_1 + 1.8e-6 * x_1 * x_1 * x_1 + 6.7e-7 * m * x_3 * x_3 * x_3 * x_3);
}

//! Sample the Kalbach-Mann angular distribution
double sample_kalbach_mu(double slope, double r, uint64_t* seed)
{
  if (slope <= 1.0e-8)
    return uniform_distribution(-1.0, 1.0, seed);

  double xi = prn(seed);
  if (prn(seed) < r) {
    // forward component proportional to exp(a mu)
    double exp_neg_2a = slope < 350.0 ? std::exp(-2.0 * slope) : 0.0;
    double mu = 1.0 + std::log(xi + (1.0 - xi) * exp_neg_2a) / slope;
    return std::clamp(mu, -1.0, 1.0);
  }
  // symmetric component proportional to cosh(a mu)
  double t = (2.0 * xi - 1.0) * std::sinh(slope);
  double mu = std::log(t + std::sqrt(t * t + 1.0)) / slope;
  return std::clamp(mu, -1.0, 1.0);
}

double particle_mass_ev(ParticleType type)
{
  if (type.is_photon())
    return 0.0;
  if (std::abs(type.pdg_number()) == PDG_ELECTRON)
    return MASS_ELECTRON * AMU_EV;
  NuclearNumbers za {type.atomic_number(), type.mass_number()};
  double mass = nuclear_mass_ev(za);
  if (mass > 0.0)
    return mass;
  // Unlisted nuclide: use A u only as an inertial mass. Q-value calculations
  // call nuclear_mass_ev() directly and therefore never take this fallback.
  if (type.is_nucleus() && za.A > 0)
    return za.A * AMU_EV;
  return 0.0;
}

double kalbach_precompound_fraction(double E_cm, double E_max_shape,
  double E_in, NuclearNumbers daughter, const AngularParams& par)
{
  if (daughter.A <= 0)
    return 0.0;
  double x = E_cm / std::max(E_max_shape, 1.0);
  x = std::clamp(x, 0.0, 1.0);
  double a_third = std::pow(static_cast<double>(daughter.A), -1.0 / 3.0);
  double excess = (daughter.A - 2.0 * daughter.Z) / daughter.A;
  double u = par.c0 + par.c1 * x + par.c2 * std::log1p(E_in / 10.0e6) +
             par.c3 * a_third + par.c4 * excess;
  u = std::clamp(u, -60.0, 60.0);
  return 1.0 / (1.0 + std::exp(-u));
}

double light_ion_log_pdf(double E, double E_max, int Z_b, int A_b, int Z_d,
  int A_d, const LightIonParams& par)
{
  if (!(E > 0.0) || E >= E_max || !std::isfinite(E_max))
    return -INFTY;

  double log_p = std::log(E) + par.nu * std::log1p(-E / E_max);
  if (Z_b > 0 && Z_d > 0) {
    double radius = par.r0 * (std::cbrt(static_cast<double>(A_b)) +
                               std::cbrt(static_cast<double>(A_d)));
    double barrier = COULOMB_EV_FM * Z_b * Z_d / radius;

    // Sommerfeld parameter eta = Z_b Z_D alpha sqrt(m_red c^2 / 2E). The
    // E^-1/2 inside the exponent makes the decay constant itself grow as E
    // falls, which is the widening of the barrier at lower energy; a
    // transmission with a fixed diffuseness falls at one rate everywhere and
    // cannot reproduce it. Normalized so that T = 1/2 at E = V_C.
    // Use the A-u inertial fallback for an unlisted mass. These masses enter
    // only through the square root of the reduced mass.
    double m_b = particle_mass_ev(ion_type({Z_b, A_b})) / AMU_EV;
    double m_d = particle_mass_ev(ion_type({Z_d, A_d})) / AMU_EV;
    double m_reduced = m_b * m_d / (m_b + m_d);
    // FINE_STRUCTURE is the *inverse* fine-structure constant in OpenMC
    double pref =
      Z_b * Z_d / FINE_STRUCTURE * std::sqrt(0.5 * m_reduced * AMU_EV);
    double arg =
      par.g * 2.0 * PI *
      (pref / std::sqrt(E) - pref / std::sqrt(std::max(barrier, 1.0)));

    // log T_C = -log(1 + e^arg) = -softplus(arg). The stable form does not
    // evaluate e^arg for a large positive argument. The exponent diverges as
    // E -> 0 and must not be clamped: an artificial floor would flatten the
    // sub-barrier spectrum. Log space preserves relative probabilities even
    // when the density underflows.
    log_p -= softplus(arg);
  }
  return log_p;
}

double sample_light_ion_energy(double E_max, double E_limit, int Z_b, int A_b,
  int Z_d, int A_d, uint64_t* seed, const LightIonParams& par)
{
  if (!(E_max > 0.0) || !std::isfinite(E_max))
    return 0.0;
  if (!(E_limit > 0.0) || E_limit > E_max)
    E_limit = E_max;

  // Tabulate the spectrum in log space, then invert its cumulative. The work
  // is exactly N_TABLE evaluations regardless of how peaked the spectrum is.
  // Piecewise-linear interpolation between nodes gives a controlled error that
  // decreases quadratically with the grid spacing.
  double f[N_TABLE + 1];
  double log_peak = -INFTY;
  const double dE = E_limit / N_TABLE;
  for (int i = 0; i <= N_TABLE; ++i) {
    f[i] = light_ion_log_pdf(dE * i, E_max, Z_b, A_b, Z_d, A_d, par);
    if (f[i] > log_peak)
      log_peak = f[i];
  }
  if (log_peak == -INFTY)
    return 0.0; // the channel carries no probability at all

  // Shifting by the peak before exponentiating is what keeps a spectrum
  // spanning hundreds of decades representable.
  for (int i = 0; i <= N_TABLE; ++i) {
    f[i] = (f[i] == -INFTY) ? 0.0 : std::exp(f[i] - log_peak);
  }

  double cdf[N_TABLE + 1];
  cdf[0] = 0.0;
  for (int i = 1; i <= N_TABLE; ++i) {
    cdf[i] = cdf[i - 1] + 0.5 * (f[i - 1] + f[i]) * dE;
  }
  double total = cdf[N_TABLE];
  if (!(total > 0.0) || !std::isfinite(total))
    return 0.0;

  double xi = prn(seed) * total;
  int k = static_cast<int>(std::lower_bound(cdf, cdf + N_TABLE + 1, xi) - cdf);
  k = std::clamp(k - 1, 0, N_TABLE - 1);

  // Within the interval the density is f_k + s(E - E_k), so the remaining
  // probability r fixes the offset through 0.5 s d^2 + f_k d = r. The
  // rationalized root stays accurate when s is small, where the textbook
  // formula cancels.
  double r = xi - cdf[k];
  double slope = (f[k + 1] - f[k]) / dE;
  double disc = std::max(f[k] * f[k] + 2.0 * slope * r, 0.0);
  double denom = f[k] + std::sqrt(disc);
  double d = (denom > 0.0) ? 2.0 * r / denom : 0.0;

  return std::clamp(dE * k + d, 0.0, E_limit);
}

namespace {

//! Additive nuclear mass of everything not yet emitted
double remaining_system_mass(ParticleType recoil, int n_neutrons,
  const ChargedProducts& ions, bool include_ions)
{
  double mass = particle_mass_ev(recoil);
  assert(mass > 0.0);
  mass += n_neutrons * MASS_NEUTRON_EV;
  if (include_ions) {
    for (int i = 0; i < ions.size; ++i) {
      double ion_mass = particle_mass_ev(ion_type(ions.za[i]));
      assert(ion_mass > 0.0);
      mass += ion_mass;
    }
  }
  return mass;
}

//! Create a recoil secondary particle
bool bank_recoil(Particle& p, double weight, Direction p_recoil, double mass,
  ParticleType type)
{
  if (weight <= 0.0)
    return false;

  double p2 = p_recoil.dot(p_recoil);
  if (!std::isfinite(p2) || p2 <= 0.0)
    return false;

  if (mass <= 0.0 || !std::isfinite(mass))
    return false;

  double E = p2 / (2.0 * mass);
  if (!std::isfinite(E) || E <= 0.0)
    return false;

  double p_magnitude = std::sqrt(p2);
  Direction u_recoil = p_recoil / p_magnitude;
  return p.create_secondary(weight, u_recoil, E, type);
}

//! Emit the light charged particles that the library does not describe
//!
//! Each ion is emitted sequentially in the rest frame of the system that has
//! not yet decayed. Its center-of-mass energy comes from
//! sample_light_ion_energy() with the endpoint of the corresponding pure
//! channel, then is rejected and resampled until it also fits inside this
//! event's remaining budget. Its direction follows Kalbach-Mann systematics.
//!
//! \return false if the exit channel is kinematically impossible, in which case
//!         the caller keeps the momentum balance it already had
bool emit_light_ions(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double E_in, Direction u_in, ChargedProducts& ions, EmissionState& state,
  SampledIon* sampled, int& n_sampled)
{
  n_sampled = 0;
  if (ions.size == 0)
    return true;
  assert(state.mass > 0.0);
  if (state.energy_available <= 0.0)
    return false;

  uint64_t* seed = p.current_seed();

  // Emission order affects which ion sees the larger budget; randomize it so
  // that no ion is systematically favored.
  fisher_yates_shuffle(ions.to_span(), seed);

  for (int i = 0; i < ions.size; ++i) {
    NuclearNumbers b = ions.za[i];
    ParticleType b_type = ion_type(b);
    double m_b = particle_mass_ev(b_type);
    NuclearNumbers d = state.za - b;
    double m_d = state.mass - m_b;
    assert(m_b > 0.0);
    assert(m_d > 0.0);
    assert(d.is_valid());

    // Kinematic endpoints for this event and for the pure channel that emits
    // this ion alone. The latter sets the spectrum shape. Both use the same
    // kinematic contract and remove the entrance-channel translational energy.
    double E_max_event = two_body_endpoint(state.energy_available, m_b, m_d);
    if (E_max_event <= 0.0 || !std::isfinite(E_max_event))
      return false;
    double E_max_shape = shape_endpoint(E_in, {nuc.Z_, nuc.A_}, b);
    // The single-ion channel has the larger Q, so its endpoint normally
    // dominates; the guard covers a missing mass, which returns zero.
    E_max_shape = std::max(E_max_event, E_max_shape);

    double E_cm;
    if (is_discrete_charged_level(rx.mt_)) {
      // MT=600-849 name a single recoil level, so the exit channel is exactly
      // two-body and the light-ion energy is fixed by the reaction Q value.
      E_cm = E_max_event;
    } else {
      E_cm = sample_light_ion_energy(
        E_max_shape, E_max_event, b.Z, b.A, d.Z, d.A, seed);
    }
    if (E_cm <= 0.0 || E_cm > E_max_event || !std::isfinite(E_cm))
      return false;

    double mu;
    if (is_discrete_charged_level(rx.mt_)) {
      // Kalbach's systematics describe continuum pre-equilibrium emission, not
      // a named residual level. With no evaluated charged-particle angle for
      // these two-body channels, sample the center-of-mass direction
      // isotropically.
      mu = 2.0 * prn(seed) - 1.0;
    } else {
      AngularParams ang {};
      double slope =
        ang.slope_scale * kalbach_slope(E_in, E_cm, b, nuc.Z_, nuc.A_);
      double r = kalbach_precompound_fraction(E_cm, E_max_shape, E_in, d, ang);
      mu = sample_kalbach_mu(slope, r, seed);
    }
    Direction u_cm = rotate_angle(u_in, mu, nullptr, seed);

    // Boost from the rest frame of the decaying system into the laboratory
    Direction v_parent = state.momentum / state.mass;
    Direction p_lab = m_b * v_parent + momentum_from_energy(m_b, E_cm, u_cm);
    double p2 = p_lab.dot(p_lab);
    if (p2 <= 0.0 || !std::isfinite(p2))
      return false;
    double E_lab = p2 / (2.0 * m_b);

    sampled[n_sampled++] = {b_type, p_lab / std::sqrt(p2), E_lab};

    state.momentum -= p_lab;
    state.mass = m_d;
    state.za = d;
    state.energy_emitted += E_lab;
    // The two-body decay took E_cm from the ion and E_cm*m_b/m_d from the
    // daughter out of the internal energy budget. E_cm was drawn no larger
    // than the endpoint, so this cannot go negative by more than rounding.
    state.energy_available =
      residual_excitation(state.energy_available, E_cm, m_b, m_d);
    if (state.energy_available < 0.0)
      state.energy_available = 0.0;
  }
  return true;
}

//! Sample an outgoing neutron from a reaction, converting CM to LAB if needed
// Keep this transformation in functional parity with the transported-neutron
// path in inelastic_scatter().
void sample_reaction_neutron(const Nuclide& nuc, const Reaction& rx,
  double E_in, Direction u_in, uint64_t* seed, double& E_out, Direction& u_out)
{
  double mu;
  rx.products_[0].sample(E_in, E_out, mu, seed);

  if (rx.scatter_in_cm_) {
    double E_cm = E_out;
    double A = nuc.awr_;
    E_out = E_cm + (E_in + 2.0 * mu * (A + 1.0) * std::sqrt(E_in * E_cm)) /
                     ((A + 1.0) * (A + 1.0));
    mu = mu * std::sqrt(E_cm / E_out) + std::sqrt(E_in / E_out) / (A + 1.0);
  }
  if (std::abs(mu) > 1.0)
    mu = std::copysign(1.0, mu);
  u_out = rotate_angle(u_in, mu, nullptr, seed);
}

//! Sample the photons of a reaction and return their total momentum and energy
struct PhotonKick {
  Direction momentum {};
  double energy {0.0};
};

//! Scale a sampled capture cascade onto its energy budget
//!
//! Processed reaction data contain inclusive photon distributions rather than
//! an eventwise correlated cascade, so independently sampled photons do not
//! generally close the event energy balance. Apply one common scale factor so
//! that the constructed cascade and recoil satisfy
//!
//! \f[ \sum_i E_i + E_R = E_\text{in} + Q, \f]
//!
//! while preserving the sampled directions and multiplicity. These photons are
//! used only to determine recoil momentum; transported photons are sampled
//! separately. See the recoil methods documentation for the full rationale.
//!
//! \param[in,out] kick    Cascade momentum and energy, scaled in place
//! \param[in] p_in        Momentum carried into the reaction in [eV]
//! \param[in] mass        Mass of the recoiling compound nucleus in [eV]
//! \param[in] budget      Energy available to the exit channel in [eV]
void constrain_photon_kick(
  PhotonKick& kick, Direction p_in, double mass, double budget)
{
  if (kick.energy <= 0.0 || mass <= 0.0 || budget <= 0.0)
    return;

  // Solve lambda^2 |P|^2/(2m) + lambda (S - p.P/m) + |p|^2/(2m) - budget = 0
  double a = kick.momentum.dot(kick.momentum) / (2.0 * mass);
  double b = kick.energy - p_in.dot(kick.momentum) / mass;
  double c = p_in.dot(p_in) / (2.0 * mass) - budget;

  double lambda;
  if (a > 0.0) {
    // The quadratic coefficient is of order the recoil energy and the linear
    // one of order the cascade energy, so the textbook root formula loses
    // essentially all of its significant digits here. Use the factored form.
    double disc = std::max(b * b - 4.0 * a * c, 0.0);
    double q = -0.5 * (b + std::copysign(std::sqrt(disc), b));
    lambda = (q != 0.0) ? c / q : 0.0;
  } else {
    lambda = (b > 0.0) ? -c / b : 0.0;
  }

  if (!(lambda > 0.0) || !std::isfinite(lambda))
    return;
  kick.momentum *= lambda;
  kick.energy *= lambda;
}

PhotonKick sample_photon_kick(
  const Reaction& rx, double E_in, Direction u_in, uint64_t* seed)
{
  PhotonKick kick;
  for (const auto& product : rx.products_) {
    if (!product.particle_.is_photon())
      continue;

    double y = (*product.yield_)(E_in);
    if (y <= 0.0)
      continue;

    int n = static_cast<int>(y);
    if (prn(seed) < y - n)
      ++n;

    for (int i = 0; i < n; ++i) {
      double E_gamma;
      double mu_gamma;
      product.sample(E_in, E_gamma, mu_gamma, seed);
      Direction u_gamma = rotate_angle(u_in, mu_gamma, nullptr, seed);
      if (E_gamma > 0.0)
        kick.momentum += E_gamma * u_gamma;
      kick.energy += E_gamma;
    }
  }
  return kick;
}

//! Finish an event: model the missing light ions and bank all products
//!
//! Either every product of the exit channel is accounted for or nothing is
//! banked. The recoil's identity asserts that a particular set of particles
//! left; banking it while silently omitting one of them makes the secondary
//! particle's mass, charge, and momentum disagree with its own label.
void finish_event(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double weight, double E_in, Direction u_in, ChargedProducts& ions,
  EmissionState& state, ParticleType recoil, double budget,
  bool validate_closure)
{
  SampledIon sampled[ChargedProducts::MAX];
  int n_sampled = 0;

  if (ions.size > 0 &&
      settings::recoil.light_ion_model == RecoilLightIonModel::statistical) {
    EmissionState trial = state;
    if (!emit_light_ions(
          p, nuc, rx, E_in, u_in, ions, trial, sampled, n_sampled)) {
      // The exit channel could not be completed. A partial event would assign
      // the recoil an identity inconsistent with its mass and momentum.
      return;
    }
    state = trial;
  }

  // Pin reconstructed events to the mass and momentum actually banked. An
  // evaluated one-neutron law is exempt because its sampled neutron remains
  // authoritative and the evaluation does not promise eventwise closure.
  if (validate_closure && !update_internal_energy(state, budget)) {
    return;
  }

  if (settings::recoil.emitted_ions) {
    for (int i = 0; i < n_sampled; ++i) {
      p.create_secondary(
        weight, sampled[i].direction, sampled[i].energy, sampled[i].type);
    }
  }
  bank_recoil(p, weight, state.momentum, state.mass, recoil);
}

} // namespace

//==============================================================================
// Entry points
//==============================================================================

void from_elastic(Particle& p, const Nuclide& nuc, double E_in, Direction u_in,
  double E_out, Direction u_out)
{
  // The neutron momentum transfer is the recoil momentum in the incident
  // target's rest frame. The free-gas target velocity shapes this transfer
  // through the sampled outgoing neutron, but the secondary particle is
  // neither the target's final laboratory kinetic energy nor its signed
  // laboratory kinetic-energy change.
  Direction p_recoil =
    neutron_momentum(E_in, u_in) - neutron_momentum(E_out, u_out);
  bank_recoil(
    p, p.wgt(), p_recoil, nuc.awr_ * MASS_NEUTRON_EV, nuc.particle_type());
}

void from_inelastic(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double wgt, double E_in, Direction u_in, double E_out, Direction u_out,
  double yield)
{
  if (!rx.recoil_.supported)
    return;
  const auto& emitted = rx.recoil_.emitted;
  double q = rx.recoil_.q_value;
  ParticleType recoil = rx.recoil_.residual;

  ChargedProducts ions {emitted};

  // Additional neutrons of a multiplicity > 1 channel. OpenMC transports only
  // one sampled neutron, so the others are sampled independently from the same
  // evaluated distribution for recoil kinematics only.
  int n_extra = 0;
  if (emitted.neutron > 1) {
    n_extra = emitted.neutron - 1;
  } else if (yield > 1.0 && std::floor(yield) == yield) {
    n_extra = static_cast<int>(std::round(yield)) - 1;
  }

  // Start with the incident momentum, then remove the transported neutron.
  // The inertial mass is not the entrance compound mass: in this
  // nonrelativistic construction Q is energy, while the mass of the undecayed
  // system is the additive ground-state nuclear mass of its final products.
  EmissionState state;
  state.momentum = neutron_momentum(E_in, u_in);
  state.za = {nuc.Z_, nuc.A_ + 1};

  // Lab-frame kinetic energy release. The closure calculation subtracts the
  // translational energy of the final-state mass that has not yet decayed.
  double budget = E_in + q;

  state.momentum -= neutron_momentum(E_out, u_out);
  state.za -= NEUTRON_NUMBERS;
  state.energy_emitted += E_out;

  bool include_ions = ions.size > 0 && settings::recoil.light_ion_model ==
                                         RecoilLightIonModel::statistical;
  // A one-neutron inelastic law is sampled with the evaluation's target AWR,
  // just like elastic scattering. Preserve that processed-law mass so the
  // recoil remains the exact complement of the evaluated neutron. Reconstructed
  // multiparticle states use the additive nuclear-mass inventory instead.
  bool evaluated_neutron_only =
    emitted.neutron == 1 && n_extra == 0 && ions.size == 0;
  if (evaluated_neutron_only) {
    state.mass = nuc.awr_ * MASS_NEUTRON_EV;
  } else {
    state.mass = remaining_system_mass(recoil, n_extra, ions, include_ions);
  }
  if (!evaluated_neutron_only && !update_internal_energy(state, budget)) {
    // The evaluated first-neutron marginal left no kinematically possible
    // remainder. It is already committed to transport and must not be
    // resampled merely to manufacture a recoil secondary particle.
    return;
  }

  for (int i = 0; i < n_extra; ++i) {
    double E_extra;
    Direction u_extra;
    bool accepted = false;
    for (int attempt = 0; attempt < MAX_BUDGET_TRIES; ++attempt) {
      sample_reaction_neutron(
        nuc, rx, E_in, u_in, p.current_seed(), E_extra, u_extra);
      EmissionState trial = state;
      trial.momentum -= neutron_momentum(E_extra, u_extra);
      trial.mass -= MASS_NEUTRON_EV;
      trial.za -= NEUTRON_NUMBERS;
      trial.energy_emitted += E_extra;
      if (update_internal_energy(trial, budget)) {
        state = trial;
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      // The recoil identity requires this neutron. Do not bank a partial event
      // whose mass number and momentum describe different exit channels.
      return;
    }
  }

  finish_event(p, nuc, rx, wgt, E_in, u_in, ions, state, recoil, budget,
    !evaluated_neutron_only);
}

void from_absorption(Particle& p, int i_nuclide, double weight, double E_in,
  Direction u_in, const Reaction& rx)
{
  if (weight <= 0.0 || !rx.recoil_.supported)
    return;

  const auto& nuc {data::nuclides[i_nuclide]};
  const auto& emitted = rx.recoil_.emitted;
  double q = rx.recoil_.q_value;
  ParticleType recoil = rx.recoil_.residual;

  ChargedProducts ions {emitted};

  EmissionState state;
  state.momentum = neutron_momentum(E_in, u_in);
  state.za = {nuc->Z_, nuc->A_ + 1};

  bool include_ions = ions.size > 0 && settings::recoil.light_ion_model ==
                                         RecoilLightIonModel::statistical;
  state.mass = remaining_system_mass(recoil, 0, ions, include_ions);
  assert(state.mass > 0.0 && std::isfinite(state.mass));

  double budget = E_in + q;

  // Radiative capture: independently sample photon multiplicity, energies, and
  // directions from this reaction's inclusive photon products, then constrain
  // the constructed cascade to the event energy budget.
  if (rx.mt_ == N_GAMMA) {
    PhotonKick kick = sample_photon_kick(rx, E_in, u_in, p.current_seed());
    constrain_photon_kick(kick, state.momentum, state.mass, budget);
    state.momentum -= kick.momentum;
    state.energy_emitted += kick.energy;
  }

  if (!update_internal_energy(state, budget)) {
    return;
  }

  finish_event(p, *nuc, rx, weight, E_in, u_in, ions, state, recoil, budget,
    include_ions || ions.size == 0);
}

} // namespace recoil
} // namespace openmc
