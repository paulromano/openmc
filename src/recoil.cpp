#include "openmc/recoil.h"

#include "openmc/atomic_mass.h"
#include "openmc/constants.h"
#include "openmc/distribution_multi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/reaction.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"

#include <algorithm> // for max, min, swap, lower_bound
#include <cctype>    // for tolower, isdigit
#include <cmath>     // for sqrt, exp, cbrt, log
#include <limits>    // for numeric_limits
#include <string>

namespace openmc {
namespace recoil {

namespace {

//==============================================================================
// Tunable constants of the light-ion emission model
//
// The calibration itself lives in the defaults of LightIonParams and
// AngularParams, so that a candidate parameter set can be scored through this
// same code path. It was fitted to 3708 evaluated centre-of-mass spectra drawn
// from six libraries, 29 nuclides, five light ions and incident energies from
// 2 to 20 MeV.
//==============================================================================

//! Coulomb constant e^2/(4 pi eps0) in [eV fm]
constexpr double COULOMB_EV_FM = 1.44e6;

//! Event counters, so that a skipped event is visible rather than merely absent
//!
//! Every path that declines to bank a recoil increments one of these. The
//! feature's failure mode used to be a record that looked ordinary and was
//! wrong; the failure mode now is a record that is missing, which is only an
//! improvement if it can be counted.
int64_t COUNTERS[static_cast<int>(RecoilCounter::size)] = {};

//! Points in the tabulated cumulative the light-ion sampler inverts
//!
//! Sixty-four intervals put the sampled mean within 0.1% of the analytic one
//! for every ion, daughter and truncation tested, at a fixed cost the previous
//! rejection sampler only matched in its best case.
constexpr int N_TABLE = 64;

//! Negative infinity, for a log density that is identically zero
constexpr double INFTY = std::numeric_limits<double>::infinity();

//! log(1 + e^x), evaluated so that neither limb overflows
double softplus(double x)
{
  return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x)));
}

//! Maximum attempts to resample a modelled product inside the energy budget
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

Direction photon_momentum(double E, Direction u)
{
  return E <= 0.0 ? Direction {} : E * u;
}

//==============================================================================
// Exit-channel bookkeeping
//==============================================================================

//! Number of each light particle emitted by a reaction
struct EmittedParticles {
  int neutron {0};
  int proton {0};
  int deuteron {0};
  int triton {0};
  int he3 {0};
  int alpha {0};
};

bool parse_channel(std::string channel, EmittedParticles& out)
{
  for (auto& c : channel) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  if (channel == "gamma")
    return true;

  size_t i = 0;
  while (i < channel.size()) {
    // "3he" has to be recognized before the multiplicity digits, or the
    // leading 3 of helium-3 is consumed as a multiplicity and the rest of the
    // channel no longer parses. Helium-3 never appears with a multiplicity
    // above one in the ENDF reaction names.
    if (channel.compare(i, 3, "3he") == 0) {
      ++out.he3;
      i += 3;
      continue;
    }

    int multiplicity = 0;
    while (i < channel.size() &&
           std::isdigit(static_cast<unsigned char>(channel[i]))) {
      multiplicity = 10 * multiplicity + (channel[i] - '0');
      ++i;
    }
    if (multiplicity == 0)
      multiplicity = 1;
    if (i >= channel.size())
      return false;

    if (channel.compare(i, 3, "3he") == 0) {
      out.he3 += multiplicity;
      i += 3;
      continue;
    }
    if (channel.compare(i, 5, "gamma") == 0) {
      i += 5;
      continue;
    }

    switch (channel[i]) {
    case 'n':
      out.neutron += multiplicity;
      break;
    case 'p':
      out.proton += multiplicity;
      break;
    case 'd':
      out.deuteron += multiplicity;
      break;
    case 't':
      out.triton += multiplicity;
      break;
    case 'a':
      out.alpha += multiplicity;
      break;
    default:
      return false;
    }
    ++i;
  }
  return true;
}

//! Determine the light particles emitted by a reaction from its MT number
bool emitted_particles(int mt, EmittedParticles& out)
{
  out = EmittedParticles {};

  if (mt == ELASTIC || (mt >= N_N1 && mt <= N_NC)) {
    out.neutron = 1;
    return true;
  }
  if (mt >= N_P0 && mt <= N_PC) {
    out.proton = 1;
    return true;
  }
  if (mt >= N_D0 && mt <= N_DC) {
    out.deuteron = 1;
    return true;
  }
  if (mt >= N_T0 && mt <= N_TC) {
    out.triton = 1;
    return true;
  }
  if (mt >= N_3HE0 && mt <= N_3HEC) {
    out.he3 = 1;
    return true;
  }
  if (mt >= N_A0 && mt <= N_AC) {
    out.alpha = 1;
    return true;
  }
  if (mt >= N_2N0 && mt <= N_2NC) {
    out.neutron = 2;
    return true;
  }

  // Fall back on parsing the "(n,...)" reaction name
  std::string name = reaction_name(mt);
  if (name.size() < 5 || name.compare(0, 3, "(n,") != 0 || name.back() != ')')
    return false;
  return parse_channel(name.substr(3, name.size() - 4), out);
}

AtomicNumbers particle_za(ParticleType type)
{
  switch (type.pdg_number()) {
  case PDG_NEUTRON:
    return {0, 1};
  case PDG_PROTON:
    return {1, 1};
  case PDG_DEUTERON:
    return {1, 2};
  case PDG_TRITON:
    return {1, 3};
  case PDG_ALPHA:
    return {2, 4};
  default:
    if (type.is_nucleus()) {
      int pdg = type.pdg_number();
      return {(pdg / 10000) % 1000, (pdg / 10) % 1000};
    }
    return {};
  }
}

ParticleType ion_type(AtomicNumbers za)
{
  if (za.Z == 0 && za.A == 1)
    return ParticleType::neutron();
  if (za.Z == 1 && za.A == 1)
    return ParticleType::proton();
  if (za.Z == 1 && za.A == 2)
    return ParticleType::deuteron();
  if (za.Z == 1 && za.A == 3)
    return ParticleType::triton();
  return ParticleType {za.Z, za.A, 0};
}

//! Push the charged ions of an exit channel onto a fixed-capacity list
struct ChargedProducts {
  static constexpr int MAX = 8;
  AtomicNumbers za[MAX];
  int n {0};

  void add(AtomicNumbers value, int count)
  {
    for (int i = 0; i < count && n < MAX; ++i) {
      za[n++] = value;
    }
  }

  void fill(const EmittedParticles& e)
  {
    add({1, 1}, e.proton);
    add({1, 2}, e.deuteron);
    add({1, 3}, e.triton);
    add({2, 3}, e.he3);
    add({2, 4}, e.alpha);
  }
};

//! Fill \p out with the emitted particles of a channel, light ions and all
//!
//! The complement of ChargedProducts::fill(), which keeps only the charged
//! ones. A mass budget needs every particle that leaves.
int all_products(const EmittedParticles& e, AtomicNumbers* out, int capacity)
{
  const AtomicNumbers kinds[6] = {
    {0, 1}, {1, 1}, {1, 2}, {1, 3}, {2, 3}, {2, 4}};
  const int counts[6] = {
    e.neutron, e.proton, e.deuteron, e.triton, e.he3, e.alpha};
  int n = 0;
  for (int k = 0; k < 6; ++k) {
    for (int i = 0; i < counts[k] && n < capacity; ++i) {
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
// reproduces the qualitative behaviour of evaluated r values, which rise from
// nearly zero at low outgoing energy to 0.5-0.9 near the kinematic maximum.
//==============================================================================

double breakup_energy(AtomicNumbers p)
{
  if (p.Z == 1 && p.A == 2)
    return 2.224566;
  if (p.Z == 1 && p.A == 3)
    return 8.481798;
  if (p.Z == 2 && p.A == 3)
    return 7.718043;
  if (p.Z == 2 && p.A == 4)
    return 28.29566;
  return 0.0;
}

//! Kalbach's semi-empirical separation energy in [MeV]
double separation_energy(
  AtomicNumbers compound, AtomicNumbers nucleus, AtomicNumbers particle)
{
  double A_c = compound.A;
  double Z_c = compound.Z;
  double N_c = compound.A - compound.Z;
  double A_a = nucleus.A;
  double Z_a = nucleus.Z;
  double N_a = nucleus.A - nucleus.Z;

  return 15.68 * (A_c - A_a) -
         28.07 *
           ((N_c - Z_c) * (N_c - Z_c) / A_c - (N_a - Z_a) * (N_a - Z_a) / A_a) -
         18.56 * (std::pow(A_c, 2.0 / 3.0) - std::pow(A_a, 2.0 / 3.0)) +
         33.22 * ((N_c - Z_c) * (N_c - Z_c) / std::pow(A_c, 4.0 / 3.0) -
                   (N_a - Z_a) * (N_a - Z_a) / std::pow(A_a, 4.0 / 3.0)) -
         0.717 * (Z_c * Z_c / std::cbrt(A_c) - Z_a * Z_a / std::cbrt(A_a)) +
         1.211 * (Z_c * Z_c / A_c - Z_a * Z_a / A_a) - breakup_energy(particle);
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
//! sum of the modelled kinetic energies can never exceed the budget.
struct EmissionState {
  Direction momentum {};    //!< lab momentum of the undecayed system [eV]
  double mass {0.0};        //!< mass of the undecayed system [eV]
  double internal {0.0};    //!< kinetic energy available in its rest frame [eV]
  AtomicNumbers za {};      //!< charge and mass number of the undecayed system
  double emitted_kin {0.0}; //!< lab kinetic energy already given to products
};

//! Record of a modelled light ion, kept until we know whether to bank it
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
  double scale =
    std::max({1.0, std::abs(budget), std::abs(state.emitted_kin), translation});
  double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  double internal = remaining_internal_energy(
    budget, state.emitted_kin, state.momentum, state.mass);
  if (!std::isfinite(internal) || internal < -tolerance)
    return false;
  state.internal = std::max(0.0, internal);
  return true;
}

//! True for the discrete charged-particle level MTs, which are exactly two-body
bool is_discrete_charged_level(int mt)
{
  return (mt >= N_P0 && mt < N_PC) || (mt >= N_D0 && mt < N_DC) ||
         (mt >= N_T0 && mt < N_TC) || (mt >= N_3HE0 && mt < N_3HEC) ||
         (mt >= N_A0 && mt < N_AC);
}

//! True for every MT that names one residual level, charged or neutron
//!
//! These are the MTs whose evaluated \c QI is a trustworthy energy release,
//! because the MT identifies the state the residual is left in. For anything
//! else \c QI is a threshold-setting value; see event_q().
bool names_one_level(int mt)
{
  return (mt >= N_N1 && mt < N_NC) || is_discrete_charged_level(mt);
}

//! How far an evaluated level Q may sit above the ground-state mass budget
//!
//! A level Q above the ground-state budget means a negative excitation, which
//! is unphysical. Small excesses are mass-table disagreement: across 418
//! evaluated channels in four libraries, MF=3 \c QM departs from the AME2020
//! mass difference by at most 90 keV. This threshold is comfortably above that
//! and far below a real inconsistency, so an excess inside it is absorbed by
//! capping and an excess beyond it means the evaluation cannot be reconciled.
constexpr double Q_LEVEL_TOLERANCE = 0.25e6;

//! Rest-mass energy release available to an event, in [eV]
//!
//! ENDF stores two Q values and OpenMC keeps only one of them. \c Reaction::
//! q_value_ is MF=3 \c QI, the Q of the lowest state the MT represents, or an
//! effective value chosen to put the threshold in the right place when the MT
//! names no unique state. For a discrete level that is exactly the budget
//! wanted. For a continuum, level range or summation channel it is not: over
//! the 190 split-representation channels of the calibration libraries it lies
//! a median 3.5 keV and up to 7.0 MeV *below* the mass-difference Q, and using
//! it truncates the modelled spectrum inside the evaluated one.
//!
//! \param[in] nuc      Target nuclide
//! \param[in] rx       Reaction that was sampled
//! \param[in] emitted  Light particles the exit channel produces
//! \param[out] ok      False when no trustworthy budget can be formed, in
//!                     which case the caller must produce nothing rather than
//!                     bank a recoil built on a budget it does not believe
double event_q(const Nuclide& nuc, const Reaction& rx,
  const EmittedParticles& emitted, bool& ok)
{
  ok = true;
  AtomicNumbers target {nuc.Z_, nuc.A_};
  AtomicNumbers products[ChargedProducts::MAX + 4];
  int n = all_products(emitted, products, ChargedProducts::MAX + 4);

  bool have_masses = false;
  double q_mass = mass_difference_q(target, products, n, have_masses);

  if (!names_one_level(rx.mt_)) {
    if (have_masses)
      return q_mass;
    // No tabulated mass for this daughter. QI is then the only budget on
    // offer; it is the right one for a lumped channel, where QI and QM agree
    // in every evaluation examined, and too small for a split continuum.
    return rx.q_value_;
  }

  // A named level: QI is the level-specific Q, which is what a two-body
  // channel needs. Check it against the ground-state budget, since the two
  // must differ by the level excitation and that cannot be negative.
  if (!have_masses)
    return rx.q_value_;
  if (rx.q_value_ <= q_mass)
    return rx.q_value_;
  if (rx.q_value_ <= q_mass + Q_LEVEL_TOLERANCE)
    return q_mass;
  ok = false;
  return 0.0;
}

void count(RecoilCounter which)
{
  int64_t& slot = COUNTERS[static_cast<int>(which)];
#pragma omp atomic update
  ++slot;
}

} // namespace

int64_t counter(RecoilCounter which)
{
  return COUNTERS[static_cast<int>(which)];
}

void reset_counters()
{
  for (auto& c : COUNTERS)
    c = 0;
}

//==============================================================================
// Kinematic contract
//
// Shared expressions, documented in recoil.h, mirrored by recoil.kinematics in
// the analysis repository and checked against a shared fixture.
//==============================================================================

double nuclear_mass_ev(AtomicNumbers za)
{
  // The five light ions come from their CODATA values. Reading them out of
  // ATOMIC_MASS would work for four of them and quietly fail for the fifth:
  // that table holds bare nuclear masses for the proton, deuteron, helion and
  // alpha, which are CODATA particle constants, but the *atomic* H-3 mass at
  // PDG 1000010030, which is 0.55 mu heavier than the triton.
  if (za.Z == 0 && za.A == 1)
    return MASS_NEUTRON * AMU_EV;
  if (za.Z == 1 && za.A == 1)
    return MASS_PROTON * AMU_EV;
  if (za.Z == 1 && za.A == 2)
    return MASS_DEUTRON * AMU_EV;
  if (za.Z == 1 && za.A == 3)
    return 3.01550071621 * AMU_EV; // triton, nuclear
  if (za.Z == 2 && za.A == 3)
    return MASS_HELION * AMU_EV;
  if (za.Z == 2 && za.A == 4)
    return MASS_ALPHA * AMU_EV;

  if (za.Z < 0 || za.A <= 0 || za.Z > za.A)
    return 0.0;
  int32_t pdg = 1000000000 + za.Z * 10000 + za.A * 10;
  auto it = ATOMIC_MASS.find(pdg);
  if (it == ATOMIC_MASS.end())
    return 0.0;
  // Atomic minus Z electrons. Electron *binding* energy is neglected; it does
  // not cancel exactly in a charged-particle Q value, but the residue is a few
  // keV on a mid-mass target against a budget of MeV.
  return (it->second - za.Z * MASS_ELECTRON) * AMU_EV;
}

double mass_excess_ev(AtomicNumbers za)
{
  double m = nuclear_mass_ev(za);
  return m <= 0.0 ? 0.0 : m - za.A * AMU_EV;
}

double mass_difference_q(
  AtomicNumbers target, const AtomicNumbers* emitted, int n_emitted, bool& ok)
{
  ok = false;
  AtomicNumbers daughter {target.Z, target.A + 1};
  // Mass excesses rather than masses. A W-184 channel differences four numbers
  // of order 1.7e11 eV to reach one of order 1e6, and a double carries about
  // sixteen digits, so the direct subtraction returns a Q good to only ten.
  // The mass numbers cancel identically because A_T + 1 = A_D + sum A_j, so
  // working in excesses of order 1e7 eV gives the same Q to fifteen.
  double exit_excess = 0.0;
  for (int i = 0; i < n_emitted; ++i) {
    daughter.Z -= emitted[i].Z;
    daughter.A -= emitted[i].A;
    if (nuclear_mass_ev(emitted[i]) <= 0.0)
      return 0.0;
    exit_excess += mass_excess_ev(emitted[i]);
  }
  if (daughter.Z < 0 || daughter.A <= 0 || daughter.Z > daughter.A)
    return 0.0;
  if (nuclear_mass_ev(target) <= 0.0 || nuclear_mass_ev(daughter) <= 0.0)
    return 0.0;

  ok = true;
  return mass_excess_ev(target) + mass_excess_ev({0, 1}) -
         mass_excess_ev(daughter) - exit_excess;
}

double entrance_internal_energy(double E_in, double m_target, double q)
{
  double m_n = MASS_NEUTRON * AMU_EV;
  if (m_target <= 0.0)
    return E_in + q;
  return E_in * m_target / (m_target + m_n) + q;
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
  double budget, double emitted_kin, Direction momentum, double mass)
{
  if (mass <= 0.0 || !std::isfinite(mass))
    return -INFTY;
  return budget - emitted_kin - momentum.dot(momentum) / (2.0 * mass);
}

double shape_endpoint(double E_in, AtomicNumbers target, AtomicNumbers ion)
{
  bool ok = false;
  double q = mass_difference_q(target, &ion, 1, ok);
  if (!ok)
    return 0.0;
  AtomicNumbers daughter {target.Z + 0 - ion.Z, target.A + 1 - ion.A};
  double u = entrance_internal_energy(E_in, nuclear_mass_ev(target), q);
  return two_body_endpoint(u, nuclear_mass_ev(ion), nuclear_mass_ev(daughter));
}

//==============================================================================
// Public helpers
//==============================================================================

double kalbach_slope(
  double E_in, double E_cm, AtomicNumbers emitted, int Z_t, int A_t)
{
  AtomicNumbers target {Z_t, A_t};
  AtomicNumbers compound {target.Z, target.A + 1};
  AtomicNumbers recoil {compound.Z - emitted.Z, compound.A - emitted.A};
  if (recoil.Z < 0 || recoil.A <= 0 || recoil.Z > recoil.A)
    return 0.0;

  double epsilon_a = E_in * target.A / (target.A + 1.0) / 1.0e6;
  double epsilon_b = E_cm * (recoil.A + emitted.A) / (recoil.A * 1.0e6);
  double e_a = epsilon_a + separation_energy(compound, target, {0, 1});
  double e_b = epsilon_b + separation_energy(compound, recoil, emitted);
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

double kalbach_slope(
  double E_in, double E_cm, AtomicNumbers emitted, const Nuclide& nuc)
{
  return kalbach_slope(E_in, E_cm, emitted, nuc.Z_, nuc.A_);
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
    return std::min(1.0, std::max(-1.0, mu));
  }
  // symmetric component proportional to cosh(a mu)
  double t = (2.0 * xi - 1.0) * std::sinh(slope);
  double mu = std::log(t + std::sqrt(t * t + 1.0)) / slope;
  return std::min(1.0, std::max(-1.0, mu));
}

double particle_mass_ev(ParticleType type)
{
  if (type.is_photon())
    return 0.0;
  int32_t pdg = std::abs(type.pdg_number());
  auto it = ATOMIC_MASS.find(pdg);
  if (it != ATOMIC_MASS.end())
    return it->second * AMU_EV;
  // Unlisted nuclide: fall back on the mass number
  if (type.is_nucleus())
    return particle_za(type).A * AMU_EV;
  return 0.0;
}

ParticleType recoil_particle_type(const Nuclide& nuc, int mt)
{
  EmittedParticles e;
  if (!emitted_particles(mt, e))
    return nuc.particle_type();

  int emitted_A = e.neutron + e.proton + 2 * e.deuteron + 3 * e.triton +
                  3 * e.he3 + 4 * e.alpha;
  int emitted_Z = e.proton + e.deuteron + e.triton + 2 * e.he3 + 2 * e.alpha;

  int Z_res = nuc.Z_ - emitted_Z;
  int A_res = nuc.A_ + 1 - emitted_A;
  if (Z_res <= 0 || A_res <= 0 || Z_res > A_res)
    return nuc.particle_type();
  return ParticleType {Z_res, A_res, 0};
}

double kalbach_precompound_fraction(double E_cm, double E_max_shape,
  double E_in, AtomicNumbers daughter, const AngularParams& par)
{
  if (daughter.A <= 0)
    return 0.0;
  double x = E_cm / std::max(E_max_shape, 1.0);
  x = std::min(1.0, std::max(0.0, x));
  double a_third = std::pow(static_cast<double>(daughter.A), -1.0 / 3.0);
  double excess = (daughter.A - 2.0 * daughter.Z) / daughter.A;
  double u = par.c0 + par.c1 * x + par.c2 * std::log1p(E_in / 10.0e6) +
             par.c3 * a_third + par.c4 * excess;
  u = std::min(60.0, std::max(-60.0, u));
  return 1.0 / (1.0 + std::exp(-u));
}

//! Nuclear mass in [amu] for the Gamow reduced mass
//!
//! Deferred to the kinematic contract rather than tabulated separately here,
//! which is how this came to hold a proton mass differing from the contract's
//! in its eighth digit and a daughter mass of A neutron masses. Neither
//! mattered physically -- the reduced mass moved by 1e-4 -- but two mass
//! tables in one file is one too many, and the cross-language fixture caught
//! the disagreement.
//!
//! Falls back on the mass number when the nuclide is not tabulated, which is
//! acceptable here and only here: the reduced mass enters through a square
//! root and a heavy daughter's contribution to it is already saturated.
double gamow_mass_amu(AtomicNumbers za)
{
  double m = nuclear_mass_ev(za);
  return m > 0.0 ? m / AMU_EV : static_cast<double>(za.A);
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

    // Sommerfeld parameter eta = Z_b Z_D alpha sqrt(mu c^2 / 2E). The E^-1/2
    // inside the exponent makes the decay constant itself grow as E falls,
    // which is the widening of the barrier at lower energy; a transmission
    // with a fixed diffuseness falls at one rate everywhere and cannot
    // reproduce it. Normalized so that T = 1/2 at E = V_C.
    double m_b = gamow_mass_amu({Z_b, A_b});
    double m_d = gamow_mass_amu({Z_d, A_d});
    double mu = m_b * m_d / (m_b + m_d);
    // FINE_STRUCTURE is the *inverse* fine-structure constant in OpenMC
    double pref = Z_b * Z_d / FINE_STRUCTURE * std::sqrt(0.5 * mu * AMU_EV);
    double arg =
      par.g * 2.0 * PI *
      (pref / std::sqrt(E) - pref / std::sqrt(std::max(barrier, 1.0)));

    // log T_C = -log(1 + e^arg) = -softplus(arg), evaluated so that neither
    // limb overflows. There is no bound on the exponent here, and its absence
    // is the point: the exponent diverges as E -> 0, and clamping it at 60 --
    // which the direct form needed, because 1/(1+e^arg) underflows to exactly
    // zero and leaves nothing to normalize -- floors the transmission at a
    // constant and flattens the spectrum to E (1-E/E_max)^nu wherever a
    // channel lies deep below the barrier. That erased the barrier shape from
    // the He-3 channels almost entirely. In log space the relative
    // probabilities survive to any depth.
    log_p -= softplus(arg);
  }
  return log_p;
}

double light_ion_log_pdf(
  double E, double E_max, int Z_b, int A_b, int Z_d, int A_d)
{
  return light_ion_log_pdf(E, E_max, Z_b, A_b, Z_d, A_d, LightIonParams {});
}

double light_ion_pdf(double E, double E_max, int Z_b, int A_b, int Z_d, int A_d,
  const LightIonParams& par)
{
  double log_p = light_ion_log_pdf(E, E_max, Z_b, A_b, Z_d, A_d, par);
  return log_p == -INFTY ? 0.0 : std::exp(log_p);
}

double light_ion_pdf(double E, double E_max, int Z_b, int A_b, int Z_d, int A_d)
{
  return light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d, LightIonParams {});
}

double sample_light_ion_energy(double E_max, double E_limit, int Z_b, int A_b,
  int Z_d, int A_d, uint64_t* seed, const LightIonParams& par)
{
  if (!(E_max > 0.0) || !std::isfinite(E_max))
    return 0.0;
  if (!(E_limit > 0.0) || E_limit > E_max)
    E_limit = E_max;

  // Tabulate the spectrum in log space, then invert its cumulative.
  //
  // The previous sampler scanned 33 points for a maximum, inflated it by 10%,
  // and rejected up to 200 times; if all attempts failed it returned the scan
  // maximizer, putting a point mass into the distribution. Nothing guaranteed
  // that the inflated scan maximum bounded the true one, and a broad check
  // that it happened to for the deployed parameters is not a property of the
  // method. Inversion has neither problem: the work is exactly N_TABLE
  // evaluations however peaked the spectrum is, there is no envelope to be
  // wrong about, and there is no fallback.
  //
  // What it approximates instead is the shape between table points, by a
  // straight line. That error is controlled and shrinks as N^-2, where a point
  // mass is not controlled at all.
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
  k = std::min(std::max(k - 1, 0), N_TABLE - 1);

  // Within the interval the density is f_k + s(E - E_k), so the remaining
  // probability r fixes the offset through 0.5 s d^2 + f_k d = r. The
  // rationalized root stays accurate when s is small, where the textbook
  // formula cancels.
  double r = xi - cdf[k];
  double slope = (f[k + 1] - f[k]) / dE;
  double disc = std::max(f[k] * f[k] + 2.0 * slope * r, 0.0);
  double denom = f[k] + std::sqrt(disc);
  double d = (denom > 0.0) ? 2.0 * r / denom : 0.0;

  return std::min(std::max(dE * k + d, 0.0), E_limit);
}

double sample_light_ion_energy(double E_max, double E_limit, int Z_b, int A_b,
  int Z_d, int A_d, uint64_t* seed)
{
  return sample_light_ion_energy(
    E_max, E_limit, Z_b, A_b, Z_d, A_d, seed, LightIonParams {});
}

namespace {

//! Create the recoil production record
bool bank_recoil(Particle& p, const Nuclide& nuc, double weight,
  Direction p_recoil, ParticleType type)
{
  if (weight <= 0.0)
    return false;

  double p2 = p_recoil.dot(p_recoil);
  if (!std::isfinite(p2) || p2 <= 0.0)
    return false;

  double mass = particle_mass_ev(type);
  if (mass <= 0.0 || !std::isfinite(mass))
    mass = nuc.awr_ * MASS_NEUTRON_EV;
  if (mass <= 0.0 || !std::isfinite(mass))
    return false;

  double E = p2 / (2.0 * mass);
  if (!std::isfinite(E) || E <= 0.0)
    return false;

  return p.create_secondary(weight, p_recoil / std::sqrt(p2), E, type);
}

//! Emit the light charged particles that the library does not describe
//!
//! Each ion is emitted sequentially in the rest frame of the system that has
//! not yet decayed. Its centre-of-mass energy comes from
//! sample_light_ion_energy() with the endpoint of the corresponding pure
//! channel, then is rejected and resampled until it also fits inside this
//! event's remaining budget. Its direction follows Kalbach-Mann systematics.
//!
//! \return false if the exit channel is kinematically impossible, in which case
//!         the caller keeps the momentum balance it already had
bool emit_light_ions(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double E_in, Direction u_in, const ChargedProducts& ions,
  EmissionState& state, SampledIon* sampled, int& n_sampled)
{
  n_sampled = 0;
  if (ions.n == 0)
    return true;
  if (state.mass <= 0.0 || state.internal <= 0.0)
    return false;

  uint64_t* seed = p.current_seed();

  // Emission order affects which ion sees the larger budget; randomize it so
  // that no ion is systematically favoured.
  AtomicNumbers order[ChargedProducts::MAX];
  for (int i = 0; i < ions.n; ++i)
    order[i] = ions.za[i];
  for (int i = ions.n - 1; i > 0; --i) {
    int j = static_cast<int>(uniform_int_distribution(0, i, seed));
    std::swap(order[i], order[j]);
  }

  for (int i = 0; i < ions.n; ++i) {
    AtomicNumbers b = order[i];
    ParticleType b_type = ion_type(b);
    double m_b = particle_mass_ev(b_type);
    AtomicNumbers d {state.za.Z - b.Z, state.za.A - b.A};
    double m_d = state.mass - m_b;
    if (m_b <= 0.0 || m_d <= 0.0 || d.A <= 0 || d.Z < 0 || d.Z > d.A)
      return false;

    // Kinematic endpoint allowed by this event and by the pure channel that
    // emits this ion alone. The latter sets the spectrum shape. Both go
    // through the kinematic contract, so both remove the translational energy
    // of the entrance channel; the shape endpoint used to keep it, which made
    // the distribution the transport kernel sampled differ from the one the
    // calibration was fitted to.
    double E_max_event = two_body_endpoint(state.internal, m_b, m_d);
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
      // Kalbach's systematics describe a continuum channel fed by
      // pre-equilibrium emission, and carrying them onto a named level makes
      // the recoil distribution worse than assuming nothing: against 25,000
      // evaluated discrete distributions the systematics double the recoil
      // Wasserstein error and bias the first Legendre moment by +0.12, while
      // isotropy leaves it at -0.04. Nothing fitted on top of isotropy
      // survived being held out.
      mu = 2.0 * prn(seed) - 1.0;
    } else {
      AngularParams ang {};
      double slope = ang.slope_scale * kalbach_slope(E_in, E_cm, b, nuc);
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
    state.emitted_kin += E_lab;
    // The two-body decay took E_cm from the ion and E_cm*m_b/m_d from the
    // daughter out of the internal energy budget. E_cm was drawn no larger
    // than the endpoint, so this cannot go negative by more than rounding.
    state.internal = residual_excitation(state.internal, E_cm, m_b, m_d);
    if (state.internal < 0.0)
      state.internal = 0.0;
  }
  return true;
}

//! Sample an outgoing neutron from a reaction, converting CM to LAB if needed
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
//! Radiative capture leaves no massive light ion, so the recoil is kicked only
//! by the photons, and the compound nucleus de-excites all the way to the
//! ground state. The cascade therefore carries the *whole* excitation energy,
//!
//! \f[ \sum_i E_i + E_R = E_\text{in} + Q, \f]
//!
//! an equality rather than a bound. The evaluated files supply only an average
//! photon spectrum and an average multiplicity, never the joint distribution of
//! a cascade, so photons drawn independently from that spectrum do not obey it:
//! their sum spans more than an order of magnitude about the budget and exceeds
//! it in roughly 45% of events at 14 MeV, which inflates the width of the
//! recoil spectrum by half.
//!
//! The cascade is therefore drawn first and then scaled by the single factor
//! that satisfies the equality. Scaling rather than rejecting keeps the
//! evaluated multiplicity exactly -- rejecting whole cascades would bias it
//! low, because one with more photons is likelier to overshoot -- and keeps the
//! mean photon energy within a couple of percent, where rejection would lose
//! nearly half of it. What it gives up is the spread of individual photon
//! energies, which is the right thing to give up here: these photons exist only
//! to build the recoil momentum, are never transported or tallied, and the
//! cascade that is actually transported is sampled separately in
//! sample_secondary_photons().
//!
//! \param[in,out] kick    Cascade momentum and energy, scaled in place
//! \param[in] p_in        Momentum carried into the reaction in [sqrt(amu eV)]
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
      kick.momentum += photon_momentum(E_gamma, u_gamma);
      kick.energy += E_gamma;
    }
  }
  return kick;
}

//! Finish an event: model the missing light ions and bank all products
//!
//! Either every product of the exit channel is accounted for or nothing is
//! banked. The recoil's identity asserts that a particular set of particles
//! left; banking it while silently omitting one of them makes the record's
//! mass, charge and momentum disagree with its own label.
void finish_event(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double weight, double E_in, Direction u_in, const ChargedProducts& ions,
  EmissionState& state, ParticleType recoil, double budget,
  bool enforce_closure)
{
  SampledIon sampled[ChargedProducts::MAX];
  int n_sampled = 0;

  if (ions.n > 0 &&
      settings::recoil.light_ion_model == RecoilLightIonModel::statistical) {
    EmissionState trial = state;
    if (!emit_light_ions(
          p, nuc, rx, E_in, u_in, ions, trial, sampled, n_sampled)) {
      // The exit channel could not be completed. Banking the recoil anyway --
      // which is what this did -- leaves a record labelled with a nuclide that
      // is lighter than the momentum balance it carries by exactly the ions
      // that were dropped.
      count(RecoilCounter::incomplete_emission);
      return;
    }
    state = trial;
  }

  // Multi-neutron candidates are constrained after every emission, and this
  // final check pins the invariant to the mass and momentum actually banked.
  if (enforce_closure && !update_internal_energy(state, budget)) {
    count(RecoilCounter::incomplete_emission);
    return;
  }

  if (settings::recoil.emitted_ions) {
    for (int i = 0; i < n_sampled; ++i) {
      p.create_secondary(
        weight, sampled[i].direction, sampled[i].energy, sampled[i].type);
    }
  }
  count(bank_recoil(p, nuc, weight, state.momentum, recoil)
          ? RecoilCounter::banked
          : RecoilCounter::unbankable);
}

} // namespace

//==============================================================================
// Entry points
//==============================================================================

void from_elastic(Particle& p, const Nuclide& nuc, double E_in, Direction u_in,
  double E_out, Direction u_out)
{
  if (!settings::recoil_production)
    return;

  // The momentum the collision transferred to the target, deliberately not
  // including the target's own thermal momentum. A PKA energy is the energy a
  // collision *imparts* to an atom, which is the quantity displacement models
  // and NJOY recoil matrices are built on; adding the target's pre-collision
  // momentum would instead report its total kinetic energy, which below a few
  // eV is dominated by thermal motion and has nothing to do with damage. The
  // free-gas treatment still shapes the result through the sampled outgoing
  // neutron.
  Direction p_recoil =
    neutron_momentum(E_in, u_in) - neutron_momentum(E_out, u_out);
  count(bank_recoil(p, nuc, p.wgt(), p_recoil, nuc.particle_type())
          ? RecoilCounter::banked
          : RecoilCounter::unbankable);
}

void from_inelastic(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double wgt, double E_in, Direction u_in, double E_out, Direction u_out,
  double yield)
{
  if (!settings::recoil_production)
    return;

  EmittedParticles emitted;
  if (!emitted_particles(rx.mt_, emitted)) {
    // The exit channel cannot be determined from the MT number, so neither can
    // the identity of the recoil. MT=5 (n,misc) is the case that matters: it
    // is a catch-all with inclusive product yields and no single recoil, and
    // some evaluations put a substantial part of the charged-particle
    // production there. Producing nothing is better than producing a record
    // labeled with the wrong nuclide.
    count(RecoilCounter::unknown_channel);
    return;
  }

  bool budget_ok = false;
  double q = event_q(nuc, rx, emitted, budget_ok);
  if (!budget_ok) {
    // The evaluated level Q is above the ground-state mass budget by more than
    // any mass table disagrees, so the event has no trustworthy energy
    // release. Producing nothing is better than banking a recoil built on a
    // budget we do not believe.
    count(RecoilCounter::no_budget);
    return;
  }

  ParticleType recoil = recoil_particle_type(nuc, rx.mt_);

  ChargedProducts ions;
  ions.fill(emitted);

  // Additional neutrons of a multiplicity > 1 channel. OpenMC transports only
  // one sampled neutron, so the others are sampled independently from the same
  // evaluated distribution for recoil kinematics only.
  int n_extra = 0;
  if (emitted.neutron > 1) {
    n_extra = emitted.neutron - 1;
  } else if (yield > 1.0 && std::floor(yield) == yield) {
    n_extra = static_cast<int>(std::round(yield)) - 1;
  }

  // Start from the compound system and remove the transported neutron
  EmissionState state;
  state.momentum = neutron_momentum(E_in, u_in);
  state.mass = MASS_NEUTRON_EV + nuc.awr_ * MASS_NEUTRON_EV;
  state.za = {nuc.Z_, nuc.A_ + 1};
  state.emitted_kin = 0.0;

  // Lab-frame kinetic energy release. The closure calculation subtracts the
  // translational energy of whatever has not decayed, so its first evaluation
  // is exactly the entrance internal energy of the contract,
  // E_in M_T/(M_T + m_n) + Q.
  double budget = E_in + q;

  state.momentum -= neutron_momentum(E_out, u_out);
  state.mass -= MASS_NEUTRON_EV;
  state.za.A -= 1;
  state.emitted_kin += E_out;

  if (n_extra > 0) {
    // Use the same residual mass here and in bank_recoil(). The old
    // compound-minus-products mass differed by binding energies, allowing the
    // acceptance predicate and the banked recoil energy to disagree.
    state.mass = particle_mass_ev(recoil) + n_extra * MASS_NEUTRON_EV;
    for (int i = 0; i < ions.n; ++i) {
      state.mass += particle_mass_ev(ion_type(ions.za[i]));
    }
    if (!update_internal_energy(state, budget)) {
      // The evaluated first-neutron marginal left no kinematically possible
      // remainder. It is already committed to transport and must not be
      // resampled merely to manufacture a recoil record.
      count(RecoilCounter::infeasible_primary);
      count(RecoilCounter::incomplete_emission);
      return;
    }
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
      trial.za.A -= 1;
      trial.emitted_kin += E_extra;
      if (update_internal_energy(trial, budget)) {
        state = trial;
        accepted = true;
        break;
      }
      count(RecoilCounter::budget_rejection);
    }
    if (!accepted) {
      // The recoil's identity says this neutron left. Dropping it and banking
      // the recoil anyway -- which is what a `continue` here did -- makes the
      // record's mass number one higher than the momentum it carries.
      count(RecoilCounter::budget_exhaustion);
      count(RecoilCounter::incomplete_emission);
      return;
    }
  }

  if (n_extra == 0) {
    state.internal = remaining_internal_energy(
      budget, state.emitted_kin, state.momentum, state.mass);
  }

  finish_event(
    p, nuc, rx, wgt, E_in, u_in, ions, state, recoil, budget, n_extra > 0);
}

void from_absorption(Particle& p, int i_nuclide, double weight, double E_in,
  Direction u_in, const Reaction* rx)
{
  if (!settings::recoil_production || weight <= 0.0)
    return;

  if (!rx)
    return;
  const auto& nuc {data::nuclides[i_nuclide]};

  EmittedParticles emitted;
  if (!emitted_particles(rx->mt_, emitted)) {
    return; // unknown exit channel; see the note in from_inelastic()
  }
  bool budget_ok = false;
  double q = event_q(*nuc, *rx, emitted, budget_ok);
  if (!budget_ok)
    return; // see the note in from_inelastic()

  ParticleType recoil = recoil_particle_type(*nuc, rx->mt_);

  EmissionState state;
  state.momentum = neutron_momentum(E_in, u_in);
  state.mass = MASS_NEUTRON_EV + nuc->awr_ * MASS_NEUTRON_EV;
  state.za = {nuc->Z_, nuc->A_ + 1};

  double budget = E_in + q;

  // Radiative capture: the recoil is kicked by the emitted photons. Their
  // momenta are resampled from this reaction's own photon distribution so that
  // the recoil belongs to the reaction the ReactionFilter reports.
  if (rx->mt_ == N_GAMMA) {
    PhotonKick kick = sample_photon_kick(*rx, E_in, u_in, p.current_seed());
    constrain_photon_kick(kick, state.momentum, state.mass, budget);
    state.momentum -= kick.momentum;
    state.emitted_kin += kick.energy;
  }

  ChargedProducts ions;
  ions.fill(emitted);

  state.internal = remaining_internal_energy(
    budget, state.emitted_kin, state.momentum, state.mass);

  finish_event(
    p, *nuc, *rx, weight, E_in, u_in, ions, state, recoil, budget, false);
}

} // namespace recoil
} // namespace openmc
