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

#include <algorithm> // for max, min, swap
#include <cctype>    // for tolower, isdigit
#include <cmath>     // for sqrt, exp, cbrt, log
#include <string>

namespace openmc {
namespace recoil {

namespace {

//==============================================================================
// Tunable constants of the light-ion emission model
//
// These two numbers set the effective Coulomb barrier felt by an emitted light
// ion. They are not optical-model quantities: the reduced radius is larger and
// the diffuseness softer than a geometric barrier so that the single smooth
// transmission factor also stands in for sub-barrier tunnelling. Both were
// calibrated by matching the mean centre-of-mass light-ion energy of this model
// against the evaluated ENDF MF=6 spectra of MT=103-107 for 13 nuclides from
// Be-9 to Ta-181 between 5 and 20 MeV.
//==============================================================================

//! Effective barrier reduced radius in [fm]
constexpr double BARRIER_RADIUS = 1.7537;

//! Effective barrier diffuseness in [eV]
constexpr double BARRIER_DIFFUSENESS = 0.8e6;

//! Coulomb constant e^2/(4 pi eps0) in [eV fm]
constexpr double COULOMB_EV_FM = 1.44e6;

//! Maximum rejection attempts before falling back to a tabulated CDF
constexpr int MAX_REJECTION = 200;

//! Maximum attempts to resample a modelled product inside the energy budget
constexpr int MAX_BUDGET_TRIES = 20;

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

//! Q value of the channel that emits a single ion \p b from \p nuc
//!
//! Uses the evaluated Q value of the corresponding pure channel (MT=103-107)
//! when the library provides it. This is the energy that would be available to
//! the ion if it were the first particle emitted by the compound nucleus, and
//! it sets the shape of the modelled spectrum even in channels where other
//! particles have already taken part of the budget.
double single_channel_q(const Nuclide& nuc, AtomicNumbers b, double fallback)
{
  int mt = 0;
  if (b.Z == 1 && b.A == 1)
    mt = N_P;
  else if (b.Z == 1 && b.A == 2)
    mt = N_D;
  else if (b.Z == 1 && b.A == 3)
    mt = N_T;
  else if (b.Z == 2 && b.A == 3)
    mt = N_3HE;
  else if (b.Z == 2 && b.A == 4)
    mt = N_A;

  if (mt > 0) {
    // reaction_index_ is filled with C_NONE, which wraps to a huge value in
    // its unsigned element type; the bounds check below covers both cases.
    size_t i = nuc.reaction_index_[mt];
    if (i < nuc.reactions_.size()) {
      const auto& rx = nuc.reactions_[i];
      if (rx && rx->mt_ == mt)
        return rx->q_value_;
    }
  }
  return fallback;
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

//! True for the discrete charged-particle level MTs, which are exactly two-body
bool is_discrete_charged_level(int mt)
{
  return (mt >= N_P0 && mt < N_PC) || (mt >= N_D0 && mt < N_DC) ||
         (mt >= N_T0 && mt < N_TC) || (mt >= N_3HE0 && mt < N_3HEC) ||
         (mt >= N_A0 && mt < N_AC);
}

//! Internal energy of a system with the given lab momentum, mass, and budget
double internal_energy(double budget, Direction momentum, double mass)
{
  if (mass <= 0.0)
    return 0.0;
  return budget - momentum.dot(momentum) / (2.0 * mass);
}

} // namespace

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

double light_ion_pdf(double E, double E_max, int Z_b, int A_b, int Z_d, int A_d,
  const LightIonParams& par)
{
  if (E <= 0.0 || E >= E_max)
    return 0.0;

  double transmission = 1.0;
  if (Z_b > 0 && Z_d > 0) {
    double radius = par.r0 * (std::cbrt(static_cast<double>(A_b)) +
                               std::cbrt(static_cast<double>(A_d)));
    double barrier = COULOMB_EV_FM * Z_b * Z_d / radius;
    double arg = (barrier - E) / par.delta;
    if (arg > 60.0) {
      transmission = std::exp(-arg);
    } else if (arg > -60.0) {
      transmission = 1.0 / (1.0 + std::exp(arg));
    }
  }
  double endpoint = 1.0 - E / E_max;
  // nu = 1/2 is by far the common case and std::sqrt is both faster and more
  // accurate than std::pow for it
  double level_density =
    (par.nu == 0.5) ? std::sqrt(endpoint) : std::pow(endpoint, par.nu);
  return E * transmission * level_density;
}

double light_ion_pdf(double E, double E_max, int Z_b, int A_b, int Z_d, int A_d)
{
  return light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d, LightIonParams {});
}

double sample_light_ion_energy(double E_max, double E_limit, int Z_b, int A_b,
  int Z_d, int A_d, uint64_t* seed, const LightIonParams& par)
{
  if (E_max <= 0.0 || !std::isfinite(E_max))
    return 0.0;
  if (!(E_limit > 0.0) || E_limit > E_max)
    E_limit = E_max;

  // The spectrum rises from the Coulomb barrier and falls to zero at E_max, so
  // rejection against its maximum over [0, E_limit] converges quickly even when
  // the budget cuts deep into the sub-barrier tail.
  constexpr int N_SCAN = 33;
  double peak = 0.0;
  double best = 0.5 * E_limit;
  for (int i = 0; i < N_SCAN; ++i) {
    double E = E_limit * i / (N_SCAN - 1);
    double pdf = light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d, par);
    if (pdf > peak) {
      peak = pdf;
      best = E;
    }
  }
  if (peak <= 0.0)
    return 0.5 * E_limit;
  peak *= 1.1; // guard against the true maximum falling between scan points

  for (int i = 0; i < MAX_REJECTION; ++i) {
    double E = E_limit * prn(seed);
    if (prn(seed) * peak <= light_ion_pdf(E, E_max, Z_b, A_b, Z_d, A_d, par))
      return E;
  }
  // Spectrum too peaked to sample by rejection: use the scan maximum
  return best;
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
    // emits this ion alone. The latter sets the spectrum shape.
    double E_max_event = state.internal * m_d / (m_b + m_d);
    if (E_max_event <= 0.0 || !std::isfinite(E_max_event))
      return false;
    double q_single = single_channel_q(nuc, b, rx.q_value_);
    double E_max_shape =
      std::max(E_max_event, (E_in + q_single) * m_d / (m_b + m_d));

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

    AngularParams ang {};
    double slope = ang.slope_scale * kalbach_slope(E_in, E_cm, b, nuc);
    double r = kalbach_precompound_fraction(E_cm, E_max_shape, E_in, d, ang);
    double mu = sample_kalbach_mu(slope, r, seed);
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
    // daughter out of the internal energy budget.
    state.internal -= E_cm * (1.0 + m_b / m_d);
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
void finish_event(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double weight, double E_in, Direction u_in, const ChargedProducts& ions,
  EmissionState& state, ParticleType recoil)
{
  SampledIon sampled[ChargedProducts::MAX];
  int n_sampled = 0;

  if (ions.n > 0 &&
      settings::recoil.light_ion_model == RecoilLightIonModel::statistical) {
    EmissionState trial = state;
    if (emit_light_ions(
          p, nuc, rx, E_in, u_in, ions, trial, sampled, n_sampled)) {
      state = trial;
    } else {
      // Kinematically impossible exit channel for this event; keep the
      // momentum balance built from the evaluated products only.
      n_sampled = 0;
    }
  }

  if (settings::recoil.emitted_ions) {
    for (int i = 0; i < n_sampled; ++i) {
      p.create_secondary(
        weight, sampled[i].direction, sampled[i].energy, sampled[i].type);
    }
  }
  bank_recoil(p, nuc, weight, state.momentum, recoil);
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
  bank_recoil(p, nuc, p.wgt(), p_recoil, nuc.particle_type());
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
    // labelled with the wrong nuclide.
    return;
  }

  ParticleType recoil = recoil_particle_type(nuc, rx.mt_);

  ChargedProducts ions;
  ions.fill(emitted);

  // Start from the compound system and remove the transported neutron
  EmissionState state;
  state.momentum = neutron_momentum(E_in, u_in);
  state.mass = MASS_NEUTRON_EV + nuc.awr_ * MASS_NEUTRON_EV;
  state.za = {nuc.Z_, nuc.A_ + 1};
  state.emitted_kin = 0.0;

  double budget = E_in + rx.q_value_;

  state.momentum -= neutron_momentum(E_out, u_out);
  state.mass -= MASS_NEUTRON_EV;
  state.za.A -= 1;
  state.emitted_kin += E_out;

  // Additional neutrons of a multiplicity > 1 channel. OpenMC transports only
  // one sampled neutron, so the others are sampled independently from the same
  // evaluated distribution. The ENDF distribution is inclusive and carries no
  // joint final state, so a sample is rejected when it would overrun the
  // event's energy budget; that keeps every event kinematically possible while
  // leaving the marginal spectrum close to the evaluated one.
  int n_extra = 0;
  if (emitted.neutron > 1) {
    n_extra = emitted.neutron - 1;
  } else if (yield > 1.0 && std::floor(yield) == yield) {
    n_extra = static_cast<int>(std::round(yield)) - 1;
  }

  for (int i = 0; i < n_extra; ++i) {
    double E_extra;
    Direction u_extra;
    bool accepted = false;
    for (int attempt = 0; attempt < MAX_BUDGET_TRIES; ++attempt) {
      sample_reaction_neutron(
        nuc, rx, E_in, u_in, p.current_seed(), E_extra, u_extra);
      if (state.emitted_kin + E_extra <= budget) {
        accepted = true;
        break;
      }
    }
    if (!accepted)
      continue;
    state.momentum -= neutron_momentum(E_extra, u_extra);
    state.mass -= MASS_NEUTRON_EV;
    state.za.A -= 1;
    state.emitted_kin += E_extra;
  }

  state.internal =
    internal_energy(budget - state.emitted_kin, state.momentum, state.mass);

  finish_event(p, nuc, rx, wgt, E_in, u_in, ions, state, recoil);
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
  ParticleType recoil = recoil_particle_type(*nuc, rx->mt_);

  EmissionState state;
  state.momentum = neutron_momentum(E_in, u_in);
  state.mass = MASS_NEUTRON_EV + nuc->awr_ * MASS_NEUTRON_EV;
  state.za = {nuc->Z_, nuc->A_ + 1};

  double budget = E_in + rx->q_value_;

  // Radiative capture: the recoil is kicked by the emitted photons. Their
  // momenta are resampled from this reaction's own photon distribution so that
  // the recoil belongs to the reaction the ReactionFilter reports.
  if (rx->mt_ == N_GAMMA) {
    PhotonKick kick = sample_photon_kick(*rx, E_in, u_in, p.current_seed());
    state.momentum -= kick.momentum;
    state.emitted_kin += kick.energy;
  }

  ChargedProducts ions;
  ions.fill(emitted);

  state.internal =
    internal_energy(budget - state.emitted_kin, state.momentum, state.mass);

  finish_event(p, *nuc, *rx, weight, E_in, u_in, ions, state, recoil);
}

} // namespace recoil
} // namespace openmc
