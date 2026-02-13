#include "openmc/physics.h"

#include "openmc/bank.h"
#include "openmc/bremsstrahlung.h"
#include "openmc/chain.h"
#include "openmc/constants.h"
#include "openmc/distribution_multi.h"
#include "openmc/eigenvalue.h"
#include "openmc/endf.h"
#include "openmc/error.h"
#include "openmc/ifp.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/ncrystal_interface.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/physics_common.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/reaction.h"
#include "openmc/search.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"
#include "openmc/tallies/tally.h"
#include "openmc/thermal.h"
#include "openmc/weight_windows.h"

#include <fmt/core.h>

#include "openmc/tensor.h"
#include <algorithm> // for max, min, max_element
#include <cctype>    // for tolower, isdigit
#include <cmath>     // for sqrt, exp, log, abs, copysign

namespace openmc {

namespace {

double particle_mass_ev(ParticleType type)
{
  switch (type.pdg_number()) {
  case PDG_NEUTRON:
  case PDG_PROTON:
    return MASS_NEUTRON_EV;
  case PDG_DEUTERON:
    return 2.0 * MASS_NEUTRON_EV;
  case PDG_TRITON:
    return 3.0 * MASS_NEUTRON_EV;
  case PDG_ALPHA:
    return 4.0 * MASS_NEUTRON_EV;
  case PDG_PHOTON:
    return 0.0;
  default:
    if (type.is_nucleus()) {
      int A = (type.pdg_number() / 10) % 1000;
      if (A > 0) {
        return A * MASS_NEUTRON_EV;
      }
    }
    return MASS_NEUTRON_EV;
  }
}

Direction momentum_from_kinetic_energy(double mass, double E, Direction u)
{
  if (mass <= 0.0 || E <= 0.0) {
    return {};
  }
  return std::sqrt(2.0 * mass * E) * u;
}

Direction neutron_momentum(double E, Direction u)
{
  return momentum_from_kinetic_energy(MASS_NEUTRON_EV, E, u);
}

Direction photon_momentum(double E, Direction u)
{
  if (E <= 0.0) {
    return {};
  }
  return E * u;
}

struct PhotonMomentumInfo {
  Direction momentum {};
  double energy {0.0};
};

bool parse_emitted_particles_from_channel(std::string channel, int& n_neutron,
  int& n_proton, int& n_deuteron, int& n_triton, int& n_he3, int& n_alpha)
{
  n_neutron = 0;
  n_proton = 0;
  n_deuteron = 0;
  n_triton = 0;
  n_he3 = 0;
  n_alpha = 0;

  for (auto& c : channel) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (channel == "gamma") {
    return true;
  }

  int i = 0;
  while (i < channel.size()) {
    int multiplicity = 0;
    while (i < channel.size() &&
           std::isdigit(static_cast<unsigned char>(channel[i]))) {
      multiplicity = 10 * multiplicity + (channel[i] - '0');
      ++i;
    }
    if (multiplicity == 0) {
      multiplicity = 1;
    }
    if (i >= channel.size()) {
      return false;
    }

    if (i + 2 < channel.size() && channel.substr(i, 3) == "3he") {
      n_he3 += multiplicity;
      i += 3;
      continue;
    }
    if (i + 4 < channel.size() && channel.substr(i, 5) == "gamma") {
      i += 5;
      continue;
    }

    switch (channel[i]) {
    case 'n':
      n_neutron += multiplicity;
      ++i;
      break;
    case 'p':
      n_proton += multiplicity;
      ++i;
      break;
    case 'd':
      n_deuteron += multiplicity;
      ++i;
      break;
    case 't':
      n_triton += multiplicity;
      ++i;
      break;
    case 'a':
      n_alpha += multiplicity;
      ++i;
      break;
    default:
      return false;
    }
  }
  return true;
}

bool emitted_particle_counts(int mt, int& n_neutron, int& n_proton,
  int& n_deuteron, int& n_triton, int& n_he3, int& n_alpha)
{
  n_neutron = 0;
  n_proton = 0;
  n_deuteron = 0;
  n_triton = 0;
  n_he3 = 0;
  n_alpha = 0;

  if (mt == ELASTIC || (mt >= N_N1 && mt <= N_NC)) {
    n_neutron = 1;
    return true;
  }
  if (mt >= N_P0 && mt <= N_PC) {
    n_proton = 1;
    return true;
  }
  if (mt >= N_D0 && mt <= N_DC) {
    n_deuteron = 1;
    return true;
  }
  if (mt >= N_T0 && mt <= N_TC) {
    n_triton = 1;
    return true;
  }
  if (mt >= N_3HE0 && mt <= N_3HEC) {
    n_he3 = 1;
    return true;
  }
  if (mt >= N_A0 && mt <= N_AC) {
    n_alpha = 1;
    return true;
  }
  if (mt >= N_2N0 && mt <= N_2NC) {
    n_neutron = 2;
    return true;
  }

  std::string name = reaction_name(mt);
  if (name.size() < 5 || name[0] != '(' || name[1] != 'n' || name[2] != ',') {
    return false;
  }
  if (name.back() != ')') {
    return false;
  }

  std::string channel = name.substr(3, name.size() - 4);
  return parse_emitted_particles_from_channel(
    channel, n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha);
}

ParticleType residual_particle_type(const Nuclide& nuc, int mt)
{
  int n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha;
  if (!emitted_particle_counts(
        mt, n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha)) {
    return nuc.particle_type();
  }

  int emitted_A = n_neutron + n_proton + 2 * n_deuteron + 3 * n_triton +
                  3 * n_he3 + 4 * n_alpha;
  int emitted_Z = n_proton + n_deuteron + n_triton + 2 * n_he3 + 2 * n_alpha;

  int Z_res = nuc.Z_ - emitted_Z;
  int A_res = nuc.A_ + 1 - emitted_A;
  if (Z_res <= 0 || A_res <= 0) {
    return nuc.particle_type();
  }
  return ParticleType {Z_res, A_res, 0};
}

bool create_recoil_secondary(Particle& p, const Nuclide& nuc, double weight,
  Direction p_recoil, ParticleType type)
{
  if (!settings::recoil_production || !settings::recoil.bank_residual ||
      weight <= 0.0) {
    return false;
  }

  double p2 = p_recoil.dot(p_recoil);
  if (!std::isfinite(p2) || p2 <= 0.0) {
    return false;
  }

  double mass = particle_mass_ev(type);
  if (mass <= 0.0 || !std::isfinite(mass)) {
    mass = nuc.awr_ * MASS_NEUTRON_EV;
  }
  if (mass <= 0.0 || !std::isfinite(mass)) {
    return false;
  }

  double E_recoil = p2 / (2.0 * mass);
  if (!std::isfinite(E_recoil) || E_recoil <= 0.0) {
    return false;
  }

  Direction u_recoil;
  if (settings::recoil.direction == RecoilDirection::momentum) {
    u_recoil = p_recoil / std::sqrt(p2);
  } else {
    u_recoil = isotropic_direction(p.current_seed());
  }

  return p.create_secondary(weight, u_recoil, E_recoil, type);
}

void check_recoil_sanity(
  const Nuclide& nuc, const Reaction& rx, double E_in, double emitted_kinetic)
{
  if (!settings::recoil.q_sanity_check) {
    return;
  }

  double available = E_in + rx.q_value_;
  if (!std::isfinite(available) || available <= 0.0) {
    return;
  }

  if (emitted_kinetic > available * (1.0 + 1e-6)) {
    static int n_warnings = 0;
    if (n_warnings < 10) {
      warning(fmt::format(
        "Recoil Q check: emitted kinetic energy exceeds E + Q for {} MT={} "
        "(E={} eV, Q={} eV, emitted={} eV).",
        nuc.name_, rx.mt_, E_in, rx.q_value_, emitted_kinetic));
      ++n_warnings;
    }
    if (settings::recoil.fail_on_nonphysical) {
      fatal_error(
        "Recoil nonphysical kinematics encountered with fail_on_nonphysical.");
    }
  }
}

const Reaction* sample_disappearance_reaction(int i_nuclide, Particle& p)
{
  const auto& nuc {data::nuclides[i_nuclide]};
  const auto& micro {p.neutron_xs(i_nuclide)};

  double total = 0.0;
  for (const auto& rx : nuc->reactions_) {
    if (rx->redundant_) {
      continue;
    }
    if (!is_disappearance(rx->mt_) || is_fission(rx->mt_)) {
      continue;
    }
    total += rx->xs(micro);
  }
  if (total <= 0.0) {
    return nullptr;
  }

  double cutoff = prn(p.current_seed()) * total;
  double prob = 0.0;
  for (const auto& rx : nuc->reactions_) {
    if (rx->redundant_) {
      continue;
    }
    if (!is_disappearance(rx->mt_) || is_fission(rx->mt_)) {
      continue;
    }
    prob += rx->xs(micro);
    if (prob > cutoff) {
      return rx.get();
    }
  }

  return nullptr;
}

bool choose_charged_absorption_product(
  int mt, ParticleType& emitted_type, bool& has_multiple_products)
{
  int n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha;
  has_multiple_products = false;
  if (!emitted_particle_counts(
        mt, n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha)) {
    return false;
  }
  if (n_neutron != 0) {
    return false;
  }

  int n_charged = n_proton + n_deuteron + n_triton + n_he3 + n_alpha;
  if (n_charged == 0) {
    return false;
  }
  if (n_charged != 1) {
    has_multiple_products = true;
    return false;
  }

  if (n_proton == 1) {
    emitted_type = ParticleType::proton();
  } else if (n_deuteron == 1) {
    emitted_type = ParticleType::deuteron();
  } else if (n_triton == 1) {
    emitted_type = ParticleType::triton();
  } else if (n_he3 == 1) {
    emitted_type = ParticleType {2, 3, 0};
  } else if (n_alpha == 1) {
    emitted_type = ParticleType::alpha();
  } else {
    return false;
  }
  return true;
}

bool sample_two_body_charged_product(const Nuclide& nuc, const Reaction& rx,
  ParticleType emitted_type, double E_in, Direction u_in, uint64_t* seed,
  Direction& p_emitted, double& E_emitted, Direction& u_emitted)
{
  double m_n = MASS_NEUTRON_EV;
  double m_target = nuc.awr_ * MASS_NEUTRON_EV;
  double m_emitted = particle_mass_ev(emitted_type);
  double m_residual = particle_mass_ev(residual_particle_type(nuc, rx.mt_));

  if (m_target <= 0.0 || m_emitted <= 0.0 || m_residual <= 0.0) {
    return false;
  }

  double E_cm = E_in * m_target / (m_n + m_target);
  double available = E_cm + rx.q_value_;
  if (available <= 0.0) {
    return false;
  }

  double mu_final = m_emitted * m_residual / (m_emitted + m_residual);
  double p_cm = std::sqrt(2.0 * mu_final * available);
  Direction u_cm = isotropic_direction(seed);

  Direction p_in = neutron_momentum(E_in, u_in);
  Direction v_cm = p_in / (m_n + m_target);
  Direction v_emitted = v_cm + (p_cm / m_emitted) * u_cm;
  p_emitted = m_emitted * v_emitted;

  double p2 = p_emitted.dot(p_emitted);
  if (p2 <= 0.0) {
    return false;
  }
  E_emitted = p2 / (2.0 * m_emitted);
  if (!std::isfinite(E_emitted) || E_emitted <= 0.0) {
    return false;
  }

  u_emitted = p_emitted / std::sqrt(p2);
  return true;
}

bool recoil_includes_photon_momentum(int mt)
{
  switch (settings::recoil.include_photon_momentum) {
  case RecoilIncludePhotonMomentum::none:
    return false;
  case RecoilIncludePhotonMomentum::all:
    return true;
  case RecoilIncludePhotonMomentum::capture_only:
    return mt == N_GAMMA;
  }
  return false;
}

PhotonMomentumInfo sample_capture_photon_momentum(
  const Reaction& rx, double E_in, Direction u_in, uint64_t* seed)
{
  PhotonMomentumInfo info;
  for (const auto& product : rx.products_) {
    if (!product.particle_.is_photon()) {
      continue;
    }

    double y = (*product.yield_)(E_in);
    if (y <= 0.0) {
      continue;
    }

    int n = static_cast<int>(y);
    if (prn(seed) < y - n) {
      ++n;
    }

    for (int i = 0; i < n; ++i) {
      double E_gamma;
      double mu_gamma;
      product.sample(E_in, E_gamma, mu_gamma, seed);
      Direction u_gamma = rotate_angle(u_in, mu_gamma, nullptr, seed);
      info.momentum += photon_momentum(E_gamma, u_gamma);
      info.energy += E_gamma;
    }
  }
  return info;
}

PhotonMomentumInfo banked_capture_photon_momentum(
  const Particle& p, double wgt_in)
{
  PhotonMomentumInfo info;
  int start = p.secondary_bank_index();
  int end = start + p.n_secondaries();

  for (int i = start; i < end; ++i) {
    const auto& site = p.secondary_bank(i);
    if (!site.particle.is_photon()) {
      continue;
    }

    double multiplicity = 1.0;
    if (wgt_in > 0.0) {
      multiplicity = site.wgt / wgt_in;
      if (settings::run_mode == RunMode::EIGENVALUE && simulation::keff > 0.0) {
        multiplicity /= simulation::keff;
      }
    }
    multiplicity = std::max(0.0, multiplicity);
    info.momentum += multiplicity * photon_momentum(site.E, site.u);
    info.energy += multiplicity * site.E;
  }
  return info;
}

void create_absorption_recoil(Particle& p, int i_nuclide, double recoil_weight,
  double E_in, Direction u_in, double incident_wgt, const Reaction* rx)
{
  if (!settings::recoil_production || recoil_weight <= 0.0) {
    return;
  }

  const auto& nuc {data::nuclides[i_nuclide]};
  if (!rx) {
    rx = sample_disappearance_reaction(i_nuclide, p);
  }
  if (!rx) {
    return;
  }

  Direction p_in = neutron_momentum(E_in, u_in);
  Direction p_recoil = p_in;
  double emitted_kinetic = 0.0;
  ParticleType recoil_type = residual_particle_type(*nuc, rx->mt_);

  if (rx->mt_ == N_GAMMA) {
    if (recoil_includes_photon_momentum(rx->mt_)) {
      PhotonMomentumInfo gamma_info;
      if (settings::recoil.capture_photons == RecoilCapturePhotons::phantom) {
        gamma_info =
          sample_capture_photon_momentum(*rx, E_in, u_in, p.current_seed());
      } else {
        gamma_info = banked_capture_photon_momentum(p, incident_wgt);
      }
      p_recoil -= gamma_info.momentum;
      emitted_kinetic = gamma_info.energy;
    }
  } else {
    ParticleType emitted_type;
    bool has_multiple_products = false;
    if (choose_charged_absorption_product(
          rx->mt_, emitted_type, has_multiple_products)) {
      Direction p_emitted;
      Direction u_emitted;
      double E_emitted;
      if (sample_two_body_charged_product(*nuc, *rx, emitted_type, E_in, u_in,
            p.current_seed(), p_emitted, E_emitted, u_emitted)) {
        p_recoil -= p_emitted;
        emitted_kinetic = E_emitted;
        if (settings::recoil.bank_emitted_ions) {
          p.create_secondary(recoil_weight, u_emitted, E_emitted, emitted_type);
        }
      } else if (settings::recoil.fail_on_nonphysical) {
        fatal_error(
          "Failed to sample two-body charged absorption recoil kinematics.");
      }
    } else if (has_multiple_products) {
      static bool warned_multicharged_absorption {false};
      if (!warned_multicharged_absorption) {
        warning("Multi-particle charged absorption recoil fallback: using "
                "compound recoil from incident neutron momentum only.");
        warned_multicharged_absorption = true;
      }
    }
  }

  check_recoil_sanity(*nuc, *rx, E_in, emitted_kinetic);
  create_recoil_secondary(p, *nuc, recoil_weight, p_recoil, recoil_type);
}

} // namespace

//==============================================================================
// Non-member functions
//==============================================================================

void collision(Particle& p)
{
  // Add to collision counter for particle
  ++(p.n_collision());
  p.secondary_bank_index() = p.local_secondary_bank().size();

  // Sample reaction for the material the particle is in
  switch (p.type().pdg_number()) {
  case PDG_NEUTRON:
    sample_neutron_reaction(p);
    break;
  case PDG_PHOTON:
    sample_photon_reaction(p);
    break;
  case PDG_ELECTRON:
    sample_electron_reaction(p);
    break;
  case PDG_POSITRON:
    sample_positron_reaction(p);
    break;
  default:
    sample_other_reaction(p);
    break;
  }

  if (settings::weight_windows_on) {
    auto [ww_found, ww] = search_weight_window(p);
    if (!ww_found && p.type() == ParticleType::neutron()) {
      // if the weight window is not valid, apply russian roulette for neutrons
      // (regardless of weight window collision checkpoint setting)
      apply_russian_roulette(p);
    } else if (settings::weight_window_checkpoint_collision) {
      // if collision checkpointing is on, apply weight window
      apply_weight_window(p, ww);
    }
  }

  // Kill particle if energy falls below cutoff
  int type = p.type().transport_index();
  if (type == C_NONE || p.E() < settings::energy_cutoff[type]) {
    p.wgt() = 0.0;
  }

  // Display information about collision
  if (settings::verbosity >= 10 || p.trace()) {
    std::string msg;
    if (p.event() == TallyEvent::KILL) {
      msg = fmt::format("    Killed. Energy = {} eV.", p.E());
    } else if (p.type().is_neutron()) {
      msg = fmt::format("    {} with {}. Energy = {} eV.",
        reaction_name(p.event_mt()), data::nuclides[p.event_nuclide()]->name_,
        p.E());
    } else if (p.type().is_photon()) {
      msg = fmt::format("    {} with {}. Energy = {} eV.",
        reaction_name(p.event_mt()),
        to_element(data::nuclides[p.event_nuclide()]->name_), p.E());
    } else {
      msg = fmt::format("    Disappeared. Energy = {} eV.", p.E());
    }
    write_message(msg, 1);
  }
}

void sample_neutron_reaction(Particle& p)
{
  // Sample a nuclide within the material
  int i_nuclide = sample_nuclide(p);

  // Save which nuclide particle had collision with
  p.event_nuclide() = i_nuclide;

  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  const auto& nuc {data::nuclides[i_nuclide]};

  if (nuc->fissionable_ && p.neutron_xs(i_nuclide).fission > 0.0) {
    auto& rx = sample_fission(i_nuclide, p);
    if (settings::run_mode == RunMode::EIGENVALUE) {
      create_fission_sites(p, i_nuclide, rx);
    } else if (settings::run_mode == RunMode::FIXED_SOURCE &&
               settings::create_fission_neutrons) {
      create_fission_sites(p, i_nuclide, rx);

      // Make sure particle population doesn't grow out of control for
      // subcritical multiplication problems.
      if (p.local_secondary_bank().size() >= settings::max_secondaries &&
          !settings::use_shared_secondary_bank) {
        fatal_error(
          "The secondary particle bank appears to be growing without "
          "bound. You are likely running a subcritical multiplication problem "
          "with k-effective close to or greater than one.");
      }
    }
    p.event_mt() = rx.mt_;
  }

  // Create secondary photons
  if (settings::photon_transport) {
    sample_secondary_photons(p, i_nuclide);
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs

  if (p.neutron_xs(i_nuclide).absorption > 0.0) {
    absorption(p, i_nuclide);
  }
  if (!p.alive())
    return;

  // Sample a scattering reaction and determine the secondary energy of the
  // exiting neutron
  const auto& ncrystal_mat = model::materials[p.material()]->ncrystal_mat();
  if (ncrystal_mat && p.E() < NCRYSTAL_MAX_ENERGY) {
    ncrystal_mat.scatter(p);
  } else {
    scatter(p, i_nuclide);
  }

  // Advance URR seed stream 'N' times after energy changes
  if (p.E() != p.E_last()) {
    advance_prn_seed(data::nuclides.size(), &p.seeds(STREAM_URR_PTABLE));
  }

  // Play russian roulette if there are no weight windows
  if (!settings::weight_windows_on)
    apply_russian_roulette(p);
}

void create_fission_sites(Particle& p, int i_nuclide, const Reaction& rx)
{
  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p.wgt() / simulation::keff * weight *
                p.neutron_xs(i_nuclide).nu_fission /
                p.neutron_xs(i_nuclide).total;

  // Sample the number of neutrons produced
  int nu = static_cast<int>(nu_t);
  if (prn(p.current_seed()) <= (nu_t - nu))
    ++nu;

  // If no neutrons were produced then don't continue
  if (nu == 0)
    return;

  // Initialize the counter of delayed neutrons encountered for each delayed
  // group.
  double nu_d[MAX_DELAYED_GROUPS] = {0.};

  // Clear out particle's nu fission bank
  p.nu_bank().clear();

  p.fission() = true;

  // Determine whether to place fission sites into the shared fission bank
  // or the secondary particle bank.
  bool use_fission_bank = (settings::run_mode == RunMode::EIGENVALUE);

  // Counter for the number of fission sites successfully stored to the shared
  // fission bank or the secondary particle bank
  int n_sites_stored;

  for (n_sites_stored = 0; n_sites_stored < nu; n_sites_stored++) {
    // Initialize fission site object with particle data
    SourceSite site;
    site.r = p.r();
    site.particle = ParticleType::neutron();
    site.time = p.time();
    site.wgt = 1. / weight;
    site.surf_id = 0;

    // Sample delayed group and angle/energy for fission reaction
    sample_fission_neutron(i_nuclide, rx, &site, p);

    // Reject site if it exceeds time cutoff
    if (site.delayed_group > 0) {
      double t_cutoff = settings::time_cutoff[site.particle.transport_index()];
      if (site.time > t_cutoff) {
        continue;
      }
    }

    // Set parent and progeny IDs
    site.parent_id = p.current_work();
    site.progeny_id = p.n_progeny()++;

    // Store fission site in bank
    if (use_fission_bank) {
      int64_t idx = simulation::fission_bank.thread_safe_append(site);
      if (idx == -1) {
        warning(
          "The shared fission bank is full. Additional fission sites created "
          "in this generation will not be banked. Results may be "
          "non-deterministic.");

        // Decrement number of particle progeny as storage was unsuccessful.
        // This step is needed so that the sum of all progeny is equal to the
        // size of the shared fission bank.
        p.n_progeny()--;

        // Break out of loop as no more sites can be added to fission bank
        break;
      }
      // Iterated Fission Probability (IFP) method
      if (settings::ifp_on) {
        ifp(p, idx);
      }
    } else {
      site.wgt_born = p.wgt_born();
      site.wgt_ww_born = p.wgt_ww_born();
      site.n_split = p.n_split();
      p.local_secondary_bank().push_back(site);
      p.n_secondaries()++;
    }

    // Increment the number of neutrons born delayed
    if (site.delayed_group > 0) {
      nu_d[site.delayed_group - 1]++;
    }

    // Write fission particles to nuBank
    NuBank& nu_bank_entry = p.nu_bank().emplace_back();
    nu_bank_entry.wgt = site.wgt;
    nu_bank_entry.E = site.E;
    nu_bank_entry.delayed_group = site.delayed_group;
  }

  // If shared fission bank was full, and no fissions could be added,
  // set the particle fission flag to false.
  if (n_sites_stored == 0) {
    p.fission() = false;
    return;
  }

  // Set nu to the number of fission sites successfully stored. If the fission
  // bank was not found to be full then these values are already equivalent.
  nu = n_sites_stored;

  // Store the total weight banked for analog fission tallies
  p.n_bank() = nu;
  p.wgt_bank() = nu / weight;
  for (size_t d = 0; d < MAX_DELAYED_GROUPS; d++) {
    p.n_delayed_bank(d) = nu_d[d];
  }
}

void sample_photon_reaction(Particle& p)
{
  // Kill photon if below energy cutoff -- an extra check is made here because
  // photons with energy below the cutoff may have been produced by neutrons
  // reactions or atomic relaxation
  int photon = ParticleType::photon().transport_index();
  if (p.E() < settings::energy_cutoff[photon]) {
    p.E() = 0.0;
    p.wgt() = 0.0;
    return;
  }

  // Sample element within material
  int i_element = sample_element(p);
  const auto& micro {p.photon_xs(i_element)};
  const auto& element {*data::elements[i_element]};

  // Calculate photon energy over electron rest mass equivalent
  double alpha = p.E() / MASS_ELECTRON_EV;

  // For tallying purposes, this routine might be called directly. In that
  // case, we need to sample a reaction via the cutoff variable
  double prob = 0.0;
  double cutoff = prn(p.current_seed()) * micro.total;

  // Coherent (Rayleigh) scattering
  prob += micro.coherent;
  if (prob > cutoff) {
    p.mu() = element.rayleigh_scatter(alpha, p.current_seed());
    p.u() = rotate_angle(p.u(), p.mu(), nullptr, p.current_seed());
    p.event() = TallyEvent::SCATTER;
    p.event_mt() = COHERENT;
    return;
  }

  // Incoherent (Compton) scattering
  prob += micro.incoherent;
  if (prob > cutoff) {
    double alpha_out;
    int i_shell;
    element.compton_scatter(
      alpha, true, &alpha_out, &p.mu(), &i_shell, p.current_seed());

    // Determine binding energy of shell. The binding energy is 0.0 if
    // doppler broadening is not used.
    double e_b;
    if (i_shell == -1) {
      e_b = 0.0;
    } else {
      e_b = element.binding_energy_[i_shell];
    }

    // Create Compton electron
    double phi = uniform_distribution(0., 2.0 * PI, p.current_seed());
    double E_electron = (alpha - alpha_out) * MASS_ELECTRON_EV - e_b;
    int electron = ParticleType::electron().transport_index();
    if (E_electron >= settings::energy_cutoff[electron]) {
      double mu_electron = (alpha - alpha_out * p.mu()) /
                           std::sqrt(alpha * alpha + alpha_out * alpha_out -
                                     2.0 * alpha * alpha_out * p.mu());
      Direction u = rotate_angle(p.u(), mu_electron, &phi, p.current_seed());
      process_charged_secondary(p, u, E_electron, ParticleType::electron());
    }

    // Allow electrons to fill orbital and produce Auger electrons and
    // fluorescent photons. Since Compton subshell data does not match atomic
    // relaxation data, use the mapping between the data to find the subshell
    if (settings::atomic_relaxation && element.has_atomic_relaxation_ &&
        i_shell >= 0 && element.subshell_map_[i_shell] >= 0) {
      element.atomic_relaxation(element.subshell_map_[i_shell], p);
    }

    phi += PI;
    p.E() = alpha_out * MASS_ELECTRON_EV;
    p.u() = rotate_angle(p.u(), p.mu(), &phi, p.current_seed());
    p.event() = TallyEvent::SCATTER;
    p.event_mt() = INCOHERENT;
    return;
  }

  // Photoelectric effect
  double prob_after = prob + micro.photoelectric;

  if (prob_after > cutoff) {
    // Get grid index, interpolation factor, and bounding subshell
    // cross sections
    int i_grid = micro.index_grid;
    double f = micro.interp_factor;
    tensor::View<const double> xs_lower = element.cross_sections_.slice(i_grid);
    tensor::View<const double> xs_upper =
      element.cross_sections_.slice(i_grid + 1);

    for (int i_shell = 0; i_shell < element.shells_.size(); ++i_shell) {
      const auto& shell {element.shells_[i_shell]};

      // Check threshold of reaction
      if (xs_lower(i_shell) == 0)
        continue;

      //  Evaluation subshell photoionization cross section
      prob += std::exp(
        xs_lower(i_shell) + f * (xs_upper(i_shell) - xs_lower(i_shell)));

      if (prob > cutoff) {
        // Determine binding energy based on whether atomic relaxation data is
        // present (if not, use value from Compton profile data)
        double binding_energy = element.has_atomic_relaxation_
                                  ? shell.binding_energy
                                  : element.binding_energy_[i_shell];

        // Determine energy of secondary electron
        double E_electron = p.E() - binding_energy;

        // Sample mu using non-relativistic Sauter distribution.
        // See Eqns 3.19 and 3.20 in "Implementing a photon physics
        // model in Serpent 2" by Toni Kaltiaisenaho
        double mu;
        while (true) {
          double r = prn(p.current_seed());
          if (4.0 * (1.0 - r) * r >= prn(p.current_seed())) {
            double rel_vel =
              std::sqrt(E_electron * (E_electron + 2.0 * MASS_ELECTRON_EV)) /
              (E_electron + MASS_ELECTRON_EV);
            mu =
              (2.0 * r + rel_vel - 1.0) / (2.0 * rel_vel * r - rel_vel + 1.0);
            break;
          }
        }

        double phi = uniform_distribution(0., 2.0 * PI, p.current_seed());
        Direction u;
        u.x = mu;
        u.y = std::sqrt(1.0 - mu * mu) * std::cos(phi);
        u.z = std::sqrt(1.0 - mu * mu) * std::sin(phi);

        // Process secondary electron at the photon collision site.
        process_charged_secondary(p, u, E_electron, ParticleType::electron());

        // Allow electrons to fill orbital and produce auger electrons
        // and fluorescent photons
        if (settings::atomic_relaxation) {
          element.atomic_relaxation(i_shell, p);
        }
        p.event() = TallyEvent::ABSORB;
        p.event_mt() = 533 + shell.index_subshell;
        p.wgt() = 0.0;
        p.E() = 0.0;
        return;
      }
    }
  }
  prob = prob_after;

  // Pair production
  prob += micro.pair_production;
  if (prob > cutoff) {
    double E_electron, E_positron;
    double mu_electron, mu_positron;
    element.pair_production(alpha, &E_electron, &E_positron, &mu_electron,
      &mu_positron, p.current_seed());

    // Process secondary electron at the photon collision site.
    Direction u = rotate_angle(p.u(), mu_electron, nullptr, p.current_seed());
    process_charged_secondary(p, u, E_electron, ParticleType::electron());

    // Process secondary positron at the photon collision site.
    u = rotate_angle(p.u(), mu_positron, nullptr, p.current_seed());
    process_charged_secondary(p, u, E_positron, ParticleType::positron());
    p.event() = TallyEvent::ABSORB;
    p.event_mt() = PAIR_PROD;
    p.wgt() = 0.0;
    p.E() = 0.0;
  }
}

void process_charged_secondary(
  Particle& p, Direction u, double E, ParticleType type)
{
  int idx = type.transport_index();
  if (idx == C_NONE || E < settings::energy_cutoff[idx])
    return;

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    thick_target_bremsstrahlung(p, type, u, E);
  }

  if (type == ParticleType::positron()) {
    Direction photon_u = isotropic_direction(p.current_seed());
    p.create_secondary(
      p.wgt(), photon_u, MASS_ELECTRON_EV, ParticleType::photon());
    p.create_secondary(
      p.wgt(), -photon_u, MASS_ELECTRON_EV, ParticleType::photon());

    // The annihilation photons are now emitted during the parent photon
    // collision. Offset the pair-production Q value in the energy balance so
    // heating matches the prior explicit positron slowing-down sequence.
    p.bank_second_E() -= 2 * MASS_ELECTRON_EV;
  }
}

void sample_electron_reaction(Particle& p)
{
  // TODO: create reaction types

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    thick_target_bremsstrahlung(p);
  }

  p.E() = 0.0;
  p.wgt() = 0.0;
  p.event() = TallyEvent::ABSORB;
}

void sample_positron_reaction(Particle& p)
{
  // TODO: create reaction types

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    thick_target_bremsstrahlung(p);
  }

  // Sample angle isotropically
  Direction u = isotropic_direction(p.current_seed());

  // Create annihilation photon pair traveling in opposite directions
  p.create_secondary(p.wgt(), u, MASS_ELECTRON_EV, ParticleType::photon());
  p.create_secondary(p.wgt(), -u, MASS_ELECTRON_EV, ParticleType::photon());

  p.E() = 0.0;
  p.wgt() = 0.0;
  p.event() = TallyEvent::ABSORB;
}

void sample_other_reaction(Particle& p)
{
  // For particles we don't handle, just kill the particle
  p.E() = 0.0;
  p.wgt() = 0.0;
  p.event() = TallyEvent::ABSORB;
}

int sample_nuclide(Particle& p)
{
  // Sample cumulative distribution function
  double cutoff = prn(p.current_seed()) * p.macro_xs().total;

  // Get pointers to nuclide/density arrays
  const auto& mat {model::materials[p.material()]};
  int n = mat->nuclide_.size();

  double prob = 0.0;
  for (int i = 0; i < n; ++i) {
    // Get atom density
    int i_nuclide = mat->nuclide_[i];
    double atom_density = mat->atom_density(i, p.density_mult());

    // Increment probability to compare to cutoff
    prob += atom_density * p.neutron_xs(i_nuclide).total;
    if (prob >= cutoff)
      return i_nuclide;
  }

  // If we reach here, no nuclide was sampled
  p.write_restart();
  throw std::runtime_error {"Did not sample any nuclide during collision."};
}

int sample_element(Particle& p)
{
  // Sample cumulative distribution function
  double cutoff = prn(p.current_seed()) * p.macro_xs().total;

  // Get pointers to elements, densities
  const auto& mat {model::materials[p.material()]};

  double prob = 0.0;
  for (int i = 0; i < mat->element_.size(); ++i) {
    // Find atom density
    int i_element = mat->element_[i];
    double atom_density = mat->atom_density(i, p.density_mult());

    // Determine microscopic cross section
    double sigma = atom_density * p.photon_xs(i_element).total;

    // Increment probability to compare to cutoff
    prob += sigma;
    if (prob > cutoff) {
      // Save which nuclide particle had collision with for tally purpose
      p.event_nuclide() = mat->nuclide_[i];

      return i_element;
    }
  }

  // If we made it here, no element was sampled
  p.write_restart();
  fatal_error("Did not sample any element during collision.");
}

Reaction& sample_fission(int i_nuclide, Particle& p)
{
  // Get pointer to nuclide
  const auto& nuc {data::nuclides[i_nuclide]};

  // If we're in the URR, by default use the first fission reaction. We also
  // default to the first reaction if we know that there are no partial fission
  // reactions
  if (p.neutron_xs(i_nuclide).use_ptable || !nuc->has_partial_fission_) {
    return *nuc->fission_rx_[0];
  }

  // Check to see if we are in a windowed multipole range.  WMP only supports
  // the first fission reaction.
  if (nuc->multipole_) {
    if (p.E() >= nuc->multipole_->E_min_ && p.E() <= nuc->multipole_->E_max_) {
      return *nuc->fission_rx_[0];
    }
  }

  // Get grid index and interpolation factor and sample fission cdf
  const auto& micro = p.neutron_xs(i_nuclide);
  double cutoff = prn(p.current_seed()) * p.neutron_xs(i_nuclide).fission;
  double prob = 0.0;

  // Loop through each partial fission reaction type
  for (auto& rx : nuc->fission_rx_) {
    // add to cumulative probability
    prob += rx->xs(micro);

    // Create fission bank sites if fission occurs
    if (prob > cutoff)
      return *rx;
  }

  // If we reached here, no reaction was sampled
  throw std::runtime_error {
    "No fission reaction was sampled for " + nuc->name_};
}

void sample_photon_product(
  int i_nuclide, Particle& p, int* i_rx, int* i_product)
{
  // Get grid index and interpolation factor and sample photon production cdf
  const auto& micro = p.neutron_xs(i_nuclide);
  double cutoff = prn(p.current_seed()) * micro.photon_prod;
  double prob = 0.0;

  // Loop through each reaction type
  const auto& nuc {data::nuclides[i_nuclide]};
  for (int i = 0; i < nuc->reactions_.size(); ++i) {
    // Evaluate neutron cross section
    const auto& rx = nuc->reactions_[i];
    double xs = rx->xs(micro);

    // if cross section is zero for this reaction, skip it
    if (xs == 0.0)
      continue;

    for (int j = 0; j < rx->products_.size(); ++j) {
      if (rx->products_[j].particle_.is_photon()) {
        // For fission, artificially increase the photon yield to account
        // for delayed photons
        double f = 1.0;
        if (settings::delayed_photon_scaling) {
          if (is_fission(rx->mt_)) {
            if (nuc->prompt_photons_ && nuc->delayed_photons_) {
              double energy_prompt = (*nuc->prompt_photons_)(p.E());
              double energy_delayed = (*nuc->delayed_photons_)(p.E());
              f = (energy_prompt + energy_delayed) / (energy_prompt);
            }
          }
        }

        // add to cumulative probability
        prob += f * (*rx->products_[j].yield_)(p.E()) * xs;

        *i_rx = i;
        *i_product = j;
        if (prob > cutoff)
          return;
      }
    }
  }
}

void absorption(Particle& p, int i_nuclide)
{
  double E_in = p.E();
  Direction u_in = p.u();
  double wgt_in = p.wgt();
  const Reaction* absorption_rx = nullptr;

  if (settings::survival_biasing) {
    // Determine weight absorbed in survival biasing
    const double wgt_absorb = p.wgt() * p.neutron_xs(i_nuclide).absorption /
                              p.neutron_xs(i_nuclide).total;

    if (wgt_absorb > 0.0) {
      absorption_rx = sample_disappearance_reaction(i_nuclide, p);
    }

    // Generate recoil with absorbed weight in survival biasing mode
    if (wgt_absorb > 0.0 && !p.fission()) {
      create_absorption_recoil(
        p, i_nuclide, wgt_absorb, E_in, u_in, wgt_in, absorption_rx);
    }

    // Adjust weight of particle by probability of absorption
    p.wgt() -= wgt_absorb;

    // Score implicit absorption estimate of keff
    if (settings::run_mode == RunMode::EIGENVALUE) {
      p.keff_tally_absorption() += wgt_absorb *
                                   p.neutron_xs(i_nuclide).nu_fission /
                                   p.neutron_xs(i_nuclide).absorption;
    }
  } else {
    // See if disappearance reaction happens
    if (p.neutron_xs(i_nuclide).absorption >
        prn(p.current_seed()) * p.neutron_xs(i_nuclide).total) {
      absorption_rx = sample_disappearance_reaction(i_nuclide, p);
      if (!p.fission() && absorption_rx) {
        p.event_mt() = absorption_rx->mt_;
      }

      // Generate recoil for explicit absorption event
      if (!p.fission()) {
        create_absorption_recoil(
          p, i_nuclide, p.wgt(), E_in, u_in, wgt_in, absorption_rx);
      }

      // Score absorption estimate of keff
      if (settings::run_mode == RunMode::EIGENVALUE) {
        p.keff_tally_absorption() += p.wgt() *
                                     p.neutron_xs(i_nuclide).nu_fission /
                                     p.neutron_xs(i_nuclide).absorption;
      }

      p.wgt() = 0.0;
      p.event() = TallyEvent::ABSORB;
      if (!p.fission() && !absorption_rx) {
        p.event_mt() = N_DISAPPEAR;
      }
    }
  }
}

void scatter(Particle& p, int i_nuclide)
{
  // copy incoming direction
  Direction u_old {p.u()};

  // Get pointer to nuclide and grid index/interpolation factor
  const auto& nuc {data::nuclides[i_nuclide]};
  const auto& micro {p.neutron_xs(i_nuclide)};
  int i_temp = micro.index_temp;

  // For tallying purposes, this routine might be called directly. In that
  // case, we need to sample a reaction via the cutoff variable
  double cutoff = prn(p.current_seed()) * (micro.total - micro.absorption);
  bool sampled = false;

  // Calculate elastic cross section if it wasn't precalculated
  if (micro.elastic == CACHE_INVALID) {
    nuc->calculate_elastic_xs(p);
  }

  double prob = micro.elastic - micro.thermal;
  if (prob > cutoff) {
    // =======================================================================
    // NON-S(A,B) ELASTIC SCATTERING

    // Determine temperature
    double kT = nuc->multipole_ ? p.sqrtkT() * p.sqrtkT() : nuc->kTs_[i_temp];

    // Perform collision physics for elastic scattering
    elastic_scatter(i_nuclide, *nuc->reactions_[0], kT, p);

    p.event_mt() = ELASTIC;
    sampled = true;
  }

  prob = micro.elastic;
  if (prob > cutoff && !sampled) {
    // =======================================================================
    // S(A,B) SCATTERING

    sab_scatter(i_nuclide, micro.index_sab, p);

    p.event_mt() = ELASTIC;
    sampled = true;
  }

  if (!sampled) {
    // =======================================================================
    // INELASTIC SCATTERING

    int n = nuc->index_inelastic_scatter_.size();
    int i = 0;
    for (int j = 0; j < n && prob < cutoff; ++j) {
      i = nuc->index_inelastic_scatter_[j];

      // add to cumulative probability
      prob += nuc->reactions_[i]->xs(micro);
    }

    // Perform collision physics for inelastic scattering
    const auto& rx {nuc->reactions_[i]};
    inelastic_scatter(*nuc, *rx, p);
    p.event_mt() = rx->mt_;
  }

  // Set event component
  p.event() = TallyEvent::SCATTER;

  // Sample new outgoing angle for isotropic-in-lab scattering
  const auto& mat {model::materials[p.material()]};
  if (!mat->p0_.empty()) {
    int i_nuc_mat = mat->mat_nuclide_index_[i_nuclide];
    if (mat->p0_[i_nuc_mat]) {
      // Sample isotropic-in-lab outgoing direction
      p.u() = isotropic_direction(p.current_seed());
      p.mu() = u_old.dot(p.u());
    }
  }
}

void elastic_scatter(int i_nuclide, const Reaction& rx, double kT, Particle& p)
{
  // get pointer to nuclide
  const auto& nuc {data::nuclides[i_nuclide]};
  Direction u_in = p.u();

  double vel = std::sqrt(p.E());
  double awr = nuc->awr_;

  // Neutron velocity in LAB
  Direction v_n = vel * p.u();

  // Sample velocity of target nucleus
  Direction v_t {};
  if (!p.neutron_xs(i_nuclide).use_ptable) {
    v_t = sample_target_velocity(*nuc, p.E(), p.u(), v_n,
      p.neutron_xs(i_nuclide).elastic, kT, p.current_seed());
  }

  // Velocity of center-of-mass
  Direction v_cm = (v_n + awr * v_t) / (awr + 1.0);

  // Transform to CM frame
  v_n -= v_cm;

  // Find speed of neutron in CM
  vel = v_n.norm();

  // Sample scattering angle, checking if angle distribution is present (assume
  // isotropic otherwise)
  double mu_cm;
  auto& d = rx.products_[0].distribution_[0];
  auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
  if (!d_->angle().empty()) {
    mu_cm = d_->angle().sample(p.E(), p.current_seed());
  } else {
    mu_cm = uniform_distribution(-1., 1., p.current_seed());
  }

  // Determine direction cosines in CM
  Direction u_cm = v_n / vel;

  // Rotate neutron velocity vector to new angle -- note that the speed of the
  // neutron in CM does not change in elastic scattering. However, the speed
  // will change when we convert back to LAB
  v_n = vel * rotate_angle(u_cm, mu_cm, nullptr, p.current_seed());

  // Transform back to LAB frame
  v_n += v_cm;

  double E_in = p.E();
  p.E() = v_n.dot(v_n);
  vel = std::sqrt(p.E());

  // compute cosine of scattering angle in LAB frame by taking dot product of
  // neutron's pre- and post-collision angle
  p.mu() = p.u().dot(v_n) / vel;

  // Set energy and direction of particle in LAB frame
  p.u() = v_n / vel;

  // Because of floating-point roundoff, it may be possible for mu_lab to be
  // outside of the range [-1,1). In these cases, we just set mu_lab to exactly
  // -1 or 1
  if (std::abs(p.mu()) > 1.0)
    p.mu() = std::copysign(1.0, p.mu());

  // Generate recoil
  if (settings::recoil_production) {
    Direction p_recoil =
      neutron_momentum(E_in, u_in) - neutron_momentum(p.E(), p.u());
    create_recoil_secondary(p, *nuc, p.wgt(), p_recoil, nuc->particle_type());
  }
}

void sab_scatter(int i_nuclide, int i_sab, Particle& p)
{
  // Determine temperature index
  const auto& micro {p.neutron_xs(i_nuclide)};
  int i_temp = micro.index_temp_sab;

  // Sample energy and angle
  double E_out;
  data::thermal_scatt[i_sab]->data_[i_temp].sample(
    micro, p.E(), &E_out, &p.mu(), p.current_seed());

  // Set energy to outgoing, change direction of particle
  p.E() = E_out;
  p.u() = rotate_angle(p.u(), p.mu(), nullptr, p.current_seed());
}

Direction sample_target_velocity(const Nuclide& nuc, double E, Direction u,
  Direction v_neut, double xs_eff, double kT, uint64_t* seed)
{
  // check if nuclide is a resonant scatterer
  ResScatMethod sampling_method;
  if (nuc.resonant_) {

    // sampling method to use
    sampling_method = settings::res_scat_method;

    // upper resonance scattering energy bound (target is at rest above this E)
    if (E > settings::res_scat_energy_max) {
      return {};

      // lower resonance scattering energy bound (should be no resonances below)
    } else if (E < settings::res_scat_energy_min) {
      sampling_method = ResScatMethod::cxs;
    }

    // otherwise, use free gas model
  } else {
    if (E >= settings::free_gas_threshold * kT && nuc.awr_ > 1.0) {
      return {};
    } else {
      sampling_method = ResScatMethod::cxs;
    }
  }

  // use appropriate target velocity sampling method
  switch (sampling_method) {
  case ResScatMethod::cxs:

    // sample target velocity with the constant cross section (cxs) approx.
    return sample_cxs_target_velocity(nuc.awr_, E, u, kT, seed);

  case ResScatMethod::dbrc:
  case ResScatMethod::rvs: {
    double E_red = std::sqrt(nuc.awr_ * E / kT);
    double E_low = std::pow(std::max(0.0, E_red - 4.0), 2) * kT / nuc.awr_;
    double E_up = (E_red + 4.0) * (E_red + 4.0) * kT / nuc.awr_;

    // find lower and upper energy bound indices
    // lower index
    int i_E_low;
    if (E_low < nuc.energy_0K_.front()) {
      i_E_low = 0;
    } else if (E_low > nuc.energy_0K_.back()) {
      i_E_low = nuc.energy_0K_.size() - 2;
    } else {
      i_E_low =
        lower_bound_index(nuc.energy_0K_.begin(), nuc.energy_0K_.end(), E_low);
    }

    // upper index
    int i_E_up;
    if (E_up < nuc.energy_0K_.front()) {
      i_E_up = 0;
    } else if (E_up > nuc.energy_0K_.back()) {
      i_E_up = nuc.energy_0K_.size() - 2;
    } else {
      i_E_up =
        lower_bound_index(nuc.energy_0K_.begin(), nuc.energy_0K_.end(), E_up);
    }

    if (i_E_up == i_E_low) {
      // Handle degenerate case -- if the upper/lower bounds occur for the same
      // index, then using cxs is probably a good approximation
      return sample_cxs_target_velocity(nuc.awr_, E, u, kT, seed);
    }

    if (sampling_method == ResScatMethod::dbrc) {
      // interpolate xs since we're not exactly at the energy indices
      double xs_low = nuc.elastic_0K_[i_E_low];
      double m = (nuc.elastic_0K_[i_E_low + 1] - xs_low) /
                 (nuc.energy_0K_[i_E_low + 1] - nuc.energy_0K_[i_E_low]);
      xs_low += m * (E_low - nuc.energy_0K_[i_E_low]);
      double xs_up = nuc.elastic_0K_[i_E_up];
      m = (nuc.elastic_0K_[i_E_up + 1] - xs_up) /
          (nuc.energy_0K_[i_E_up + 1] - nuc.energy_0K_[i_E_up]);
      xs_up += m * (E_up - nuc.energy_0K_[i_E_up]);

      // get max 0K xs value over range of practical relative energies
      double xs_max = *std::max_element(
        &nuc.elastic_0K_[i_E_low + 1], &nuc.elastic_0K_[i_E_up + 1]);
      xs_max = std::max({xs_low, xs_max, xs_up});

      while (true) {
        double E_rel;
        Direction v_target;
        while (true) {
          // sample target velocity with the constant cross section (cxs)
          // approx.
          v_target = sample_cxs_target_velocity(nuc.awr_, E, u, kT, seed);
          Direction v_rel = v_neut - v_target;
          E_rel = v_rel.dot(v_rel);
          if (E_rel < E_up)
            break;
        }

        // perform Doppler broadening rejection correction (dbrc)
        double xs_0K = nuc.elastic_xs_0K(E_rel);
        double R = xs_0K / xs_max;
        if (prn(seed) < R)
          return v_target;
      }

    } else if (sampling_method == ResScatMethod::rvs) {
      // interpolate xs CDF since we're not exactly at the energy indices
      // cdf value at lower bound attainable energy
      double cdf_low = 0.0;
      if (E_low > nuc.energy_0K_.front()) {
        double m = (nuc.xs_cdf_[i_E_low + 1] - nuc.xs_cdf_[i_E_low]) /
                   (nuc.energy_0K_[i_E_low + 1] - nuc.energy_0K_[i_E_low]);
        cdf_low = nuc.xs_cdf_[i_E_low] + m * (E_low - nuc.energy_0K_[i_E_low]);
      }

      // cdf value at upper bound attainable energy
      double m = (nuc.xs_cdf_[i_E_up + 1] - nuc.xs_cdf_[i_E_up]) /
                 (nuc.energy_0K_[i_E_up + 1] - nuc.energy_0K_[i_E_up]);
      double cdf_up = nuc.xs_cdf_[i_E_up] + m * (E_up - nuc.energy_0K_[i_E_up]);

      while (true) {
        // directly sample Maxwellian
        double E_t = -kT * std::log(prn(seed));

        // sample a relative energy using the xs cdf
        double cdf_rel = cdf_low + prn(seed) * (cdf_up - cdf_low);
        int i_E_rel = lower_bound_index(nuc.xs_cdf_.begin() + i_E_low,
          nuc.xs_cdf_.begin() + i_E_up + 2, cdf_rel);
        double E_rel = nuc.energy_0K_[i_E_low + i_E_rel];
        double m = (nuc.xs_cdf_[i_E_low + i_E_rel + 1] -
                     nuc.xs_cdf_[i_E_low + i_E_rel]) /
                   (nuc.energy_0K_[i_E_low + i_E_rel + 1] -
                     nuc.energy_0K_[i_E_low + i_E_rel]);
        E_rel += (cdf_rel - nuc.xs_cdf_[i_E_low + i_E_rel]) / m;

        // perform rejection sampling on cosine between
        // neutron and target velocities
        double mu = (E_t + nuc.awr_ * (E - E_rel)) /
                    (2.0 * std::sqrt(nuc.awr_ * E * E_t));

        if (std::abs(mu) < 1.0) {
          // set and accept target velocity
          E_t /= nuc.awr_;
          return std::sqrt(E_t) * rotate_angle(u, mu, nullptr, seed);
        }
      }
    }
  } // case RVS, DBRC
  } // switch (sampling_method)

  UNREACHABLE();
}

Direction sample_cxs_target_velocity(
  double awr, double E, Direction u, double kT, uint64_t* seed)
{
  double beta_vn = std::sqrt(awr * E / kT);
  double alpha = 1.0 / (1.0 + std::sqrt(PI) * beta_vn / 2.0);

  double beta_vt_sq;
  double mu;
  while (true) {
    // Sample two random numbers
    double r1 = prn(seed);
    double r2 = prn(seed);

    if (prn(seed) < alpha) {
      // With probability alpha, we sample the distribution p(y) =
      // y*e^(-y). This can be done with sampling scheme C45 from the Monte
      // Carlo sampler

      beta_vt_sq = -std::log(r1 * r2);

    } else {
      // With probability 1-alpha, we sample the distribution p(y) = y^2 *
      // e^(-y^2). This can be done with sampling scheme C61 from the Monte
      // Carlo sampler

      double c = std::cos(PI / 2.0 * prn(seed));
      beta_vt_sq = -std::log(r1) - std::log(r2) * c * c;
    }

    // Determine beta * vt
    double beta_vt = std::sqrt(beta_vt_sq);

    // Sample cosine of angle between neutron and target velocity
    mu = uniform_distribution(-1., 1., seed);

    // Determine rejection probability
    double accept_prob =
      std::sqrt(beta_vn * beta_vn + beta_vt_sq - 2 * beta_vn * beta_vt * mu) /
      (beta_vn + beta_vt);

    // Perform rejection sampling on vt and mu
    if (prn(seed) < accept_prob)
      break;
  }

  // Determine speed of target nucleus
  double vt = std::sqrt(beta_vt_sq * kT / awr);

  // Determine velocity vector of target nucleus based on neutron's velocity
  // and the sampled angle between them
  return vt * rotate_angle(u, mu, nullptr, seed);
}

void sample_fission_neutron(
  int i_nuclide, const Reaction& rx, SourceSite* site, Particle& p)
{
  // Get attributes of particle
  double E_in = p.E();
  uint64_t* seed = p.current_seed();

  // Determine total nu, delayed nu, and delayed neutron fraction
  const auto& nuc {data::nuclides[i_nuclide]};
  double nu_t = nuc->nu(E_in, Nuclide::EmissionMode::total);
  double nu_d = nuc->nu(E_in, Nuclide::EmissionMode::delayed);
  double beta = nu_d / nu_t;

  if (prn(seed) < beta) {
    // ====================================================================
    // DELAYED NEUTRON SAMPLED

    // sampled delayed precursor group
    double xi = prn(seed) * nu_d;
    double prob = 0.0;
    int group;
    for (group = 1; group < nuc->n_precursor_; ++group) {
      // determine delayed neutron precursor yield for group j
      double yield = (*rx.products_[group].yield_)(E_in);

      // Check if this group is sampled
      prob += yield;
      if (xi < prob)
        break;
    }

    // if the sum of the probabilities is slightly less than one and the
    // random number is greater, j will be greater than nuc %
    // n_precursor -- check for this condition
    group = std::min(group, nuc->n_precursor_);

    // set the delayed group for the particle born from fission
    site->delayed_group = group;

    // Sample time of emission based on decay constant of precursor
    double decay_rate = rx.products_[site->delayed_group].decay_rate_;
    site->time -= std::log(prn(p.current_seed())) / decay_rate;

  } else {
    // ====================================================================
    // PROMPT NEUTRON SAMPLED

    // set the delayed group for the particle born from fission to 0
    site->delayed_group = 0;
  }

  // sample from prompt neutron energy distribution
  int n_sample = 0;
  double mu;
  while (true) {
    rx.products_[site->delayed_group].sample(E_in, site->E, mu, seed);

    // resample if energy is greater than maximum neutron energy
    int neutron = ParticleType::neutron().transport_index();
    if (site->E < data::energy_max[neutron])
      break;

    // check for large number of resamples
    ++n_sample;
    if (n_sample == MAX_SAMPLE) {
      // particle_write_restart(p)
      fatal_error("Resampled energy distribution maximum number of times "
                  "for nuclide " +
                  nuc->name_);
    }
  }

  // Sample azimuthal angle uniformly in [0, 2*pi) and assign angle
  site->u = rotate_angle(p.u(), mu, nullptr, seed);
}

void inelastic_scatter(const Nuclide& nuc, const Reaction& rx, Particle& p)
{
  Direction u_in = p.u();
  double E_in = p.E();

  auto sample_neutron_out = [&](
                              double& E_out, double& mu_out, Direction& u_out) {
    rx.products_[0].sample(E_in, E_out, mu_out, p.current_seed());

    if (rx.scatter_in_cm_) {
      double E_cm = E_out;
      double A = nuc.awr_;
      E_out =
        E_cm + (E_in + 2.0 * mu_out * (A + 1.0) * std::sqrt(E_in * E_cm)) /
                 ((A + 1.0) * (A + 1.0));
      mu_out = mu_out * std::sqrt(E_cm / E_out) +
               1.0 / (A + 1.0) * std::sqrt(E_in / E_out);
    }

    if (std::abs(mu_out) > 1.0) {
      mu_out = std::copysign(1.0, mu_out);
    }

    u_out = rotate_angle(u_in, mu_out, nullptr, p.current_seed());
  };

  double E_out;
  double mu_out;
  Direction u_out;
  sample_neutron_out(E_out, mu_out, u_out);

  p.E() = E_out;
  p.mu() = mu_out;
  p.u() = u_out;

  // evaluate yield
  double yield = (*rx.products_[0].yield_)(E_in);
  if (std::floor(yield) == yield && yield > 0) {
    // If yield is integral, create exactly that many secondary particles
    for (int i = 0; i < static_cast<int>(std::round(yield)) - 1; ++i) {
      p.create_secondary(p.wgt(), p.u(), p.E(), ParticleType::neutron());
    }
  } else {
    // Otherwise, change weight of particle based on yield
    p.wgt() *= yield;
  }

  if (!settings::recoil_production) {
    return;
  }

  Direction p_recoil = neutron_momentum(E_in, u_in);
  double emitted_kinetic = 0.0;

  int n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha;
  bool have_counts = emitted_particle_counts(
    rx.mt_, n_neutron, n_proton, n_deuteron, n_triton, n_he3, n_alpha);

  bool has_unmodeled_charged_products =
    have_counts && n_neutron > 0 &&
    (n_proton + n_deuteron + n_triton + n_he3 + n_alpha > 0);
  if (has_unmodeled_charged_products && settings::recoil.missing_products !=
                                          RecoilMissingProducts::neutron_only) {
    static bool warned_missing_products {false};
    if (!warned_missing_products) {
      warning("Recoil setting missing_products mode is not implemented for "
              "charged-particle channels without MF=6 in this build; "
              "falling back to neutron_only.");
      warned_missing_products = true;
    }
  }

  if (settings::recoil.multi_neutron_mode ==
      RecoilMultiNeutronMode::one_particle) {
    p_recoil -= neutron_momentum(E_out, u_out);
    emitted_kinetic += E_out;
  } else if (settings::recoil.multi_neutron_mode ==
             RecoilMultiNeutronMode::duplicate_as_transport) {
    double multiplicity = std::max(0.0, yield);
    p_recoil -= multiplicity * neutron_momentum(E_out, u_out);
    emitted_kinetic += multiplicity * E_out;
  } else {
    int n_secondary = 0;
    if (std::floor(yield) == yield && yield > 0.0) {
      n_secondary = static_cast<int>(std::round(yield)) - 1;
    }
    p_recoil -= neutron_momentum(E_out, u_out);
    emitted_kinetic += E_out;

    if (n_secondary > 0) {
      for (int i = 0; i < n_secondary; ++i) {
        double E_extra;
        double mu_extra;
        Direction u_extra;
        sample_neutron_out(E_extra, mu_extra, u_extra);
        p_recoil -= neutron_momentum(E_extra, u_extra);
        emitted_kinetic += E_extra;
      }
    } else if (yield > 1.0) {
      double extra = yield - 1.0;
      p_recoil -= extra * neutron_momentum(E_out, u_out);
      emitted_kinetic += extra * E_out;
    }
  }

  check_recoil_sanity(nuc, rx, E_in, emitted_kinetic);
  create_recoil_secondary(
    p, nuc, p.wgt(), p_recoil, residual_particle_type(nuc, rx.mt_));
}

void sample_secondary_photons(Particle& p, int i_nuclide)
{
  // Sample the number of photons produced
  double y_t =
    p.neutron_xs(i_nuclide).photon_prod / p.neutron_xs(i_nuclide).total;
  double photon_wgt = p.wgt();
  int y = 1;

  if (settings::use_decay_photons) {
    // For decay photons, sample a single photon and modify the weight
    if (y_t <= 0.0)
      return;
    photon_wgt *= y_t;
  } else {
    // For prompt photons, sample an integral number of photons with weight
    // equal to the neutron's weight
    y = static_cast<int>(y_t);
    if (prn(p.current_seed()) <= y_t - y)
      ++y;
  }

  // Sample each secondary photon
  for (int i = 0; i < y; ++i) {
    // Sample the reaction and product
    int i_rx;
    int i_product;
    sample_photon_product(i_nuclide, p, &i_rx, &i_product);

    // Sample the outgoing energy and angle
    auto& rx = data::nuclides[i_nuclide]->reactions_[i_rx];
    double E;
    double mu;
    rx->products_[i_product].sample(p.E(), E, mu, p.current_seed());

    // Sample the new direction
    Direction u = rotate_angle(p.u(), mu, nullptr, p.current_seed());

    // In a k-eigenvalue simulation, it's necessary to provide higher weight to
    // secondary photons from non-fission reactions to properly balance energy
    // release and deposition. See D. P. Griesheimer, S. J. Douglass, and M. H.
    // Stedry, "Self-consistent energy normalization for quasistatic reactor
    // calculations", Proc. PHYSOR, Cambridge, UK, Mar 29-Apr 2, 2020.
    double wgt = photon_wgt;
    if (settings::run_mode == RunMode::EIGENVALUE && !is_fission(rx->mt_)) {
      wgt *= simulation::keff;
    }

    // Create the secondary photon
    bool created_photon = p.create_secondary(wgt, u, E, ParticleType::photon());

    // Pre-add photon energy to pht_storage so pht_secondary_particles()
    // subtraction results in net zero
    if (created_photon && !model::active_pulse_height_tallies.empty()) {
      auto it = std::find(model::pulse_height_cells.begin(),
        model::pulse_height_cells.end(), p.lowest_coord().cell());
      if (it != model::pulse_height_cells.end()) {
        int index = std::distance(model::pulse_height_cells.begin(), it);
        p.pht_storage()[index] += E;
      }
    }

    // Tag secondary particle with parent nuclide
    if (created_photon && settings::use_decay_photons) {
      p.local_secondary_bank().back().parent_nuclide =
        rx->products_[i_product].parent_nuclide_;
    }
  }
}

} // namespace openmc
