//! \file recoil.h
//! \brief Production of recoil nuclei and emitted light ions

#ifndef OPENMC_RECOIL_H
#define OPENMC_RECOIL_H

#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/reaction.h"

namespace openmc {

//==============================================================================
//! Recoil (primary knock-on atom) production
//!
//! When \ref settings::recoil_production is enabled, each continuous-energy
//! neutron collision creates additional entries in the local secondary bank
//! describing the recoil nucleus and, optionally, the light ions that
//! the reaction emits. These entries are *production records*: they are scored
//! by ParticleProductionFilter and are not transported.
//!
//! Every model here is built on the same nonrelativistic momentum balance,
//!
//! \f[ \mathbf{p}_R = \mathbf{p}_{n,\text{in}} + \mathbf{p}_{T,\text{in}}
//!                  - \sum_i \mathbf{p}_i, \qquad
//!     E_R = \frac{|\mathbf{p}_R|^2}{2 M_R}, \f]
//!
//! where the sum runs over every emitted particle. Products that OpenMC
//! samples from evaluated data (the outgoing neutron, capture photons) are used
//! exactly as sampled; products for which the library carries no distribution
//! (light charged particles, additional neutrons of a multiplicity > 1 channel)
//! are modelled subject to the event's remaining energy budget so that
//! \f$\sum_i E_i + E_R \le E_\text{in} + Q\f$ holds by construction.
//==============================================================================

namespace recoil {

//! Charge and mass number of a nuclide or light ion
struct AtomicNumbers {
  int Z {0};
  int A {0};
};

//==============================================================================
//! \name Event diagnostics
//!
//! Every path that declines to produce a recoil increments one of these. The
//! feature's old failure mode was a record that looked ordinary and was wrong
//! -- a recoil labelled with a nuclide whose momentum balance omitted a
//! product that could not be sampled. The new failure mode is a record that is
//! absent, which is only an improvement if it can be counted.
//! @{
//==============================================================================

enum class RecoilCounter {
  banked,              //!< a complete, consistent production record was made
  incomplete_emission, //!< an exit-channel product could not be emitted
  no_budget,           //!< no trustworthy energy release for the event
  unknown_channel,     //!< the exit channel does not follow from the MT
  unbankable,          //!< the recoil momentum or mass was unusable
  size
};

//! Number of events counted in one bin since the last reset
int64_t counter(RecoilCounter which);

//! Zero every counter. Intended for tests and for the start of a run.
void reset_counters();

//! @}

//==============================================================================
//! \name Kinematic contract
//!
//! Every energy budget, endpoint and excitation in this file comes from these
//! five expressions, so that the transport kernel and the offline calibration
//! cannot disagree about what "the energy available to this channel" means.
//! The Python half is \c recoil.kinematics in the analysis repository and both
//! are checked against \c tests/data/kinematics_fixture.json.
//!
//! Masses are **nuclear**, obtained from the atomic mass as
//! \f$M_\text{nuc}(Z,A) = M_\text{atom}(Z,A) - Z m_e\f$ with the five light
//! ions taken from their CODATA values. Mixing conventions is the error to
//! avoid: a Q value built from atomic targets and nuclear light ions is wrong
//! by \f$Z_b m_e c^2\f$, which is 1.02 MeV for an alpha channel.
//! @{
//==============================================================================

//! Nuclear rest mass of a nuclide or light ion in [eV]
//!
//! \return zero when the mass is not tabulated. Callers must check: a
//!         fabricated mass produces a Q value wrong by tens of MeV, which is
//!         worse than producing no recoil at all.
double nuclear_mass_ev(AtomicNumbers za);

//! Nuclear mass less \f$A\f$ mass units, in [eV]
//!
//! Q values are built from these rather than from the masses themselves. A
//! W-184 channel differences four numbers of order 1.7e11 eV to reach one of
//! order 1e6, and a double carries about sixteen digits, so the direct
//! subtraction returns a Q good to only ten. The mass numbers cancel
//! identically, so working in excesses of order 1e7 eV recovers five digits.
double mass_excess_ev(AtomicNumbers za);

//! Ground-state mass budget of an exit channel in [eV]
//!
//! \f[ Q_M = \left[M_T + m_n - M_D - \sum_j m_j\right] c^2 \f]
//!
//! where the daughter \f$D\f$ is whatever the emitted particles leave behind.
//! This is the rest-mass energy release, which is what an event budget needs.
//! It is *not* what ENDF MF=3 \c QI holds for a continuum, level-range or
//! summation channel: there \c QI is a threshold-setting value, and across
//! 418 evaluated channels in four libraries it sits up to 7.0 MeV below
//! \f$Q_M\f$.
//!
//! \param[in] target     Charge and mass number of the target
//! \param[in] emitted    Charge and mass number of each emitted particle
//! \param[in] n_emitted  Number of entries in \p emitted
//! \param[out] ok        False if any required mass is missing or the channel
//!                       leaves an impossible daughter
double mass_difference_q(
  AtomicNumbers target, const AtomicNumbers* emitted, int n_emitted, bool& ok);

//! Energy available in the compound system's rest frame in [eV]
//!
//! \f[ U_0 = E_\text{in}\frac{M_T}{M_T+m_n} + Q \f]
//!
//! The centre of mass carries \f$E_\text{in} m_n/(M_T+m_n)\f$ that no exit
//! channel can spend, so \f$E_\text{in}+Q\f$ overstates the budget. The excess
//! is a few percent for a mid-mass target and unbounded near threshold.
double entrance_internal_energy(double E_in, double m_target, double q);

//! Largest centre-of-mass energy ion \p m_b can take from internal energy \p u
//!
//! \f[ E_{b,\max} = u\,\frac{M_D}{m_b + M_D} \f]
double two_body_endpoint(double u, double m_b, double m_d);

//! Internal energy left in the daughter after a two-body emission in [eV]
//!
//! \f[ E_x = u - E_b^\text{cm}\left(1 + \frac{m_b}{M_D}\right) \f]
//!
//! Negative only if \p e_cm exceeded two_body_endpoint(), which every caller
//! should assert rather than clamp.
double residual_excitation(double u, double e_cm, double m_b, double m_d);

//! Endpoint of the channel that emits \p ion alone from \p target, in [eV]
//!
//! Sets the shape of the modelled spectrum: the \f$(1-E/E_\text{max})^\nu\f$
//! factor is the level density of the daughter that ion would leave if nothing
//! else were emitted. Built from masses and the entrance CM energy, never from
//! an evaluated Q, because the MT one would reach for is often a level range.
//!
//! \return zero if a required mass is missing
double shape_endpoint(double E_in, AtomicNumbers target, AtomicNumbers ion);

//! @}

//! Recoil products following an elastic scattering event
//!
//! The recoil energy is the energy the collision transfers to the target, so
//! the target's own thermal momentum is not included; see the implementation
//! for why.
//!
//! \param[in,out] p        Colliding neutron, after the outgoing state is set
//! \param[in] nuc          Target nuclide
//! \param[in] E_in         Incident neutron energy in [eV]
//! \param[in] u_in         Incident neutron direction
//! \param[in] E_out        Outgoing neutron energy in [eV]
//! \param[in] u_out        Outgoing neutron direction
void from_elastic(Particle& p, const Nuclide& nuc, double E_in, Direction u_in,
  double E_out, Direction u_out);

//! Recoil products following a reaction that emits neutrons
//!
//! The sampled outgoing neutron is used exactly as transported. Additional
//! neutrons required by the reaction multiplicity are sampled independently
//! from the same evaluated distribution, and light charged particles in the
//! exit channel are modelled, both within the event's energy budget.
//!
//! \param[in,out] p    Colliding neutron, after the outgoing state is set
//! \param[in] nuc      Target nuclide
//! \param[in] rx       Reaction that was sampled
//! \param[in] wgt      Collision weight before any neutron-yield scaling
//! \param[in] E_in     Incident neutron energy in [eV]
//! \param[in] u_in     Incident neutron direction
//! \param[in] E_out    Sampled outgoing neutron energy in [eV]
//! \param[in] u_out    Sampled outgoing neutron direction
//! \param[in] yield    Neutron yield of the reaction at \p E_in
void from_inelastic(Particle& p, const Nuclide& nuc, const Reaction& rx,
  double wgt, double E_in, Direction u_in, double E_out, Direction u_out,
  double yield);

//! Recoil products following a disappearance (absorption) reaction
//!
//! \param[in,out] p         Colliding neutron
//! \param[in] i_nuclide     Index of the target nuclide
//! \param[in] weight        Weight to assign to the recoil products
//! \param[in] E_in          Incident neutron energy in [eV]
//! \param[in] u_in          Incident neutron direction
//! \param[in] rx            Reaction that was sampled; no products are made
//!                          if this is nullptr
void from_absorption(Particle& p, int i_nuclide, double weight, double E_in,
  Direction u_in, const Reaction* rx);

//! Identity of the recoil left by a reaction
//!
//! \param[in] nuc  Target nuclide
//! \param[in] mt   Reaction MT number
//! \return Recoil nucleus in its ground state; the target itself if the exit
//!         channel cannot be determined from \p mt
ParticleType recoil_particle_type(const Nuclide& nuc, int mt);

//! Rest mass of a particle in [eV]
//!
//! Returns zero for photons and falls back to \f$A\f$ atomic mass units for
//! nuclides that are absent from the tabulated mass table.
double particle_mass_ev(ParticleType type);

//! Centre-of-mass kinetic energy of a light ion emitted from an excited system
//!
//! Samples the empirical evaporation spectrum described in
//! \ref light_ion_pdf() by rejection.
//!
//! \param[in] E_max     Endpoint that sets the spectrum shape in [eV], namely
//!                      the kinematic maximum of the channel that emits this
//!                      ion alone
//! \param[in] E_limit   Largest energy this event can afford in [eV]; the
//!                      spectrum is sampled truncated to [0, E_limit]
//! \param[in] Z_b,A_b   Charge and mass number of the emitted ion
//! \param[in] Z_d,A_d   Charge and mass number of the daughter nucleus
//! \param[in] seed      Pseudorandom number seed pointer
double sample_light_ion_energy(double E_max, double E_limit, int Z_b, int A_b,
  int Z_d, int A_d, uint64_t* seed);

//! Kalbach-Mann slope parameter \f$a\f$ from the 1988 systematics
//!
//! Exposed so that angular candidates can be scored through this exact code
//! path rather than a reimplementation of it.
//!
//! \param[in] E_in     Incident neutron energy in [eV]
//! \param[in] E_cm     Emitted particle centre-of-mass energy in [eV]
//! \param[in] emitted  Charge and mass number of the emitted particle
//! \param[in] Z_t,A_t  Charge and mass number of the target
double kalbach_slope(
  double E_in, double E_cm, AtomicNumbers emitted, int Z_t, int A_t);

//! Kalbach-Mann slope parameter, taking the target from a Nuclide
double kalbach_slope(
  double E_in, double E_cm, AtomicNumbers emitted, const Nuclide& nuc);

//! Sample \f$f(\mu) \propto \cosh(a\mu) + r\sinh(a\mu)\f$ by exact inversion
//!
//! \param[in] slope  Kalbach slope \f$a\f$
//! \param[in] r      Pre-equilibrium fraction, in [0, 1]
//! \param[in] seed   Pseudorandom number seed pointer
double sample_kalbach_mu(double slope, double r, uint64_t* seed);

//! Unnormalized light-ion emission spectrum
//!
//! \f[ P(E) \propto E\, T_C(E) \left(1 - E/E_\text{max}\right)^{\nu}, \qquad
//!     T_C(E) = \left[1 + e^{2\pi g (\eta(E) - \eta(V_C))}\right]^{-1}, \f]
//!
//! with the Sommerfeld parameter and Coulomb barrier
//!
//! \f[ \eta(E) = \frac{Z_b Z_d}{137.036}\sqrt{\frac{\mu c^2}{2E}}, \qquad
//!     V_C = \frac{1.44\,\text{MeV fm}\, Z_b Z_d}
//!                {r_0 (A_b^{1/3} + A_d^{1/3})}. \f]
//!
//! The factor \f$E\,T_C(E)\f$ stands in for the inverse-reaction cross section
//! of a Weisskopf-Ewing evaporation spectrum and
//! \f$(1-E/E_\text{max})^{\nu}\f$ for the level density of the recoil, with
//! \f$\nu = n - 2\f$ for an \f$n\f$-exciton state.
//!
//! The transmission is a WKB Coulomb penetrability normalized so that
//! \f$T_C = 1/2\f$ at the barrier top. Because \f$\eta \propto E^{-1/2}\f$ the
//! decay constant itself grows as \f$E\f$ falls, which is the widening of the
//! barrier at lower energy; a transmission with a fixed diffuseness falls at
//! one rate everywhere and cannot reproduce it. The effective radius
//! \f$r_0\f$, the WKB strength \f$g\f$ and the exponent \f$\nu\f$ are
//! calibrated against evaluated ENDF MF=6 charged-particle spectra rather than
//! derived from an optical model; see the recoil section of the user
//! documentation.
//!
//! \note The transmission is evaluated as \f$-\operatorname{softplus}\f$ of
//!       the exponent, with no bound on it. An earlier implementation clamped
//!       the exponent at \f$\pm 60\f$ because \f$1/(1+e^x)\f$ underflows to
//!       exactly zero below the barrier and leaves nothing to normalize. That
//!       clamp floored the transmission at a constant and flattened the
//!       spectrum to \f$E(1-E/E_\text{max})^\nu\f$ wherever a channel lay
//!       deep below the barrier, which erased the barrier shape from the
//!       helium-3 channels almost entirely.
double light_ion_pdf(
  double E, double E_max, int Z_b, int A_b, int Z_d, int A_d);

//! Calibration constants of the light-ion spectrum
//!
//! There is one model and these are its parameters; the struct groups them so
//! that they are documented and tested in one place, and so that a candidate
//! set can be scored through the overloads below using this exact code path
//! rather than a reimplementation of it. Nothing in the transport kernel passes
//! anything but the defaults.
struct LightIonParams {
  double r0 {1.36093}; //!< effective barrier radius in [fm]
  double g {0.48566};  //!< scales the WKB barrier exponent
  double nu {1.18221}; //!< endpoint exponent of the level-density factor
};

//! Calibration constants of the light-ion angular distribution
//!
//! The pre-equilibrium fraction of the Kalbach form. The deployed prescription
//! was \f$r = E/E_{\max,\text{shape}}\f$, which is identically 1 at a
//! ground-state channel because the shape endpoint is *defined* from that
//! channel's Q value, so it asserted purely direct emission regardless of the
//! reaction dynamics. Replaced by a logistic in quantities the transport kernel
//! already has; see the recoil section of the methods documentation.
//!
//! These constants describe a **continuum** channel and are not applied to a
//! named level, which is sampled isotropically in the centre of mass. Nor were
//! they refitted when the corpus they were calibrated against was found to
//! carry an exactly-zero pre-equilibrium fraction over 79% of its probability:
//! measured on the nodes that are not that fill, they already agree with the
//! evaluations to three decimal places in the bands where the fill dominates,
//! and refitting without it moves them by less than the disagreement between
//! libraries. The apparent discrepancy against evaluated MF=6 is documented in
//! the methods section so that it is not mistaken for a defect here.
struct AngularParams {
  double c0 {-8.42909};         //!< constant
  double c1 {4.89255};          //!< coefficient of E/E_max_shape
  double c2 {9.85086};          //!< coefficient of log(1 + E_in / 10 MeV)
  double c3 {-7.85305};         //!< coefficient of A_D^(-1/3)
  double c4 {5.69894};          //!< coefficient of the daughter neutron excess
  double slope_scale {0.99102}; //!< multiplies the Kalbach 1988 slope
};

//! Pre-equilibrium fraction \f$r\f$ of the Kalbach angular distribution
//!
//! \f[ r = \left[1 + e^{-u}\right]^{-1}, \qquad
//!     u = c_0 + c_1\frac{E}{E_{\max,\text{shape}}}
//!         + c_2\ln\!\left(1 + \frac{E_\text{in}}{10\,\text{MeV}}\right)
//!         + c_3 A_D^{-1/3} + c_4\frac{N_D - Z_D}{A_D}. \f]
//!
//! \param[in] E_cm          Emitted ion centre-of-mass energy in [eV]
//! \param[in] E_max_shape   Shape endpoint of the ground-state channel in [eV]
//! \param[in] E_in          Incident neutron energy in [eV]
//! \param[in] daughter      Charge and mass number of the recoil nucleus
//! \param[in] par           Calibration parameters
double kalbach_precompound_fraction(double E_cm, double E_max_shape,
  double E_in, AtomicNumbers daughter, const AngularParams& par = {});

//! Natural logarithm of the unnormalized light-ion emission spectrum
//!
//! \f[ \ln P(E) = \ln E - \operatorname{softplus}
//!     \big(2\pi g[\eta(E)-\eta(V_C)]\big)
//!     + \nu \ln\!\left(1 - E/E_\text{max}\right) \f]
//!
//! This is the primitive; light_ion_pdf() exponentiates it. Sub-barrier the
//! spectrum spans hundreds of decades, so a code that needs relative
//! probabilities there -- the sampler, or a fit -- must work here instead.
//!
//! \return \f$-\infty\f$ outside \f$(0, E_\text{max})\f$
double light_ion_log_pdf(
  double E, double E_max, int Z_b, int A_b, int Z_d, int A_d);

//! Natural logarithm of the spectrum with explicit parameters
double light_ion_log_pdf(double E, double E_max, int Z_b, int A_b, int Z_d,
  int A_d, const LightIonParams& par);

//! Unnormalized light-ion emission spectrum with explicit parameters
//!
//! \f$ P(E) \propto E\,T_C(E)\,(1 - E/E_\text{max})^\nu \f$. Calling this
//! with a default-constructed LightIonParams is identical to the fixed-constant
//! overload above; that equality is asserted in the C++ unit tests.
double light_ion_pdf(double E, double E_max, int Z_b, int A_b, int Z_d, int A_d,
  const LightIonParams& par);

//! Sample the light-ion spectrum with explicit parameters
double sample_light_ion_energy(double E_max, double E_limit, int Z_b, int A_b,
  int Z_d, int A_d, uint64_t* seed, const LightIonParams& par);

} // namespace recoil
} // namespace openmc

#endif // OPENMC_RECOIL_H
