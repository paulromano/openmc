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
//! path rather than a reimplementation of it. See
//! `jnm-recoil/cpp/verify_angular.cpp`.
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
//! \f[ P(E) \propto E\, T_C(E) \sqrt{1 - E/E_\text{max}}, \qquad
//!     T_C(E) = \left[1 + e^{(V_C - E)/\Delta}\right]^{-1}, \qquad
//!     V_C = \frac{1.44\,\text{MeV fm}\, Z_b Z_d}
//!                {r_0 (A_b^{1/3} + A_d^{1/3})}. \f]
//!
//! The factor \f$E\,T_C(E)\f$ stands in for the inverse-reaction cross section
//! of a Weisskopf-Ewing evaporation spectrum and \f$\sqrt{1-E/E_\text{max}}\f$
//! for the level density of the recoil. The effective barrier radius
//! \f$r_0\f$ and diffuseness \f$\Delta\f$ are calibrated against evaluated
//! ENDF MF=6 charged-particle spectra rather than derived from an optical
//! model; see the recoil section of the user documentation.
double light_ion_pdf(
  double E, double E_max, int Z_b, int A_b, int Z_d, int A_d);

//! Calibration parameters of the light-ion spectrum
//!
//! The deployed values are the defaults. The parameterized overloads below
//! exist so that a candidate parameter set can be scored against evaluated data
//! using this exact code path, rather than against a reimplementation of it in
//! the fitting scripts. See `jnm-recoil/cpp/verify_light_ion.cpp`.
struct LightIonParams {
  double r0 {1.7537};   //!< effective barrier radius in [fm]
  double delta {0.8e6}; //!< barrier diffuseness in [eV]
  double nu {1.2110};   //!< endpoint exponent of the level-density factor
};

//! Calibration parameters of the light-ion angular distribution
//!
//! The pre-equilibrium fraction of the Kalbach form. The deployed prescription
//! was \f$r = E/E_{\max,\text{shape}}\f$, which is identically 1 at a
//! ground-state channel because the shape endpoint is *defined* from that
//! channel's Q value, so it asserted purely direct emission regardless of the
//! reaction dynamics. Replaced by a logistic in quantities the transport kernel
//! already has; see the recoil section of the methods documentation.
struct AngularParams {
  double c0 {-7.72232};         //!< constant
  double c1 {3.78801};          //!< coefficient of E/E_max_shape
  double c2 {10.59511};         //!< coefficient of log(1 + E_in / 10 MeV)
  double c3 {-18.91239};        //!< coefficient of A_D^(-1/3)
  double c4 {8.76885};          //!< coefficient of the daughter neutron excess
  double slope_scale {1.01913}; //!< multiplies the Kalbach 1988 slope
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
