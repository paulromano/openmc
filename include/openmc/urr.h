//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"

#include "xtensor/xtensor.hpp"

#include <unordered_map>
#include <utility>

namespace openmc {

//==============================================================================
//! UrrData contains probability tables for the unresolved resonance range.
//==============================================================================

class UrrData{
public:
  Interpolation interp_;          //!< interpolation type
  int inelastic_flag_;            //!< inelastic competition flag
  int absorption_flag_;           //!< other absorption flag
  bool multiply_smooth_;          //!< multiply by smooth cross section?
  int n_energy_;                  //!< number of energy points
  xt::xtensor<double, 1> energy_; //!< incident energies
  xt::xtensor<double, 3> prob_;   //!< Actual probability tables

  //! \brief Load the URR data from the provided HDF5 group
  explicit UrrData(hid_t group_id);
};

//==============================================================================
// Cross section evaluated at a particular energy and temperature
//==============================================================================

struct URRXS {
  double elastic {0.0};
  double capture {0.0};
  double fission {0.0};
  double competitive {0.0};
  double total {0.0};
};

//==============================================================================
// Resonance ladder spanning the entire unresolved resonance range
//==============================================================================

class ResonanceLadder {
public:
  // Types
  struct Resonance {
    double E; //!< Energy
    int l; //!< Neutron orbital angular momentum
    int j; //!< Total angular momentum
    double gt; //!< Total width
    double gn; //!< Energy-dependent neutron width
    double gg; //!< Radiation width
    double gf; //!< Fission width
    double gx; //!< Competitive width
    double p; //!< Penetration factor
    double s; //!< Shift factor
  };

  // Methods
  URRXS evaluate(double E, double sqrtkT, double target_spin, double awr,
    const Function1D& channel_radius, const Function1D& scattering_radius) const;

  // Data members
  std::vector<Resonance> res_; //!< Sampled resonance parameters
  std::unordered_map<int, std::vector<int>> l_values_;
};

//==============================================================================
// Unresolved resonance parameters
//==============================================================================

class Unresolved {
public:
  // Types, enums
  enum class Case {
    A, B, C
  };

  struct URParameters {
    double E;       //!< Energy
    double avg_d;   //!< Average level spacing
    double df_x;    //!< Degrees of freedom in competetive width distribution
    double df_n;    //!< Degrees of freedom in neutron width distribution
    double df_f;    //!< Degrees of freedom in fission width distribution
    double avg_gx;  //!< Average competitive reaction width
    double avg_gn0; //!< Average reduced neutron width
    double avg_gg;  //!< Average radiation width
    double avg_gf;  //!< Average fission width
  };

  struct SpinSequence {
    int l; //!< Neutron orbital angular momentum
    int j; //!< Total angular momentum
    std::vector<URParameters> params; //!< Unresolved resonance parameters
  };

  // Constructors
  Unresolved(hid_t group);

  // Methods
  ResonanceLadder sample_full_ladder(uint64_t* seed) const;
  ResonanceLadder sample_ladder(double energy, uint64_t* seed) const;

  // Data members
  Case case_; //!< Which of 3 cases
  bool add_to_background_; //!< Whether to add File 3 cross sections
  double energy_min_; //!< Minimum energy of the unresolved resonance range in eV
  double energy_max_; //!< Maximum energy of the unresolved resonance range in eV
  double target_spin_; //!< Intrinsic spin of the target nuclide
  double awr_; //!< Atomic weight ratio
  std::unique_ptr<Function1D> channel_radius_;
  std::unique_ptr<Function1D> scattering_radius_;
  std::vector<double> energy_; //!< Energy at which parameters are tabulated
  std::vector<SpinSequence> ljs_; //!< Unresolved resonance parameters at each (l,j)

private:
  //! Sample a resonance to add to a ladder
  //
  //! \param[in] E  Energy of the resonance
  //! \param[in] E_neutron  Energy at which penetration/shift are calculated at
  //! \param[in] l  Neutron orbital angular momentum
  //! \param[in] j  Total angular momentum
  //! \param[in] p  Average unresolved resonance parameters
  //! \param[in,out] seed  PRNG seed
  //! \return Sampled resonance
  ResonanceLadder::Resonance sample_resonance(double E, double E_neutron,
    int l, int j, const URParameters& p, uint64_t* seed) const;

  //! Linearly interpolate average resonance parameters
  //
  //! \param[in] left  Parameters to the left of energy E
  //! \param[in] right  Parameters to the right of energy E
  //! \param[in] E  Energy to interpolate to
  //! \return Linearly interpolated parameters
  URParameters interpolate_parameters(const URParameters& left,
    const URParameters& right, double E) const;
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Calculate the neutron wave number in center-of-mass system.
//!
//! \param[in] A Ratio of target mass to neutron mass
//! \param[in] E Energy in eV
//! \return Neutron wave number in b^-0.5
double wave_number(double A, double E);

//! Calculate hardsphere phase shift
//!
//! \param[in] l Angular momentum quantum number
//! \param[in] rho Product of the wave number and the channel radius
//! \return Hardsphere phase shift
double phase_shift(int l, double rho);

//! Calculate shift and penetration factors
//!
//! \param[in] l Angular momentum quantum number
//! \param[in] rho Product of the wave number and the channel radius
//! \return Penetration factor and shift factor for given l
std::pair<double, double> penetration_shift(int l, double rho);

//! Sample from a chi-square distribution
//!
//! \param[in] df Number of degrees of freedom
//! \param[in,out] seed PRNG seed
//! \return Sample drawn from a chi-square distribution
double chi_square(int df, uint64_t* seed);

//! Sample resonance width from chi-squared distribution
//
//! \param[in] avg_width Average resonance width
//! \param[in] df Number of degrees of freedom
//! \param[in,out] seed PRNG seed
//! \return Sampled resonance width
double sample_width(double avg_width, double df, uint64_t* seed);

//! Sample resonance spacing from Wigner distribution
//
//! \param[in] avg_spacing Average resonance spacing
//! \param[in,out] seed PRNG seed
//! \return Sampled resonance spacing
double sample_spacing(double avg_spacing, uint64_t* seed);

} // namespace openmc

#endif // OPENMC_URR_H
