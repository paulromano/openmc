//! \brief UrrData information for the unresolved resonance treatment

#ifndef OPENMC_URR_H
#define OPENMC_URR_H

#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/hdf5_interface.h"

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
//! Resonance parameters
//==============================================================================

struct Resonance {
  double energy;
  int l; //!< Neutron orbital angular momentum
  int j; //!< Total angular momentum
  double gt; //!< Total width
  double gn; //!< Energy-dependent neutron width
  double gg; //!< Radiation width
  double gf; //!< Fission width
  double gx; //!< Competitive width
};

//==============================================================================
// Resonance ladder spanning the entire unresolved resonance range
//==============================================================================

class ResonanceLadder {
public:
  // Methods
  double evaluate(double E, double T) const;

  // Data members
  std::vector<Resonance> res_; //!< Sampled resonance parameters
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

  // Constructors
  Unresolved(hid_t group);

  // Methods
  ResonanceLadder sample(double E) const;

  // Data members
  Case case_; //!< Which of 3 cases
  bool add_to_background_; //!< Whether to add File 3 cross sections
  double energy_min_; //!< Minimum energy of the unresolved resonance range in eV
  double energy_max_; //!< Maximum energy of the unresolved resonance range in eV
  double target_spin_; //!< Intrinsic spin of the target nuclide
  std::unique_ptr<Function1D> channel_radius_;
  std::unique_ptr<Function1D> scattering_radius_;
  xt::xtensor<double, 1> energy_; //!< Energy at which parameters are tabulated
  xt::xtensor<double, 2> params_; //!< Unresolved resonance parameters
};

} // namespace openmc

#endif // OPENMC_URR_H
