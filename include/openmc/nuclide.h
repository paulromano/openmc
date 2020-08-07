//! \file nuclide.h
//! \brief Nuclide type and other associated types/data

#ifndef OPENMC_NUCLIDE_H
#define OPENMC_NUCLIDE_H

#include <array>
#include <memory> // for unique_ptr
#include <unordered_map>
#include <vector>

#include <gsl/gsl>
#include <hdf5.h>

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/particle.h"
#include "openmc/reaction.h"
#include "openmc/reaction_product.h"
#include "openmc/urr.h"
#include "openmc/wmp.h"

namespace openmc {

//==============================================================================
// Data for a nuclide
//==============================================================================

class Nuclide {
public:
  // Types, aliases
  using EmissionMode = ReactionProduct::EmissionMode;
  struct EnergyGrid {
    std::vector<int> grid_index;
    std::vector<double> energy;
  };

  struct InterpResult {
    gsl::index i;  //!< Index in tabulated data
    double f;      //!< Interpolation factor between i and i+1
  };

  // Constructors/destructors
  Nuclide(hid_t group, const std::vector<double>& temperature);
  ~Nuclide();

  //! Initialize logarithmic grid for energy searches
  void init_grid();

  void calculate_xs(int i_sab, int i_log_union, double sab_frac, Particle& p);

  void calculate_sab_xs(int i_sab, double sab_frac, Particle& p);

  // Methods
  double nu(double E, EmissionMode mode, int group=0) const;
  void calculate_elastic_xs(Particle& p) const;

  //! Determines the microscopic 0K elastic cross section at a trial relative
  //! energy used in resonance scattering
  double elastic_xs_0K(double E) const;

  //! \brief Determines cross sections in the unresolved resonance range
  //! from probability tables.
  void calculate_urr_xs(int i_temp, Particle& p) const;

  //! \brief Calculate reaction rate based on group-wise flux distribution
  //
  //! \param[in] MT ENDF MT value for desired reaction
  //! \param[in] temperature Temperature in [K]
  //! \param[in] energy Energy group boundaries in [eV]
  //! \param[in] flux Flux in each energy group (not normalized per eV)
  //! \return Reaction rate
  double collapse_rate(int MT, double temperature, gsl::span<const double> energy,
    gsl::span<const double> flux) const;

  // Data members
  std::string name_; //!< Name of nuclide, e.g. "U235"
  int Z_; //!< Atomic number
  int A_; //!< Mass number
  int metastable_; //!< Metastable state
  double awr_; //!< Atomic weight ratio
  gsl::index index_; //!< Index in the nuclides array

  // Temperature dependent cross section data
  std::vector<double> kTs_; //!< temperatures in eV (k*T)
  std::vector<EnergyGrid> grid_; //!< Energy grid at each temperature
  std::vector<xt::xtensor<double, 2>> xs_; //!< Cross sections at each temperature

  // Multipole data
  std::unique_ptr<WindowedMultipole> multipole_;

  // Fission data
  bool fissionable_ {false}; //!< Whether nuclide is fissionable
  bool has_partial_fission_ {false}; //!< has partial fission reactions?
  std::vector<Reaction*> fission_rx_; //!< Fission reactions
  int n_precursor_ {0}; //!< Number of delayed neutron precursors
  std::unique_ptr<Function1D> total_nu_; //!< Total neutron yield
  std::unique_ptr<Function1D> fission_q_prompt_; //!< Prompt fission energy release
  std::unique_ptr<Function1D> fission_q_recov_; //!< Recoverable fission energy release
  std::unique_ptr<Function1D> prompt_photons_; //!< Prompt photon energy release
  std::unique_ptr<Function1D> delayed_photons_; //!< Delayed photon energy release
  std::unique_ptr<Function1D> fragments_; //!< Fission fragment energy release
  std::unique_ptr<Function1D> betas_; //!< Delayed beta energy release

  // Resonance scattering information
  bool resonant_ {false};
  std::vector<double> energy_0K_;
  std::vector<double> elastic_0K_;
  std::vector<double> xs_cdf_;

  // Unresolved resonance range information
  bool urr_present_ {false};
  int urr_inelastic_ {C_NONE};
  std::vector<UrrData> urr_data_;

  std::vector<std::unique_ptr<Reaction>> reactions_; //!< Reactions
  std::array<size_t, 902> reaction_index_; //!< Index of each reaction
  std::vector<int> index_inelastic_scatter_;

private:
  void create_derived(const Function1D* prompt_photons, const Function1D* delayed_photons);

  InterpResult find_temperature(double T) const;

  static int XS_TOTAL;
  static int XS_ABSORPTION;
  static int XS_FISSION;
  static int XS_NU_FISSION;
  static int XS_PHOTON_PROD;
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Checks for the right version of nuclear data within HDF5 files
void check_data_version(hid_t file_id);

bool multipole_in_range(const Nuclide& nuc, double E);

//==============================================================================
// Global variables
//==============================================================================

namespace data {

// Minimum/maximum transport energy for each particle type. Order corresponds to
// that of the ParticleType enum
extern std::array<double, 2> energy_min;
extern std::array<double, 2> energy_max;

//! Minimum temperature in [K] that nuclide data is available at
extern double temperature_min;

//! Maximum temperature in [K] that nuclide data is available at
extern double temperature_max;

extern std::unordered_map<std::string, int> nuclide_map;
extern std::vector<std::unique_ptr<Nuclide>> nuclides;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void nuclides_clear();

} // namespace openmc

#endif // OPENMC_NUCLIDE_H
