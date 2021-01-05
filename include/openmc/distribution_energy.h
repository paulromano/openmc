//! \file distribution_energy.h
//! Energy distributions that depend on incident particle energy

#ifndef OPENMC_DISTRIBUTION_ENERGY_H
#define OPENMC_DISTRIBUTION_ENERGY_H

#include <vector>

#include "xtensor/xtensor.hpp"
#include "hdf5.h"

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/serialize.h"

namespace openmc {

enum class EnergyDistType {
  DISCRETE_PHOTON,
  LEVEL_INELASTIC,
  CONTINUOUS_TABULAR,
  EVAPORATION,
  MAXWELL,
  WATT
};

class EnergyDistributionFlat {
public:
  EnergyDistributionFlat(EnergyDistType type, const uint8_t* data)
    : type_(type), data_(data) { }

  double sample(double E, uint64_t* seed) const;
private:
  EnergyDistType type_;
  const uint8_t* data_;
};

//===============================================================================
//! Abstract class defining an energy distribution that is a function of the
//! incident energy of a projectile. Each derived type must implement a sample()
//! function that returns a sampled outgoing energy given an incoming energy
//===============================================================================

class EnergyDistribution {
public:
  virtual double sample(double E, uint64_t* seed) const = 0;
  virtual size_t nbytes() const = 0;
  virtual void serialize(DataBuffer& buffer) const = 0;
  virtual ~EnergyDistribution() = default;
};

//===============================================================================
//! Discrete photon energy distribution
//===============================================================================

class DiscretePhoton : public EnergyDistribution {
public:
  explicit DiscretePhoton(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;

  size_t nbytes() const { return 24; }

  void serialize(DataBuffer& buffer) const;
private:
  int primary_flag_; //!< Indicator of whether the photon is a primary or
                     //!< non-primary photon.
  double energy_; //!< Photon energy or binding energy
  double A_; //!< Atomic weight ratio of the target nuclide
};

class DiscretePhotonFlat {
public:
  explicit DiscretePhotonFlat(const uint8_t* data) : data_(data) { }

  double sample(double E, uint64_t* seed) const;
private:
  int primary_flag() const { return *reinterpret_cast<const int*>(data_); }
  double energy() const { return *reinterpret_cast<const double*>(data_ + 4); }
  double A() const { return *reinterpret_cast<const double*>(data_ + 12); }

  const uint8_t* data_;
};

//===============================================================================
//! Level inelastic scattering distribution
//===============================================================================

class LevelInelastic : public EnergyDistribution {
public:
  explicit LevelInelastic(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;

  size_t nbytes() const { return 20; }

  void serialize(DataBuffer& buffer) const;
private:
  double threshold_; //!< Energy threshold in lab, (A + 1)/A * |Q|
  double mass_ratio_; //!< (A/(A+1))^2
};

class LevelInelasticFlat {
public:
  explicit LevelInelasticFlat(const uint8_t* data) : data_(data) { }

  double sample(double E, uint64_t* seed) const;
private:
  double threshold() const { return *reinterpret_cast<const double*>(data_); }
  double mass_ratio() const { return *reinterpret_cast<const double*>(data_ + 8); }

  const uint8_t* data_;
};

//===============================================================================
//! An energy distribution represented as a tabular distribution with histogram
//! or linear-linear interpolation. This corresponds to ACE law 4, which NJOY
//! produces for a number of ENDF energy distributions.
//===============================================================================

class ContinuousTabular : public EnergyDistribution {
public:
  explicit ContinuousTabular(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;

  size_t nbytes() const;

  void serialize(DataBuffer& buffer) const;
private:
  //! Outgoing energy for a single incoming energy
  struct CTTable {
    Interpolation interpolation; //!< Interpolation law
    int n_discrete; //!< Number of of discrete energies
    xt::xtensor<double, 1> e_out; //!< Outgoing energies in [eV]
    xt::xtensor<double, 1> p; //!< Probability density
    xt::xtensor<double, 1> c; //!< Cumulative distribution
  };

  int n_region_; //!< Number of inteprolation regions
  std::vector<int> breakpoints_; //!< Breakpoints between regions
  std::vector<Interpolation> interpolation_; //!< Interpolation laws
  std::vector<double> energy_; //!< Incident energy in [eV]
  std::vector<CTTable> distribution_; //!< Distributions for each incident energy
};

class CTTableFlat {
public:
  explicit CTTableFlat(const uint8_t* data);

  Interpolation interpolation() const;
  int n_discrete() const;
  gsl::span<const double> e_out() const;
  gsl::span<const double> p() const;
  gsl::span<const double> c() const;
private:
  const uint8_t* data_;
  size_t n_eout_;
};

class ContinuousTabularFlat {
public:
  explicit ContinuousTabularFlat(const uint8_t* data);

  double sample(double E, uint64_t* seed) const;
private:
  gsl::span<const int> breakpoints() const;
  Interpolation interpolation(gsl::index i) const;
  gsl::span<const double> energy() const;
  CTTableFlat distribution(gsl::index i) const;

  const uint8_t* data_;
  size_t n_region_;
  size_t n_energy_;
};

//===============================================================================
//! Evaporation spectrum corresponding to ACE law 9 and ENDF File 5, LF=9.
//===============================================================================

class Evaporation : public EnergyDistribution {
public:
  explicit Evaporation(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;

  size_t nbytes() const { return 12 + theta_.nbytes(); }

  void serialize(DataBuffer& buffer) const;
private:
  Tabulated1D theta_; //!< Incoming energy dependent parameter
  double u_; //!< Restriction energy
};

class EvaporationFlat {
public:
  explicit EvaporationFlat(const uint8_t* data) : data_(data) { }

  double sample(double E, uint64_t* seed) const;

private:
  double u() const;
  Tabulated1DFlat theta() const;

  const uint8_t* data_;
};

//===============================================================================
//! Energy distribution of neutrons emitted from a Maxwell fission spectrum.
//! This corresponds to ACE law 7 and ENDF File 5, LF=7.
//===============================================================================

class MaxwellEnergy : public EnergyDistribution {
public:
  explicit MaxwellEnergy(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;

  size_t nbytes() const { return 12 + theta_.nbytes(); }

  void serialize(DataBuffer& buffer) const;
private:
  Tabulated1D theta_; //!< Incoming energy dependent parameter
  double u_; //!< Restriction energy
};

class MaxwellFlat {
public:
  explicit MaxwellFlat(const uint8_t* data) : data_(data) { }

  double sample(double E, uint64_t* seed) const;

private:
  double u() const;
  Tabulated1DFlat theta() const;

  const uint8_t* data_;
};

//===============================================================================
//! Energy distribution of neutrons emitted from a Watt fission spectrum. This
//! corresponds to ACE law 11 and ENDF File 5, LF=11.
//===============================================================================

class WattEnergy : public EnergyDistribution {
public:
  explicit WattEnergy(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;

  size_t nbytes() const;

  void serialize(DataBuffer& buffer) const;
private:
  Tabulated1D a_; //!< Energy-dependent 'a' parameter
  Tabulated1D b_; //!< Energy-dependent 'b' parameter
  double u_; //!< Restriction energy
};

class WattFlat {
public:
  explicit WattFlat(const uint8_t* data) : data_(data) { }

  double sample(double E, uint64_t* seed) const;

private:
  Tabulated1DFlat a() const;
  Tabulated1DFlat b() const;
  double u() const;

  const uint8_t* data_;
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_ENERGY_H
