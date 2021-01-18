#ifndef OPENMC_WMP_H
#define OPENMC_WMP_H

#include "hdf5.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "xtensor/xtensor.hpp"

#include <array>
#include <complex>
#include <memory>
#include <string>
#include <tuple>

namespace openmc {

//========================================================================
// Constants
//========================================================================

// Constants that determine which value to access
constexpr int MP_EA {0}; // Pole
constexpr int MP_RS {1}; // Residue scattering
constexpr int MP_RA {2}; // Residue absorption
constexpr int MP_RF {3}; // Residue fission

// Polynomial fit indices
constexpr int FIT_S {0}; // Scattering
constexpr int FIT_A {1}; // Absorption
constexpr int FIT_F {2}; // Fission

// Multipole HDF5 file version
constexpr std::array<int, 2> WMP_VERSION {1, 1};

//========================================================================
// Windowed multipole data
//========================================================================

class WindowedMultipole {
public:
  WindowedMultipole(const std::string& name, double const& E_min,
    double const& E_max, double const& sqrt_awr, double const& inv_spacing)
    : name_(name), E_min_(E_min), E_max_(E_max), sqrt_awr_(sqrt_awr),
      inv_spacing_(inv_spacing)
  {}
  virtual ~WindowedMultipole() = default;

  //! \brief Evaluate the windowed multipole equations for cross sections in the
  //! resolved resonance regions
  //!
  //! \param E Incident neutron energy in [eV]
  //! \param sqrtkT Square root of temperature times Boltzmann constant
  //! \return Tuple of elastic scattering, absorption, and fission cross sections in [b]
  std::tuple<double, double, double> virtual evaluate(
    double E, double sqrtkT) const = 0;

  //! \brief Evaluates the windowed multipole equations for the derivative of
  //! cross sections in the resolved resonance regions with respect to
  //! temperature.
  //!
  //! \param E Incident neutron energy in [eV]
  //! \param sqrtkT Square root of temperature times Boltzmann constant
  //! \return Tuple of derivatives of elastic scattering, absorption, and
  //!         fission cross sections in [b/K]
  std::tuple<double, double, double> virtual evaluate_deriv(
    double E, double sqrtkT) const = 0;

  const std::string name_;
  const double E_min_;       //!< Minimum energy in [eV]
  const double E_max_;       //!< Maximum energy in [eV]
  const double sqrt_awr_;    //!< Square root of atomic weight ratio
  const double inv_spacing_; //!< 1 / spacing in sqrt(E) space
};

template<unsigned FitOrder, bool Fissionable>
class WindowedMultipoleImplementation : public WindowedMultipole {
public:
  // Types
  struct WindowInfo {
    int index_start;   // Index of starting pole
    int index_end;     // Index of ending pole
    bool broaden_poly; // Whether to broaden polynomial curvefit
  };
  static constexpr unsigned POLYNOMIAL_DATA_LENGTH = 2 + Fissionable;
  static constexpr unsigned POLE_DATA_LENGTH = 3 + Fissionable;
  static constexpr unsigned NUM_POLY_COEFF = FitOrder + 1;
  typedef std::array<std::array<double, POLYNOMIAL_DATA_LENGTH>, NUM_POLY_COEFF>
    PolynomialCoefficients;
  typedef std::array<std::complex<double>, POLE_DATA_LENGTH> PoleAndResidues;

  WindowedMultipoleImplementation(std::string const& name, double const& E_min,
    double const& E_max, double const& sqrt_awr, double const& inv_spacing,
    xt::xtensor<int, 2> const& windows,
    xt::xtensor<bool, 1> const& broaden_poly,
    xt::xtensor<std::complex<double>, 2> const& data,
    xt::xtensor<double, 3> const& curve_fit)
    : WindowedMultipole(name, E_min, E_max, sqrt_awr, inv_spacing),
      window_info_(windows.shape()[0]), curve_fit_(windows.shape()[0]),
      data_(data.shape()[0])
  {
    const unsigned n_windows = window_info_.size();

    // If we looked specifically at the column or row-majorness of the
    // xtensors here, we could just use std::copy. But direct indexing
    // guarantees we're copying the right order.

    // Copy window info
    for (unsigned w = 0; w < n_windows; ++w) {
      WindowInfo& window = window_info_[w];
      window.index_start = windows(w, 0);
      window.index_end = windows(w, 1);
      window.broaden_poly = broaden_poly[w];
    }

    // Copy polynomial curve fits
    for (unsigned w = 0; w < n_windows; ++w) {
#pragma unroll
      for (unsigned coeff = 0; coeff < NUM_POLY_COEFF; ++coeff) {
#pragma unroll
        for (unsigned reaction = 0; reaction < POLYNOMIAL_DATA_LENGTH;
             ++reaction) {
          curve_fit_[w][coeff][reaction] = curve_fit(w, coeff, reaction);
        }
      }
    }

    // Copy poles and residues
    for (unsigned pole = 0; pole < data.shape()[0]; ++pole) {
#pragma unroll
      for (unsigned pole_data_i = 0; pole_data_i < POLE_DATA_LENGTH;
           ++pole_data_i) {
        data_[pole][pole_data_i] = data(pole, pole_data_i);
      }
    }
  }

  std::tuple<double, double, double> virtual evaluate(
    double E, double sqrtkT) const override;
  std::tuple<double, double, double> virtual evaluate_deriv(
    double E, double sqrtkT) const override;

  std::vector<WindowInfo> window_info_; // Information about a window
  std::vector<PolynomialCoefficients> curve_fit_;
  std::vector<PoleAndResidues> data_;
};

//========================================================================
// Non-member functions
//========================================================================

//! Check to make sure WMP library data version matches
//!
//! \param[in] file  HDF5 file object
void check_wmp_version(hid_t file);

//! \brief Checks for the existence of a multipole library in the directory and
//! loads it
//!
//! \param[in] i_nuclide  Index in global nuclides array
void read_multipole_data(int i_nuclide);

std::unique_ptr<WindowedMultipole> make_multipole_implementation(hid_t group);

//==============================================================================
//! Doppler broadens the windowed multipole curvefit.
//!
//! The curvefit is a polynomial of the form a/E + b/sqrt(E) + c + d sqrt(E)...
//! This version is kept for C API compatibility. The templated version does
//! unrolling.
//!
//! \param E       The energy to evaluate the broadening at
//! \param dopp    sqrt(atomic weight ratio / kT) with kT given in eV
//! \param n       The number of components to the polynomial
//! \param factors The output leading coefficient
//==============================================================================

extern "C" void broaden_wmp_polynomials(double E, double dopp, int n, double factors[]);

template<unsigned FitOrder>
inline void broaden_wmp_polynomials(
  double const& E, double const& dopp, double* factors)
{
  double sqrtE = std::sqrt(E);
  double beta = sqrtE * dopp;
  double half_inv_dopp2 = 0.5 / (dopp * dopp);
  double quarter_inv_dopp4 = half_inv_dopp2 * half_inv_dopp2;

  double erf_beta;    // error function of beta
  double exp_m_beta2; // exp(-beta**2)
  if (beta > 6.0) {
    // Save time, ERF(6) is 1 to machine precision.
    // beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
    erf_beta = 1.;
    exp_m_beta2 = 0.;
  } else {
    erf_beta = std::erf(beta);
    exp_m_beta2 = std::exp(-beta * beta);
  }

  // Assume that, for sure, we'll use a second order (1/E, 1/V, const)
  // fit, and no less.

  factors[0] = erf_beta / E;
  factors[1] = 1. / sqrtE;
  factors[2] =
    factors[0] * (half_inv_dopp2 + E) + exp_m_beta2 / (beta * SQRT_PI);

  // Perform recursive broadening of high order components
  constexpr int n = FitOrder + 1;
  constexpr int n_higher_order = n - 3;

  if (FitOrder > 2)
    factors[3] = factors[1] * (E + 3. * half_inv_dopp2);
#pragma unroll
  for (int i = 1; i < n_higher_order; i++) {
    double ip1_dbl = i + 1;
    factors[i + 3] =
      -factors[i - 1] * (ip1_dbl - 1.) * ip1_dbl * quarter_inv_dopp4 +
      factors[i + 1] * (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
  }
}

template<unsigned FitOrder, bool Fissionable>
std::tuple<double, double, double>
WindowedMultipoleImplementation<FitOrder, Fissionable>::evaluate(
  double E, double sqrtkT) const
{
  using namespace std::complex_literals;

  // Define some frequently used variables.
  const double sqrtE = std::sqrt(E);
  const double invE = 1.0 / E;

  // Locate window containing energy
  const int i_window = std::min(window_info_.size() - 1,
    static_cast<size_t>((sqrtE - std::sqrt(E_min_)) * inv_spacing_));
  const auto& window {window_info_[i_window]};
  const int& startw = window.index_start;
  const int& endw = window.index_end;

  // Initialize the ouptut cross sections
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && window.broaden_poly) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    std::array<double, NUM_POLY_COEFF> broadened_polynomials;
    broaden_wmp_polynomials<FitOrder>(E, dopp, broadened_polynomials.data());
#pragma unroll
    for (int i_poly = 0; i_poly < NUM_POLY_COEFF; ++i_poly) {
      sig_s +=
        curve_fit_[i_window][i_poly][FIT_S] * broadened_polynomials[i_poly];
      sig_a +=
        curve_fit_[i_window][i_poly][FIT_A] * broadened_polynomials[i_poly];
      if (Fissionable) {
        sig_f +=
          curve_fit_[i_window][i_poly][FIT_F] * broadened_polynomials[i_poly];
      }
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
#pragma unroll
    for (int i_poly = 0; i_poly < NUM_POLY_COEFF; ++i_poly) {
      sig_s += curve_fit_[i_window][i_poly][FIT_S] * temp;
      sig_a += curve_fit_[i_window][i_poly][FIT_A] * temp;
      if (Fissionable) {
        sig_f += curve_fit_[i_window][i_poly][FIT_F] * temp;
      }
      temp *= sqrtE;
    }
  }

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = startw; i_pole <= endw; ++i_pole) {
      std::complex<double> psi_chi = -1.0i / (data_[i_pole][MP_EA] - sqrtE);
      std::complex<double> c_temp = psi_chi * invE;
      sig_s += (data_[i_pole][MP_RS] * c_temp).real();
      sig_a += (data_[i_pole][MP_RA] * c_temp).real();
      if (Fissionable) {
        sig_f += (data_[i_pole][MP_RF] * c_temp).real();
      }
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    for (int i_pole = startw; i_pole <= endw; ++i_pole) {
      std::complex<double> z = (sqrtE - data_[i_pole][MP_EA]) * dopp;
      std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
      sig_s += (data_[i_pole][MP_RS] * w_val).real();
      sig_a += (data_[i_pole][MP_RA] * w_val).real();
      if (Fissionable) {
        sig_f += (data_[i_pole][MP_RF] * w_val).real();
      }
    }
  }
  return std::make_tuple(sig_s, sig_a, sig_f);
}

template<unsigned FitOrder, bool Fissionable>
std::tuple<double, double, double>
WindowedMultipoleImplementation<FitOrder, Fissionable>::evaluate_deriv(
  double E, double sqrtkT) const
{
  // Define some frequently used variables.
  const double sqrtE = std::sqrt(E);
  const double invE = 1.0 / E;
  double T = sqrtkT * sqrtkT / K_BOLTZMANN;

  if (sqrtkT == 0.0) {
    fatal_error("Windowed multipole temperature derivatives are not implemented"
                " for 0 Kelvin cross sections.");
  }

  // Locate us
  int i_window = (sqrtE - std::sqrt(E_min_)) * inv_spacing_;
  const auto& window {window_info_[i_window]};
  int startw = window.index_start;
  int endw = window.index_end;

  // Initialize the ouptut cross sections.
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // TODO Polynomials: Some of the curvefit polynomials Doppler broaden so
  // rigorously we should be computing the derivative of those.  But in
  // practice, those derivatives are only large at very low energy and they
  // have no effect on reactor calculations.

  // ==========================================================================
  // Add the contribution from the poles in this window.

  double dopp = sqrt_awr_ / sqrtkT;
  for (int i_pole = startw; i_pole <= endw; ++i_pole) {
    std::complex<double> z = (sqrtE - data_[i_pole][MP_EA]) * dopp;
    std::complex<double> w_val = -invE * SQRT_PI * 0.5 * w_derivative(z, 2);
    sig_s += (data_[i_pole][MP_RS] * w_val).real();
    sig_a += (data_[i_pole][MP_RA] * w_val).real();
    if (Fissionable) {
      sig_f += (data_[i_pole][MP_RF] * w_val).real();
    }
  }
  double norm = -0.5 * sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
  sig_s *= norm;
  sig_a *= norm;
  sig_f *= norm;

  return std::make_tuple(sig_s, sig_a, sig_f);
}

} // namespace openmc

#endif // OPENMC_WMP_H
