#include "openmc/wmp.h"

#include "openmc/constants.h"
#include "openmc/cross_sections.h"
#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"

#include <fmt/core.h>

#include <cmath>
#include <algorithm> // for min

namespace openmc {

void check_wmp_version(hid_t file)
{
  if (attribute_exists(file, "version")) {
    std::array<int, 2> version;
    read_attribute(file, "version", version);
    if (version[0] != WMP_VERSION[0]) {
      fatal_error(fmt::format(
        "WMP data format uses version {}.{} whereas your installation of "
        "OpenMC expects version {}.x data.",
        version[0], version[1], WMP_VERSION[0]));
    }
  } else {
    fatal_error(fmt::format("WMP data does not indicate a version. Your "
      "installation of OpenMC expects version {}x data.", WMP_VERSION[0]));
  }
}

void read_multipole_data(int i_nuclide)
{
  // Look for WMP data in cross_sections.xml
  const auto& nuc {data::nuclides[i_nuclide]};
  auto it = data::library_map.find({Library::Type::wmp, nuc->name_});

  // If no WMP library for this nuclide, just return
  if (it == data::library_map.end()) return;

  // Check if WMP library exists
  int idx = it->second;
  std::string& filename = data::libraries[idx].path_;

  // Display message
  write_message(6, "Reading {} WMP data from {}", nuc->name_, filename);

  // Open file and make sure version is sufficient
  hid_t file = file_open(filename, 'r');
  check_wmp_version(file);

  // Read nuclide data from HDF5
  hid_t group = open_group(file, nuc->name_.c_str());
  nuc->multipole_ = make_multipole_implementation(group);
  close_group(group);
  file_close(file);
}

void broaden_wmp_polynomials(double E, double dopp, int n, double factors[])
{
  // Factors is already pre-allocated
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
  factors[2] = factors[0] * (half_inv_dopp2 + E) + exp_m_beta2 /
       (beta * SQRT_PI);

  // Perform recursive broadening of high order components
  for (int i = 0; i < n - 3; i++) {
    double ip1_dbl = i + 1;
    if (i != 0) {
      factors[i + 3] = -factors[i - 1] * (ip1_dbl - 1.) * ip1_dbl *
           quarter_inv_dopp4 + factors[i + 1] *
           (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
    } else {
      // Although it's mathematically identical, factors[0] will contain
      // nothing, and we don't want to have to worry about memory.
      factors[i + 3] = factors[i + 1] *
           (E + (1. + 2. * ip1_dbl) * half_inv_dopp2);
    }
  }
}

std::unique_ptr<WindowedMultipole> make_multipole_implementation(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  std::string name = object_name(group).substr(1);

  // Read scalar values.
  double inv_spacing, sqrt_awr, E_min, E_max;
  read_dataset(group, "spacing", inv_spacing);
  inv_spacing = 1.0 / inv_spacing;
  read_dataset(group, "sqrtAWR", sqrt_awr);
  read_dataset(group, "E_min", E_min);
  read_dataset(group, "E_max", E_max);

  // Read the "data" array.  Use its shape to figure out the number of poles
  // and residue types in this data.
  xt::xtensor<std::complex<double>, 2> data;
  read_dataset(group, "data", data);
  int n_residues = data.shape()[1] - 1;

  // Check to see if this data includes fission residues.
  bool fissionable = (n_residues == 3);

  // Read the "windows" array and use its shape to figure out the number of
  // windows.
  xt::xtensor<int, 2> windows;
  read_dataset(group, "windows", windows);
  int n_windows = windows.shape()[0];
  windows -= 1; // Adjust to 0-based indices

  // Read the "broaden_poly" arrays.
  xt::xtensor<bool, 1> broaden_poly;
  read_dataset(group, "broaden_poly", broaden_poly);
  if (n_windows != broaden_poly.shape()[0]) {
    fatal_error("broaden_poly array shape is not consistent with the windows "
                "array shape in WMP library for " +
                name + ".");
  }

  // Read the "curvefit" array.
  xt::xtensor<double, 3> curvefit;
  read_dataset(group, "curvefit", curvefit);
  if (n_windows != curvefit.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
                "array shape in WMP library for " +
                name + ".");
  }
  unsigned fit_order = curvefit.shape()[1] - 1;

  // Check the code is compiling to work with sufficiently high fit order
  constexpr unsigned MAX_POLY_COEFFICIENTS = 11;
  if (fit_order + 1 > MAX_POLY_COEFFICIENTS) {
    fatal_error(fmt::format(
      "Need to compile with WindowedMultipole::MAX_POLY_COEFFICIENTS = {} and "
      "update switch table in make_multipole_implementation",
      fit_order + 1));
  }

  if (fissionable) {
    constexpr bool Fissionable = true;
    switch (fit_order) {
    case 0: {
      auto ptr =
        new WindowedMultipoleImplementation<0, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 1: {
      auto ptr =
        new WindowedMultipoleImplementation<1, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 2: {
      auto ptr =
        new WindowedMultipoleImplementation<2, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 3: {
      auto ptr =
        new WindowedMultipoleImplementation<3, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 4: {
      auto ptr =
        new WindowedMultipoleImplementation<4, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 5: {
      auto ptr =
        new WindowedMultipoleImplementation<5, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 6: {
      auto ptr =
        new WindowedMultipoleImplementation<6, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 7: {
      auto ptr =
        new WindowedMultipoleImplementation<7, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 8: {
      auto ptr =
        new WindowedMultipoleImplementation<8, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 9: {
      auto ptr =
        new WindowedMultipoleImplementation<9, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 10: {
      auto ptr =
        new WindowedMultipoleImplementation<10, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 11: {
      auto ptr =
        new WindowedMultipoleImplementation<11, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    default: {
      fatal_error("Fit order too high...");
      std::unique_ptr<WindowedMultipole> ptr(nullptr);
      return ptr;
    }
    }
  } else {
    constexpr bool Fissionable = false;
    switch (fit_order) {
    case 0: {
      auto ptr =
        new WindowedMultipoleImplementation<0, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 1: {
      auto ptr =
        new WindowedMultipoleImplementation<1, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 2: {
      auto ptr =
        new WindowedMultipoleImplementation<2, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 3: {
      auto ptr =
        new WindowedMultipoleImplementation<3, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 4: {
      auto ptr =
        new WindowedMultipoleImplementation<4, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 5: {
      auto ptr =
        new WindowedMultipoleImplementation<5, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 6: {
      auto ptr =
        new WindowedMultipoleImplementation<6, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 7: {
      auto ptr =
        new WindowedMultipoleImplementation<7, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 8: {
      auto ptr =
        new WindowedMultipoleImplementation<8, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 9: {
      auto ptr =
        new WindowedMultipoleImplementation<9, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 10: {
      auto ptr =
        new WindowedMultipoleImplementation<10, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    case 11: {
      auto ptr =
        new WindowedMultipoleImplementation<11, Fissionable>(name, E_min, E_max,
          sqrt_awr, inv_spacing, windows, broaden_poly, data, curvefit);
      return std::unique_ptr<WindowedMultipole>(ptr);
    }
    default: {
      fatal_error("Fit order too high...");
      std::unique_ptr<WindowedMultipole> ptr(nullptr);
      return ptr;
    }
    }
  }
}

} // namespace openmc
