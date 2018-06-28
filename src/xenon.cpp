#include <cmath>
#include <cstring>
#include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <iostream>

#include "openmc.h"
#include "error.h"
#include "pugixml.hpp"


void CHECK(int err) {
  if (err) openmc::fatal_error(openmc_err_msg);
}

namespace xenon {

// Global variables

std::unordered_map<std::string, double> i135_yield;
std::unordered_map<std::string, double> xe135_yield;
std::unordered_map<int, std::string> nucname;
double i135_decay;
double xe135_decay;

std::vector<std::string> split(const char* s) {
  std::istringstream iss {s};
  return {std::istream_iterator<std::string>{iss},
          std::istream_iterator<std::string>{}};
}

void get_chain_data() {
  using namespace pugi;

  const char* chain_file = std::getenv("OPENMC_DEPLETE_CHAIN");
  if (!chain_file)
    openmc::fatal_error("Must have a depletion chain file.");

  xml_document doc;
  xml_parse_result result = doc.load_file(chain_file);
  xml_node root = doc.document_element();

  for (xml_node nuclide : root.children("nuclide")) {
    //////////////////////////////////////////////////////////////////////////
    // Determine decay constants for Xe135 and I135
    const char* name = nuclide.attribute("name").value();
    if (std::strcmp(name, "Xe135") == 0) {
      double half_life = nuclide.attribute("half_life").as_double();
      xe135_decay = std::log(2.0) / half_life;
    } else if (std::strcmp(name, "I135") == 0) {
      double half_life = nuclide.attribute("half_life").as_double();
      i135_decay = std::log(2.0) / half_life;
    }

    //////////////////////////////////////////////////////////////////////////
    // Determine Xe135 and I135 yields for each fissionable nuclide
    if (xml_node nfy = nuclide.child("neutron_fission_yields")) {
      // Get the first <fission_yields> node. This usually corresponds to the
      // yields at thermal energies but NOT ALWAYS. So, it is an approximation
      // for now.
      xml_node yields = nfy.child("fission_yields");

      auto products = split(yields.child("products").text().get());
      auto data = split(yields.child("data").text().get());
      // determine number of products
      for (int i = 0; i < products.size(); ++i) {
        if (products[i] == "I135") {
          i135_yield[name] = std::stod(data[i]);
        } else if (products[i] == "Xe135") {
          xe135_yield[name] = std::stod(data[i]);
        }
      } // <products>
    } // <neutron_fission_yields>
  } // <nuclide>
}

void get_nuclide_names () {
  char* name;
  for (int i = 1; i <= n_nuclides; ++i) {
    CHECK(openmc_nuclide_name(i, &name));
    for (size_t j = 0; j < 20; ++j) {
      if (name[j] == ' ') {
        nucname[i] = std::string{name, j};
        break;
      }
    }
  }
}

std::tuple<std::vector<int32_t>, std::vector<int32_t>> fissionable_materials() {
  std::vector<int32_t> material_indices;

  int* nuclides;
  double* densities;
  int n;
  for (int i = 1; i <= n_materials; ++i) {
    // Get arrays of nuclide indices and densities
    CHECK(openmc_material_get_densities(i, &nuclides, &densities, &n));

    // If material is fissionable add it to vector
    char* name;
    for (int j = 0; j < n; ++j) {
    }
  }
}

} // namespace xenon

int main(int argc, char* argv[]) {
  // Initialize OpenMC
  CHECK(openmc_init(argc, argv, nullptr));

  xenon::get_chain_data();
  std::cout << "I135 decay constant: " << xenon::i135_decay << '\n';
  std::cout << "Xe135 decay constant: " << xenon::xe135_decay << '\n';
  for (const auto& kv : xenon::i135_yield) {
    std::cout << "I135 yield for " << kv.first << ": " << kv.second << '\n';
  }
  for (const auto& kv : xenon::xe135_yield) {
    std::cout << "Xe135 yield for " << kv.first << ": " << kv.second << '\n';
  }

  xenon::get_nuclide_names();
  for (const auto& kv : xenon::nucname) {
    std::cout << "Nuclide " << kv.first << ": " << kv.second << '\n';
  }

  // Run
  // CHECK(openmc_run());

  // Finalize simulation
  CHECK(openmc_finalize());
}
