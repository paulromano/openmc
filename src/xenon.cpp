#include <cmath>
#include <cstring>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <iostream>

#include "openmc.h"
#include "error.h"
#include "nuclide.h"
#include "pugixml.hpp"


void CHECK(int err) {
  if (err) openmc::fatal_error(openmc_err_msg);
}

namespace xenon {

// Global variables

std::unordered_map<std::string, double> i135_yield;
std::unordered_map<std::string, double> xe135_yield;
std::unordered_map<std::string, double> fission_q;
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
    const char* name = nuclide.attribute("name").value();

    //////////////////////////////////////////////////////////////////////////
    // Determine fission Q values
    for (xml_node reaction : nuclide.children("reaction")) {
      const char* rx_type = reaction.attribute("type").value();
      if (std::strcmp(rx_type, "fission") == 0) {
        fission_q[name] = reaction.attribute("Q").as_double();
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Determine decay constants for Xe135 and I135
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
  for (int i = 1; i <= n_nuclides; ++i) {
    nucname[i] = openmc::nuclide_name(i);
  }
}

std::pair<std::set<int>, std::vector<int32_t>> fissionable_materials() {
  std::set<int> fissionable_nuclides;
  std::vector<int32_t> material_indices;

  int* nuclides;
  double* densities;
  int n;
  for (int i = 1; i <= n_materials; ++i) {
    // Get arrays of nuclide indices and densities
    CHECK(openmc_material_get_densities(i, &nuclides, &densities, &n));

    // loop over nuclides
    bool added = false;
    for (int j = 0; j < n; ++j) {
      // get name of j-th nuclide
      std::string name = nucname[nuclides[j]];

      // if nuclide is fissionable, add it to set of fissionable nuclides and
      // add material to list of fuel materials
      if (fission_q.find(name) != fission_q.end()) {
        fissionable_nuclides.insert(nuclides[j]);
        if (!added) {
          material_indices.push_back(i);
          added = true;
        }
      }
    }
  }

  return {fissionable_nuclides, material_indices};
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
  for (const auto& kv : xenon::fission_q) {
    std::cout << "fission Q for " << kv.first << ": " << kv.second << '\n';
  }

  xenon::get_nuclide_names();
  for (const auto& kv : xenon::nucname) {
    std::cout << "Nuclide " << kv.first << ": " << kv.second << '\n';
  }

  std::set<int> nucs;
  std::vector<int32_t> mats;
  std::tie(nucs, mats) = xenon::fissionable_materials();

  std::cout << "Fissionable nuclides (in model): \n";
  for (int i : nucs) {
    std::cout << "  " << xenon::nucname[i] << '\n';
  }

  std::cout << "Fissionable materials:\n";
  for (int32_t idx : mats) {
    std::cout << idx << '\n';
  }

  // Create filter
  int32_t f_idx;
  CHECK(openmc_extend_filters(1, &f_idx, nullptr));
  CHECK(openmc_filter_set_type(f_idx, "material"));
  CHECK(openmc_material_filter_set_bins(f_idx, mats.size(), mats.data()));

  // Create tally for fission
  int32_t t_idx;
  CHECK(openmc_extend_tallies(1, &t_idx, nullptr));
  CHECK(openmc_tally_set_type(t_idx, "generic"));
  const int32_t indices[] {f_idx};
  CHECK(openmc_tally_set_filters(t_idx, 1, indices));

  size_t num_nucs = nucs.size();
  const char* nuclides[num_nucs];
  int j = 0;
  for (int i : nucs) {
    nuclides[j++] = xenon::nucname[i].c_str();
  }
  CHECK(openmc_tally_set_nuclides(t_idx, num_nucs, nuclides));

  char score_array[][20]{"fission"};
  const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
  CHECK(openmc_tally_set_scores(t_idx, 1, scores));

  // Run
  // CHECK(openmc_run());

  // Finalize simulation
  CHECK(openmc_finalize());
}
