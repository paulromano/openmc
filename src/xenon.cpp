#include <algorithm>
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

int32_t mat_filter_idx;
int32_t fission_tally_idx;
int32_t absorb_tally_idx;

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

  /*
  std::cout << "I135 decay constant: " << i135_decay << '\n';
  std::cout << "Xe135 decay constant: " << xe135_decay << '\n';
  for (const auto& kv : i135_yield) {
    std::cout << "I135 yield for " << kv.first << ": " << kv.second << '\n';
  }
  for (const auto& kv : xe135_yield) {
    std::cout << "Xe135 yield for " << kv.first << ": " << kv.second << '\n';
  }
  for (const auto& kv : fission_q) {
    std::cout << "fission Q for " << kv.first << ": " << kv.second << '\n';
  }
  */
}

void get_nuclide_names () {
  for (int i = 1; i <= n_nuclides; ++i) {
    nucname[i] = openmc::nuclide_name(i);
  }
}

void create_tallies(const std::set<int>& nucs, const std::vector<int32_t>& mats) {
  // Determine maximum tally/filter ID used so far
  int32_t max_filter_id = 0;
  int32_t filter_id;
  for (int i = 1; i <= n_filters; ++i) {
    CHECK(openmc_filter_get_id(i, &filter_id));
    max_filter_id = std::max(max_filter_id, filter_id);
  }

  int32_t max_tally_id = 0;
  int32_t tally_id;
  for (int i = 1; i <= n_tallies; ++i) {
    CHECK(openmc_tally_get_id(i, &tally_id));
    max_tally_id = std::max(max_tally_id, tally_id);
  }

  // Create filter for each material with fissionable material
  CHECK(openmc_extend_filters(1, &mat_filter_idx, nullptr));
  CHECK(openmc_filter_set_type(mat_filter_idx, "material"));
  CHECK(openmc_filter_set_id(mat_filter_idx, max_filter_id + 1));
  CHECK(openmc_material_filter_set_bins(mat_filter_idx, mats.size(), mats.data()));

  //////////////////////////////////////////////////////////////////////////////
  // Create tally for fission in each fissionable nuclide
  CHECK(openmc_extend_tallies(1, &fission_tally_idx, nullptr));
  CHECK(openmc_tally_set_type(fission_tally_idx, "generic"));
  CHECK(openmc_tally_set_id(fission_tally_idx, max_tally_id + 1));
  const int32_t indices[] {mat_filter_idx};
  CHECK(openmc_tally_set_filters(fission_tally_idx, 1, indices));

  // Create array of pointers to nuclide C-strings
  size_t num_nucs = nucs.size();
  const char* nuclides[num_nucs];
  int j = 0;
  for (int i : nucs) {
    nuclides[j++] = nucname[i].c_str();
  }
  // Set nuclides
  CHECK(openmc_tally_set_nuclides(fission_tally_idx, num_nucs, nuclides));

  char score_array[][20]{"fission"};
  const char *scores[]{score_array[0]}; // OpenMC expects a const char**, ugh
  CHECK(openmc_tally_set_scores(fission_tally_idx, 1, scores));

  //////////////////////////////////////////////////////////////////////////////
  // Create tally for absorption in Xe135

  CHECK(openmc_extend_tallies(1, &absorb_tally_idx, nullptr));
  CHECK(openmc_tally_set_type(absorb_tally_idx, "generic"));
  CHECK(openmc_tally_set_id(absorb_tally_idx, max_tally_id + 2));
  CHECK(openmc_tally_set_filters(absorb_tally_idx, 1, indices));

  // Load Xe135
  CHECK(openmc_load_nuclide("Xe135"));

  // Set nuclides
  nuclides[0] = "Xe135";
  CHECK(openmc_tally_set_nuclides(absorb_tally_idx, 1, nuclides));

  std::strncpy(score_array[0], "absorption", 20);
  CHECK(openmc_tally_set_scores(absorb_tally_idx, 1, scores));
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

void init() {
  // Get yields/decay constants from depletion chain file
  xenon::get_chain_data();

  // Determine fissionable materials and nuclides
  std::set<int> nucs;
  std::vector<int32_t> mats;
  std::tie(nucs, mats) = xenon::fissionable_materials();

  // Create tallies needed for Xenon feedback
  xenon::create_tallies(nucs, mats);
}

} // namespace xenon

int main(int argc, char* argv[]) {
  // Initialize OpenMC
  CHECK(openmc_init(argc, argv, nullptr));

  // Xe feedback initialization
  xenon::init();

  // Run
  CHECK(openmc_run());

  // Finalize simulation
  CHECK(openmc_finalize());
}
