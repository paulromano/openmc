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
std::unordered_map<int32_t, double> volume;
double i135_decay;
double xe135_decay;
int32_t mat_filter_idx;
int32_t fission_tally_idx;
int32_t absorb_tally_idx;


double actual_power = 20.0;
int batch_size = 5;
int restart_batches = 4;

// Helper function

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
  CHECK(openmc_tally_set_active(fission_tally_idx, true));
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
  CHECK(openmc_tally_set_active(absorb_tally_idx, true));
  CHECK(openmc_tally_set_filters(absorb_tally_idx, 1, indices));

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
  for (int32_t i = 1; i <= n_materials; ++i) {
    // Get arrays of nuclide indices and densities
    CHECK(openmc_material_get_densities(i, &nuclides, &densities, &n));

    // loop over nuclides
    bool added = false;
    bool has_xe135 = false;
    bool has_i135 = false;
    for (int j = 0; j < n; ++j) {
      // get name of j-th nuclide
      std::string name = openmc::nuclide_name(nuclides[j]);

      // Check if Xe135/I135 is present
      if (name == "Xe135") has_xe135 = true;
      if (name == "I135") has_i135 = true;

      // if nuclide is fissionable, add it to set of fissionable nuclides and
      // add material to list of fuel materials
      if (fission_q.find(name) != fission_q.end()) {
        fissionable_nuclides.insert(nuclides[j]);
        if (!added) {
          // Determine volume -- note that material_get_volume will fail if
          // volume is not set
          double vol;
          CHECK(openmc_material_get_volume(i, &vol));
          volume[i] = vol;

          material_indices.push_back(i);
          added = true;
        }
      }
    } // nuclides

    if (!has_xe135) CHECK(openmc_material_add_nuclide(i, "Xe135", 1e-14));
    if (!has_i135) CHECK(openmc_material_add_nuclide(i, "I135", 1e-14));
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

  // Save nuclide names for quick access
  xenon::get_nuclide_names();

  // Create tallies needed for Xenon feedback
  xenon::create_tallies(nucs, mats);
}

void update() {
  // Get material bins
  int32_t* mats;
  int32_t n_mats;
  CHECK(openmc_material_filter_get_bins(mat_filter_idx, &mats, &n_mats));

  // Get nuclides
  int* nuclides;
  int n_nucs;
  CHECK(openmc_tally_get_nuclides(fission_tally_idx, &nuclides, &n_nucs));
  double* fission;
  int shape[3];
  CHECK(openmc_tally_results(fission_tally_idx, &fission, shape));
  int32_t m;
  CHECK(openmc_tally_get_n_realizations(fission_tally_idx, &m));

  double* absorb;
  CHECK(openmc_tally_results(absorb_tally_idx, &absorb, shape));

  // Create vectors for I135/Xe135 production, Xe135 absorption
  std::vector<double> i135_prod(n_mats);
  std::vector<double> xe135_prod(n_mats);
  std::vector<double> xe135_abs(n_mats);

  double power = 0.0;

  for (int i = 0; i < n_mats; ++i) {
    // Determine I135 and Xe135 production rates. These have to be summed over
    // individual nuclides because each nuclide has a different fission product
    // yield.
    for (int j = 0; j < n_nucs; ++j) {
      double fission_rate = fission[3*n_nucs*i + 3*j + 1] / m;
      std::string nuc = nucname[nuclides[j]];
      i135_prod[i] += i135_yield[nuc] * fission_rate;
      xe135_prod[i] += xe135_yield[nuc] * fission_rate;
      power += fission_q[nuc] * fission_rate;
    }

    // Determine Xe135 absorption rate
    xe135_abs[i] = absorb[3*i + 1] / m;
  }

  // Now that the raw production/absorption rates and total power have been
  // determined, we need to normalize the values by the total power and update
  // densities

  for (int i = 0; i < n_mats; ++i) {
    // Normalize production rates by specified power and volume
    double V = volume[mats[i]];
    double normalization = actual_power / (power * V);
    i135_prod[i] *= normalization;
    xe135_prod[i] *= normalization;
    xe135_abs[i] *= normalization;

    // Get array of densities in current material
    int* nucs_in_mat;
    double* densities;
    int n_nucs_in_mat;
    CHECK(openmc_material_get_densities(mats[i], &nucs_in_mat, &densities,
                                        &n_nucs_in_mat));

    // Divide Xe135 absorption rate by current Xe135 density
    int idx_i135;
    int idx_xe135;
    for (int j = 0; j < n_nucs_in_mat; ++j) {
      if (nucname[nucs_in_mat[j]] == "Xe135") {
        xe135_abs[i] /= densities[j];
        idx_xe135 = j;
      } else if (nucname[nucs_in_mat[j]] == "I135") {
        idx_i135 = j;
      }
    }

    // solve equation (9) and (10) for equilibrium concentrations
    double i135_eq = i135_prod[i] / i135_decay;
    double xe135_eq = (i135_prod[i] + xe135_prod[i])/(xe135_decay + xe135_abs[i]);

    // Create array of pointers to nuclide C-strings
    const char* new_nucnames[n_nucs_in_mat];
    for (int j = 0; j < n_nucs_in_mat; ++j) {
      new_nucnames[j] = nucname[nucs_in_mat[j]].c_str();
    }

    // Copy existing densities and update I135 and Xe135
    double new_densities[n_nucs_in_mat];
    std::copy(densities, densities + n_nucs_in_mat, new_densities);
    new_densities[idx_i135] = i135_eq;
    new_densities[idx_xe135] = xe135_eq;

    // Set densities for material
    CHECK(openmc_material_set_densities(mats[i], n_nucs_in_mat, new_nucnames,
                                        new_densities));
  }
}

int run() {
  openmc_simulation_init();

  int err = 0;
  int status = 0;
  int generation = 1;
  int batch = 1;
  while (true) {
    err = openmc_next_batch(&status);
    if (status != 0 || err != 0) break;

    if (batch_size == generation) {
      std::cout << "Updating Xenon densities\n";
      xenon::update();
      if (batch <= restart_batches) {
        std::cout << "Resetting Xenon tallies\n";
        CHECK(openmc_tally_reset(fission_tally_idx));
        CHECK(openmc_tally_reset(absorb_tally_idx));
      }
      ++batch;
      generation = 1;
    } else {
      ++generation;
    }
  }

  openmc_simulation_finalize();
  return err;
}

} // namespace xenon

int main(int argc, char* argv[]) {
  // Initialize OpenMC
  CHECK(openmc_init(argc, argv, nullptr));

  // Xe feedback initialization
  xenon::init();

  // Run
  CHECK(xenon::run());

  // Finalize simulation
  CHECK(openmc_finalize());
}
