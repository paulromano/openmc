#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/string_utils.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/tally.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"

#include "pugixml.hpp"

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


void CHECK(int err) {
  if (err) openmc::fatal_error(openmc_err_msg);
}

namespace openmc{
namespace xenon {

// Global variables

struct NuclideProperties {
  double yield_te135;
  double yield_i135;
  double yield_xe135;
  double yield_xe135m;
  double fission_q;
};

std::unordered_map<std::string, NuclideProperties> properties;
std::vector<std::string> nucname;
std::unordered_map<int32_t, double> volume;
double i135_decay;
double xe135_decay;
double b_m;
double b_g;
double b_it;
MaterialFilter* mat_filter;
Tally* fission_tally;
Tally* absorb_tally;


double actual_power = 104.5*PI*0.4096*0.4096;
int batch_size = 5;
int restart_batches = 4;


void get_chain_data() {
  using namespace pugi;

  // TODO: get chain from cross_sections.xml
  const char* chain_file = std::getenv("OPENMC_DEPLETE_CHAIN");
  if (!chain_file)
    throw std::runtime_error{"Must have a depletion chain file."};

  xml_document doc;
  xml_parse_result result = doc.load_file(chain_file);
  xml_node root = doc.document_element();

  for (xml_node nuclide : root.children("nuclide")) {
    std::string name = nuclide.attribute("name").value();

    //////////////////////////////////////////////////////////////////////////
    // Determine fission Q values
    for (xml_node reaction : nuclide.children("reaction")) {
      std::string rx_type = reaction.attribute("type").value();
      if (rx_type == "fission") {
        properties[name].fission_q = openmc::JOULE_PER_EV * reaction.attribute("Q").as_double();
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Determine branching ratios for I135 to Xe135 and Xe135m and Xe135m to
    // Xe135
    for (xml_node decay : nuclide.children("decay")) {
      std::string target = decay.attribute("target").value();
      if (name == "I135") {
        if (target == "Xe135") {
          b_g = decay.attribute("branching_ratio").as_double();
        } else if (target == "Xe135_m1") {
          b_m = decay.attribute("branching_ratio").as_double();
        }
      } else if (name == "Xe135_m1" && target == "Xe135") {
        b_it = decay.attribute("branching_ratio").as_double();
      }
    }

    //////////////////////////////////////////////////////////////////////////
    // Determine decay constants for Xe135 and I135
    if (name == "Xe135") {
      double half_life = nuclide.attribute("half_life").as_double();
      xe135_decay = std::log(2.0) / half_life;
    } else if (name == "I135") {
      double half_life = nuclide.attribute("half_life").as_double();
      i135_decay = std::log(2.0) / half_life;
    }

    //////////////////////////////////////////////////////////////////////////
    // Determine Xe135, Xe135m, Te135, and I135 yields for each fissionable
    // nuclide
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
          properties[name].yield_i135 = std::stod(data[i]);
        } else if (products[i] == "Te135") {
          properties[name].yield_te135 = std::stod(data[i]);
        } else if (products[i] == "Xe135") {
          properties[name].yield_xe135 = std::stod(data[i]);
        } else if (products[i] == "Xe135_m1") {
          properties[name].yield_xe135m = std::stod(data[i]);
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
  for (int i = 0; i < data::nuclides.size(); ++i) {
    nucname[i] = data::nuclides[i]->name_;
  }
}

void create_tallies(const std::vector<std::string>& nucs, const std::vector<int32_t>& mats) {
  // Create filter for each material with fissionable material
  mat_filter = dynamic_cast<MaterialFilter*>(Filter::create("material"));
  mat_filter->set_materials(mats);

  //////////////////////////////////////////////////////////////////////////////
  // Create tally for fission in each fissionable nuclide
  fission_tally = Tally::create();
  fission_tally->set_active(true);
  std::vector<Filter*> filters = {mat_filter};
  fission_tally->set_filters(filters);
  fission_tally->set_nuclides(nucs);
  fission_tally->set_scores({"fission"});

  //////////////////////////////////////////////////////////////////////////////
  // Create tally for absorption in Xe135

  absorb_tally = Tally::create();
  absorb_tally->set_active(true);
  absorb_tally->set_filters(filters);
  absorb_tally->set_nuclides({"Xe135"});
  absorb_tally->set_scores({"absorption"});
}

std::pair<std::vector<std::string>, std::vector<int32_t>>
fissionable_materials()
{
  std::vector<std::string> fissionable_nuclides;
  std::vector<int32_t> material_indices;

  int* nuclides;
  double* densities;
  int n;
  for (int32_t i = 0; i < model::materials.size(); ++i) {
    const auto& mat {model::materials[i]};

    // loop over nuclides
    bool added = false;
    bool has_xe135 = false;
    bool has_i135 = false;
    for (int j = 0; j < mat->nuclide_.size(); ++j) {
      // get name of j-th nuclide
      int i_nuc = mat->nuclide_[j];
      auto name = data::nuclides[i_nuc]->name_;

      // Check if Xe135/I135 is present
      if (name == "Xe135") has_xe135 = true;
      if (name == "I135") has_i135 = true;

      // if nuclide is fissionable, add it to set of fissionable nuclides and
      // add material to list of fuel materials
      if (properties.find(name) != properties.end()) {
        fissionable_nuclides.push_back(name);
        if (!added) {
          volume[i] = mat->volume();
          material_indices.push_back(i);
          added = true;
        }
      }
    } // nuclides

    if (added && !has_xe135) mat->add_nuclide("Xe135", 1e-14);
    if (added && !has_i135) mat->add_nuclide("I135", 1e-14);
  }

  return {fissionable_nuclides, material_indices};
}

void init() {
  // Get yields/decay constants from depletion chain file
  xenon::get_chain_data();

  // Determine fissionable materials and nuclides
  std::vector<std::string> nucs;
  std::vector<int32_t> mats;
  std::tie(nucs, mats) = xenon::fissionable_materials();

  // Save nuclide names for quick access
  xenon::get_nuclide_names();

  // Create tallies needed for Xenon feedback
  xenon::create_tallies(nucs, mats);
}

void update() {
  // Get material bins
  const auto& mats {mat_filter->materials()};
  auto n_mats = mats.size();

  // Create vectors for Te135/I135/Xe135/Xe135m production, Xe135 absorption
  std::vector<double> te135_prod(n_mats);
  std::vector<double> i135_prod(n_mats);
  std::vector<double> xe135_prod(n_mats);
  std::vector<double> xe135m_prod(n_mats);
  std::vector<double> xe135_abs(n_mats);

  double power = 0.0;

  for (gsl::index i = 0; i < n_mats; ++i) {
    // Determine Te135, I135 and Xe135 production rates. These have to be
    // summed over individual nuclides because each nuclide has a different
    // fission product yield.
    for (int j = 0; j < fission_tally->nuclides_.size(); ++j) {
      double fission_rate = fission_tally->results_(i, j, 1) / fission_tally->n_realizations_;
      std::string nuc = nucname[fission_tally->nuclides_[j]];
      const auto& props = properties[nuc];
      i135_prod[i] += props.yield_i135 * fission_rate;
      te135_prod[i] += props.yield_te135 * fission_rate;
      xe135_prod[i] += props.yield_xe135 * fission_rate;
      xe135m_prod[i] += props.yield_xe135m * fission_rate;
      power += props.fission_q * fission_rate;
    }

    // Determine Xe135 absorption rate
    xe135_abs[i] = absorb_tally->results_(i, 0, 1) / absorb_tally->n_realizations_;
  }

  // Now that the raw production/absorption rates and total power have been
  // determined, we need to normalize the values by the total power and update
  // densities

  for (gsl::index i = 0; i < n_mats; ++i) {
    // Normalize production rates by specified power and volume. Note that
    // production terms need to be in units of [1/(s-b-cm)] since nuclide
    // densities are in units of [atom/b-cm].
    constexpr double BARN_PER_CM_SQ {1.0e24};
    double V = volume[mats[i]];
    double normalization = actual_power / (power * BARN_PER_CM_SQ * V);
    i135_prod[i] *= normalization;
    te135_prod[i] *= normalization;
    xe135_prod[i] *= normalization;
    xe135m_prod[i] *= normalization;
    xe135_abs[i] *= normalization;

    // Get array of densities in current material
    const auto& mat {model::materials[mats[i]]};
    const auto& nucs_in_mat {mat->nuclides()};
    const auto& densities {mat->densities()};

    // Divide Xe135 absorption rate by current Xe135 density to get absorption
    // coefficient in units of [1/s]
    int idx_i135;
    int idx_xe135;
    for (gsl::index j = 0; j < nucs_in_mat.size(); ++j) {
      if (nucname[nucs_in_mat[j]] == "Xe135") {
        xe135_abs[i] /= densities[j];
        idx_xe135 = j;
      } else if (nucname[nucs_in_mat[j]] == "I135") {
        idx_i135 = j;
      }
    }

    // solve for equilibrium concentrations. These are analogous to eq (9) and
    // (10) in Griesheimer, but they explicitly account for Te135 and Xe135_m1.
    double i135_eq = (i135_prod[i] + te135_prod[i]) / i135_decay;
    double xe135_eq = (xe135_prod[i] + b_g*(i135_prod[i] + te135_prod[i]) +
                       b_it*(xe135m_prod[i] + b_m*(i135_prod[i] +
                       te135_prod[i]))) / (xe135_decay + xe135_abs[i]);

    // Create array of pointers to nuclide C-strings
    std::vector<std::string> new_nucnames;
    for (auto index : nucs_in_mat) {
      new_nucnames.push_back(nucname[index]);
    }

    // Copy existing densities and update I135 and Xe135
    std::vector<double> new_densities{densities.begin(), densities.end()};
    new_densities[idx_i135] = i135_eq;
    new_densities[idx_xe135] = xe135_eq;

    // Set densities for material
    mat->set_densities(new_nucnames, new_densities);
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
        fission_tally->reset();
        absorb_tally->reset();
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
} // namespace openmc

int main(int argc, char* argv[]) {
  // Initialize OpenMC
  CHECK(openmc_init(argc, argv, nullptr));

  // Xe feedback initialization
  openmc::xenon::init();

  // Run
  CHECK(openmc::xenon::run());

  // Finalize simulation
  CHECK(openmc_finalize());
}
