#include "openmc/urr.h"

#include "openmc/endf.h"

#include <iostream>

namespace openmc {

UrrData::UrrData(hid_t group_id)
{
  // Read interpolation and other flags
  int interp_temp;
  read_attribute(group_id, "interpolation", interp_temp);
  interp_ = static_cast<Interpolation>(interp_temp);

  // read the metadata
  read_attribute(group_id, "inelastic", inelastic_flag_);
  read_attribute(group_id, "absorption", absorption_flag_);
  int temp_multiply_smooth;
  read_attribute(group_id, "multiply_smooth", temp_multiply_smooth);
  multiply_smooth_ = (temp_multiply_smooth == 1);

  // read the energies at which tables exist
  read_dataset(group_id, "energy", energy_);

  // Set n_energy_
  n_energy_ = energy_.shape()[0];

  // Read URR tables
  read_dataset(group_id, "table", prob_);
}

Unresolved::Unresolved(hid_t group)
{
  // Read the metadata
  char c;
  read_attribute(group, "case", c);
  if (c == 'A') {
    case_ = Case::A;
  } else if (c == 'B') {
    case_ = Case::B;
  } else {
    case_ = Case::C;
  }
  read_attribute(group, "add_to_background", add_to_background_);
  read_attribute(group, "energy_min", energy_min_);
  read_attribute(group, "energy_max", energy_max_);
  read_attribute(group, "target_spin", target_spin_);

  // Read channel and scattering radii
  channel_radius_ = read_function(group, "channel_radius");
  scattering_radius_ = read_function(group, "scattering_radius");

  // Read energies for energy-dependent cases
  if (case_ != Case::A) {
    read_dataset(group, "energies", energy_);
  }

  // Read parameters
  read_dataset(group, "parameters", params_);
}

ResonanceLadder Unresolved::sample(double E) const
{
  // TODO
}

double ResonanceLadder::evaluate(double E, double T) const
{
  // TODO
}

}
