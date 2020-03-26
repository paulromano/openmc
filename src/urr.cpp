#include "openmc/math_functions.h"
#include "openmc/urr.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"

#include <iostream>

namespace openmc {

//==============================================================================
// URRData implementation
//==============================================================================

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
  n_energy_ = energy_.size();

  // Read URR tables
  read_dataset(group_id, "table", prob_);
}

//==============================================================================
// Unresolved implementation
//==============================================================================

Unresolved::Unresolved(hid_t group)
{
  // Read the metadata
  std::string tmp;
  read_attribute(group, "case", tmp);
  if (tmp == "A") {
    case_ = Case::A;
  } else if (tmp == "B") {
    case_ = Case::B;
  } else {
    case_ = Case::C;
  }
  read_attribute(group, "add_to_background", add_to_background_);
  read_attribute(group, "energy_min", energy_min_);
  read_attribute(group, "energy_max", energy_max_);
  read_attribute(group, "target_spin", target_spin_);
  read_attribute(group, "atomic_weight_ratio", awr_);

  // Read channel and scattering radii
  channel_radius_ = read_function(group, "channel_radius");
  scattering_radius_ = read_function(group, "scattering_radius");

  // Read energies for energy-dependent cases
  if (case_ != Case::A) {
    read_dataset(group, "energies", energy_);
  }

  // Read unresolved resonance parameters
  xt::xtensor<double, 2> matrix;
  read_dataset(group, "parameters", matrix);
  auto n_rows = matrix.shape()[0];
  for (int i = 0; i < n_rows; ++i) {
    int l = matrix(i, 0);
    int j = matrix(i, 1);

    // New spin sequence (l, j)
    if (ljs_.empty() || j != ljs_.back().j || l != ljs_.back().l) {
      ljs_.emplace_back();
      ljs_.back().l = l;
      ljs_.back().j = j;
    }

    // New set of parameters for this spin sequence/energy
    ljs_.back().params.emplace_back();
    auto& p {ljs_.back().params.back()};
    p.E = matrix(i, 2);
    p.avg_d = matrix(i, 3);
    p.df_x = matrix(i, 4);
    p.df_n = matrix(i, 5);
    p.df_f = matrix(i, 6);
    p.avg_gx = matrix(i, 7);
    p.avg_gn0 = matrix(i, 8);
    p.avg_gg = matrix(i, 9);
    p.avg_gf = matrix(i, 10);
  }
}

void Unresolved::sample_full_ladder(ResonanceLadder* ladder, uint64_t* seed)
  const
{
  ladder->res_.clear();
  ladder->l_values_.clear();

  int i_res = 0;

  for (auto& spin_seq : ljs_) {
    // Keep about 300 resonances on either end of the unresolved resonance
    // range to avoid truncation effects
    double d_min = spin_seq.params.front().avg_d;
    double d_max = spin_seq.params.back().avg_d;
    double E_min = energy_min_ - 300 * d_min;
    double E_max = energy_max_ + 300 * d_max;

    // Select a starting energy with random offset for this spin sequence
    double E = E_min + prn(seed) * d_min;

    // Average parameters interpolated at incident neutron energy
    URParameters p, p_l, p_r;

    // Loop over energies
    for (int i = 0; i < spin_seq.params.size(); ++i) {
      p_l = spin_seq.params[i];

      if (case_ == Case::A) {
        p = p_l;
      } else {
        // Parameters are energy-dependent, so get next values to interpolate
        if (E < energy_min_ || p_l.E >= energy_max_) {
          p_r = p_l;
        } else {
          p_r = spin_seq.params[i+1];
        }
      }

      while (E < E_max) {
        // Interpolate average parameters
        if (case_ != Case::A) {
          double f = p_l.E == p_r.E ? 0 : (E - p_l.E) / (p_r.E - p_l.E);
          p.avg_d = p_l.avg_d + f * (p_r.avg_d - p_l.avg_d);
          p.df_x = p_l.df_x + f * (p_r.df_x - p_l.df_x);
          p.df_n = p_l.df_n + f * (p_r.df_n - p_l.df_n);
          p.df_f = p_l.df_f + f * (p_r.df_f - p_l.df_f);
          p.avg_gx = p_l.avg_gx + f * (p_r.avg_gx - p_l.avg_gx);
          p.avg_gn0 = p_l.avg_gn0 + f * (p_r.avg_gn0 - p_l.avg_gn0);
          p.avg_gg = p_l.avg_gg + f * (p_r.avg_gg - p_l.avg_gg);
          p.avg_gf = p_l.avg_gf + f * (p_r.avg_gf - p_l.avg_gf);
        }

        // Create resonance
        ladder->res_.emplace_back();
        auto& res {ladder->res_.back()};
        res.E = E;
        res.l = spin_seq.l;
        res.j = spin_seq.j;

        // Sample fission width
        if (p.avg_gf == 0) {
          res.gf = 0;
        } else {
          double xf = chi_square(p.df_f, seed);
          res.gf = xf * p.avg_gf / p.df_f;
        }

        // Sample competetive width
        if (p.avg_gx == 0) {
          res.gx = 0;
        } else {
          double xx = chi_square(p.df_x, seed);
          res.gx = xx * p.avg_gx / p.df_x;
        }

        // Calculate energy-dependent neutron width
        double xn0 = chi_square(p.df_n, seed);
        double k = wave_number(awr_, E);
        double rho = k * (*channel_radius_)(E);
        std::tie(res.p, res.s) = penetration_shift(res.l, rho);
        res.gn = res.p / rho * std::sqrt(E) * xn0 * p.avg_gn0;

        // Calculate total width
        res.gt = res.gn + p.avg_gg + res.gf + res.gx;

        // Sample level spacing and update energy
        double d = p.avg_d * std::sqrt(-4.*std::log(prn(seed)) / PI);
        E += d;

        // Add the index of this resonance to the map of l-values
        ladder->l_values_[res.l].push_back(i_res);
        i_res++;

        // If the parameters are energy-dependent (Case C) or fission widths
        // are energy-dependent (Case B), get the parameters for the next
        // energy bin for this spin sequence
        if (case_ != Case::A && E > p_r.E && p_l.E < energy_max_) break;
      }
    }
  }
}

void Unresolved::sample_ladder(double energy, ResonanceLadder* ladder,
  uint64_t* seed) const
{
  // Number of resonances to sample
  int n_res = 100;

  ladder->res_.clear();
  ladder->l_values_.clear();

  // Find the energy bin
  int i_grid;
  if (energy < energy_.front()) {
    i_grid = 0;
  } else if (energy >= energy_.back()) {
    i_grid = energy_.size() - 1;
  } else {
    i_grid = lower_bound_index(energy_.begin(), energy_.end(), energy);
  }

  int i_res = 0;
  for (auto& spin_seq : ljs_) {
    // Get the average parameters, interpolated at incident neutron energy if
    // the parameters are energy-dependent
    URParameters p;
    auto p_l = spin_seq.params[i_grid];
    if (case_ == Case::A) {
      p = p_l;
    } else {
      // Parameters are energy-dependent, so get next values to interpolate
      URParameters p_r;
      if (energy <= energy_min_ || energy >= energy_max_) {
        p_r = p_l;
      } else {
        p_r = spin_seq.params[i_grid + 1];
      }
      double f = p_l.E == p_r.E ? 0 : (energy - p_l.E) / (p_r.E - p_l.E);
      p.avg_d = p_l.avg_d + f * (p_r.avg_d - p_l.avg_d);
      p.df_x = p_l.df_x + f * (p_r.df_x - p_l.df_x);
      p.df_n = p_l.df_n + f * (p_r.df_n - p_l.df_n);
      p.df_f = p_l.df_f + f * (p_r.df_f - p_l.df_f);
      p.avg_gx = p_l.avg_gx + f * (p_r.avg_gx - p_l.avg_gx);
      p.avg_gn0 = p_l.avg_gn0 + f * (p_r.avg_gn0 - p_l.avg_gn0);
      p.avg_gg = p_l.avg_gg + f * (p_r.avg_gg - p_l.avg_gg);
      p.avg_gf = p_l.avg_gf + f * (p_r.avg_gf - p_l.avg_gf);
    }

    // Select a starting energy with random offset for this spin sequence
    double E_start = energy + prn(seed) * p.avg_d;
    double E = E_start;

    int i_mid = n_res/2;
    for (int i = 0; i < n_res; ++i) {
      // Create resonance
      ResonanceLadder::Resonance res;
      res.l = spin_seq.l;
      res.j = spin_seq.j;

      // Finished sampling all the resonances to the right; now go to the left
      if (i == i_mid) {
        double d = p.avg_d * std::sqrt(-4.*std::log(prn(seed)) / PI);
        E = E_start - d;
      }

      // Sample fission width
      if (p.avg_gf == 0) {
        res.gf = 0;
      } else {
        double xf = chi_square(p.df_f, seed);
        res.gf = xf * p.avg_gf / p.df_f;
      }

      // Sample competetive width
      if (p.avg_gx == 0) {
        res.gx = 0;
      } else {
        double xx = chi_square(p.df_x, seed);
        res.gx = xx * p.avg_gx / p.df_x;
      }

      // Calculate energy-dependent neutron width
      double xn0 = chi_square(p.df_n, seed);
      double k = wave_number(awr_, energy);
      double rho = k * (*channel_radius_)(energy);
      std::tie(res.p, res.s) = penetration_shift(res.l, rho);
      res.gn = res.p / rho * std::sqrt(energy) * xn0 * p.avg_gn0;

      // Calculate total width
      res.gt = res.gn + p.avg_gg + res.gf + res.gx;

      // Sample level spacing
      double d = p.avg_d * std::sqrt(-4.*std::log(prn(seed)) / PI);

      // Update resonance parameters and energy
      res.E = E;
      ladder->res_.push_back(res);
      if (i < i_mid) {
        E += d;
      } else {
        E -= d;
      }

      // Add the index of this resonance to the map of l-values
      ladder->l_values_[res.l].push_back(i_res);
      i_res++;
    }
  }
}

//==============================================================================
// ResonanceLadder implementation
//==============================================================================

void ResonanceLadder::evaluate(double E, double sqrtkT, double target_spin, double
  awr, Function1D& channel_radius, Function1D& scattering_radius, URRXS& xs) const
{
  using namespace std::complex_literals;

  xs.elastic = 0.;
  xs.capture = 0.;
  xs.fission = 0.;
  xs.competitive = 0.;

  double k = wave_number(awr, E);
  double rho = k * channel_radius(E);
  double rhohat = k * scattering_radius(E);

  // Doppler width
  double d = 2*sqrtkT * std::sqrt(E/awr);

  for (auto& it : l_values_) {
    int l = it.first;
    double p, s;
    std::tie(p, s) = penetration_shift(l, rho);
    double phi = phase_shift(l, rhohat);
    double cos2phi = std::cos(2*phi);
    double sin2phi = std::sin(2*phi);
    double sinphi2 = std::pow(std::sin(phi), 2);

    // Add potential scattering
    xs.elastic += 4*PI / (k*k) * (2*l + 1) * sinphi2;

    for (auto& i_res : it.second) {
      // Get the resonance parameters
      auto& res = res_[i_res];

      // Calculate neutron and total width at energy E
      double gnE = p * res.gn / res.p;
      double gtE = gnE + res.gg + res.gf + res.gx;

      double Eprime = res.E + (res.s - s) / (2*res.p) * res.gn;
      double gJ = (2*res.j + 1) / (4*target_spin + 2);

      // Calculate the symmetric and asymmetric Breit-Wigner line shape
      // functions
      double theta = gtE/d;
      double x = 2*(E - Eprime)/gtE;
      std::complex<double> z = theta*x/2 + theta/2 * 1.0i;
      std::complex<double> w = theta * std::sqrt(PI/2) * faddeeva(z);
      double psi = w.real();
      double chi = w.imag();

      // Calculate common factor for elastic, capture, and fission cross
      // sections
      double f = (4*PI)/(k*k) * gJ * res.gn/gtE;

      // Add contribution to elastic
      xs.elastic += f * (psi * (cos2phi - (1 - gnE/gtE)) + chi * sin2phi);

      // Add contribution to capture
      xs.capture += f * res.gg / gtE * psi;

      // Add contribution to fission
      xs.fission += f * res.gf / gtE * psi;

      // Add contribution to competitive
      xs.competitive += f * res.gx / gtE * psi;
    }
  }
  xs.total = xs.elastic + xs.capture + xs.fission;
}

//==============================================================================
// Non-member functions
//==============================================================================

double wave_number(double A, double E)
{
  // Reduced Planck constant times c in eV-b^0.5
  double hbar_c = 1e4 * PLANCK_C / (2 * PI);
  return A / (A + 1) * std::sqrt(2*MASS_NEUTRON_EV * std::abs(E)) / hbar_c;
}

double phase_shift(int l, double rho)
{
  double phi;

  if (l == 0) {
    phi = rho;
  } else if  (l == 1) {
    phi = rho - std::atan(rho);
  } else if (l == 2) {
    phi = rho - std::atan(3*rho/(3 - rho*rho));
  } else if (l == 3) {
    phi = rho - std::atan((15*rho - std::pow(rho, 3))/(15 - 6*rho*rho));
  } else if (l == 4) {
    phi = rho - std::atan((105*rho - 10*std::pow(rho, 3))/(105 - 45*rho*rho +
      std::pow(rho, 4)));
  }

  return phi;
}

std::pair<double, double> penetration_shift(int l, double rho)
{
  double p, s;

  if (l == 0) {
    return {rho, 0.};
  } else if  (l == 1) {
    double den = 1 + rho*rho;
    p = std::pow(rho, 3)/den;
    s = -1/den;
  } else if (l == 2) {
    double den = 9 + 3*rho*rho + std::pow(rho, 4);
    p = std::pow(rho, 5)/den;
    s = -(18 + 3*rho*rho)/den;
  } else if (l == 3) {
    double den = 225 + 45*rho*rho + 6*std::pow(rho, 4) + std::pow(rho, 6);
    p = std::pow(rho, 7)/den;
    s = -(675 + 90*rho*rho + 6*std::pow(rho, 4))/den;
  } else if (l == 4) {
    double den = 11025 + 1575*rho*rho + 135*std::pow(rho, 4) + 10*std::pow(rho,6)
      + std::pow(rho, 8);
    p = std::pow(rho, 9)/den;
    s = -(44100 + 4725*rho*rho + 270*std::pow(rho, 4) + 10*std::pow(rho, 6))/den;
  }

  return {p, s};
}

double chi_square(int df, uint64_t* seed)
{
  double q = 0.;

  // Sample df random numbers from the standard Normal distribution using
  // Box-Muller and sum their squares to get a sample from the chi-square
  // distribution
  for (int i = 1; i <= df; i += 2) {
    double u1 = prn(seed);
    double u2 = prn(seed);
    double r = std::sqrt(-2*std::log(u1));
    double theta = 2*PI*u2;
    double z1 = r * std::cos(theta);
    q += z1*z1;
    if (i + 1 > df) break;
    double z2 = r * std::sin(theta);
    q += z2*z2;
  }

  return q;
}

}
