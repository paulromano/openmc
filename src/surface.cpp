#include "openmc/surface.h"
#include <chrono>
#include <array>
#include <cmath>
#include <utility>

#include <fmt/core.h>
#include <gsl/gsl>

#include <iomanip>
#include "openmc/error.h"
#include "openmc/dagmc.h"
#include "openmc/hdf5_interface.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/math_functions.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::vector<std::unique_ptr<Surface>> surfaces;
  std::unordered_map<int, int> surface_map;
} // namespace model

//==============================================================================
// Helper functions for reading the "coeffs" node of an XML surface element
//==============================================================================

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 1) {
    fatal_error(fmt::format("Surface {} expects 1 coeff but was given {}",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf", &c1);
  if (stat != 1) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 3) {
    fatal_error(fmt::format("Surface {} expects 3 coeffs but was given {}",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf", &c1, &c2, &c3);
  if (stat != 3) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3, double &c4)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 4) {
    fatal_error(fmt::format("Surface {} expects 4 coeffs but was given ",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf", &c1, &c2, &c3, &c4);
  if (stat != 4) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
		 double &c3, double &c4, double &c5, bool strict_count)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);

  if (n_words != 5 && strict_count ) {
    std::stringstream err_msg;
    err_msg << "Surface " << surf_id << " expects 5 coeffs but was given "
	    << n_words;
    fatal_error(err_msg);
  } else if (n_words != 5 && !strict_count ) {
    // fall back to reading the 4 type
    read_coeffs(surf_node, surf_id, c1, c2, c3, c4);
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf %lf", &c1, &c2, &c3, &c4, &c5);
  if (stat != 5 && strict_count ) {
      fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3, double &c4, double &c5, double &c6)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 6) {
    std::stringstream err_msg;
    err_msg << "Surface " << surf_id << " expects 6 coeffs but was given "
            << n_words;
    fatal_error(err_msg);
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf %lf %lf", &c1, &c2, &c3, &c4, &c5, &c6);
  if (stat != 6) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, int surf_id, double &c1, double &c2,
                 double &c3, double &c4, double &c5, double &c6, double &c7,
                 double &c8, double &c9, double &c10)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 10) {
    fatal_error(fmt::format("Surface {} expects 10 coeffs but was given {}",
      surf_id, n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10);
  if (stat != 10) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

//==============================================================================
// Surface implementation
//==============================================================================

Surface::Surface() {} // empty constructor

Surface::Surface(pugi::xml_node surf_node)
{
  if (check_for_node(surf_node, "id")) {
    id_ = std::stoi(get_node_value(surf_node, "id"));
  } else {
    fatal_error("Must specify id of surface in geometry XML file.");
  }

  if (check_for_node(surf_node, "name")) {
    name_ = get_node_value(surf_node, "name", false);
  }

  if (check_for_node(surf_node, "boundary")) {
    std::string surf_bc = get_node_value(surf_node, "boundary", true, true);

    if (surf_bc == "transmission" || surf_bc == "transmit" ||surf_bc.empty()) {
      bc_ = BoundaryType::TRANSMIT;

    } else if (surf_bc == "vacuum") {
      bc_ = BoundaryType::VACUUM;

    } else if (surf_bc == "reflective" || surf_bc == "reflect"
               || surf_bc == "reflecting") {
      bc_ = BoundaryType::REFLECT;
    } else if (surf_bc == "white") {
      bc_ = BoundaryType::WHITE;
    } else if (surf_bc == "periodic") {
      bc_ = BoundaryType::PERIODIC;
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
        "on surface {}", surf_bc, id_));
    }

  } else {
    bc_ = BoundaryType::TRANSMIT;
  }

}

bool
Surface::sense(Position r, Direction u) const
{
  // Evaluate the surface equation at the particle's coordinates to determine
  // which side the particle is on.
  const double f = evaluate(r);

  // Check which side of surface the point is on.
  if (std::abs(f) < FP_COINCIDENT) {
    // Particle may be coincident with this surface. To determine the sense, we
    // look at the direction of the particle relative to the surface normal (by
    // default in the positive direction) via their dot product.
    return u.dot(normal(r)) > 0.0;
  }
  return f > 0.0;
}

Direction
Surface::reflect(Position r, Direction u, Particle* p) const
{
  // Determine projection of direction onto normal and squared magnitude of
  // normal.
  Direction n = normal(r);
  const double projection = n.dot(u);
  const double magnitude = n.dot(n);

  // Reflect direction according to normal.
  return u -= (2.0 * projection / magnitude) * n;
}

Direction
Surface::diffuse_reflect(Position r, Direction u, uint64_t* seed) const
{
  // Diffuse reflect direction according to the normal.
  // cosine distribution

  Direction n = this->normal(r);
  n /= n.norm();
  const double projection = n.dot(u);

  // sample from inverse function, u=sqrt(rand) since p(u)=2u, so F(u)=u^2
  const double mu = (projection>=0.0) ?
                  -std::sqrt(prn(seed)) : std::sqrt(prn(seed));

  // sample azimuthal distribution uniformly
  u = rotate_angle(n, mu, nullptr, seed);

  // normalize the direction
  return u/u.norm();
}

CSGSurface::CSGSurface() : Surface{} {};
CSGSurface::CSGSurface(pugi::xml_node surf_node) : Surface{surf_node} {};

void
CSGSurface::to_hdf5(hid_t group_id) const
{
  std::string group_name {"surface "};
  group_name += std::to_string(id_);

  hid_t surf_group = create_group(group_id, group_name);

  switch(bc_) {
    case BoundaryType::TRANSMIT :
      write_string(surf_group, "boundary_type", "transmission", false);
      break;
    case BoundaryType::VACUUM :
      write_string(surf_group, "boundary_type", "vacuum", false);
      break;
    case BoundaryType::REFLECT :
      write_string(surf_group, "boundary_type", "reflective", false);
      break;
    case BoundaryType::WHITE :
      write_string(surf_group, "boundary_type", "white", false);
      break;
    case BoundaryType::PERIODIC :
      write_string(surf_group, "boundary_type", "periodic", false);
      break;
  }

  if (!name_.empty()) {
    write_string(surf_group, "name", name_, false);
  }

  to_hdf5_inner(surf_group);

  close_group(surf_group);
}

//==============================================================================
// DAGSurface implementation
//==============================================================================
#ifdef DAGMC
DAGSurface::DAGSurface() : Surface{} {} // empty constructor

double DAGSurface::evaluate(Position r) const
{
  return 0.0;
}

double
DAGSurface::distance(Position r, Direction u, bool coincident) const
{
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  moab::EntityHandle hit_surf;
  double dist;
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  rval = dagmc_ptr_->ray_fire(surf, pnt, dir, hit_surf, dist, NULL, 0, 0);
  MB_CHK_ERR_CONT(rval);
  if (dist < 0.0) dist = INFTY;
  return dist;
}

Direction DAGSurface::normal(Position r) const
{
  moab::ErrorCode rval;
  moab::EntityHandle surf = dagmc_ptr_->entity_by_index(2, dag_index_);
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3];
  rval = dagmc_ptr_->get_angle(surf, pnt, dir);
  MB_CHK_ERR_CONT(rval);
  return dir;
}

Direction DAGSurface::reflect(Position r, Direction u, Particle* p) const
{
  Expects(p);
  p->history_.reset_to_last_intersection();
  p->last_dir_ = Surface::reflect(r, u, p);
  return p->last_dir_;
}

void DAGSurface::to_hdf5(hid_t group_id) const {}

#endif
//==============================================================================
// PeriodicSurface implementation
//==============================================================================

PeriodicSurface::PeriodicSurface(pugi::xml_node surf_node)
  : CSGSurface {surf_node}
{
  if (check_for_node(surf_node, "periodic_surface_id")) {
    i_periodic_ = std::stoi(get_node_value(surf_node, "periodic_surface_id"));
  }
}

//==============================================================================
// Generic functions for x-, y-, and z-, planes.
//==============================================================================

// The template parameter indicates the axis normal to the plane.
template<int i> double
axis_aligned_plane_distance(Position r, Direction u, bool coincident, double offset)
{
  const double f = offset - r[i];
  if (coincident || std::abs(f) < FP_COINCIDENT || u[i] == 0.0) return INFTY;
  const double d = f / u[i];
  if (d < 0.0) return INFTY;
  return d;
}

//==============================================================================
// SurfaceXPlane implementation
//==============================================================================

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_);
}

double SurfaceXPlane::evaluate(Position r) const
{
  return r.x - x0_;
}

double SurfaceXPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<0>(r, u, coincident, x0_);
}

Direction SurfaceXPlane::normal(Position r) const
{
  return {1., 0., 0.};
}

void SurfaceXPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-plane", false);
  std::array<double, 1> coeffs {{x0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfaceXPlane::periodic_translate(const PeriodicSurface* other,
                                       Position& r,  Direction& u) const
{
  Direction other_n = other->normal(r);
  if (other_n.x == 1 && other_n.y == 0 && other_n.z == 0) {
    r.x = x0_;
    return false;
  } else {
    // Assume the partner is an YPlane (the only supported partner).  Use the
    // evaluate function to find y0_, then adjust position/Direction for
    // rotational symmetry.
    double y0_ = -other->evaluate({0., 0., 0.});
    r.y = r.x - x0_ + y0_;
    r.x = x0_;

    double ux = u.x;
    u.x = -u.y;
    u.y = ux;

    return true;
  }
}

BoundingBox
SurfaceXPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {x0_, INFTY, -INFTY, INFTY, -INFTY, INFTY};
  } else {
    return {-INFTY, x0_, -INFTY, INFTY, -INFTY, INFTY};
  }
}

//==============================================================================
// SurfaceYPlane implementation
//==============================================================================

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id_, y0_);
}

double SurfaceYPlane::evaluate(Position r) const
{
  return r.y - y0_;
}

double SurfaceYPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<1>(r, u, coincident, y0_);
}

Direction SurfaceYPlane::normal(Position r) const
{
  return {0., 1., 0.};
}

void SurfaceYPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-plane", false);
  std::array<double, 1> coeffs {{y0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfaceYPlane::periodic_translate(const PeriodicSurface* other,
                                       Position& r, Direction& u) const
{
  Direction other_n = other->normal(r);
  if (other_n.x == 0 && other_n.y == 1 && other_n.z == 0) {
    // The periodic partner is also aligned along y.  Just change the y coord.
    r.y = y0_;
    return false;
  } else {
    // Assume the partner is an XPlane (the only supported partner).  Use the
    // evaluate function to find x0_, then adjust position/Direction for rotational
    // symmetry.
    double x0_ = -other->evaluate({0., 0., 0.});
    r.x = r.y - y0_ + x0_;
    r.y = y0_;

    double ux = u.x;
    u.x = u.y;
    u.y = -ux;

    return true;
  }
}

BoundingBox
SurfaceYPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {-INFTY, INFTY, y0_, INFTY, -INFTY, INFTY};
  } else {
    return {-INFTY, INFTY, -INFTY, y0_, -INFTY, INFTY};
  }
}

//==============================================================================
// SurfaceZPlane implementation
//==============================================================================

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id_, z0_);
}

double SurfaceZPlane::evaluate(Position r) const
{
  return r.z - z0_;
}

double SurfaceZPlane::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_plane_distance<2>(r, u, coincident, z0_);
}

Direction SurfaceZPlane::normal(Position r) const
{
  return {0., 0., 1.};
}

void SurfaceZPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-plane", false);
  std::array<double, 1> coeffs {{z0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfaceZPlane::periodic_translate(const PeriodicSurface* other,
                                       Position& r, Direction& u) const
{
  // Assume the other plane is aligned along z.  Just change the z coord.
  r.z = z0_;
  return false;
}

BoundingBox
SurfaceZPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {-INFTY, INFTY, -INFTY, INFTY, z0_, INFTY};
  } else {
    return {-INFTY, INFTY, -INFTY, INFTY, -INFTY, z0_};
  }
}

//==============================================================================
// SurfacePlane implementation
//==============================================================================

SurfacePlane::SurfacePlane(pugi::xml_node surf_node)
  : PeriodicSurface(surf_node)
{
  read_coeffs(surf_node, id_, A_, B_, C_, D_);
}

double
SurfacePlane::evaluate(Position r) const
{
  return A_*r.x + B_*r.y + C_*r.z - D_;
}

double
SurfacePlane::distance(Position r, Direction u, bool coincident) const
{
  const double f = A_*r.x + B_*r.y + C_*r.z - D_;
  const double projection = A_*u.x + B_*u.y + C_*u.z;
  if (coincident || std::abs(f) < FP_COINCIDENT || projection == 0.0) {
    return INFTY;
  } else {
    const double d = -f / projection;
    if (d < 0.0) return INFTY;
    return d;
  }
}

Direction
SurfacePlane::normal(Position r) const
{
  return {A_, B_, C_};
}

void SurfacePlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "plane", false);
  std::array<double, 4> coeffs {{A_, B_, C_, D_}};
  write_dataset(group_id, "coefficients", coeffs);
}

bool SurfacePlane::periodic_translate(const PeriodicSurface* other, Position& r,
                                      Direction& u) const
{
  // This function assumes the other plane shares this plane's normal direction.

  // Determine the distance to intersection.
  double d = evaluate(r) / (A_*A_ + B_*B_ + C_*C_);

  // Move the particle that distance along the normal vector.
  r.x -= d * A_;
  r.y -= d * B_;
  r.z -= d * C_;

  return false;
}

//==============================================================================
// Generic functions for x-, y-, and z-, cylinders
//==============================================================================

// The template parameters indicate the axes perpendicular to the axis of the
// cylinder.  offset1 and offset2 should correspond with i1 and i2,
// respectively.
template<int i1, int i2> double
axis_aligned_cylinder_evaluate(Position r, double offset1,
                               double offset2, double radius)
{
  const double r1 = r[i1] - offset1;
  const double r2 = r[i2] - offset2;
  return r1*r1 + r2*r2 - radius*radius;
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cylinder_distance(Position r, Direction u,
     bool coincident, double offset1, double offset2, double radius)
{
  const double a = 1.0 - u[i1]*u[i1];  // u^2 + v^2
  if (a == 0.0) return INFTY;

  const double r2 = r[i2] - offset1;
  const double r3 = r[i3] - offset2;
  const double k = r2 * u[i2] + r3 * u[i3];
  const double c = r2*r2 + r3*r3 - radius*radius;
  const double quad = k*k - a*c;

  if (quad < 0.0) {
    // No intersection with cylinder.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the cylinder, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) {
      return INFTY;
    } else {
      return (-k + sqrt(quad)) / a;
    }

  } else if (c < 0.0) {
    // Particle is inside the cylinder, thus one distance must be negative
    // and one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad).
    return (-k + sqrt(quad)) / a;

  } else {
    // Particle is outside the cylinder, thus both distances are either
    // positive or negative. If positive, the smaller distance is the one
    // with positive sign on sqrt(quad).
    const double d = (-k - sqrt(quad)) / a;
    if (d < 0.0) return INFTY;
    return d;
  }
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3> Direction
axis_aligned_cylinder_normal(Position r, double offset1, double offset2)
{
  Direction u;
  u[i2] = 2.0 * (r[i2] - offset1);
  u[i3] = 2.0 * (r[i3] - offset2);
  u[i1] = 0.0;
  return u;
}

//==============================================================================
// SurfaceXCylinder implementation
//==============================================================================

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, y0_, z0_, radius_);
}

double SurfaceXCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<1, 2>(r, y0_, z0_, radius_);
}

double SurfaceXCylinder::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<0, 1, 2>(r, u, coincident, y0_, z0_,
                                                 radius_);
}

Direction SurfaceXCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<0, 1, 2>(r, y0_, z0_);
}

void SurfaceXCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cylinder", false);
  std::array<double, 3> coeffs {{y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceXCylinder::bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {-INFTY, INFTY, y0_ - radius_, y0_ + radius_, z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}
//==============================================================================
// SurfaceYCylinder implementation
//==============================================================================

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, z0_, radius_);
}

double SurfaceYCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 2>(r, x0_, z0_, radius_);
}

double SurfaceYCylinder::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<1, 0, 2>(r, u, coincident, x0_, z0_,
                                                 radius_);
}

Direction SurfaceYCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<1, 0, 2>(r, x0_, z0_);
}

void SurfaceYCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cylinder", false);
  std::array<double, 3> coeffs {{x0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceYCylinder::bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, -INFTY, INFTY, z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}

//==============================================================================
// SurfaceZCylinder implementation
//==============================================================================

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, radius_);
}

double SurfaceZCylinder::evaluate(Position r) const
{
  return axis_aligned_cylinder_evaluate<0, 1>(r, x0_, y0_, radius_);
}

double SurfaceZCylinder::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cylinder_distance<2, 0, 1>(r, u, coincident, x0_, y0_,
                                                 radius_);
}

Direction SurfaceZCylinder::normal(Position r) const
{
  return axis_aligned_cylinder_normal<2, 0, 1>(r, x0_, y0_);
}

void SurfaceZCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cylinder", false);
  std::array<double, 3> coeffs {{x0_, y0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceZCylinder::bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, y0_ - radius_, y0_ + radius_, -INFTY, INFTY};
  } else {
    return {};
  }
}


//==============================================================================
// SurfaceSphere implementation
//==============================================================================

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_);
}

double SurfaceSphere::evaluate(Position r) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  return x*x + y*y + z*z - radius_*radius_;
}

double SurfaceSphere::distance(Position r, Direction u, bool coincident) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  const double k = x*u.x + y*u.y + z*u.z;
  const double c = x*x + y*y + z*z - radius_*radius_;
  const double quad = k*k - c;

  if (quad < 0.0) {
    // No intersection with sphere.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the sphere, thus one distance is positive/negative and
    // the other is zero. The sign of k determines if we are facing in or out.
    if (k >= 0.0) {
      return INFTY;
    } else {
      return -k + sqrt(quad);
    }

  } else if (c < 0.0) {
    // Particle is inside the sphere, thus one distance must be negative and
    // one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad)
    return -k + sqrt(quad);

  } else {
    // Particle is outside the sphere, thus both distances are either positive
    // or negative. If positive, the smaller distance is the one with positive
    // sign on sqrt(quad).
    const double d = -k - sqrt(quad);
    if (d < 0.0) return INFTY;
    return d;
  }
}

Direction SurfaceSphere::normal(Position r) const
{
  return {2.0*(r.x - x0_), 2.0*(r.y - y0_), 2.0*(r.z - z0_)};
}

void SurfaceSphere::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "sphere", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox SurfaceSphere::bounding_box(bool pos_side) const {
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_,
            y0_ - radius_, y0_ + radius_,
            z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}

//==============================================================================
// Generic functions for x-, y-, and z-, cones
//==============================================================================

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cone_evaluate(Position r, double offset1,
                           double offset2, double offset3, double radius_sq)
{
  const double r1 = r[i1] - offset1;
  const double r2 = r[i2] - offset2;
  const double r3 = r[i3] - offset3;
  return r2*r2 + r3*r3 - radius_sq*r1*r1;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> double
axis_aligned_cone_distance(Position r, Direction u,
     bool coincident, double offset1, double offset2, double offset3,
     double radius_sq)
{
  const double r1 = r[i1] - offset1;
  const double r2 = r[i2] - offset2;
  const double r3 = r[i3] - offset3;
  const double a = u[i2]*u[i2] + u[i3]*u[i3]
                   - radius_sq*u[i1]*u[i1];
  const double k = r2*u[i2] + r3*u[i3] - radius_sq*r1*u[i1];
  const double c = r2*r2 + r3*r3 - radius_sq*r1*r1;
  double quad = k*k - a*c;

  double d;

  if (quad < 0.0) {
    // No intersection with cone.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the cone, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    const double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0) d = b;
    } else {
      if (b > 0.0) {
        if (b < d) d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0) return INFTY;
  return d;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3> Direction
axis_aligned_cone_normal(Position r, double offset1, double offset2,
                         double offset3, double radius_sq)
{
  Direction u;
  u[i1] = -2.0 * radius_sq * (r[i1] - offset1);
  u[i2] = 2.0 * (r[i2] - offset2);
  u[i3] = 2.0 * (r[i3] - offset3);
  return u;
}

//==============================================================================
// SurfaceXCone implementation
//==============================================================================

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_);
}

double SurfaceXCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

double SurfaceXCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<0, 1, 2>(r, u, coincident, x0_, y0_, z0_,
                                             radius_sq_);
}

Direction SurfaceXCone::normal(Position r) const
{
  return axis_aligned_cone_normal<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

void SurfaceXCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCone implementation
//==============================================================================

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_);
}

double SurfaceYCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

double SurfaceYCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<1, 0, 2>(r, u, coincident, y0_, x0_, z0_,
                                             radius_sq_);
}

Direction SurfaceYCone::normal(Position r) const
{
  return axis_aligned_cone_normal<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

void SurfaceYCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCone implementation
//==============================================================================

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, radius_sq_);
}

double SurfaceZCone::evaluate(Position r) const
{
  return axis_aligned_cone_evaluate<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

double SurfaceZCone::distance(Position r, Direction u, bool coincident) const
{
  return axis_aligned_cone_distance<2, 0, 1>(r, u, coincident, z0_, x0_, y0_,
                                             radius_sq_);
}

Direction SurfaceZCone::normal(Position r) const
{
  return axis_aligned_cone_normal<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

void SurfaceZCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceQuadric implementation
//==============================================================================

SurfaceQuadric::SurfaceQuadric(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, A_, B_, C_, D_, E_, F_, G_, H_, J_, K_);
}

double
SurfaceQuadric::evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;
  return x*(A_*x + D_*y + G_) +
         y*(B_*y + E_*z + H_) +
         z*(C_*z + F_*x + J_) + K_;
}

double
SurfaceQuadric::distance(Position r, Direction ang, bool coincident) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  const double &u = ang.x;
  const double &v = ang.y;
  const double &w = ang.z;

  const double a = A_*u*u + B_*v*v + C_*w*w + D_*u*v + E_*v*w + F_*u*w;
  const double k = A_*u*x + B_*v*y + C_*w*z + 0.5*(D_*(u*y + v*x)
                   + E_*(v*z + w*y) + F_*(w*x + u*z) + G_*u + H_*v + J_*w);
  const double c = A_*x*x + B_*y*y + C_*z*z + D_*x*y + E_*y*z +  F_*x*z + G_*x
                   + H_*y + J_*z + K_;
  double quad = k*k - a*c;

  double d;

  if (quad < 0.0) {
    // No intersection with surface.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the surface, thus one distance is positive/negative and
    // the other is zero. The sign of k determines which distance is zero and
    // which is not.
    if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0) d = b;
    } else {
      if (b > 0.0) {
        if (b < d) d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0) return INFTY;
  return d;
}

Direction
SurfaceQuadric::normal(Position r) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  return {2.0*A_*x + D_*y + F_*z + G_,
          2.0*B_*y + D_*x + E_*z + H_,
          2.0*C_*z + E_*y + F_*x + J_};
}

void SurfaceQuadric::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "quadric", false);
  std::array<double, 10> coeffs {{A_, B_, C_, D_, E_, F_, G_, H_, J_, K_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// Generic functions for cubic & quartic solver
//==============================================================================

int cubic_solve(double a, double b, double c, double d, std::array<double,3> &x) {

  double a2 = a*a;
  double q  = (a2 - 3*b)/9;
	double r  = (a*(2*a2-9*b) + 27*c)/54;
  double r2 = r*r;
	double q3 = q*q*q;
	double A,B;
  if(r2<q3) {
    double t=r/sqrt(q3);
    if( t<-1) t=-1;
    if( t> 1) t= 1;
    t=acos(t);
    a/=3; q=-2*sqrt(q);
    x[0]=q*std::cos(t/3)-a;
    x[1]=q*std::cos((t+2*PI)/3)-a;
    x[2]=q*std::cos((t-2*PI)/3)-a;
    return 3;
  } else {
    A =-std::pow(std::fabs(r)+std::sqrt(r2-q3),1./3);
    if( r<0 ) A=-A;
    B = (0==A ? 0 : q/A);
          
		a/=3;
		x[0] =(A+B)-a;
		x[1] =-0.5*(A+B)-a;
		x[2] = 0.5*sqrt(3.)*(A-B);
		if(std::fabs(x[2])<1e-32) { x[2]=x[1]; return 2; }
		
		return 1;
  }
}


void quartic_solve(double a, double b, double c, double d, double e, std::array<double,4> &roots) {
  
  double descrim = 256*std::pow(a,3)*std::pow(e,3) - 192*std::pow(a,2)*b*d*std::pow(e,2) - 128*std::pow(a,2)*std::pow(c,2)*std::pow(e,2) + 144*std::pow(a,2)*c*std::pow(d,2)*e -27*std::pow(a,2)*std::pow(d,4) + 144*a*std::pow(b,2)*c*std::pow(e,2) - 6*a*std::pow(b,2)*std::pow(d,2)*e - 80*a*b*std::pow(c,2)*d*e + 18*a*b*c*std::pow(d,3) + 16*a*std::pow(c,4)*e - 4*a*std::pow(c,3)*std::pow(d,2) - 27*std::pow(b,4)*std::pow(e,2) + 18*std::pow(b,3)*c*d*e - 4*std::pow(b,3)*std::pow(d,3) - 4*std::pow(b,2)*std::pow(c,3)*e + std::pow(b,2)*std::pow(c,2)*std::pow(d,2);
  double P = 8*a*c - 3*std::pow(b,2);
  double R = std::pow(b,3) + 8*d*std::pow(a,2) - 4*a*b*c;
  double d_0 = std::pow(c,2) - 3*b*d + 12*a*e;
  double D = 64*std::pow(a,3)*e - 16*std::pow(a,2)*std::pow(c,2) + 16*a*std::pow(b,2)*c - 16*std::pow(a,2)*b*d - 3*std::pow(b,4);
  if(descrim > 0 ) { // 4 real roots
    double p = P/(8*std::pow(a,2));
    double q = R/(8*std::pow(a,3));

    double delta0 = std::pow(c,2) - 3*b*d + 12*a*e;
    double delta1 = 2*std::pow(c,3) - 9*b*c*d + 27*std::pow(b,2)*e + 27*a*std::pow(d,2) - 72*a*c*e;

    double phi = std::acos(delta1/(2.*std::sqrt(std::pow(delta0,3))));
    double S = 0.5*sqrt(-2./3.*p + 2./(3.*a)*std::sqrt(delta0)*std::cos(phi/3.));
    double coeff = std::sqrt( -4*std::pow(S,2) - 2*p + q/S);
    roots[0] = -b/(4*a) - S + 0.5*coeff;
    roots[1] = -b/(4*a) - S - 0.5*coeff;
    roots[2] = -b/(4*a) + S + 0.5*coeff;
    roots[3] = -b/(4*a) + S - 0.5*coeff;
    std::sort(roots.begin(),roots.end());
  }

  
  // depressed quartic coefficients 
  double f = c - (3*std::pow(b,2)/8);
  double g = d + std::pow(b,3)/8 - (b*c/2);
  double h = e - (3*std::pow(b,4)/256) + (std::pow(b,2)*c/16) - (b*d/4);


  double x1 = f/2;
  double x2 = std::pow(f,2)-4*h/16;
  double x3 = std::pow(g,2)/64;

  std::array<double,3> cubic_roots;
  int n_roots = cubic_solve(1,x1,x2,x3,cubic_roots);
  std::cout << cubic_roots[0] << " ";
  std::cout << cubic_roots[1] << " ";
  std::cout << cubic_roots[2] << std::endl;  

  double r1,r2;

  if ( n_roots == 3 ) 

  if(cubic_roots[0] > 0 && cubic_roots[1] > 0 ) {
    r1 = cubic_roots[0];
    r2 = cubic_roots[1];
  } else if (cubic_roots[1] > 0 && cubic_roots[2] > 0  ) {
    r1 = cubic_roots[1];
    r2 = cubic_roots[2];
  } else if (cubic_roots[0] > 0 && cubic_roots[2] > 0  ) {
    r1 = cubic_roots[0];
    r2 = cubic_roots[2];
  } 
    
  double p = std::sqrt(r1);
  double q = std::sqrt(r2);
  double r = -g/(8*p*q);
  double s = b/(4*a);

  std::cout << p << " " << q << " " << r << " " << s << std::endl;

    roots[0] =  p + q + r - s;
    roots[1] =  p - q - r - s;
    roots[2] = -p + q - r - s; 
    roots[3] = -p - q + r - s;
  
  return;
}

//==============================================================================
// Generic functions for x-, y-, and z-, tori
//==============================================================================

SurfaceXTorus::SurfaceXTorus(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, A_, B_, C_);
}

// todo should do some alebgra first to get rid of the sqrt

double
SurfaceXTorus::evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;

  const double x_coeff2 = std::pow(x-x0_,2);
  const double y_coeff2 = std::pow(y-y0_,2);
  const double z_coeff2 = std::pow(z-z0_,2);

  const double B2 = B_*B_;
  const double C2 = C_*C_;

  return x_coeff2/B2 + std::pow(std::sqrt(y_coeff2 + z_coeff2)-A_,2)/C2 - 1.;
}

double
SurfaceXTorus::distance(Position r, Direction ang, bool coincident) const
{
  return 0;
}

Direction
SurfaceXTorus::normal(Position r) const
{
  return {0,0,0};
}

SurfaceYTorus::SurfaceYTorus(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, A_, B_, C_);
}

void SurfaceXTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_ }};
  write_dataset(group_id, "coefficients", coeffs);
}

double
SurfaceYTorus::evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;

  const double x_coeff2 = std::pow(x-x0_,2);
  const double y_coeff2 = std::pow(y-y0_,2);
  const double z_coeff2 = std::pow(z-z0_,2);

  const double B2 = B_*B_;
  const double C2 = C_*C_;

  return y_coeff2/B2 + std::pow(std::sqrt(x_coeff2 + z_coeff2)-A_,2)/C2 - 1.;
}

double
SurfaceYTorus::distance(Position r, Direction ang, bool coincident) const
{
  return 0;
}

Direction
SurfaceYTorus::normal(Position r) const
{
  return {0,0,0};
}

void SurfaceYTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_ }};
  write_dataset(group_id, "coefficients", coeffs);
}

SurfaceZTorus::SurfaceZTorus(pugi::xml_node surf_node)
  : CSGSurface(surf_node)
{
  read_coeffs(surf_node, id_, x0_, y0_, z0_, A_, B_, C_);
}

double
SurfaceZTorus::evaluate(Position r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;

  const double x_coeff2 = std::pow(x-x0_,2);
  const double y_coeff2 = std::pow(y-y0_,2);
  const double z_coeff2 = std::pow(z-z0_,2);

  const double B2 = B_*B_;
  const double C2 = C_*C_;

  return z_coeff2/B2 + std::pow(std::sqrt(x_coeff2 + y_coeff2)-A_,2)/C2 - 1.;
}

double
SurfaceZTorus::distance(Position r, Direction ang, bool coincident) const
{
  double u = ang.x;
  double v = ang.y;
  double w = ang.z;

  double xr = r.x;
  double yr = r.y;
  double zr = r.z;

  double A = A_;
  double B = B_;
  double C = C_;

  double x0 = x0_;
  double y0 = y0_;
  double z0 = z0_;

  double a,b,c,d,e; // coefficients of quartic

  // 4th order coefficient
  a = std::pow(u,4) + 2*std::pow(u,2)*std::pow(v,2) + std::pow(v,4) + 2*std::pow(C,2)*std::pow(u,2)*std::pow(w,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(v,2)*std::pow(w,2)/std::pow(B,2) + std::pow(C,4)*std::pow(w,4)/std::pow(B,4);

  // 3rd order coefficient
  b = -4*std::pow(u,3)*x0 + 4*std::pow(u,3)*xr - 4*std::pow(u,2)*v*y0 + 4*std::pow(u,2)*v*yr - 4*u*std::pow(v,2)*x0 + 4*u*std::pow(v,2)*xr - 4*std::pow(v,3)*y0 + 4*std::pow(v,3)*yr - 4*std::pow(C,2)*std::pow(u,2)*w*z0/std::pow(B,2) + 4*std::pow(C,2)*std::pow(u,2)*w*zr/std::pow(B,2) - 4*std::pow(C,2)*u*std::pow(w,2)*x0/std::pow(B,2) + 4*std::pow(C,2)*u*std::pow(w,2)*xr/std::pow(B,2) - 4*std::pow(C,2)*std::pow(v,2)*w*z0/std::pow(B,2) + 4*std::pow(C,2)*std::pow(v,2)*w*zr/std::pow(B,2) - 4*std::pow(C,2)*v*std::pow(w,2)*y0/std::pow(B,2) + 4*std::pow(C,2)*v*std::pow(w,2)*yr/std::pow(B,2) - 4*std::pow(C,4)*std::pow(w,3)*z0/std::pow(B,4) + 4*std::pow(C,4)*std::pow(w,3)*zr/std::pow(B,4);
  // 2nd order coefficient
  c= -2*std::pow(A,2)*std::pow(u,2) - 2*std::pow(A,2)*std::pow(v,2) + 2*std::pow(A,2)*std::pow(C,2)*std::pow(w,2)/std::pow(B,2) - 2*std::pow(C,2)*std::pow(u,2) - 2*std::pow(C,2)*std::pow(v,2) + 6*std::pow(u,2)*std::pow(x0,2) - 12*std::pow(u,2)*x0*xr + 6*std::pow(u,2)*std::pow(xr,2) + 2*std::pow(u,2)*std::pow(y0,2) - 4*std::pow(u,2)*y0*yr + 2*std::pow(u,2)*std::pow(yr,2) + 8*u*v*x0*y0 - 8*u*v*x0*yr - 8*u*v*xr*y0 + 8*u*v*xr*yr + 2*std::pow(v,2)*std::pow(x0,2) - 4*std::pow(v,2)*x0*xr + 2*std::pow(v,2)*std::pow(xr,2) + 6*std::pow(v,2)*std::pow(y0,2) - 12*std::pow(v,2)*y0*yr + 6*std::pow(v,2)*std::pow(yr,2) - 2*std::pow(C,4)*std::pow(w,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(u,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(u,2)*z0*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(u,2)*std::pow(zr,2)/std::pow(B,2) + 8*std::pow(C,2)*u*w*x0*z0/std::pow(B,2) - 8*std::pow(C,2)*u*w*x0*zr/std::pow(B,2) - 8*std::pow(C,2)*u*w*xr*z0/std::pow(B,2) + 8*std::pow(C,2)*u*w*xr*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(v,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(v,2)*z0*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(v,2)*std::pow(zr,2)/std::pow(B,2) + 8*std::pow(C,2)*v*w*y0*z0/std::pow(B,2) - 8*std::pow(C,2)*v*w*y0*zr/std::pow(B,2) - 8*std::pow(C,2)*v*w*yr*z0/std::pow(B,2) + 8*std::pow(C,2)*v*w*yr*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(w,2)*std::pow(x0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(w,2)*x0*xr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(w,2)*std::pow(xr,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(w,2)*std::pow(y0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(w,2)*y0*yr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(w,2)*std::pow(yr,2)/std::pow(B,2) + 6*std::pow(C,4)*std::pow(w,2)*std::pow(z0,2)/std::pow(B,4) - 12*std::pow(C,4)*std::pow(w,2)*z0*zr/std::pow(B,4) + 6*std::pow(C,4)*std::pow(w,2)*std::pow(zr,2)/std::pow(B,4);
  // the 1st order terms
  d = 4*std::pow(A,2)*u*x0 - 4*std::pow(A,2)*u*xr + 4*std::pow(A,2)*v*y0 - 4*std::pow(A,2)*v*yr - 4*std::pow(A,2)*std::pow(C,2)*w*z0/std::pow(B,2) + 4*std::pow(A,2)*std::pow(C,2)*w*zr/std::pow(B,2) + 4*std::pow(C,2)*u*x0 - 4*std::pow(C,2)*u*xr + 4*std::pow(C,2)*v*y0 - 4*std::pow(C,2)*v*yr - 4*u*std::pow(x0,3) + 12*u*std::pow(x0,2)*xr - 12*u*x0*std::pow(xr,2) - 4*u*x0*std::pow(y0,2) + 8*u*x0*y0*yr - 4*u*x0*std::pow(yr,2) + 4*u*std::pow(xr,3) + 4*u*xr*std::pow(y0,2) - 8*u*xr*y0*yr + 4*u*xr*std::pow(yr,2) - 4*v*std::pow(x0,2)*y0 + 4*v*std::pow(x0,2)*yr + 8*v*x0*xr*y0 - 8*v*x0*xr*yr - 4*v*std::pow(xr,2)*y0 + 4*v*std::pow(xr,2)*yr - 4*v*std::pow(y0,3) + 12*v*std::pow(y0,2)*yr - 12*v*y0*std::pow(yr,2) + 4*v*std::pow(yr,3) + 4*std::pow(C,4)*w*z0/std::pow(B,2) - 4*std::pow(C,4)*w*zr/std::pow(B,2) - 4*std::pow(C,2)*u*x0*std::pow(z0,2)/std::pow(B,2) + 8*std::pow(C,2)*u*x0*z0*zr/std::pow(B,2) - 4*std::pow(C,2)*u*x0*std::pow(zr,2)/std::pow(B,2) + 4*std::pow(C,2)*u*xr*std::pow(z0,2)/std::pow(B,2) - 8*std::pow(C,2)*u*xr*z0*zr/std::pow(B,2) + 4*std::pow(C,2)*u*xr*std::pow(zr,2)/std::pow(B,2) - 4*std::pow(C,2)*v*y0*std::pow(z0,2)/std::pow(B,2) + 8*std::pow(C,2)*v*y0*z0*zr/std::pow(B,2) - 4*std::pow(C,2)*v*y0*std::pow(zr,2)/std::pow(B,2) + 4*std::pow(C,2)*v*yr*std::pow(z0,2)/std::pow(B,2) - 8*std::pow(C,2)*v*yr*z0*zr/std::pow(B,2) + 4*std::pow(C,2)*v*yr*std::pow(zr,2)/std::pow(B,2) - 4*std::pow(C,2)*w*std::pow(x0,2)*z0/std::pow(B,2) + 4*std::pow(C,2)*w*std::pow(x0,2)*zr/std::pow(B,2) + 8*std::pow(C,2)*w*x0*xr*z0/std::pow(B,2) - 8*std::pow(C,2)*w*x0*xr*zr/std::pow(B,2) - 4*std::pow(C,2)*w*std::pow(xr,2)*z0/std::pow(B,2) + 4*std::pow(C,2)*w*std::pow(xr,2)*zr/std::pow(B,2) - 4*std::pow(C,2)*w*std::pow(y0,2)*z0/std::pow(B,2) + 4*std::pow(C,2)*w*std::pow(y0,2)*zr/std::pow(B,2) + 8*std::pow(C,2)*w*y0*yr*z0/std::pow(B,2) - 8*std::pow(C,2)*w*y0*yr*zr/std::pow(B,2) - 4*std::pow(C,2)*w*std::pow(yr,2)*z0/std::pow(B,2) + 4*std::pow(C,2)*w*std::pow(yr,2)*zr/std::pow(B,2) - 4*std::pow(C,4)*w*std::pow(z0,3)/std::pow(B,4) + 12*std::pow(C,4)*w*std::pow(z0,2)*zr/std::pow(B,4) - 12*std::pow(C,4)*w*z0*std::pow(zr,2)/std::pow(B,4) + 4*std::pow(C,4)*w*std::pow(zr,3)/std::pow(B,4);
  // the 0th order terms
  e = std::pow(A,4) - 2*std::pow(A,2)*std::pow(C,2) - 2*std::pow(A,2)*std::pow(x0,2) + 4*std::pow(A,2)*x0*xr - 2*std::pow(A,2)*std::pow(xr,2) - 2*std::pow(A,2)*std::pow(y0,2) + 4*std::pow(A,2)*y0*yr - 2*std::pow(A,2)*std::pow(yr,2) + 2*std::pow(A,2)*std::pow(C,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(A,2)*std::pow(C,2)*z0*zr/std::pow(B,2) + 2*std::pow(A,2)*std::pow(C,2)*std::pow(zr,2)/std::pow(B,2) + std::pow(C,4) - 2*std::pow(C,2)*std::pow(x0,2) + 4*std::pow(C,2)*x0*xr - 2*std::pow(C,2)*std::pow(xr,2) - 2*std::pow(C,2)*std::pow(y0,2) + 4*std::pow(C,2)*y0*yr - 2*std::pow(C,2)*std::pow(yr,2) + std::pow(x0,4) - 4*std::pow(x0,3)*xr + 6*std::pow(x0,2)*std::pow(xr,2) + 2*std::pow(x0,2)*std::pow(y0,2) - 4*std::pow(x0,2)*y0*yr + 2*std::pow(x0,2)*std::pow(yr,2) - 4*x0*std::pow(xr,3) - 4*x0*xr*std::pow(y0,2) + 8*x0*xr*y0*yr - 4*x0*xr*std::pow(yr,2) + std::pow(xr,4) + 2*std::pow(xr,2)*std::pow(y0,2) - 4*std::pow(xr,2)*y0*yr + 2*std::pow(xr,2)*std::pow(yr,2) + std::pow(y0,4) - 4*std::pow(y0,3)*yr + 6*std::pow(y0,2)*std::pow(yr,2) - 4*y0*std::pow(yr,3) + std::pow(yr,4) - 2*std::pow(C,4)*std::pow(z0,2)/std::pow(B,2) + 4*std::pow(C,4)*z0*zr/std::pow(B,2) - 2*std::pow(C,4)*std::pow(zr,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(x0,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(x0,2)*z0*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(x0,2)*std::pow(zr,2)/std::pow(B,2) - 4*std::pow(C,2)*x0*xr*std::pow(z0,2)/std::pow(B,2) + 8*std::pow(C,2)*x0*xr*z0*zr/std::pow(B,2) - 4*std::pow(C,2)*x0*xr*std::pow(zr,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(xr,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(xr,2)*z0*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(xr,2)*std::pow(zr,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(y0,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(y0,2)*z0*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(y0,2)*std::pow(zr,2)/std::pow(B,2) - 4*std::pow(C,2)*y0*yr*std::pow(z0,2)/std::pow(B,2) + 8*std::pow(C,2)*y0*yr*z0*zr/std::pow(B,2) - 4*std::pow(C,2)*y0*yr*std::pow(zr,2)/std::pow(B,2) + 2*std::pow(C,2)*std::pow(yr,2)*std::pow(z0,2)/std::pow(B,2) - 4*std::pow(C,2)*std::pow(yr,2)*z0*zr/std::pow(B,2) + 2*std::pow(C,2)*std::pow(yr,2)*std::pow(zr,2)/std::pow(B,2) + std::pow(C,4)*std::pow(z0,4)/std::pow(B,4) - 4*std::pow(C,4)*std::pow(z0,3)*zr/std::pow(B,4) + 6*std::pow(C,4)*std::pow(z0,2)*std::pow(zr,2)/std::pow(B,4) - 4*std::pow(C,4)*z0*std::pow(zr,3)/std::pow(B,4) + std::pow(C,4)*std::pow(zr,4)/std::pow(B,4);

  std::cout << std::scientific;
  std::cout << std::setprecision(9);
  std::cout << a << " " << b << " " << c << " " << d << " " << e << std::endl;
  // solve it
  std::array<double,4> roots;
  auto start = std::chrono::high_resolution_clock::now();
  quartic_solve(a,b,c,d,e,roots);
  auto finish = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration <double, std::micro>(finish-start).count() << " us" << std::endl;
  std::cout << roots[0] << " " << roots[1] << " " << roots[2] << " " << roots[3] << std::endl;
  return 0;
}

Direction
SurfaceZTorus::normal(Position r) const
{
  return {0,0,0};
}

void SurfaceZTorus::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-torus", false);
  std::array<double, 6> coeffs {{x0_, y0_, z0_, A_, B_, C_ }};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================

void read_surfaces(pugi::xml_node node)
{
  // Count the number of surfaces.
  int n_surfaces = 0;
  for (pugi::xml_node surf_node : node.children("surface")) {n_surfaces++;}
  if (n_surfaces == 0) {
    fatal_error("No surfaces found in geometry.xml!");
  }

  // Loop over XML surface elements and populate the array.
  model::surfaces.reserve(n_surfaces);
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node.child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      std::string surf_type = get_node_value(surf_node, "type", true, true);

      if (surf_type == "x-plane") {
        model::surfaces.push_back(std::make_unique<SurfaceXPlane>(surf_node));

      } else if (surf_type == "y-plane") {
        model::surfaces.push_back(std::make_unique<SurfaceYPlane>(surf_node));

      } else if (surf_type == "z-plane") {
        model::surfaces.push_back(std::make_unique<SurfaceZPlane>(surf_node));

      } else if (surf_type == "plane") {
        model::surfaces.push_back(std::make_unique<SurfacePlane>(surf_node));

      } else if (surf_type == "x-cylinder") {
        model::surfaces.push_back(std::make_unique<SurfaceXCylinder>(surf_node));

      } else if (surf_type == "y-cylinder") {
        model::surfaces.push_back(std::make_unique<SurfaceYCylinder>(surf_node));

      } else if (surf_type == "z-cylinder") {
        model::surfaces.push_back(std::make_unique<SurfaceZCylinder>(surf_node));

      } else if (surf_type == "sphere") {
        model::surfaces.push_back(std::make_unique<SurfaceSphere>(surf_node));

      } else if (surf_type == "x-cone") {
        model::surfaces.push_back(std::make_unique<SurfaceXCone>(surf_node));

      } else if (surf_type == "y-cone") {
        model::surfaces.push_back(std::make_unique<SurfaceYCone>(surf_node));

      } else if (surf_type == "z-cone") {
        model::surfaces.push_back(std::make_unique<SurfaceZCone>(surf_node));

      } else if (surf_type == "quadric") {
        model::surfaces.push_back(std::make_unique<SurfaceQuadric>(surf_node));

      } else if (surf_type == "x-torus") {
        model::surfaces.push_back(std::make_unique<SurfaceXTorus>(surf_node));

      } else if (surf_type == "y-torus") {
        model::surfaces.push_back(std::make_unique<SurfaceYTorus>(surf_node));

      } else if (surf_type == "z-torus") {
        model::surfaces.push_back(std::make_unique<SurfaceZTorus>(surf_node));

      } else {
        fatal_error(fmt::format("Invalid surface type, \"{}\"", surf_type));
      }
    }
  }

  // Fill the surface map.
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    int id = model::surfaces[i_surf]->id_;
    auto in_map = model::surface_map.find(id);
    if (in_map == model::surface_map.end()) {
      model::surface_map[id] = i_surf;
    } else {
      fatal_error(fmt::format(
        "Two or more surfaces use the same unique ID: {}", id));
    }
  }

  // Find the global bounding box (of periodic BC surfaces).
  double xmin {INFTY}, xmax {-INFTY}, ymin {INFTY}, ymax {-INFTY},
         zmin {INFTY}, zmax {-INFTY};
  int i_xmin, i_xmax, i_ymin, i_ymax, i_zmin, i_zmax;
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    if (model::surfaces[i_surf]->bc_ == Surface::BoundaryType::PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface* surf_base = model::surfaces[i_surf].get();
      auto surf = dynamic_cast<PeriodicSurface*>(surf_base);

      // Make sure this surface inherits from PeriodicSurface.
      if (!surf) {
        fatal_error(fmt::format(
          "Periodic boundary condition not supported for surface {}. Periodic "
          "BCs are only supported for planar surfaces.", surf_base->id_));
      }

      // See if this surface makes part of the global bounding box.
      auto bb = surf->bounding_box(true) & surf->bounding_box(false);
      if (bb.xmin > -INFTY && bb.xmin < xmin) {
        xmin = bb.xmin;
        i_xmin = i_surf;
      }
      if (bb.xmax < INFTY && bb.xmax > xmax) {
        xmax = bb.xmax;
        i_xmax = i_surf;
      }
      if (bb.ymin > -INFTY && bb.ymin < ymin) {
        ymin = bb.ymin;
        i_ymin = i_surf;
      }
      if (bb.ymax < INFTY && bb.ymax > ymax) {
        ymax = bb.ymax;
        i_ymax = i_surf;
      }
      if (bb.zmin > -INFTY && bb.zmin < zmin) {
        zmin = bb.zmin;
        i_zmin = i_surf;
      }
      if (bb.zmax < INFTY && bb.zmax > zmax) {
        zmax = bb.zmax;
        i_zmax = i_surf;
      }
    }
  }

  // Set i_periodic for periodic BC surfaces.
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    if (model::surfaces[i_surf]->bc_ == Surface::BoundaryType::PERIODIC) {
      // Downcast to the PeriodicSurface type.
      Surface* surf_base = model::surfaces[i_surf].get();
      auto surf = dynamic_cast<PeriodicSurface*>(surf_base);

      // Also try downcasting to the SurfacePlane type (which must be handled
      // differently).
      SurfacePlane* surf_p = dynamic_cast<SurfacePlane*>(surf);

      if (!surf_p) {
        // This is not a SurfacePlane.
        if (surf->i_periodic_ == C_NONE) {
          // The user did not specify the matching periodic surface.  See if we
          // can find the partnered surface from the bounding box information.
          if (i_surf == i_xmin) {
            surf->i_periodic_ = i_xmax;
          } else if (i_surf == i_xmax) {
            surf->i_periodic_ = i_xmin;
          } else if (i_surf == i_ymin) {
            surf->i_periodic_ = i_ymax;
          } else if (i_surf == i_ymax) {
            surf->i_periodic_ = i_ymin;
          } else if (i_surf == i_zmin) {
            surf->i_periodic_ = i_zmax;
          } else if (i_surf == i_zmax) {
            surf->i_periodic_ = i_zmin;
          } else {
            fatal_error("Periodic boundary condition applied to interior "
                        "surface");
          }
        } else {
          // Convert the surface id to an index.
          surf->i_periodic_ = model::surface_map[surf->i_periodic_];
        }
      } else {
        // This is a SurfacePlane.  We won't try to find it's partner if the
        // user didn't specify one.
        if (surf->i_periodic_ == C_NONE) {
          fatal_error(fmt::format("No matching periodic surface specified for "
            "periodic boundary condition on surface {}", surf->id_));
        } else {
          // Convert the surface id to an index.
          surf->i_periodic_ = model::surface_map[surf->i_periodic_];
        }
      }

      // Make sure the opposite surface is also periodic.
      if (model::surfaces[surf->i_periodic_]->bc_ != Surface::BoundaryType::PERIODIC) {
        fatal_error(fmt::format("Could not find matching surface for periodic "
          "boundary condition on surface {}", surf->id_));
      }
    }
  }

  // Check to make sure a boundary condition was applied to at least one
  // surface
  bool boundary_exists = false;
  for (const auto& surf : model::surfaces) {
    if (surf->bc_ != Surface::BoundaryType::TRANSMIT) {
      boundary_exists = true;
      break;
    }
  }
  if (settings::run_mode != RunMode::PLOTTING && !boundary_exists) {
    fatal_error("No boundary conditions were applied to any surfaces!");
  }
}

void free_memory_surfaces()
{
  model::surfaces.clear();
  model::surface_map.clear();
}

} // namespace openmc
