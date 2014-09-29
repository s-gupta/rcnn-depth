/*
 * Point 3D.
 */
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/streams/ostream.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/exact.hh"
#include "math/geometry/point_3D.hh"
#include "math/math.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::streams::ostream;
using math::exact;
using lang::pointers::auto_ptr;

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/
 
/*
 * Default constructor.
 * Return the origin.
 */
point_3D::point_3D()
 : _x(0), _y(0), _z(0)
{ }

/*
 * Constructor.
 * Return the point with the given coordinates.
 */
point_3D::point_3D(const double& x, const double& y, const double& z)
 : _x(x), _y(y), _z(z)
{ }

/*
 * Copy constructor.
 */
point_3D::point_3D(const point_3D& p)
 : _x(p._x), _y(p._y), _z(p._z)
{ }

/*
 * Destructor.
 */
point_3D::~point_3D() {
   /* do nothing */
}

/***************************************************************************
 * Named constructors.
 ***************************************************************************/

/*
 * Polar form.
 * Return the point at the given radius from the origin, theta radians from
 * the positive x-axis, and psi radians above the xy-plane.
 */
point_3D point_3D::polar(
   const double& r, const double& theta, const double& psi)
{
   double r_xy = r*math::cos(psi);
   return point_3D(
      r_xy*math::cos(theta), r_xy*math::sin(theta), r*math::sin(psi)
   );
}

/***************************************************************************
 * Serialization.
 ***************************************************************************/

/*
 * Serialize.
 */
void point_3D::serialize(serial_output_stream& s) const {
   s << _x << _y << _z;
}

/*
 * Deserialize.
 */
auto_ptr<point_3D> point_3D::deserialize(serial_input_stream& s) {
   auto_ptr<point_3D> p(new point_3D());
   s >> p->_x >> p->_y >> p->_z;
   return p;
}

/***************************************************************************
 * I/O.
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
ostream& operator<<(ostream& os, const point_3D& p) {
   os << "(" << p._x << ", " << p._y << ", " << p._z << ")";
   return os;
}

/***************************************************************************
 * Assignment operator.
 ***************************************************************************/

point_3D& point_3D::operator=(const point_3D& p) {
   _x = p._x;
   _y = p._y;
   _z = p._z;
   return *this;
}

/***************************************************************************
 * Binary operators.
 ***************************************************************************/

point_3D operator*(const double& s, const point_3D& p) {
   return point_3D(s*p._x, s*p._y, s*p._z);
}

point_3D operator*(const point_3D& p, const double& s) {
   return point_3D(p._x*s, p._y*s, p._z*s);
}

point_3D operator/(const point_3D& p, const double& s) {
   return point_3D(p._x/s, p._y/s, p._z/s);
}

point_3D operator+(const point_3D& p0, const point_3D& p1) {
   return point_3D(p0._x + p1._x, p0._y + p1._y, p0._z + p1._z);
}

point_3D operator-(const point_3D& p0, const point_3D& p1) {
   return point_3D(p0._x - p1._x, p0._y - p1._y, p0._z - p1._z);
}

/***************************************************************************
 * Unary operators.
 ***************************************************************************/

point_3D point_3D::operator+() const {
   return *this;
}

point_3D point_3D::operator-() const {
   return point_3D(-_x, -_y, -_z);
}

/***************************************************************************
 * Binary equality tests.
 ***************************************************************************/

bool operator==(const point_3D& p0, const point_3D& p1) {
   return ((p0._x == p1._x) && (p0._y == p1._y) && (p0._z == p1._z));
}

bool operator!=(const point_3D& p0, const point_3D& p1) {
   return ((p0._x != p1._x) || (p0._y != p1._y) || (p0._z != p1._z));
}

/***************************************************************************
 * Comparators.
 ***************************************************************************/

/*
 * Binary comparators.
 * These comparators enforce a lexicographic ordering on points.
 */
bool operator<(const point_3D& p0, const point_3D& p1) {
   return (
      (p0._x < p1._x) ||
      ((p0._x == p1._x) &&
         ((p0._y < p1._y) || ((p0._y == p1._y) && (p0._z < p1._z))))
   );
}

bool operator>(const point_3D& p0, const point_3D& p1) {
   return (
      (p0._x > p1._x) ||
      ((p0._x == p1._x) &&
         ((p0._y > p1._y) || ((p0._y == p1._y) && (p0._z > p1._z))))
   );
}

bool operator<=(const point_3D& p0, const point_3D& p1) {
   return (
      (p0._x < p1._x) ||
      ((p0._x == p1._x) &&
         ((p0._y < p1._y) || ((p0._y == p1._y) && (p0._z <= p1._z))))
   );
}

bool operator>=(const point_3D& p0, const point_3D& p1) {
   return (
      (p0._x > p1._x) ||
      ((p0._x == p1._x) &&
         ((p0._y > p1._y) || ((p0._y == p1._y) && (p0._z >= p1._z))))
   );
}

/*
 * Comparison method (for lexicographic ordering on points).
 */
int point_3D::compare_to(const point_3D& p) const {
   int cmp_x = (_x < p._x) ? -1 : ((_x > p._x) ? 1 : 0);
   int cmp_y = (_y < p._y) ? -1 : ((_y > p._y) ? 1 : 0);
   int cmp_z = (_z < p._z) ? -1 : ((_z > p._z) ? 1 : 0);
   return ((cmp_x == 0) ? ((cmp_y == 0) ? cmp_z : cmp_y) : cmp_x);
}

/***************************************************************************
 * Polar coordinates.
 ***************************************************************************/

/*
 * Magnitude (distance from origin).
 */
double abs(const point_3D& p) {
   return math::sqrt(p._x*p._x + p._y*p._y + p._z*p._z);
}

/*
 * Argument (angle in radians from positive x-axis).
 */
double arg(const point_3D& p) {
   return math::atan2(p._y, p._x);
}

/*
 * Elevation (angle in radians from the xy-plane).
 */
double elev(const point_3D& p) {
   double r_xy = math::sqrt(p._x*p._x + p._y*p._y);
   return math::atan2(p._z, r_xy);
}

/***************************************************************************
 * Dot product.
 ***************************************************************************/

/*
 * Dot product.
 */
double dot(const point_3D& p0, const point_3D& p1) {
   return (p0._x * p1._x + p0._y * p1._y + p0._z * p1._z);
}

/***************************************************************************
 * Distance from point to line/segment.
 ***************************************************************************/

/*
 * Distance from the point to the line passing through the specified points.
 */
double point_3D::distance_to_line(
   const point_3D& p, const point_3D& q) const
{
   /* compute vectors between points */
   point_3D a = q - p;
   point_3D b = *this - p;
   /* compute distance */
   double mag_a_sq = a._x*a._x + a._y*a._y + a._z*a._z;
   double mag_b_sq = b._x*b._x + b._y*b._y + b._z*b._z;
   double a_dot_b = dot(a, b);
   return (mag_a_sq*mag_b_sq - a_dot_b*a_dot_b)/mag_a_sq;
}

/*
 * Distance from the point to the line segment with the specified endpoints.
 */
double point_3D::distance_to_segment(
   const point_3D& p, const point_3D& q) const
{
   /* compute vectors between points */
   point_3D a = q - p;
   point_3D b = *this - p;
   point_3D c = *this - q;
   /* compute distance to segment */
   if (dot(a, c) >= 0) {
      /* q is closest point */
      return abs(c);
   } else if (dot(-a, b) >= 0) {
      /* p is closest point */
      return abs(b);
   } else {
      /* compute distance from point to line */
      return this->distance_to_line(p, q);
   }
}

/***************************************************************************
 * Geometric predicates. 
 ***************************************************************************/

/*
 * Robust geometric orientation test (using exact arithmetic).
 * Return -1 if the point is above plane pqr,
 *         0 if the points are coplanar, and
 *         1 if the point is below plane pqr.
 *
 * Note that "below" and "above" are defined so that p, q, and r appear
 * in counterclockwise order when viewed from above plane pqr.
 */
int point_3D::orientation(
   const point_3D& p, const point_3D& q, const point_3D& r) const
{
   /* compute error bound beyond which exact arithmetic is needed */
   static const double error_bound = 
      (7.0 + 56.0 * exact<>::epsilon()) * exact<>::epsilon();
   /* perform approximate orientation test */
   double ax = p.x() - _x;  double ay = p.y() - _y;  double az = p.z() - _z;
   double bx = q.x() - _x;  double by = q.y() - _y;  double bz = q.z() - _z;
   double cx = r.x() - _x;  double cy = r.y() - _y;  double cz = r.z() - _z;
   double ax_by = ax * by;  double bx_ay = bx * ay;
   double ax_cy = ax * cy;  double cx_ay = cx * ay;
   double bx_cy = bx * cy;  double cx_by = cx * by;
   double det =
      az * (bx_cy - cx_by) + bz * (cx_ay - ax_cy) + cz * (ax_by - bx_ay);
   /* check if error cannot change result */
   double permanent = 
      math::abs(az) * (math::abs(bx_cy) + math::abs(cx_by)) +
      math::abs(bz) * (math::abs(cx_ay) + math::abs(ax_cy)) +
      math::abs(cz) * (math::abs(ax_by) + math::abs(bx_ay));
   double max_error = error_bound * permanent;
   if ((-det) < max_error) {
      return -1;
   } else if (det > max_error) {
      return 1;
   } else {
      /* perform exact orientation test */
      exact<> ax_e = exact<>(p.x()) - _x;
      exact<> ay_e = exact<>(p.y()) - _y;
      exact<> az_e = exact<>(p.z()) - _z;
      exact<> bx_e = exact<>(q.x()) - _x;
      exact<> by_e = exact<>(q.y()) - _y;
      exact<> bz_e = exact<>(q.z()) - _z;
      exact<> cx_e = exact<>(r.x()) - _x;
      exact<> cy_e = exact<>(r.y()) - _y;
      exact<> cz_e = exact<>(r.z()) - _z;
      exact<> ax_by_e = ax_e * by_e;  exact<> bx_ay_e = bx_e * ay_e;
      exact<> ax_cy_e = ax_e * cy_e;  exact<> cx_ay_e = cx_e * ay_e;
      exact<> bx_cy_e = bx_e * cy_e;  exact<> cx_by_e = cx_e * by_e;
      exact<> det_exact = 
         az_e * (bx_cy_e - cx_by_e) +
         bz_e * (cx_ay_e - ax_cy_e) +
         cz_e * (ax_by_e - bx_ay_e);
      return det_exact.sign();
   }
}

/*
 * Robust geometric in-sphere test (using exact arithmetic).
 * Return -1 if the point is outside the sphere,
 *         0 if the point is on the sphere, and 
 *         1 if the point is inside the sphere.
 *
 * The four points defining the sphere must appear in order so that the 
 * fourth point is below the plane defined by the first three (positive 
 * orientation), or the sign of the result will be reversed.
 */
int point_3D::in_sphere(
   const point_3D& p0,
   const point_3D& p1,
   const point_3D& p2,
   const point_3D& p3) const
{
   /* compute error bound beyond which exact arithmetic is needed */
   static const double error_bound =
      (16.0 + 224.0 * exact<>::epsilon()) * exact<>::epsilon();
   /* perform approximate in-sphere test */
   double ax = p0.x() - _x;  double ay = p0.y() - _y;  double az = p0.z() - _z;
   double bx = p1.x() - _x;  double by = p1.y() - _y;  double bz = p1.z() - _z;
   double cx = p2.x() - _x;  double cy = p2.y() - _y;  double cz = p2.z() - _z;
   double dx = p3.x() - _x;  double dy = p3.y() - _y;  double dz = p3.z() - _z;
   double ax_by = ax * by;  double bx_ay = bx * ay;
   double bx_cy = bx * cy;  double cx_by = cx * by;
   double cx_dy = cx * dy;  double dx_cy = dx * cy;
   double dx_ay = dx * ay;  double ax_dy = ax * dy;
   double ax_cy = ax * cy;  double cx_ay = cx * ay;
   double bx_dy = bx * dy;  double dx_by = dx * by;
   double ab_diff = ax_by - bx_ay;
   double bc_diff = bx_cy - cx_by;
   double cd_diff = cx_dy - dx_cy;
   double da_diff = dx_ay - ax_dy;
   double ac_diff = ax_cy - cx_ay;
   double bd_diff = bx_dy - dx_by;
   double abc = az * bc_diff - bz * ac_diff + cz * ab_diff;
   double bcd = bz * cd_diff - cz * bd_diff + dz * bc_diff;
   double cda = cz * da_diff + dz * ac_diff + az * cd_diff;
   double dab = dz * ab_diff + az * bd_diff + bz * da_diff;
   double a_sq = ax * ax + ay * ay + az * az;
   double b_sq = bx * bx + by * by + bz * bz;
   double c_sq = cx * cx + cy * cy + cz * cz;
   double d_sq = dx * dx + dy * dy + dz * dz;
   double det = (d_sq * abc - c_sq * dab) + (b_sq * cda - a_sq * bcd);
   /* check if error cannot change result */
   double permanent = 
      ((math::abs(cx_dy) + math::abs(dx_cy)) * math::abs(bz) +
       (math::abs(dx_by) + math::abs(bx_dy)) * math::abs(cz) +
       (math::abs(bx_cy) + math::abs(cx_by)) * math::abs(dz)) * a_sq +
      ((math::abs(dx_ay) + math::abs(ax_dy)) * math::abs(cz) +
       (math::abs(ax_cy) + math::abs(cx_ay)) * math::abs(dz) +
       (math::abs(cx_dy) + math::abs(dx_cy)) * math::abs(az)) * b_sq + 
      ((math::abs(ax_by) + math::abs(bx_ay)) * math::abs(dz) +
       (math::abs(bx_dy) + math::abs(dx_by)) * math::abs(az) +
       (math::abs(dx_ay) + math::abs(ax_dy)) * math::abs(bz)) * c_sq +
      ((math::abs(bx_cy) + math::abs(cx_by)) * math::abs(az) +
       (math::abs(cx_ay) + math::abs(ax_cy)) * math::abs(bz) +
       (math::abs(ax_by) + math::abs(bx_ay)) * math::abs(cz)) * d_sq;
   double max_error = error_bound * permanent;
   if ((-det) < max_error) {
      return -1;
   } else if (det > max_error) {
      return 1;
   } else {
      /* perform exact in-sphere test */ 
      exact<> ax_e = exact<>(p0.x()) - _x;
      exact<> ay_e = exact<>(p0.y()) - _y;
      exact<> az_e = exact<>(p0.z()) - _z;
      exact<> bx_e = exact<>(p1.x()) - _x;
      exact<> by_e = exact<>(p1.y()) - _y;
      exact<> bz_e = exact<>(p1.z()) - _z;
      exact<> cx_e = exact<>(p2.x()) - _x;
      exact<> cy_e = exact<>(p2.y()) - _y;
      exact<> cz_e = exact<>(p2.z()) - _z;
      exact<> dx_e = exact<>(p3.x()) - _x;
      exact<> dy_e = exact<>(p3.y()) - _y;
      exact<> dz_e = exact<>(p3.z()) - _z;
      exact<> ab_diff_e = ax_e * by_e - bx_e * ay_e;
      exact<> bc_diff_e = bx_e * cy_e - cx_e * by_e;
      exact<> cd_diff_e = cx_e * dy_e - dx_e * cy_e;
      exact<> da_diff_e = dx_e * ay_e - ax_e * dy_e;
      exact<> ac_diff_e = ax_e * cy_e - cx_e * ay_e;
      exact<> bd_diff_e = bx_e * dy_e - dx_e * by_e;
      exact<> abc_e = az_e * bc_diff_e - bz_e * ac_diff_e + cz_e * ab_diff_e;
      exact<> bcd_e = bz_e * cd_diff_e - cz_e * bd_diff_e + dz_e * bc_diff_e;
      exact<> cda_e = cz_e * da_diff_e + dz_e * ac_diff_e + az_e * cd_diff_e;
      exact<> dab_e = dz_e * ab_diff_e + az_e * bd_diff_e + bz_e * da_diff_e; 
      exact<> a_sq_e = ax_e * ax_e + ay_e * ay_e + az_e * az_e;
      exact<> b_sq_e = bx_e * bx_e + by_e * by_e + bz_e * bz_e;
      exact<> c_sq_e = cx_e * cx_e + cy_e * cy_e + cz_e * cz_e;
      exact<> d_sq_e = dx_e * dx_e + dy_e * dy_e + dz_e * dz_e;
      exact<> det_exact = 
         (d_sq_e * abc_e - c_sq_e * dab_e) + (b_sq_e * cda_e - a_sq_e * bcd_e);
      return det_exact.sign();
   }
}

} /* namespace geometry */
} /* namespace math */
