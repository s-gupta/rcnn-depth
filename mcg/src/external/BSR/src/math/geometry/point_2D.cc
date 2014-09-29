/*
 * Point 2D.
 */
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/streams/ostream.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/exact.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/triangle_2D.hh"
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
point_2D::point_2D()
 : _x(0), _y(0)
{ }

/*
 * Constructor.
 * Return the point with the given coordinates.
 */
point_2D::point_2D(const double& x, const double& y)
 : _x(x), _y(y)
{ }

/*
 * Copy constructor.
 */
point_2D::point_2D(const point_2D& p)
 : _x(p._x), _y(p._y)
{ }

/*
 * Destructor.
 */
point_2D::~point_2D() {
   /* do nothing */
}

/***************************************************************************
 * Named constructors.
 ***************************************************************************/

/*
 * Polar form.
 */
point_2D point_2D::polar(const double& r, const double& theta) {
   return point_2D(r*math::cos(theta), r*math::sin(theta));
}

/***************************************************************************
 * Serialization.
 ***************************************************************************/

/*
 * Serialize.
 */
void point_2D::serialize(serial_output_stream& s) const {
   s << _x << _y;
}

/*
 * Deserialize.
 */
auto_ptr<point_2D> point_2D::deserialize(serial_input_stream& s) {
   auto_ptr<point_2D> p(new point_2D());
   s >> p->_x >> p->_y;
   return p;
}

/***************************************************************************
 * I/O.
 ***************************************************************************/

/*
 * Formatted output to stream.
 */
ostream& operator<<(ostream& os, const point_2D& p) {
   os << "(" << p._x << ", " << p._y << ")";
   return os;
}

/***************************************************************************
 * Assignment operator.
 ***************************************************************************/

point_2D& point_2D::operator=(const point_2D& p) {
   _x = p._x;
   _y = p._y;
   return *this;
}

/***************************************************************************
 * Binary operators.
 ***************************************************************************/

point_2D operator*(const double& s, const point_2D& p) {
   return point_2D(s*p._x, s*p._y);
}

point_2D operator*(const point_2D& p, const double& s) {
   return point_2D(p._x*s, p._y*s);
}

point_2D operator/(const point_2D& p, const double& s) {
   return point_2D(p._x/s, p._y/s);
}

point_2D operator+(const point_2D& p0, const point_2D& p1) {
   return point_2D(p0._x + p1._x, p0._y + p1._y);
}

point_2D operator-(const point_2D& p0, const point_2D& p1) {
   return point_2D(p0._x - p1._x, p0._y - p1._y);
}

/***************************************************************************
 * Unary operators.
 ***************************************************************************/

point_2D point_2D::operator+() const {
   return *this;
}

point_2D point_2D::operator-() const {
   return point_2D(-_x, -_y);
}

/***************************************************************************
 * Binary equality tests.
 ***************************************************************************/

bool operator==(const point_2D& p0, const point_2D& p1) {
   return ((p0._x == p1._x) && (p0._y == p1._y));
}

bool operator!=(const point_2D& p0, const point_2D& p1) {
   return ((p0._x != p1._x) || (p0._y != p1._y));
}

/***************************************************************************
 * Comparators.
 ***************************************************************************/

/*
 * Binary comparators.
 * These comparators enforce a lexicographic ordering on points.
 */
bool operator<(const point_2D& p0, const point_2D& p1) {
   return ((p0._x < p1._x) || ((p0._x == p1._x) && (p0._y < p1._y)));
}

bool operator>(const point_2D& p0, const point_2D& p1) {
   return ((p0._x > p1._x) || ((p0._x == p1._x) && (p0._y > p1._y)));
}

bool operator<=(const point_2D& p0, const point_2D& p1) {
   return ((p0._x < p1._x) || ((p0._x == p1._x) && (p0._y <= p1._y)));
}

bool operator>=(const point_2D& p0, const point_2D& p1) {
   return ((p0._x > p1._x) || ((p0._x == p1._x) && (p0._y >= p1._y)));
}

/*
 * Comparison method (for lexicographic ordering on points).
 */
int point_2D::compare_to(const point_2D& p) const {
   int cmp_x = (_x < p._x) ? -1 : ((_x > p._x) ? 1 : 0);
   int cmp_y = (_y < p._y) ? -1 : ((_y > p._y) ? 1 : 0);
   return ((cmp_x == 0) ? cmp_y : cmp_x);
}

/***************************************************************************
 * Polar coordinates.
 ***************************************************************************/

/*
 * Magnitude (distance from origin).
 */
double abs(const point_2D& p) {
   return math::sqrt(p._x*p._x + p._y*p._y);
}

/*
 * Argument (angle in radians from positive x-axis).
 */
double arg(const point_2D& p) {
   return math::atan2(p._y, p._x);
}

/***************************************************************************
 * Dot product.
 ***************************************************************************/

/*
 * Dot product.
 */
double dot(const point_2D& p0, const point_2D& p1) {
   return (p0._x * p1._x + p0._y * p1._y);
}

/***************************************************************************
 * Distance from point to line/segment.
 ***************************************************************************/

/*
 * Distance from the point to the line passing through the specified points.
 */
double point_2D::distance_to_line(
   const point_2D& p, const point_2D& q) const
{
   double xa = q._x - p._x;  double ya = q._y - p._y;
   double xb = _x   - p._x;  double yb = _y   - p._y;
   return (math::abs(xa*yb - ya*xb))/(math::sqrt(xa*xa + ya*ya));
}

/*
 * Distance from the point to the line segment with the specified endpoints.
 */
double point_2D::distance_to_segment(
   const point_2D& p, const point_2D& q) const
{
   /* compute vectors between points */
   double xa = q._x - p._x;  double ya = q._y - p._y;
   double xb = _x   - p._x;  double yb = _y   - p._y;
   double xc = _x   - q._x;  double yc = _y   - q._y;
   /* compute distance to segment */
   if ((xa*xc + ya*yc) >= 0) {   
      /* q is the closest point */
      return math::sqrt(xc*xc + yc*yc);
   } else if ((-xa*xb - ya*yb) >= 0) {
      /* p is the closest point */
      return math::sqrt(xb*xb + yb*yb);
   } else {
      /* compute distance from point to line */
      return (math::abs(xa*yb - ya*xb))/(math::sqrt(xa*xa + ya*ya));
   }
}

/***************************************************************************
 * Geometric predicates. 
 ***************************************************************************/

/*
 * Robust geometric orientation test (using exact arithmetic).
 * Return -1 if the point is to the right of ray pq (clockwise),
 *         0 if the point is on the line pq, and
 *         1 if the point is to the left of ray pq (counterclockwise).
 */
int point_2D::orientation(const point_2D& p, const point_2D& q) const {
   /* compute error bound beyond which exact arithmetic is needed */
   static const double error_bound = 
      (3.0 + 16.0 * exact<>::epsilon()) * exact<>::epsilon();
   /* perform approximate orientation test */
   double det_left  = (p.x() - _x) * (q.y() - _y);
   double det_right = (q.x() - _x) * (p.y() - _y);
   /* return if result is certain */
   double det_sum = 0;
   if (det_left < 0) {
      if (det_right >= 0)
         return -1;
      else
         det_sum = -det_left - det_right;
   } else if (det_left > 0) {
      if (det_right <= 0)
         return 1;
      else
         det_sum = det_left + det_right;
   } else {
      return (det_right > 0) ? -1 : ((det_right < 0) ? 1 : 0);
   }
   /* check if error cannot change result */
   double max_error = error_bound * det_sum;
   double det = det_left - det_right;
   if ((-det) >= max_error) {
      return -1;
   } else if (det >= max_error) {
      return 1;
   } else {
      /* perform exact orientation test */
      exact<> det_exact = 
         (exact<>(p.x()) - _x) * (exact<>(q.y()) - _y) - 
         (exact<>(q.x()) - _x) * (exact<>(p.y()) - _y);
      return det_exact.sign();
   }
}

/*
 * Robust geometric in-circle test (using exact arithmetic).
 * Return -1 if the point is outside the circle,
 *         0 if the point is on the circle, and 
 *         1 if the point is inside the circle.
 *
 * The three points defining the circle must appear in counterclockwise 
 * order, or the sign of the result will be reversed.
 */
int point_2D::in_circle(
   const point_2D& p0, const point_2D& p1, const point_2D& p2) const
{
   /* compute error bound beyond which exact arithmetic is needed */
   static const double error_bound =
      (10.0 + 96.0 * exact<>::epsilon()) * exact<>::epsilon();
   /* perform approximate in-circle test */
   double ax = p0.x() - _x;  double ay = p0.y() - _y;
   double bx = p1.x() - _x;  double by = p1.y() - _y;
   double cx = p2.x() - _x;  double cy = p2.y() - _y;
   double a_sq = ax * ax + ay * ay;
   double b_sq = bx * bx + by * by;
   double c_sq = cx * cx + cy * cy;
   double ax_by = ax * by;  double bx_ay = bx * ay;
   double ax_cy = ax * cy;  double cx_ay = cx * ay;
   double bx_cy = bx * cy;  double cx_by = cx * by;
   double det =
      a_sq * (bx_cy - cx_by) + b_sq * (cx_ay - ax_cy) + c_sq * (ax_by - bx_ay);
   /* check if error cannot change result */
   double permanent = 
      a_sq * (math::abs(bx_cy) + math::abs(cx_by)) + 
      b_sq * (math::abs(cx_ay) + math::abs(ax_cy)) +
      c_sq * (math::abs(ax_by) + math::abs(bx_ay));
   double max_error = error_bound * permanent;
   if ((-det) > max_error) {
      return -1;
   } else if (det > max_error) {
      return 1;
   } else {
      /* perform exact in-circle test */
      exact<> p0x(p0.x());  exact<> p0y(p0.y());
      exact<> p1x(p1.x());  exact<> p1y(p1.y());
      exact<> p2x(p2.x());  exact<> p2y(p2.y());
      exact<> ax_e = p0x - _x;  exact<> ay_e = p0y - _y;
      exact<> bx_e = p1x - _x;  exact<> by_e = p1y - _y;
      exact<> cx_e = p2x - _x;  exact<> cy_e = p2y - _y;
      exact<> a_sq_e = ax_e * ax_e + ay_e * ay_e;
      exact<> b_sq_e = bx_e * bx_e + by_e * by_e;
      exact<> c_sq_e = cx_e * cx_e + cy_e * cy_e;
      exact<> ax_by_e = ax_e * by_e;  exact<> bx_ay_e = bx_e * ay_e;
      exact<> ax_cy_e = ax_e * cy_e;  exact<> cx_ay_e = cx_e * ay_e;
      exact<> bx_cy_e = bx_e * cy_e;  exact<> cx_by_e = cx_e * by_e;
      exact<> det_exact =
         a_sq_e * (bx_cy_e - cx_by_e) + 
         b_sq_e * (cx_ay_e - ax_cy_e) +
         c_sq_e * (ax_by_e - bx_ay_e);
      return det_exact.sign();
   }
}

/*
 * Robust geometric in-triangle test (using exact arithmetic).
 * Return -1 if the point is outside the triangle,
 *         0 if the point is on the triangle boundary, and
 *         1 if the point is inside the triangle.
 *
 * The order of the vertices defining the triangle does not matter.
 */
int point_2D::in_triangle(const triangle_2D& t) const {
   return this->in_triangle(t.vertex(0), t.vertex(1), t.vertex(2));
}

int point_2D::in_triangle(
   const point_2D& p0, const point_2D& p1, const point_2D& p2) const
{
   /* compute initial orientation factor */
   int orient = p2.orientation(p0, p1);
   /* handle degenerate triangles */
   if (orient == 0) {
      const point_2D& min01 = (p0 < p1) ? p0 : p1;
      const point_2D& max01 = (p0 > p1) ? p0 : p1;
      const point_2D& p_min = (min01 < p2) ? min01 : p2;
      const point_2D& p_max = (max01 > p2) ? max01 : p2;
      return (
         (*this >= p_min) &&
         (*this <= p_max) &&
         (this->orientation(p_min, p_max) == 0)
      ) ? 0 : -1;
   }  
   /* perform orientation test for each side */
   int cmp0 = orient * this->orientation(p0, p1); if (cmp0 < 0) { return -1; }
   int cmp1 = orient * this->orientation(p1, p2); if (cmp1 < 0) { return -1; }
   int cmp2 = orient * this->orientation(p2, p0); if (cmp2 < 0) { return -1; }
   return ((cmp0 == 0) || (cmp1 == 0) || (cmp2 == 0)) ? 0 : 1;
}

} /* namespace geometry */
} /* namespace math */
