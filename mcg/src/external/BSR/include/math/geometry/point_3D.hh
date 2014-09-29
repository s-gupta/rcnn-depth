/*
 * Point 3D.
 */
#ifndef MATH__GEOMETRY__POINT_3D_HH
#define MATH__GEOMETRY__POINT_3D_HH

#include "interfaces/comparable.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/streams/ostream.hh"
#include "lang/pointers/auto_ptr.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using interfaces::comparable;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::streams::ostream;
using lang::pointers::auto_ptr;

/*
 * Point 3D.
 */
class point_3D : public comparable<point_3D> {
public:
   /*
    * Constructors.
    */
   point_3D();
   explicit point_3D(
      const double& /* x */, const double& /* y */, const double& /* z */
   );

   /*
    * Copy constructor.
    */
   point_3D(const point_3D&);

   /*
    * Destructor.
    */
   virtual ~point_3D();

   /*
    * Polar form.
    * Return the point at the given radius from the origin, theta radians from
    * the positive x-axis, and psi radians above the xy-plane.
    */
   static point_3D polar(
      const double& /* r */, const double& /* theta */, const double& /* psi */
   );

   /*
    * Coordinate reference.
    */
   inline double& x() { return _x; }
   inline double& y() { return _y; }
   inline double& z() { return _z; }
   inline const double& x() const { return _x; }
   inline const double& y() const { return _y; }
   inline const double& z() const { return _z; }

   /*
    * Serialize.
    */
   void serialize(serial_output_stream&) const;

   /*
    * Deserialize.
    */
   static auto_ptr<point_3D> deserialize(serial_input_stream&);

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<<(ostream&, const point_3D&);

   /*
    * Assignment operator.
    */
   point_3D& operator=(const point_3D&);

   /*
    * Binary operators.
    */
   friend point_3D operator*(const double&, const point_3D&);
   friend point_3D operator*(const point_3D&, const double&);
   friend point_3D operator/(const point_3D&, const double&);

   friend point_3D operator+(const point_3D&, const point_3D&);
   friend point_3D operator-(const point_3D&, const point_3D&);

   /*
    * Unary operators.
    */
   point_3D operator+() const;
   point_3D operator-() const;

   /*
    * Binary equality tests.
    */
   friend bool operator==(const point_3D&, const point_3D&);
   friend bool operator!=(const point_3D&, const point_3D&);

   /*
    * Binary comparators.
    * These comparators enforce a lexicographic ordering on points.
    */
   friend bool operator< (const point_3D&, const point_3D&);
   friend bool operator> (const point_3D&, const point_3D&);
   friend bool operator<=(const point_3D&, const point_3D&);
   friend bool operator>=(const point_3D&, const point_3D&);

   /*
    * Comparison method (for lexicographic ordering on points).
    */
   int compare_to(const point_3D&) const;
    
   /*
    * Magnitude (distance from origin).
    */
   friend double abs(const point_3D&);

   /*
    * Argument (angle in radians from positive x-axis).
    */
   friend double arg(const point_3D&);

   /*
    * Elevation (angle in radians from the xy-plane).
    */
   friend double elev(const point_3D&);

   /*
    * Dot product.
    */
   friend double dot(const point_3D&, const point_3D&);

   /*
    * Distance from the point to the line passing through the specified points.
    */
   double distance_to_line(const point_3D&, const point_3D&) const;

   /*
    * Distance from the point to the line segment with the specified endpoints.
    */
   double distance_to_segment(const point_3D&, const point_3D&) const;
   
   /*
    * Robust geometric orientation test (using exact arithmetic).
    * Return -1 if the point is above plane pqr,
    *         0 if the points are coplanar, and
    *         1 if the point is below plane pqr.
    *
    * Note that "below" and "above" are defined so that p, q, and r appear
    * in counterclockwise order when viewed from above plane pqr.
    */
   int orientation(
      const point_3D& /* p */, const point_3D& /* q */, const point_3D& /* r */
   ) const;

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
   int in_sphere(
      const point_3D&, const point_3D&, const point_3D&, const point_3D&
   ) const;
   
protected:
   double _x;  /* x-coordinate */
   double _y;  /* y-coordinate */
   double _z;  /* z-coordinate */
};

} /* namespace geometry */
} /* namespace math */

#endif
