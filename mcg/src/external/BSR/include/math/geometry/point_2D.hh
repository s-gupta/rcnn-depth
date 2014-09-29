/*
 * Point 2D.
 */
#ifndef MATH__GEOMETRY__POINT_2D_HH
#define MATH__GEOMETRY__POINT_2D_HH

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
 * Declare existence of triangle_2D class.
 */
class triangle_2D;

/*
 * Point 2D.
 */
class point_2D : public comparable<point_2D> {
public:
   /*
    * Constructors.
    */
   point_2D();
   explicit point_2D(const double& /* x */, const double& /* y */);

   /*
    * Copy constructor.
    */
   point_2D(const point_2D&);

   /*
    * Destructor.
    */
   virtual ~point_2D();

   /*
    * Polar form.
    */
   static point_2D polar(const double& /* r */, const double& /* theta */);

   /*
    * Coordinate reference.
    */
   inline double& x() { return _x; }
   inline double& y() { return _y; }
   inline const double& x() const { return _x; }
   inline const double& y() const { return _y; }

   /*
    * Serialize.
    */
   void serialize(serial_output_stream&) const;

   /*
    * Deserialize.
    */
   static auto_ptr<point_2D> deserialize(serial_input_stream&);

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<<(ostream&, const point_2D&);

   /*
    * Assignment operator.
    */
   point_2D& operator=(const point_2D&);

   /*
    * Binary operators.
    */
   friend point_2D operator*(const double&, const point_2D&);
   friend point_2D operator*(const point_2D&, const double&);
   friend point_2D operator/(const point_2D&, const double&);

   friend point_2D operator+(const point_2D&, const point_2D&);
   friend point_2D operator-(const point_2D&, const point_2D&);

   /*
    * Unary operators.
    */
   point_2D operator+() const;
   point_2D operator-() const;

   /*
    * Binary equality tests.
    */
   friend bool operator==(const point_2D&, const point_2D&);
   friend bool operator!=(const point_2D&, const point_2D&);

   /*
    * Binary comparators.
    * These comparators enforce a lexicographic ordering on points.
    */
   friend bool operator< (const point_2D&, const point_2D&);
   friend bool operator> (const point_2D&, const point_2D&);
   friend bool operator<=(const point_2D&, const point_2D&);
   friend bool operator>=(const point_2D&, const point_2D&);

   /*
    * Comparison method (for lexicographic ordering on points).
    */
   int compare_to(const point_2D&) const;
   
   /*
    * Magnitude (distance from origin).
    */
   friend double abs(const point_2D&);

   /*
    * Argument (angle in radians from positive x-axis).
    */
   friend double arg(const point_2D&);

   /*
    * Dot product.
    */
   friend double dot(const point_2D&, const point_2D&);

   /*
    * Distance from the point to the line passing through the specified points.
    */
   double distance_to_line(const point_2D&, const point_2D&) const;

   /*
    * Distance from the point to the line segment with the specified endpoints.
    */
   double distance_to_segment(const point_2D&, const point_2D&) const;
    
   /*
    * Robust geometric orientation test (using exact arithmetic).
    * Return -1 if the point is to the right of ray pq (clockwise),
    *         0 if the point is on the line pq, and
    *         1 if the point is to the left of ray pq (counterclockwise).
    */
   int orientation(const point_2D& /* p */, const point_2D& /* q */) const;

   /*
    * Robust geometric in-circle test (using exact arithmetic).
    * Return -1 if the point is outside the circle,
    *         0 if the point is on the circle, and 
    *         1 if the point is inside the circle.
    *
    * The three points defining the circle must appear in counterclockwise 
    * order, or the sign of the result will be reversed.
    */
   int in_circle(const point_2D&, const point_2D&, const point_2D&) const;
   
   /*
    * Robust geometric in-triangle test (using exact arithmetic).
    * Return -1 if the point is outside the triangle,
    *         0 if the point is on the triangle boundary, and
    *         1 if the point is inside the triangle.
    *
    * The order of the vertices defining the triangle does not matter.
    */
   int in_triangle(const triangle_2D&) const;
   
   int in_triangle(const point_2D&, const point_2D&, const point_2D&) const;

protected:
   double _x;  /* x-coordinate */
   double _y;  /* y-coordinate */
};

} /* namespace geometry */
} /* namespace math */

#endif
