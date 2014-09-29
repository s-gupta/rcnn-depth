/*
 * Triangle with 2D vertex coordinates.
 */
#ifndef MATH__GEOMETRY__TRIANGLE_2D_HH
#define MATH__GEOMETRY__TRIANGLE_2D_HH

#include "io/streams/ostream.hh"
#include "math/geometry/point_2D.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using io::streams::ostream;

/*
 * Triangle with 2D vertex coordinates.
 * Note that vertices are stored by reference.
 */
class triangle_2D {
public:
   /*
    * Constructor.
    */
   explicit triangle_2D(point_2D&, point_2D&, point_2D&);
   
   /*
    * Copy constructor.
    */
   triangle_2D(const triangle_2D&);

   /*
    * Destructor.
    */
   virtual ~triangle_2D();

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<<(ostream&, const triangle_2D&);

   /*
    * Vertex reference.
    */
   point_2D& vertex(unsigned long /* vertex # (0, 1, or 2) */) const;

   /*
    * Angle at vertex.
    */
   double angle(unsigned long /* vertex # (0, 1, or 2) */) const;

   /*
    * Altitude from vertex to opposite side.
    */
   double altitude(unsigned long /* vertex # (0, 1, or 2) */) const;

   /*
    * Length of side opposite vertex.
    */
   double side_length(unsigned long /* vertex # (0, 1, or 2) */) const;
   
   /*
    * Area of triangle.
    */
   double area() const;

   /*
    * Centroid.
    */
   point_2D centroid() const;

protected:
   point_2D& _v0;    /* vertex 0 */
   point_2D& _v1;    /* vertex 1 */
   point_2D& _v2;    /* vertex 2 */
};

} /* namespace geometry */
} /* namespace math */

#endif
