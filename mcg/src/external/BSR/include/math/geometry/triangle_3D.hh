/*
 * Triangle with 3D vertex coordinates.
 */
#ifndef MATH__GEOMETRY__TRIANGLE_3D_HH
#define MATH__GEOMETRY__TRIANGLE_3D_HH

#include "io/streams/ostream.hh"
#include "math/geometry/point_3D.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using io::streams::ostream;

/*
 * Triangle with 3D vertex coordinates.
 * Note that vertices are stored by reference.
 */
class triangle_3D {
public:
   /*
    * Constructor.
    */
   explicit triangle_3D(point_3D&, point_3D&, point_3D&);
   
   /*
    * Copy constructor.
    */
   triangle_3D(const triangle_3D&);

   /*
    * Destructor.
    */
   virtual ~triangle_3D();

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<<(ostream&, const triangle_3D&);

   /*
    * Vertex reference.
    */
   point_3D& vertex(unsigned long /* vertex # (0, 1, or 2) */) const;

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
   point_3D centroid() const;

protected:
   point_3D& _v0;    /* vertex 0 */
   point_3D& _v1;    /* vertex 1 */
   point_3D& _v2;    /* vertex 2 */
};

} /* namespace geometry */
} /* namespace math */

#endif
