/*
 * Boundary element.
 */
#ifndef VISION__SEGMENTATION__BOUNDARY_ELEMENT_HH
#define VISION__SEGMENTATION__BOUNDARY_ELEMENT_HH

#include "collections/array_list.hh"
#include "collections/pointers/auto_collection.hh"
#include "lang/array.hh"
#include "math/libraries/lib_image.hh"
#include "math/matrices/matrix.hh"

namespace vision {
namespace segmentation {
/*
 * Imports.
 */
using collections::array_list;
using collections::pointers::auto_collection;
using lang::array;
using math::libraries::lib_image;
using math::matrices::matrix;

/*
 * Boundary element.
 */
class boundary_element {
public:
   /*
    * Create boundary elements for the given collection of contours.
    */
   static auto_collection< boundary_element, array_list<boundary_element> > 
      create_boundary_elements(
         const matrix<unsigned long>&,                 /* contour assignments */
         const array_list<lib_image::contour_vertex>&, /* contour vertices */
         const array_list<lib_image::contour_edge>&    /* contour edges */
      );

   /*
    * Copy constructor.
    */
   boundary_element(const boundary_element&);

   /*
    * Destructor.
    */
   virtual ~boundary_element();

   /*
    * Get the id of the smooth contour represented by the boundary element.
    */
   unsigned long contour_id() const;

   /*
    * Get the id of the contour set to which the boundary element belongs.
    *
    * Each contour set is a longest possible sequence of connected contours 
    * such that each connection vertex joins two constraint boundaries and 
    * all other boundaries intersecting it are completion edges.
    *
    * Contour sets correspond to the image contours that existed prior to
    * subdivision into approximate straight line segments.  Each completion
    * edge is in its own separate contour set.
    */
   unsigned long contour_set_id() const;
   
   /*
    * Check whether the boundary element is a completion edge.
    */
   bool is_completion() const;
   
   /*
    * Check whether the boundary element is an intruding edge.
    * An intruding edge is a completion edge which interrupts a contour set.
    */
   bool is_intruder() const;
   
   /*
    * Check whether the boundary element is an interior edge to two constraint
    * edges.  An interior edge is connected to at least one constraint edge at
    * each of its vertices.
    */
   bool is_interior() const;
   
   /*
    * Check whether the boundary element lies entirely along the image border.
    */
   bool is_border() const;

   /*
    * Get the number of pixels along the boundary element's smooth contour.
    */
   unsigned long size() const;

   /*
    * Get the length of the straight line segment approxmation to the boundary.
    */
   double length() const;

   /*
    * Get the orientation of the normal to the boundary.
    * The returned orientation is in the range [0, pi).
    */
   double orientation() const;

   /*
    * Get x-coordinates of points along the boundary element.
    */
   const array<unsigned long>& x_coords() const;

   /*
    * Get y-coordinates of points along the boundary element.
    */
   const array<unsigned long>& y_coords() const;
   
protected:
   /*
    * Protected constructor.
    * Create a boundary element corresponding to the given contour edge.
    * Initialize its status flags to false.
    */
   explicit boundary_element(const lib_image::contour_edge& /* contour edge */);
   
   /*
    * Boundary element data.
    */
   const lib_image::contour_edge& _edge;  /* corresponding contour edge */
   bool _is_intruder;                     /* is intruding edge? */
   bool _is_interior;                     /* is interior edge? */
   bool _is_border;                       /* is border edge? */
}; 
 
} /* namespace segmentation */
} /* namespace vision */

#endif
