/*
 * Line segment intersection.
 */
#ifndef MATH__GEOMETRY__SEGMENT_INTERSECTION_HH
#define MATH__GEOMETRY__SEGMENT_INTERSECTION_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/splay_set.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using collections::splay_set;
using lang::array;
using lang::pointers::auto_ptr;

/*
 * Line segment intersection.
 * Points are stored by reference in the segment intersection data structure.
 */
class seg_intersect {
public:
   /*
    * Constructor.
    *
    * Compute the intersection points of the given set of unique line segments.
    * For segments that intersect on an interval, only the endpoints of the
    * interval are reported as intersection points.
    *
    * Intersections are found using the O(n*log(n) + I*log(n)) sweep line 
    * algorithm, where n is the number of input segments and I is the number
    * of intersection points.
    *
    * Intersection points in the interior of segments cannot always be computed
    * exactly.  Intersections located within the optionally specified tolerance
    * distance from one another are considered to occur at the same point.
    
    * A tolerance of zero guarantees all endpoint and pairwise segment interior
    * intersections are detected.  However, correctly identifying a common
    * interior intersection of three or more segments may require a tolerance
    * greater than zero.
    */
   explicit seg_intersect(
      const collection<point_2D>&,  /* vertices of line segments */
      const array<unsigned long>&,  /* segment first endpoint ids */
      const array<unsigned long>&,  /* segment second endpoint ids */
      double = 0                    /* intersection distance tolerance */
   );

   /*
    * Copy constructor.
    */
   seg_intersect(const seg_intersect&);

   /*
    * Destructor.
    */
   virtual ~seg_intersect();
   
   /*
    * Vertex and segment lookup.
    *
    * Vertices and segments are numbered according to the order in which they
    * appear in the collection and arrays specified at construction time.
    *
    * Additional intersection vertices are numbered consecutively, starting 
    * with the next available id after the input vertices.
    */
   
   /*
    * Get number of vertices (including intersection points).
    */
   unsigned long vertices_size() const;

   /*
    * Get number of original segments.
    */
   unsigned long segments_size() const;
  
   /*
    * Get number of segments containing the specified vertex.
    */
   unsigned long intersection_size(unsigned long /* vertex id */) const;

   /*
    * Get the ids of segments that intersect at the specified vertex.
    */
   array<unsigned long> intersection(unsigned long /* vertex id */) const;

   /*
    * Check whether a vertex is in the interior of at least one segment.
    */
   bool is_interior_intersection(unsigned long /* vertex id */) const;
   
   /*
    * Get the ids of segments whos interior contains the specified vertex.
    */
   array<unsigned long> interior_intersection(
      unsigned long /* vertex id */
   ) const;
   
   /*
    * Get the ids of segments with an endpoint at the specified vertex.
    */
   array<unsigned long> endpoint_intersection(
      unsigned long /* vertex id */
   ) const;
   
   /*
    * Get the id of the specified endpoint of the segment.
    * The endpoints are numbered in lexicographic order.
    */
   unsigned long segment_vertex_id(
      unsigned long /* segment id */, unsigned long /* vertex # (0 or 1) */
   ) const;

   /*
    * Get the ids of both endpoints of the segment.
    * The endpoints appear in lexicographic order.
    */
   array<unsigned long> segment_vertex_ids(
      unsigned long /* segment id */
   ) const;

   /*
    * Return the specified vertex (by reference).
    */
   point_2D& vertex(unsigned long /* vertex id */) const;

protected:
   /************************************************************************
    * Line segment intersection data structure.
    ************************************************************************/
   
   /*
    * Declare vertex and segment types.
    */
   class vrtx;
   class seg;

   /*
    * Vertex.
    */
   class vrtx {
   public:
      /* constructor */
      explicit vrtx(unsigned long /* id */, point_2D& /* point */);

      /* data */
      unsigned long id;       /* vertex id */
      point_2D&     p;        /* point at vertex */
      list<seg>     lower;    /* segments whos lower endpoint is the vertex */
      list<seg>     upper;    /* segments whos upper endpoint is the vertex */
      list<seg>     inter;    /* segments whos interior intersects vertex */
   };
   
   /*
    * Segment.
    */
   class seg {
   public:
      /* constructor */
      explicit seg(
         unsigned long /* id */, vrtx& /* lower */, vrtx& /* upper */
      );

      /* check if segment is vertical */
      bool is_vertical() const;

      /* data */
      unsigned long id;       /* segment id */
      vrtx& v_lower;          /* lower endpoint */
      vrtx& v_upper;          /* upper endpoint */
   };
  
   /*
    * Intersection points that are not segment endpoints.
    */
   auto_collection< point_2D, array_list<point_2D> > _p_int; /* intersections */
   
   /*
    * Vertex and segment arrays.
    */
   auto_collection< vrtx, array_list<vrtx> > _vertices;  /* vertices */
   auto_collection< seg,  array_list<seg>  > _segments;  /* segments */

   /************************************************************************
    * Argument checking helper functions.
    ************************************************************************/

   /*
    * Check that the given vertex or segment id is valid.
    * Throw an exception (ex_index_out_of_bounds) if it is invalid.
    */
   void check_vertex_id(unsigned long /* vertex id */) const;
   void check_segment_id(unsigned long /* segment id */) const;

   /************************************************************************
    * Segment intersection helper functions.
    ************************************************************************/

   /*
    * Determine (exactly) whether the line segment contains the point.
    */
   static bool segment_contains_point(
      const seg& /* segment */, const point_2D& /* point */
   );

   /*
    * Determine (approximately) whether the line segment contains the point.
    * The segment contains the point iff the distance from the point to the 
    * segment is less than or equal to the specified tolerance.
    */
   static bool segment_contains_point(
      const seg& /* segment */, const point_2D& /* point */, double /* tol */
   );

   /*
    * Compute (approximately) the lexicographically smallest intersection point
    * of the two line segments if the line segments intersect.
    *
    * Otherwise, return a NULL pointer.
    *
    * In addition, if the segments intersect, update the distance tolerance 
    * to be at least the required distance tolerance in the neighborhood of
    * that intersection.
    *
    * The following numerical robustness guarantees are made:
    *
    * (1) Determination of whether or not the segments intersect is exact.
    *
    * (2) If the lexicographically smallest intersection point is an endpoint
    *     of one of the segments, that exact endpoint is returned.
    *
    * (3) If the lexicographically smallest intersection point lies in the
    *     interior of both segments, the approximate intersection is returned.
    */
   static auto_ptr<point_2D> compute_segment_intersection(
      const seg&, const seg&, double& /* global distance tolerance */
   );

   /*
    * Compute (approximately) the lexicographically smallest intersection point
    * of the given segment with the vertical sweep line through the specified
    * event point.
    *
    * If the segment does not intersect the sweep line, this functions throws
    * an exception (ex_invalid_argument).
    *
    * The following numerical robustness properties are guaranteed:
    *
    * (1) If the intersection point is a segment endpoint, that exact endpoint
    *     is returned.
    *
    * (2) In all other cases, an approximate intersection point is returned.
    *
    * (3) The intersection point has the sweep line's exact x-coordinate.
    */
   static point_2D compute_sweep_line_intersection(
      const seg& /* segment */, const point_2D& /* sweep line event point */
   );

   /*
    * Compute the angle of the given segment about its intersection with the
    * sweep line as viewed from the current sweep line event point.
    */
   static double compute_sweep_line_angle(
      const seg&,       /* segment */
      const point_2D&,  /* intersection point */
      const point_2D&   /* sweep line event point */
   );

   /*
    * Find the event (if it exists) corresponding to the intersection of the
    * two line segments.
    *
    * If this event vertex is both lexicographically greater than the current
    * event vertex and not yet contained in the event queue, then it is a new 
    * intersection and is enqueued into the event queue.
    *
    * In addition, upon discovering a new intersection point that is not a
    * segment endpoint, that point and its corresponding event vertex are added
    * to the segment intersection data structure.
    *
    * Also, update the global distance tolerance to reflect that required for 
    * any new intersection.
    */
   void find_event(
      seg&,             /* left segment */
      seg&,             /* right segment */
      const vrtx&,      /* current event vertex */
      splay_set<vrtx>&, /* event queue */
      double&           /* global distance tolerance */
   );
};

} /* namespace geometry */
} /* namespace math */

#endif
