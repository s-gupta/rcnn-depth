/*
 * Line segment intersection.
 */
#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/splay_set.hh"
#include "functors/comparable_functors.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/exact.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/seg_intersect.hh"
#include "math/math.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::splay_set;
using functors::comparable_functor;
using lang::array;
using lang::exceptions::ex_index_out_of_bounds;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::exact;

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

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
seg_intersect::seg_intersect(
   const collection<point_2D>& points,
   const array<unsigned long>& segment_start_ids,
   const array<unsigned long>& segment_end_ids,
   double                      tol)
 : _p_int(new array_list<point_2D>()),
   _vertices(new array_list<vrtx>()),
   _segments(new array_list<seg>())
{
   /* create vertices for points, compute sweep line starting point */
   unsigned long n_p = 0;
   point_2D p_min;
   for (auto_ptr< iterator<point_2D> > i = points.iter_create();
      i->has_next(); )
   {
      /* update sweep line starting point */
      point_2D& p = i->next();
      if ((p < p_min) || (n_p == 0))
         p_min = p;
      /* create vertex */
      auto_ptr<vrtx> v(new vrtx(n_p++, p));
      _vertices->add(*v);
      v.release();
   }
   /* create segments */
   unsigned long n_s = segment_start_ids.size();
   if (n_s != segment_end_ids.size())
      throw ex_invalid_argument("mismatch in segment endpoint array sizes");
   for (unsigned long n = 0; n < n_s; n++) {
      vrtx& v_start = (*_vertices)[segment_start_ids[n]];
      vrtx& v_end   = (*_vertices)[segment_end_ids[n]];
      bool v_cmp = (v_start.p < v_end.p);
      vrtx& v_lower = v_cmp ? v_start : v_end;
      vrtx& v_upper = v_cmp ? v_end : v_start;
      auto_ptr<seg> s(new seg(n, v_lower, v_upper));
      v_lower.lower.add(*s);
      v_upper.upper.add(*s);
      _segments->add(*s);
      s.release();
   }
   /* check that each vertex lies on at least one segment */
   for (unsigned long n = 0; n < n_p; n++) {
      vrtx& v = (*_vertices)[n];
      if (v.lower.is_empty() && v.upper.is_empty())
         throw ex_invalid_argument(
            "vertex specified which is not a segment endpoint"
         );
   }
   /* create comparison functor for event queue */
   class event_compare : public comparable_functor<vrtx> {
   public:
      event_compare(){} // EDIT: Jordi Pont-Tuset <jordi.pont@upc.edu> to compile on OSX Mavericks
      int operator()(const vrtx& v0, const vrtx& v1) const {
         return v0.p.compare_to(v1.p);
      }
   };
   static const event_compare event_cmp;
   /* create comparison functor for segment tree */
   class seg_compare : public comparable_functor<seg> {
   public:
      /* constructor - set sweep line event point */
      explicit seg_compare(point_2D p, double& tol) : _p(p), _tol(tol) { }
      
      /* access the sweep line event point */
      point_2D& p() { return _p; }
      const point_2D& p() const { return _p; }

      /* comparison function - order by intersection with sweep line */
      int operator()(const seg& s0, const seg& s1) const {
         /* compute intersection points with sweep line */
         point_2D p0 = seg_intersect::compute_sweep_line_intersection(s0, _p);
         point_2D p1 = seg_intersect::compute_sweep_line_intersection(s1, _p);
         int cmp = p0.compare_to(p1);
         /* if intersection points are identical, order by angle with line */
         if ((cmp == 0) || (math::abs(p0.y() - p1.y()) <= _tol)) {
            /* pick reference intersection point */
            point_2D& pr = (cmp < 0) ? p0 : p1;
            /* compute angle with sweep line */
            double ang0 = seg_intersect::compute_sweep_line_angle(s0, pr, _p);
            double ang1 = seg_intersect::compute_sweep_line_angle(s1, pr, _p);
            cmp = (ang0 < ang1) ? -1 : ((ang0 > ang1) ? 1 : 0);
            /* if angle with line is identical, order by segment id */
            static const double angle_tolerance = 16.0 * exact<>::epsilon();
            if ((cmp == 0) || (math::abs(ang0 - ang1) <= angle_tolerance))
               cmp = (s0.id < s1.id) ? -1 : ((s0.id > s1.id) ? 1 : 0);
         }
         return cmp;
      }
   protected:
      point_2D _p;      /* current sweep line event point */
      double&  _tol;    /* global distance tolerance */
   };
   seg_compare seg_cmp(p_min, tol);
   /* initialize event queue - enqueue all vertices */
   splay_set<vrtx> event_q(event_cmp);
   for (unsigned long n = 0; n < n_p; n++)
      event_q.add((*_vertices)[n]);
   /* initialize empty segment tree */
   splay_set<seg> seg_tree(seg_cmp);
   /* initialize vertical segment list */
   list<seg> vertical_segments;
   /* repeatedly dequeue and handle events */
   while (!(event_q.is_empty())) {
      /* dequeue event point */
      vrtx& v = event_q.find_min();
      event_q.remove(v);
      /* remove segments whos upper endpoint or interior intersects v */
      seg_tree.remove(v.upper);
      seg_tree.remove(v.inter);
      /* remove remaining vertical segments (their interior intersects v) */
      while (!(vertical_segments.is_empty())) {
         seg& s_vertical = vertical_segments.remove_head();
         if (seg_tree.contains(s_vertical)) {
            seg_tree.remove(s_vertical);
            v.inter.add(s_vertical);
         }
      }
      /* move the sweep line */
      seg_cmp.p() = v.p;
      /* grab a segment containing the event point */
      unsigned long n_lower = v.lower.size();
      unsigned long n_upper = v.upper.size();
      unsigned long n_inter = v.inter.size();
      seg& s =
         (n_lower == 0) ?
            ((n_upper == 0) ? v.inter.head() : v.upper.head())
          : (v.lower.head());
      /* find additional segments whos interior intersects v */
      while (seg_tree.contains_prev(s)) {
         seg& s_prev = seg_tree.find_prev(s);
         if (seg_intersect::segment_contains_point(s_prev, v.p, tol)) {
            seg_tree.remove(s_prev);
            v.inter.add(s_prev);
         } else {
            break;
         }
      }
      while (seg_tree.contains_next(s)) {
         seg& s_next = seg_tree.find_next(s);
         if (seg_intersect::segment_contains_point(s_next, v.p, tol)) {
            seg_tree.remove(s_next);
            v.inter.add(s_next);
         } else {
            break;
         }
      }
      /* update number of interior intersections */
      n_inter = v.inter.size();
      /* add segments whos lower endpoint or interior intersects v */
      seg_tree.add(v.lower);
      seg_tree.add(v.inter);
      /* sort lists of segments with lower endpoint or interior intersections */
      v.lower.sort(seg_cmp);
      v.inter.sort(seg_cmp);
      /* update list of vertical segments in tree */
      for (auto_ptr< iterator<seg> > i = v.lower.iter_create(); i->has_next(); )
      {
         seg& s_lower = i->next();
         if (s_lower.is_vertical())
            vertical_segments.add(s_lower);
      }
      for (auto_ptr< iterator<seg> > i = v.inter.iter_create(); i->has_next(); )
      {
         seg& s_inter = i->next();
         if (s_inter.is_vertical())
            vertical_segments.add(s_inter);
      }
      /* find new event points */
      if ((n_lower + n_inter) == 0) {
         /* check for left/right neighbors of segment that ended at v */
         if (seg_tree.contains_prev(s) && seg_tree.contains_next(s)) {
            seg& s_left  = seg_tree.find_prev(s);
            seg& s_right = seg_tree.find_next(s);
            this->find_event(s_left, s_right, v, event_q, tol);
         }
      } else {
         /* extract leftmost segment through event point */
         seg& s_lower_left = (n_lower == 0) ? v.inter.head() : v.lower.head();
         seg& s_inter_left = (n_inter == 0) ? v.lower.head() : v.inter.head();
         seg& s_left =
            (seg_cmp(s_lower_left, s_inter_left) < 0) ?
               s_lower_left : s_inter_left;
         /* extract rightmost segment through event point */
         seg& s_lower_right = (n_lower == 0) ? v.inter.tail() : v.lower.tail();
         seg& s_inter_right = (n_inter == 0) ? v.lower.tail() : v.inter.tail();
         seg& s_right =
            (seg_cmp(s_lower_right, s_inter_right) > 0) ?
               s_lower_right : s_inter_right;
         /* find new event on left side */
         if (seg_tree.contains_prev(s_left)) {
            seg& s_left_left = seg_tree.find_prev(s_left);
            this->find_event(s_left_left, s_left, v, event_q, tol);
         }
         /* find new event on right side */
         if (seg_tree.contains_next(s_right)) {
            seg& s_right_right = seg_tree.find_next(s_right);
            this->find_event(s_right, s_right_right, v, event_q, tol);
         }
      }
   }
}

/*
 * Copy constructor.
 */
seg_intersect::seg_intersect(const seg_intersect& seg_int)
 : _p_int(new array_list<point_2D>()),
   _vertices(new array_list<vrtx>()),
   _segments(new array_list<seg>())
{
   /* copy non-endpoint intersection points */
   unsigned long n_p = seg_int._p_int->size();
   for (unsigned long n = 0; n < n_p; n++) {
      const point_2D& p = (*(seg_int._p_int))[n];
      auto_ptr<point_2D> p_copy(new point_2D(p));
      _p_int->add(*p_copy);
      p_copy.release();
   }
   /* copy vertices */
   unsigned long n_v = seg_int._vertices->size();
   unsigned long n_v_endpoint = n_v - n_p;
   for (unsigned long n = 0; n < n_v_endpoint; n++) {
      const vrtx& v = (*(seg_int._vertices))[n];
      auto_ptr<vrtx> v_copy(new vrtx(v.id, v.p));
      _vertices->add(*v_copy);
      v_copy.release();
   }
   for (unsigned long n = 0; n < n_p; n++) {
      const vrtx& v = (*(seg_int._vertices))[n + n_v_endpoint];
      auto_ptr<vrtx> v_copy(new vrtx(v.id, (*_p_int)[n]));
      _vertices->add(*v_copy);
      v_copy.release();
   }
   /* copy segments */
   unsigned long n_s = seg_int._segments->size();
   for (unsigned long n = 0; n < n_s; n++) {
      const seg& s = (*(seg_int._segments))[n];
      vrtx& v_lower = (*_vertices)[s.v_lower.id];
      vrtx& v_upper = (*_vertices)[s.v_upper.id];
      auto_ptr<seg> s_copy(new seg(s.id, v_lower, v_upper));
      _segments->add(*s_copy);
      s_copy.release();
   }
   /* build vertex segment lists */
   for (unsigned long n = 0; n < n_v; n++) {
      const vrtx& v = (*(seg_int._vertices))[n];
      vrtx& v_copy = (*_vertices)[n];
      auto_ptr< iterator<seg> > i = v.lower.iter_create();
      while (i->has_next())
         v_copy.lower.add((*_segments)[i->next().id]);
      i = v.upper.iter_create();
      while (i->has_next())
         v_copy.upper.add((*_segments)[i->next().id]);
      i = v.inter.iter_create();
      while (i->has_next())
         v_copy.inter.add((*_segments)[i->next().id]);
   }
}

/*
 * Destructor.
 */
seg_intersect::~seg_intersect() {
   /* do nothing */
}

/***************************************************************************
 * Size.
 ***************************************************************************/

/*
 * Get number of vertices (including intersection points).
 */
unsigned long seg_intersect::vertices_size() const {
   return _vertices->size();
}

/*
 * Get number of original segments.
 */
unsigned long seg_intersect::segments_size() const {
   return _segments->size();
}

/*
 * Get number of segments containing the specified vertex.
 */
unsigned long seg_intersect::intersection_size(unsigned long v_id) const {
   this->check_vertex_id(v_id);
   vrtx& v = (*_vertices)[v_id];
   return (v.lower.size() + v.upper.size() + v.inter.size());
}

/***************************************************************************
 * Vertex and segment lookup.
 ***************************************************************************/

/*
 * Get the ids of segments that intersect at the specified vertex.
 */
array<unsigned long> seg_intersect::intersection(unsigned long v_id) const {
   this->check_vertex_id(v_id);
   vrtx& v = (*_vertices)[v_id];
   unsigned long n_lower = v.lower.size();
   unsigned long n_upper = v.upper.size();
   unsigned long n_inter = v.inter.size();
   unsigned long n_seg = n_lower + n_upper + n_inter;
   array<unsigned long> seg_ids(n_seg);
   unsigned long n = 0;
   for (auto_ptr< iterator<seg> > i = v.lower.iter_create(); i->has_next(); ) {
      seg& s = i->next();
      seg_ids[n++] = s.id;
   }
   for (auto_ptr< iterator<seg> > i = v.upper.iter_create(); i->has_next(); ) {
      seg& s = i->next();
      seg_ids[n++] = s.id;
   }
   for (auto_ptr< iterator<seg> > i = v.inter.iter_create(); i->has_next(); ) {
      seg& s = i->next();
      seg_ids[n++] = s.id;
   }
   return seg_ids;
}

/*
 * Check whether a vertex is in the interior of at least one segment.
 */
bool seg_intersect::is_interior_intersection(unsigned long v_id) const {
   this->check_vertex_id(v_id);
   vrtx& v = (*_vertices)[v_id];
   return (v.inter.size() > 0);
}

/*
 * Get the ids of segments whos interior contains the specified vertex.
 */
array<unsigned long> seg_intersect::interior_intersection(
   unsigned long v_id) const
{
   this->check_vertex_id(v_id);
   vrtx& v = (*_vertices)[v_id];
   unsigned long n_inter = v.inter.size();
   array<unsigned long> seg_ids(n_inter);
   unsigned long n = 0;
   for (auto_ptr< iterator<seg> > i = v.inter.iter_create(); i->has_next(); ) {
      seg& s = i->next();
      seg_ids[n++] = s.id;
   }
   return seg_ids;
}

/*
 * Get the ids of segments with an endpoint at the specified vertex.
 */
array<unsigned long> seg_intersect::endpoint_intersection(
   unsigned long v_id) const
{
   this->check_vertex_id(v_id);
   vrtx& v = (*_vertices)[v_id];
   unsigned long n_lower = v.lower.size();
   unsigned long n_upper = v.upper.size();
   unsigned long n_seg = n_lower + n_upper;
   array<unsigned long> seg_ids(n_seg);
   unsigned long n = 0;
   for (auto_ptr< iterator<seg> > i = v.lower.iter_create(); i->has_next(); ) {
      seg& s = i->next();
      seg_ids[n++] = s.id;
   }
   for (auto_ptr< iterator<seg> > i = v.upper.iter_create(); i->has_next(); ) {
      seg& s = i->next();
      seg_ids[n++] = s.id;
   }
   return seg_ids;
}
   
/*
 * Get the id of the specified endpoint of the segment.
 * The endpoints are numbered in lexicographic order.
 */
unsigned long seg_intersect::segment_vertex_id(
   unsigned long s_id, unsigned long v_num) const
{
   this->check_segment_id(s_id);
   if (v_num == 0)
      return (*_segments)[s_id].v_lower.id;
   else if (v_num == 1)
      return (*_segments)[s_id].v_upper.id;
   else
      throw ex_invalid_argument("vertex # must be 0 or 1");
}
   
/*
 * Get the ids of both endpoints of the segment.
 * The endpoints appear in lexicographic order.
 */
array<unsigned long> seg_intersect::segment_vertex_ids(
   unsigned long s_id) const
{
   this->check_segment_id(s_id);
   array<unsigned long> v_ids(2);
   v_ids[0] = (*_segments)[s_id].v_lower.id;
   v_ids[1] = (*_segments)[s_id].v_upper.id;
   return v_ids;
}

/*
 * Return the specified vertex (by reference).
 */
point_2D& seg_intersect::vertex(unsigned long v_id) const {
   this->check_vertex_id(v_id);
   vrtx& v = (*_vertices)[v_id];
   return v.p;
}

/***************************************************************************
 * Argument checking helper functions.
 ***************************************************************************/

/* 
 * Check that the given vertex id is valid.
 * Throw an exception (ex_index_out_of_bounds) if it is invalid.
 */
void seg_intersect::check_vertex_id(unsigned long v_id) const {
   if (v_id >= _vertices->size())
      throw ex_index_out_of_bounds("invalid vertex id", v_id);
}

/* 
 * Check that the given segment id is valid.
 * Throw an exception (ex_index_out_of_bounds) if it is invalid.
 */
void seg_intersect::check_segment_id(unsigned long s_id) const {
   if (s_id >= _segments->size())
      throw ex_index_out_of_bounds("invalid segment id", s_id);
}

/***************************************************************************
 * Vertex and segment classes.
 ***************************************************************************/

/*
 * Vertex constructor.
 */
seg_intersect::vrtx::vrtx(unsigned long id, point_2D& p)
 : id(id), p(p), lower(), upper(), inter()
{ }

/*
 * Segment constructor.
 */
seg_intersect::seg::seg(unsigned long id, vrtx& v_lower, vrtx& v_upper)
 : id(id), v_lower(v_lower), v_upper(v_upper)
{ }

/*
 * Check if segment is vertical.
 */
bool seg_intersect::seg::is_vertical() const {
   return (v_lower.p.x() == v_upper.p.x());
}

/***************************************************************************
 * Segment intersection helper functions.
 ***************************************************************************/

/*
 * Determine (exactly) whether the line segment contains the point.
 */
bool seg_intersect::segment_contains_point(
   const seg& s, const point_2D& p)
{
   /* extract segment endpoints */
   point_2D& p_lower = s.v_lower.p;
   point_2D& p_upper = s.v_upper.p;
   /* determine if p is between lower and upper endpoints and on the line */
   return (
      (p_lower <= p) && (p <= p_upper) && (p.orientation(p_lower, p_upper) == 0)
   );
}

/*
 * Determine (approximately) whether the line segment contains the point.
 * The segment contains the point iff the distance from the point to the 
 * segment is less than or equal to the specified tolerance.
 */
bool seg_intersect::segment_contains_point(
   const seg& s, const point_2D& p, double tol)
{
   /* check arguments */
   if (tol < 0)
      throw ex_invalid_argument("tolerance must be >= 0");
   /* check both approximate and exact conditions for containment   */
   /* (note that the approximation may report a nonzero distance    */
   /*  even if the exact distance is zero, so both must be checked) */
   return (
      (p.distance_to_segment(s.v_lower.p, s.v_upper.p) <= tol) ||
      (seg_intersect::segment_contains_point(s, p))
   );
}

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
auto_ptr<point_2D> seg_intersect::compute_segment_intersection(
   const seg& s0, const seg& s1, double& tol)
{
   /* get segment endpoints */
   point_2D& pa = s0.v_lower.p;
   point_2D& pb = s0.v_upper.p;
   point_2D& pc = s1.v_lower.p;
   point_2D& pd = s1.v_upper.p;
   /* perform orientation tests relative to segment ab */
   int ab_c = pc.orientation(pa, pb);
   int ab_d = pd.orientation(pa, pb);
   /* check if segments lie on the same line */
   if ((ab_c == 0) && (ab_d == 0)) {
      if ((pd < pa) || (pc > pb)) {
         /* segments do not interset */
         return auto_ptr<point_2D>(NULL);
      } else if (pa > pc) {
         /* vertex a is the lexicographically smallest intersection */
         return auto_ptr<point_2D>(new point_2D(pa));
      } else {
         /* vertex c is the lexicographically smallest intersection */
         return auto_ptr<point_2D>(new point_2D(pc));
      }
   }
   /* check if segment ab contains c */
   if ((ab_c == 0) && (pa <= pc) && (pc <= pb))
      return auto_ptr<point_2D>(new point_2D(pc));
   /* check if segment ab contains d */
   if ((ab_d == 0) && (pa <= pd) && (pd <= pb))
      return auto_ptr<point_2D>(new point_2D(pd));
   /* perform orientation tests relative to segment cd */
   int cd_a = pa.orientation(pc, pd);
   int cd_b = pb.orientation(pc, pd);
   /* check if segment cd contains a */
   if ((cd_a == 0) && (pc <= pa) && (pa <= pd)) 
      return auto_ptr<point_2D>(new point_2D(pa));
   /* check if segment cd contains b */
   if ((cd_b == 0) && (pc <= pb) && (pb <= pd)) 
      return auto_ptr<point_2D>(new point_2D(pb));
   /* check if segment interiors intersect */
   if (((ab_c * ab_d) < 0) && ((cd_a * cd_b) < 0)) {
      /* get endpoint coordinates */
      double xa = pa.x();  double ya = pa.y();
      double xb = pb.x();  double yb = pb.y();
      double xc = pc.x();  double yc = pc.y();
      double xd = pd.x();  double yd = pd.y();
      /* check for vertical segments */
      if (xa == xb) {
         /* first segment is vertical - compute sweep line intersection */
         return auto_ptr<point_2D>(
            new point_2D(
               seg_intersect::compute_sweep_line_intersection(s1, pa)
            )
         );
      } else if (xc == xd) {
         /* second segment is vertical - compute sweep line intersection */
         return auto_ptr<point_2D>(
            new point_2D(
               seg_intersect::compute_sweep_line_intersection(s0, pc)
            )
         );
      } else {
         /* compute fractional distance along segment to intersection */
         double alpha =
            ((xd - xc)*(ya - yc) - (yd - yc)*(xa - xc)) /
            ((yd - yc)*(xb - xa) - (xd - xc)*(yb - ya));
         /* compute approximate x-coordinate of intersection */
         double x = (xa + alpha*(xb - xa));
         /* compute approximate y-coordinates of intersection */
         double y_s0_approx = ya + (yb - ya)/(xb - xa)*(x - xa);
         double y_s1_approx = yc + (yd - yc)/(xd - xc)*(x - xc);
         /* choose largest y-coordinate to be consistent with sweep logic */
         double y = (y_s0_approx > y_s1_approx) ? y_s0_approx : y_s1_approx;
         /* update global distance tolerance */
         double min_tol = 2.0 * math::abs(y_s0_approx - y_s1_approx);
         if (tol < min_tol)
            tol = min_tol;
         /* return intersection point */
         return auto_ptr<point_2D>(new point_2D(x, y));
      }
   } else {
      /* segments do not intersect */
      return auto_ptr<point_2D>(NULL);
   }
}

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
point_2D seg_intersect::compute_sweep_line_intersection(
   const seg& s, const point_2D& p)
{
   /* get segment endpoints */
   point_2D& pa = s.v_lower.p;
   point_2D& pb = s.v_upper.p;
   /* get endpoint coordinates */
   double xa = pa.x();  double ya = pa.y();
   double xb = pb.x();  double yb = pb.y();
   /* get sweep line event point coordinates */
   double x = p.x();  double y = p.y();
   /* determine intersection point */
   if ((xa > x) || (xb < x)) {
      /* segment does not intersect sweep line */
      throw ex_invalid_argument(
         "segment does not intersect sweep line"
      );
   } else if (xa == xb) {
      /* vertical segment - intersection must be on or above event point */
      if (ya >= y)
         return pa;
      else if (y <= yb)
         return p;
      else
         throw ex_invalid_argument(
            "vertical segment does not intersect sweep line"
         );
   } else if (xa == x) {
      /* lexicographically smallest endpoint intersects sweep line */
      return pa;
   } else if (xb == x) {
      /* lexicographically largest endpoint intersects sweep line */
      return pb;
   } else {
      /* compute approximate intersection point with segment interior */
      return point_2D(x, ya + (yb - ya)/(xb - xa)*(x - xa));
   }
}

/*
 * Compute the angle of the given segment about its intersection with the
 * sweep line as viewed from the current sweep line event point.
 */
double seg_intersect::compute_sweep_line_angle(
   const seg&      s,
   const point_2D& p_intersect,
   const point_2D& p_event)
{
   if (p_event < p_intersect) {
      /* compute angle using the part of the segment left of the sweep line */
      const point_2D& p = s.v_lower.p;
      return math::atan2(p_intersect.x() - p.x(), p_intersect.y() - p.y());
   } else {
      /* compute angle using the part of the segment right of the sweep line */
      const point_2D& p = s.v_upper.p;
      return math::atan2(p.x() - p_intersect.x(), p_intersect.y() - p.y());
   }
}

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
void seg_intersect::find_event(
   seg&             s_left,
   seg&             s_right,
   const vrtx&      v,
   splay_set<vrtx>& event_q,
   double&          tol)
{
   /* compute segment intersection (if it exists) */
   auto_ptr<point_2D> auto_p_event =
      seg_intersect::compute_segment_intersection(s_left, s_right, tol);
   point_2D* p_event = auto_p_event.get();
   /* check if intersection exists and appears after current event vertex */
   if ((p_event != NULL) && (*p_event > v.p)) {
      /* initialize vertex for event point */
      auto_ptr<vrtx> auto_v_event(NULL);
      vrtx* v_event = NULL;
      /* check if vertex is endpoint of left segment */
      bool is_left_endpoint = true;
      if (*p_event == s_left.v_lower.p) {
         v_event = &(s_left.v_lower);
      } else if (*p_event == s_left.v_upper.p) {
         v_event = &(s_left.v_upper);
      } else {
         is_left_endpoint = false;
      }
      /* check if vertex is endpoint of right segment */
      bool is_right_endpoint = true;
      if (*p_event == s_right.v_lower.p) {
         v_event = &(s_right.v_lower);
      } else if (*p_event == s_right.v_upper.p) {
         v_event = &(s_right.v_upper);
      } else {
         is_right_endpoint = false;
      }
      /* create vertex if it is not an endpoint */
      bool is_endpoint = is_left_endpoint || is_right_endpoint;
      if (!is_endpoint) {
         unsigned long n_v = _vertices->size();
         auto_v_event.reset(new vrtx(n_v, *p_event));
         v_event = auto_v_event.get();
      }
      /* if vertex is in interior of either segment, record that fact */
      if (!is_left_endpoint)
         v_event->inter.add(s_left);
      if (!is_right_endpoint)
         v_event->inter.add(s_right);
      /* check if vertex is contained in event queue */
      if (!event_q.contains(*v_event)) {
         /* not in event queue - add vertex to it */
         event_q.add(*v_event);
         /* if new vertex (not an endpoint), update point and vertex arrays */
         if (!is_endpoint) {
            _p_int->add(*p_event);     auto_p_event.release();
            _vertices->add(*v_event);  auto_v_event.release();
         }
      }
   }
}

} /* namespace geometry */
} /* namespace math */
