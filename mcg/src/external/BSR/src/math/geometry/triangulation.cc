/*
 * Triangulation.
 */
#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/child_thread.hh"
#include "functors/comparable_functors.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/sym_edge.hh"
#include "math/geometry/triangle_2D.hh"
#include "math/geometry/triangulation.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using concurrent::threads::child_thread;
using functors::comparable_functor;
using lang::array;
using lang::exceptions::ex_index_out_of_bounds;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

/*
 * Constructor.
 * Create an empty triangulation.
 */
triangulation::triangulation()
 : _vertices(new array_list<vrtx>()),
   _edges(new array_list< sym_edge<vrtx,edge> >()),
   _edge_data(new list<edge>()),
   _triangles(new array_list<tri>())
{ }

/*
 * Copy constructor.
 */
triangulation::triangulation(const triangulation& t_map)
 : _vertices(new array_list<vrtx>()),
   _edges(new array_list< sym_edge<vrtx,edge> >()),
   _edge_data(new list<edge>()),
   _triangles(new array_list<tri>())
{
   /* copy vertices */
   unsigned long n_v = t_map._vertices->size();
   for (unsigned long n = 0; n < n_v; n++) {
      vrtx& v = (*(t_map._vertices))[n];
      auto_ptr<vrtx> v_copy(new vrtx(v.id, v.p));
      _vertices->add(*v_copy);
      v_copy.release();
   }
   /* copy edges and edge data */
   unsigned long n_e = t_map._edges->size();
   for (unsigned long n = 0; n < n_e; n++) {
      sym_edge<vrtx,edge>& e = (*(t_map._edges))[n];
      unsigned long e_id = e.data().id;
      auto_ptr<edge> e_data(new edge(e_id));
      auto_ptr<edge> e_sym_data(new edge(e_id));
      unsigned long v_org_id  = e.origin().id;
      unsigned long v_dest_id = e.destination().id;
      auto_ptr< sym_edge<vrtx,edge> > e_copy(new sym_edge<vrtx,edge>(
         (*_vertices)[v_org_id], (*_vertices)[v_dest_id], *e_data, *e_sym_data
      ));
      _edge_data->add(*e_data);      e_data.release();
      _edge_data->add(*e_sym_data);  e_sym_data.release();
      _edges->add(*e_copy);          e_copy.release();
   }
   /* copy triangles */
   unsigned long n_t = t_map._triangles->size();
   for (unsigned long n = 0; n < n_t; n++) {
      tri& t = (*(t_map._triangles))[n];
      auto_ptr<tri> t_copy(new tri(t.id));
      _triangles->add(*t_copy);
      t_copy.release();
   }
   /* assemble edge rings */
   for (unsigned long n = 0; n < n_v; n++) {
      vrtx& v_src  = (*(t_map._vertices))[n];
      vrtx& v_copy = (*_vertices)[n];
      if (v_src.e != NULL) {
         /* link vertex -> edge out of vertex */
         unsigned long e_id = v_src.e->data().id;
         sym_edge<vrtx,edge>* e = &((*_edges)[e_id]);
         if (&(e->origin()) != &v_copy)
            e = &(e->sym());
         v_copy.e = e;
         /* splice other edges into ring */
         sym_edge<vrtx,edge>* e_src_stop = v_src.e;
         sym_edge<vrtx,edge>* e_src  = &(e_src_stop->origin_next());
         while (e_src != e_src_stop) {
            sym_edge<vrtx,edge>* e_copy = e;
            e_id = e_src->data().id;
            e = &((*_edges)[e_id]);
            if (&(e->origin()) != &v_copy)
               e = &(e->sym());
            sym_edge<vrtx,edge>::splice(*e_copy, *e);
            e_src = &(e_src->origin_next());
         }
      }
   }
   /* link edges <-> triangles */
   for (unsigned long n = 0; n < n_e; n++) {
      /* get edge */
      sym_edge<vrtx,edge>& e      = (*(t_map._edges))[n];
      sym_edge<vrtx,edge>& e_copy = (*_edges)[n];
      /* link left triangle */
      tri* e_t_left = e.data().t_left;
      if (e_t_left != NULL) {
         unsigned long t_id_left = e_t_left->id;
         tri& t = (*_triangles)[t_id_left];
         e_copy.data().t_left = &t;
         if ((*(t_map._triangles))[t_id_left].e == &e)
            t.e = &e_copy;
      }
      /* link right triangle */
      tri* e_t_right = e.sym().data().t_left;
      if (e_t_right != NULL) {
         unsigned long t_id_right = e_t_right->id;
         tri& t = (*_triangles)[t_id_right];
         e_copy.sym().data().t_left = &t;
         if ((*(t_map._triangles))[t_id_right].e == &(e.sym()))
            t.e = &(e_copy.sym());
      }
   }
}

/*
 * Destructor.
 */
triangulation::~triangulation() {
   /* do nothing */
}

/***************************************************************************
 * Named constructors.
 ***************************************************************************/

/*
 * Delaunay triangulation.
 * Create a Delaunay triangulation of the given set of unique points.
 *
 * The O(n + (n/p)*log(n/p)) divide and conquer algorithm is used to
 * construct the Delaunay triangulation, where n is the number of points
 * and p is the number of available processors.
 */
auto_ptr<triangulation> triangulation::delaunay(
   const collection<point_2D>& points)
{
   array<unsigned long> constraint_start_ids;
   array<unsigned long> constraint_end_ids;
   return triangulation::delaunay(
      points, constraint_start_ids, constraint_end_ids
   );
}

/*
 * Constrained Delaunay triangulation (CDT).
 *
 * Points are numbered according to the order in which they appear in the
 * specified collection.  The constraint arrays specify, for each constraint
 * edge, the ids of the two endpoints of the edge.
 *
 * The constraints must be consistent, in that no two constraint edges 
 * cross.  Otherwise, an exception (ex_invalid_argument) is throw.
 */
auto_ptr<triangulation> triangulation::delaunay(
   const collection<point_2D>& points,
   const array<unsigned long>& constraint_start_ids,
   const array<unsigned long>& constraint_end_ids)
{
   /* build vertices */
   auto_ptr<triangulation> t_map(new triangulation());
   t_map->build_vertices(points);
   /* sort vertices in lexicographic order */
   class vrtx_compare : public comparable_functor<vrtx> {
   public:
       vrtx_compare(){} // EDIT: Jordi Pont-Tuset <jordi.pont@upc.edu> to compile on OSX Mavericks
      int operator()(const vrtx& v0, const vrtx& v1) const {
         return v0.p.compare_to(v1.p);
      }
   };
   static const vrtx_compare f_compare;
   t_map->_vertices->sort(f_compare);
   /* check sorted vertices for duplicates */
   unsigned long n_v = t_map->_vertices->size();
   for (unsigned long n = 1; n < n_v; n++) {
      if ((*(t_map->_vertices))[n-1].p == (*(t_map->_vertices))[n].p)
         throw ex_invalid_argument("points in set must be unique");
   }
   /* check constraints arrays */
   unsigned long n_constraints = constraint_start_ids.size();
   if (n_constraints != constraint_end_ids.size())
      throw ex_invalid_argument("mismatch in constraint array sizes");
   for (unsigned long n = 0; n < n_constraints; n++) {
      if ((constraint_start_ids[n] > n_v) || (constraint_end_ids[n] > n_v))
         throw ex_invalid_argument(
            "constraints specified for vertices that don't exist"
         );
   }
   /* build delaunay triangulation using divide and conquer algorithm */
   if (n_v > 0) {
      /* compute delaunay triangulation */
      auto_collection< sym_edge<vrtx,edge>, list< sym_edge<vrtx,edge> > > edges;
      delaunay_triangulator dt(*(t_map->_vertices), 0, n_v - 1, edges);
      dt.run();
      /* place edges into edge array */
      unsigned long n_e = edges->size();
      for (unsigned long n = 0; n < n_e; n++) {
         auto_ptr< sym_edge<vrtx,edge> > e(&(edges->remove_head()));
         t_map->_edges->add(*e);
         e.release();
      }
   }
   /* put vertices back in original order */
   for (unsigned long n = 0; n < n_v; ) {
      vrtx& v = (*(t_map->_vertices))[n];
      if (v.id != n) {
         /* swap vertex into correct position */
         vrtx& v_temp = t_map->_vertices->replace(v.id, v);
         t_map->_vertices->replace(n, v_temp);
      } else {
         n++;
      }
   }
   /* link vertices to edge rings */
   t_map->build_vertex_edge_links();
   /* label all edges with the maximum id */
   unsigned long n_e = t_map->_edges->size();
   for (unsigned long n = 0; n < n_e; n++) {
      sym_edge<vrtx,edge>& e = (*(t_map->_edges))[n];
      auto_ptr<edge> e_data(new edge(n_e));
      auto_ptr<edge> e_sym_data(new edge(n_e));
      t_map->_edge_data->add(*e_data);
      t_map->_edge_data->add(*e_sym_data);
      e.data(*(e_data.release()));
      e.sym().data(*(e_sym_data.release()));
   }
   /* insert constraint edges */
   for (unsigned long n = 0; n < n_constraints; n++)
      t_map->insert_edge(constraint_start_ids[n], constraint_end_ids[n], n);
   /* label non-constraint edges */
   for (unsigned long n = 0, id = n_constraints; n < n_e; n++) {
      sym_edge<vrtx,edge>& e = (*(t_map->_edges))[n];
      unsigned long e_id = e.data().id;
      if (e_id == n_e) {
         e.data().id       = id;
         e.sym().data().id = id;
         id++;
      }
   }
   /* move edges into position by id (constraint edges to front of array) */
   for (unsigned long n = 0; n < n_e; ) {
      sym_edge<vrtx,edge>& e = (*(t_map->_edges))[n];
      unsigned long e_id = e.data().id;
      if (e_id != n) {
         /* swap edge into correct position */
         sym_edge<vrtx,edge>& e_temp = t_map->_edges->replace(e_id, e);
         t_map->_edges->replace(n, e_temp);
      } else {
         n++;
      }
   }
   /* create triangle nodes */
   t_map->build_triangles();
   return t_map;
}

/***************************************************************************
 * Size.
 ***************************************************************************/

/*
 * Get number of vertices in triangulation.
 */
unsigned long triangulation::vertices_size() const {
   return _vertices->size();
}

/*
 * Get number of edges in triangulation.
 */
unsigned long triangulation::edges_size() const {
   return _edges->size();
}

/*
 * Get number of triangles in triangulation.
 */
unsigned long triangulation::triangles_size() const {
   return _triangles->size();
}

/*
 * Get the number of edges containing the specified vertex.
 */
unsigned long triangulation::vertex_edges_size(unsigned long v_id) const {
   /* check argument */
   this->check_vertex_id(v_id);
   /* get edge out of vertex */
   sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
   /* count edges */
   unsigned long n_e = 0;
   if (e != NULL) {
      sym_edge<vrtx,edge>* e_curr = e;
      do {
         n_e++;
         e_curr = &(e_curr->origin_next());
      } while (e_curr != e);
   }
   return n_e;
}

/*
 * Get the number of triangles containing the specified vertex.
 */
unsigned long triangulation::vertex_triangles_size(unsigned long v_id) const {
   /* check argument */
   this->check_vertex_id(v_id);
   /* get edge out of vertex */
   sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
   /* count triangles */
   unsigned long n_t = 0;
   if (e != NULL) {
      sym_edge<vrtx,edge>* e_curr = e;
      do {
         tri* t = e_curr->data().t_left;
         if (t != NULL)
            n_t++;
         e_curr = &(e_curr->origin_next());
      } while (e_curr != e);
   }
   return n_t;
}

/***************************************************************************
 * Vertex -> edge, triangle lookup.
 ***************************************************************************/

/*
 * Get the id of the specified edge containing the vertex.
 */
unsigned long triangulation::vertex_edge_id(
   unsigned long v_id, unsigned long e_num) const
{
   /* check argument */
   this->check_vertex_id(v_id);
   /* get edge out of vertex */
   sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
   if (e == NULL)
      throw ex_index_out_of_bounds("specified edge does not exist", e_num);   
   /* move to requested edge */
   sym_edge<vrtx,edge>* e_curr = e;
   unsigned long n = 0;
   while (n < e_num) {
      e_curr = &(e_curr->origin_next());
      if (e_curr == e)
         throw ex_index_out_of_bounds("specified edge does not exist", e_num);
      n++;
   }
   return e_curr->data().id;
}

/*
 * Get the id of the specified triangle containing the vertex.
 */
unsigned long triangulation::vertex_triangle_id(
   unsigned long v_id, unsigned long t_num) const
{
   /* check argument */
   this->check_vertex_id(v_id);
   /* get edge out of vertex */
   sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
   if (e == NULL)
      throw ex_index_out_of_bounds("specified triangle does not exist", t_num);
   /* move to requested triangle */
   sym_edge<vrtx,edge>* e_curr = e;
   tri* t = e->data().t_left;
   unsigned long n = (t != NULL) ? 1 : 0;
   while (n <= t_num) {
      e_curr = &(e_curr->origin_next());
      if (e_curr == e)
         throw ex_index_out_of_bounds(
            "specified triangle does not exist", t_num
         );
      t = e_curr->data().t_left;
      if (t != NULL)
         n++;
   }
   return t->id;
}

/*
 * Get the ids of the edges containing the vertex.
 */
array<unsigned long> triangulation::vertex_edge_ids(
   unsigned long v_id) const
{
   /* get number of edges containing vertex */
   unsigned long n_e = this->vertex_edges_size(v_id);
   /* build array of edge ids */
   array<unsigned long> e_ids(n_e);
   sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
   unsigned long n = 0;
   while (n < n_e) {
      e_ids[n] = e->data().id;
      n++;
      e = &(e->origin_next());
   }
   return e_ids;
}

/*
 * Get the ids of the triangles containing the vertex.
 */
array<unsigned long> triangulation::vertex_triangle_ids(
   unsigned long v_id) const
{
   /* get number of triangles containing vertex */
   unsigned long n_t = this->vertex_triangles_size(v_id);
   /* build array of triangle ids */
   array<unsigned long> t_ids(n_t);
   sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
   unsigned long n = 0;
   while (n < n_t) {
      tri* t = e->data().t_left;
      if (t != NULL) { 
         t_ids[n] = t->id;
         n++;
      }
      e = &(e->origin_next());
   }
   return t_ids;
}

/***************************************************************************
 * Edge -> vertex, triangle lookup.
 ***************************************************************************/

/*
 * Get the id of the specified vertex of the edge.
 */
unsigned long triangulation::edge_vertex_id(
   unsigned long e_id, unsigned long v_num) const
{
   this->check_edge_id(e_id);
   sym_edge<vrtx,edge>& e = (*_edges)[e_id];
   if (v_num == 0)
      return e.origin().id;
   else if (v_num == 1)
      return e.destination().id;
   else
      throw ex_index_out_of_bounds(
         "edge vertex number must be 0 or 1", v_num
      );
}

/*
 * Get the id of the specified triangle adjacent to the edge.
 */
unsigned long triangulation::edge_triangle_id(
   unsigned long e_id, unsigned long t_num) const
{
   this->check_edge_id(e_id);
   sym_edge<vrtx,edge>& e = (*_edges)[e_id];
   tri* t = NULL;
   if (t_num == 0)
      t = e.data().t_left;
   else if (t_num == 1)
      t = e.sym().data().t_left;
   else
      throw ex_index_out_of_bounds(
         "edge triangle number must be 0 or 1", t_num
      );
   return (t != NULL) ? (t->id) : (_triangles->size());
}

/*
 * Get the ids of both vertices of the edge.
 */
array<unsigned long> triangulation::edge_vertex_ids(
   unsigned long e_id) const
{
   this->check_edge_id(e_id);
   sym_edge<vrtx,edge>& e = (*_edges)[e_id];
   array<unsigned long> v_ids(2);
   v_ids[0] = e.origin().id;
   v_ids[1] = e.destination().id;
   return v_ids;
}

/*
 * Get the ids of both triangles adjacent to the edge.
 */
array<unsigned long> triangulation::edge_triangle_ids(
   unsigned long e_id) const
{
   this->check_edge_id(e_id);
   sym_edge<vrtx,edge>& e = (*_edges)[e_id];
   tri* t_left  = e.data().t_left;
   tri* t_right = e.sym().data().t_left;
   array<unsigned long> t_ids(2);
   t_ids[0] = (t_left != NULL)  ? (t_left->id)  : (_triangles->size());
   t_ids[1] = (t_right != NULL) ? (t_right->id) : (_triangles->size());
   return t_ids;
}

/***************************************************************************
 * Triangle -> vertex, edge lookup.
 ***************************************************************************/

/*
 * Get the id of the specified vertex of the triangle.
 */
unsigned long triangulation::triangle_vertex_id(
   unsigned long t_id, unsigned long v_num) const
{
   this->check_triangle_id(t_id);
   sym_edge<vrtx,edge>* e = (*_triangles)[t_id].e;
   if (v_num == 0)
      return e->origin().id;
   else if (v_num == 1)
      return e->left_next().origin().id;
   else if (v_num == 2)
      return e->left_prev().origin().id;
   else
      throw ex_index_out_of_bounds(
         "triangle vertex number must be 0, 1, or 2", v_num
      );
}

/*
 * Get the id of the specified edge of the triangle.
 */
unsigned long triangulation::triangle_edge_id(
   unsigned long t_id, unsigned long e_num) const
{
   this->check_triangle_id(t_id);
   sym_edge<vrtx,edge>* e = (*_triangles)[t_id].e;
   if (e_num == 0)
      return e->left_next().data().id;
   else if (e_num == 1)
      return e->left_prev().data().id;
   else if (e_num == 2)
      return e->data().id;
   else
      throw ex_index_out_of_bounds(
         "triangle edge number must be 0, 1, or 2", e_num
      );
}

/*
 * Get the ids of the vertices of the triangle.
 */
array<unsigned long> triangulation::triangle_vertex_ids(
   unsigned long t_id) const
{
   this->check_triangle_id(t_id);
   sym_edge<vrtx,edge>* e = (*_triangles)[t_id].e;
   array<unsigned long> v_ids(3);
   v_ids[0] = e->origin().id;
   v_ids[1] = e->left_next().origin().id;
   v_ids[2] = e->left_prev().origin().id;
   return v_ids;
}
   
/*
 * Get the ids of the edges of the triangle.
 */
array<unsigned long> triangulation::triangle_edge_ids(
   unsigned long t_id) const
{
   this->check_triangle_id(t_id);
   sym_edge<vrtx,edge>* e = (*_triangles)[t_id].e;
   array<unsigned long> e_ids(3);
   e_ids[0] = e->left_next().data().id;
   e_ids[1] = e->left_prev().data().id;
   e_ids[2] = e->data().id;
   return e_ids;
}

/***************************************************************************
 * Adjacent triangle lookup.
 ***************************************************************************/

/*
 * Get the id of the specified adjacent (sharing an edge) triangle.
 */
unsigned long triangulation::adjacent_triangle_id(
   unsigned long t_id, unsigned long t_num) const
{
   this->check_triangle_id(t_id);
   sym_edge<vrtx,edge>* e = (*_triangles)[t_id].e;
   if (t_num == 0)
      e = &(e->left_next());
   else if (t_num == 1)
      e = &(e->left_prev());
   else if (t_num != 2)
      throw ex_index_out_of_bounds(
         "adjacent triangle number must be 0, 1, or 2", t_num
      );
   tri* adjacent_t = e->sym().data().t_left;
   return (adjacent_t != NULL) ? (adjacent_t->id) : (_triangles->size());
}

/*
 * Get the ids of all adjacent (sharing an edge) triangles.
 */
array<unsigned long> triangulation::adjacent_triangle_ids(
   unsigned long t_id) const
{
   this->check_triangle_id(t_id);
   sym_edge<vrtx,edge>* e = (*_triangles)[t_id].e;
   tri* adjacent_t0 = e->left_next().sym().data().t_left;
   tri* adjacent_t1 = e->left_prev().sym().data().t_left;
   tri* adjacent_t2 = e->sym().data().t_left;
   array<unsigned long> t_ids(3);
   t_ids[0] = (adjacent_t0 != NULL) ? (adjacent_t0->id) : (_triangles->size());
   t_ids[1] = (adjacent_t1 != NULL) ? (adjacent_t1->id) : (_triangles->size());
   t_ids[2] = (adjacent_t2 != NULL) ? (adjacent_t2->id) : (_triangles->size());
   return t_ids;
}

/***************************************************************************
 * Vertex and triangle reference.
 ***************************************************************************/

/*
 * Return the specified vertex (by reference).
 */
point_2D& triangulation::vertex(unsigned long v_id) const {
   this->check_vertex_id(v_id);
   return (*_vertices)[v_id].p;
}

/*
 * Return the specified triangle (containing vertices by reference).
 */
triangle_2D triangulation::triangle(unsigned long t_id) const {
   array<unsigned long> v_ids = this->triangle_vertex_ids(t_id);
   return triangle_2D(
      (*_vertices)[v_ids[0]].p,
      (*_vertices)[v_ids[1]].p,
      (*_vertices)[v_ids[2]].p
   );
}

/***************************************************************************
 * Primary and dual graphs.
 ***************************************************************************/

/*
 * Return the undirected graph corresponding to the triangulation.
 *
 * For each vertex, return an array indicating the vertices of lower id
 * to which an edge connects it.
 */
auto_collection< array<unsigned long>, array_list< array<unsigned long> > >
   triangulation::graph() const
{
   auto_collection< array<unsigned long>, array_list< array<unsigned long> > >
      g(new array_list< array<unsigned long> >());
   unsigned long n_v = _vertices->size();
   for (unsigned long v_id = 0; v_id < n_v; v_id++) {
      /* get edge out of vertex */
      sym_edge<vrtx,edge>* e = (*_vertices)[v_id].e;
      /* count connected vertices with lower id */
      unsigned long n_e = 0;
      if (e != NULL) {
         sym_edge<vrtx,edge>* e_curr = e;
         do {
            if (e_curr->destination().id < v_id)
               n_e++;
            e_curr = &(e_curr->origin_next());
         } while (e_curr != e);
      }
      /* retrieve connected vertices with lower id */
      auto_ptr< array<unsigned long> > e_v_ids(new array<unsigned long>(n_e));
      n_e = 0;
      if (e != NULL) {
         sym_edge<vrtx,edge>* e_curr = e;
         do {
            unsigned long e_v_id = e_curr->destination().id;
            if (e_v_id < v_id)
               (*e_v_ids)[n_e++] = e_v_id;
            e_curr = &(e_curr->origin_next());
         } while (e_curr != e);
      }
      /* store vertex id array */
      g->add(*e_v_ids);
      e_v_ids.release();
   }
   return g;
}

/*
 * Return the undirected dual graph of the triangulation.
 *
 * For each triangle, return an array indicating the triangles of lower id
 * with which it shares an edge.
 */
auto_collection< array<unsigned long>, array_list< array<unsigned long> > >
   triangulation::graph_dual() const
{
   auto_collection< array<unsigned long>, array_list< array<unsigned long> > >
      g(new array_list< array<unsigned long> >());
   unsigned long n_t = _triangles->size();
   for (unsigned long t_id = 0; t_id < n_t; t_id++) {
      /* count adjacent triangles with lower id */
      array<unsigned long> t_ids_adjacent = this->adjacent_triangle_ids(t_id);
      unsigned long n_adjacent = 0;
      for (unsigned long n = 0; n < 3; n++) {
         if (t_ids_adjacent[n] < t_id)
            n_adjacent++;
      }
      /* retrieve adjacent triangles with lower id*/
      auto_ptr< array<unsigned long> > t_ids(
         new array<unsigned long>(n_adjacent)
      );
      n_adjacent = 0;
      for (unsigned long n = 0; n < 3; n++) {
         unsigned long id = t_ids_adjacent[n];
         if (id < t_id)
            (*t_ids)[n_adjacent++] = id;
      }
      /* store triangle id array */
      g->add(*t_ids);
      t_ids.release();
   }
   return g;
}  

/***************************************************************************
 * Argument checking helper functions.
 ***************************************************************************/

/*
 * Check that the given vertex id is valid.
 * Throw an exception (ex_index_out_of_bounds) if it is invalid.
 */
void triangulation::check_vertex_id(unsigned long v_id) const {
   if (v_id >= _vertices->size())
      throw ex_index_out_of_bounds("invalid vertex id", v_id);
}

/*
 * Check that the given edge id is valid.
 * Throw an exception (ex_index_out_of_bounds) if it is invalid.
 */
void triangulation::check_edge_id(unsigned long e_id) const {
   if (e_id >= _edges->size())
      throw ex_index_out_of_bounds("invalid edge id", e_id);
}

/*
 * Check that the given triangle id is valid.
 * Throw an exception (ex_index_out_of_bounds) if it is invalid.
 */
void triangulation::check_triangle_id(unsigned long t_id) const {
   if (t_id >= _triangles->size())
      throw ex_index_out_of_bounds("invalid triangle id", t_id);
}

/***************************************************************************
 * Vertex, edge, and triangle classes.
 ***************************************************************************/

/*
 * Vertex constructor.
 */
triangulation::vrtx::vrtx(unsigned long id, point_2D& p)
 : id(id), p(p), e(NULL)
{ }

/*
 * Edge constructors.
 */
triangulation::edge::edge(unsigned long id)
 : id(id), t_left(NULL)
{ }

triangulation::edge::edge(unsigned long id, tri* t_left)
 : id(id), t_left(t_left)
{ }

/*
 * Triangle node constructors.
 */
triangulation::tri::tri(unsigned long id)
 : id(id), e(NULL)
{ }

triangulation::tri::tri(unsigned long id, sym_edge<vrtx,edge>* e)
 : id(id), e(e)
{ }

/***************************************************************************
 * Triangulation helper functions.
 ***************************************************************************/
 
/*
 * Return true iff the point is to the left of the edge.
 */
bool triangulation::left_of(const point_2D& p, const sym_edge<vrtx,edge>& e) {
   return (p.orientation(e.origin().p, e.destination().p) > 0);
}

/*
 * Return true iff the vertex is to the left of the edge.
 */
bool triangulation::left_of(const vrtx& v, const sym_edge<vrtx,edge>& e) {
   return left_of(v.p, e);
}

/*
 * Return true iff the point is to the right of the edge.
 */
bool triangulation::right_of(const point_2D& p, const sym_edge<vrtx,edge>& e) {
   return (p.orientation(e.origin().p, e.destination().p) < 0);
}

/*
 * Return true iff the vertex is to the right of the edge.
 */
bool triangulation::right_of(const vrtx& v, const sym_edge<vrtx,edge>& e) {
   return right_of(v.p, e);
}

/*
 * Return true iff the first edge is above the second base edge.
 */
bool triangulation::above_base_edge(
   const sym_edge<vrtx,edge>& e, const sym_edge<vrtx,edge>& base)
{
   return (
      e.destination().p.orientation(base.destination().p, base.origin().p) > 0
   );
}

/*
 * Create the constrained Delaunay triangulation of the polygon located to 
 * the left of the given edge.  Consume the given edges in forming the
 * triangulation.
 */
void triangulation::triangulate_polygon(
   sym_edge<vrtx,edge>&         e_base,
   list< sym_edge<vrtx,edge> >& edges)
{
   /* get endpoints of base edge */
   vrtx& v_start = e_base.origin();
   vrtx& v_end   = e_base.destination();
   /* find first polygon vertex hit by circle containing base edge as chord */   
   sym_edge<vrtx,edge>* e_hit_r = &(e_base.left_next());
   sym_edge<vrtx,edge>* e = &(e_hit_r->left_next());
   vrtx* v_hit = &(e_hit_r->destination());
   vrtx* v = &(e->destination());
   while (v != &v_start) {
      if (v->p.in_circle(v_start.p, v_end.p, v_hit->p) > 0) {
         e_hit_r = e;
         v_hit = v;
      }
      e = &(e->left_next());
      v = &(e->destination());
   }
   sym_edge<vrtx,edge>* e_hit_l = &(e_hit_r->left_next());
   /* update vertex -> edge pointer */
   v_hit->e = e_hit_l;
   /* add edge from v_hit to v_end, retriangulate its polygon */
   if (&(e_hit_r->origin()) != &v_end) {
      /* add edge */
      sym_edge<vrtx,edge>& e_r = edges.remove_head();
      e_r.reconnect(*e_hit_r, e_base.left_next());
      /* retriangulate */
      triangulation::triangulate_polygon(e_r, edges);
   }
   /* add edge from v_start to v_hit, retriangulate its polygon */
   if (&(e_hit_l->destination()) != &v_start) {
      /* add edge */
      sym_edge<vrtx,edge>& e_l = edges.remove_head();
      e_l.reconnect(e_base.left_prev(), *e_hit_l);
      /* retriangulate */
      triangulation::triangulate_polygon(e_l, edges);
   }
}

/*
 * Insert a constraint edge between the vertices with the given ids.
 */
void triangulation::insert_edge(
   unsigned long v_start_id, unsigned long v_end_id, unsigned long e_id)
{
   /* get vertices at endpoints of edge */
   vrtx& v_start = (*_vertices)[v_start_id];
   vrtx& v_end   = (*_vertices)[v_end_id];
   /* find edge out of v_start which is left of v_end (if it exists) */ 
   sym_edge<vrtx,edge>* e_left = v_start.e;
   while (triangulation::right_of(v_end, *e_left)) {
      e_left = &(e_left->origin_prev());
      if (e_left == v_start.e)
         break;
   }
   while (!(triangulation::right_of(v_end, *e_left))) {
      e_left = &(e_left->origin_next());
      if (e_left == v_start.e)
         break;
   }
   /* find edge out of v_start which is right of v_end (if it exists) */
   sym_edge<vrtx,edge>* e_right = v_start.e;
   while (triangulation::left_of(v_end, *e_right)) {
      e_right = &(e_right->origin_next());
      if (e_right == v_start.e)
         break;
   }
   while (!(triangulation::left_of(v_end, *e_right))) {
      e_right = &(e_right->origin_prev());
      if (e_right == v_start.e)
         break;
   }
   /* check left/right edges */
   sym_edge<vrtx,edge>* e = NULL;
   if (triangulation::right_of(v_end, *e_left)) {
      e_right = &(e_left->origin_prev());
      e = e_right;
   } else if (triangulation::left_of(v_end, *e_right)) {
      e_left = &(e_right->origin_next());
      e = e_left;
   } else {
      /* degenerate triangulation - all vertices are colinear */
      e = e_left;
      if (&(e->destination()) != &v_end)
         e = &(e->origin_next());
   }
   /* check if constraint edges already exists */
   if (&(e->destination()) == &v_end) {
      /* check for duplicate constraint */
      unsigned long n_e = _edges->size();
      if (e->data().id < n_e)
         throw ex_invalid_argument("duplicate constraint edges specified");
      /* mark edge as constrained */
      e->data().id       = e_id;
      e->sym().data().id = e_id;
   } else if (v_end.p.orientation(e->origin().p, e->destination().p) == 0) {
       /* v_end is collinear with current edge - invalid constraint */
      throw ex_invalid_argument(
         "constraint edge specified to pass through 3 collinear vertices"
      );
   } else {
      /* initialize list of removed edges */
      list< sym_edge<vrtx,edge> > removed_edges;
      /* remove edges that cross the constraint edge */
      sym_edge<vrtx,edge>* e_l = e_left;
      sym_edge<vrtx,edge>* e_r = e_right;
      unsigned long n_e = _edges->size();
      while (true) {
         /* remove edge */
         e = &(e_r->left_next());
         if (e->data().id < n_e)
            throw ex_invalid_argument("constraint edges cross");
         e->disconnect();
         removed_edges.add(*e);
         /* update left/right bounding edges */
         sym_edge<vrtx,edge>* e_l_next = &(e_l->right_prev());
         sym_edge<vrtx,edge>* e_r_next = &(e_r->left_next());
         vrtx& v_next = e_l_next->destination();
         int orient = v_next.p.orientation(v_start.p, v_end.p);
         if (orient > 0) {
            /* move the left edge */
            e_l = e_l_next;
         } else if (orient < 0) {
            /* move the right edge */
            e_r = e_r_next;
         } else if (&v_next == &v_end) {
            /* done - move both edges */
            e_l = e_l_next;
            e_r = e_r_next;
            break;
         } else {
            /* v_end is collinear with current edge - invalid constraint */
            throw ex_invalid_argument(
               "constraint edge specified to pass through 3 collinear vertices"
            );
         }
      }
      /* insert the constraint edge */
      sym_edge<vrtx,edge>& e_base = removed_edges.remove_head();
      e_base.reconnect(e_left->sym(), e_l->sym());
      /* mark edge as constrained */
      e_base.data().id       = e_id;
      e_base.sym().data().id = e_id;
      /* update vertex -> edge pointers */
      v_start.e = &e_base;
      v_end.e = &(e_base.sym());
      /* retriangulate left/right bounding polygons */
      triangulation::triangulate_polygon(e_base, removed_edges);
      triangulation::triangulate_polygon(e_base.sym(), removed_edges);
   }
}

/*
 * Build the vertices from the given collection of points.
 */
void triangulation::build_vertices(const collection<point_2D>& points) {
   auto_ptr< iterator<point_2D> > i = points.iter_create();
   unsigned long n_v = 0;
   while (i->has_next()) {
      point_2D& p = i->next();
      auto_ptr<vrtx> v(new vrtx(n_v++, p));
      _vertices->add(*v);
      v.release();
   }
}

/*
 * Link each vertex to an edge out of it if such an edge exists.
 * The vertices and edges must already exist.
 */
void triangulation::build_vertex_edge_links() {
   unsigned long n_e = _edges->size();
   for (unsigned long n = 0; n < n_e; n++) {
      /* get both sides of edge */
      sym_edge<vrtx,edge>& e     = (*_edges)[n];
      sym_edge<vrtx,edge>& e_sym = e.sym();
      /* link currently unlinked origin vertices */
      vrtx& e_org     = e.origin();
      vrtx& e_sym_org = e_sym.origin();
      if (e_org.e == NULL)
         e_org.e = &e;
      if (e_sym_org.e == NULL)
         e_sym_org.e = &e_sym;
   }
}

/*
 * Build triangle node for the triangle to the left of the given edge.
 * Link the triangle to the edge structure.
 */
void triangulation::build_triangle(sym_edge<vrtx,edge>& e) {
   /* skip edge if already linked */
   edge& e_data = e.data();
   if (e_data.t_left == NULL) {
      /* check that a valid triangle lies left of the edge */
      sym_edge<vrtx,edge>& e_onext = e.origin_next();
      sym_edge<vrtx,edge>& e_lnext = e.left_next();
      vrtx& v = e_onext.destination();
      if ((&v == &(e_lnext.destination())) && triangulation::left_of(v, e)) {
         /* get edge data for other sides of triangle */
         edge& e_lnext_data = e_lnext.data();
         edge& e_lprev_data = e.left_prev().data();
         /* create triangle node */
         auto_ptr<tri> t(new tri(_triangles->size(), &e));
         _triangles->add(*t);
         tri* t_ptr = t.release();
         /* link edges to triangle */
         e_data.t_left       = t_ptr;
         e_lnext_data.t_left = t_ptr;
         e_lprev_data.t_left = t_ptr;
      }
   }
}

/*
 * Build triangle nodes and link them to the edge structure.
 * The vertices, edges, and edge data must already exist.
 */
void triangulation::build_triangles() {
   unsigned long n_e = _edges->size();
   for (unsigned long n = 0; n < n_e; n++) {
      /* get both sides of edge */
      sym_edge<vrtx,edge>& e     = (*_edges)[n];
      sym_edge<vrtx,edge>& e_sym = e.sym();
      /* link both sides to triangles (if needed) */
      this->build_triangle(e);
      this->build_triangle(e_sym);
   }
}

/***************************************************************************
 * Delaunay triangulation.
 ***************************************************************************/

/*
 * Constructor.
 * The input vertex set must be sorted lexicographically.
 */
triangulation::delaunay_triangulator::delaunay_triangulator(
   const array_list<vrtx>&           vertices,
   unsigned long                     start,
   unsigned long                     end,
   auto_collection< 
      sym_edge<vrtx,edge>,
      list< sym_edge<vrtx,edge> > >& edges)
 : _vertices(vertices),
   _start(start),
   _end(end),
   _edges(edges),
   _le(NULL),
   _re(NULL)
{ }

/*
 * Destructor.
 */
triangulation::delaunay_triangulator::~delaunay_triangulator() {
   /* do nothing */
}

/*
 * Compute Delaunay triangulation of vertices in subrange.
 */
void triangulation::delaunay_triangulator::run() {
   /* check size of vertices set */
   unsigned long n_vertices = (_start <= _end) ? (_end - _start + 1) : 0;
   if (n_vertices == 2) {
      /* initialize edge list */
      _edges.reset(new list< sym_edge<vrtx,edge> >());
      /* create an edge */
      auto_ptr< sym_edge<vrtx,edge> > e(
         new sym_edge<vrtx,edge>(_vertices[_start], _vertices[_end])
      );
      _edges->add(*e);
      _le = e.release();
      _re = &(_le->sym());
   } else if (n_vertices == 3) {
      /* initialize edge list */
      _edges.reset(new list< sym_edge<vrtx,edge> >());
      /* create a triangle */
      vrtx& v0 = _vertices[_start];
      vrtx& v1 = _vertices[_start+1];
      vrtx& v2 = _vertices[_end];
      /* create edge a connecting v0 to v1, and b connecting v1 to v2 */
      auto_ptr< sym_edge<vrtx,edge> > a(new sym_edge<vrtx,edge>(v0, v1));
      auto_ptr< sym_edge<vrtx,edge> > b(new sym_edge<vrtx,edge>(v1, v2));
      sym_edge<vrtx,edge>::splice(a->sym(), *b);
      /* close the triangle */
      if (v2.p.orientation(v0.p, v1.p) > 0) {
         /* v0, v1, v2 are in counterclockwise order */
         auto_ptr< sym_edge<vrtx,edge> > c =
            sym_edge<vrtx,edge>::connect(*b, *a);
         _edges->add(*c);
         c.release();
         _le = a.get();
         _re = &(b->sym());
      } else if (v1.p.orientation(v0.p, v2.p) > 0) {
         /* v0, v2, v1 are in counterclockwise order */
         auto_ptr< sym_edge<vrtx,edge> > c =
            sym_edge<vrtx,edge>::connect(*b, *a);
         _edges->add(*c);
         _re = c.release();
         _le = &(_re->sym());
      } else {
         /* v0, v1, v2 are collinear */
         _le = a.get();
         _re = &(b->sym());
      }
      _edges->add(*a);  a.release();
      _edges->add(*b);  b.release();
   } else if (n_vertices > 3) {
      /* recursively compute triangulation of each half */
      auto_collection< sym_edge<vrtx,edge>, list < sym_edge<vrtx,edge> > > erht;
      unsigned long mid = (_start + _end)/2;
      delaunay_triangulator dt_left(_vertices, _start, mid, _edges);
      delaunay_triangulator dt_right(_vertices, mid+1, _end, erht);
      child_thread::run(dt_left, dt_right);
      /* combine edge lists */
      unsigned long n_erht = erht->size();
      for (unsigned long n = 0; n < n_erht; n++) {
         auto_ptr< sym_edge<vrtx,edge> > e(&(erht->remove_head()));
         _edges->add(*e);
         e.release();
      }
      /* initialize edge pointers */
      sym_edge<vrtx,edge>* ldo = dt_left._le;
      sym_edge<vrtx,edge>* ldi = dt_left._re;
      sym_edge<vrtx,edge>* rdi = dt_right._le;
      sym_edge<vrtx,edge>* rdo = dt_right._re;
      /* compute the lower common tangent of left and right sets */
      while (true) {
         if (triangulation::left_of(rdi->origin(), *ldi))
            ldi = &(ldi->left_next());
         else if (triangulation::right_of(ldi->origin(), *rdi))
            rdi = &(rdi->right_prev());
         else
            break;
      }
      /* create a cross edge basel from rdi->origin() to ldi->origin() */
      auto_ptr< sym_edge<vrtx,edge> > auto_basel =
         sym_edge<vrtx,edge>::connect(rdi->sym(), *ldi);
      _edges->add(*auto_basel);
      sym_edge<vrtx,edge>* basel = auto_basel.release();
      if (&(ldi->origin()) == &(ldo->origin())) { ldo = &(basel->sym()); }
      if (&(rdi->origin()) == &(rdo->origin())) { rdo = basel; }
      /* create edge data object used to mark deleted edges */
      auto_ptr<edge> e_data_deleted(new edge(0));
      /* merge loop */
      while (true) {
         /*
          * locate the first left point (lcand->destination()) to be 
          * encountered by the rising bubble and delete left edges out of
          * (basel->destination()) that fail the in_circle test
          */
         sym_edge<vrtx,edge>* lcand = &(basel->sym().origin_next());
         if (triangulation::above_base_edge(*lcand, *basel)) {
            while (lcand->origin_next().destination().p.in_circle(
               basel->destination().p,
               basel->origin().p,
               lcand->destination().p) > 0)
            {
               sym_edge<vrtx,edge>* t = &(lcand->origin_next());
               lcand->disconnect();
               lcand->data(*e_data_deleted);
               lcand->sym().data(*e_data_deleted);
               lcand = t;
            }
         }
         /*
          * symmetrically, locate the first right point to be hit, and delete
          * right edges that fail the in_circle test
          */
         sym_edge<vrtx,edge>* rcand = &(basel->origin_prev());
         if (triangulation::above_base_edge(*rcand, *basel)) {
            while (rcand->origin_prev().destination().p.in_circle(
               basel->destination().p,
               basel->origin().p,
               rcand->destination().p) > 0)
            {
               sym_edge<vrtx,edge>* t = &(rcand->origin_prev());
               rcand->disconnect();
               rcand->data(*e_data_deleted);
               rcand->sym().data(*e_data_deleted);
               rcand = t;
            }
         }
         /* check if basel is the upper common tangent */
         bool l_valid = triangulation::above_base_edge(*lcand, *basel);
         bool r_valid = triangulation::above_base_edge(*rcand, *basel);
         if ((!l_valid) && (!r_valid))
            break;
         /*
          * the next cross edge is connected to either:
          *    lcand->destination() or rcand->destination()
          * if both are valid, use the in_circle test to choose one
          */
         if ((!l_valid) || 
             (r_valid && 
             (rcand->destination().p.in_circle(
               lcand->destination().p,
               lcand->origin().p,
               rcand->origin().p) > 0)))
         {
            /* add edge from rcand->destination() to basel->destination() */
            auto_basel = sym_edge<vrtx,edge>::connect(
               *rcand, basel->sym()
            );
            _edges->add(*auto_basel);
            basel = auto_basel.release();
         } else {   
            /* add edge from basel->origin() to lcand->destination() */
            auto_basel = sym_edge<vrtx,edge>::connect(
               basel->sym(), lcand->sym()
            );
            _edges->add(*auto_basel);
            basel = auto_basel.release();
         }
      } /* end of merge loop */
      _le = ldo;
      _re = rdo;
      /* remove deleted edges from edge list */
      unsigned long n_edges = _edges->size();
      for (unsigned long n = 0; n < n_edges; n++) {
         auto_ptr< sym_edge<vrtx,edge> > e(&(_edges->remove_head()));
         if (!(e->has_data())) {
            _edges->add(*e);
            e.release();
         }
      }
   }
}

} /* namespace geometry */
} /* namespace math */
