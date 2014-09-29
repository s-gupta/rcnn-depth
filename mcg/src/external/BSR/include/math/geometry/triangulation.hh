/*
 * Triangulation.
 */
#ifndef MATH__GEOMETRY__TRIANGULATION_HH
#define MATH__GEOMETRY__TRIANGULATION_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/runnable.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/sym_edge.hh"
#include "math/geometry/triangle_2D.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using concurrent::threads::runnable;
using lang::array;
using lang::pointers::auto_ptr;

/*
 * Triangulation of a set of points in the plane.
 * Points are stored by reference in the triangulation.
 */
class triangulation {
public:
   /*
    * Constructor.
    * Create an empty triangulation.
    */
   triangulation();

   /*
    * Copy constructor.
    */
   triangulation(const triangulation&);

   /*
    * Destructor.
    */
   virtual ~triangulation();

   /*
    * Delaunay triangulation.
    * Create a Delaunay triangulation of the given set of unique points.
    *
    * The O(n + (n/p)*log(n/p)) divide and conquer algorithm is used to
    * construct the Delaunay triangulation, where n is the number of points
    * and p is the number of available processors.
    */
   static auto_ptr<triangulation> delaunay(const collection<point_2D>&);

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
   static auto_ptr<triangulation> delaunay(
      const collection<point_2D>&,  /* points */
      const array<unsigned long>&,  /* constraint edge first endpoint ids */
      const array<unsigned long>&   /* constraint edge second endpoint ids */
   );

   /*
    * Get number of vertices in triangulation.
    */
   unsigned long vertices_size() const;
   
   /*
    * Get number of edges in triangulation.
    */
   unsigned long edges_size() const;

   /*
    * Get number of triangles in triangulation.
    */
   unsigned long triangles_size() const;

   /*
    * Get the number of edges/triangles containing the specified vertex.
    */
   unsigned long vertex_edges_size(unsigned long /* vertex id */) const;

   unsigned long vertex_triangles_size(unsigned long /* vertex id */) const;
   
   /*
    * Vertex, edge, and triangle ids.
    *
    * Let Nv, Ne, and Nt be the number of vertices, edges, and triangles
    * respectively.  Vertex, edge, and triangle ids are integers in the
    * range [0,Nv-1], [0,Ne-1], and [0,Nt-1], respectively.
    *
    * In addition, vertices are numbered in the same order as they appear in
    * the initial collection of points.  Constaint edges appear first in the
    * set of edges and are numbered in the same order as they are found in the
    * initial set of constraints.
    *
    * The exterior region (not contained in any triangle) is assigned the
    * special triangle id Nt.  This id is only returned when accessing the
    * triangle on the exterior side of a boundary edge.  This id cannot be 
    * used as an argument to any method requiring a triangle id.
    *
    * Within a triangle, vertices and edges are numbered so that vertex n is
    * the vertex opposite edge n, for n = 0, 1, 2.
    */
    
   /*
    * Get the id(s) of the specified edge(s)/triangle(s) containing the vertex.
    */
   unsigned long vertex_edge_id(
      unsigned long,    /* vertex id */
      unsigned long     /* edge # */
   ) const;

   unsigned long vertex_triangle_id(
      unsigned long,    /* vertex id */
      unsigned long     /* triangle # */
   ) const;
   
   array<unsigned long> vertex_edge_ids(
      unsigned long     /* vertex id */
   ) const;

   array<unsigned long> vertex_triangle_ids(
      unsigned long     /* vertex id */
   ) const;
    
   /*
    * Get the id(s) of the specified component(s) of an edge.
    */
   unsigned long edge_vertex_id(
      unsigned long,    /* edge id */
      unsigned long     /* vertex # (0 or 1) */
   ) const;

   unsigned long edge_triangle_id(
      unsigned long,    /* edge id */
      unsigned long     /* triangle # (0 or 1) */
   ) const;

   array<unsigned long> edge_vertex_ids(
      unsigned long     /* edge id */
   ) const;

   array<unsigned long> edge_triangle_ids(
      unsigned long     /* edge id */
   ) const;

   /*
    * Get the id(s) of the specified component(s) of a triangle.
    */
   unsigned long triangle_vertex_id(
      unsigned long,    /* triangle id */
      unsigned long     /* vertex # (0, 1, or 2) */
   ) const;

   unsigned long triangle_edge_id(
      unsigned long,    /* triangle id */
      unsigned long     /* edge # (0, 1, or 2) */
   ) const;

   array<unsigned long> triangle_vertex_ids(
      unsigned long     /* triangle id */
   ) const;
      
   array<unsigned long> triangle_edge_ids(
      unsigned long     /* triangle id */
   ) const;
   
   /*
    * Get the id of the specified adjacent (sharing an edge) triangle.
    */
   unsigned long adjacent_triangle_id(
      unsigned long,    /* triangle id */
      unsigned long     /* triangle # (0, 1, or 2) */
   ) const;

   /*
    * Get the ids of all adjacent (sharing an edge) triangles.
    */
   array<unsigned long> adjacent_triangle_ids(
      unsigned long     /* triangle id */
   ) const;
   
   /*
    * Return the specified vertex (by reference).
    */
   point_2D& vertex(unsigned long /* vertex id */) const;

   /*
    * Return the specified triangle (containing vertices by reference).
    */
   triangle_2D triangle(unsigned long /* triangle id */) const;

   /*
    * Return the undirected graph corresponding to the triangulation.
    *
    * For each vertex, return an array indicating the vertices of lower id
    * to which an edge connects it.
    */
   auto_collection< array<unsigned long>,
      array_list< array<unsigned long> > > graph() const;

   /*
    * Return the undirected dual graph of the triangulation.
    *
    * For each triangle, return an array indicating the triangles of lower id
    * with which it shares an edge.
    */
   auto_collection< array<unsigned long>,
      array_list< array<unsigned long> > > graph_dual() const;
   
protected:
   /************************************************************************
    * Triangulation data structures.
    ************************************************************************/
  
   /*
    * Declare vertex, edge, and triangle types.
    */
   class vrtx;
   class edge;
   class tri;
    
   /*
    * Vertex.
    */
   class vrtx {
   public:
      /* constructor */
      explicit vrtx(unsigned long /* id */, point_2D& /* point */);

      /* data */
      unsigned long        id;   /* vertex id */
      point_2D&            p;    /* point at vertex */
      sym_edge<vrtx,edge>* e;    /* edge out of vertex (if it exists) */
   };

   /*
    * Edge data.
    */
   class edge {
   public:
      /* constructors */
      explicit edge(unsigned long /* id */);
      explicit edge(unsigned long /* id */, tri* /* triangle to left */);

      /* data */
      unsigned long id;          /* undirected edge id */
      tri*          t_left;      /* triangle to left of edge (if it exists) */
   };

   /*
    * Triangle node.
    */
   class tri {
   public:
      /* constructors */
      explicit tri(unsigned long /* id */);
      explicit tri(unsigned long /* id */, sym_edge<vrtx,edge>* /* edge */);
      
      /* data */
      unsigned long id;          /* triangle id */
      sym_edge<vrtx,edge>* e;    /* edge which triangle is left of */
   };
      
   /*
    * Vertex, edge, and triangle arrays.
    */
   auto_collection< vrtx, array_list<vrtx> > _vertices;     /* vertices */
   auto_collection< 
      sym_edge<vrtx,edge>,
      array_list< sym_edge<vrtx,edge> > >    _edges;        /* edges */
   auto_collection< edge, list<edge> >       _edge_data;    /* edge data */
   auto_collection< tri, array_list<tri> >   _triangles;    /* triangles */

   /************************************************************************
    * Argument checking helper functions.
    ************************************************************************/

   /*
    * Check that the given vertex, edge, or triangle id is valid.
    * Throw an exception (ex_index_out_of_bounds) if it is invalid.
    */
   void check_vertex_id(unsigned long /* vertex id */) const;
   void check_edge_id(unsigned long /* edge id */) const;
   void check_triangle_id(unsigned long /* triangle id */) const;

   /************************************************************************
    * Triangulation helper functions.
    ************************************************************************/

   /*
    * Return true iff the point/vertex is to the left of the edge.
    */
   static bool left_of(const point_2D&, const sym_edge<vrtx,edge>&);
   static bool left_of(const vrtx&,     const sym_edge<vrtx,edge>&);

   /*
    * Return true iff the point/vertex is to the right of the edge.
    */
   static bool right_of(const point_2D&, const sym_edge<vrtx,edge>&);
   static bool right_of(const vrtx&,     const sym_edge<vrtx,edge>&);

   /*
    * Return true iff the first edge is above the second base edge.
    */
   static bool above_base_edge(
      const sym_edge<vrtx,edge>&, const sym_edge<vrtx,edge>&
   );

   /*
    * Create the constrained Delaunay triangulation of the polygon located to 
    * the left of the given edge.  Consume the given edges in forming the
    * triangulation.
    */
   static void triangulate_polygon(
      sym_edge<vrtx,edge>&,            /* edge which polygon is left of */
      list< sym_edge<vrtx,edge> >&     /* edges to use in triangulation */
   );
 
   /*
    * Insert a constraint edge between the vertices with the given ids.
    */
   void insert_edge(
      unsigned long,    /* vertex id */
      unsigned long,    /* vertex id */
      unsigned long     /* constraint edge id */
   );
  
   /*
    * Build the vertices from the given collection of points.
    */
   void build_vertices(const collection<point_2D>&);
  
   /*
    * Link each vertex to an edge out of it if such an edge exists.
    * The vertices and edges must already exist.
    */
   void build_vertex_edge_links();

   /*
    * Build triangle node for the triangle to the left of the given edge.
    * Link the triangle to the edge structure.
    */
   void build_triangle(sym_edge<vrtx,edge>& e);
   
   /*
    * Build triangle nodes and link them to the edge structure.
    * The vertices, edges, and edge data must already exist.
    */
   void build_triangles();
   
   /************************************************************************
    * Delaunay triangulation.
    ************************************************************************/
   
   /*
    * Runnable object for computing the Delaunary triangulation using the 
    * O(n + (n/p)*log(n/p)) time divide and conquer algorithm, where n is the 
    * number of vertices and p is the number of available processors.
    */
   class delaunay_triangulator : public runnable {
   public:
      /*
       * Constructor.
       * The input vertex set must be sorted lexicographically.
       */
      explicit delaunay_triangulator(
         const array_list<vrtx>&,                  /* vertex set */
         unsigned long,                            /* start of subrange */
         unsigned long,                            /* end of subrange */
         auto_collection< 
            sym_edge<vrtx,edge>,
            list< sym_edge<vrtx,edge> > >&         /* returned edges */
      );

      /*
       * Destructor.
       */
      virtual ~delaunay_triangulator();

      /*
       * Compute Delaunay triangulation of vertices in subrange.
       */
      virtual void run();

   protected:
      const array_list<vrtx>&           _vertices; /* vertex set */
      unsigned long                     _start;    /* start of subrange */
      unsigned long                     _end;      /* end of subrange */
      auto_collection< 
         sym_edge<vrtx,edge>,
         list< sym_edge<vrtx,edge> > >& _edges;    /* returned edges */
      sym_edge<vrtx,edge>* _le;  /* ccw hull edge out of leftmost vertex */
      sym_edge<vrtx,edge>* _re;  /* cw hull edge out of rightmost vertex */
   };
};

} /* namespace geometry */
} /* namespace math */

#endif
