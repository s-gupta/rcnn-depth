/*
 * Agglomerative cluterer on undirected graphs.
 *
 * Each step of the agglomeration merges the vertices at the ends of the
 * highest priority edge and updates the costs of the edges connected to 
 * the resulting vertex.
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__GRAPH__BASIC_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__GRAPH__BASIC_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/queue_set.hh"
#include "collections/set.hh"
#include "functors/comparable_functors.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "mlearning/clustering/clusterers/abstract/clusterer.hh"

namespace mlearning {
namespace clustering {
namespace clusterers {
namespace graph {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using collections::queue_set;
using collections::set;
using functors::comparable_functor;
using functors::compare_functors;
using lang::array;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using mlearning::clustering::clusterers::abstract::clusterer;

/*
 * Agglomerative cluterer on undirected graphs.
 *
 * V is the vertex type for the graph.
 * E is the edge type for the graph.
 * 
 * Note that V must be copy-constructible.
 */
template <typename V, typename E>
class basic_clusterer : public clusterer<V> {
public:
   /*
    * Abstract base class for agglomeration functionality.
    */
   class agglomerator {
   public:
      /*
       * Destructor.
       */
      virtual ~agglomerator() { }

      /*
       * Return whether the vertices at the end of the edge can be merged.
       * Agglomeration stops when the highest priority edge is unmergeable.
       */
      virtual bool is_mergeable(const E&) const = 0;
      
      /*
       * Merge two vertices and return the merged vertex.
       */
      virtual auto_ptr<V> merge(const V&, const V&) const = 0;

      /*
       * Compute the edge between the two vertices in the graph.
       */
      virtual auto_ptr<E> update(const V&, const V&) const = 0;
   };
   
   /*
    * Constructor.
    * Specify the agglomerator to use for merge and update operations.
    * Optionally specify the comparison functor for prioritizing edges.
    */
   explicit basic_clusterer(
      const agglomerator&,
      const comparable_functor<E>& = compare_functors<E>::f_compare()
   );

   /*
    * Copy constructor.
    * Create a clusterer that uses the same agglomerator and priority functor.
    */
   basic_clusterer(const basic_clusterer<V,E>&);

   /*
    * Destructor.
    */
   virtual ~basic_clusterer();

   /*
    * Clustering.
    * Assume a fully connected undirected graph.
    * Return the cluster assignments.
    */
   virtual array<unsigned long> cluster(
      const collection<V>&                         /* items to cluster */
   ) const;
  
   /*
    * Clustering.
    * Assume a fully connected undirected graph.
    * Return both the cluster assignments and the clustering result (the set
    * of unmerged vertices, one for each cluster).
    */
   virtual array<unsigned long> cluster(
      const collection<V>&,                        /* items to cluster */
      auto_collection< V, array_list<V> >&         /* clustering result */
   ) const;

   /*
    * Clustering.
    * Specify the edges in the undirected graph.
    * Return the cluster assignments.
    *
    * Vertices are numbered according to the order in which they appear in
    * the collection.  For each vertex, an array must be given that lists 
    * (by number) the vertices to which it is connected.
    *
    * Vertex numbering starts at zero.
    *
    * Only one direction of each undirected edge need be listed.
    */
   virtual array<unsigned long> cluster(
      const collection<V>&,                        /* items to cluster */
      const collection< array<unsigned long> >&    /* edges in graph */
   ) const;

   /*
    * Clustering.
    * Specify the edges in the undirected graph.
    * Return both the cluster assignments and the clustering result.
    */
   virtual array<unsigned long> cluster(
      const collection<V>&,                        /* items to cluster */
      const collection< array<unsigned long> >&,   /* edges in graph */
      auto_collection< V, array_list<V> >&         /* clustering result */
   ) const;

protected:
   /*
    * List of merged vertex numbers.
    */
   class v_list {
   public:
      /*
       * Constructor.
       * Create a one-element list containing the given vertex number.
       */
      explicit v_list(unsigned long num) : v_num(num), next() { }

      /*
       * Destructor.
       */
      ~v_list() { /* do nothing */ }

      /*
       * Vertex list data.
       */
      unsigned long    v_num;    /* vertex number */
      auto_ptr<v_list> next;     /* next list element */
   private:
      /*
       * Private copy constructor.
       * Vertex lists should not be copied.
       */
      v_list(const v_list& v_lst) : v_num(v_lst.v_num), next(v_lst.next) { }
   };
   
   /*
    * Declare node and edge classes.
    */
   class node;
   class edge;
   
   /*
    * Node in the graph.
    */
   class node {
   public:
      /*
       * Constructor.
       * Specify the vertex data.
       */
      explicit node(auto_ptr<V> v_data)
       : data(v_data),
         v_list_start(),
         v_list_end(NULL),
         edges_start(compare_functors<edge>::f_compare_memloc()),
         edges_end(compare_functors<edge>::f_compare_memloc())
      { }

      /*
       * Destructor.
       */
      ~node() { /* do nothing */ }
   
      /*
       * Node data.
       */
      auto_ptr<V>      data;           /* vertex data */
      auto_ptr<v_list> v_list_start;   /* list of merged original vertices */
      v_list*          v_list_end;     /* (start and end of list)          */
      set<edge>        edges_start;    /* edges that start at node */
      set<edge>        edges_end;      /* edges that end at node */
      
   private:
      /*
       * Private copy constructor.
       * Nodes should not be copied.
       */
      explicit node(const node& n)
       : data(n.data),
         v_list_start(n.v_list_start), 
         v_list_end(n.v_list_end),
         edges_start(n.edges_start),
         edges_end(n.edges_end)
      { }
   };
   
   /*
    * Edge in the graph.
    */
   class edge {
   public:
      /*
       * Constructor.
       * Specify vertex connectivity.
       */
      explicit edge(node* node_start, node* node_end)
       : data(), start(node_start), end(node_end) { }
      
      /*
       * Destructor.
       */
      ~edge() { /* do nothing */ }

      /*
       * Edge data.
       */
      auto_ptr<E> data;    /* edge data */
      node*       start;   /* node at which edge starts */
      node*       end;     /* node at which edge ends */

   private:
      /*
       * Private copy constructor.
       * Edges should not be copied.
       */
      explicit edge(const edge& e)
       : data(e.data), start(e.start), end(e.end) { }
   };

   /*
    * Priority comparison functor on edges.
    * Compare the edge data using the given comparison functor.
    */
   class edge_prio_functor : public comparable_functor<edge> {
   public:
      /*
       * Constructor.
       */
      explicit edge_prio_functor(const comparable_functor<E>& f)
       : _f(f) { }

      /*
       * Copy constructor.
       */
      explicit edge_prio_functor(const edge_prio_functor& f)
       : _f(f._f) { }

      /*
       * Priority comparison function.
       */
      int operator()(const edge& e0, const edge& e1) const {
         return _f(*(e0.data), *(e1.data));
      }

   protected:
      const comparable_functor<E>& _f;
   };
  
   /*
    * Search comparison functor on edges.
    * Compared undirected edges based on uniqueness of their vertices.
    */
   class edge_search_functor : public comparable_functor<edge> {
   public:
      int operator()(const edge& e0, const edge& e1) const {
         bool e0_cmp = (e0.start < e0.end);
         bool e1_cmp = (e1.start < e1.end);
         node* e0_min = (e0_cmp) ? (e0.start) : (e0.end);
         node* e1_min = (e1_cmp) ? (e1.start) : (e1.end);
         if (e0_min < e1_min) {
            return -1;
         } else if (e0_min > e1_min) {
            return 1;
         } else {
            node* e0_max = (e0_cmp) ? (e0.end) : (e0.start);
            node* e1_max = (e1_cmp) ? (e1.end) : (e1.start);
            return (e0_max < e1_max) ? (-1) : ((e0_max > e1_max) ? 1 : 0);
         }
      }
   };

   /*
    * Static search comparison functor on edges.
    */
   static const edge_search_functor& f_edge_search() {
      static const edge_search_functor* f = new edge_search_functor();
      return *f;
   }
   
   /*
    * Clusterer data.
    */
   const agglomerator&     _agglm;        /* agglomerator for merge, update */
   const edge_prio_functor _f_edge_prio;  /* edge priority functor */
};

/***************************************************************************
 * Clusterer implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Specify the agglomerator to use for merge and update operations.
 * Specify the comparison functor for prioritizing edges.
 */
template <typename V, typename E>
basic_clusterer<V,E>::basic_clusterer(
   const agglomerator&          agglm,
   const comparable_functor<E>& f)
 : _agglm(agglm),
   _f_edge_prio(f)
{ }

/*
 * Copy constructor.
 * Create a clusterer that uses the same agglomerator and priority functor.
 */
template <typename V, typename E>
basic_clusterer<V,E>::basic_clusterer(const basic_clusterer<V,E>& c)
 : _agglm(c._agglm),
   _f_edge_prio(c._f_edge_prio)
{ }

/*
 * Destructor.
 */
template <typename V, typename E>
basic_clusterer<V,E>::~basic_clusterer() {
   /* do nothing */
}

/*
 * Clustering.
 * Assume a fully connected undirected graph.
 * Return the cluster assignments.
 */
template <typename V, typename E>
array<unsigned long> basic_clusterer<V,E>::cluster(
   const collection<V>& items) const
{
   /* call the version that also returns vertices */
   auto_collection< V, array_list<V> > vertices;
   return this->cluster(items, vertices);
}

/*
 * Clustering.
 * Assume a fully connected undirected graph.
 * Return both the cluster assignments and the clustering result (the set
 * of unmerged vertices, one for each cluster).
 */
template <typename V, typename E>
array<unsigned long> basic_clusterer<V,E>::cluster(
   const collection<V>&                 items,
   auto_collection< V, array_list<V> >& vertices) const
{
   /* create edge indicators for fully connected undirected graph */
   auto_collection< array<unsigned long>, list< array<unsigned long> > > 
      edges(new list< array<unsigned long> >());
   unsigned long n_items = items.size();
   for (unsigned long n = 0; n < n_items; n++) {
      auto_ptr< array<unsigned long> > edge_inds(
         new array<unsigned long>(n_items - n - 1)
      );
      for (unsigned long nn = 0; nn < (n_items - n - 1); nn++)
         (*edge_inds)[nn] = n + nn + 1;
      edges->add(*edge_inds);
      edge_inds.release();
   }
   /* cluster */
   return this->cluster(items, *edges, vertices);
}

/*
 * Clustering.
 * Specify the edges in the undirected graph.
 * Return the cluster assignments.
 */
template <typename V, typename E>
array<unsigned long> basic_clusterer<V,E>::cluster(
   const collection<V>&                      items,
   const collection< array<unsigned long> >& edges) const
{
   /* call the version that also returns vertices */
   auto_collection< V, array_list<V> > vertices;
   return this->cluster(items, edges, vertices);
}

/*
 * Clustering.
 * Specify the edges in the undirected graph.
 * Return both the cluster assignments and the clustering result.
 */
template <typename V, typename E>
array<unsigned long> basic_clusterer<V,E>::cluster(
   const collection<V>&                      items,
   const collection< array<unsigned long> >& edges,
   auto_collection< V, array_list<V> >&      vertices) const
{
   /* allocate data structures to hold graph nodes and edges */
   auto_collection< node, set<node> > node_set(
      new set<node>(compare_functors<node>::f_compare_memloc())
   );
   auto_collection<edge, queue_set<edge> > edge_q(
      new queue_set<edge>(_f_edge_prio, basic_clusterer<V,E>::f_edge_search())
   );
   /* create nodes */
   array_list<node> node_list;
   unsigned long n_items = 0;
   for (auto_ptr< iterator<V> > i = items.iter_create(); i->has_next(); ) {
      V& item = i->next();
      auto_ptr<V>      item_copy(new V(item));
      auto_ptr<v_list> v_lst(new v_list(n_items));
      auto_ptr<node>   v(new node(item_copy));
      v->v_list_start = v_lst;
      v->v_list_end = v->v_list_start.get();
      node_list.add(*v);
      node_set->add(*v);
      v.release();
      n_items++;
   }
   /* create edges */
   {
      auto_ptr< iterator< array<unsigned long> > > i = edges.iter_create();
      for (unsigned long n = 0; n < n_items; n++) {
         array<unsigned long>& edge_inds = i->next();
         node& v0 = node_list[n];
         unsigned long n_edges = edge_inds.size();
         for (unsigned long n_e = 0; n_e < n_edges; n_e++) {
            node& v1 = node_list[edge_inds[n_e]];
            auto_ptr<edge> e(new edge(&v0, &v1));
            if (!(edge_q->contains(*e))) {
               e->data = _agglm.update(*(v0.data), *(v1.data));
               v0.edges_start.add(*e);
               v1.edges_end.add(*e);
               edge_q->add(*e);
               e.release();
            }
         }
      }
   }
   /* iteratively merge nodes */
   while (!(edge_q->is_empty())) {
      /* dequeue highest priority edge */
      auto_ptr<edge> e_best(&(edge_q->dequeue()));
      /* check that edge is valid merge candidate */
      if (!(_agglm.is_mergeable(*(e_best->data))))
         break;
      /* merge nodes - merge vertex data */
      node* v0 = e_best->start;
      node* v1 = e_best->end;
      auto_ptr<V> v_data = _agglm.merge(*(v0->data), *(v1->data));
      /* merge nodes - create merged node */
      auto_ptr<node> v(new node(v_data));
      /* merge nodes - merge vertex lists */
      v0->v_list_end->next = v1->v_list_start;
      v->v_list_start = v0->v_list_start;
      v->v_list_end   = v1->v_list_end;
      /* merge nodes - merge edge start sets */
      list<edge> v_edges_start;
      for (typename set<edge>::iterator_t i(v0->edges_start); i.has_next(); ) {
         edge& e = i.next();
         if (e.end != v1)
            v_edges_start.add(e);
      }
      v_edges_start.add(v1->edges_start);
      /* merge nodes - merge edge end sets */
      list<edge> v_edges_end;
      for (typename set<edge>::iterator_t i(v1->edges_end); i.has_next(); ) {
         edge& e = i.next();
         if (e.start != v0)
            v_edges_end.add(e);
      }
      v_edges_end.add(v0->edges_end);
      /* remove original nodes */
      node_set->remove(*v0);
      delete v0;
      node_set->remove(*v1);
      delete v1;
      /* update edges containing merged node */
      for (typename list<edge>::iterator_t i(v_edges_start); i.has_next(); ) {
         edge& e = i.next();
         edge_q->remove(e);
         auto_ptr<edge> e_ptr(&e);  
         e.start = v.get();
         if (edge_q->contains(e)) {
            e.end->edges_end.remove(e);
         } else {
            e.data = _agglm.update(*(e.start->data), *(e.end->data));
            v->edges_start.add(e);
            edge_q->add(e);
            e_ptr.release();
         }
      }
      for (typename list<edge>::iterator_t i(v_edges_end); i.has_next(); ) {
         edge& e = i.next();
         edge_q->remove(e);
         auto_ptr<edge> e_ptr(&e);
         e.end = v.get();
         if (edge_q->contains(e)) {
            e.start->edges_start.remove(e);
         } else {
            e.data = _agglm.update(*(e.start->data), *(e.end->data));
            v->edges_end.add(e);
            edge_q->add(e);
            e_ptr.release();
         }
      }
      /* add merged node to node set */
      node_set->add(*v);
      v.release();
   }
   /* number and return the final set of vertices */
   array<unsigned long> assignments(n_items);
   vertices.reset(new array_list<V>());
   unsigned long n = 0;
   for (typename set<node>::iterator_t i(*node_set); i.has_next(); n++) {
      node& v = i.next();
      vertices->add(*(v.data));
      v.data.release();
      for (v_list* v_lst = v.v_list_start.get();
           v_lst != NULL;
           v_lst = v_lst->next.get())
         assignments[v_lst->v_num] = n;
   }
   return assignments;         
}

} /* namespace graph */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
