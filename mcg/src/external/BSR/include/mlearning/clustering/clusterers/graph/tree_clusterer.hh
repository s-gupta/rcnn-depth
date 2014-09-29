/*
 * Clusterer on undirected graphs that produces an agglomeration tree.
 *
 * Graph vertices are initialized to each be a single element tree.  Each 
 * step of the agglomeration merges the vertices at the ends of the highest 
 * priority edge and updates the costs of the edges connected to the merged 
 * vertex.  The merged vertex is a tree containing the original vertices as
 * subtrees.  
 */
#ifndef MLEARNING__CLUSTERING__CLUSTERERS__GRAPH__TREE_CLUSTERER_HH
#define MLEARNING__CLUSTERING__CLUSTERERS__GRAPH__TREE_CLUSTERER_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "functors/comparable_functors.hh"
#include "lang/array.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "mlearning/clustering/clusterers/abstract/clusterer.hh"
#include "mlearning/clustering/clusterers/graph/basic_clusterer.hh"

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
using functors::comparable_functor;
using functors::compare_functors;
using lang::array;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using mlearning::clustering::clusterers::abstract::clusterer;

/*
 * Clusterer on undirected graphs that produces an agglomeration tree.
 * 
 * T is the type of data contained in each node of the tree.
 * E is the edge type for the graph.
 *
 * Note that T must be copy-constructible.
 *
 * Vertices in the graph are binary trees whos nodes are elements of type T.
 */
template <typename T, typename E>
class tree_clusterer : public clusterer<T> {
public:
   /*
    * Agglomeration tree.
    */
   class tree {
   public:
      /*
       * Constructors.
       */
      explicit tree(auto_ptr<T> t)
       : data(t), left(), right() { }
      
      explicit tree(auto_ptr<T> t, auto_ptr<tree> l, auto_ptr<tree> r) 
       : data(t), left(l), right(r) { }

      /*
       * Copy constructor.
       * Transfer ownership of tree's elements to the copy.
       */
      explicit tree(const tree& t) 
       : data(t.data), left(t.left), right(t.right) { }
       
      /*
       * Destructor.
       */
      ~tree() { /* do nothing */ }

      /*
       * Tree data.
       */
      mutable auto_ptr<T>    data;     /* data */
      mutable auto_ptr<tree> left;     /* left subtree */
      mutable auto_ptr<tree> right;    /* right subtree */
   };

   /*
    * Abstract base class for agglomeration functionality on trees.
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
       *
       * Default to always allowing merges.  In this case, the tree clusterer
       * returns a single tree rather than a forest of trees.
       */
      virtual bool is_mergeable(const E&) const {
         return true;
      }
      
      /*
       * Merge two vertices and return the data for the merged vertex.
       */
      virtual auto_ptr<T> merge(const tree&, const tree&) const = 0;

      /*
       * Compute the edge between the two vertices in the graph.
       */
      virtual auto_ptr<E> update(const tree&, const tree&) const = 0;
   };

   /*
    * Constructor. 
    * Specify the agglomerator to use for merge and update operations.
    * Optionally specify the comparison functor for prioritizing edges.
    */
   explicit tree_clusterer(
      const agglomerator&,
      const comparable_functor<E>& = compare_functors<E>::f_compare()
   );
   
   /*
    * Copy constructor.
    * Create a clusterer that uses the same agglomerator and priority functor.
    */
   tree_clusterer(const tree_clusterer<T,E>&);

   /*
    * Destructor.
    */
   virtual ~tree_clusterer();

   /*
    * Clustering.
    * Assume a fully connected undirected graph.
    * Return the cluster assignments.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&                         /* items to cluster */
   ) const;
   
   /*
    * Clustering.
    * Assume a fully connected undirected graph.
    * Return both the cluster assignments and the clustering result.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&,                        /* items to cluster */
      auto_collection< tree, array_list<tree> >&   /* clustering result */
   ) const;

   /*
    * Clustering.
    * Specify the edges in the undirected graph.
    * Return the cluster assignments.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&,                        /* items to cluster */
      const collection< array<unsigned long> >&    /* edges in graph */
   ) const;
   
   /*
    * Clustering.
    * Specify the edges in the undirected graph.
    * Return both the cluster assignments and the clustering result.
    */
   virtual array<unsigned long> cluster(
      const collection<T>&,                        /* items to cluster */
      const collection< array<unsigned long> >&,   /* edges in graph */
      auto_collection< tree, array_list<tree> >&   /* clustering result */
   ) const;

protected:
   /*
    * Agglomerator which translates a tree agglomerator into one 
    * that can be used with the basic graph clustering algorithm.
    */
   class basic_agglomerator
    : public basic_clusterer<tree,E>::agglomerator {
   public:
      /*
       * Constructor.
       */
      explicit basic_agglomerator(const agglomerator& agglm)
       : _agglm(agglm) { }

      /*
       * Copy constructor.
       */
      basic_agglomerator(const basic_agglomerator& agglm)
       : _agglm(agglm._agglm) { }
   
      /*
       * Return whether the vertices at the end of the edge can be merged.
       * Agglomeration stops when the highest priority edge is unmergeable.
       */
      bool is_mergeable(const E& e) const {
         return _agglm.is_mergeable(e);
      }
      
      /*
       * Merge two vertices and return the merged vertex.
       */
      auto_ptr<tree> merge(const tree& t0, const tree& t1) const {
         /* compute data for merged node */
         auto_ptr<T> data = _agglm.merge(t0, t1);
         /* construct tree that owns t0, t1 subtrees */
         auto_ptr<tree> t0_copy(new tree(t0));
         auto_ptr<tree> t1_copy(new tree(t1));
         auto_ptr<tree> t(new tree(data, t0_copy, t1_copy));
         return t;
      }

      /*
       * Compute the edge between the two vertices in the graph.
       */
      auto_ptr<E> update(const tree& t0, const tree& t1) const {
         return _agglm.update(t0, t1);
      }

   protected:
      const agglomerator& _agglm;   /* agglomerator on trees */
   };

   /*
    * Helper function.
    * Return a collection of single element trees given a collection of 
    * elements.
    */
   static auto_collection<tree> make_trees(const collection<T>&);

   /*
    * Tree clusterer data.
    */
   const basic_agglomerator      _agglm; /* agglomerator for basic clusterer */
   const comparable_functor<E>&  _f;     /* edge priority functor */
   const basic_clusterer<tree,E> _c;     /* basic graph clusterer */
};

/***************************************************************************
 * Tree clusterer implementation.
 ***************************************************************************/

/*
 * Constructor. 
 * Specify the agglomerator to use for merge and update operations.
 * Specify the comparison functor for prioritizing edges.
 */
template <typename T, typename E>
tree_clusterer<T,E>::tree_clusterer(
   const agglomerator&          agglm,
   const comparable_functor<E>& f)
 : _agglm(agglm),
   _f(f),
   _c(_agglm, _f)
{ }

/*
 * Copy constructor.
 * Create a clusterer that uses the same agglomerator and priority functor.
 */
template <typename T, typename E>
tree_clusterer<T,E>::tree_clusterer(const tree_clusterer<T,E>& c)
 : _agglm(c._agglm),
   _f(c._f),
   _c(_agglm, _f)
{ }

/*
 * Destructor.
 */
template <typename T, typename E>
tree_clusterer<T,E>::~tree_clusterer() {
   /* do nothing */
}

/*
 * Clustering.
 * Assume a fully connected undirected graph.
 * Return the cluster assignments.
 */
template <typename T, typename E>
array<unsigned long> tree_clusterer<T,E>::cluster(
   const collection<T>& items) const
{
   auto_collection<tree> trees = tree_clusterer<T,E>::make_trees(items);
   return _c.cluster(*trees);
}

/*
 * Clustering.
 * Assume a fully connected undirected graph.
 * Return both the cluster assignments and the clustering result.
 */
template <typename T, typename E>
array<unsigned long> tree_clusterer<T,E>::cluster(
   const collection<T>&                       items,
   auto_collection< tree, array_list<tree> >& vertices) const
{
   auto_collection<tree> trees = tree_clusterer<T,E>::make_trees(items);
   return _c.cluster(*trees, vertices);
}

/*
 * Clustering.
 * Specify the edges in the undirected graph.
 * Return the cluster assignments.
 */
template <typename T, typename E>
array<unsigned long> tree_clusterer<T,E>::cluster(
   const collection<T>&                      items,
   const collection< array<unsigned long> >& edges) const
{
   auto_collection<tree> trees = tree_clusterer<T,E>::make_trees(items);
   return _c.cluster(*trees, edges);
}

/*
 * Clustering.
 * Specify the edges in the undirected graph.
 * Return both the cluster assignments and the clustering result.
 */
template <typename T, typename E>
array<unsigned long> tree_clusterer<T,E>::cluster(
   const collection<T>&                       items,
   const collection< array<unsigned long> >&  edges,
   auto_collection< tree, array_list<tree> >& vertices) const
{
   auto_collection<tree> trees = tree_clusterer<T,E>::make_trees(items);
   return _c.cluster(*trees, edges, vertices);
}

/*
 * Helper function.
 * Return a collection of single element trees given a collection of 
 * elements.
 */
template <typename T, typename E>
auto_collection<typename tree_clusterer<T,E>::tree> 
   tree_clusterer<T,E>::make_trees(const collection<T>& items)
{
   auto_collection<tree> trees(new list<tree>);
   for (auto_ptr< iterator<T> > i = items.iter_create(); i->has_next(); ) {
      T& item = i->next();
      auto_ptr<T> item_copy(new T(item));
      auto_ptr<tree> t(new tree(item_copy));
      trees->add(*t);
      t.release();
   }
   return trees;
}

} /* namespace graph */
} /* namesapce clusterers */
} /* namespace clustering */
} /* namespace mlearning */

#endif
