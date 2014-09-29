/*
 * kd-trees (thread-safe).
 *
 * A kd-tree is a key tree whos internal keys consist of a split dimension and
 * a split value along that dimension and whos elements are members of some
 * D-dimensional space such that both the distance between elements and the 
 * distance from any element to an axis-aligned hyperplane can be measured.
 *
 * kd-trees enable efficient retrieval of approximate nearest neighbors for 
 * high dimensional data.
 */
#ifndef COLLECTIONS__KD_TREE_HH
#define COLLECTIONS__KD_TREE_HH

#include "collections/abstract/collection.hh"
#include "collections/key_tree.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "collections/queue.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
#include "functors/distanceable_functors.hh"
#include "functors/equalable_functors.hh"
#include "functors/filterable_functors.hh"
#include "functors/keyable_functors.hh"
#include "functors/kd_treeable_functors.hh"
#include "interfaces/comparable.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::pointers::auto_collection;
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using functors::compare_functors;
using functors::distanceable_functor;
using functors::distance_functors;
using functors::equalable_functor;
using functors::equal_functors;
using functors::filterable_functor;
using functors::filter_functors;
using functors::keyable_compare_functor;
using functors::keyable_split_functor;
using functors::key_compare_functors;
using functors::key_split_functors;
using functors::kd_treeable;
using functors::kd_tree_key;
using functors::kd_treeable_key_distance_functor;
using functors::kd_tree_key_distance_functors;
using interfaces::comparable;
using lang::exceptions::ex_not_found;
using lang::pointers::auto_ptr;

/*
 * Options for finding nearest neighbors in kd-trees.
 *
 * The search options returned by the constructor are those for finding exact
 * nearest neighbors.  However, note that finding exact nearest neighbors in
 * kd-trees for high-dimensional data can be an expensive operation.  The 
 * options below allow one to decrease the search time at the cost of 
 * returning only approximate, rather than exact, nearest neighbors.
 */
template <typename V = double>
class kd_tree_search_options {
public:
   /*
    * Constructor.
    * Return the options for finding the specified number 
    * of exact nearest neighbors of an item in a kd-tree.
    */
   explicit kd_tree_search_options(unsigned long = 1);

   /*
    * Copy constructor.
    */
   kd_tree_search_options(const kd_tree_search_options<V>&);
   
   /*
    * Destructor.
    */
   virtual ~kd_tree_search_options();

   /*
    * Set the number of nearest neighbors requested.
    */
   void set_num_nn(unsigned long);

   /*
    * Optionally limit the number of items examined during a search.
    */
   void set_item_limit(unsigned long);
   void unset_item_limit();
   
   /*
    * Optionally limit the maximum distance of any neighbor returned by
    * the search.
    */
   void set_distance_limit(V);
   void unset_distance_limit();

protected:
   unsigned long _num_nn;     /* number of nearest neighbors requested */
   bool _use_item_limit;      /* limit number of items search examines? */
   bool _use_dist_limit;      /* limit distance from search item? */
   unsigned long _item_limit; /* maximum number of items examined (if limited) */
   V _dist_limit;             /* maximum distance of neighbors from search item (if limited) */
};

/*
 * kd-tree class.
 */
template <typename T, typename V = double, typename Syn = unsynchronized>
class kd_tree : public key_tree<T, kd_tree_key<V>, Syn> {
public:
   /*
    * Constructor.
    * Return an empty kd-tree.
    * Optionally specify custom functors for tree operations.
    */
   explicit kd_tree(
      const equalable_functor<T>&                         = equal_functors<T>::f_equal(),
      const keyable_compare_functor< T, kd_tree_key<V> >& = (key_compare_functors< T, kd_tree_key<V> >::f_key_compare()),
      const keyable_split_functor< T, kd_tree_key<V> >&   = (key_split_functors<T, kd_tree_key<V> >::f_key_split()),
      const kd_treeable_key_distance_functor<T,V>&        = (kd_tree_key_distance_functors<T,V>::f_key_distance()),
      const distanceable_functor<T,V>&                    = (distance_functors<T,V>::f_distance()),
      unsigned long /* dimensionality of items */         = (kd_treeable<T,V>::dimensionality(list<T>()))
   );

   /*
    * Constructor.
    * Create a kd-tree from the given collection.
    * Optionally specify custom functors for tree operations.
    * Automatically determine the item dimensionality.
    */
   explicit kd_tree(
      const collection<T>&,
      const equalable_functor<T>&                         = equal_functors<T>::f_equal(),
      const keyable_compare_functor< T, kd_tree_key<V> >& = (key_compare_functors< T, kd_tree_key<V> >::f_key_compare()),
      const keyable_split_functor< T, kd_tree_key<V> >&   = (key_split_functors<T, kd_tree_key<V> >::f_key_split()),
      const kd_treeable_key_distance_functor<T,V>&        = (kd_tree_key_distance_functors<T,V>::f_key_distance()),
      const distanceable_functor<T,V>&                    = (distance_functors<T,V>::f_distance())
   );

   /*
    * Constructor.
    * Create a kd-tree from the given collection.
    * Specify custom functors for tree operations.
    * Also specify a custom item dimensionality.
    */
   explicit kd_tree(
      const collection<T>&,
      const equalable_functor<T>&,
      const keyable_compare_functor< T, kd_tree_key<V> >&,
      const keyable_split_functor< T, kd_tree_key<V> >&,
      const kd_treeable_key_distance_functor<T,V>&,
      const distanceable_functor<T,V>&,
      unsigned long /* dimensionality of items */
   );

   /*
    * Copy constructor.
    */
   kd_tree(const kd_tree<T,V,Syn>&);

   /*
    * Destructor.
    */
   virtual ~kd_tree();

   /*
    * Add element(s) after tree construction.
    * WARNING: This may cause the tree to become unbalanced.
    */
   kd_tree<T,V,Syn>& add(T&);
   kd_tree<T,V,Syn>& add(const collection<T>&);
   
   /*
    * Remove element(s) after tree construction.
    * WARNING: This may cause the tree to become unbalanced.
    */
   kd_tree<T,V,Syn>& remove(const T&);
   kd_tree<T,V,Syn>& remove(const collection<T>&);
   
   /*
    * Remove all instances of the element(s).
    * WARNING: This may cause the tree to become unbalanced.
    */
   kd_tree<T,V,Syn>& remove_all(const T&);
   kd_tree<T,V,Syn>& remove_all(const collection<T>&);
   
   /*
    * Remove all element(s) from the tree.
    * Return a reference to the tree.
    */
   kd_tree<T,V,Syn>& clear();

   /*
    * Find (exact/approximate) nearest neighbor.
    * Search for the given item in the tree and return the closest neighbor
    * found.  Throw an exception (ex_not_found) if no items satisfy the search
    * constraints.
    */
   T& find_nn(
      const T&, 
      const kd_tree_search_options<V>&   = kd_tree_search_options<V>(),
      const filterable_functor<const T>& = filter_functors<const T>::f_true()
   ) const;

   /*
    * Find multiple (exact/approximate) nearest neighbors.
    * Search for the given item in the tree and return the closest neighbor(s)
    * found.  An empty list is returned if no items satisfy the search 
    * constraints.
    */
   list<T> find_nns(
      const T&, 
      const kd_tree_search_options<V>&   = kd_tree_search_options<V>(),
      const filterable_functor<const T>& = filter_functors<const T>::f_true()
   ) const;
   
protected:
   /*
    * kd-tree data.
    */
   const kd_treeable_key_distance_functor<T,V>& _f_key_distance;  /* key distance function */
   const distanceable_functor<T,V>&             _f_distance;      /* item distance function */
   const unsigned long                          _n_dims;          /* dimensionality of items */

   /*
    * Declare tree branch and leaf classes.
    */
   class tree_branch;
   class tree_leaf;
   
   /*
    * Search status.
    */
   class search_status : public kd_tree_search_options<V> {
   public:
      /*
       * Constructor.
       * Initialize search status for the given search options.
       */
      explicit search_status(
         const kd_tree_search_options<V>&,
         const filterable_functor<const T>& 
      );
     
      /*
       * Copy constructor.
       */
      search_status(const search_status&);
      
      /*
       * Destructor.
       */
      virtual ~search_status();

      /*
       * Enqueue an internal branch during a search (if allowed by search
       * options).  Update the search status.
       */
      void enqueue_branch_q(queue<tree_branch>&, auto_ptr<tree_branch>);

      /*
       * Enqueue a leaf during a search (if allowed by search options).
       * Update the search status.
       */
      void enqueue_leaf_q(queue<tree_leaf>&, auto_ptr<tree_leaf>);
      
      /*
       * Check whether search has completed.
       */
      bool is_complete() const;

   protected:
      const filterable_functor<const T>& _f_filter;   /* filter for items found during search */
      V _dist_sq_limit;                               /* square of distance limit */
      bool _complete;                                 /* has search completed? */
   };
     
   /*
    * Branch of the tree (as viewed during a search).
    *
    * At a branch during the search we store:
    * (1) its corresponding root node in the tree
    * (2) coordinates of the closest point in the branch to the search item
    * (3) distance from the closest point in the branch to the search item
    *     where "closest point" means the closest point in the space covered
    *     by the branch (not the closest item the branch actually contains).
    */
   class tree_branch : public comparable<const tree_branch> {
   public:
      /*
       * Friend classes.
       */
      friend class search_status;

      /*
       * Constructor.
       * Create a branch representing the tree root.
       */
      explicit tree_branch(const kd_tree<T,V,Syn>&);

      /*
       * Destructor.
       */
      virtual ~tree_branch();

      /*
       * Descend a branch searching for the given item.
       *
       * All internal branches touched, but not taken, during the descent 
       * are enqueued into the first queue, while all leaf branches seen 
       * are enqueued into the second.
       *
       * The original branch ends up pointing to a leaf.
       */
      void descend(
         const T&, 
         queue<tree_branch>&, 
         queue<tree_leaf>&,
         const keyable_compare_functor<T, kd_tree_key<V> >&,
         const kd_treeable_key_distance_functor<T,V>&,
         const distanceable_functor<T,V>&,
         search_status&
      );
      
      /*
       * Comparison (for best-bin first priority search).
       */
      int compare_to(const tree_branch&) const;
     
   protected:
      /*
       * Copy constructor.
       */
      explicit tree_branch(const tree_branch&);
      
      /*
       * Branch data.
       */
      typename kd_tree<T,V,Syn>::tree_node* _node; /* node of the tree */
      lang::array<V> _dists;                       /* distance along each dimension to closest point in branch */
      V _dist_sq;                                  /* squared distance of closest point in branch */
   };

   /*
    * Leaf of the tree (as viewed during a search).
    */
   class tree_leaf : public comparable<const tree_leaf> {
   public:
      /*
       * Friend classes.
       */
      friend class search_status;
      friend class kd_tree<T,V,Syn>;
      
      /*
       * Constructor.
       */
      tree_leaf(T&, V);
      
      /*
       * Destructor.
       */
      virtual ~tree_leaf();
      
      /*
       * Comparison (for best-bin first priority search).
       */
      int compare_to(const tree_leaf&) const;
      
   protected:
      /*
       * Leaf data.
       */
      T& _t;         /* item contained in leaf */
      V  _dist_sq;   /* squared distance of leaf from search item */
   };
};

/***************************************************************************
 * Kd-tree search options implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Return the options for finding the specified number 
 * of exact nearest neighbors of an item in a kd-tree.
 */
template <typename V>
kd_tree_search_options<V>::kd_tree_search_options(unsigned long n)
 : _num_nn(n),
   _use_item_limit(false),
   _use_dist_limit(false),
   _item_limit(0),
   _dist_limit(V()) 
{ }

/*
 * Copy constructor.
 */
template <typename V>
kd_tree_search_options<V>::kd_tree_search_options(
   const kd_tree_search_options<V>& s_opt)
 : _num_nn(s_opt._num_nn), 
   _use_item_limit(s_opt._use_item_limit),
   _use_dist_limit(s_opt._use_dist_limit), 
   _item_limit(s_opt._item_limit), 
   _dist_limit(s_opt._dist_limit) 
{ }

/*
 * Destructor.
 */
template <typename V>
kd_tree_search_options<V>::~kd_tree_search_options() {
   /* do nothing */
}

/*
 * Set the number of nearest neighbors requested.
 */
template <typename V>
void kd_tree_search_options<V>::set_num_nn(unsigned long n) {
   _num_nn = n;
}

/*
 * Limit the number of items examined during a search.
 */
template <typename V>
void kd_tree_search_options<V>::set_item_limit(unsigned long n) {
   _use_item_limit = true;
   _item_limit = n;
}

/*
 * Remove limit on the number of items during a search.
 */
template <typename V>
void kd_tree_search_options<V>::unset_item_limit() {
   _use_item_limit = false;
}

/*
 * Limit the maximum distance of any neighbor returned by the search.
 */
template <typename V>
void kd_tree_search_options<V>::set_distance_limit(V dist) {
   _use_dist_limit = true;
   _dist_limit = dist;
}

/*
 * Remove limit on the maximum distance of any neighbor returned by the search.
 */
template <typename V>
void kd_tree_search_options<V>::unset_distance_limit() {
   _use_dist_limit = false;
}

/***************************************************************************
 * Kd-tree search status helper class implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Initialize search status for the given search options.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::search_status::search_status(
   const kd_tree_search_options<V>& s_opt,
   const filterable_functor<const T>& f_filter)
 : kd_tree_search_options<V>(s_opt),
   _f_filter(f_filter),
   _dist_sq_limit((this->_dist_limit)*(this->_dist_limit)),
   _complete(
      (this->_num_nn == 0) || 
      ((this->_use_item_limit) && (this->_item_limit == 0))
   )
{ }

/*
 * Copy constructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::search_status::search_status(const search_status& s)
 : kd_tree_search_options<V>(s),
   _f_filter(s._f_filter),
   _dist_sq_limit(s._dist_sq_limit),
   _complete(s._complete)
{ }
     
/*
 * Destructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::search_status::~search_status() {
   /* do nothing */
}

/*
 * Enqueue an internal branch during a search (if allowed by search options).
 * Update the search status.
 */
template <typename T, typename V, typename Syn>
void kd_tree<T,V,Syn>::search_status::enqueue_branch_q(
   queue<tree_branch>& branch_q, 
   auto_ptr<tree_branch> b) 
{
   /* check that branch does not exceed maximum allowable distance */
   if ((!(this->_use_dist_limit)) || (b->_dist_sq < this->_dist_sq_limit)) {
      /* branch is candidate for future search */
      branch_q.enqueue(*b);
      b.release();
   }
}

/*
 * Enqueue a leaf branch during a search (if allowed by search options).
 * Update the search status.
 */
template <typename T, typename V, typename Syn>
void kd_tree<T,V,Syn>::search_status::enqueue_leaf_q(
   queue<tree_leaf>& leaf_q, 
   auto_ptr<tree_leaf> l)
{
   /* check that search is not complete */
   if (!(_complete)) {
      /* check that leaf should not be ignored */
      if (_f_filter(l->_t)) {
         /* decrement item limit */
         if (this->_use_item_limit) {
            this->_item_limit--;
            _complete = (this->_item_limit == 0);
         }
         /* check distance limit */
         if ((!(this->_use_dist_limit)) || (l->_dist_sq < _dist_sq_limit)) {
            /* enqueue the current leaf */
            leaf_q.enqueue(*l);
            l.release();
            /* don't keep more leaves than necessary */
            if (leaf_q.size() > this->_num_nn) {
               /* delete the worst leaf */
               delete &(leaf_q.dequeue());
            }
            /* place more restrictive distance limit if possible */
            if (leaf_q.size() == this->_num_nn) {
               /* we are now only interested in neighbors */
               /* closer than those we have already found */
               this->_use_dist_limit = true;
               _dist_sq_limit = leaf_q.head()._dist_sq;
            }
         }
      }
   }
}

/*
 * Check whether search has completed.
 */
template <typename T, typename V, typename Syn>
bool kd_tree<T,V,Syn>::search_status::is_complete() const {
   return _complete;
}

/***************************************************************************
 * Kd-tree branch helper class implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Create a branch representing the tree root.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::tree_branch::tree_branch(const kd_tree<T,V,Syn>& k_tree)
 : _node(k_tree._root), 
   _dists(k_tree._n_dims), 
   _dist_sq(V()) 
{ }

/*
 * Copy constructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::tree_branch::tree_branch(const tree_branch& branch)
 : _node(branch._node), 
   _dists(branch._dists), 
   _dist_sq(branch._dist_sq) 
{ }

/*
 * Destructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::tree_branch::~tree_branch() {
   /* do nothing */
}

/*
 * Descend a branch searching for the given item.
 *
 * All internal branches touched, but not taken, during the descent 
 * are enqueued into the first queue, while all leaf branches seen 
 * are enqueued into the second.
 *      
 * The original branch ends up pointing at a leaf.
 */
template <typename T, typename V, typename Syn>
void kd_tree<T,V,Syn>::tree_branch::descend(
   const T&                                           t, 
   queue<tree_branch>&                                branch_q, 
   queue<tree_leaf>&                                  leaf_q,
   const keyable_compare_functor<T, kd_tree_key<V> >& f_key_compare,
   const kd_treeable_key_distance_functor<T,V>&       f_key_distance,
   const distanceable_functor<T,V>&                   f_distance,
   search_status&                                     s_status) 
{
   auto_ptr<tree_branch> this_branch(this);
   if (_node != NULL) {
      /* descend the current branch */
      while (_node->key.get() != NULL) {
         /* determine branch farther from item */
         kd_tree_key<V>& k = _node->key->k;
         auto_ptr<tree_branch> b(new tree_branch(*this));
         if (f_key_compare(t, k) < 0) {
            /* right branch is farther */
            b->_node = _node->right;
            _node    = _node->left;
         } else {
            /* left branch is farther */
            b->_node = _node->left;
            _node    = _node->right;
         }
         /* enqueue branch not taken */
         if (b->_node != NULL) {
            if (b->_node->key.get() != NULL) {
               /* internal branch - update distance to it */
               unsigned long dim = k.split_dimension();
               V dist           = f_key_distance(t, k);
               V dist_old       = b->_dists[dim];
               b->_dists[dim]   = dist;
               b->_dist_sq      = (dist*dist) - (dist_old*dist_old);
               /* enqueue into branch queue */
               s_status.enqueue_branch_q(branch_q, b);
            } else {
               /* leaf branch - enqueue all items into leaf queue */
               for (typename list<T>::iterator_t i(b->_node->items);
                    i.has_next(); )
               {
                  T& t_leaf = i.next();
                  V dist = f_distance(t, t_leaf);
                  auto_ptr<tree_leaf> l(new tree_leaf(t_leaf, dist*dist));
                  s_status.enqueue_leaf_q(leaf_q, l);
               }
            }
         }
      }
      /* current branch is now a leaf - enqueue all items in it */
      for (typename list<T>::iterator_t i(_node->items); i.has_next(); ) {
         T& t_leaf = i.next();
         V dist = f_distance(t, t_leaf);
         auto_ptr<tree_leaf> l(new tree_leaf(t_leaf, dist*dist));
         s_status.enqueue_leaf_q(leaf_q, l);
      }
   }
}

/*
 * Comparison (for best-bin first priority search).
 * The tree branch closer to the search item is given priority.
 */
template <typename T, typename V, typename Syn>
int kd_tree<T,V,Syn>::tree_branch::compare_to(const tree_branch& b) const {
   return (_dist_sq < b._dist_sq) ? (-1) : ((_dist_sq > b._dist_sq) ? 1 : 0);
}

/***************************************************************************
 * Kd-tree leaf helper class implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::tree_leaf::tree_leaf(T& t, V v)
 : _t(t),
   _dist_sq(v)
{ }

/*
 * Destructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::tree_leaf::~tree_leaf() {
   /* do nothing */
}

/*
 * Comparison (for best-bin first priority search).
 */
template <typename T, typename V, typename Syn>
int kd_tree<T,V,Syn>::tree_leaf::compare_to(const tree_leaf& l) const {
   return (_dist_sq < l._dist_sq) ? (-1) : ((_dist_sq > l._dist_sq) ? 1 : 0);
}

/***************************************************************************
 * Kd-tree implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Return an empty kd-tree.
 * Use the specified functors for tree operations.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::kd_tree(
   const equalable_functor<T>&                         f_eq,
   const keyable_compare_functor< T, kd_tree_key<V> >& f_key_c,
   const keyable_split_functor< T, kd_tree_key<V> >&   f_key_s,
   const kd_treeable_key_distance_functor<T,V>&        f_key_dist,
   const distanceable_functor<T,V>&                    f_dist,
   unsigned long                                       n_dims)
 : key_tree<T,kd_tree_key<V>,Syn>(
      f_eq, 
      f_key_c, 
      f_key_s), 
   _f_key_distance(f_key_dist), 
   _f_distance(f_dist), 
   _n_dims(n_dims) 
{ }

/*
 * Constructor.
 * Create a kd-tree from the given collection.
 * Use the specified functors for tree operations.
 * Automatically determine the item dimensionality.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::kd_tree(
   const collection<T>&                                c,
   const equalable_functor<T>&                         f_eq,
   const keyable_compare_functor< T, kd_tree_key<V> >& f_key_c,
   const keyable_split_functor< T, kd_tree_key<V> >&   f_key_s,
   const kd_treeable_key_distance_functor<T,V>&        f_key_dist,
   const distanceable_functor<T,V>&                    f_dist)
 : key_tree<T,kd_tree_key<V>,Syn>(
      c, 
      f_eq, 
      f_key_c, 
      f_key_s), 
   _f_key_distance(f_key_dist), 
   _f_distance(f_dist), 
   _n_dims(kd_treeable<T,V>::dimensionality(c)) 
{ }

/*
 * Constructor.
 * Create a kd-tree from the given collection.
 * Use the specified functors for tree operations.
 * Use the specified item dimensionality.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::kd_tree(
   const collection<T>&                                c,
   const equalable_functor<T>&                         f_eq,
   const keyable_compare_functor< T, kd_tree_key<V> >& f_key_c,
   const keyable_split_functor< T, kd_tree_key<V> >&   f_key_s,
   const kd_treeable_key_distance_functor<T,V>&        f_key_dist,
   const distanceable_functor<T,V>&                    f_dist,
   unsigned long                                       n_dims)
 : key_tree<T,kd_tree_key<V>,Syn>(
      c, 
      f_eq, 
      f_key_c, 
      f_key_s), 
   _f_key_distance(f_key_dist), 
   _f_distance(f_dist), 
   _n_dims(n_dims) 
{ }

/*
 * Copy constructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::kd_tree(const kd_tree<T,V,Syn>& k_tree)
 : key_tree<T,kd_tree_key<V>,Syn>(
      k_tree, 
      k_tree._f_equal, 
      k_tree._f_key_compare, 
      k_tree._f_key_split), 
   _f_key_distance(k_tree._f_key_distance), 
   _f_distance(k_tree._f_distance), 
   _n_dims(k_tree._n_dims)
{ }
   
/*
 * Destructor.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>::~kd_tree() {
   /* do nothing - superclass destructor should handle everything */
}

/*
 * Add element(s) after tree construction.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::add(T& t) {
   this->key_tree<T,kd_tree_key<V>,Syn>::add(t);
   return *this;
}

template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::add(const collection<T>& c) {
   this->key_tree<T,kd_tree_key<V>,Syn>::add(c);
   return *this;
}

/*
 * Remove element(s) after tree construction.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::remove(const T& t) {
   this->key_tree<T,kd_tree_key<V>,Syn>::remove(t);
   return *this;
}

template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::remove(const collection<T>& c) {
   this->key_tree<T,kd_tree_key<V>,Syn>::remove(c);
   return *this;
}

/*
 * Remove all instances of the element(s).
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::remove_all(const T& t) {
   this->key_tree<T,kd_tree_key<V>,Syn>::remove_all(t);
   return *this;
}

template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::remove_all(const collection<T>& c) {
   this->key_tree<T,kd_tree_key<V>,Syn>::remove_all(c);
   return *this;
}

/*
 * Remove all element(s) from the tree.
 * Return a reference to the tree.
 */
template <typename T, typename V, typename Syn>
kd_tree<T,V,Syn>& kd_tree<T,V,Syn>::clear() {
   this->key_tree<T,kd_tree_key<V>,Syn>::clear();
   return *this;
}

/*
 * Find (exact/approximate) nearest neighbor.
 * Search for the given item in the tree and return the closest neighbor
 * found.  Throw an exception (ex_not_found) if no items satisfy the search
 * constraints.
 */
template <typename T, typename V, typename Syn>
T& kd_tree<T,V,Syn>::find_nn(
   const T& t, 
   const kd_tree_search_options<V>& s_opt,
   const filterable_functor<const T>& f_filter) const 
{
   list<T> nns = this->find_nns(t, s_opt, f_filter);
   if (nns.is_empty())
      throw ex_not_found(
         "item not found in kd_tree"
      );
   else
      return nns.head();
}

/*
 * Find multiple (exact/approximate) nearest neighbors.
 * Search for the given item in the tree and return the closest neighbor(s)
 * found.  An empty list is returned if no items satisfy the search 
 * constraints.
 */
template <typename T, typename V, typename Syn>
list<T> kd_tree<T,V,Syn>::find_nns(
   const T& t, 
   const kd_tree_search_options<V>& s_opt,
   const filterable_functor<const T>& f_filter) const 
{
   /* lock tree */
   auto_read_lock<const Syn> rlock(*this);
   /* initialize search status */
   search_status s_status(s_opt, f_filter);
   /* initialize search queues */
   auto_collection< tree_leaf, queue<tree_leaf> > leaf_q(
      new queue<tree_leaf>(
         compare_functors<tree_leaf>::f_compare_reverse()
      )
   );
   auto_collection< tree_branch, queue<tree_branch> > branch_q(
      new queue<tree_branch>()
   );
   /* enqueue the root */
   auto_ptr<tree_branch> root_branch(new tree_branch(*this));
   branch_q->enqueue(*root_branch);
   root_branch.release();
   /* search */
   do {
      /* explore the closest branch */
      branch_q->dequeue().descend(
         t,
         *branch_q,
         *leaf_q,
         this->_f_key_compare,
         _f_key_distance,
         _f_distance,
         s_status
      );
   } while (!(s_status.is_complete() || branch_q->is_empty()));
   /* place all neighbors found into list ordered by distance */
   list<T> nns;
   while (!(leaf_q->is_empty())) {
      auto_ptr<tree_leaf> l(&(leaf_q->dequeue()));
      nns.prepend(l->_t);
   }
   /* return neighbors */
   return nns;
}

} /* namespace collections */

#endif
