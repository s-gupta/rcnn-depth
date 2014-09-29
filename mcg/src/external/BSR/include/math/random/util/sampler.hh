/*
 * Sampler (thread-safe).
 */
#ifndef MATH__RANDOM__UTIL__SAMPLER_HH
#define MATH__RANDOM__UTIL__SAMPLER_HH

#include "collections/abstract/collection.hh"
#include "collections/list.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/random/generators/rand_gen_uniform.hh"

namespace math {
namespace random {
namespace util {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::list;
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using lang::exceptions::ex_not_found;
using lang::pointers::auto_ptr;
using math::random::generators::rand_gen_uniform;

/*
 * Sampler.
 */
template <typename T, typename Syn = unsynchronized>
class sampler : protected Syn {
public:
   /*
    * Constructor.
    */
   sampler(
      const collection<T>&,      /* items to sample */
      const collection<double>&, /* weight of each item */
      bool = true                /* is sampling with replacement? */
   );

   /*
    * Copy constructor.
    */
   sampler(const sampler<T,Syn>&);
   
   /*
    * Destructor.
    */
   virtual ~sampler();
   
   /*
    * Sample an item.
    * Throw an ex_not_found exception if there are no more samples to draw.
    */
   T& sample();

   /*
    * Size.
    * The size is the number of distinct items that may be sampled.
    */
   bool is_empty() const;
   unsigned long size() const;

protected:
   /************************************************************************
    * Sampler data structures.
    ************************************************************************/
    
   /*
    * Item stored at a leaf node.
    */
   class tree_item {
   public:
      /*
       * Constructors.
       */
      explicit tree_item(T&);
      explicit tree_item(const tree_item&);

      /*
       * Destructor.
       */
      ~tree_item();
      
      /*
       * Reference to item.
       */
      T& t;
   };

   /*
    * Node in tree.
    */
   class tree_node {
   public:
      /*
       * Constructor.
       * Create an empty node.
       */
      tree_node();
      
      /*
       * Constructor.
       * Create a node containing the given item with the given weight.
       */
      explicit tree_node(auto_ptr<tree_item>, double);


      /*
       * Copy constructor.
       * Create a node with the same contents.
       * Initialize parent, left, and right to NULL.
       */
      explicit tree_node(const tree_node&);

      /*
       * Destructor.
       */
      ~tree_node();
      
      /*
       * Node data.
       */
      auto_ptr<tree_item> item;     /* item (if leaf node) */
      double              weight;   /* total weight of subtree */
      tree_node*          parent;   /* parent node */
      tree_node*          left;     /* left subtree  */
      tree_node*          right;    /* right subtree */
      unsigned long       size;     /* total # items in subtree */
   };

   /*
    * Sampler data.
    */
   tree_node*         _root;     /* root node of tree */
   bool               _replace;  /* is sampling with replacement? */
   rand_gen_uniform<> _rand;     /* random number generator for sampling */

   /************************************************************************
    * Sampler helper functions.
    ************************************************************************/

   /*
    * Swap the root nodes of two samplers.
    * Note: Samplers are not locked by this function.
    */
   static void swap_root(sampler<T,Syn>&, sampler<T,Syn>&);

   /*
    * Constructor returning empty sampler.
    */
   sampler();
   
   /*
    * Build a sampler over the given items with the given weights.
    */
   void build(const collection<T>&, const collection<double>&);

   /*
    * Runnable object for recursively constructing a sampling tree.
    *
    * Building the tree takes O(n + (n/p)*log(n/p)) time, where 
    * n is the tree size and p is the number of available processors.
    */
   class builder : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit builder(
         tree_node*&,               /* root of subtree */
         auto_ptr< list<T> >,       /* items to place in subtree */
         auto_ptr< list<double> >   /* weight of items */
      );

      /*
       * Destructor.
       */
      virtual ~builder();

      /*
       * Build the tree.
       */
      virtual void run();

   protected:
      tree_node*&              _root;
      auto_ptr< list<T> >      _items;
      auto_ptr< list<double> > _weights;
   };

   /*
    * Runnable object for recursively copying a sampling tree.
    * Copying the tree takes O(n/p + log(n)) time.
    */
   class copier : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit copier(
         tree_node*,       /* source */
         tree_node*&       /* destination */
      );

      /*
       * Destructor.
       */
      virtual ~copier();

      /*
       * Copy the tree.
       */
      virtual void run();

   protected:
      tree_node*  _src;
      tree_node*& _dest;
   };
   
   /*
    * Runnable object for recursively destructing a sampling tree.
    * Destroying the tree takes O(n/p + log(n)) time.
    */
   class destroyer : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit destroyer(tree_node*);

      /*
       * Destructor.
       */
      virtual ~destroyer();

      /*
       * Destroy the tree.
       */
      virtual void run();
      
   protected:
      tree_node* _root;
   };
};
 
/***************************************************************************
 * Sampler helper class implementation - tree_item.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_item::tree_item(T& item)
 : t(item)
{ }

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_item::tree_item(const tree_item& t_item)
 : t(t_item.t)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_item::~tree_item() { 
   /* do nothing */
}

/***************************************************************************
 * Sampler helper class implementation - tree_node.
 ***************************************************************************/

/*
 * Constructor.
 * Create an empty node.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_node::tree_node()
 : item(NULL),
   weight(0),
   parent(NULL),
   left(NULL),
   right(NULL),
   size(0)
{ }

/*
 * Constructor.
 * Create a node containing the given item with the given weight.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_node::tree_node(
   auto_ptr<tree_item> t_item,
   double              t_weight)
 : item(t_item),
   weight(t_weight),
   parent(NULL),
   left(NULL),
   right(NULL),
   size((item.get() != NULL) ? 1 : 0)
{ }

/*
 * Copy constructor.
 * Create a node with the same contents.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_node::tree_node(const tree_node& t_node)
 : item((t_node.item.get() == NULL) ? NULL : (new tree_item(*(t_node.item)))),
   weight(t_node.weight),
   parent(NULL),
   left(NULL),
   right(NULL),
   size((item.get() != NULL) ? 1 : 0)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::tree_node::~tree_node() {
   /* do nothing */
}

/***************************************************************************
 * Sampler implementation.
 ***************************************************************************/

/*
 * Swap the root nodes of two samplers.
 * Note: Samplers are not locked by this function.
 */
template <typename T, typename Syn>
void sampler<T,Syn>::swap_root(sampler<T,Syn>& s0, sampler<T,Syn>& s1) {
   tree_node* temp_root = s0._root;
   s0._root = s1._root;
   s1._root = temp_root;
}

/*
 * Constructor returning empty sampler.
 */
template <typename T, typename Syn>
sampler<T,Syn>::sampler()
 : Syn(),
   _root(NULL),
   _replace(true),
   _rand()
{ }

/*
 * Constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::sampler(
   const collection<T>&      items,
   const collection<double>& weights,
   bool                      replace)
 : Syn(),
   _root(NULL),
   _replace(replace),
   _rand()
{
   /* build sampling tree */
   this->build(items, weights);
}

/*
 * Build a sampler from the given collection.
 */
template <typename T, typename Syn>
void sampler<T,Syn>::build(
   const collection<T>&      items,
   const collection<double>& weights)
{
   /* create temporary sampler */
   sampler<T,Syn> s;
   /* build tree for sampling */
   auto_ptr< list<T> >      items_list(new list<T>(items));
   auto_ptr< list<double> > weights_list(new list<double>(weights));
   builder b(s._root, items_list, weights_list);
   b.run();
   /* take ownership of tree */
   sampler<T,Syn>::swap_root(*this, s);
}

/*
 * Builder - Constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::builder::builder(
   tree_node*&              root,
   auto_ptr< list<T> >      items,
   auto_ptr< list<double> > weights)
 : _root(root),
   _items(items),
   _weights(weights)
{ }

/*
 * Builder - Destructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::builder::~builder() {
   /* do nothing */
}

/*
 * Builder - Build the tree.
 */
template <typename T, typename Syn>
void sampler<T,Syn>::builder::run() {
   unsigned long n_items = _items->size();
   if (n_items == 0) {
      /* empty tree */
      _root = NULL;
   } else if (n_items == 1) {
      /* leaf node containing single element */
      auto_ptr<tree_item> item(new tree_item(_items->head()));
      _root = new tree_node(item, _weights->head());
   } else {
      /* break items and weights into two */
      auto_ptr< list<T> >      items_list(new list<T>());
      auto_ptr< list<double> > weights_list(new list<double>());
      for (unsigned long n = 0; n < (n_items/2); n++) {
         items_list->add(_items->remove_head());
         weights_list->add(_weights->remove_head());
      }
      /* build the root node */
      _root = new tree_node();
      /* recursively build the subtrees */
      builder b_left( _root->left, items_list, weights_list);
      builder b_right(_root->right, _items, _weights);
      child_thread::run(b_left, b_right);
      /* link subtrees to root */
      if (_root->left != NULL) {
         _root->left->parent = _root;
         _root->size   += _root->left->size;
         _root->weight += _root->left->weight;
      }
      if (_root->right != NULL) {
         _root->right->parent = _root;
         _root->size   += _root->right->size;
         _root->weight += _root->right->weight;
      }
   }
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::sampler(const sampler<T,Syn>& s)
 : Syn(),
   _root(NULL),
   _replace(true),
   _rand()
{
   /* lock source sampler */
   auto_read_lock<const Syn> rlock(s);
   /* set replacement option */
   _replace = s._replace;
   /* create temporary sampler */
   sampler<T,Syn> s_copy;
   /* copy tree */
   copier c(s._root, s_copy._root);
   c.run();
   /* take ownership of copy */
   sampler<T,Syn>::swap_root(*this, s_copy);
}

/*
 * Copier - Constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::copier::copier(tree_node* src, tree_node*& dest)
 : _src(src),
   _dest(dest)
{ }

/*
 * Copier - Destructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::copier::~copier() {
   /* do nothing */
}

/*
 * Copier - Copy the tree.
 */
template <typename T, typename Syn>
void sampler<T,Syn>::copier::run() {
   if (_src != NULL) {
      /* copy current node */
      _dest = new tree_node(*_src);
      /* recursively copy the subtrees */
      copier c_left(_src->left, _dest->left);
      copier c_right(_src->right, _dest->right);
      child_thread::run(c_left, c_right);
      /* link subtrees to parent */
      if (_dest->left != NULL) {
         _dest->left->parent = _dest;
         _dest->size += _dest->left->size;
      }
      if (_dest->right != NULL) {
         _dest->right->parent = _dest;
         _dest->size += _dest->right->size;
      }
   }
}
   
/*
 * Destructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::~sampler() {
   destroyer d(_root);
   d.run();
}

/*
 * Destroyer - Constructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::destroyer::destroyer(tree_node* root)
 : _root(root)
{ }

/*
 * Destroyer - Destructor.
 */
template <typename T, typename Syn>
sampler<T,Syn>::destroyer::~destroyer() {
   /* do nothing */
}

/*
 * Destroyer - Destroy the tree.
 */
template <typename T, typename Syn>
void sampler<T,Syn>::destroyer::run() {
   if (_root != NULL) {
      /* create destroyers for each subtree */
      destroyer d_left(_root->left);
      destroyer d_right(_root->right);
      /* delete the root */
      delete _root;
      /* recursively delete the subtrees */
      child_thread::run(d_left, d_right);
   }
}

/*
 * Sample an item.
 * Throw an ex_not_found exception if there are no more samples to draw.
 */
template <typename T, typename Syn>
T& sampler<T,Syn>::sample() {
   auto_write_lock<const Syn> wlock(*this);
   /* check that samples exist */
   if (_root == NULL)
      throw ex_not_found(
         "attempt to sample from empty sampler"
      );
   /* generate random number in range [0, total weight] */
   double r = _rand.generate() * (_root->weight);
   /* descend tree to leftmost leaf whos interval contains r */
   tree_node* curr = _root;
   while (curr->size > 1) {
      tree_node* left = curr->left;
      double w_left   = left->weight;
      if (r <= w_left) {
         /* branch left */
         curr = left;
      } else {
         /* branch right */
         curr = curr->right;
         r -= w_left;
      }
   }
   /* retreive leaf item */
   T& t = curr->item->t;
   /* check if sampling without replacement */
   if (!_replace) {
      /* remove leaf */
      tree_node* parent = curr->parent;
      if (parent == NULL) {
         /* current node was root */
         _root = NULL;
      } else {
         /* move sibling up the tree */
         tree_node* sib = 
            (curr == parent->left) ? 
               (parent->right) : 
               (parent->left);
         parent->item  = sib->item;
         parent->left  = sib->left;
         parent->right = sib->right;
         delete sib;
         if (parent->left != NULL)
            parent->left->parent = parent;
         if (parent->right != NULL)
            parent->right->parent = parent;
         /* update tree weights and sizes */
         double w = curr->weight;
         while (parent != NULL) {
            parent->weight -= w;
            parent->size--;
            parent = parent->parent;
         }
      }
      delete curr;
   }
   return t;
}

/*
 * Return whether there are any items left to be sampled.
 */
template <typename T, typename Syn>
bool sampler<T,Syn>::is_empty() const {
   return (this->size() == 0);
}

/*
 * Return the number of distinct items that may be sampled.
 */
template <typename T, typename Syn>
unsigned long sampler<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return ((_root == NULL) ? 0 : (_root->size));
}

} /* namespace util */
} /* namespace random */
} /* namespace math */

#endif
