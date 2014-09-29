/*
 * Splay sets (thread-safe).
 *
 * Splay sets are self-balancing binary search search trees which perform
 * basic insertion, search, and removal operations in O(log(n)) amortized time.
 * They support some additional search operations in addition to the standard
 * operations defined in the abstract::set class.
 *
 * Splay sets do not allow duplicate items.
 * Adding an item already in the set has no effect.
 */
#ifndef COLLECTIONS__SPLAY_SET_HH
#define COLLECTIONS__SPLAY_SET_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/set.hh"
#include "collections/list.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "functors/comparable_functors.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
/*
 * Imports.
 */
using collections::abstract::collection;
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using functors::comparable_functor;
using functors::compare_functors;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over splay sets.
 */
template <typename T, typename Syn>
class splay_set_iterator;

template <typename T, typename Syn>
class splay_set_iterator_reverse;

/*
 * Splay sets.
 */
template <typename T, typename Syn = unsynchronized>
class splay_set : public abstract::set<T>, 
                  protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class splay_set_iterator<T,Syn>;
   friend class splay_set_iterator_reverse<T,Syn>;
   
   /*
    * Define the iterator types.
    */
   typedef splay_set_iterator<T,Syn> iterator_t;
   typedef splay_set_iterator_reverse<T,Syn> iterator_reverse_t;

   /*
    * Constructor.
    */
   splay_set(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a splay set from a collection.
    */
   explicit splay_set(
      const collection<T>&, 
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );
   
   /*
    * Copy constructor.
    */
   splay_set(const splay_set<T,Syn>&);
      
   /*
    * Destructor.
    */
   virtual ~splay_set();

   /*
    * Add element(s) to the set.
    * Return a reference to the set.
    */
   splay_set<T,Syn>& add(T&);
   splay_set<T,Syn>& add(const collection<T>&);

   /*
    * Remove element(s) from the set.
    * Return a reference to the set.
    */
   splay_set<T,Syn>& remove(const T&);
   splay_set<T,Syn>& remove(const collection<T>&);

   /*
    * Remove all element(s) from the set.
    * Return a reference to the set.
    */
   splay_set<T,Syn>& clear();

   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not 
    * contained in the set.
    */
   bool contains(const T&) const;
   T& find(const T&) const;

   /*
    * Search.
    * Find the element in the set immediately previous/next of the given 
    * element (which itself may or may not be in the set).  Throw an exception
    * (ex_not_found) when attempting to retrieve a previous/next element which 
    * does not exist.
    */
   bool contains_prev(const T&) const;
   bool contains_next(const T&) const;
   T& find_prev(const T&) const;
   T& find_next(const T&) const;

   /*
    * Search.
    * Find the minium/maximum element in the set (as defined by the comparison
    * functor).  Throw an exception (ex_not_found) if the set is empty.
    */
   T& find_min() const;
   T& find_max() const;
  
   /*
    * Size.
    */
   unsigned long size() const;

   /*
    * Return iterators over elements.
    */
   auto_ptr< iterator<T> > iter_create() const;
   auto_ptr< iterator<T> > iter_reverse_create() const;

protected:
   /************************************************************************
    * Splay set data structures.
    ************************************************************************/

   /*
    * Node in the splay tree.
    */
   class splay_node {
   public:
      /*
       * Constructors.
       */
      explicit splay_node(T&);
      explicit splay_node(const splay_node&);

      /*
       * Destructor.
       */
      ~splay_node();

      /*
       * Node data.
       */
      T& t;                   /* item stored in the node */
      splay_node* parent;     /* parent node */
      splay_node* left;       /* left subtree */
      splay_node* right;      /* right subtree */
   };

   /*
    * Splay set data structure.
    */
   const comparable_functor<T>& _f_compare;  /* comparison function */
   unsigned long                _size;       /* number of elements in splay set */
   mutable splay_node*          _root;       /* root of splay tree */

   /************************************************************************
    * Splay set helper functions.
    ************************************************************************/
    
   /*
    * Swap the contents of two splay sets.
    * The splay sets must have the same comparison functor.
    * Note: Splay sets are not locked by this function.
    */
   static void swap(splay_set<T,Syn>&, splay_set<T,Syn>&);

   /*
    * Runnable object for recursively copying a splay set.
    * Copying the set takes O(n/p + log(n)) time where n is the
    * size of the set and p is the number of available processors.
    */
   class copier : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit copier(
         splay_node*,      /* source */
         splay_node*&      /* destination */
      );

      /*
       * Destructor.
       */
      virtual ~copier();

      /*
       * Copy.
       */
      virtual void run();

   protected:
      splay_node*  _src;
      splay_node*& _dest;
   };

   /*
    * Runnable object for recursively destructing a splay set.
    * Destroying the set takes O(n/p + log(n)) time.
    */
   class destroyer : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit destroyer(splay_node*);

      /*
       * Destructor.
       */
      virtual ~destroyer();

      /*
       * Destroy the set.
       */
      virtual void run();

   protected:
      splay_node* _root;
   };
   
   /*
    * Add an element to the set.
    */
   void add_item(T&);

   /*
    * Remove an element from the set.
    */
   void remove_item(const T&);
    
   /*
    * Rotate.
    *
    * The rotate operation promotes a node one level in the 
    * tree, preserving the inorder numbering of the tree. 
    * If the element is already the root, nothing is done.
    * Otherwise, the following rotation patterns are used:
    *
    *     y                              x
    *    / \     ---- rotate x --->     / \
    *   x   c                          a   y 
    *  / \       <--- rotate y ----       / \
    * a   b                              b   c
    */
   static void rotate(splay_node*);

   /*
    * Splay.
    *
    * The splay(t) operation rearranges the binary tree, preserving 
    * its inorder property, such that the new root of the tree is:
    *
    *    t (if t is in the tree) and otherwise is either
    *    min { k in tree | k > t } or 
    *    max { k in tree | k < t }
    *
    * This operation tends to yield a more balanced tree.
    */
   void splay(const T&) const;

   /*
    * Splay minimum.
    *
    * The splay_min() operation acts as if splay() were called on an 
    * element less than all items in the tree.  The new root is the 
    * minimum element in the tree.
    */
   void splay_min() const;

   /*
    * Splay maximum.
    *
    * The splay_max() operation acts as if splay() were called on an 
    * element greater than all items in the tree.  The new root is the 
    * maximum element in the tree.
    */
   void splay_max() const;
   
   /*
    * Recursive implementation of the splay operation.
    */
   static splay_node* splay(
      const T&,
      splay_node*, 
      const comparable_functor<T>&
   );

   /*
    * Recursive implementation of the splay_min operation.
    */
   static splay_node* splay_min(splay_node*);

   /*
    * Recursive implementation of the splay_max operation.
    */
   static splay_node* splay_max(splay_node*);
};

/*
 * Splay set iterator.
 */
template <typename T, typename Syn = unsynchronized>
class splay_set_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit splay_set_iterator(const splay_set<T,Syn>&);

   /*
    * Copy constructor.
    */
   splay_set_iterator(const splay_set_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~splay_set_iterator();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (ex_not_found) if there are no more items.
    */
   T& next();

protected:
   const splay_set<T,Syn>&                _set;    /* set being iterator over */
   typename splay_set<T,Syn>::splay_node* _curr;   /* current node in set */
};

/*
 * Splay set reverse iterator.
 */
template <typename T, typename Syn = unsynchronized>
class splay_set_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit splay_set_iterator_reverse(const splay_set<T,Syn>&);

   /*
    * Copy constructor.
    */
   splay_set_iterator_reverse(const splay_set_iterator_reverse<T,Syn>&);
   
   /*
    * Destructor.
    */
   virtual ~splay_set_iterator_reverse();

   /*
    * Check if there is another item available.
    */
   bool has_next() const;

   /*
    * Return the next item.
    * Throw an exception (ex_not_found) if there are no more items.
    */
   T& next();

protected:
   const splay_set<T,Syn>&                _set;    /* set being iterator over */
   typename splay_set<T,Syn>::splay_node* _curr;   /* current node in set */
};

/***************************************************************************
 * Splay set iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Lock splay set and initialize position.
 */
template <typename T, typename Syn>
splay_set_iterator<T,Syn>::splay_set_iterator(
   const splay_set<T,Syn>& s)
 : _set(s) 
{
   _set.read_lock();
   _curr = _set._root;
   if (_curr != NULL) {
      while (_curr->left != NULL)
         _curr = _curr->left;
   }
}

template <typename T, typename Syn>
splay_set_iterator_reverse<T,Syn>::splay_set_iterator_reverse(
   const splay_set<T,Syn>& s)
 : _set(s)
{
   _set.read_lock();
   _curr = _set._root;
   if (_curr != NULL) {
      while (_curr->right != NULL)
         _curr = _curr->right;
   }
}

/*
 * Copy constructors.
 */
template <typename T, typename Syn>
splay_set_iterator<T,Syn>::splay_set_iterator(
   const splay_set_iterator<T,Syn>& i)
 : _set(i._set),
   _curr(i._curr)
{
   _set.read_lock();
}

template <typename T, typename Syn>
splay_set_iterator_reverse<T,Syn>::splay_set_iterator_reverse(
   const splay_set_iterator_reverse<T,Syn>& i)
 : _set(i._set),
   _curr(i._curr)
{
   _set.read_lock();
}

/*
 * Destructors.
 * Unlock the splay set.
 */
template <typename T, typename Syn>
splay_set_iterator<T,Syn>::~splay_set_iterator() {
   _set.read_unlock();
}

template <typename T, typename Syn>
splay_set_iterator_reverse<T,Syn>::~splay_set_iterator_reverse() {
   _set.read_unlock();
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool splay_set_iterator<T,Syn>::has_next() const {
   return (_curr != NULL);
}

template <typename T, typename Syn>
bool splay_set_iterator_reverse<T,Syn>::has_next() const {
   return (_curr != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& splay_set_iterator<T,Syn>::next() {
   if (_curr != NULL) {
      /* get item at current node */
      typename splay_set<T,Syn>::splay_node* node = _curr;
      T& t = node->t;
      /* move to next node on inorder traversal */
      if (node->right != NULL) {
         /* right subtree unvisited -> move to leftmost node in right subtree */
         node = node->right;
         while (node->left != NULL)
            node = node->left;
      } else {
         /* both subtrees visited -> move up until we encounter parent from the left */
         typename splay_set<T,Syn>::splay_node* node_prev = node;
         node = node->parent;
         while (node != NULL) {
            if (node->left == node_prev)
               break;
            node_prev = node;
            node = node->parent;
         }
      }
      /* now at the next node */
      _curr = node;
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of splay_set"
      );
   }
}
 
template <typename T, typename Syn>
T& splay_set_iterator_reverse<T,Syn>::next() {
   if (_curr != NULL) {
      /* get item at current node */
      typename splay_set<T,Syn>::splay_node* node = _curr;
      T& t = node->t;
      /* move to next node on reverse inorder traversal */
      if (node->left != NULL) {
         /* left subtree unvisited -> move to rightmost node in left subtree */
         node = node->left;
         while (node->right != NULL)
            node = node->right;
      } else {
         /* both subtrees visited -> move up until we encounter parent from the right */
         typename splay_set<T,Syn>::splay_node* node_prev = node;
         node = node->parent;
         while (node != NULL) {
            if (node->right == node_prev)
               break;
            node_prev = node;
            node = node->parent;
         }
      }
      /* now at the next node */
      _curr = node;
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of splay_set"
      );
   }
}

/***************************************************************************
 * Splay node helper class implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Create a node containing the given item.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::splay_node::splay_node(T& item)
 : t(item),
   parent(NULL),
   left(NULL),
   right(NULL)
{ }

/*
 * Copy constructor.
 * Create a node identical to the given node.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::splay_node::splay_node(const splay_node& node)
 : t(node.t),
   parent(node.parent),
   left(node.left),
   right(node.right)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::splay_node::~splay_node() {
   /* do nothing */
}

/***************************************************************************
 * Splay set implementation.
 ***************************************************************************/

/*
 * Swap the contents of two splay sets.
 * The splay sets must have the same comparison functor.
 * Note: Splay sets are not locked by this function.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::swap(splay_set<T,Syn>& s0, splay_set<T,Syn>& s1) {
   unsigned long temp_size = s0._size;
   splay_node*   temp_root = s0._root;
   s0._size = s1._size;
   s0._root = s1._root;
   s1._size = temp_size;
   s1._root = temp_root;
}

/*
 * Constructor.
 * Specify the comparison function to use when ordering set elements.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::splay_set(const comparable_functor<T>& f)
 : Syn(), 
   _f_compare(f),
   _size(0),
   _root(NULL)
{ }

/*
 * Constructor.
 * Create a splay set from a collection.
 * Specify the comparison function to use when ordering set elements.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::splay_set(
   const collection<T>& c,
   const comparable_functor<T>& f)
 : Syn(), 
   _f_compare(f),
   _size(0),
   _root(NULL)
   
{
   /* add elements to set */
   splay_set<T,Syn> s(_f_compare);
   for (auto_ptr< iterator<T> > i = c.iter_create(); i->has_next(); )
      s.add_item(i->next());
   /* take ownership of set */
   splay_set<T,Syn>::swap(*this, s);
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::splay_set(const splay_set<T,Syn>& s)
 : Syn(), 
   _f_compare(s._f_compare)
{
   /* lock source set */
   auto_read_lock<const Syn> rlock(s);
   /* create copy */
   splay_set<T,Syn> s_copy(_f_compare);
   copier c(s._root, s_copy._root);
   c.run();
   /* take ownership of copy */
   splay_set<T,Syn>::swap(*this, s_copy);
}

/*
 * Copier - Constructor.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::copier::copier(splay_node* src, splay_node*& dest)
 : _src(src),
   _dest(dest)
{ }

/*
 * Copier - Destructor.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::copier::~copier() {
   /* do nothing */
}

/*
 * Copier - Copy.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::copier::run() {
   if (_src != NULL) {
      /* copy current node */
      _dest = new splay_node(_src->t);
      /* create copiers for each subtree */
      copier c_left(_src->left, _dest->left);
      copier c_right(_src->right, _dest->right);
      /* recursively copy the subtrees */
      child_thread::run(c_left, c_right);
      /* link subtrees to parent */
      if (_dest->left != NULL)
         _dest->left->parent = _dest;
      if (_dest->right != NULL)
         _dest->right->parent = _dest;
   }
}

/*
 * Destructor.
 * Delete the tree nodes.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::~splay_set() {
   destroyer d(_root);
   d.run();
}

/*
 * Destroyer - Constructor.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::destroyer::destroyer(splay_node* root) : _root(root) { }

/*
 * Destroyer - Destructor.
 */
template <typename T, typename Syn>
splay_set<T,Syn>::destroyer::~destroyer() {
   /* do nothing */
}

/*
 * Destroyer - Destroy the set.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::destroyer::run() {
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
 * Add an element to the set.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::add_item(T& t) {
   /* promote predecessor or successor to top */
   this->splay(t);
   splay_node* x = _root;
   /* check for case of empty tree */
   if (x == NULL) {
      /* empty tree - allocate root */
      _root = new splay_node(t);
      _size = 1;
   } else {
      /* compare root item to t */
      int x_i = _f_compare(x->t, t);
      if (x_i != 0) {
         /* t not in set - add it */
         splay_node* t_node = new splay_node(t);
         if (x_i < 0) {
            t_node->left  = x;
            t_node->right = x->right;
            x->right = NULL;
         } else {
            t_node->right = x;
            t_node->left  = x->left;
            x->left = NULL;
         }
         if (t_node->left != NULL)
            t_node->left->parent = t_node;
         if (t_node->right != NULL)
            t_node->right->parent = t_node;
         _root = t_node;
         _size++;
      }
   }
}

/*
 * Add an element to the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
splay_set<T,Syn>& splay_set<T,Syn>::add(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->add_item(t);
   return *this;
}

/*
 * Add multiple elements to the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
splay_set<T,Syn>& splay_set<T,Syn>::add(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->add_item(i.next());
   return *this;
}

/*
 * Remove an element from the set.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::remove_item(const T& t) {
   /* check that set is nonempty */
   if (_root != NULL) {
      /* promote item to top */
      this->splay(t);
      /* check that root is desired item */
      splay_node* x = _root;
      if (_f_compare(t, x->t) == 0) {
         /* remove item */
         splay_node* x_left  = x->left;
         splay_node* x_right = x->right;
         delete x;
         if (x_left != NULL) {
            /* meld left and right branches */
            x_left->parent = NULL;
            _root = x_left;
            this->splay(t);
            _root->right = x_right;
            if (x_right != NULL)
               x_right->parent = _root;
         } else {
            /* right branch is entire tree */
            _root = x_right;
            if (x_right != NULL)
               x_right->parent = NULL;
         }
         /* decrement size */
         _size--;
      }
   }
}

/*
 * Remove an element from the set.
 * Return a reference to the set. 
 */
template <typename T, typename Syn>
splay_set<T,Syn>& splay_set<T,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

/*
 * Remove multiple elements from the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
splay_set<T,Syn>& splay_set<T,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Remove all element(s) from the set.
 * Return a reference to the set.
 */
template <typename T, typename Syn>
splay_set<T,Syn>& splay_set<T,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   splay_set<T,Syn> s(_f_compare);
   splay_set<T,Syn>::swap(*this, s);
   return *this;
}

/*
 * Check if an element is in the set.
 */
template <typename T, typename Syn>
bool splay_set<T,Syn>::contains(const T& t) const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      return false;
   /* bring element to root */
   this->splay(t);
   return (_f_compare(t, _root->t) == 0);
}

/*
 * Find the element in the set that is equal to the given element.
 * Throw an exception (ex_not_found) if the element is not in the set.
 */
template <typename T, typename Syn>
T& splay_set<T,Syn>::find(const T& t) const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      throw ex_not_found("item not found in set");
   /* bring element to root */
   this->splay(t);
   T& t_found = _root->t;
   /* check equality */
   if (_f_compare(t, t_found) != 0)
      throw ex_not_found("item not found in set");
   return t_found;
}

/*
 * Check if the set contains an element that comes before the given element.
 */
template <typename T, typename Syn>
bool splay_set<T,Syn>::contains_prev(const T& t) const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      return false;
   /* bring element to root */
   this->splay(t);
   /* check for previous element */
   if (_f_compare(t, _root->t) > 0) {
      /* root element is previous */
      return true;
   } else {
      /* root element is min { k in tree | k >= t } */
      return (_root->left != NULL);
   }
}

/*
 * Check if the set contains an element that comes after the given element.
 */
template <typename T, typename Syn>
bool splay_set<T,Syn>::contains_next(const T& t) const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      return false;
   /* bring element to root */
   this->splay(t);
   /* check for next element */
   if (_f_compare(t, _root->t) < 0) {
      /* root element is next */
      return true;
   } else {
      /* root element is max { k in tree | k <= t } */
      return (_root->right != NULL);
   }
}

/*
 * Find the element in the set immediately previous of the given element (which
 * itself may or may not be in the set).  Throw an exception (ex_not_found)
 * if no previous element exists.
 */
template <typename T, typename Syn>
T& splay_set<T,Syn>::find_prev(const T& t) const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      throw ex_not_found("item not found in set");
   /* bring element to root */
   this->splay(t);
   /* check for previous element */
   if (_f_compare(t, _root->t) > 0) {
      /* root element is previous */
      return _root->t;
   } else {
      /* root element is min { k in tree | k >= t } */
      /* previous element is rightmost element in left subtree */
      splay_node* curr = _root->left;
      if (curr == NULL)
         throw ex_not_found("item not found in set");
      while (curr->right != NULL)
         curr = curr->right;
      return curr->t;
   }
}

/*
 * Find the element in the set immediately next of the given element (which
 * itself may or may not be in the set).  Throw an exception (ex_not_found)
 * if no next element exists.
 */
template <typename T, typename Syn>
T& splay_set<T,Syn>::find_next(const T& t) const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      throw ex_not_found("item not found in set");
   /* bring element to root */
   this->splay(t);
   /* check for next element */
   if (_f_compare(t, _root->t) < 0) {
      /* root element is next */
      return _root->t;
   } else {
      /* root element is max { k in tree | k <= t } */
      /* next element is leftmost element in right subtree */
      splay_node* curr = _root->right;
      if (curr == NULL)
         throw ex_not_found("item not found in set");
      while (curr->left != NULL)
         curr = curr->left;
      return curr->t;
   }
}

/*
 * Find the minium element in the set (as defined by the comparison functor).
 * Throw an exception (ex_not_found) if the set is empty.
 */
template <typename T, typename Syn>
T& splay_set<T,Syn>::find_min() const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      throw ex_not_found("item not found in set");
   /* bring minimum element to root */
   this->splay_min();
   return _root->t;
}

/*
 * Find the maximum element in the set (as defined by the comparison functor).
 * Throw an exception (ex_not_found) if the set is empty.
 */
template <typename T, typename Syn>
T& splay_set<T,Syn>::find_max() const {
   auto_write_lock<const Syn> wlock(*this);
   /* check that set is nonempty */
   if (_root == NULL)
      throw ex_not_found("item not found in set");
   /* bring maximum element to root */
   this->splay_max();
   return _root->t;
}

/*
 * Get number of elements in set.
 */
template <typename T, typename Syn>
unsigned long splay_set<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > splay_set<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(
      new splay_set_iterator<T,Syn>(*this)
   );
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > splay_set<T,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(
      new splay_set_iterator_reverse<T,Syn>(*this)
   );
}

/*
 * Rotate.
 *
 * The rotate operation promotes a node one level in the 
 * tree, preserving the inorder numbering of the tree. 
 * If the element is already the root, nothing is done.
 * Otherwise, the following rotation patterns are used:
 *
 *     y                              x
 *    / \     ---- rotate x --->     / \
 *   x   c                          a   y 
 *  / \       <--- rotate y ----       / \
 * a   b                              b   c
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::rotate(splay_node* x) {
   if (x != NULL) {
      splay_node* y = x->parent;
      if (y != NULL) {
         /* update y's parent to point at x */
         splay_node* y_parent = y->parent;
         if (y_parent != NULL) {
            if (y == y_parent->left)
               y_parent->left = x;
            else /* (y == y_parent->right) */
               y_parent->right = x;
         }
         x->parent = y_parent;
         /* rotate */
         splay_node* temp = NULL;
         if (x == y->left) {
            temp = x->right;
            x->right = y;
            y->left = temp;
         } else /* (x == y->right) */ {
            temp = x->left;
            x->left = y;
            y->right = temp;
         }
         /* update parent pointers */
         y->parent = x;
         if (temp != NULL)
            temp->parent = y;
      }
   }
}

/*
 * Splay.
 *
 * The splay(t) operation rearranges the binary tree, preserving 
 * its inorder property, such that the new root of the tree is:
 *
 *    t (if t is in the tree) and otherwise is either
 *    min { k in tree | k > t } or 
 *    max { k in tree | k < t }
 *
 * This operation tends to yield a more balanced tree.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::splay(const T& t) const {
   _root = splay_set<T,Syn>::splay(t, _root, _f_compare);
}

/*
 * Splay minimum.
 *
 * The splay_min() operation acts as if splay() were called on an 
 * element less than all items in the tree.  The new root is the 
 * minimum element in the tree.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::splay_min() const {
   _root = splay_set<T,Syn>::splay_min(_root);
}

/*
 * Splay maximum.
 *
 * The splay_max() operation acts as if splay() were called on an 
 * element greater than all items in the tree.  The new root is the 
 * maximum element in the tree.
 */
template <typename T, typename Syn>
void splay_set<T,Syn>::splay_max() const {
   _root = splay_set<T,Syn>::splay_max(_root);
}

/*
 * Recursive implementation of the splay operation.
 *
 * This function does the actual work for the splay
 * operation, promoting a given node to the root
 * of an ordered binary tree.
 *                                                   
 * If a node has no grandparent, it is rotated to the
 * root as follows:
 *
 *     y                               x
 *    / \     ---- promote x --->     / \
 *   x   c                           a   y 
 *  / \       <--- promote y ----       / \
 * a   b                               b   c
 *
 * Otherwise, if a node has a grandparent, it is rotated
 * up two levels in the tree according to the following
 * configurations:
 *
 *       z                           x
 *      / \   ---- promote x --->   / \
 *     y   d                       a   y
 *    / \                             / \
 *   x   c    <--- promote z ----    b   z
 *  / \                                 / \
 * a   b                               c   d
 *
 *     z                               x
 *    / \                             / \
 *   y   d    ---- promote x --->    /   \                  
 *  / \                             y     z
 * a   x                           / \   / \   
 *    / \                         a   b c   d     
 *   b   c                               
 *
 *     z                               x
 *    / \                             / \
 *   a   y    ---- promote x --->    /   \                  
 *      / \                         z     y
 *     x   d                       / \   / \   
 *    / \                         a   b c   d     
 *   b   c                               
 *
 * These promotion operations are applied until the node
 * is at the root through recursive calls to splay.
 *
 * splay takes 3 arguments:
 *
 * t          - the value to promote to the top of the tree
 *              (if t is not in the tree, either
 *               min { k in tree | k > t } or 
 *               max { k in tree | k < t }
 *              is promoted to become the root node)
 *
 * node       - the root node of the current subtree
 *
 * f_compare  - comparison functor for the set
 *
 * and returns the new root of the subtree.
 */
template <typename T, typename Syn>
typename splay_set<T,Syn>::splay_node* splay_set<T,Syn>::splay(
   const T& t, 
   splay_node* z, 
   const comparable_functor<T>& f_compare)
{
   /* check for case of empty tree */
   if (z == NULL)
      return NULL;
   /* explore branch one level below root */
   int z_i = f_compare(t, z->t);
   splay_node* y = 
      (z_i == 0) ? NULL :
      ((z_i < 0) ? 
         (splay_set<T,Syn>::splay(t, z->left, f_compare)) : 
         (splay_set<T,Syn>::splay(t, z->right, f_compare)));
   /* check if z is desired item */
   if (y == NULL)
      return z;
   /* explore branch two levels below root */
   int y_i = f_compare(t, y->t);
   splay_node* x = 
      (y_i == 0) ? NULL : 
      ((y_i < 0) ? (y->left) : (y->right));
   /* check if y is desired item */
   if (x == NULL) {
      if (z->parent != NULL) {
         /* recursive call from z's parent will promote y */
         return z;
      } else {
         /* promote y to root */
         splay_set<T,Syn>::rotate(y);
         return y;
      }
   }
   /* promote x */
   bool same_branch_direction = (z_i < 0) ? (y_i < 0) : (y_i > 0);
   if (same_branch_direction)
      splay_set<T,Syn>::rotate(y);
   else
      splay_set<T,Syn>::rotate(x);
   splay_set<T,Syn>::rotate(x);
   return x;
}

/*
 * Recursive implementation of the splay_min operation.
 */
template <typename T, typename Syn>
typename splay_set<T,Syn>::splay_node* splay_set<T,Syn>::splay_min(
   splay_node* z)
{
   /* check for case of empty tree */
   if (z == NULL)
      return NULL;
   /* explore branch one level below root */
   splay_node* y = splay_set<T,Syn>::splay_min(z->left);
   /* check if z is desired item */
   if (y == NULL)
      return z;
   /* explore branch two levels below root */
   splay_node* x = y->left;
   /* check if y is desired item */
   if (x == NULL) {
      if (z->parent != NULL) {
         /* recursive call from z's parent will promote y */
         return z;
      } else {
         /* promote y to root */
         splay_set<T,Syn>::rotate(y);
         return y;
      }
   }
   /* promote x */
   splay_set<T,Syn>::rotate(y);
   splay_set<T,Syn>::rotate(x);
   return x;
}

/*
 * Recursive implementation of the splay_max operation.
 */
template <typename T, typename Syn>
typename splay_set<T,Syn>::splay_node* splay_set<T,Syn>::splay_max(
   splay_node* z)
{
   /* check for case of empty tree */
   if (z == NULL)
      return NULL;
   /* explore branch one level below root */
   splay_node* y = splay_set<T,Syn>::splay_max(z->right);
   /* check if z is desired item */
   if (y == NULL)
      return z;
   /* explore branch two levels below root */
   splay_node* x = y->right;
   /* check if y is desired item */
   if (x == NULL) {
      if (z->parent != NULL) {
         /* recursive call from z's parent will promote y */
         return z;
      } else {
         /* promote y to root */
         splay_set<T,Syn>::rotate(y);
         return y;
      }
   }
   /* promote x */
   splay_set<T,Syn>::rotate(y);
   splay_set<T,Syn>::rotate(x);
   return x;
}

} /* namespace collections */

#endif
