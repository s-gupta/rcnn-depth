/*
 * Key trees (thread-safe).
 *
 * A key tree is a binary tree whos internal nodes store keys and whos leaf 
 * nodes store elements of a collection.  Each key separates (by means of 
 * comparison) the elements in its two subtrees.  Leaves contain collections
 * of element(s) which cannot be split further using a key.
 *
 * Searching for an element is equivalent to descending the tree by comparing
 * the element to successive keys and finally checking equality between the 
 * element and those items stored at the resulting leaf node.
 *
 * Key trees are multisets.  They allow duplicate elements.
 *
 * Key trees can be automatically created for object which implement the 
 * equalable and keyable interfaces.  Key trees for other classes may be 
 * constructed by specifying appropriate functors at creation time.
 */
#ifndef COLLECTIONS__KEY_TREE_HH
#define COLLECTIONS__KEY_TREE_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/multiset.hh"
#include "collections/list.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "functors/equalable_functors.hh"
#include "functors/keyable_functors.hh"
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
using functors::equalable_functor;
using functors::equal_functors;
using functors::keyable_compare_functor;
using functors::keyable_split_functor;
using functors::key_compare_functors;
using functors::key_split_functors;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over key trees.
 */
template <typename T, typename K, typename Syn>
class key_tree_iterator;

template <typename T, typename K, typename Syn>
class key_tree_iterator_reverse;

/*
 * Key trees.
 */
template <typename T, typename K, typename Syn = unsynchronized>
class key_tree : public abstract::multiset<T>,
                 protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class key_tree_iterator<T,K,Syn>;
   friend class key_tree_iterator_reverse<T,K,Syn>;

   /*
    * Define the iterator types.
    */
   typedef key_tree_iterator<T,K,Syn> iterator_t;
   typedef key_tree_iterator_reverse<T,K,Syn> iterator_reverse_t;
   
   /*
    * Constructor.
    * Return an empty tree.
    * Optionally specify custom functors for tree operations.
    */
   explicit key_tree(
      const equalable_functor<T>&         = equal_functors<T>::f_equal(),
      const keyable_compare_functor<T,K>& = (key_compare_functors<T,K>::f_key_compare()),
      const keyable_split_functor<T,K>&   = (key_split_functors<T,K>::f_key_split())
   );

   /*
    * Constructor.
    * Create a tree from the given collection.
    * Optionally specify custom functors for tree operations.
    */
   explicit key_tree(
      const collection<T>&,
      const equalable_functor<T>&         = equal_functors<T>::f_equal(),
      const keyable_compare_functor<T,K>& = (key_compare_functors<T,K>::f_key_compare()),
      const keyable_split_functor<T,K>&   = (key_split_functors<T,K>::f_key_split())
   );

   /*
    * Copy constructor.
    */
   key_tree(const key_tree<T,K,Syn>&);
   
   /*
    * Destructor.
    */
   virtual ~key_tree();

   /*
    * Add element(s) after tree construction.
    * WARNING: This may cause the tree to become unbalanced.
    */
   virtual key_tree<T,K,Syn>& add(T&);
   virtual key_tree<T,K,Syn>& add(const collection<T>&);
   
   /*
    * Remove element(s) after tree construction.
    * WARNING: This may cause the tree to become unbalanced.
    */
   virtual key_tree<T,K,Syn>& remove(const T&);
   virtual key_tree<T,K,Syn>& remove(const collection<T>&);
   
   /*
    * Remove all instances of the element(s).
    * WARNING: This may cause the tree to become unbalanced.
    */
   virtual key_tree<T,K,Syn>& remove_all(const T&);
   virtual key_tree<T,K,Syn>& remove_all(const collection<T>&);
   
   /*
    * Remove all element(s) from the tree.
    * Return a reference to the tree.
    */
   virtual key_tree<T,K,Syn>& clear();

   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not
    * contained in the tree.
    */
   bool contains(const T&) const;
   T& find(const T&) const;

   /*
    * Search.
    * Return a list of all element(s) matching the given element.
    */
   list<T> find_all(const T&) const;
   
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
    * Key tree data structures.
    ************************************************************************/

   /*
    * Key data stored at an internal node.
    */
   class tree_key {
   public:
      /*
       * Constructors.
       */
      explicit tree_key(const K& key) : k(key) { }
      explicit tree_key(const tree_key& t_key) : k(t_key.k) { }

      /*
       * Destructor.
       */
      ~tree_key() { /* do nothing */ }
      
      /*
       * Key.
       */
      K k;
   };
   
   /*
    * Node in the tree.
    */
   class tree_node {
   public:
      /*
       * Constructor.
       * Create an internal node containing the given key.
       */
      explicit tree_node(auto_ptr<tree_key> t_key)
       : key(t_key), items(), parent(NULL), left(NULL), right(NULL), size(0)
      { }

      /*
       * Constructor.
       * Create a leaf node containing the given item.
       */
      explicit tree_node(T& t)
       : key(), items(), parent(NULL), left(NULL), right(NULL), size(1)
      { items.add(t); }
      
      /*
       * Constructor.
       * Create a leaf node containing the given items.
       */
      explicit tree_node(const list<T>& t_items)
       : key(), items(t_items), parent(NULL), left(NULL), right(NULL),
         size(items.size())
      { }
      
      /*
       * Copy constructor.
       * Create a node containing the same key or items.
       */
      explicit tree_node(const tree_node& t_node)
       : key((t_node.key.get() == NULL) ?
            NULL : (new tree_key(*(t_node.key)))
         ),
         items(t_node.items),
         parent(NULL),
         left(NULL),
         right(NULL), 
         size(items.size())
      { }

      /*
       * Destructor.
       */
      ~tree_node() { /* do nothing */ }
      
      /*
       * Node data.
       */
      auto_ptr<tree_key> key;       /* key (if internal node) */
      list<T>            items;     /* items (if leaf node) */
      tree_node*         parent;    /* parent node */
      tree_node*         left;      /* left subtree  */
      tree_node*         right;     /* right subtree */
      unsigned long      size;      /* total # of items in node and subtrees */
   };
   
   /*
    * Tree data.
    */
   const equalable_functor<T>&         _f_equal;         /* item equality test function */
   const keyable_compare_functor<T,K>& _f_key_compare;   /* key comparison function */
   const keyable_split_functor<T,K>&   _f_key_split;     /* split function */
   tree_node*                          _root;            /* root node of tree */

   /************************************************************************
    * Key tree helper functions.
    ************************************************************************/

   /*
    * Swap the contents of two key trees.
    * The trees should have identical functors.
    * Note: Key trees are not locked by this function.
    */
   static void swap(key_tree<T,K,Syn>&, key_tree<T,K,Syn>&);

   /*
    * Build a key tree from the given collection.
    */
   void build(const collection<T>& c);

   /*
    * Runnable object for recursively constructing a key tree.
    *
    * Building the key tree takes O(n + (n/p)*log(n/p)) time, where 
    * n is the tree size and p is the number of available processors.
    */
   class builder : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit builder(
         tree_node*&,                           /* root of subtree */
         auto_ptr< list<T> >,                   /* items to place in subtree */
         const keyable_compare_functor<T,K>&,   /* key comparison function */
         const keyable_split_functor<T,K>&      /* split function */
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
      tree_node*&                         _root;
      auto_ptr< list<T> >                 _items;
      const keyable_compare_functor<T,K>& _f_key_compare;
      const keyable_split_functor<T,K>&   _f_key_split;
   };

   /*
    * Runnable object for recursively copying a key tree.
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
    * Runnable object for recursively destructing a key tree.
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

   /*
    * Add an item to the key tree.
    */
   void add_item(T&);

   /*
    * Remove an item from the key tree.
    */
   void remove_item(const T&);
    
   /*
    * Remove all items matching the given item from the key tree.
    */
   void remove_item_all(const T&);
  
   /*
    * Find node in the given tree that potentially contains the specified item.
    * Return NULL if there is no such node.
    */
   static tree_node* find_node(
      tree_node*, 
      const T&, 
      const equalable_functor<T>&, 
      const keyable_compare_functor<T,K>&
   );
};

/*
 * Key tree iterator.
 */
template <typename T, typename K, typename Syn = unsynchronized>
class key_tree_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit key_tree_iterator(const key_tree<T,K,Syn>&);

   /*
    * Copy constructor.
    */
   key_tree_iterator(const key_tree_iterator<T,K,Syn>&);

   /*
    * Destructor.
    */
   virtual ~key_tree_iterator();

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
   const key_tree<T,K,Syn>&               _tree; /* tree being iterated over */
   typename key_tree<T,K,Syn>::tree_node* _curr; /* current leaf of tree */
   auto_ptr< list_iterator<T> >      _item_iter; /* iterator over leaf items */
};

/*
 * Key tree reverse iterator.
 */
template <typename T, typename K, typename Syn = unsynchronized>
class key_tree_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit key_tree_iterator_reverse(const key_tree<T,K,Syn>&);

   /*
    * Copy constructor.
    */
   key_tree_iterator_reverse(const key_tree_iterator_reverse<T,K,Syn>&);

   /*
    * Destructor.
    */
   virtual ~key_tree_iterator_reverse();

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
   const key_tree<T,K,Syn>&                  _tree; /* tree being iterated over */
   typename key_tree<T,K,Syn>::tree_node*    _curr; /* current leaf of tree */
   auto_ptr< list_iterator_reverse<T> > _item_iter; /* iterator over leaf items */
};

/***************************************************************************
 * Key tree iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Lock tree and initialize position.
 */
template <typename T, typename K, typename Syn>
key_tree_iterator<T,K,Syn>::key_tree_iterator(
   const key_tree<T,K,Syn>& k_tree)
 : _tree(k_tree),
   _curr(NULL),
   _item_iter()
{
   _tree.read_lock();
   _curr = _tree._root;
   if (_curr != NULL) {
      while (_curr->left != NULL) 
         _curr = _curr->left;
      _item_iter.reset(new list_iterator<T>(_curr->items));
   }
}

template <typename T, typename K, typename Syn>
key_tree_iterator_reverse<T,K,Syn>::key_tree_iterator_reverse(
   const key_tree<T,K,Syn>& k_tree)
 : _tree(k_tree),
   _curr(NULL),
   _item_iter()
{
   _tree.read_lock();
   _curr = _tree._root;
   if (_curr != NULL) {
      while (_curr->right != NULL) 
         _curr = _curr->right;
      _item_iter.reset(new list_iterator_reverse<T>(_curr->items));
   }
}

/*
 * Copy constructors.
 */
template <typename T, typename K, typename Syn>
key_tree_iterator<T,K,Syn>::key_tree_iterator(
   const key_tree_iterator<T,K,Syn>& i)
 : _tree(i._tree),
   _curr(i._curr),
   _item_iter(i._item_iter)
{
   _tree.read_lock();
}

template <typename T, typename K, typename Syn>
key_tree_iterator_reverse<T,K,Syn>::key_tree_iterator_reverse(
   const key_tree_iterator_reverse<T,K,Syn>& i)
 : _tree(i._tree),
   _curr(i._curr),
   _item_iter(i._item_iter)
{
   _tree.read_lock();
}

/*
 * Destructors.
 * Unlock the tree.
 */
template <typename T, typename K, typename Syn>
key_tree_iterator<T,K,Syn>::~key_tree_iterator() {
   _tree.read_unlock();
}

template <typename T, typename K, typename Syn>
key_tree_iterator_reverse<T,K,Syn>::~key_tree_iterator_reverse() {
   _tree.read_unlock();
}

/*
 * Check if there is another item available.
 */
template <typename T, typename K, typename Syn>
bool key_tree_iterator<T,K,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}

template <typename T, typename K, typename Syn>
bool key_tree_iterator_reverse<T,K,Syn>::has_next() const {
   return (_item_iter.get() != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename K, typename Syn>
T& key_tree_iterator<T,K,Syn>::next() {
   if (_item_iter.get() != NULL) {
      /* get current item */
      T& t = _item_iter->next();
      /* move to next item */
      if (!(_item_iter->has_next())) {
         /* move to the next leaf node - start by moving up once */
         typename key_tree<T,K,Syn>::tree_node* node_prev = _curr;
         typename key_tree<T,K,Syn>::tree_node* node = node_prev->parent;
         /* move up until we can move to a new right node */
         while (node != NULL) {
            if (node->right != node_prev)
               break;
            node_prev = node;
            node = node->parent;
         }
         /* move right once */
         if (node != NULL)
            node = node->right;
         /* move left as many times as possible */
         if (node != NULL) {
            while (node->left != NULL)
               node = node->left;
         }
         /* now at the next leaf */
         _curr = node;
         _item_iter.reset(
            (_curr == NULL) ?
               NULL : (new list_iterator<T>(_curr->items))
         );
      }
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of key_tree"
      );
   }
}

template <typename T, typename K, typename Syn>
T& key_tree_iterator_reverse<T,K,Syn>::next() {
   if (_item_iter.get() != NULL) {
      /* get current item */
      T& t = _item_iter->next();
      /* move to next item */
      if (!(_item_iter->has_next())) {
         /* move to the next leaf node - start by moving up once */
         typename key_tree<T,K,Syn>::tree_node* node_prev = _curr;
         typename key_tree<T,K,Syn>::tree_node* node = node_prev->parent;
         /* move up until we can move to a new left node */
         while (node != NULL) {
            if (node->left != node_prev)
               break;
            node_prev = node;
            node = node->parent;
         }
         /* move left once */
         if (node != NULL)
            node = node->left;
         /* move right as many times as possible */
         if (node != NULL) {
            while (node->right != NULL)
               node = node->right;
         }
         /* now at the next leaf */
         _curr = node;
         _item_iter.reset(
            (_curr == NULL) ?
               NULL : (new list_iterator_reverse<T>(_curr->items))
         );
      }
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past start of key_tree"
      );
   }
}

/***************************************************************************
 * Key tree implementation.
 ***************************************************************************/

/*
 * Swap the contents of two key trees.
 * The trees should have identical functors.
 * Note: Key trees are not locked by this function.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::swap(
   key_tree<T,K,Syn>& k_tree0,
   key_tree<T,K,Syn>& k_tree1)
{
   tree_node* temp_root = k_tree0._root;
   k_tree0._root = k_tree1._root;
   k_tree1._root = temp_root;
}

/*
 * Constructor.
 * Return an empty tree.
 * Use the specified functors for tree operations.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::key_tree(
   const equalable_functor<T>&         f_eq,
   const keyable_compare_functor<T,K>& f_key_c,
   const keyable_split_functor<T,K>&   f_key_s) 
 : Syn(), 
   _f_equal(f_eq), 
   _f_key_compare(f_key_c), 
   _f_key_split(f_key_s),
   _root(NULL)
{ }

/*
 * Constructor.
 * Create a tree from the given collection.
 * Use the specified functors for tree operations.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::key_tree(
   const collection<T>&                c,
   const equalable_functor<T>&         f_eq,
   const keyable_compare_functor<T,K>& f_key_c,
   const keyable_split_functor<T,K>&   f_key_s)
 : Syn(), 
   _f_equal(f_eq), 
   _f_key_compare(f_key_c), 
   _f_key_split(f_key_s), 
   _root(NULL)
{
   /* build tree */
   this->build(c);
}

/*
 * Copy constructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::key_tree(const key_tree<T,K,Syn>& k_tree) 
 : Syn(),
   _f_equal(k_tree._f_equal),
   _f_key_compare(k_tree._f_key_compare),
   _f_key_split(k_tree._f_key_split),
   _root(NULL)
{
   /* lock source tree */
   auto_read_lock<const Syn> rlock(k_tree);
   /* create temporary tree */
   key_tree<T,K,Syn> k_tree_copy(
      _f_equal,
      _f_key_compare,
      _f_key_split
   );
   /* copy tree */
   copier c(k_tree._root, k_tree_copy._root);
   c.run();
   /* take ownership of copy */
   key_tree<T,K,Syn>::swap(*this, k_tree_copy);
}

/*
 * Build a key tree from the given collection.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::build(const collection<T>& c) {
   /* create temporary tree */
   key_tree<T,K,Syn> k_tree(
      _f_equal,
      _f_key_compare,
      _f_key_split
   );
   /* build tree */
   auto_ptr< list<T> > items(new list<T>(c));
   builder b(
      k_tree._root,
      items,
      k_tree._f_key_compare,
      k_tree._f_key_split
   );
   b.run();
   /* take ownership of tree */
   key_tree<T,K,Syn>::swap(*this, k_tree);
}

/*
 * Builder - Constructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::builder::builder(
   tree_node*&                         root,
   auto_ptr< list<T> >                 items,
   const keyable_compare_functor<T,K>& f_key_compare,
   const keyable_split_functor<T,K>&   f_key_split)
 : _root(root),
   _items(items),
   _f_key_compare(f_key_compare),
   _f_key_split(f_key_split)
{ }

/*
 * Builder - Destructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::builder::~builder() {
   /* do nothing */
}

/*
 * Builder - Build the tree.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::builder::run() {
   unsigned long n_items = _items->size();
   if (n_items == 0) {
      /* empty tree */
      _root = NULL;
   } else if (n_items == 1) {
      /* leaf node containing single element */
      _root = new tree_node(*_items);
   } else {
      /* compute the split key */
      auto_ptr<tree_key> key(new tree_key(_f_key_split(*_items)));
      /* split items using key */
      auto_ptr< list<T> > items_left(new list<T>());
      auto_ptr< list<T> > items_right(new list<T>());
      while (!(_items->is_empty())) {
         T& t = _items->remove_head();
         if (_f_key_compare(t, key->k) < 0)
            items_left->add(t);
         else
            items_right->add(t);
      }
      /* delete original list */
      _items.reset();
      /* check if items were split */
      if (items_left->is_empty()) {
         /* items could not be split */
         _root = new tree_node(*items_right);
      } else if (items_right->is_empty()) {
         /* items could not be split */
         _root = new tree_node(*items_left);
      } else {
         /* build the split node */
         _root = new tree_node(key);
         /* recursively build the subtrees */
         builder b_left(
            _root->left,
            items_left,
            _f_key_compare,
            _f_key_split
         );
         builder b_right(
            _root->right,
            items_right,
            _f_key_compare,
            _f_key_split
         );
         child_thread::run(b_left, b_right);
         /* link subtrees to root */
         if (_root->left != NULL) {
            _root->left->parent = _root;
            _root->size += _root->left->size;
         }
         if (_root->right != NULL) {
            _root->right->parent = _root;
            _root->size += _root->right->size;
         }
      }
   }
}

/*
 * Copier - Constructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::copier::copier(tree_node* src, tree_node*& dest)
 : _src(src),
   _dest(dest)
{ }

/*
 * Copier - Destructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::copier::~copier() {
   /* do nothing */
}

/*
 * Copier - Copy the tree.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::copier::run() {
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
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::~key_tree() {
   destroyer d(_root);
   d.run();
}

/*
 * Destroyer - Constructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::destroyer::destroyer(tree_node* root)
 : _root(root)
{ }

/*
 * Destroyer - Destructor.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>::destroyer::~destroyer() {
   /* do nothing */
}

/*
 * Destroyer - Destroy the tree.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::destroyer::run() {
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
 * Add element after tree construction.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::add(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->add_item(t);
   return *this;
}

/*
 * Add element(s) after tree construction.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::add(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->add_item(i.next());
   return *this;
}

/*
 * Remove element after tree construction.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

/*
 * Remove element(s) after tree construction.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Remove all instances of the element.
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::remove_all(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item_all(t);
   return *this;
}

/*
 * Remove all instances of the element(s).
 * WARNING: This may cause the tree to become unbalanced.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::remove_all(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item_all(i.next());
   return *this;
}

/*
 * Remove all element(s) from the tree.
 * Return a reference to the tree.
 */
template <typename T, typename K, typename Syn>
key_tree<T,K,Syn>& key_tree<T,K,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   key_tree<T,K,Syn> k_tree(
      _f_equal,
      _f_key_compare,
      _f_key_split
   );
   key_tree<T,K,Syn>::swap(*this, k_tree);
   return *this;
}

/*
 * Add an item to the key tree.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::add_item(T& t) {
   if (_root == NULL) {
      /* create new leaf node for item */
      _root = new tree_node(t);
   } else {
      /* descend to leaf node closest to item */
      tree_node* curr = _root;
      while (curr->key.get() != NULL) {
         if (_f_key_compare(t, curr->key->k) < 0) 
            curr = curr->left;
         else
            curr = curr->right;
      }
      /* gather items at current node */
      list<T> items;
      items.add(curr->items);
      items.add(t);
      /* build tree for items */
      key_tree k_tree(
         items,
         _f_equal,
         _f_key_compare,
         _f_key_split
      );
      /* swap current node with tree */
      tree_node* parent = curr->parent;
      if (parent == NULL) {
         /* current node is root node */
         _root = k_tree._root;
      } else if (curr == parent->left) {
         /* replace left node */
         parent->left = k_tree._root;
         parent->left->parent = parent;
      } else {
         /* replace right node */
         parent->right = k_tree._root;
         parent->right->parent = parent;
      }
      curr->parent = NULL;
      k_tree._root = curr;
      /* update size along path to root */
      while (parent != NULL) {
         parent->size++;
         parent = parent->parent;
      }
   }
}

/*
 * Remove an item from the key tree.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::remove_item(const T& t) {
   /* find node potentially containing item */
   tree_node* t_node = key_tree<T,K,Syn>::find_node(
      _root, 
      t, 
      _f_equal, 
      _f_key_compare
   );
   /* remove the item */
   if (t_node != NULL) {
      /* remove item from leaf node */
      bool found = false;
      unsigned long n_items = t_node->items.size();
      for (unsigned long n = 0; (n < n_items) && (!found); n++) {
         T& item = t_node->items.remove_head();
         if (_f_equal(t, item)) {
            found = true;
         } else {
            t_node->items.append(item);
         }
      }
      /* remove leaf node if empty */
      tree_node* t_node_parent = t_node->parent;
      if (t_node->items.is_empty()) {
         if (t_node_parent != NULL) {
            /* move sibling up the tree */
            tree_node* t_node_sib = 
               (t_node_parent->left == t_node) ? 
                  (t_node_parent->right) : 
                  (t_node_parent->left);
            t_node_parent->key   = t_node_sib->key;
            t_node_parent->items.add(t_node_sib->items);
            t_node_parent->left  = t_node_sib->left;
            t_node_parent->right = t_node_sib->right;
            delete t_node_sib;
            if (t_node_parent->left != NULL)
               t_node_parent->left->parent = t_node_parent;
            if (t_node_parent->right != NULL)
               t_node_parent->right->parent = t_node_parent;
         } else {
            /* single element tree - node was root */
            _root = NULL;
         }
         /* delete the node */
         delete t_node;
      }
      /* update size along path to root */
      while (t_node_parent != NULL) {
         t_node_parent->size--;
         t_node_parent = t_node_parent->parent;
      }
   }
}

/*
 * Remove all items matching the given item from the key tree.
 */
template <typename T, typename K, typename Syn>
void key_tree<T,K,Syn>::remove_item_all(const T& t) {
   /* find node potentially containing item */
   tree_node* t_node = key_tree<T,K,Syn>::find_node(
      _root, 
      t, 
      _f_equal, 
      _f_key_compare
   );
   /* remove the item(s) */
   if (t_node != NULL) {
      /* remove item(s) from leaf node */
      unsigned long n_items = t_node->items.size();
      unsigned long n_removed = 0;
      for (unsigned long n = 0; n < n_items; n++) {
         T& item = t_node->items.remove_head();
         if (_f_equal(t, item)) {
            n_removed++;
         } else {
            t_node->items.append(item);
         }
      }
      /* remove leaf node if empty */
      tree_node* t_node_parent = t_node->parent;
      if (t_node->items.is_empty()) {
         if (t_node_parent != NULL) {
            /* move sibling up the tree */
            tree_node* t_node_sib = 
               (t_node_parent->left == t_node) ? 
                  (t_node_parent->right) : 
                  (t_node_parent->left);
            t_node_parent->key   = t_node_sib->key;
            t_node_parent->items.add(t_node_sib->items);
            t_node_parent->left  = t_node_sib->left;
            t_node_parent->right = t_node_sib->right;
            delete t_node_sib;
            if (t_node_parent->left != NULL)
               t_node_parent->left->parent = t_node_parent;
            if (t_node_parent->right != NULL)
               t_node_parent->right->parent = t_node_parent;
         } else {
            /* single element tree - node was root */
            _root = NULL;
         }
         /* delete the node */
         delete t_node;
      }
      /* update size along path to root */
      while (t_node_parent != NULL) {
         t_node_parent->size -= n_removed;
         t_node_parent = t_node_parent->parent;
      }
   }
}

/*
 * Find node in the given tree that potentially contains the specified item.
 * Return NULL if there is no such node.
 */
template <typename T, typename K, typename Syn>
typename key_tree<T,K,Syn>::tree_node* key_tree<T,K,Syn>::find_node(
   tree_node* t_node, 
   const T& t, 
   const equalable_functor<T>& f_equal, 
   const keyable_compare_functor<T,K>& f_key_compare)
{
   if (t_node == NULL) {
      /* empty tree - no leaf nodes */
      return NULL;
   } else if (t_node->key.get() != NULL) {
      /* internal node - recurse on subtree */
      if (f_key_compare(t, t_node->key->k) < 0) { 
         return (key_tree<T,K,Syn>::find_node(
            t_node->left,  t, f_equal, f_key_compare));
      } else {
         return (key_tree<T,K,Syn>::find_node(
            t_node->right, t, f_equal, f_key_compare));
      }
   } else {
      /* leaf node */
      return t_node;
   }
}

/*
 * Check if the tree contains the specified item.
 */
template <typename T, typename K, typename Syn>
bool key_tree<T,K,Syn>::contains(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   /* find node potentially containing item */
   tree_node* t_node = key_tree<T,K,Syn>::find_node(
      _root, 
      t, 
      _f_equal, 
      _f_key_compare
   );
   /* check if node actually contains item */
   bool found = false;
   if (t_node != NULL) {
      for (typename list<T>::iterator_t i(t_node->items);
           (i.has_next()) && (!found); )
         found = _f_equal(t, i.next());
   }
   return found;
}

/*
 * Find the specified item in the tree.
 * Throw an exception (ex_not_found) when attempting to find an element not
 * contained in the tree.
 */
template <typename T, typename K, typename Syn>
T& key_tree<T,K,Syn>::find(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   /* find node potentially containing item */
   tree_node* t_node = key_tree<T,K,Syn>::find_node(
      _root, 
      t, 
      _f_equal, 
      _f_key_compare
   );
   /* find item in node */
   if (t_node != NULL) {
      for (typename list<T>::iterator_t i(t_node->items); i.has_next(); ) {
         T& item = i.next();
         if (_f_equal(t, item))
            return item;
      }
   }
   throw ex_not_found(
      "item not found in key_tree"
   );
}

/*
 * Return a list of all element(s) matching the given element.
 */
template <typename T, typename K, typename Syn>
list<T> key_tree<T,K,Syn>::find_all(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   /* find node potentially containing item */
   tree_node* t_node = key_tree<T,K,Syn>::find_node(
      _root, 
      t, 
      _f_equal, 
      _f_key_compare
   );
   /* find item in node */
   list<T> items;
   if (t_node != NULL) {
      for (typename list<T>::iterator_t i(t_node->items); i.has_next(); ) {
         T& item = i.next();
         if (_f_equal(t, item))
            items.add(item);
      }
   }
   return items;
}

/*
 * Get number of elements in tree.
 */
template <typename T, typename K, typename Syn>
unsigned long key_tree<T,K,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return (_root != NULL) ? (_root->size) : 0;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename K, typename Syn>
auto_ptr< iterator<T> >key_tree<T,K,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(
      new key_tree_iterator<T,K,Syn>(*this)
   );
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename K, typename Syn>
auto_ptr< iterator<T> > key_tree<T,K,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(
      new key_tree_iterator_reverse<T,K,Syn>(*this)
   );
}

} /* namespace collections */

#endif
