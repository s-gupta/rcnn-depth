/*
 * Priority queues (thread-safe).
 *
 * Smaller priorities are placed at the front of the queue.
 */
#ifndef COLLECTIONS__QUEUE_HH
#define COLLECTIONS__QUEUE_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/queue.hh"
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
 * Declare class for iterator over queues.
 */
template <typename T, typename Syn>
class queue_iterator;

/*
 * Priority queues.
 */
template <typename T, typename Syn = unsynchronized>
class queue : public abstract::queue<T>,
              protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class queue_iterator<T,Syn>;

   /*
    * Define the iterator type.
    */
   typedef queue_iterator<T,Syn> iterator_t;
   
   /*
    * Constructor.
    */
   explicit queue(
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Constructor.
    * Create a queue from a collection.
    */
   explicit queue(
      const collection<T>&, 
      const comparable_functor<T>& = compare_functors<T>::f_compare()
   );

   /*
    * Copy constructor.
    */
   queue(const queue<T,Syn>&);
   
   /*
    * Destructor.
    */ 
   virtual ~queue();

   /*
    * Collection interface for enqueueing elements.
    * Return a reference to the queue.
    */
   queue<T,Syn>& add(T&);
   queue<T,Syn>& add(const collection<T>&);
   
   /*
    * Enqueue.
    * Return a reference to the queue.
    */
   queue<T,Syn>& enqueue(T&);
   queue<T,Syn>& enqueue(const collection<T>&);

   /*
    * Dequeue.
    * Throw an exception (ex_not_found) if the queue is empty.
    */
   T& dequeue();

   /*
    * Return top element of queue without removing it.
    * Throw an exception (ex_not_found) if the queue is empty.
    */
   T& head() const;

   /*
    * Remove all element(s) from the queue.
    * Return a reference to the queue.
    */
   queue<T,Syn>& clear();

   /*
    * Size.
    */
   unsigned long size() const;
   
   /*
    * Return iterator over elements.
    */
   auto_ptr< iterator<T> > iter_create() const;
 
protected:
   /************************************************************************
    * Queue data structures.
    ************************************************************************/

   /*
    * Declare node class.
    */
   class queue_node;
   
   /*
    * Item in the heap.
    */
   class queue_item {
   public:
      /*
       * Constructors.
       */
      explicit queue_item(T&);
      explicit queue_item(const queue_item&);

      /*
       * Destructor.
       */
      ~queue_item();

      /*
       * Item data.
       */
      T& t;                /* reference to item */
      queue_node* q_node;  /* node containing queue_item */
   };
  
   /* 
    * Node in the heap.
    */
   class queue_node {
   public:
      /* 
       * Constructors.
       */
      explicit queue_node(auto_ptr<queue_item>);
      explicit queue_node(const queue_node&);

      /*
       * Destructor.
       */
      ~queue_node();
      
      /*
       * Node data.
       */
      auto_ptr<queue_item> q_item;  /* item stored in the node */
      queue_node* parent;           /* tree structure */
      queue_node* left;
      queue_node* right;
      queue_node* prev;             /* list structure */
      queue_node* next;
   };

   /*
    * Queue data structure.
    */
   const comparable_functor<T>& _f_compare;  /* comparison function */
   unsigned long _size;                      /* number of elements in queue */
   queue_node* _root;                        /* heap root */
   queue_node* _leaf;                        /* bottom, rightmost leaf node */

   /************************************************************************
    * Queue helper functions.
    ************************************************************************/
    
   /*
    * Swap the contents of two queues.
    * The queues must have the same comparison functor.
    * Note: queues are not locked by this function.
    */
   static void swap(queue<T,Syn>&, queue<T,Syn>&);
    
   /*
    * Runnable object for recursively destructing a queue.
    *
    * Destroying the queue takes O(n/p + log(n)) time where n is the queue 
    * size and p is the number of available processors.
    */
   class destroyer : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit destroyer(queue_node*);

      /*
       * Destructor.
       */
      virtual ~destroyer();

      /*
       * Destroy the queue.
       */
      virtual void run();
      
   protected:
      queue_node* _root;
   };
    
   /*
    * Heapify routines (swap up/down).
    * Swap up/down starting from the given node.
    * Return the final location of the start item.
    */
   queue_node* heapify_up(queue_node*);
   queue_node* heapify_down(queue_node*);

   /*
    * Enqueue an item into the heap.
    * Return a pointer to heap node containing the item.
    */
   queue_node* enqueue_heap(auto_ptr<queue_item>);

   /*
    * Remove an arbitrary node from the heap.
    * Heapify up/down to maintain the heap invariant.
    */
   void remove_heap_node(queue_node*);
};

/*
 * Queue iterator.
 */
template <typename T, typename Syn = unsynchronized>
class queue_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    * Seek to the head of the queue.
    */
   explicit queue_iterator(const queue<T,Syn>&);

   /*
    * Copy constructor.
    */
   queue_iterator(const queue_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~queue_iterator();

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
   const queue<T,Syn>& _q;                   /* queue being iterated over */
   typename queue<T,Syn>::queue_node* _curr; /* current position in queue */
};

/************************************************************************
 * Queue iterator implementation.
 ************************************************************************/

/*
 * Constructor.
 * Lock queue and initialize position.
 */
template <typename T, typename Syn>
queue_iterator<T,Syn>::queue_iterator(const queue<T,Syn>& q)
 : _q(q)
{
   _q.read_lock();
   _curr = _q._root;
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
queue_iterator<T,Syn>::queue_iterator(const queue_iterator<T,Syn>& i)
 : _q(i._q),
   _curr(i._curr)
{
   _q.read_lock();
}

/*
 * Destructor.
 * Unlock the queue.
 */
template <typename T, typename Syn>
queue_iterator<T,Syn>::~queue_iterator() {
   _q.read_unlock();
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool queue_iterator<T,Syn>::has_next() const {
   return (_curr != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& queue_iterator<T,Syn>::next() {
   if (_curr != NULL) {
      T& t = _curr->q_item->t;
      _curr = _curr->next;
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past end of queue"
      );
   }
}

/************************************************************************
 * Queue item helper class implementation.
 ************************************************************************/

/*
 * Constructor.
 * Create an item referencing the specified object.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue_item::queue_item(T& item)
 : t(item),
   q_node(NULL)
{ }

/*
 * Copy constructor.
 * Create an item referencing the same object as the original item.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue_item::queue_item(const queue_item& q_item)
 : t(q_item.t),
   q_node(NULL)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue_item::~queue_item() { 
   /* do nothing */
}

/************************************************************************
 * Queue node helper class implementation.
 ************************************************************************/

/*
 * Constructor.
 * Create a node containing the given item.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue_node::queue_node(auto_ptr<queue_item> item)
 : q_item(item),
   parent(NULL),
   left(NULL),
   right(NULL),
   prev(NULL),
   next(NULL)
{ 
   /* store backpointer from item to node */
   if (q_item.get() != NULL) 
      q_item->q_node = this;
}

/*
 * Copy constructor.
 * Create a node identical to the given node.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue_node::queue_node(const queue_node& node)
 : q_item(
      (node.q_item.get() == NULL) ? 
         NULL : (new queue_item(*(node.q_item)))
   ), 
   parent(node.parent),
   left(node.left),
   right(node.right),
   prev(node.prev),
   next(node.next)
{
   /* store backpointer from item to node */
   if (q_item.get() != NULL) 
      q_item->q_node = this;
}

/*
 * Destructor.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue_node::~queue_node() {
   /* do nothing */
}

/************************************************************************
 * Queue implementation.
 ************************************************************************/

/*
 * Swap the contents of two queues.
 * The queues must have the same comparison functor.
 * Note: queues are not locked by this function.
 */
template <typename T, typename Syn>
void queue<T,Syn>::swap(queue<T,Syn>& q0, queue<T,Syn>& q1) {
   unsigned long temp_size = q0._size;
   queue_node*   temp_root = q0._root;
   queue_node*   temp_leaf = q0._leaf;
   q0._size = q1._size;
   q0._root = q1._root;
   q0._leaf = q1._leaf;
   q1._size = temp_size;
   q1._root = temp_root;
   q1._leaf = temp_leaf;
}

/*
 * Constructor.
 * Specify the comparison function to use when prioritizing queue elements.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue(const comparable_functor<T>& f)
 : Syn(), 
   _f_compare(f),
   _size(0),
   _root(NULL),
   _leaf(NULL)
{ }

/*
 * Constructor.
 * Create a queue from a collection.
 * Specify the comparison function to use when prioritizing queue elements.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue(const collection<T>& c, const comparable_functor<T>& f)
 : Syn(),
   _f_compare(f),
   _size(0),
   _root(NULL),
   _leaf(NULL)
{
   /* enqueue elements into queue */
   queue<T,Syn> q(_f_compare);
   for (auto_ptr< iterator<T> > i = c.iter_create(); i->has_next(); ) {
      auto_ptr<queue_item> q_item(new queue_item(i->next()));
      q.enqueue_heap(q_item);
   }
   /* take ownership of queue */
   queue<T,Syn>::swap(*this, q);
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
queue<T,Syn>::queue(const queue<T,Syn>& q)
 : Syn(), 
   _f_compare(q._f_compare), 
   _size(0),
   _root(NULL),
   _leaf(NULL)
{
   /* initialize copy */
   queue<T,Syn> q_copy(_f_compare);
   /* lock source queue */
   auto_read_lock<const Syn> rlock(q);
   /*
    * Simply enqueue everything from the old queue.
    * This is an O(n) procedure as no swaps occur.
    */
   queue_node* curr = q._root;
   while (curr != NULL) {
      auto_ptr<queue_item> q_item(new queue_item(curr->q_item->t));
      q_copy.enqueue_heap(q_item);
      curr = curr->next;
   }
   /* take ownership of queue */
   queue<T,Syn>::swap(*this, q_copy);
}

/*
 * Destructor.
 * Delete heap nodes.
 */
template <typename T, typename Syn>
queue<T,Syn>::~queue() {
   destroyer d(_root);
   d.run();
}

/*
 * Destroyer - Constructor.
 */
template <typename T, typename Syn>
queue<T,Syn>::destroyer::destroyer(queue_node* root)
 : _root(root)
{ }

/*
 * Destroyer - Destructor.
 */
template <typename T, typename Syn>
queue<T,Syn>::destroyer::~destroyer() {
   /* do nothing */
}

/*
 * Destroyer - Destroy the queue.
 */
template <typename T, typename Syn>
void queue<T,Syn>::destroyer::run() {
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
 * Heapify up starting from the specified node.
 * Return the final location of the start item.
 */
template <typename T, typename Syn>
typename queue<T,Syn>::queue_node* queue<T,Syn>::heapify_up(queue_node* curr)
{
   /* get item at starting node */
   auto_ptr<queue_item> q_item = curr->q_item;
   T& t = q_item->t;
   /* swap upwards */
   queue_node* parent = curr->parent;
   while (parent != NULL) {
      if (_f_compare(t, parent->q_item->t) < 0) {
         curr->q_item = parent->q_item;
         curr->q_item->q_node = curr;
         curr = parent;
         parent = curr->parent;
      } else {
         break;
      }
   }
   /* insert starting item */
   curr->q_item = q_item;
   curr->q_item->q_node = curr;
   return curr;
}

/*
 * Heapify down starting from the specified node.
 * Return the final location of the start item.
 */
template <typename T, typename Syn>
typename queue<T,Syn>::queue_node* queue<T,Syn>::heapify_down(queue_node* curr)
{
   /* get item at starting node */
   auto_ptr<queue_item> q_item = curr->q_item;
   T& t = q_item->t;
   /* initialize children */
   queue_node* left  = curr->left;
   queue_node* right = curr->right;
   /* handle all two-child nodes */
   while (right != NULL) {
      queue_item* left_item  = left->q_item.get();
      queue_item* right_item = right->q_item.get();
      bool l_less_than_r = (_f_compare(left_item->t, right_item->t) < 0);
      if (l_less_than_r && (_f_compare(left_item->t, t) < 0)) {
         /* swap left */
         curr->q_item = left->q_item;
         curr->q_item->q_node = curr;
         curr = left;
     } else if ((!l_less_than_r) && (_f_compare(right_item->t, t) < 0)) {
         /* swap right */
         curr->q_item = right->q_item;
         curr->q_item->q_node = curr;
         curr = right;
      } else {
         break;
      }
      left  = curr->left;
      right = curr->right;
   }
   /* handle the one-child node if it exists */
   if (left != NULL) {
      queue_item* left_item = left->q_item.get();
      if (_f_compare(left_item->t, t) < 0) {
         /* swap left */
         curr->q_item = left->q_item;
         curr->q_item->q_node = curr;
         curr = left;
      }
   }
   /* insert starting item */
   curr->q_item = q_item;
   curr->q_item->q_node = curr;
   return curr;
}

/*
 * Enqueue an item into the heap.
 * Return a pointer to heap node containing the item.
 */
template <typename T, typename Syn>
typename queue<T,Syn>::queue_node* queue<T,Syn>::enqueue_heap(
   auto_ptr<queue_item> q_item)
{
   /* allocate the new leaf node */
   queue_node* leaf_new = new queue_node(q_item);
   leaf_new->prev = _leaf;
   /* initialize pointer to new leaf's parent */
   queue_node* leaf_parent = NULL;
   /* check whether heap is empty */
   if (_leaf != NULL) {
      /* non-empty heap - add end node */
      _leaf->next = leaf_new;
      leaf_parent = _leaf->parent;
      if (leaf_parent != NULL) {
         /* old leaf is not the root */
         if (leaf_parent->right != NULL) {
            /* no more room at parent */
            leaf_parent = leaf_parent->next;
            leaf_parent->left = leaf_new;
         } else {
            /* more room at parent */
            leaf_parent->right = leaf_new;
         }
         leaf_new->parent = leaf_parent;
      } else {
         /* leaf is the root node */
         _leaf->left = leaf_new;
         leaf_new->parent = _leaf;
      }
   } else {
      /* empty heap - add first node */
      _root = leaf_new;
   }
   /* update heap to reflect new leaf */
   _leaf = leaf_new;
   _size++;
   /* swap upwards as needed */
   return (this->heapify_up(leaf_new));
}

/*
 * Remove an arbitrary node from the heap.
 * Heapify up/down to maintain the heap invariant.
 */
template <typename T, typename Syn>
void queue<T,Syn>::remove_heap_node(queue_node* curr) {
   /* check if queue becomes empty */
   if (_size == 1) {
      /* current node must be root and leaf */
      delete curr;
      _size = 0;
      _root = NULL;
      _leaf = NULL;
   } else {
      /* check if current node is leaf */
      bool is_leaf = (curr == _leaf);
      /* swap leaf item into current position */
      if (!is_leaf) {
         curr->q_item = _leaf->q_item;
         curr->q_item->q_node = curr;
      }
      /* remove the leaf */
      queue_node* leaf_parent = _leaf->parent;
      queue_node* leaf_prev   = _leaf->prev;
      if (leaf_parent->right != NULL)
         leaf_parent->right = NULL;
      else
         leaf_parent->left = NULL;
      leaf_prev->next = NULL;
      delete _leaf;
      _leaf = leaf_prev;
      /* adjust queue size */
      _size--;
      /* heapify */
      if (!is_leaf) {
         this->heapify_up(curr);
         this->heapify_down(curr);
      }
   }
}

/*
 * Collection interface for enqueueing an element.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue<T,Syn>& queue<T,Syn>::add(T& t) {
   return this->enqueue(t);
}

/*
 * Collection interface for enqueueing elements.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue<T,Syn>& queue<T,Syn>::add(const collection<T>& c) {
   return this->enqueue(c);
}

/*
 * Enqueue.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue<T,Syn>& queue<T,Syn>::enqueue(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   auto_ptr<queue_item> q_item(new queue_item(t));
   this->enqueue_heap(q_item);
   return *this;
}

/*
 * Enqueue all elements in collection.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue<T,Syn>& queue<T,Syn>::enqueue(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); ) {
      auto_ptr<queue_item> q_item(new queue_item(i.next()));
      this->enqueue_heap(q_item);
   }
   return *this;
}
 
/*
 * Dequeue.
 * Throw an exception (ex_not_found) if the queue is empty.
 */
template <typename T, typename Syn>
T& queue<T,Syn>::dequeue() {
   auto_write_lock<const Syn> wlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to dequeue from empty queue"
      );
   T& t = _root->q_item->t;
   this->remove_heap_node(_root);
   return t;
}

/*
 * Return top element of queue without removing it.
 * Throw an exception (ex_not_found) if the queue is empty.
 */
template <typename T, typename Syn>
T& queue<T,Syn>::head() const {
   auto_read_lock<const Syn> rlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to grab head of empty queue"
      );
   return _root->q_item->t;
}

/*
 * Remove all element(s) from the queue.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue<T,Syn>& queue<T,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   queue<T,Syn> q(_f_compare);
   queue<T,Syn>::swap(*this, q);
   return *this;
}

/*
 * Get number of elements in queue.
 */
template <typename T, typename Syn>
unsigned long queue<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > queue<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(new queue_iterator<T,Syn>(*this));
}

} /* namespace collections */

#endif
