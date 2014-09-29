/*
 * Priority queue set (thread-safe).
 *
 * A queue_set is a priority queue that also supports set operations.  In 
 * particular, it is possible to remove an arbitrary element from the queue.
 *
 * Smaller priorities are placed at the front of the queue.
 *
 * Note that two comparison functors are used in construction of a queue set:
 *    (1) a priority comparison functor for the queue and 
 *    (2) a search comparison functor for the set
 *
 * Both queue and set data structures are maintained.
 *
 * Note that queue sets do not allow duplicate items (as defined by the search 
 * comparison functor for the set).  Adding (enqueueing) an item already in 
 * the set has no effect.
 */
#ifndef COLLECTIONS__QUEUE_SET_HH
#define COLLECTIONS__QUEUE_SET_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/queue.hh"
#include "collections/abstract/set.hh"
#include "collections/list.hh"
#include "collections/queue.hh"
#include "collections/set.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
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
using functors::comparable_functor;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;


/*
 * Declare class for iterator over queue sets.
 */
template <typename T, typename Syn>
class queue_set_iterator;

/*
 * Priority queue sets.
 */
template <typename T, typename Syn = unsynchronized>
class queue_set : public abstract::queue<T>,
                  public abstract::set<T>,
                  protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class queue_set_iterator<T,Syn>;

   /*
    * Define the iterator type.
    */
   typedef queue_set_iterator<T,Syn> iterator_t;

   /*
    * Constructor.
    */
   explicit queue_set(
      const comparable_functor<T>&, /* priority comparison functor */
      const comparable_functor<T>&  /* search comparison functor */
   );

   /*
    * Constructor.
    * Create a queue set from a collection.
    */
   explicit queue_set(
      const collection<T>&, 
      const comparable_functor<T>&, /* priority comparison functor */
      const comparable_functor<T>&  /* search comparison functor */
   );

   /*
    * Copy constructor.
    */
   queue_set(const queue_set<T,Syn>&);
   
   /*
    * Destructor.
    */ 
   virtual ~queue_set();

   /*
    * Collection interface for enqueueing elements.
    * Return a reference to the queue.
    */
   queue_set<T,Syn>& add(T&);
   queue_set<T,Syn>& add(const collection<T>&);
   
   /*
    * Enqueue.
    * Return a reference to the queue.
    */
   queue_set<T,Syn>& enqueue(T&);
   queue_set<T,Syn>& enqueue(const collection<T>&);

   /*
    * Update.
    * Update the queue to reflect an element whos priority has changed.
    */
   queue_set<T,Syn>& update(const T&);
   
   /*
    * Remove element(s) from the queue.
    * Return a reference to the queue.
    */
   queue_set<T,Syn>& remove(const T&);
   queue_set<T,Syn>& remove(const collection<T>&);

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
   queue_set<T,Syn>& clear();

   /*
    * Search.
    * Throw an exception (ex_not_found) when attempting to find an element not 
    * contained in the queue.
    */
   bool contains(const T&) const;
   T& find(const T&) const;

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
    * Queue set data structures.
    ************************************************************************/
    
   /*
    * Declare class for set search comparison functor.
    */
   class queue_item_compare_functor;

   /*
    * Declare class that gives queue sets access to internal queue data 
    * structures found in the collections::queue class.
    */
   class queue_set_q : public queue<T> {
   public:
      /*
       * Friend classes.
       */
      friend class queue_set<T,Syn>;
      friend class queue_item_compare_functor;

      /*
       * Constructor.
       */
      explicit queue_set_q(
         const comparable_functor<T>& f)
       : queue<T>(f)
      { }

      /*
       * Constructor.
       * Create a queue from a collection.
       */
      explicit queue_set_q(
         const collection<T>& c,
         const comparable_functor<T>& f)
       : queue<T>(c,f)
      { }

      /*
       * Copy constructor.
       */
      queue_set_q(const queue_set_q& q)
       : queue<T>(q)
      { }
      
      /*
       * Destructor.
       */ 
      virtual ~queue_set_q() { /* do nothing */ }
   };

   /*
    * Set search comparison functor on queue items.
    * This calls the given comparison functor on the elements contained by 
    * the queue items.
    */
   class queue_item_compare_functor
    : public comparable_functor<typename queue_set_q::queue_item> {
   public:
      /*
       * Constructor.
       */
      explicit queue_item_compare_functor(const comparable_functor<T>& f)
       : _f(f) { }
   
      /*
       * Copy constructor.
       */
      explicit queue_item_compare_functor(const queue_item_compare_functor& f)
       : _f(f._f) { }
       
      /*
       * Comparison function.
       */
      int operator()(
         const typename queue_set_q::queue_item& i0,
         const typename queue_set_q::queue_item& i1) const
      {
         return _f(i0.t, i1.t);
      }

   protected:
      const comparable_functor<T>& _f;
   };

   /*
    * Queue set data.
    */
   const comparable_functor<T>&       _f_prio; /* queue priority comparison */
   const queue_item_compare_functor _f_search; /* set search comparison */
   queue_set_q                             _q; /* queue of elements */
   set<typename queue_set_q::queue_item>   _s; /* set of items in queue */

   /************************************************************************
    * Queue set helper functions.
    ************************************************************************/

   /*
    * Enqueue an item.
    */
   void enqueue_item(T&);
   
   /*
    * Remove an item.
    */
   void remove_item(const T&);
};

/*
 * Queue set iterator.
 */
template <typename T, typename Syn = unsynchronized>
class queue_set_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit queue_set_iterator(const queue_set<T,Syn>&);

   /*
    * Copy constructor.
    */
   queue_set_iterator(const queue_set_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~queue_set_iterator();

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
   const queue_set<T,Syn>&   _q_set;   /* queue set being iterated over */
   auto_read_lock<const Syn> _rlock;   /* read lock on queue set */
   queue_iterator<T>         _iter;    /* iterator over elements in queue */
};

/***************************************************************************
 * Queue set iterator implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename Syn>
queue_set_iterator<T,Syn>::queue_set_iterator(const queue_set<T,Syn>& q_set)
 : _q_set(q_set),
   _rlock(_q_set),
   _iter(_q_set._q)
{ }

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
queue_set_iterator<T,Syn>::queue_set_iterator(
   const queue_set_iterator<T,Syn>& i)
 : _q_set(i._q_set),
   _rlock(_q_set),
   _iter(i._iter)
{ }

/*
 * Destructor.
 * Do nothing as the queue set is automatically unlocked upon destruction of 
 * the read lock.
 */
template <typename T, typename Syn>
queue_set_iterator<T,Syn>::~queue_set_iterator() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool queue_set_iterator<T,Syn>::has_next() const {
   return _iter.has_next();
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& queue_set_iterator<T,Syn>::next() {
   return _iter.next();
}

/***************************************************************************
 * Queue set implementation.
 ***************************************************************************/

/*
 * Constructor.
 * Specify the priority and search comparison functors.
 */
template <typename T, typename Syn>
queue_set<T,Syn>::queue_set(
   const comparable_functor<T>& f_prio,
   const comparable_functor<T>& f_search)
 : Syn(),
   _f_prio(f_prio),
   _f_search(f_search),
   _q(_f_prio),
   _s(_f_search)
{ }

/*
 * Constructor.
 * Create a queue set from a collection.
 * Specify the priority and search comparison functors.
 */
template <typename T, typename Syn>
queue_set<T,Syn>::queue_set(
   const collection<T>& c, 
   const comparable_functor<T>& f_prio,
   const comparable_functor<T>& f_search)
 : Syn(),
   _f_prio(f_prio),
   _f_search(f_search),
   _q(c, _f_prio),
   _s(_f_search)
{
   /* add all queue_items (containers of elements in queue) to the set */
   typename queue_set_q::queue_node* curr = _q._root;
   while (curr != NULL) {
      _s.add(*(curr->q_item));
      curr = curr->next;
   }
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
queue_set<T,Syn>::queue_set(const queue_set<T,Syn>& q_set)
 : Syn(),
   _f_prio(q_set._f_prio),
   _f_search(q_set._f_search),
   _q(q_set, _f_prio),
   _s(_f_search)
{
   /* add all queue_items (containers of elements in queue) to the set */
   typename queue_set_q::queue_node* curr = _q._root;
   while (curr != NULL) {
      _s.add(*(curr->q_item));
      curr = curr->next;
   }
}

/*
 * Destructor.
 */ 
template <typename T, typename Syn>
queue_set<T,Syn>::~queue_set() {
   /* do nothing */
}

/*
 * Collection interface for enqueueing elements.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::add(T& t) {
   return this->enqueue(t);
}

template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::add(const collection<T>& c) {
   return this->enqueue(c);
}

/*
 * Enqueue an item.
 */
template <typename T, typename Syn>
void queue_set<T,Syn>::enqueue_item(T& t) {
   auto_ptr<typename queue_set_q::queue_item> q_item(
      new typename queue_set_q::queue_item(t)
   );
   if (!(_s.contains(*q_item))) {
      _s.add(*q_item);
      _q.enqueue_heap(q_item);
   }
}  

/*
 * Enqueue.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::enqueue(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->enqueue_item(t);
   return *this;
}

template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::enqueue(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->enqueue_item(i.next());
   return *this;
}

/*
 * Update.
 * Update the queue to reflect an element whos priority has changed.
 */
template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::update(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   /* note: const_cast is safe (and required) below */
   typename queue_set_q::queue_item q_item(const_cast<T&>(t));
   typename queue_set_q::queue_item& q_item_match = _s.find(q_item);
   typename queue_set_q::queue_node* q_node = q_item_match.q_node;
   _q.heapify_up(q_node);
   _q.heapify_down(q_node);
   return *this;
}

/*
 * Remove an item.
 */
template <typename T, typename Syn>
void queue_set<T,Syn>::remove_item(const T& t) {
   /* note: const_cast is safe (and required) below */
   typename queue_set_q::queue_item q_item(const_cast<T&>(t));
   if (_s.contains(q_item)) {
      typename queue_set_q::queue_item& q_item_match = _s.find(q_item);
      typename queue_set_q::queue_node* q_node = q_item_match.q_node;
      _s.remove(q_item_match);
      _q.remove_heap_node(q_node);
   }
}

/*
 * Remove element(s) from the queue.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::remove(const T& t) {
   auto_write_lock<const Syn> wlock(*this);
   this->remove_item(t);
   return *this;
}

template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::remove(const collection<T>& c) {
   list<T> lst(c);
   auto_write_lock<const Syn> wlock(*this);
   for (typename list<T>::iterator_t i(lst); i.has_next(); )
      this->remove_item(i.next());
   return *this;
}

/*
 * Dequeue.
 * Throw an exception (ex_not_found) if the queue is empty.
 */
template <typename T, typename Syn>
T& queue_set<T,Syn>::dequeue() {
   auto_write_lock<const Syn> wlock(*this);
   if (_q._size == 0)
      throw ex_not_found(
         "attempt to dequeue from empty queue_set"
      );
   T& t = _q._root->q_item->t;
   _s.remove(*(_q._root->q_item));
   _q.remove_heap_node(_q._root);
   return t;
}
     
/*
 * Return top element of queue without removing it.
 * Throw an exception (ex_not_found) if the queue is empty.
 */
template <typename T, typename Syn>
T& queue_set<T,Syn>::head() const {
   auto_read_lock<const Syn> rlock(*this);
   return _q.head();
}

/*
 * Remove all element(s) from the queue.
 * Return a reference to the queue.
 */
template <typename T, typename Syn>
queue_set<T,Syn>& queue_set<T,Syn>::clear() {
   auto_write_lock<const Syn> wlock(*this);
   _s.clear();
   _q.clear();
   return *this;
}

/*
 * Search.
 * Throw an exception (ex_not_found) when attempting to find an element not 
 * contained in the queue.
 */
template <typename T, typename Syn>
bool queue_set<T,Syn>::contains(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   /* note: const_cast is safe (and required) below */
   typename queue_set_q::queue_item q_item(const_cast<T&>(t));
   return _s.contains(q_item);
}

template <typename T, typename Syn>
T& queue_set<T,Syn>::find(const T& t) const {
   auto_read_lock<const Syn> rlock(*this);
   /* note: const_cast is safe (and required) below */
   typename queue_set_q::queue_item q_item(const_cast<T&>(t));
   typename queue_set_q::queue_item& q_item_match = _s.find(q_item);
   return q_item_match.t;
}

/*
 * Size.
 */
template <typename T, typename Syn>
unsigned long queue_set<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _q._size;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > queue_set<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(new queue_set_iterator<T,Syn>(*this));
}

} /* namespace collections */

#endif
