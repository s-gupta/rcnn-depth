/*
 * Abstract queue.
 */
#ifndef COLLECTIONS__ABSTRACT__QUEUE_HH
#define COLLECTIONS__ABSTRACT__QUEUE_HH

#include "collections/abstract/collection.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
namespace abstract {
/*
 * Imports.
 */
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Abstract base class for queues.
 */
template <typename T>
class queue : virtual public collection<T> {
public:
   /*
    * Destructor.
    */
   virtual ~queue() = 0;

   /*
    * Collection interface for adding element(s) to queue.
    * Enqueue element(s).
    * Return a reference to the queue.
    */
   virtual queue<T>& add(T&);
   virtual queue<T>& add(const collection<T>&);

   /*
    * Enqueue.
    * Return a reference to the queue.
    */
   virtual queue<T>& enqueue(T&) = 0;
   virtual queue<T>& enqueue(const collection<T>&) = 0;

   /*
    * Dequeue.
    * Throw an exception (ex_not_found) if the queue is empty.
    */
   virtual T& dequeue() = 0;

   /*
    * Return top element of queue without removing it.
    * Throw an exception (ex_not_found) if the queue is empty.
    */
   virtual T& head() const = 0;

   /*
    * Remove all element(s) from the queue.
    * Return a reference to the queue.
    */
   virtual queue<T>& clear() = 0;

   /*
    * Size.
    */
   virtual unsigned long size() const = 0;
   
   /*
    * Return iterator over elements.
    */
   virtual auto_ptr< iterator<T> > iter_create() const = 0;
};

/*
 * Pure virtual destructor.
 */
template <typename T>
queue<T>::~queue() { }

/*
 * Collection interface to add an element to the queue.
 * Enqueue the element.
 * Return a reference to the queue.
 */
template <typename T>
queue<T>& queue<T>::add(T& t) {
   return this->enqueue(t);
}

/*
 * Collection interface to add elements to the queue.
 * Enqueue the elements.
 * Return a reference to the queue.
 */
template <typename T>
queue<T>& queue<T>::add(const collection<T>& c) {
   return this->enqueue(c);
}

} /* namespace abstract */
} /* namespace collections */

#endif
