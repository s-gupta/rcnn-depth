/*
 * Thread-safe iterators.
 */
#ifndef LANG__ITERATORS__ITERATOR_SYN_HH
#define LANG__ITERATORS__ITERATOR_SYN_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "lang/iterators/iterator.hh"

namespace lang {
namespace iterators {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;

/*
 * A synchronized iterator wraps synchronization locks around another iterator.
 */
template <typename T, typename Syn = unsynchronized>
class iterator_syn : public iterator<T>,
                     protected Syn {
public:
   /*
    * Constructor.
    */
   iterator_syn(iterator<T>&);

   /*
    * Copy constructor.
    */
   explicit iterator_syn(const iterator_syn<T>&);

    /*
    * Destructor.
    */
   virtual ~iterator_syn();

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
   iterator<T>& _iter;   /* contained iterator */
};

/*
 * Constructor.
 */
template <typename T, typename Syn>
iterator_syn<T,Syn>::iterator_syn(iterator<T>& i)
 : Syn(),
   _iter(i)
{ }

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
iterator_syn<T,Syn>::iterator_syn(const iterator_syn<T>& i)
 : Syn(), 
   _iter(i._iter)
{ }

 /*
 * Destructor.
 */
template <typename T, typename Syn>
iterator_syn<T,Syn>::~iterator_syn() {
   /* do nothing */
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool iterator_syn<T,Syn>::has_next() const {
   auto_read_lock<const Syn> rlock(*this);
   return _iter.has_next();
}

/*
 * Return the next item.
 * Throw an exception (not found) if there are no more items.
 */
template <typename T, typename Syn>
T& iterator_syn<T,Syn>::next() {
   auto_write_lock<const Syn> wlock(*this);
   return _iter.next();
}

} /* namespace iterators */
} /* namespace lang */

#endif
