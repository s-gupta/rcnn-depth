/*
 * Iterators.
 */
#ifndef LANG__ITERATORS__ITERATOR_HH
#define LANG__ITERATORS__ITERATOR_HH

namespace lang {
namespace iterators {

/*
 * Abstract base class for iterators.
 */
template <typename T>
class iterator {
public:
   /*
    * Destructor.
    */
   virtual ~iterator() = 0;

   /*
    * Check if there is another item available.
    */
   virtual bool has_next() const = 0;

   /*
    * Return the next item.
    * Throw an exception (ex_not_found) if there are no more items.
    */
   virtual T& next() = 0;
};

/*
 * Iterator pure virtual destructor.
 */
template <typename T>
iterator<T>::~iterator() { }

} /* namespace iterators */
} /* namespace lang */

#endif
