/*
 * Linked lists (thread-safe).
 */
#ifndef COLLECTIONS__LIST_HH
#define COLLECTIONS__LIST_HH

#include "collections/abstract/collection.hh"
#include "collections/abstract/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "functors/comparable_functors.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/serialization/serializers.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"

namespace collections {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::pointers::auto_collection;
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using functors::comparable_functor;
using functors::compare_functors;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::serialization::serializer;
using io::serialization::serializers;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;

/*
 * Declare classes for iterators over linked lists.
 */
template <typename T, typename Syn>
class list_iterator;

template <typename T, typename Syn>
class list_iterator_reverse;

/*
 * Linked lists.
 */ 
template <typename T, typename Syn = unsynchronized>
class list : public abstract::list<T>, 
             protected Syn {
public:
   /*
    * Friend classes.
    */
   friend class list_iterator<T,Syn>;
   friend class list_iterator_reverse<T,Syn>;
   
   /*
    * Define the iterator types.
    */
   typedef list_iterator<T,Syn> iterator_t;
   typedef list_iterator_reverse<T,Syn> iterator_reverse_t;

   /*
    * Constructor.
    */
   list();

   /*
    * Constructor.
    * Create a list from a collection.
    */
   explicit list(const collection<T>&);

   /*
    * Copy constructor.
    */
   list(const list<T,Syn>&);
   
   /*
    * Destructor.
    */
   virtual ~list();

   /*
    * Serialize.
    */
   void serialize(
      serial_output_stream&,
      const serializer<T>& = serializers<T>::s_default()
   ) const;

   /*
    * Deserialize.
    */
   static auto_collection< T, list<T,Syn> > deserialize(
      serial_input_stream&,
      const serializer<T>& = serializers<T>::s_default()
   );

   /*
    * Return head element(s).
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   T& head() const;
   void head(unsigned long, collection<T>&) const;

   /*
    * Return tail element(s). 
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   T& tail() const;
   void tail(unsigned long, collection<T>&) const;

   /*
    * Collection interface for adding element(s) to tail of list.
    * Return a reference to the list.
    */
   list<T,Syn>& add(T&);
   list<T,Syn>& add(const collection<T>&);
    
   /*
    * Addition of element(s) to head of list.
    * Return a reference to the list.
    */
   list<T,Syn>& prepend(T&);
   list<T,Syn>& prepend(const collection<T>&);
   
   /*
    * Addition of element(s) to tail of list.
    * Return a reference to the list.
    */
   list<T,Syn>& append(T&);
   list<T,Syn>& append(const collection<T>&);

   /*
    * Remove and return head element(s).
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   T& remove_head();
   void remove_head(unsigned long, collection<T>&);

   /*
    * Remove and return tail element(s).
    * Throw an exception (ex_not_found) if the list doesn't contain enough
    * elements.
    */
   T& remove_tail();
   void remove_tail(unsigned long, collection<T>&);
   
   /*
    * Remove all element(s) from the list.
    * Return a reference to the list.
    */
   list<T,Syn>& clear();

   /*
    * Reverse list.
    * Return a reference to the list.
    */
   list<T,Syn>& reverse();
   
   /*
    * Size.
    */
   unsigned long size() const;

   /*
    * Return iterators over elements.
    */
   auto_ptr< iterator<T> > iter_create() const;
   auto_ptr< iterator<T> > iter_reverse_create() const;
 
   /*
    * Sort elements in ascending order according to the given comparison
    * functor.
    *
    * Sorting uses an O(n + (n/p)*log(n/p)) time algorithm, where
    * n is the list size and p is the number of available processors.
    */
   void sort(const comparable_functor<T>& = compare_functors<T>::f_compare());
   
   /*
    * Remove duplicate elements (as defined by the given comparison functor)
    * and place the remaining unique elements in sorted order.
    */
   void unique(const comparable_functor<T>& = compare_functors<T>::f_compare());

protected:
   /************************************************************************
    * List data structures.
    ************************************************************************/
   
   /*
    * Node in the list.
    */
   class list_node {
   public:
      /*
       * Constructors.
       */
      explicit list_node(T& item)
       : t(item),
         prev(NULL),
         next(NULL)
      { }
      
      explicit list_node(const list_node& node)
       : t(node.t),
         prev(node.prev),
         next(node.next)
      { }

      /*
       * Destructor.
       */
      ~list_node() { /* do nothing */ }

      /*
       * Node data.
       */
      T& t;                /* item stored in the node */
      list_node* prev;     /* previous node in list */
      list_node* next;     /* next node in list */
   };

   /*
    * List data structure.
    */
   unsigned long _size;    /* number of elements in the list */
   list_node* _head;       /* beginning of list */
   list_node* _tail;       /* end of list */

   /************************************************************************
    * List helper functions.
    ************************************************************************/

   /*
    * Swap the contents of two lists.
    * Note: Lists are not locked by this function.
    */
   static void swap(list<T,Syn>&, list<T,Syn>&);

   /*
    * Runnable object for merge sorting a list in parallel given its head node,
    * tail node, size, and a comparison functor.
    *
    * This is an O(n + (n/p)*log(n/p)) time algorithm, where n is the list size
    * and p is the number of available processors.
    */
   class sorter : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit sorter(
         list_node*&,                           /* head node */
         list_node*&,                           /* tail node */
         unsigned long,                         /* size */
         const comparable_functor<T>&           /* comparison functor */
      );
      
      /*
       * Destructor.
       */
      virtual ~sorter();
      
      /*
       * Run the sort.
       */
      virtual void run();
      
   protected:
      list_node*& _head;                        /* head node */
      list_node*& _tail;                        /* tail node */
      unsigned long _size;                      /* size */
      const comparable_functor<T>& _f_compare;  /* comparison functor */
   };
};

/*
 * List iterator.
 */
template <typename T, typename Syn = unsynchronized>
class list_iterator : public iterator<T> {
public:
   /*
    * Constructor.
    * Seek to the head of the list.
    */
   explicit list_iterator(const list<T,Syn>&);

   /*
    * Copy constructor.
    */
   list_iterator(const list_iterator<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~list_iterator();

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
   const list<T,Syn>& _lst;                  /* list being iterated over */
   typename list<T,Syn>::list_node* _curr;   /* current position in list */
};

/*
 * List reverse iterator.
 * Behavior is what the list_iterator class would do on the reverse list.
 */
template <typename T, typename Syn = unsynchronized>
class list_iterator_reverse : public iterator<T> {
public:
   /*
    * Constructor.
    */
   explicit list_iterator_reverse(const list<T,Syn>&);

   /*
    * Copy constructor.
    */
   list_iterator_reverse(const list_iterator_reverse<T,Syn>&);

   /*
    * Destructor.
    */
   virtual ~list_iterator_reverse();
   
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
   const list<T,Syn>& _lst;                  /* list being iterated over */
   typename list<T,Syn>::list_node* _curr;   /* current position in list */
};

/***************************************************************************
 * List iterator implementation.
 ***************************************************************************/

/*
 * Constructors.
 * Lock list and initialize position.
 */
template <typename T, typename Syn>
list_iterator<T,Syn>::list_iterator(const list<T,Syn>& l)
 : _lst(l)
{
   _lst.read_lock();
   _curr = _lst._head;
}

template <typename T, typename Syn>
list_iterator_reverse<T,Syn>::list_iterator_reverse(const list<T,Syn>& l)
 : _lst(l)
{
   _lst.read_lock();
   _curr = _lst._tail;
}

/*
 * Copy constructors.
 * Lock list again and copy position.
 */
template <typename T, typename Syn>
list_iterator<T,Syn>::list_iterator(
   const list_iterator<T,Syn>& i)
 : _lst(i._lst),
   _curr(i._curr)
{
   _lst.read_lock();
}

template <typename T, typename Syn>
list_iterator_reverse<T,Syn>::list_iterator_reverse(
   const list_iterator_reverse<T,Syn>& i)
 : _lst(i._lst),
   _curr(i._curr)
{
   _lst.read_lock();
}

/*
 * Destructors.
 * Unlock the list.
 */
template <typename T, typename Syn>
list_iterator<T,Syn>::~list_iterator() {
   _lst.read_unlock();
}

template <typename T, typename Syn>
list_iterator_reverse<T,Syn>::~list_iterator_reverse() {
   _lst.read_unlock();
}

/*
 * Check if there is another item available.
 */
template <typename T, typename Syn>
bool list_iterator<T,Syn>::has_next() const {
   return (_curr != NULL);
}

template <typename T, typename Syn>
bool list_iterator_reverse<T,Syn>::has_next() const {
   return (_curr != NULL);
}

/*
 * Return the next item.
 * Throw an exception (ex_not_found) if there are no more items.
 */
template <typename T, typename Syn>
T& list_iterator<T,Syn>::next() {
   if (_curr != NULL) {
      T& t = _curr->t;
      _curr = _curr->next;
      return t;
   } else {
      throw ex_not_found(
         "iterator attempted to read past tail of list"
      );
   }
}

template <typename T, typename Syn>
T& list_iterator_reverse<T,Syn>::next() {
   if (_curr != NULL) {
      T& t = _curr->t;
      _curr = _curr->prev;
      return t;
   } else {
      throw ex_not_found(
         "reverse iterator attempted to read past head of list"
      );
   }
}

/***************************************************************************
 * List implementation.
 ***************************************************************************/

/*
 * Swap the contents of two lists.
 * Note: Lists are not locked by this function.
 */
template <typename T, typename Syn>
void list<T,Syn>::swap(list<T,Syn>& l0, list<T,Syn>& l1) {
   unsigned long temp_size = l0._size;
   list_node*    temp_head = l0._head;
   list_node*    temp_tail = l0._tail;
   l0._size = l1._size;
   l0._head = l1._head;
   l0._tail = l1._tail;
   l1._size = temp_size;
   l1._head = temp_head;
   l1._tail = temp_tail;
}

/*
 * Default constructor.
 * Create an empty list.
 */
template <typename T, typename Syn>
list<T,Syn>::list()
 : Syn(),
   _size(0),
   _head(NULL),
   _tail(NULL)
{ }

/*
 * Constructor.
 * Create a list from a collection.
 */
template <typename T, typename Syn>
list<T,Syn>::list(const collection<T>& c)
 : Syn(),
   _size(0),
   _head(NULL),
   _tail(NULL)
{
   /* add elements to empty list */
   list<T,Syn> l;
   auto_ptr< iterator<T> > i = c.iter_create();
   if (i->has_next()) {
      l._head = new list_node(i->next());
      l._tail = l._head;
      l._size++;
   }
   while (i->has_next()) {
      list_node* node = new list_node(i->next());
      node->prev = l._tail;
      l._tail->next = node;
      l._tail = node;
      l._size++;
   }
   /* take ownership of list */
   list<T,Syn>::swap(*this, l);
}

/*
 * Copy constructor.
 */
template <typename T, typename Syn>
list<T,Syn>::list(const list<T,Syn>& l)
 : Syn(),
   _size(0),
   _head(NULL),
   _tail(NULL)
{
   /* copy source list */
   auto_read_lock<const Syn> rlock(l);
   list<T,Syn> l_copy;
   list_node* curr = l._head;
   if (curr != NULL) {
      l_copy._head = new list_node(curr->t);
      l_copy._tail = l_copy._head;
      l_copy._size++;
      curr = curr->next;
   }
   while (curr != NULL) {
      list_node* node = new list_node(curr->t);
      node->prev = l_copy._tail;
      l_copy._tail->next = node;
      l_copy._tail = node;
      l_copy._size++;
      curr = curr->next;
   }
   /* take ownership of copy */
   list<T,Syn>::swap(*this, l_copy);
}

/*
 * Destructor.
 * Delete the list nodes.
 */
template <typename T, typename Syn>
list<T,Syn>::~list() {
   list_node* curr = _head;
   while (curr != NULL) {
      list_node* temp = curr;
      curr = curr->next;
      delete temp;
   }
}

/*
 * Serialize.
 */
template <typename T, typename Syn>
void list<T,Syn>::serialize(
   serial_output_stream& s, const serializer<T>& slzr) const
{
   auto_read_lock<const Syn> rlock(*this);
   s << _size;
   list_node* curr = _head;
   while (curr != NULL) {
      slzr.serialize(s, curr->t);
      curr = curr->next;
   }
}

/*
 * Deserialize.
 */
template <typename T, typename Syn>
auto_collection< T, list<T,Syn> > list<T,Syn>::deserialize(
   serial_input_stream& s, const serializer<T>& slzr)
{
   auto_collection< T, list<T,Syn> > l(new list<T,Syn>());
   unsigned long size = 0;
   s >> size;
   for (unsigned long n = 0; n < size; n++) {
      auto_ptr<T> t = slzr.deserialize(s);
      l->add(*t);
      t.release();
   }
   return l;
}

/*
 * Return head element.
 * Throw an exception (ex_not_found) if the list is empty.
 */
template <typename T, typename Syn>
T& list<T,Syn>::head() const {
   auto_read_lock<const Syn> rlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to fetch head element of empty list"
      );
   return _head->t;
}

/*
 * Add the first n elements of the list to the given collection.
 * Throw an exception (ex_not_found) if the list doesn't contain at least
 * n elements.
 */
template <typename T, typename Syn>
void list<T,Syn>::head(unsigned long n, collection<T>& c) const {
   list<T> l;
   {
      /* lock list and obtain first n elements */
      auto_read_lock<const Syn> rlock(*this);
      if (_size < n)
         throw ex_not_found(
            "attempt to fetch too many head elements of list"
         );
      list_node* curr = _head;
      for (unsigned long i = 0; i < n; i++) {
         l.append(curr->t);
         curr = curr->next;
      }
   }
   c.add(l);
}

/*
 * Return tail element.
 * Throw an exception (ex_not_found) if the list is empty.
 */
template <typename T, typename Syn>
T& list<T,Syn>::tail() const {
   auto_read_lock<const Syn> rlock(*this);
   if (_size == 0)
      throw ex_not_found(
         "attempt to fetch tail element of empty list"
      );
   return _tail->t;
}

/*
 * Add the last n elements of the list to the given collection.
 * Throw an exception (ex_not_found) if the list doesn't contain at least
 * n elements.
 */
template <typename T, typename Syn>
void list<T,Syn>::tail(unsigned long n, collection<T>& c) const {
   list<T> l;
   {
      /* lock list and obtain first n elements */
      auto_read_lock<const Syn> rlock(*this);
      if (_size < n)
         throw ex_not_found(
            "attempt to fetch too many tail elements of list"
         );
      list_node* curr = _tail;
      for (unsigned long i = 0; i < n; i++) {
         l.prepend(curr->t);
         curr = curr->prev;
      }
   }
   c.add(l);
}

/*
 * Collection interface to add an element to the tail of the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::add(T& t) {
   return this->append(t);
}

/*
 * Add all elements in a collection to the tail of the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::add(const collection<T>& c) {
   return this->append(c);
}

/*
 * Add an element to the head of the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::prepend(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   list_node* curr = new list_node(t);
   curr->next = _head;
   if (_head != NULL)
      _head->prev = curr;
   else
      _tail = curr;
   _head = curr;
   _size++;
   return *this;
}

/*
 * Add multiple elements to the head of the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::prepend(const collection<T>& c) {
   /* create list nodes for collection elements */
   list<T,Syn> l(c);
   /* link nodes */
   auto_write_lock<const Syn> wlock(*this);
   if (_size == 0) {
      list<T,Syn>::swap(*this, l);
   } else if (l._size != 0) {
      l._tail->next = _head;
      _head->prev = l._tail;
      _head = l._head;
      _size += l._size;
      l._size = 0;
      l._head = NULL;
      l._tail = NULL;
   }
   return *this;
}

/*
 * Add an element to the tail of the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::append(T& t) {
   auto_write_lock<const Syn> wlock(*this);
   list_node* curr = new list_node(t);
   curr->prev = _tail;
   if (_tail != NULL) 
      _tail->next = curr;
   else
      _head = curr;
   _tail = curr;
   _size++;
   return *this;
}

/*
 * Add multiple elements to the tail of the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::append(const collection<T>& c) {
   /* create list nodes for collection elements */
   list<T,Syn> l(c);
   /* link nodes */
   auto_write_lock<const Syn> wlock(*this);
   if (_size == 0) {
      list<T,Syn>::swap(*this, l);
   } else if (l._size != 0) {
      l._head->prev = _tail;
      _tail->next = l._head;
      _tail = l._tail;
      _size += l._size;
      l._size = 0;
      l._head = NULL;
      l._tail = NULL;
   }
   return *this;
}

/*
 * Remove and return the head element of the list.
 * Throw an exception (ex_not_found) if the list is empty.
 */
template <typename T, typename Syn>
T& list<T,Syn>::remove_head() {
   auto_write_lock<const Syn> wlock(*this);
   /* check size */
   if (_size == 0)
      throw ex_not_found(
         "attempt to remove head element of empty list"
      );
   /* remove head element */
   T& t = _head->t;
   list_node* temp = _head;
   _head = _head->next;
   if (_head != NULL)
      _head->prev = NULL;
   else
      _tail = NULL;
   _size--;
   delete temp;
   return t;
}

/*
 * Remove and return the first n elements of the list.
 * Throw an exception (ex_not_found) if the list contains less than n elements.
 */
template <typename T, typename Syn>
void list<T,Syn>::remove_head(unsigned long n, collection<T>& c) {
   list<T,Syn> l;
   {
      /* lock list, obtain and remove first n elements */
      auto_write_lock<const Syn> wlock(*this);
      if (_size < n) {
         throw ex_not_found(
            "attempt to remove too many head elements of list"
         );
      } else if (_size == n) {
         list<T,Syn>::swap(*this, l);
      } else if (n > 0) {
         l._size = n;
         l._head = _head;
         for (unsigned long i = 0; i < n; i++)
            _head = _head->next;
         l._tail = _head->prev;
         l._tail->next = NULL;
         _head->prev = NULL;
         _size -= n;
      }
   }
   c.add(l);
}

/*
 * Remove and return the tail element of the list.
 * Throw an exception (ex_not_found) if the list is empty.
 */
template <typename T, typename Syn>
T& list<T,Syn>::remove_tail() {
   auto_write_lock<const Syn> wlock(*this);
   /* check size */
   if (_size == 0)
      throw ex_not_found(
         "attempt to remove tail element of empty list"
      );
   /* remove tail element */
   T& t = _tail->t;
   list_node* temp = _tail;
   _tail = _tail->prev;
   if (_tail != NULL)
      _tail->next = NULL;
   else
      _head = NULL;
   _size--;
   delete temp;
   return t;
}

/*
 * Remove and return the last n elements of the list.
 * Throw an exception (ex_not_found) if the list contains less than n elements.
 */
template <typename T, typename Syn>
void list<T,Syn>::remove_tail(unsigned long n, collection<T>& c) {
   list<T,Syn> l;
   {
      /* lock list, obtain and remove last n elements */
      auto_write_lock<const Syn> wlock(*this);
      if (_size < n) {
         throw ex_not_found(
            "attempt to remove too many tail elements of list"
         );
      } else if (_size == n) {
         list<T,Syn>::swap(*this, l);
      } else if (n > 0) {
         l._size = n;
         l._tail = _tail;
         for (unsigned long i = 0; i < n; i++)
            _tail = _tail->prev;
         l._head = _tail->next;
         l._head->prev = NULL;
         _tail->next = NULL;
         _size -= n;
      }
   }
   c.add(l);
}

/*
 * Remove all element(s) from the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::clear() {
   list<T,Syn> l;
   auto_write_lock<const Syn> wlock(*this);
   list<T,Syn>::swap(*this, l);
   return *this;
}

/*
 * Reverse the order of elements in the list.
 * Return a reference to the list.
 */
template <typename T, typename Syn>
list<T,Syn>& list<T,Syn>::reverse() {
   auto_write_lock<const Syn> wlock(*this);
   list_node* curr = _head;
   while (curr != NULL) {
      list_node* temp = curr->next;
      curr->next = curr->prev;
      curr->prev = temp;
      curr = temp;
   }
   list_node* temp = _tail;
   _tail = _head;
   _head = temp;
   return *this;
}

/*
 * Get number of elements in list.
 */
template <typename T, typename Syn>
unsigned long list<T,Syn>::size() const {
   auto_read_lock<const Syn> rlock(*this);
   return _size;
}

/*
 * Return iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > list<T,Syn>::iter_create() const {
   return auto_ptr< iterator<T> >(new list_iterator<T,Syn>(*this));
}

/*
 * Return reverse iterator over elements.
 */
template <typename T, typename Syn>
auto_ptr< iterator<T> > list<T,Syn>::iter_reverse_create() const {
   return auto_ptr< iterator<T> >(new list_iterator_reverse<T,Syn>(*this));
}

/*
 * Sort elements in ascending order according to the given comparison functor.
 *
 * Sorting uses the O(n + (n/p)*log(n/p)) time merge sort algorithm, where
 * n is the list size and p is the number of available processors.
 */
template <typename T, typename Syn>
void list<T,Syn>::sort(const comparable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   sorter s(_head, _tail, _size, f);
   s.run();
}

/*
 * Remove duplicate elements (as defined by the given comparison functor)
 * and place the remaining unique elements in sorted order.
 */
template <typename T, typename Syn>
void list<T,Syn>::unique(const comparable_functor<T>& f) {
   auto_write_lock<const Syn> wlock(*this);
   /* sort the list */
   sorter s(_head, _tail, _size, f);
   s.run();
   /* remove duplicate elements */
   list_node* curr = _head;
   list_node* temp = (_head == NULL) ? NULL : (_head->next);
   while (temp != NULL) {
      /* check if duplicate */
      if (f(temp->t, curr->t) == 0) {
         /* remove duplicate */
         list_node* temp_next = temp->next;
         curr->next = temp_next;
         if (temp_next != NULL)
            temp_next->prev = curr;
         delete temp;
      } else {
         /* move to next element in list */
         curr = temp;
      }
      temp = curr->next;
   }
   /* set list tail */
   _tail = curr;
}

/***************************************************************************
 * List sorter implementation.
 ***************************************************************************/

/*
 * Constructor.
 */
template <typename T, typename Syn>
list<T,Syn>::sorter::sorter(
   list_node*& head, 
   list_node*& tail, 
   unsigned long size, 
   const comparable_functor<T>& f)
 : _head(head),
   _tail(tail),
   _size(size),
   _f_compare(f)
{ }

/*
 * Destructor.
 */
template <typename T, typename Syn>
list<T,Syn>::sorter::~sorter() {
   /* do nothing */
}

/*
 * Run the sort.
 */
template <typename T, typename Syn>
void list<T,Syn>::sorter::run() {
   /* recursively sort each half of the list and merge the result */
   if (_size > 1) {
      /* split list - compute size of each half */
      unsigned long size0 = _size/2;
      unsigned long size1 = _size - size0;
      /* split list - locate middle nodes */
      list_node* mid0 = _head;
      for (unsigned long n = 1; n < size0; n++)
         mid0 = mid0->next;
      list_node* mid1 = mid0->next;
      /* split list - break middle node links */
      mid0->next = NULL;
      mid1->prev = NULL;
      /* sort each half */
      sorter s0(_head,  mid0, size0, _f_compare);
      sorter s1( mid1, _tail, size1, _f_compare);
      child_thread::run(s0, s1);
      /* initialize head of sorted list */
      if (_f_compare(mid1->t, _head->t) < 0) {
         list_node* temp = _head;
         _head = mid1;
         mid1 = temp;
      }
      list_node* curr      = _head;
      list_node* curr_next = _head->next;
      list_node* other     = mid1;
      /* merge sorted halves */
      while (curr_next != NULL) {
         /* compare next items from each half */
         if (_f_compare(other->t, curr_next->t) < 0) {
            /* switch halves */
            curr->next  = other;
            other->prev = curr;
            other = curr_next;
         }
         /* continue along the current half */
         curr      = curr->next;
         curr_next = curr->next;
      }
      /* link rest of other list to end of current list */
      curr->next  = other;
      other->prev = curr;
      /* determine tail of sorted list */
      if (_tail->next != NULL)
         _tail = mid0;
   }
}

} /* namespace collections */

#endif
