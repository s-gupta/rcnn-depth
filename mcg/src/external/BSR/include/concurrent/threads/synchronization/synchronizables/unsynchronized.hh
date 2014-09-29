/*
 * The unsynchronized class is a dummy synchronizable interface for use in
 * cases where actual synchronization is not needed.
 *
 * We provide versions of some of the synchronizable static member functions
 * specialized to the unsynchronized class in order to increase speed when
 * unsynchronized objects are used.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__UNSYNCHRONIZED_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__UNSYNCHRONIZED_HH

#include "concurrent/threads/synchronization/synchronizables/synchronizable.hh"

namespace concurrent {
namespace threads {
namespace synchronization {
namespace synchronizables {

class unsynchronized : public synchronizable {
public:
   /*
    * Constructor.
    */
   unsynchronized();
   
   /*
    * Destructor.
    */
   virtual ~unsynchronized();

   /*
    * Claim exclusive access.
    */
   virtual void abstract_lock() const;
   inline void lock() const;
   
   /*
    * Release exclusive access.
    */
   virtual void abstract_unlock() const; 
   inline void unlock() const;

   /*
    * Claim read access.
    */
   virtual void abstract_read_lock() const;
   inline void read_lock() const;

   /*
    * Release read access.
    */
   virtual void abstract_read_unlock() const;
   inline void read_unlock() const;

   /*
    * Claim write access.
    */
   virtual void abstract_write_lock() const;
   inline void write_lock() const;

   /*
    * Release write access.
    */
   virtual void abstract_write_unlock() const;
   inline void write_unlock() const;

   /*
    * Claim/release exclusive access to two unsynchronized objects.
    */
   inline static void lock(
      const unsynchronized&,
      const unsynchronized&
   );
   
   inline static void unlock(
      const unsynchronized&,
      const unsynchronized&
   );

   /*
    * Claim/release read access to two unsynchronized objects.
    */
   inline static void read_lock(
      const unsynchronized&, 
      const unsynchronized&
   );
   
   inline static void read_unlock(
      const unsynchronized&,
      const unsynchronized&
   );

   /*
    * Claim/release write access to two unsynchronized objects.
    */
   inline static void write_lock(
      const unsynchronized&,
      const unsynchronized&
   );
   
   inline static void write_unlock(
      const unsynchronized&,
      const unsynchronized&
   );

   /*
    * Claim/release read access to the first object and write access to the
    * second.
    */
   inline static void read_write_lock(
      const unsynchronized&,
      const unsynchronized&
   );
   
   inline static void read_write_unlock(
      const unsynchronized&,
      const unsynchronized&
   );

private:
   /*
    * Private copy constructor.
    * (unsynchronized objects should not be copied)
    */
   explicit unsynchronized(const unsynchronized&);
};
   
/*
 * Claim exclusive access (inline version).
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::lock() const { }

/*
 * Release exclusive access (inline version).
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::unlock() const { }

/*
 * Claim read access (inline version).
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::read_lock() const { }

/*
 * Release read access (inline version).
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::read_unlock() const { }

/*
 * Claim write access (inline version).
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::write_lock() const { }

/*
 * Release write access (inline version).
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::write_unlock() const { }

/*
 * Claim/release exclusive access to two unsynchronized objects.
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::lock(
   const unsynchronized&,
   const unsynchronized&)
{ }

inline void unsynchronized::unlock(
   const unsynchronized&,
   const unsynchronized&)
{ }

/*
 * Claim/release read access to two unsynchronized objects.
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::read_lock(
   const unsynchronized&,
   const unsynchronized&)
{ }

inline void unsynchronized::read_unlock(
   const unsynchronized&,
   const unsynchronized&)
{ }

/*
 * Claim/release write access to two unsynchronized objects.
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::write_lock(
   const unsynchronized&,
   const unsynchronized&)
{ }

inline void unsynchronized::write_unlock(
   const unsynchronized&,
   const unsynchronized&)
{ }

/*
 * Claim/release read access to the first object and write access to the second.
 * Do nothing as this is a dummy interface.
 */
inline void unsynchronized::read_write_lock(
   const unsynchronized&,
   const unsynchronized&)
{ }

inline void unsynchronized::read_write_unlock(
   const unsynchronized&,
   const unsynchronized&)
{ }

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif
