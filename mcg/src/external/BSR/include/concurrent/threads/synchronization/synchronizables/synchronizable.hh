/*
 * Synchronizable interface.
 * 
 * This abstract base class provides an interface for controlling read and 
 * write access to objects by different threads.  A thread requesting a 
 * particular type of access is suspended until that request is granted.
 */
#ifndef CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__SYNCHRONIZABLE_HH
#define CONCURRENT__THREADS__SYNCHRONIZATION__SYNCHRONIZABLES__SYNCHRONIZABLE_HH

namespace concurrent {
namespace threads {
namespace synchronization {
namespace synchronizables {

/*
 * Abstract base class for synchronization capabilities.
 */
class synchronizable {
public:
   /*
    * Constructor.
    */
   synchronizable();
   
   /*
    * Destructor.
    */
   virtual ~synchronizable();

   /*
    * Claim/release exclusive access (abstract interface).
    */
   virtual void abstract_lock() const = 0;
   virtual void abstract_unlock() const = 0;

   /*
    * Claim/release read access (abstract interface).
    */
   virtual void abstract_read_lock() const = 0;
   virtual void abstract_read_unlock() const = 0;

   /*
    * Claim/release write access (abstract interface).
    */
   virtual void abstract_write_lock() const = 0;
   virtual void abstract_write_unlock() const = 0;

   /*
    * Claim/release exclusive access to two synchronizable objects.
    */
   static void lock(const synchronizable&, const synchronizable&);
   static void unlock(const synchronizable&, const synchronizable&);

   /*
    * Claim/release read access to two synchronizable objects.
    */
   static void read_lock(const synchronizable&, const synchronizable&);
   static void read_unlock(const synchronizable&, const synchronizable&);

   /*
    * Claim/release write access to two synchronizable objects.
    */
   static void write_lock(const synchronizable&, const synchronizable&);
   static void write_unlock(const synchronizable&, const synchronizable&);

   /*
    * Claim/release read access to the first object and write access to the 
    * second.
    */
   static void read_write_lock(const synchronizable&, const synchronizable&);
   static void read_write_unlock(const synchronizable&, const synchronizable&);

protected:
   /*
    * Constructor.
    * Initialize the synchronizable object to have the given id.
    */
   explicit synchronizable(unsigned long);

   /*
    * Get identity of synchronizable object.
    */
   inline unsigned long syn_id() const;

private:
   /*
    * Identity of synchronizable object.
    */
   unsigned long _id;

   /*
    * Issue the next available identity.
    * The zero id is reserved for unsynchronized objects.
    */
   static unsigned long id_next();
};

/*
 * Get identity of synchronizable object.
 */
inline unsigned long synchronizable::syn_id() const {
   return _id;
}

} /* namespace synchronizables */
} /* namespace synchronization */
} /* namespace threads */
} /* namespace concurrent */

#endif 
