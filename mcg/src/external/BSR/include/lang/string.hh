/*
 * Strings (thread-safe).
 *
 * Strings are the same as arrays of characters, except that an additional 
 * hidden '\0' terminator is added to the end of the string for backwards 
 * compatibility with the char* representation of strings.
 */
#ifndef LANG__STRING_HH
#define LANG__STRING_HH

#include "concurrent/threads/synchronization/locks/auto_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_read_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_read_write_lock.hh"
#include "concurrent/threads/synchronization/locks/auto_write_lock.hh"
#include "concurrent/threads/synchronization/synchronizables/unsynchronized.hh"
#include "functors/comparable_functors.hh"
#include "functors/filterable_functors.hh"
#include "interfaces/comparable.hh"
#include "io/serialization/serial_input_stream.hh"
#include "io/serialization/serial_output_stream.hh"
#include "io/streams/ostream.hh"
#include "lang/array.hh"
#include "lang/null.hh"
#include "lang/pointers/auto_ptr.hh"

#include <cstring>

namespace lang {
/*
 * Imports.
 */
using concurrent::threads::synchronization::locks::auto_read_lock;
using concurrent::threads::synchronization::locks::auto_read_read_lock;
using concurrent::threads::synchronization::locks::auto_read_write_lock;
using concurrent::threads::synchronization::locks::auto_write_lock;
using concurrent::threads::synchronization::synchronizables::unsynchronized;
using functors::comparable_functor;
using functors::compare_functors;
using functors::filterable_functor;
using interfaces::comparable;
using io::serialization::serial_input_stream;
using io::serialization::serial_output_stream;
using io::streams::ostream;
using lang::pointers::auto_ptr;

/*
 * Declare prototypes for string template friend functions.
 */
template <typename Syn> class string;
template <typename Syn> ostream& operator<<(ostream&, const string<Syn>&);
template <typename Syn> bool operator==(const string<Syn>&, const char*);
template <typename Syn> bool operator!=(const string<Syn>&, const char*);
template <typename Syn> bool operator< (const string<Syn>&, const char*);
template <typename Syn> bool operator> (const string<Syn>&, const char*);
template <typename Syn> bool operator<=(const string<Syn>&, const char*);
template <typename Syn> bool operator>=(const string<Syn>&, const char*);
template <typename Syn> bool operator==(const char*, const string<Syn>&);
template <typename Syn> bool operator!=(const char*, const string<Syn>&);
template <typename Syn> bool operator< (const char*, const string<Syn>&);
template <typename Syn> bool operator> (const char*, const string<Syn>&);
template <typename Syn> bool operator<=(const char*, const string<Syn>&);
template <typename Syn> bool operator>=(const char*, const string<Syn>&);
template <typename Syn> bool operator==(const string<Syn>&, const string<Syn>&);
template <typename Syn> bool operator!=(const string<Syn>&, const string<Syn>&);
template <typename Syn> bool operator< (const string<Syn>&, const string<Syn>&);
template <typename Syn> bool operator> (const string<Syn>&, const string<Syn>&);
template <typename Syn> bool operator<=(const string<Syn>&, const string<Syn>&);
template <typename Syn> bool operator>=(const string<Syn>&, const string<Syn>&);

/*
 * String class.
 */
template <typename Syn = unsynchronized>
class string : public array<char,Syn>, 
               public comparable< string<Syn> > {
public:
   /*
    * Friend classes.
    */
   template <typename S> friend class string;

   /*
    * Constructors.
    */
   string();
   string(const char*);
   explicit string(unsigned long, char = '\0');

   /*
    * Copy constructors.
    */
   string(const string<Syn>&);
   template <typename S> string(const string<S>&);

   /*
    * Constructor from array.
    */
   template <typename S> explicit string(const array<char,S>&);
    
   /*
    * Destructor.
    */
   virtual ~string();

   /*
    * Assignment operator.
    */
   string<Syn>& operator=(const string<Syn>&);
   template <typename S> string<Syn>& operator=(const string<S>&);

   /*
    * Implicit conversion to char pointer.
    */
   operator const char* () const;

   /*
    * Substring operations.
    */
   string<Syn> substring(
      unsigned long,                /* start index */
      unsigned long                 /* end index */
   ) const;
   
   string<Syn> substring(
      unsigned long,                /* start index */
      unsigned long,                /* step size */
      unsigned long                 /* end index */
   ) const;
   
   string<Syn> substring(
      const array<unsigned long>&   /* indices of desired elements */
   ) const;

   /*
    * Reverse character order.
    * Return a reference to the string.
    */
   string<Syn>& reverse();
   
   /*
    * Resize.
    * Initialize new characters to null terminator (if enlarging).
    * Return a reference to the string.
    */
   string<Syn>& resize(unsigned long);
   
   /*
    * Resize.
    * Initialize new characters to the specified value (if enlarging).
    * Return a reference to the string.
    */
   string<Syn>& resize(unsigned long, const char&);

   /*
    * Conversion to lowercase/uppercase (in place).
    */
   string<Syn>& lowercase();
   string<Syn>& uppercase();

   /*
    * Conversion to lowercase/uppercase (returning new string).
    */
   static string<Syn> lowercase(const string<Syn>&);
   static string<Syn> uppercase(const string<Syn>&);

   /*
    * String concatenation (in place).
    */
   string<Syn>& concat(const string<Syn>&);
   
   /*
    * String concatenation (returning new string).
    */
   static string<Syn> concat(const string<Syn>&, const string<Syn>&);
  
   /*
    * Apply filter to the string.
    */
   string<Syn> filter(const filterable_functor<char>&);
   string<Syn> filter(const filterable_functor<const char>&) const;

   /*
    * Return the set of unique characters (as defined by the given comparison
    * functor) in the string (replacing the original string).  
    */
   void unique(
      const comparable_functor<char>& = compare_functors<char>::f_compare()
   );

   /*
    * Return the set of unique characters (as defined by the given comparison
    * functor) in the string (replacing the original string).  In addition,
    * return an index array containing the original positions of the unique
    * characters.
    */
   array<unsigned long> unique_idx(
      const comparable_functor<char>& = compare_functors<char>::f_compare()
   );

   /*
    * Serialize.
    */
   void serialize(serial_output_stream&) const;

   /*
    * Deserialize.
    */
   static auto_ptr< string<Syn> > deserialize(serial_input_stream&);

   /*
    * Formatted output to stream.
    */
   friend ostream& operator<< <Syn>(ostream&, const string<Syn>&);

   /*
    * Comparators: string-char*.
    */
   friend bool operator== <Syn>(const string<Syn>&, const char*);
   friend bool operator!= <Syn>(const string<Syn>&, const char*);
   friend bool operator<  <Syn>(const string<Syn>&, const char*);
   friend bool operator>  <Syn>(const string<Syn>&, const char*);
   friend bool operator<= <Syn>(const string<Syn>&, const char*);
   friend bool operator>= <Syn>(const string<Syn>&, const char*);

   friend bool operator== <Syn>(const char*, const string<Syn>&);
   friend bool operator!= <Syn>(const char*, const string<Syn>&);
   friend bool operator<  <Syn>(const char*, const string<Syn>&);
   friend bool operator>  <Syn>(const char*, const string<Syn>&);
   friend bool operator<= <Syn>(const char*, const string<Syn>&);
   friend bool operator>= <Syn>(const char*, const string<Syn>&);

   /*
    * Comparators: string-string.
    */
   friend bool operator== <Syn>(const string<Syn>&, const string<Syn>&);
   friend bool operator!= <Syn>(const string<Syn>&, const string<Syn>&);
   friend bool operator<  <Syn>(const string<Syn>&, const string<Syn>&);
   friend bool operator>  <Syn>(const string<Syn>&, const string<Syn>&);
   friend bool operator<= <Syn>(const string<Syn>&, const string<Syn>&);
   friend bool operator>= <Syn>(const string<Syn>&, const string<Syn>&);

   /*
    * Comparison to char*.
    */
   int compare_to(const char*) const;
   
   /*
    * Comparison to string.
    */
   int compare_to(const string<Syn>&) const;
};

/*
 * Default Constructor.
 * Create empty string.
 */
template <typename Syn>
string<Syn>::string()
 : array<char,Syn>(1,'\0') 
{
   /* hide terminator element */
   this->_size = 0;
}

/*
 * Constructor.
 * Create string from character array.
 */
template <typename Syn>
string<Syn>::string(const char* s)
 : array<char,Syn>() 
{
   this->_size = ((s == NULL) ? 0 : std::strlen(s));
   this->_data = new char[this->_size+1]();
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] = s[n];
   this->_data[this->_size] = '\0';
}

/*
 * Constructor.
 * Create a string of specified size and initialize elements.
 */
template <typename Syn>
string<Syn>::string(unsigned long size, char c)
 : array<char,Syn>(size+1, c) 
{
   /* set last element to terminator and hide it */
   this->_size = size;
   this->_data[size] = '\0';
}

/*
 * Copy constructor.
 */
template <typename Syn>
string<Syn>::string(const string<Syn>& s)
 : array<char,Syn>()
{
   /* lock source string */
   auto_read_lock<const Syn> rlock(s);
   /* copy character array */
   this->_size = s._size;
   this->_data = new char[this->_size+1]();
   for (unsigned long n = 0; n <= this->_size; n++)
      this->_data[n] = s._data[n];
}

/*
 * Copy constructor (type conversion).
 */
template <typename Syn>
template <typename S>
string<Syn>::string(const string<S>& s)
 : array<char,Syn>() 
{
   /* lock source string */
   auto_read_lock<const S> rlock(s);
   /* copy character array */
   this->_size = s._size;
   this->_data = new char[this->_size+1]();
   for (unsigned long n = 0; n <= this->_size; n++)
      this->_data[n] = s._data[n];
}

/*
 * Constructor from array.
 * Build a string from an array of characters (add null terminator).
 */
template <typename Syn>
template <typename S>
string<Syn>::string(const array<char,S>& a)
 : array<char,Syn>() 
{
   /* lock source array */
   auto_read_lock<const S> rlock(a);
   /* copy characters */
   this->_size = a._size;
   this->_data = new char[this->_size+1]();
   for (unsigned long n = 0; n < this->_size; n++)
      this->_data[n] = a._data[n];
   this->_data[this->_size] = '\0';
}

/*
 * Destructor.
 */
template <typename Syn>
string<Syn>::~string() {
   /* do nothing (array's destructor takes care of everything) */
}

/*
 * Assignment operator.
 */
template <typename Syn>
string<Syn>& string<Syn>::operator=(const string<Syn>& s) {
   /* copy source string */
   string<Syn> s_copy(s);
   /* swap copy into current string */
   auto_write_lock<const Syn> wlock(*this);
   string<Syn>::swap(*this, s_copy);
   return *this;
}

/*
 * Assignment operator (type conversion).
 */
template <typename Syn>
template <typename S>
string<Syn>& string<Syn>::operator=(const string<S>& s) {
   /* copy source string */
   string<Syn> s_copy(s);
   /* swap copy into current string */
   auto_write_lock<const Syn> wlock(*this);
   string<Syn>::swap(*this, s_copy);
   return *this;
}

/*
 * Implicit conversion to char pointer.
 */
template <typename Syn>
string<Syn>::operator const char* () const {
   auto_read_lock<const Syn> rlock(*this);
   return this->_data;
}

/*
 * Return substring of characters in range [start:end].
 */
template <typename Syn>
string<Syn> string<Syn>::substring(
   unsigned long start,
   unsigned long end) const 
{
   return this->substring(start, 1, end);
}

/*
 * Return substring of characters in range [start:step:end].
 */
template <typename Syn>
string<Syn> string<Syn>::substring(
   unsigned long start,
   unsigned long step,
   unsigned long end) const 
{
   return string<Syn>(this->subarray(start, step, end));
}

/*
 * Return the string of characters taken from the given indices.
 */
template <typename Syn>
string<Syn> string<Syn>::substring(
   const array<unsigned long>& index_arr) const
{
   return string<Syn>(this->subarray(index_arr));
} 

/*
 * Reverse character order.
 * Return a reference to the string.
 */
template <typename Syn>
string<Syn>& string<Syn>::reverse() {
   this->array<char,Syn>::reverse();
   return *this;
}

/*
 * Resize.
 * Initialize new characters to null terminator (if enlarging).
 * Return a reference to the string.
 */
template <typename Syn>
string<Syn>& string<Syn>::resize(unsigned long size) {
   /* allocate resized string */
   string<Syn> s(size);
   /* lock current string */
   auto_write_lock<const Syn> wlock(*this);
   /* copy current string into part of resized string */
   unsigned long size_min = (this->_size < s._size) ? this->_size : s._size;
   for (unsigned long n = 0; n < size_min; n++)
      s._data[n] = this->_data[n];
   /* swap resized string into current */
   string<Syn>::swap(*this, s);
   return *this;
}

/*
 * Resize.
 * Initialize new characters to the specified value (if enlarging).
 * Return a reference to the string.
 */
template <typename Syn>
string<Syn>& string<Syn>::resize(unsigned long size, const char& c) {
   /* allocate resized string */
   string<Syn> s(size, c);
   /* lock current string */
   auto_write_lock<const Syn> wlock(*this);
   /* copy current string into part of resized string */
   unsigned long size_min = (this->_size < s._size) ? this->_size : s._size;
   for (unsigned long n = 0; n < size_min; n++)
      s._data[n] = this->_data[n];
   /* swap resized string into current */
   string<Syn>::swap(*this, s);
   return *this;
}

/*
 * Conversion to lowercase (in place).
 * Return a reference to the string.
 */
template <typename Syn>
string<Syn>& string<Syn>::lowercase() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++) {
      char c = this->_data[n];
      if (('A' <= c) && (c <= 'Z'))
         this->_data[n] = c - 'A' + 'a';
   }
   return *this;
}

/*
 * Conversion to uppercase (in place).
 * Return a reference to the string.
 */
template <typename Syn>
string<Syn>& string<Syn>::uppercase() {
   auto_write_lock<const Syn> wlock(*this);
   for (unsigned long n = 0; n < this->_size; n++) {
      char c = this->_data[n];
      if (('a' <= c) && (c <= 'z'))
         this->_data[n] = c - 'a' + 'A';
   }
   return *this;
}

/*
 * Conversion to lowercase (returning new string).
 */
template <typename Syn>
string<Syn> string<Syn>::lowercase(const string<Syn>& s) {
   string<Syn> s_lower(s);
   s_lower.lowercase();
   return s_lower;
}

/*
 * Conversion to uppercase (returning new string).
 */
template <typename Syn>
string<Syn> string<Syn>::uppercase(const string<Syn>& s) {
   string<Syn> s_upper(s);
   s_upper.uppercase();
   return s_upper;
}

/*
 * String concatenation (in place).
 */
template <typename Syn>
string<Syn>& string<Syn>::concat(const string<Syn>& s) {
   auto_read_write_lock<const Syn> rwlock(s, *this);
   /* form concatenated string */
   string<Syn> s_concat(this->_size + s._size);
   for (unsigned long n = 0; n < this->_size; n++)
      s_concat._data[n] = this->_data[n];
   for (unsigned long n = 0, m = this->_size; n < s._size; n++, m++)
      s_concat._data[m] = s._data[n];
   /* swap concatented string into current */
   string<Syn>::swap(*this, s_concat);
   return *this;
}

/*
 * String concatenation (returning new string).
 */
template <typename Syn>
string<Syn> string<Syn>::concat(const string<Syn>& s0, const string<Syn>& s1) {
   auto_read_read_lock<const Syn> rrlock(s0, s1);
   string<Syn> s_concat(s0._size + s1._size);
   for (unsigned long n = 0; n < s0._size; n++)
      s_concat._data[n] = s0._data[n];
   for (unsigned long n = 0, m = s0._size; n < s1._size; n++, m++)
      s_concat._data[m] = s1._data[n];
   return s_concat;
}

/*
 * Apply filter to non-const string.
 */
template <typename Syn>
string<Syn> string<Syn>::filter(const filterable_functor<char>& f) {
   return string<Syn>(this->array<char,Syn>::filter(f));
}

/*
 * Apply filter to const string.
 */
template <typename Syn>
string<Syn> string<Syn>::filter(const filterable_functor<const char>& f) const {
   return string<Syn>(this->array<char,Syn>::filter(f));
}

/*
 * Return the set of unique characters (as defined by the given comparison
 * functor) in the string (replacing the original string).  
 */
template <typename Syn>
void string<Syn>::unique(const comparable_functor<char>& f) {
   auto_write_lock<const Syn> wlock(*this);
   /* compute unique array */
   this->make_unique(f, 1);
   /* set terminator */
   this->_data[this->_size] = '\0';
}

/*
 * Return the set of unique characters (as defined by the given comparison
 * functor) in the string (replacing the original string).  In addition,
 * return an index array containing the original positions of the unique
 * characters.
 */
template <typename Syn>
array<unsigned long> string<Syn>::unique_idx(const comparable_functor<char>& f)
{
   auto_write_lock<const Syn> wlock(*this);
   /* compute unique array and index */
   array<unsigned long> idx = this->make_unique_idx(f, 1);
   /* set terminator */
   this->_data[this->_size] = '\0';
   return idx;
}

/*
 * Serialize.
 */
template <typename Syn>
void string<Syn>::serialize(serial_output_stream& s) const {
   this->array<char,Syn>::serialize(s);
}

/*
 * Deserialize.
 */
template <typename Syn>
auto_ptr< string<Syn> > string<Syn>::deserialize(serial_input_stream& s) {
   unsigned long size = 0;
   s >> size;
   auto_ptr< string<Syn> > str(new string<Syn>(size));
   for (unsigned long n = 0; n < size; n++) {
      s >> str->_data[n];
   }
   return str;
}

/*
 * Formatted output to stream.
 */
template <typename Syn>
ostream& operator<<(ostream& os, const string<Syn>& s) {
   auto_read_lock<const Syn> rlock(s);
   os << s._data;
   return os;
}

/*
 * string == char*.
 */
template <typename Syn>
bool operator==(const string<Syn>& s, const char* c) {
   return (s.compare_to(c) == 0);
}

/*
 * string != char*.
 */
template <typename Syn>
bool operator!=(const string<Syn>& s, const char* c) {
   return (s.compare_to(c) != 0);
}

/*
 * string < char*.
 */
template <typename Syn>
bool operator<(const string<Syn>& s, const char* c) {
   return (s.compare_to(c) < 0);
}

/*
 * string > char*.
 */
template <typename Syn>
bool operator>(const string<Syn>& s, const char* c) {
   return (s.compare_to(c) > 0);
}

/*
 * string <= char*.
 */
template <typename Syn>
bool operator<=(const string<Syn>& s, const char* c) {
   return (s.compare_to(c) <= 0);
}

/*
 * string >= char*.
 */
template <typename Syn>
bool operator>=(const string<Syn>& s, const char* c) {
   return (s.compare_to(c) >= 0);
}

/*
 * char* == string.
 */
template <typename Syn>
inline bool operator==(const char* c, const string<Syn>& s) {
   return (s == c);
}

/*
 * char* != string.
 */
template <typename Syn>
inline bool operator!=(const char* c, const string<Syn>& s) {
   return (s != c);
}

/*
 * char* < string.
 */
template <typename Syn>
inline bool operator<(const char* c, const string<Syn>& s) {
   return (s > c);
}

/*
 * char* > string.
 */
template <typename Syn>
inline bool operator>(const char* c, const string<Syn>& s) {
   return (s < c);
}

/*
 * char* <= string.
 */
template <typename Syn>
inline bool operator<=(const char* c, const string<Syn>& s) {
   return (s >= c);
}

/*
 * char* >= string.
 */
template <typename Syn>
inline bool operator>=(const char* c, const string<Syn>& s) {
   return (s <= c);
}

/*
 * string == string.
 */
template <typename Syn>
bool operator==(const string<Syn>& s0, const string<Syn>& s1) {
   return (s0.compare_to(s1) == 0);
}

/*
 * string != string.
 */
template <typename Syn>
bool operator!=(const string<Syn>& s0, const string<Syn>& s1) {
   return (s0.compare_to(s1) != 0);
}

/*
 * string < string.
 */
template <typename Syn>
bool operator<(const string<Syn>& s0, const string<Syn>& s1) {
   return (s0.compare_to(s1) < 0);
}

/*
 * string > string.
 */
template <typename Syn>
bool operator>(const string<Syn>& s0, const string<Syn>& s1) {
   return (s0.compare_to(s1) > 0);
}

/*
 * string <= string.
 */
template <typename Syn>
bool operator<=(const string<Syn>& s0, const string<Syn>& s1) {
   return (s0.compare_to(s1) <= 0);
}

/*
 * string >= string.
 */
template <typename Syn>
bool operator>=(const string<Syn>& s0, const string<Syn>& s1) {
   return (s0.compare_to(s1) >= 0);
}

/*
 * Comparison to char*.
 */
template <typename Syn>
int string<Syn>::compare_to(const char* c) const {
   auto_read_lock<const Syn> rlock(*this);
   return std::strcmp(this->_data, c);
}

/*
 * Comparison to string.
 */
template <typename Syn>
int string<Syn>::compare_to(const string<Syn>& s) const {
   auto_read_read_lock<const Syn> rrlock(*this, s);
   return std::strcmp(this->_data, s._data);
}

} /* namespace lang */

#endif
