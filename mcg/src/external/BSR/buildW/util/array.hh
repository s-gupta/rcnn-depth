// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.

#ifndef ARRAY_HH
#define ARRAY_HH

#include <assert.h>
//#include "message.hh"

// TODO: 
// - replace asserts with exceptions?
// - add typing to serialization here and in matlab

// Arrays that reduce bugs by:
//  - Being allocatable on the stack, so destructors get called 
//    automatically.
//  - Doing bounds checking.
//  - Providing easy initialization.
//  - Encapsulating the address calculation.
//
// The arrays are allocated as single blocks so that all elements are
// contiguous in memory.  Latter indices change more quickly than
// former indices.  Clients can rely on this ordering.
// 
// slice returns an const array of smaller dimension.  for effeciency, it
// doesn't copy the array contents but since it's const you won't accidentaly
// resize or delete it...
//


namespace Util
{

  template < class Elem > class Array1D;
  template < class Elem > class Array2D;
  template < class Elem > class Array3D;
  template < class Elem > class Array4D;

  template < class Elem > class Array1D
  {
      public:
        Array1D ()
        {
          _alloc (0);
        }

        Array1D(const Array1D<Elem>& a)
        {
          _alloc(a._n);
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = a._array[i];
          }
        }

        Array1D (unsigned int n)
        {
          _alloc (n);
        }

        ~Array1D ()
        {
          _delete ();
        }

        void resize (unsigned int n)
        {
          if (!issize (n))
          {
            _delete ();
            _alloc (n);
            //Message::debug(String("Array1D::resized to %d",_n),5);
          }
        }

        void init (const Elem & elem)
        {
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = elem;
          }
        }

        bool issize (unsigned int n) const
        {
          return (_n == n);
        }

        int size () const
        {
          return _n;
        }

        int size (unsigned int d) const
        {
          assert (d < 1);
          return _n;
        }

        Elem *data ()
        {
          return _array;
        }

        Elem & operator()(unsigned int i)
        {
          assert (i < _n);
          return _array[i];
        }

        const Elem & operator() (unsigned int i) const
        {
          assert (i < _n);
          return _array[i];
        }

        const Array1D<Elem>& operator=(const Array1D<Elem>& rhs) 
        { 
          if (this != &rhs) 
          { 
            resize(rhs._n);
            for (unsigned int i = 0; i < _n; i++)
            {
              _array[i] = rhs._array[i];
            }
          } 
          return *this; 
        }

#define DEFOP(OP) \
    Array1D<Elem>& operator OP (const Elem& a) \
    { \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP a; \
        } \
        return *this; \
    } \
    \
    Array1D<Elem>& operator OP (const Array1D<Elem>& that)\
    { \
        assert(size(0) == that.size(0)); \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP that._array[i]; \
        } \
        return *this; \
    }
DEFOP(+=)
DEFOP(-=)
DEFOP(*=)
DEFOP(/=)
#undef DEFOP

        friend std::istream & operator >> (std::istream & in, Array1D<Elem>& a)
        {
          unsigned int dim = 0;
          in.read((char*)&dim,sizeof(int));
          assert(dim == 1);

          unsigned int d0 = 0;
          in.read((char*)&d0,sizeof(int));
          a.resize(d0);

          in.read((char*)a._array,d0*sizeof(Elem));
         
          return in;
        }

        friend std::ostream & operator << (std::ostream & out, const Array1D<Elem>& a)
        {
          const unsigned int dim = 2;
          const unsigned int d0 = a.size(0); 
          out.write((char*)&dim,sizeof(int)); 
          out.write((char*)&d0,sizeof(int));
          out.write((char*)a._array,d0*sizeof(Elem));
          return out;
        }
   
      private:

        friend class Array2D<Elem>;
        friend class Array3D<Elem>;
        friend class Array4D<Elem>;
        Array1D (const Elem* data, unsigned int size)
        {
          _n = size;
          _array = const_cast<Elem*>(data);
        }

        void _alloc (unsigned int n)
        {
          _n = n;
          _array = new Elem[_n];
        }

        void _delete ()
        {
          assert (_array != NULL);
          delete[] _array;
          _array = NULL;
        }

        unsigned int _n;
        Elem *_array;
  };                            // class Array1D

  ///////////////////////////////////////////////////////////////////////////////

  template < class Elem > class Array2D
  {
      public:

        Array2D ()
        {
          _alloc (0, 0);
        }

        Array2D(const Array2D<Elem>& a)
        {
          _alloc(a._dim[0],a._dim[1]);
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = a._array[i];
          }
        }

        Array2D (unsigned int d0, unsigned int d1)
        {
          _alloc (d0, d1);
        }

        ~Array2D ()
        {
          _delete ();
        }

        void resize (unsigned int d0, unsigned int d1)
        {
          if (!issize (d0, d1))
          {
            _delete ();
            _alloc (d0, d1);
            //Message::debug(String("Array2D::resized to %d %d",d0,d1),5);
          }
        }

        void init (const Elem & elem)
        {
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = elem;
          }
        }

        bool issize (unsigned int d0, unsigned int d1) const
        {
          return (_dim[0] == d0 && _dim[1] == d1);
        }

        int size (unsigned int d) const
        {
          assert (d < 2);
          return _dim[d];
        }

        Elem *data ()
        {
          return _array;
        }

        // return a slice wrapped in a const array
        const Array1D<Elem>* slice(unsigned int i) const
        {
          return new Array1D<Elem>( &((*this)(i,0)) , _dim[1] );
        }

        Elem & operator() (unsigned int i, unsigned int j)
        {
          assert (i < _dim[0]);
          assert (j < _dim[1]);
          unsigned int index = i * _dim[1] + j;

          assert (index < _n);
          return _array[index];
        }

        const Elem & operator() (unsigned int i, unsigned int j) const
        {
          assert (i < _dim[0]);
          assert (j < _dim[1]);
          unsigned int index = i * _dim[1] + j;
            assert (index < _n);
            return _array[index];
        }

        const Array2D<Elem>& operator=(const Array2D<Elem>& rhs) 
        { 
          if (this != &rhs) 
          { 
            resize(rhs._dim[0],rhs._dim[1]);
            for (unsigned int i = 0; i < _n; i++)
            {
              _array[i] = rhs._array[i];
            }
          } 
          return *this; 
        }

#define DEFOP(OP) \
    Array2D<Elem>& operator OP (const Elem& a) \
    { \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP a; \
        } \
        return *this; \
    } \
    \
    Array2D<Elem>& operator OP (const Array2D<Elem>& that)\
    { \
        assert(size(0) == that.size(0)); \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP that._array[i]; \
        } \
        return *this; \
    }
DEFOP(+=)
DEFOP(-=)
DEFOP(*=)
DEFOP(/=)
#undef DEFOP

        friend std::istream & operator >> (std::istream & in, Array2D<Elem>& a)
        {
          unsigned int dim = 0;
          in.read((char*)&dim,sizeof(int));
          assert(dim == 2);

          unsigned int d0 = 0;
          unsigned int d1 = 0;
          in.read((char*)&d0,sizeof(int));
          in.read((char*)&d1,sizeof(int));
          a.resize(d0,d1);
          in.read((char*)a._array,d0*d1*sizeof(Elem));
         
          return in;
        }

        friend std::ostream & operator << (std::ostream & out, const Array2D<Elem>& a)
        {
          const unsigned int dim = 2;
          const unsigned int d0 = a.size(0); 
          const unsigned int d1 = a.size(1); 
          out.write((char*)&dim,sizeof(int)); 
          out.write((char*)&d0,sizeof(int));
          out.write((char*)&d1,sizeof(int));
          out.write((char*)a._array,d0*d1*sizeof(Elem));
          return out;
        }

      private:

        friend class Array3D<Elem>;
        friend class Array4D<Elem>;
        Array2D (const Elem* data, unsigned int d0, unsigned int d1) 
        {
          _array = const_cast<Elem*>(data);
          _dim[0] = d0;
          _dim[1] = d1;
          _n = _dim[0]*_dim[1];
        }

        void _alloc (unsigned int d0, unsigned int d1)
        {
          _n = d0 * d1;
          _dim[0] = d0;
          _dim[1] = d1;
          _array = new Elem[_n];
        }

        void _delete ()
        {
          assert (_array != NULL);
          delete[]_array;
          _array = NULL;
        }

        unsigned int _dim[2];
        unsigned int _n;
        Elem *_array;

  }; // class Array2D

  ///////////////////////////////////////////////////////////////////////////////

  template < class Elem > class Array3D
  {
      public:

        Array3D ()
        {
          _alloc (0, 0, 0);
        }

        Array3D(const Array3D<Elem>& a)
        {
          _alloc(a._dim[0],a._dim[1],a._dim[2]);
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = a._array[i];
          }
        }

        Array3D (unsigned int d0, unsigned int d1, unsigned int d2)
        {
          _alloc (d0, d1, d2);
        }

        ~Array3D ()
        {
          _delete ();
        }

        void resize (unsigned int d0, unsigned int d1, unsigned int d2)
        {
          if (!issize (d0, d1, d2))
          {
            _delete ();
            _alloc (d0, d1, d2);
            //Message::debug(String("Array3D::resized to %d %d %d",d0,d1,d2),5);
          }
        }

        void init (const Elem & elem)
        {
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = elem;
          }
        }

        bool issize (unsigned int d0, unsigned int d1, unsigned int d2) const
        {
          return (_dim[0] == d0 && _dim[1] == d1 && _dim[2] == d2);
        }

        int size (unsigned int d) const
        {
          assert (d < 3);
          return _dim[d];
        }

        Elem *data ()
        {
          return _array;
        }

        // return a slice wrappedped in a const array
        const Array2D<Elem>* slice(unsigned int i) const
        {
          return new Array2D<Elem>(&((*this)(i,0,0)),_dim[1],_dim[2]);
        }

        // return a slice wrappedped in a const array
        const Array1D<Elem>* slice(unsigned int i, unsigned int j) const
        {
          return new Array1D<Elem>(&((*this)(i,j,0)),_dim[2]);
        }

        Elem & operator() (unsigned int i, unsigned int j, unsigned int k)
        {
          assert (i < _dim[0]);
          assert (j < _dim[1]);
          assert (k < _dim[2]);
          unsigned int index = (i * _dim[1] + j) * _dim[2] + k;

          assert (index < _n);
          return _array[index];
        }

        const Elem & operator() (unsigned int i, unsigned int j, unsigned int k) const
        {
          assert (i < _dim[0]);
          assert (j < _dim[1]);
          assert (k < _dim[2]);
          unsigned int index = (i * _dim[1] + j) * _dim[2] + k;
            assert (index < _n);
            return _array[index];
        }

        const Array3D<Elem>& operator=(const Array3D<Elem>& rhs) 
        { 
          if (this != &rhs) 
          { 
            resize(rhs._dim[0],rhs._dim[1],rhs._dim[2]);
            for (unsigned int i = 0; i < _n; i++)
            {
              _array[i] = rhs._array[i];
            }
          } 
          return *this; 
        }

#define DEFOP(OP) \
    Array3D<Elem>& operator OP (const Elem& a) \
    { \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP a; \
        } \
        return *this; \
    } \
    \
    Array3D<Elem>& operator OP (const Array3D<Elem>& that)\
    { \
        assert(size(0) == that.size(0)); \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP that._array[i]; \
        } \
        return *this; \
    }
DEFOP(+=)
DEFOP(-=)
DEFOP(*=)
DEFOP(/=)
#undef DEFOP

        friend std::istream & operator >> (std::istream & in, Array3D<Elem>& a)
        {
          unsigned int dim = 0;
          in.read((char*)&dim,sizeof(int));
          assert(dim == 3);

          unsigned int d0 = 0;
          unsigned int d1 = 0;
          unsigned int d2 = 0;
          in.read((char*)&d0,sizeof(int));
          in.read((char*)&d1,sizeof(int));
          in.read((char*)&d2,sizeof(int));
          a.resize(d0,d1,d2);

          in.read((char*)a._array,d0*d1*d2*sizeof(Elem));
         
          return in;
        }

        friend std::ostream & operator << (std::ostream & out, const Array3D<Elem>& a)
        {
          const unsigned int dim = 3;
          const unsigned int d0 = a.size(0); 
          const unsigned int d1 = a.size(1); 
          const unsigned int d2 = a.size(2); 
          out.write((char*)&dim,sizeof(int)); 
          out.write((char*)&d0,sizeof(int));
          out.write((char*)&d1,sizeof(int));
          out.write((char*)&d2,sizeof(int));
          out.write((char*)a._array,d0*d1*d2*sizeof(Elem));
          return out;
        }

      private:

        friend class Array4D<Elem>;
        Array3D (const Elem* data, unsigned int d0, unsigned int d1, unsigned int d2)
        {
          _array = const_cast<Elem*>(data);
          _dim[0] = d0;
          _dim[1] = d1;
          _dim[2] = d2;
          _n = _dim[0]*_dim[1]*_dim[2];
        }


        void _alloc (unsigned int d0, unsigned int d1, unsigned int d2)
        {
          _n = d0 * d1 * d2;
          _array = new Elem[_n];
          _dim[0] = d0;
          _dim[1] = d1;
          _dim[2] = d2;
        }

        void _delete ()
        {
          assert (_array != NULL);
          delete[]_array;
          _array = NULL;
        }

        unsigned int _n;
        Elem *_array;
        unsigned int _dim[3];
  };                            // class Array3D


  ///////////////////////////////////////////////////////////////////////////////


  template < class Elem > class Array4D
  {

      public:
        Array4D ()
        {
          _alloc (0, 0, 0, 0);
        }

        Array4D(const Array4D<Elem>& a)
        {
          _alloc(a._dim[0],a._dim[1],a._dim[2],a._dim[3]);
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = a._array[i];
          }
        }

        Array4D (unsigned int d0, unsigned int d1, unsigned int d2, unsigned int d3)
        {
          _alloc (d0, d1, d2, d3);
        }

        ~Array4D ()
        {
          _delete ();
        }

        void resize (unsigned int d0, unsigned int d1, unsigned int d2, unsigned int d3)
        {
          if (!issize (d0, d1, d2, d3))
          {
            _delete ();
            _alloc (d0, d1, d2, d3);
            //Message::debug(String("Array4D::resized to %d %d %d %d",d0,d1,d2,d3),5);
          }
        }

        void init (const Elem & elem)
        {
          for (unsigned int i = 0; i < _n; i++)
          {
            _array[i] = elem;
          }
        }

        bool issize (unsigned int d0, unsigned int d1, unsigned int d2, unsigned int d3) const
        {
          return (_dim[0] == d0 && _dim[1] == d1 && _dim[2] == d2
                  && _dim[3] == d3);
        }

        int size (unsigned int d) const
        {
          assert (d < 4);
          return _dim[d];
        }

        Elem *data ()
        {
          return _array;
        }

        // return a slice wrapped in a const array
        const Array3D<Elem>* slice(unsigned int i) const 
        {
          return new Array3D<Elem>(&((*this)(i,0,0,0)),_dim[1],_dim[2],_dim[3]);
        }

        // return a slice wrapped in a const array
        const Array2D<Elem>* slice(unsigned int i, unsigned int j) const 
        {
          return new Array2D<Elem>(&((*this)(i,0,0,0)),_dim[2],_dim[3]);
        }

        // return a slice wrapped in a const array
        const Array1D<Elem>* slice(unsigned int i, unsigned int j, unsigned int k) const
        {
          return new Array1D<Elem>(&((*this)(i,j,k,0)),_dim[3]);
        }

        Elem & operator() (unsigned int i, unsigned int j, unsigned int k, unsigned int m)
        {
          assert (i < _dim[0]);
          assert (j < _dim[1]);
          assert (k < _dim[2]);
          assert (m < _dim[3]);
          unsigned int index = ((i * _dim[1] + j) * _dim[2] + k) * _dim[3] + m;
          assert (index < _n);
          return _array[index];
        }

        const Elem & operator() (unsigned int i, unsigned int j, unsigned int k,
                                           unsigned int m) const
        {
          assert (i < _dim[0]);
          assert (j < _dim[1]);
          assert (k < _dim[2]);
          assert (m < _dim[3]);
          unsigned int index = ((i * _dim[1] + j) * _dim[2] + k) * _dim[3] + m;
          assert (index < _n);
          return _array[index];
        }

        //assignment
        const Array4D<Elem>& operator=(const Array4D<Elem>& rhs) 
        { 
          if (this != &rhs) 
          { 
            resize(rhs._dim[0],rhs._dim[1],rhs._dim[2],rhs._dim[3]);
            for (unsigned int i = 0; i < _n; i++)
            {
              _array[i] = rhs._array[i];
            }
          } 
          return *this; 
        }

#define DEFOP(OP) \
    Array4D<Elem>& operator OP (const Elem& a) \
    { \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP a; \
        } \
        return *this; \
    } \
    \
    Array4D<Elem>& operator OP (const Array4D<Elem>& that)\
    { \
        assert(size(0) == that.size(0)); \
        for (unsigned int i = 0; i < _n; i++) { \
            _array[i] OP that._array[i]; \
        } \
        return *this; \
    }
DEFOP(+=)
DEFOP(-=)
DEFOP(*=)
DEFOP(/=)
#undef DEFOP

        friend std::istream & operator >> (std::istream & in, Array4D<Elem>& a)
        {
          unsigned int dim = 0;
          in.read((char*)&dim,sizeof(int));
          assert(dim == 4);

          unsigned int d0 = 0;
          unsigned int d1 = 0;
          unsigned int d2 = 0;
          unsigned int d3 = 0;
          in.read((char*)&d0,sizeof(int));
          in.read((char*)&d1,sizeof(int));
          in.read((char*)&d2,sizeof(int));
          in.read((char*)&d3,sizeof(int));
          a.resize(d0,d1,d2,d3);

          in.read((char*)a._array,d0*d1*d2*d3*sizeof(Elem));
         
          return in;
        }

        friend std::ostream & operator << (std::ostream & out, const Array4D<Elem>& a)
        {
          const unsigned int dim = 4;
          const unsigned int d0 = a.size(0); 
          const unsigned int d1 = a.size(1); 
          const unsigned int d2 = a.size(2); 
          const unsigned int d3 = a.size(3); 
          out.write((char*)&dim,sizeof(int)); 
          out.write((char*)&d0,sizeof(int));
          out.write((char*)&d1,sizeof(int));
          out.write((char*)&d2,sizeof(int));
          out.write((char*)&d3,sizeof(int));
          out.write((char*)a._array,d0*d1*d2*d3*sizeof(Elem));
          return out;
        }

      private:

        Array4D (const Elem* data, unsigned int d0, unsigned int d1, unsigned int d2, unsigned int d3)
        {
          _array = const_cast<Elem*>(data);
          _dim[0] = d0;
          _dim[1] = d1;
          _dim[2] = d2;
          _dim[3] = d3;
          _n = _dim[0]*_dim[1]*_dim[2]*_dim[3];
        }

        void _alloc (unsigned int d0, unsigned int d1, unsigned int d2, unsigned int d3)
        {
          _n = d0 * d1 * d2 * d3;
          _array = new Elem[_n];
          _dim[0] = d0;
          _dim[1] = d1;
          _dim[2] = d2;
          _dim[3] = d3;
        }

        void _delete ()
        {
          assert(_array != NULL);
          delete[]_array;
          _array = NULL;
        }

        unsigned int _n;
        Elem *_array;
        unsigned int _dim[4];

  }; // class Array4D


#define DEFOP(OP,IOP,CLASS) \
    template <class Elem> CLASS<Elem> operator OP (const CLASS<Elem> a, const Elem b) \
    { \
        CLASS<Elem> c = a; \
        return c IOP b; \
    }\
    \
    template <class Elem> CLASS<Elem> operator OP (const CLASS<Elem> a, const CLASS<Elem> b) \
    { \
        CLASS<Elem> c = a; \
        return c IOP b; \
    }
DEFOP(+,+=,Array1D)
DEFOP(-,-=,Array1D)
DEFOP(*,*=,Array1D)
DEFOP(/,/=,Array1D)
DEFOP(+,+=,Array2D)
DEFOP(-,-=,Array2D)
DEFOP(*,*=,Array2D)
DEFOP(/,/=,Array2D)
DEFOP(+,+=,Array3D)
DEFOP(-,-=,Array3D)
DEFOP(*,*=,Array3D)
DEFOP(/,/=,Array3D)
DEFOP(+,+=,Array4D)
DEFOP(-,-=,Array4D)
DEFOP(*,*=,Array4D)
DEFOP(/,/=,Array4D)
#undef DEFOP

} // namespace Util

#endif // array_h
