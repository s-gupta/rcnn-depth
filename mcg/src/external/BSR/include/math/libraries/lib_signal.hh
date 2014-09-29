/*
 * Signal processing library.
 */
#ifndef MATH__LIBRARIES__LIB_SIGNAL_HH
#define MATH__LIBRARIES__LIB_SIGNAL_HH

#include "math/matrices/matrix.hh"
#include "math/matrices/cmatrix.hh"

namespace math {
namespace libraries {
/*
 * Imports.
 */
using math::matrices::matrix;
using math::matrices::cmatrix;

/*
 * Signal processing functions.
 */
class lib_signal {
public:
   /************************************************************************
    * Fourier transform.
    *
    * For an input vector x of length N, the discrete Fourier transform is 
    * the length N vector X given by
    *
    *                 N-1
    *    X(k) =       sum x(n)*exp(-i*2*pi*k*n/N), for 0 <= k <= N-1
    *                 n=0
    *
    * The inverse discrete Fourier transform is given by
    *
    *                 N-1
    *    x(n) = (1/N) sum X(k)*exp( i*2*pi*k*n/N), for 0 <= n <= N-1
    *                 k=0
    *
    * The O(n*log(n)) fast transforms (fft, ifft) only work on inputs whose 
    * length n is a power of two.
    ************************************************************************/
    
   /*
    * Discrete Fourier transform (DFT) and inverse DFT (of real-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    *
    * WARNING: This is an O(n^2) procedure.  In most cases, using the
    * O(n*log(n)) fft or ifft function instead is likely to be preferable.
    */
   static cmatrix<>  dft(const matrix<>&, unsigned long = 0);
   static cmatrix<> idft(const matrix<>&, unsigned long = 0);

   /*
    * Discrete Fourier transform (DFT) and inverse DFT (of complex-valued
    * matrix).  Compute the transform along the specified dimension of the
    * matrix.
    *
    * WARNING: This is an O(n^2) procedure.  In most cases, using the
    * O(n*log(n)) fft or ifft function instead is likely to be preferable.
    */
   static cmatrix<>  dft(const cmatrix<>&, unsigned long = 0);
   static cmatrix<> idft(const cmatrix<>&, unsigned long = 0);

   /*
    * Fast Fourier transform (FFT) and inverse FFT (of real-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    */
   static cmatrix<>  fft(const matrix<>&, unsigned long = 0);
   static cmatrix<> ifft(const matrix<>&, unsigned long = 0);

   /*
    * Fast Fourier transform (FFT) and inverse FFT (of complex-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    */
   static cmatrix<>  fft(const cmatrix<>&, unsigned long = 0);
   static cmatrix<> ifft(const cmatrix<>&, unsigned long = 0);
   
   /************************************************************************
    * Cosine transform.
    *
    * For an input vector x of length N, the discrete cosine transform is 
    * the length N vector X given by
    *
    *                 N-1
    *    X(k) =       sum x(n)*cos(pi*(n+1/2)*k/N),      for 0 <= k <= N-1
    *                 n=0
    *
    * The inverse discrete cosine transform is given by
    *
    *                 N-1
    *    x(n) = (2/N) sum w(k)*X(k)*cos(pi*(n+1/2)*k/N), for 0 <= n <= N-1
    *                 k=0
    *
    * where w(k) = 1/2 for k = 0, and w(k) = 1 for k > 0.
    *
    * The O(n*log(n)) fast transforms (fct, ifct) only work on inputs whose 
    * length n is a power of two.
    ************************************************************************/

   /*
    * Discrete cosine transform (DCT) and inverse DCT (of real-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    *
    * WARNING: This is an O(n^2) procedure.  In most cases, using the
    * O(n*log(n)) fct or ifct function instead is likely to be preferable.
    */
   static matrix<>  dct(const matrix<>&, unsigned long = 0);
   static matrix<> idct(const matrix<>&, unsigned long = 0);

   /*
    * Discrete cosine transform (DCT) and inverse DCT (of complex-valued
    * matrix).  Compute the transform along the specified dimension of the
    * matrix.
    *
    * WARNING: This is an O(n^2) procedure.  In most cases, using the
    * O(n*log(n)) fct or ifct function instead is likely to be preferable.
    */
   static cmatrix<>  dct(const cmatrix<>&, unsigned long = 0);
   static cmatrix<> idct(const cmatrix<>&, unsigned long = 0);

   /*
    * Fast cosine transform (FCT) and inverse FCT (of real-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    */
   static matrix<>  fct(const matrix<>&, unsigned long = 0);
   static matrix<> ifct(const matrix<>&, unsigned long = 0);

   /*
    * Fast cosine transform (FCT) and inverse FCT (of complex-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    */
   static cmatrix<>  fct(const cmatrix<>&, unsigned long = 0);
   static cmatrix<> ifct(const cmatrix<>&, unsigned long = 0);
   
   /************************************************************************
    * Sine transform.
    *
    * For an input vector x of length N, the discrete sine transform is 
    * the length N vector X given by
    *
    *                 N-1
    *    X(k) =       sum x(n)*sin(pi*(n+1/2)*(k+1)/N),      for 0 <= k <= N-1
    *                 n=0
    *
    * The inverse discrete sine transform is given by
    *
    *                 N-1
    *    x(n) = (2/N) sum w(k)*X(k)*sin(pi*(n+1/2)*(k+1)/N), for 0 <= n <= N-1
    *                 k=0
    *
    * where w(k) = 1 for k < N-1, and w(k) = 1/2 for k = N-1.
    *
    * The O(n*log(n)) fast transforms (fst, ifst) only work on inputs whose 
    * length n is a power of two.
    ************************************************************************/

   /*
    * Discrete sine transform (DST) and inverse DST (of real-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    *
    * WARNING: This is an O(n^2) procedure.  In most cases, using the
    * O(n*log(n)) fst or ifst function instead is likely to be preferable.
    */
   static matrix<>  dst(const matrix<>&, unsigned long = 0);
   static matrix<> idst(const matrix<>&, unsigned long = 0);

   /*
    * Discrete sine transform (DST) and inverse DST (of complex-valued
    * matrix).  Compute the transform along the specified dimension of the
    * matrix.
    *
    * WARNING: This is an O(n^2) procedure.  In most cases, using the
    * O(n*log(n)) fst or ifst function instead is likely to be preferable.
    */
   static cmatrix<>  dst(const cmatrix<>&, unsigned long = 0);
   static cmatrix<> idst(const cmatrix<>&, unsigned long = 0);

   /*
    * Fast sine transform (FST) and inverse FST (of real-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    */
   static matrix<>  fst(const matrix<>&, unsigned long = 0);
   static matrix<> ifst(const matrix<>&, unsigned long = 0);

   /*
    * Fast sine transform (FST) and inverse FST (of complex-valued matrix).
    * Compute the transform along the specified dimension of the matrix.
    */
   static cmatrix<>  fst(const cmatrix<>&, unsigned long = 0);
   static cmatrix<> ifst(const cmatrix<>&, unsigned long = 0);

   /************************************************************************
    * Hilbert transform.
    *
    * If the matrix size is a power of two, the transform is computed in
    * O(n*log(n)) time, otherwise it takes O(n^2) time.
    ************************************************************************/

   /*
    * Hilbert transform (along specified dimension).
    */
   static matrix<> hilbert(const matrix<>&, unsigned long = 0);
};

} /* namespace libraries */
} /* namespace math */

#endif
