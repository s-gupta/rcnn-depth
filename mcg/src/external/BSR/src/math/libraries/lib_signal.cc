/*
 * Signal processing library.
 */
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/complex.hh"
#include "math/libraries/lib_signal.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "math/matrices/cmatrix.hh"

/* FIXME */
#include "lang/exceptions/ex_not_implemented.hh"

namespace math {
namespace libraries {
/*
 * Imports.
 */
using lang::array;
using lang::exceptions::ex_invalid_argument;
using math::complex;
using math::matrices::matrix;
using math::matrices::cmatrix;

/* FIXME */
using lang::exceptions::ex_not_implemented;

/***************************************************************************
 * Fourier transform.
 ***************************************************************************/

namespace {
/*
 * Compute dft/idft.
 */
void compute_dft(
   const complex<>* x,  /* source vector */
   complex<>* X,        /* destination vector */
   unsigned long N,     /* vector length */
   unsigned long step,  /* step between vector elements */
   bool is_inv)         /* computing inverse transform? */
{
   double exp_factor = is_inv ?
      ( 2*M_PIl/static_cast<double>(N)) :
      (-2*M_PIl/static_cast<double>(N));
   for (unsigned long k = 0, k_ind = 0; k < N; k++, k_ind += step) {
      complex<> iw(0, exp_factor*static_cast<double>(k));
      for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step)
         X[k_ind] += x[n_ind] * math::exp(iw*static_cast<double>(n));
   }
   if (is_inv) {
      for (unsigned long k = 0, k_ind = 0; k < N; k++, k_ind += step)
         X[k_ind] /= static_cast<double>(N);
   }
}

/*
 * Check whether argument is a power of two.
 */
bool is_power_of_two(unsigned long n) {
   bool is_pwr_two = (n > 0);
   while ((n > 1) && (is_pwr_two)) {
      unsigned long half_n = n/2;
      is_pwr_two = (n == (half_n*2));
      n = half_n;
   }
   return is_pwr_two;
}

/*
 * Compute fft/ifft.
 */
void compute_fft(
   complex<>* x,        /* vector to transform (in place) */
   complex<>* temp,     /* temporary vector (for rearrangement) */ 
   unsigned long N,     /* vector length (must be a power of 2) */
   unsigned long step,  /* step between vector elements */
   complex<> omega)     /* Nth root of unity */
{
   if (N > 1) {
      /* place even elements in first half, odd elements in second */
      unsigned long half_N = N/2;
      unsigned long twice_step = step*2;
      for (unsigned long n = 0, pos = 0; n < half_N; n++, pos += twice_step)
         temp[n] = x[pos];
      for (unsigned long n = half_N, pos = step; n < N; n++, pos += twice_step)
         temp[n] = x[pos];
      for (unsigned long n = 0, pos = 0; n < N; n++, pos += step)
         x[pos] = temp[n];
      /* recursively compute fft */
      unsigned long pos_u = 0;
      unsigned long pos_v = step*half_N;
      complex<> omega_sq = omega*omega;
      compute_fft(x,           temp,            half_N, step, omega_sq);
      compute_fft(&(x[pos_v]), &(temp[half_N]), half_N, step, omega_sq);
      /* merge results */
      complex<> omega_curr(1,0);
      for (unsigned long n = 0; n < half_N; n++) {
         complex<> u = x[pos_u];
         complex<> v = x[pos_v] * omega_curr;
         x[pos_u] = u + v;
         x[pos_v] = u - v;
         pos_u += step;
         pos_v += step;
         omega_curr *= omega;
      }
   }
}
} /* namespace */

/*
 * Discrete Fourier transform (DFT) and inverse DFT (of real-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 *
 * WARNING: This is an O(n^2) procedure.  In most cases, using the
 * O(n*log(n)) fft or ifft function instead is likely to be preferable.
 */
cmatrix<> lib_signal::dft(const matrix<>& m, unsigned long d) {
   return lib_signal::dft(cmatrix<>(m), d);
}

cmatrix<> lib_signal::idft(const matrix<>& m, unsigned long d) {
   return lib_signal::idft(cmatrix<>(m), d);
}

/*
 * Discrete Fourier transform (DFT) and inverse DFT (of complex-valued
 * matrix).  Compute the transform along the specified dimension of the
 * matrix.
 *
 * WARNING: This is an O(n^2) procedure.  In most cases, using the
 * O(n*log(n)) fft or ifft function instead is likely to be preferable.
 */
cmatrix<> lib_signal::dft(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_dft(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size, false
         );
      }
   }
   return m_result;
}

cmatrix<> lib_signal::idft(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_dft(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size, true
         );
      }
   }
   return m_result;
}

/*
 * Fast Fourier transform (FFT) and inverse FFT (of real-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 */
cmatrix<> lib_signal::fft(const matrix<>& m, unsigned long d) {
   return lib_signal::fft(cmatrix<>(m), d);
}

cmatrix<> lib_signal::ifft(const matrix<>& m, unsigned long d) {
   return lib_signal::ifft(cmatrix<>(m), d);
}

/*
 * Fast Fourier transform (FFT) and inverse FFT (of complex-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 */
cmatrix<> lib_signal::fft(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m_result._dims.size();
      unsigned long n_steps = (d < n_dims) ? m_result._dims[d] : 1;
      /* check that vector length is a power of two */
      if ((n_steps > 0) && (!is_power_of_two(n_steps)))
         throw ex_invalid_argument(
            "fft called on vector with non-power of two length"
         );
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(
         m_result._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m_result._dims, d);
      /* compute nth root of unity */
      complex<> omega = math::exp(
         complex<>(0, -2*M_PIl/static_cast<double>(n_steps))
      );
      /* compute result matrix */
      cmatrix<> temp(n_steps, 1);
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_fft(
            &(m_result._data[ind]), temp._data, n_steps, step_size, omega
         );
      }
   }
   return m_result;
}

cmatrix<> lib_signal::ifft(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m_result._dims.size();
      unsigned long n_steps = (d < n_dims) ? m_result._dims[d] : 1;
      /* check that vector length is a power of two */
      if ((n_steps > 0) && (!is_power_of_two(n_steps)))
         throw ex_invalid_argument(
            "ifft called on vector with non-power of two length"
         );
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(
         m_result._dims, d
      );
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m_result._dims, d);
      /* compute nth root of unity */
      complex<> omega = math::exp(
         complex<>(0, 2*M_PIl/static_cast<double>(n_steps))
      );
      /* compute result matrix */
      cmatrix<> temp(n_steps, 1);
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         /* compute inverse transform */
         compute_fft(
            &(m_result._data[ind]), temp._data, n_steps, step_size, omega
         );
         /* scale by 1/N */
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            m_result._data[ind] /= static_cast<double>(n_steps);
            ind += step_size;
         }
      }
   }
   return m_result;
}

/***************************************************************************
 * Cosine transform.
 ***************************************************************************/

namespace {
/*
 * Compute dct.
 */
template <typename T>
void compute_dct(
   const T* x,          /* source vector */
   T* X,                /* destination vector */
   unsigned long N,     /* vector length */
   unsigned long step)  /* step between vector elements */
{
   double pi_N = M_PIl/static_cast<double>(N);
   for (unsigned long k = 0, k_ind = 0; k < N; k++, k_ind += step) {
      double pik_N  = pi_N * static_cast<double>(k);
      double pik_2N = 0.5 * pik_N;
      for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step) {
         X[k_ind] +=
            x[n_ind] * math::cos(pik_N * static_cast<double>(n) + pik_2N);
      }
   }
}

/*
 * Compute idct.
 */
template <typename T>
void compute_idct(
   const T* x,          /* source vector */
   T* X,                /* destination vector */
   unsigned long N,     /* vector length */
   unsigned long step)  /* step between vector elements */
{
   double pi_N  = M_PIl/static_cast<double>(N);
   double pi_2N = 0.5 * pi_N;
   for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step) {
      double coeff = pi_N * static_cast<double>(n) + pi_2N;
      X[n_ind] += 0.5 * x[0];
      for (unsigned long k = 1, k_ind = step; k < N; k++, k_ind += step)
         X[n_ind] += x[k_ind] * math::cos(coeff * static_cast<double>(k));
   }
   double mult_factor = 2.0/static_cast<double>(N);
   for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step)
      X[n_ind] *= mult_factor;
}
} /* namespace */

/*
 * Discrete cosine transform (DCT) and inverse DCT (of real-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 *
 * WARNING: This is an O(n^2) procedure.  In most cases, using the
 * O(n*log(n)) fct or ifct function instead is likely to be preferable.
 */
matrix<> lib_signal::dct(const matrix<>& m, unsigned long d) {
   matrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_dct(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

matrix<> lib_signal::idct(const matrix<>& m, unsigned long d) {
   matrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_idct(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

/*
 * Discrete cosine transform (DCT) and inverse DCT (of complex-valued
 * matrix).  Compute the transform along the specified dimension of the
 * matrix.
 *
 * WARNING: This is an O(n^2) procedure.  In most cases, using the
 * O(n*log(n)) fct or ifct function instead is likely to be preferable.
 */
cmatrix<> lib_signal::dct(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_dct(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

cmatrix<> lib_signal::idct(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_idct(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

/*
 * Fast cosine transform (FCT) and inverse FCT (of real-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 */
matrix<> lib_signal::fct(const matrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

matrix<> lib_signal::ifct(const matrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

/*
 * Fast cosine transform (FCT) and inverse FCT (of complex-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 */
cmatrix<> lib_signal::fct(const cmatrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

cmatrix<> lib_signal::ifct(const cmatrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

/***************************************************************************
 * Sine transform.
 ***************************************************************************/

namespace {
/*
 * Compute dst.
 */
template <typename T>
void compute_dst(
   const T* x,          /* source vector */
   T* X,                /* destination vector */
   unsigned long N,     /* vector length */
   unsigned long step)  /* step between vector elements */
{
   double pi_N = M_PIl/static_cast<double>(N);
   for (unsigned long k = 1, k_ind = 0; k <= N; k++, k_ind += step) {
      double pik_N  = pi_N * static_cast<double>(k);
      double pik_2N = 0.5 * pik_N;
      for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step) {
         X[k_ind] +=
            x[n_ind] * math::sin(pik_N * static_cast<double>(n) + pik_2N);
      }
   }
}

/*
 * Compute idst.
 */
template <typename T>
void compute_idst(
   const T* x,          /* source vector */
   T* X,                /* destination vector */
   unsigned long N,     /* vector length */
   unsigned long step)  /* step between vector elements */
{
   double pi_N  = M_PIl/static_cast<double>(N);
   double pi_2N = 0.5 * pi_N;
   for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step) {
      double coeff = pi_N * static_cast<double>(n) + pi_2N;
      unsigned long k_ind = 0;
      for (unsigned long k = 1; k < N; k++, k_ind += step)
         X[n_ind] += x[k_ind] * math::sin(coeff * static_cast<double>(k));
      X[n_ind] += 0.5 * x[k_ind] * math::sin(coeff * static_cast<double>(N));
   }
   double mult_factor = 2.0/static_cast<double>(N);
   for (unsigned long n = 0, n_ind = 0; n < N; n++, n_ind += step)
      X[n_ind] *= mult_factor;
}
} /* namespace */

/*
 * Discrete sine transform (DST) and inverse DST (of real-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 *
 * WARNING: This is an O(n^2) procedure.  In most cases, using the
 * O(n*log(n)) fst or ifst function instead is likely to be preferable.
 */
matrix<> lib_signal::dst(const matrix<>& m, unsigned long d) {
   matrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_dst(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

matrix<> lib_signal::idst(const matrix<>& m, unsigned long d) {
   matrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_idst(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

/*
 * Discrete sine transform (DST) and inverse DST (of complex-valued
 * matrix).  Compute the transform along the specified dimension of the
 * matrix.
 *
 * WARNING: This is an O(n^2) procedure.  In most cases, using the
 * O(n*log(n)) fst or ifst function instead is likely to be preferable.
 */
cmatrix<> lib_signal::dst(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_dst(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

cmatrix<> lib_signal::idst(const cmatrix<>& m, unsigned long d) {
   cmatrix<> m_result(m._dims);
   if (m_result._size > 0) {
      /* get # steps along dimension */
      unsigned long n_dims  = m._dims.size();
      unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         compute_idst(
            &(m._data[ind]), &(m_result._data[ind]), n_steps, step_size
         );
      }
   }
   return m_result;
}

/*
 * Fast sine transform (FST) and inverse FST (of real-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 */
matrix<> lib_signal::fst(const matrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

matrix<> lib_signal::ifst(const matrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

/*
 * Fast sine transform (FST) and inverse FST (of complex-valued matrix).
 * Compute the transform along the specified dimension of the matrix.
 */
cmatrix<> lib_signal::fst(const cmatrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

cmatrix<> lib_signal::ifst(const cmatrix<>& m, unsigned long d) {
   /* FIXME */
   throw ex_not_implemented();
}

/***************************************************************************
 * Hilbert transform.
 ***************************************************************************/

/*
 * Hilbert transform (along specified dimension).
 */
matrix<> lib_signal::hilbert(const matrix<>& m, unsigned long d) {
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   unsigned long half_n_steps = n_steps/2;
   /* compute fourier transform */
   bool is_pwr_two = is_power_of_two(n_steps);
   cmatrix<> m_dft =
      is_pwr_two ? lib_signal::fft(m, d) : lib_signal::dft(m, d);
   /* compute indices along first slice of dimension */
   array<unsigned long> inds = matrix<>::linear_indices_slice(
      m._dims, d
   );
   /* compute step size along dimension */
   unsigned long step_size = matrix<>::dims_size(m._dims, d);
   /* double positive frequencies   */
   /* zero out negative frequencies */
   /* leave dc component            */
   bool is_mult_two = (n_steps == (half_n_steps*2));
   unsigned long pos_freq_bound_ind = (n_steps+1)/2;
   unsigned long n_inds = inds.size();
   for (unsigned long n = 0; n < n_inds; n++) {
      unsigned long ind = inds[n];
      /* double positive frequencies */
      ind += step_size;
      for (unsigned long n_step = 1; n_step < pos_freq_bound_ind; n_step++)
      {
         m_dft._data[ind] *= 2;
         ind += step_size;
      }
      /* zero out negative frequencies */
      ind += (is_mult_two ? step_size : 0);
      for (unsigned long n_step = (half_n_steps+1); n_step < n_steps; n_step++)
      {
         m_dft._data[ind] = 0;
         ind += step_size;
      }
   }
   /* compute inverse fourier transform */
   cmatrix<> m_hilbert =
      is_pwr_two ? lib_signal::ifft(m_dft, d) : lib_signal::idft(m_dft, d);
   /* grab and return imaginary component */
   return m_hilbert.imag();
}

} /* namespace libraries */
} /* namespace math */
