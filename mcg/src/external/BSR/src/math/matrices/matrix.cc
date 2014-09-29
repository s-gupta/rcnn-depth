/*
 * Matrix.
 */
#include "lang/array.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "math/matrices/exceptions/ex_matrix_dimension_mismatch.hh"
#include "math/matrices/matrix.hh"

namespace math {
namespace matrices {
/*
 * Imports.
 */
using lang::array;
using lang::exceptions::ex_invalid_argument;
using math::matrices::exceptions::ex_matrix_dimension_mismatch;

/*
 * Pure virtual destructor.
 */
matrix_base::~matrix_base() { }

/*
 * Compute which dimension is shortest (return lower id in case of a tie).
 */
unsigned long matrix_base::dims_min(const array<unsigned long>& dims) {
   unsigned long n_dims = dims.size();
   unsigned long d      = 0;
   unsigned long d_size = 0;
   for (unsigned long n = 0; n < n_dims; n++) {
      if (dims[n] < d_size) {
         d = n;
         d_size = dims[n];
      }
   }
   return d;
}

/*
 * Compute which dimension is longest (return lower id in case of a tie).
 */
unsigned long matrix_base::dims_max(const array<unsigned long>& dims) {
   unsigned long n_dims = dims.size();
   unsigned long d      = 0;
   unsigned long d_size = (n_dims > 0) ? dims[0] : 0;
   for (unsigned long n = 1; n < n_dims; n++) {
      if (dims[n] > d_size) {
         d = n;
         d_size = dims[n];
      }
   }
   return d;
}

/*
 * Compute size of matrix with given dimensions.
 */
unsigned long matrix_base::dims_matrix_size(const array<unsigned long>& dims) {
   unsigned long n_dims = dims.size();
   unsigned long size = (n_dims > 0) ? 1 : 0;
   for (unsigned long n = 0; n < n_dims; n++)
      size *= dims[n];
   return size;
}

/*
 * Compute size of a step along the specified dimension (change in linear
 * index).
 */
unsigned long matrix_base::dims_size(
   const array<unsigned long>& dims, unsigned long d)
{
   unsigned long n_dims = dims.size();
   unsigned long size = 1;
   if (n_dims > 0) {
      for (unsigned long n = n_dims - 1; n > d; n--)
         size *= dims[n];
   }
   return size;
}

/*
 * Compute size of a step in each dimension (change in linear index).
 */
array<unsigned long> matrix_base::dims_sizes(
   const array<unsigned long>& dims)
{
   unsigned long n_dims = dims.size();
   array<unsigned long> sizes(n_dims);
   if (n_dims > 0) {
      unsigned long size = 1;
      for (unsigned long n = n_dims - 1; n > 0; n--) {
         sizes[n] = size;
         size *= dims[n];
      }
      sizes[0] = size;
   }  
   return sizes;
}

/*
 * Check if dimensions are equal.
 */
bool matrix_base::dims_equal(
   const array<unsigned long>& dims0,
   const array<unsigned long>& dims1) 
{
   unsigned long n_dims0 = dims0.size();
   unsigned long n_dims1 = dims1.size();
   if (n_dims0 == n_dims1) {
      for (unsigned long n = 0; n < n_dims0; n++) {
         if (dims0[n] != dims1[n])
            return false;
      }
      return true;
   }
   return false;
}

/*
 * Assert that dimensions are equal.
 * Throw an exception (ex_matrix_dimension_mismatch) if they are not equal.
 */
void matrix_base::assert_dims_equal(
   const array<unsigned long>& dims0, 
   const array<unsigned long>& dims1)
{
   if (!(matrix_base::dims_equal(dims0, dims1)))
      throw ex_matrix_dimension_mismatch();
}

/*
 * Check if dimensions make matrices compatible for multiplcation.
 */
bool matrix_base::dims_mult_compatible(
   const array<unsigned long>& dims_left,
   const array<unsigned long>& dims_right)
{
   unsigned long n_dims_left  = dims_left.size();
   unsigned long n_dims_right = dims_right.size();
   unsigned long last_dim_left =
      (n_dims_left < 2)   ? 1 : dims_left[n_dims_left-1];
   unsigned long first_dim_right =
      (n_dims_right == 0) ? 0 : dims_right[0];
   return (last_dim_left == first_dim_right);
}

/*
 * Compute dimensions resulting from matrix multiplication.
 * Throw an exception (ex_matrix_dimension_mismatch) if the dimensions 
 * are incompatible.
 */
array<unsigned long> matrix_base::dims_mult(
   const array<unsigned long>& dims_left,
   const array<unsigned long>& dims_right)
{
   if (matrix_base::dims_mult_compatible(dims_left, dims_right)) {
      /* get number of dimensions of left/right matrices */
      unsigned long n_dims_left  = dims_left.size();
      unsigned long n_dims_right = dims_right.size();
      /* compute number of dimensions each matrix contributes */
      unsigned long n_dims_from_left =
         (n_dims_left > 2)  ? (n_dims_left - 1)  : 1;
      unsigned long n_dims_from_right =
         (n_dims_right > 2) ? (n_dims_right - 1) : 1;
      unsigned long n_dims = n_dims_from_left + n_dims_from_right;
      array<unsigned long> dims(n_dims);
      /* set dimensions contributed by left matrix */
      if (n_dims_left == 0) {
         dims[0] = 0;
      } else {
         for (unsigned long n = 0; n < n_dims_from_left; n++)
            dims[n] = dims_left[n];
      }
      /* set dimensions contributed by right matrix */
      if (n_dims_right == 0) {
         dims[n_dims_from_left] = 1;
      } else if (n_dims_right == 1) {
         dims[n_dims_from_left] = 1;
      } else {
         for (unsigned long n = 0; n < n_dims_from_right; n++)
            dims[n_dims_from_left + n] = dims_right[n + 1];
      }
      return dims;
   } else {
      throw ex_matrix_dimension_mismatch(
         "matrix dimension mismatch during multiply"
      );
   }
}

/*
 * Extend the dimensions array to at least the given number of dimensions
 * by appending dimensions of size one (as needed).
 */
array<unsigned long> matrix_base::extend_dims(
   const array<unsigned long>& dims,
   unsigned long               n_dims_min)
{
   unsigned long n_dims = dims.size();
   array<unsigned long> dims_result(
      (n_dims > n_dims_min) ? n_dims : n_dims_min
   );
   for (unsigned long n = 0; n < n_dims; n++)
      dims_result[n] = dims[n];
   for (unsigned long n = n_dims; n < n_dims_min; n++)
      dims_result[n] = 1;
   return dims_result;
}

/*
 * Compute the range of indices in a block.
 */
array< array<unsigned long> > matrix_base::range_indices(
   const array<unsigned long>& start,
   const array<unsigned long>& end)
{
   array<unsigned long> step(start.size(), 1);
   return matrix_base::range_indices(start, step, end);
}

/*
 * Compute the range of indices in a block.
 */
array< array<unsigned long> > matrix_base::range_indices(
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end)
{
   /* check that all index arrays have the same size */
   unsigned long n_dims = start.size();
   if ((step.size() != n_dims) || (end.size() != n_dims))
      throw ex_invalid_argument(
         "start, step, and end arrays must all have the same size"
      );
   /* generate indices of desired elements */
   array< array<unsigned long> > range(n_dims);
   for (unsigned long d = 0; d < n_dims; d++) {
      /* get start, step, end along dimension d */
      unsigned long i_curr = start[d];
      unsigned long i_step = step[d];
      unsigned long i_end  = end[d];
      /* check step size */
      if (i_step == 0)
         throw ex_invalid_argument(
            "step cannot be zero along any dimension"
         );
      /* generate indices along dimension d */
      unsigned long i_size = (i_curr <= i_end) ?
         ((i_end - i_curr)/i_step + 1) : 0;
      array<unsigned long> i(i_size);
      for (unsigned long n = 0; n < i_size; n++, i_curr += i_step)
         i[n] = i_curr;
      /* store indices */
      range[d] = i;
   }
   return range;
}

/*
 * Compute linear index from two indices.
 */
unsigned long matrix_base::linear_index(
   const array<unsigned long>& dims, 
   unsigned long x, 
   unsigned long y)
{
   unsigned long n_dims = dims.size();
   unsigned long y_size = (n_dims > 1) ? dims[1] : 1;
   unsigned long index = x*y_size + y;
   for (unsigned long n = 2; n < n_dims; n++)
      index *= dims[n];
   return index;
}

/*
 * Compute linear index from multiple indices.
 */
unsigned long matrix_base::linear_index(
   const array<unsigned long>& dims,
   const array<unsigned long>& i)
{
   unsigned long n_dims     = dims.size();
   unsigned long n_dims_i   = i.size();
   unsigned long n_dims_min = (n_dims < n_dims_i) ? n_dims : n_dims_i;
   unsigned long index = 0;
   for (unsigned long n = 0; n < n_dims_min; n++)
      index = index*dims[n] + i[n];
   for (unsigned long n = n_dims_min; n < n_dims; n++)
      index *= dims[n];
   for (unsigned long n = n_dims_min; n < n_dims_i; n++)
      index += i[n];
   return index;
}

/*
 * Compute linear indices along the specified slice of the given dimenion.
 * These are the linear indices for all elements for which the coordinate 
 * along the given dimension is the given value.
 */
array<unsigned long> matrix_base::linear_indices_slice(
   const array<unsigned long>& dims,   /* dimensions */
   unsigned long               d,      /* dimension along which to slice */
   unsigned long               x)      /* coordinate of slice */
{
   array<unsigned long> dims_extended = matrix_base::extend_dims(dims, d + 1);
   unsigned long n_dims = dims_extended.size();
   array<unsigned long> start(n_dims);
   array<unsigned long> end(n_dims);
   for (unsigned long n = 0; n < n_dims; n++) {
      if (n == d) {
         start[d] = x;
         end[d] = x;
      } else {
         if (dims_extended[n] == 0)
            return array<unsigned long>();
         else
            end[n] = dims_extended[n] - 1;
      }
   }
   return matrix_base::linear_indices(dims_extended, start, end);
}

/*
 * Compute linear indices for a block of indices.
 */
array<unsigned long> matrix_base::linear_indices(
   const array<unsigned long>& dims,
   const array<unsigned long>& start,
   const array<unsigned long>& end)
{
   array< array<unsigned long> > range = range_indices(start, end);
   return matrix_base::linear_indices(dims, range);
}

/*
 * Compute linear indices for a block of indices.
 */
array<unsigned long> matrix_base::linear_indices(
   const array<unsigned long>& dims,
   const array<unsigned long>& start,
   const array<unsigned long>& step,
   const array<unsigned long>& end)
{
   array< array<unsigned long> > range = range_indices(start, step, end);
   return matrix_base::linear_indices(dims, range);
}

/*
 * Compute linear indices for an arbitrary range of indices.
 */
array<unsigned long> matrix_base::linear_indices(
   const array<unsigned long>&          dims,
   const array< array<unsigned long> >& range)
{
   array<unsigned long> step_sizes = matrix_base::dims_sizes(dims);
   return matrix_base::linear_indices_custom(step_sizes, range);
}

/*
 * Compute linear indices for an arbitrary range of indices and an
 * arbitrary matrix memory layout.  Use the specified step sizes along 
 * each dimension rather than assuming the standard memory layout.
 */
array<unsigned long> matrix_base::linear_indices_custom(
   const array<unsigned long>&          step_sizes,
   const array< array<unsigned long> >& range)
{
   /* get number of dimensions */
   unsigned long n_dims   = step_sizes.size();
   unsigned long n_dims_i = range.size();
   /* compute total number of indices */
   unsigned long n_inds = (n_dims_i > 0) ? 1 : 0;
   for (unsigned long n = 0; n < n_dims_i; n++)
      n_inds *= range[n].size();
   /* compute linear indices */
   array<unsigned long> inds(n_inds);
   if (n_inds > 0) {
      /* extend dimension size array (if needed) */
      array<unsigned long> d_sizes(step_sizes);
      d_sizes.resize(n_dims_i);
      for (unsigned long d = n_dims; d < n_dims_i; d++)
         d_sizes[d] = 1;
      /* initialize index at first position */
      array<unsigned long> pos(n_dims_i);
      unsigned long index = 0;
      for (unsigned long d = 0; d < n_dims_i; d++)
         index += (range[d])[0] * d_sizes[d];
      /* compute all linear indices in range */
      for (unsigned long n = 0; n < n_inds; n++) {
         /* store current linear index */
         inds[n] = index;
         /* compute next linear index */
         unsigned long d = n_dims_i - 1;
         index -= (range[d])[pos[d]] * d_sizes[d];
         pos[d]++;
         while ((pos[d] == range[d].size()) && (d > 0)) {
            pos[d] = 0;
            index += (range[d])[0] * d_sizes[d];
            d--;
            index -= (range[d])[pos[d]] * d_sizes[d];
            pos[d]++;
         }
         if (pos[d] < range[d].size())
            index += (range[d])[pos[d]] * d_sizes[d];
      }
   }
   return inds;
}

/*
 * Computing relative linear offsets of neighbors from the center of a
 * block of 3^d elements for a d-dimensional matrix of the given dimensions.
 *
 * Return an array of 3^d offsets (including the center of the block).
 * Note that some offsets will be negative.
 */
array<long> matrix_base::linear_neighbor_offsets(
   const array<unsigned long>& dims)
{
   /* compute indices of block of 3^d neighbors centered on (1,1,...,1) */
   unsigned long n_dims = dims.size();
   array<unsigned long> start(n_dims);
   array<unsigned long> end(n_dims, 2);
   array<unsigned long> inds = matrix_base::linear_indices(dims, start, end);
   /* compute index of (1,1,...,1) */
   array<unsigned long> center_coord(n_dims, 1);
   long center_ind = static_cast<long>(
      matrix_base::linear_index(dims, center_coord)
   );
   /* compute offsets of neighbors from center of block */
   unsigned long n_inds = inds.size();
   array<long> offsets(n_inds);
   for (unsigned long n = 0; n < n_inds; n++)
      offsets[n] = static_cast<long>(inds[n]) - center_ind;
   return offsets;
}

} /* namespace matrices */
} /* namespace math */
