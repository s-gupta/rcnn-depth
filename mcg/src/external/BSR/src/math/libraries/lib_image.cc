/*
 * Image processing library.
 */
#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/queue_set.hh"
#include "collections/pointers/auto_collection.hh"
#include "concurrent/threads/child_thread.hh"
#include "concurrent/threads/runnable.hh"
#include "concurrent/threads/thread.hh"
#include "functors/comparable_functors.hh"
#include "functors/distanceable_functors.hh"
#include "lang/array.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/exceptions/ex_not_found.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/seg_intersect.hh"
#include "math/geometry/triangulation.hh"
#include "math/libraries/lib_image.hh"
#include "math/libraries/lib_signal.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "math/matrices/cmatrix.hh"
#include "mlearning/clustering/clusterers/abstract/centroid_clusterer.hh"
#include "mlearning/clustering/clusterers/general/sample_clusterer.hh"
#include "mlearning/clustering/clusterers/kmeans/matrix_clusterer.hh"
#include "mlearning/clustering/metrics/matrix_metrics.hh"

namespace math {
namespace libraries {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::queue_set;
using collections::pointers::auto_collection;
using concurrent::threads::child_thread;
using concurrent::threads::runnable;
using concurrent::threads::thread;
using functors::comparable_functor;
using functors::compare_functor;
using functors::compare_functors;
using functors::distanceable_functor;
using lang::array;
using lang::exceptions::ex_index_out_of_bounds;
using lang::exceptions::ex_invalid_argument;
using lang::exceptions::ex_not_found;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::geometry::point_2D;
using math::geometry::seg_intersect;
using math::geometry::triangulation;
using math::matrices::matrix;
using math::matrices::cmatrix;
using mlearning::clustering::clusterers::abstract::centroid_clusterer;
using mlearning::clustering::clusterers::general::sample_clusterer;
using mlearning::clustering::metrics::matrix_metrics;
using namespace mlearning::clustering::clusterers;

/***************************************************************************
 * Image color space transforms.
 ***************************************************************************/

/*
 * Compute a grayscale image from an RGB image.
 */
matrix<> lib_image::grayscale(
   const matrix<>& r,
   const matrix<>& g,
   const matrix<>& b)
{
   /* check arguments */
   matrix<>::assert_dims_equal(r._dims, g._dims);
   matrix<>::assert_dims_equal(r._dims, b._dims);
   /* convert to grayscale */
   matrix<> m(r._dims);
   for (unsigned long n = 0; n < m._size; n++) {
      m._data[n] =
         (0.29894 * r._data[n])
       + (0.58704 * g._data[n])
       + (0.11402 * b._data[n]);
   }
   return m;
}

/*
 * Normalize a grayscale image so that intensity values lie in [0,1].
 */
void lib_image::grayscale_normalize(matrix<>& m) {
   if (m._size > 0) {
      double max_val = max(m);
      if (max_val > 0) {
         for (unsigned long n = 0; n < m._size; n++)
            m._data[n] /= max_val;
      }
   }
}

/*
 * Normalize and stretch a grayscale image so that intensity values span 
 * the full [0,1] range.
 */
void lib_image::grayscale_normalize_stretch(matrix<>& m) {
   /* subtract off minimum value */
   if (m._size > 0) {
      double min_val = min(m);
      for (unsigned long n = 0; n < m._size; n++)
         m._data[n] -= min_val;
   }
   /* normalize image */
   lib_image::grayscale_normalize(m);
}

/*
 * Gamma correct the RGB image using the given correction value.
 */
void lib_image::rgb_gamma_correct(
   matrix<>& r,
   matrix<>& g,
   matrix<>& b,
   double gamma)
{
   /* check arguments */
   matrix<>::assert_dims_equal(r._dims, g._dims);
   matrix<>::assert_dims_equal(r._dims, b._dims);
   /* gamma correct image */
   for (unsigned long n = 0; n < r._size; n++) {
      r._data[n] = math::pow(r._data[n], gamma);
      g._data[n] = math::pow(g._data[n], gamma);
      b._data[n] = math::pow(b._data[n], gamma);
   }
}

/*
 * Normalize an Lab image so that values for each channel lie in [0,1].
 */
void lib_image::lab_normalize(
   matrix<>& l,
   matrix<>& a,
   matrix<>& b)
{
   /* check arguments */
   matrix<>::assert_dims_equal(l._dims, a._dims);
   matrix<>::assert_dims_equal(l._dims, b._dims);
   /* range for a, b channels */
   const double ab_min = -73;
   const double ab_max = 95;
   const double ab_range = ab_max - ab_min;
   /* normalize Lab image */
   for (unsigned long n = 0; n < l._size; n++) {
      double l_val = l._data[n] / 100.0;
      double a_val = (a._data[n] - ab_min) / ab_range;
      double b_val = (b._data[n] - ab_min) / ab_range;
      if (l_val < 0) { l_val = 0; } else if (l_val > 1) { l_val = 1; }
      if (a_val < 0) { a_val = 0; } else if (a_val > 1) { a_val = 1; }
      if (b_val < 0) { b_val = 0; } else if (b_val > 1) { b_val = 1; }
      l._data[n] = l_val;
      a._data[n] = a_val;
      b._data[n] = b_val;
   }
}

/*
 * Convert from RGB color space to XYZ color space.
 */
void lib_image::rgb_to_xyz(matrix<>& r_x, matrix<>& g_y, matrix<>& b_z) {
   /* check arguments */
   matrix<>::assert_dims_equal(r_x._dims, g_y._dims);
   matrix<>::assert_dims_equal(r_x._dims, b_z._dims);
   /* convert RGB to XYZ */
   for (unsigned long n = 0; n < r_x._size; n++) {
      double r = r_x._data[n];
      double g = g_y._data[n];
      double b = b_z._data[n];
      r_x._data[n] = (0.412453 * r) +  (0.357580 * g) + (0.180423 * b);
      g_y._data[n] = (0.212671 * r) +  (0.715160 * g) + (0.072169 * b);
      b_z._data[n] = (0.019334 * r) +  (0.119193 * g) + (0.950227 * b);
   }
}

/*
 * Convert from RGB color space to Lab color space.
 */
void lib_image::rgb_to_lab(matrix<>& r_l, matrix<>& g_a, matrix<>& b_b) {
   lib_image::rgb_to_xyz(r_l, g_a, b_b);
   lib_image::xyz_to_lab(r_l, g_a, b_b);
}

/*
 * Convert from XYZ color space to RGB color space.
 */
void lib_image::xyz_to_rgb(matrix<>& x_r, matrix<>& y_g, matrix<>& z_b) {
   /* check arguments */
   matrix<>::assert_dims_equal(x_r._dims, y_g._dims);
   matrix<>::assert_dims_equal(x_r._dims, z_b._dims);
   /* convert XYZ to RGB */
   for (unsigned long n = 0; n < x_r._size; n++) {
      double x = x_r._data[n];
      double y = y_g._data[n];
      double z = z_b._data[n];
      x_r._data[n] = ( 3.240479 * x) +  (-1.537150 * y) + (-0.498535 * z);
      y_g._data[n] = (-0.969256 * x) +  ( 1.875992 * y) + ( 0.041556 * z);
      z_b._data[n] = ( 0.055648 * x) +  (-0.204043 * y) + ( 1.057311 * z);
   }
}

/*
 * Convert from XYZ color space to Lab color space.
 */
void lib_image::xyz_to_lab(matrix<>& x_l, matrix<>& y_a, matrix<>& z_b) {
   /* check arguments */
   matrix<>::assert_dims_equal(x_l._dims, y_a._dims);
   matrix<>::assert_dims_equal(x_l._dims, z_b._dims);
   /* D65 white point reference */
   const double x_ref = 0.950456;
   const double y_ref = 1.000000;
   const double z_ref = 1.088754;
   /* threshold value  */
   const double threshold = 0.008856;
   /* convert XYZ to Lab */
   for (unsigned long n = 0; n < x_l._size; n++) {
      /* extract xyz color value and normalize by reference point */
      double x = (x_l._data[n] / x_ref);
      double y = (y_a._data[n] / y_ref);
      double z = (z_b._data[n] / z_ref);
      /* compute fx, fy, fz */
      double fx =
         (x > threshold) ? math::pow(x,(1.0/3.0)) : (7.787*x + (16.0/116.0));
      double fy =
         (y > threshold) ? math::pow(y,(1.0/3.0)) : (7.787*y + (16.0/116.0));
      double fz =
         (z > threshold) ? math::pow(z,(1.0/3.0)) : (7.787*z + (16.0/116.0));
      /* compute Lab color value */
      x_l._data[n] =
         (y > threshold) ? (116*math::pow(y,(1.0/3.0)) - 16) : (903.3*y);
      y_a._data[n] = 500 * (fx - fy);
      z_b._data[n] = 200 * (fy - fz);
   }
}

/*
 * Convert from Lab color space to RGB color space.
 */
void lib_image::lab_to_rgb(matrix<>& l_r, matrix<>& a_g, matrix<>& b_b) {
   lib_image::lab_to_xyz(l_r, a_g, b_b);
   lib_image::xyz_to_rgb(l_r, a_g, b_b);
}

/*
 * Convert from Lab color space to XYZ color space.
 */
void lib_image::lab_to_xyz(matrix<>& l_x, matrix<>& a_y, matrix<>& b_z) {
   /* check arguments */
   matrix<>::assert_dims_equal(l_x._dims, a_y._dims);
   matrix<>::assert_dims_equal(l_x._dims, b_z._dims);
   /* D65 white point reference */
   const double x_ref = 0.950456;
   const double y_ref = 1.000000;
   const double z_ref = 1.088754;
   /* threshold value  */
   const double threshold = 0.008856;
   /* convert Lab to XYZ */
   for (unsigned long n = 0; n < l_x._size; n++) {
      /* extract Lab color value */
      double l = l_x._data[n];
      double a = a_y._data[n];
      double b = b_z._data[n];
      /* compute y (normalized by reference point) */
      double y = (l + 16.0)/116.0;
      y = y*y*y;
      if (y <= threshold)
         y = l / 903.3;
      /* compute fy, fx, fz */
      double fy =
         (y > threshold) ? math::pow(y,(1.0/3.0)) : (7.787*y + (16.0/116.0));
      double fx = (a/500) + fy;
      double fz = fy - (b/200);
      /* compute x (normalized by reference point) */
      double x = fx*fx*fx;
      if (x <= threshold)
         x = (fx - (16.0/116.0))/7.787;
      /* compute z (normalized by reference point) */
      double z = fz*fz*fz;
      if (z <= threshold)
         z = (fz - (16.0/116.0))/7.787;
      /* compute xyz color value */
      l_x._data[n] = x * x_ref;
      a_y._data[n] = y * y_ref;
      b_z._data[n] = z * z_ref;
   }
}

/***************************************************************************
 * Matrix border operations (2D).
 ***************************************************************************/

namespace {
/*
 * Compute border mirrored matrix.
 */
template <typename T>
void compute_border_mirror_2D(
   const T*      m_src,       /* source matrix */
   T*            m_dst,       /* destination matrix */
   unsigned long border_x,    /* size of border */
   unsigned long border_y,
   unsigned long size_x_src,  /* size of source */
   unsigned long size_y_src)
{
   /* compute destination size */
   unsigned long size_x_dst = size_x_src + 2*border_x;
   unsigned long size_y_dst = size_y_src + 2*border_y;
   /* compute step sizes in destination matrix (for copying interior) */
   unsigned long ind_init = border_x * size_y_dst + border_y;
   unsigned long ind_step = 2*border_y;
   /* copy interior */
   for (unsigned long n = 0, ind = ind_init, x = 0; x < size_x_src; x++) {
      for (unsigned long y = 0; y < size_y_src; y++, n++, ind++)
         m_dst[ind] = m_src[n];
      ind += ind_step;
   }
   /* mirror top */
   {
      unsigned long ind_intr = border_x * size_y_dst;
      unsigned long ind_bdr = ind_intr + size_y_dst;
      for (unsigned long x = 0; x < border_x; x++) {
         ind_bdr -= 2*size_y_dst;
         for (unsigned long y = 0; y < size_y_dst; y++)
            m_dst[ind_bdr++] = m_dst[ind_intr++];
      }
   }
   /* mirror bottom */
   {
      unsigned long ind_bdr = (size_x_src + border_x) * size_y_dst;
      unsigned long ind_intr = ind_bdr + size_y_dst;
      for (unsigned long x = 0; x < border_x; x++) {
         ind_intr -= 2*size_y_dst;
         for (unsigned long y = 0; y < size_y_dst; y++)
            m_dst[ind_bdr++] = m_dst[ind_intr++];
      }
   }
   /* mirror left */
   {
      unsigned long ind_bdr = 0;
      unsigned long ind_intr = 2*border_y;
      for (unsigned long x = 0; x < size_x_dst; x++) {
         for (unsigned long y = 0; y < border_y; y++)
            m_dst[ind_bdr++] = m_dst[--ind_intr];
         ind_bdr  += (size_y_dst - border_y);
         ind_intr += (size_y_dst + border_y);
      }
   }
   /* mirror right */
   {
      unsigned long ind_bdr  = size_y_dst - border_y;
      unsigned long ind_intr = ind_bdr;
      for (unsigned long x = 0; x < size_x_dst; x++) {
         for (unsigned long y = 0; y < border_y; y++)
            m_dst[ind_bdr++] = m_dst[--ind_intr];
         ind_bdr  += (size_y_dst - border_y);
         ind_intr += (size_y_dst + border_y);
      }
   }
}

/*
 * Compute border trimmed matrix.
 */
template <typename T>
void compute_border_trim_2D(
   const T*      m_src,       /* source matrix */
   T*            m_dst,       /* destination matrix */
   unsigned long border_x,    /* size of border */
   unsigned long border_y,
   unsigned long size_x_dst,  /* size of destination */
   unsigned long size_y_dst,
   unsigned long size_y_src)  /* size of source in y-dimension */
{
   /* compute step sizes in source matrix */
   unsigned long ind_init = border_x * size_y_src + border_y;
   unsigned long ind_step = 2*border_y;
   /* trim border */
   for (unsigned long n = 0, ind = ind_init, x = 0; x < size_x_dst; x++) {
      for (unsigned long y = 0; y < size_y_dst; y++, n++, ind++)
         m_dst[n] = m_src[ind];
      ind += ind_step;
   }
}
} /* namespace */

/*
 * Expand the border of the 2D matrix by the specified size on all sides.
 * The expanded border is filled with the mirror image of interior data.
 */
matrix<> lib_image::border_mirror_2D(const matrix<>& m, unsigned long size) {
   return lib_image::border_mirror_2D(m, size, size);
}

cmatrix<> lib_image::border_mirror_2D(const cmatrix<>& m, unsigned long size) {
   return lib_image::border_mirror_2D(m, size, size);
}

/*
 * Expand the border of the 2D matrix by the specified sizes in the x- and
 * y-dimensions.  The expanded border is filled with the mirror image of
 * interior data.
 */
matrix<> lib_image::border_mirror_2D(
   const matrix<>& m, unsigned long border_x, unsigned long border_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix size */
   unsigned long size_x_src = m._dims[0];
   unsigned long size_y_src = m._dims[1];
   /* check that interior can be mirrored */
   if ((border_x > size_x_src) || (border_y > size_y_src))
      throw ex_invalid_argument(
         "cannot create mirrored border larger than matrix interior dimensions"
      );
   /* mirror border */
   matrix<> m_dst(size_x_src + 2*border_x, size_y_src + 2*border_y);
   compute_border_mirror_2D(
      m._data, m_dst._data, border_x, border_y, size_x_src, size_y_src
   );
   return m_dst;
}

cmatrix<> lib_image::border_mirror_2D(
   const cmatrix<>& m, unsigned long border_x, unsigned long border_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix size */
   unsigned long size_x_src = m._dims[0];
   unsigned long size_y_src = m._dims[1];
   /* check that interior can be mirrored */
   if ((border_x > size_x_src) || (border_y > size_y_src))
      throw ex_invalid_argument(
         "cannot create mirrored border larger than matrix interior dimensions"
      );
   /* mirror border */
   cmatrix<> m_dst(size_x_src + 2*border_x, size_y_src + 2*border_y);
   compute_border_mirror_2D(
      m._data, m_dst._data, border_x, border_y, size_x_src, size_y_src
   );
   return m_dst;
}

/*
 * Trim the specified border size off all sides of the 2D matrix.
 */
matrix<> lib_image::border_trim_2D(const matrix<>& m, unsigned long size) {
   return lib_image::border_trim_2D(m, size, size);
}

cmatrix<> lib_image::border_trim_2D(const cmatrix<>& m, unsigned long size) {
   return lib_image::border_trim_2D(m, size, size);
}

/*
 * Trim the specified border dimensions off of the 2D matrix.
 */
matrix<> lib_image::border_trim_2D(
   const matrix<>& m, unsigned long border_x, unsigned long border_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix size */
   unsigned long size_x_src = m._dims[0];
   unsigned long size_y_src = m._dims[1];
   /* compute size of result matrix */
   unsigned long size_x_dst =
      (size_x_src > (2*border_x)) ? (size_x_src - 2*border_x) : 0;
   unsigned long size_y_dst =
      (size_y_src > (2*border_y)) ? (size_y_src - 2*border_y) : 0;
   /* trim border */
   matrix<> m_dst(size_x_dst, size_y_dst);
   compute_border_trim_2D(
      m._data, m_dst._data,
      border_x, border_y, size_x_dst, size_y_dst, size_y_src
   );
   return m_dst;
}

cmatrix<> lib_image::border_trim_2D(
   const cmatrix<>& m, unsigned long border_x, unsigned long border_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix size */
   unsigned long size_x_src = m._dims[0];
   unsigned long size_y_src = m._dims[1];
   /* compute size of result matrix */
   unsigned long size_x_dst =
      (size_x_src > (2*border_x)) ? (size_x_src - 2*border_x) : 0;
   unsigned long size_y_dst =
      (size_y_src > (2*border_y)) ? (size_y_src - 2*border_y) : 0;
   /* trim border */
   cmatrix<> m_dst(size_x_dst, size_y_dst);
   compute_border_trim_2D(
      m._data, m_dst._data,
      border_x, border_y, size_x_dst, size_y_dst, size_y_src
   );
   return m_dst;
}

/***************************************************************************
 * Matrix subdivision (2D).
 ***************************************************************************/

/*
 * Compute the optimal number of subdivisions along each dimension.
 *
 * The optimal subdivision parameters split the 2D matrix into no more than
 * the desired number of cells, while at the same time preferring nearly 
 * square cells over elongated rectangles.
 *
 * Adjacent cells overlap by the specified number of elements in each 
 * dimension.
 */
void lib_image::subdivision_params_2D(
   unsigned long  size_x,      unsigned long  size_y,
   unsigned long  n_subdiv,
   unsigned long& n_subdiv_x,  unsigned long& n_subdiv_y,
   unsigned long  overlap_x,   unsigned long  overlap_y)
{
   /* check number of subdivisions is valid */
   if (n_subdiv == 0)
      throw ex_invalid_argument("number of subdivisions must be nonzero");
   /* check matrix is not smaller than overlap */
   if (size_x < (2*overlap_x))
      throw ex_invalid_argument("x-overlap is too large");
   if (size_y < (2*overlap_y))
      throw ex_invalid_argument("y-overlap is too large");
   /* handle the case of an empty matrix */
   if ((size_x == 0) || (size_y == 0)) {
      n_subdiv_x = 1;
      n_subdiv_y = 1;
      return;
   }
   /* compute maximum number of subdivisions allowed */
   unsigned long max_subdiv_x =
      (overlap_x == 0) ? size_x : (size_x - overlap_x) / overlap_x;
   unsigned long max_subdiv_y =
      (overlap_y == 0) ? size_y : (size_y - overlap_y) / overlap_y;
   /* compute desired number of subdivisions */
   double length_x = static_cast<double>(size_x - overlap_x);
   double length_y = static_cast<double>(size_y - overlap_y);
   double num_x = math::sqrt(static_cast<double>(n_subdiv)*length_x/length_y);
   double num_y = math::sqrt(static_cast<double>(n_subdiv)*length_y/length_x);
   unsigned long nx = static_cast<unsigned long>(math::round(num_x));
   unsigned long ny = static_cast<unsigned long>(math::round(num_y));
   /* adjust to allowed subdivision limits */
   if (nx == 0) { nx = 1; } else if (nx > max_subdiv_x) { nx = max_subdiv_x; }
   if (ny == 0) { ny = 1; } else if (ny > max_subdiv_y) { ny = max_subdiv_y; }
   /* adjust to subdivision geometry */
   if ((nx > ny) || ((nx == ny) && (max_subdiv_x > max_subdiv_y))) {
      /* adjust nx, then ny */
      nx = ((n_subdiv / ny) < max_subdiv_x) ? (n_subdiv / ny) : max_subdiv_x;
      ny = ((n_subdiv / nx) < max_subdiv_y) ? (n_subdiv / nx) : max_subdiv_y;
   } else {
      /* adjust ny, then nx */
      ny = ((n_subdiv / nx) < max_subdiv_y) ? (n_subdiv / nx) : max_subdiv_y;
      nx = ((n_subdiv / ny) < max_subdiv_x) ? (n_subdiv / ny) : max_subdiv_x;
   }
   /* return subdivision parameters */
   n_subdiv_x = nx;
   n_subdiv_y = ny;
}

void lib_image::subdivision_params_2D(
   const matrix<>& m,
   unsigned long  n_subdiv,
   unsigned long& n_subdiv_x,  unsigned long& n_subdiv_y,
   unsigned long  overlap_x,   unsigned long  overlap_y)
{
   if (m.dimensionality() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   return lib_image::subdivision_params_2D(
      m.size(0), m.size(1),
      n_subdiv, n_subdiv_x, n_subdiv_y, overlap_x, overlap_y
   );
}

void lib_image::subdivision_params_2D(
   const cmatrix<>& m,
   unsigned long  n_subdiv,
   unsigned long& n_subdiv_x,  unsigned long& n_subdiv_y,
   unsigned long  overlap_x,   unsigned long  overlap_y)
{
   if (m.dimensionality() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   return lib_image::subdivision_params_2D(
      m.size(0), m.size(1),
      n_subdiv, n_subdiv_x, n_subdiv_y, overlap_x, overlap_y
   );
}

namespace {
/*
 * Compute matrix subdivision.
 */
template <typename T>
void compute_subdivide_2D(
   const matrix<T>&         m,            /* matrix */
   unsigned long            n_subdiv_x,   /* # of subdivisions in x-dimension */
   unsigned long            n_subdiv_y,   /* # of subdivisions in y-dimension */
   unsigned long            overlap_x,    /* overlap in x-dimension */
   unsigned long            overlap_y,    /* overlap in y-dimension */
   array_list< matrix<T> >& submatrices)  /* returned subdivision */
{
   /* check that matrix is 2D */
   if (m.dimensionality() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* check number of subdivisions is valid */
   if ((n_subdiv_x == 0) || (n_subdiv_y == 0))
      throw ex_invalid_argument("number of subdivisions must be nonzero");
   /* get matrix dimensions */
   unsigned long size_x = m.size(0);
   unsigned long size_y = m.size(1);
   /* compute total overlap size */
   unsigned long overlap_x_total = (n_subdiv_x + 1) * overlap_x;
   unsigned long overlap_y_total = (n_subdiv_y + 1) * overlap_y;
   /* check that all subdivisions can be filled */
   if ((overlap_x_total > size_x) || (n_subdiv_x > size_x))
      throw ex_invalid_argument("cannot fill requested # of x-subdivisions");
   if ((overlap_y_total > size_y) || (n_subdiv_y > size_y))
      throw ex_invalid_argument("cannot fill requested # of y-subdivisions");
   /* compute subdivision size (excluding overlap) */
   unsigned long interior_x_total = size_x - overlap_x_total;
   unsigned long interior_y_total = size_y - overlap_y_total;
   unsigned long interior_x = interior_x_total / n_subdiv_x;
   unsigned long interior_y = interior_y_total / n_subdiv_y;
   unsigned long extra_x = interior_x_total - (interior_x * n_subdiv_x);
   unsigned long extra_y = interior_y_total - (interior_y * n_subdiv_y);
   /* subdivide matrix */
   const T* m_data = m.data();
   for (unsigned long nx = 0, n = 0, start_x = 0; nx < n_subdiv_x; nx++) {
      /* compute submatrix x-size */
      unsigned long sx_interior = interior_x + ((nx < extra_x) ? 1 : 0);
      unsigned long sx = sx_interior + (2*overlap_x);
      for (unsigned long ny = 0, start_y = 0; ny < n_subdiv_y; ny++, n++) {
         /* compute submatrix y-size */
         unsigned long sy_interior = interior_y + ((ny < extra_y) ? 1 : 0);
         unsigned long sy = sy_interior + (2*overlap_y);
         /* copy submatrix */
         matrix<T>& m_sub = submatrices[n];
         m_sub.resize(sx, sy);
         T* m_sub_data = m_sub.data();
         unsigned long pos = start_x * size_y + start_y;
         for (unsigned long x = 0, ind = 0; x < sx; x++) {
            for (unsigned long y = 0; y < sy; y++)
               m_sub_data[ind++] = m_data[pos++];
            pos += (size_y - sy);
         }
         /* update y offset */
         start_y += (sy_interior + overlap_y);
      }
      /* update x offset */
      start_x += (sx_interior + overlap_x);
   }
}

/*
 * Compute combined matrix from subdivision.
 */
template <typename T>
void compute_combine_2D(
   const array_list< matrix<T> >& submatrices,  /* subdivision */
   unsigned long                  n_subdiv_x,   /* # of x-subdivisions */
   unsigned long                  n_subdiv_y,   /* # of y-subdivisions */
   unsigned long                  overlap_x,    /* overlap in x-dimension */
   unsigned long                  overlap_y,    /* overlap in y-dimension */
   matrix<T>&                     m)            /* returned matrix */
{
   /* check for correct number of subdivisions */
   unsigned long n_subdiv = n_subdiv_x * n_subdiv_y;
   if (submatrices.size() != n_subdiv)
      throw ex_invalid_argument("incorrect number of subdivisions specified");
   /* check that all submatrices are 2D */
   for (unsigned long n = 0; n < n_subdiv; n++) {
      if (submatrices[n].dimensionality() != 2)
         throw ex_invalid_argument("all submatrices must be 2D");
   }
   /* determine x-subdivision sizes */
   array<unsigned long> subdiv_sizes_x(n_subdiv_x);
   unsigned long size_x = overlap_x;
   for (unsigned long nx = 0, n = 0; nx < n_subdiv_x; nx++, n += n_subdiv_y) {
      unsigned long sx = submatrices[n].size(0);
      if (sx < (2*overlap_x))
         throw ex_invalid_argument("incorrect submatrix size");
      sx -= (2*overlap_x);
      size_x += (sx + overlap_x);
      subdiv_sizes_x[nx] = sx;
   }
   /* determine y-subdivision sizes */
   array<unsigned long> subdiv_sizes_y(n_subdiv_y);
   unsigned long size_y = overlap_y;
   for (unsigned long ny = 0; ny < n_subdiv_y; ny++) {
      unsigned long sy = submatrices[ny].size(1);
      if (sy < (2*overlap_y))
         throw ex_invalid_argument("incorrect submatrix size");
      sy -= (2*overlap_y);
      size_y += (sy + overlap_y);
      subdiv_sizes_y[ny] = sy;
   }     
   /* combine submatrices */
   m.resize(size_x, size_y);
   T* m_data = m.data();
   for (unsigned long nx = 0, n = 0, start_x = 0; nx < n_subdiv_x; nx++) {
      /* get submatrix x-size */
      unsigned long sx_interior = subdiv_sizes_x[nx];
      unsigned long sx = sx_interior + (2*overlap_x);
      for (unsigned long ny = 0, start_y = 0; ny < n_subdiv_y; ny++, n++) {
         /* get submatrix y-size */
         unsigned long sy_interior = subdiv_sizes_y[ny];
         unsigned long sy = sy_interior + (2*overlap_y);
         /* check submatrix size */
         const matrix<T>& m_sub = submatrices[n];
         if ((m_sub.size(0) != sx) || (m_sub.size(1) != sy))
            throw ex_invalid_argument("submatrix size disagreement");
         /* copy submatrix */
         const T* m_sub_data = m_sub.data();
         unsigned long pos = start_x * size_y + start_y;
         for (unsigned long x = 0, ind = 0; x < sx; x++) {
            for (unsigned long y = 0; y < sy; y++)
               m_data[pos++] = m_sub_data[ind++];
            pos += (size_y - sy);
         }
         /* update y offset */
         start_y += (sy_interior + overlap_y);
      }
      /* update x offset */
      start_x += (sx_interior + overlap_x);
   }
}
} /* namespace */

/*
 * Subdivide the 2D matrix into a grid of overlapping submatrices.
 * Return the collection of submatrices.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::subdivide_2D(
   const matrix<>& m,     
   unsigned long n_subdiv_x, unsigned long n_subdiv_y,
   unsigned long overlap_x,  unsigned long overlap_y)
{
   /* allocate subdivision */
   unsigned long n_subdiv = n_subdiv_x * n_subdiv_y;
   auto_collection< matrix<>, array_list< matrix<> > > submatrices(
      new array_list< matrix<> >()
   );
   for (unsigned long n = 0; n < n_subdiv; n++) {
      auto_ptr< matrix<> > m_sub(new matrix<>());
      submatrices->add(*m_sub);
      m_sub.release();
   }
   /* compute subdivision */
   compute_subdivide_2D(
      m, n_subdiv_x, n_subdiv_y, overlap_x, overlap_y, *submatrices
   );
   return submatrices;
}

auto_collection< cmatrix<>, array_list< cmatrix<> > > lib_image::subdivide_2D(
   const cmatrix<>& m,
   unsigned long n_subdiv_x, unsigned long n_subdiv_y,
   unsigned long overlap_x,  unsigned long overlap_y)
{
   /* allocate subdivision */
   unsigned long n_subdiv = n_subdiv_x * n_subdiv_y;
   auto_collection< cmatrix<>, array_list< cmatrix<> > > submatrices(
      new array_list< cmatrix<> >()
   );
   for (unsigned long n = 0; n < n_subdiv; n++) {
      auto_ptr< cmatrix<> > m_sub(new cmatrix<>());
      submatrices->add(*m_sub);
      m_sub.release();
   }
   /* create submatrix array of matrix< complex<> > type */
   array_list< matrix< complex<> > > submatrix_arr;
   for (unsigned long n = 0; n < n_subdiv; n++)
      submatrix_arr.add((*submatrices)[n]);
   /* compute subdivision */
   compute_subdivide_2D(
      m, n_subdiv_x, n_subdiv_y, overlap_x, overlap_y, submatrix_arr
   );
   return submatrices;
}

/*
 * Combine submatrices resulting from matrix subdivision.
 */
matrix<> lib_image::combine_2D(
   const array_list< matrix<> >& submatrices,
   unsigned long n_subdiv_x, unsigned long n_subdiv_y,
   unsigned long overlap_x,  unsigned long overlap_y)
{
   matrix<> m;
   compute_combine_2D(
      submatrices, n_subdiv_x, n_subdiv_y, overlap_x, overlap_y, m
   );
   return m;
}

cmatrix<> lib_image::combine_2D(
   const array_list< cmatrix<> >& submatrices,
   unsigned long n_subdiv_x, unsigned long n_subdiv_y,
   unsigned long overlap_x,  unsigned long overlap_y)
{
   /* create submatrix array of matrix< complex<> > type */
   array_list< matrix< complex<> > > submatrix_arr;
   for (unsigned long n = 0, n_subdiv = submatrices.size(); n < n_subdiv; n++)
      submatrix_arr.add(submatrices[n]);
   /* combine submatrices */
   cmatrix<> m;
   compute_combine_2D(
      submatrix_arr, n_subdiv_x, n_subdiv_y, overlap_x, overlap_y, m
   );
   return m;
}

/***************************************************************************
 * Matrix resampling (2D).
 ***************************************************************************/

namespace {
/*
 * Compute resampled 2D matrix using bilinear interpolation.
 */
template <typename T>
void compute_resample_2D(
   const T*      m_src,       /* source matrix */
   T*            m_dst,       /* destination matrix */
   unsigned long size_x_src,  /* size of source */
   unsigned long size_y_src,
   unsigned long size_x_dst,  /* size of destination */
   unsigned long size_y_dst)
{
   /* check that matrices are nonempty */
   if ((size_x_src > 0) && (size_y_src > 0) &&
       (size_x_dst > 0) && (size_y_dst > 0))
   {
      /* compute whether to anchor x, y extrema to existing extrema */
      const bool anchor_x =
         ((size_x_dst >= 3) || ((size_x_dst == 2) && (size_x_src <= 2)));
      const bool anchor_y =
         ((size_y_dst >= 3) || ((size_y_dst == 2) && (size_y_src <= 2)));
      /* compute step sizes */
      const double step_x =
         (anchor_x)
       ? (static_cast<double>(size_x_src-1)/static_cast<double>(size_x_dst-1))
       : (static_cast<double>(size_x_src-1)/static_cast<double>(size_x_dst+1));
      const double step_y =
         (anchor_y)
       ? (static_cast<double>(size_y_src-1)/static_cast<double>(size_y_dst-1))
       : (static_cast<double>(size_y_src-1)/static_cast<double>(size_y_dst+1));
      /* resample */
      double x = (anchor_x) ? 0 : step_x;
      unsigned long n = 0;
      for (unsigned long dst_x = 0; dst_x < size_x_dst; dst_x++) {
         /* compute integer coordinates bounding x-coordinate */
         unsigned long x0 = static_cast<unsigned long>(math::floor(x));
         unsigned long x1 = static_cast<unsigned long>(math::ceil(x));
         if (x1 == size_x_src) { x1--; }
         /* compute distances to x-bounds */
         const double dist_x0 = x - x0;
         const double dist_x1 = x1 - x;
         /* loop over y */
         double y = (anchor_y) ? 0 : step_y;
         for (unsigned long dst_y = 0; dst_y < size_y_dst; dst_y++) {
            /* compute integer coordinates bounding y-coordinate */
            unsigned long y0 = static_cast<unsigned long>(math::floor(y));
            unsigned long y1 = static_cast<unsigned long>(math::ceil(y));
            /* compute distances to y-bounds */
            if (y1 == size_y_src) { y1--; }
            const double dist_y0 = y - y0;
            const double dist_y1 = y1 - y;
            /* grab matrix elements */
            const T& m00 = m_src[x0*size_y_src + y0];
            const T& m01 = m_src[x0*size_y_src + y1];
            const T& m10 = m_src[x1*size_y_src + y0];
            const T& m11 = m_src[x1*size_y_src + y1];
            /* interpolate in x-direction */
            const T t0 = (x0 != x1) ? (dist_x1 * m00 + dist_x0 * m10) : m00;
            const T t1 = (x0 != x1) ? (dist_x1 * m01 + dist_x0 * m11) : m01;
            /* interpolate in y-direction */
            m_dst[n] = (y0 != y1) ? (dist_y1 * t0 + dist_y0 * t1) : t0;
            /* increment coordinate */
            n++;
            y += step_y;
         }
         x += step_x;
      }
   }
}  
} /* namespace */

/*
 * Resample the 2D matrix by the given scale factor in both the x- and y-
 * directions.  A factor less than 1 produces a result smaller than the
 * original matrix, while a factor greater than 1 produces a result larger.
 */
matrix<> lib_image::resample_2D(const matrix<>& m, double scale) {
   return lib_image::resample_2D(m, scale, scale);
}

cmatrix<> lib_image::resample_2D(const cmatrix<>& m, double scale) {
   return lib_image::resample_2D(m, scale, scale);
}

/*
 * Resample the 2D matrix by the given scale factors in the x- and y-
 * directions, respectively.
 */
matrix<> lib_image::resample_2D(
   const matrix<>& m, double scale_x, double scale_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute desired size from scale factor */
   unsigned long size_x = static_cast<unsigned long>(
      math::ceil(scale_x * static_cast<double>(m._dims[0]))
   );
   unsigned long size_y = static_cast<unsigned long>(
      math::ceil(scale_y * static_cast<double>(m._dims[1]))
   );
   /* resample */
   return lib_image::resample_2D(m, size_x, size_y);
}

cmatrix<> lib_image::resample_2D(
   const cmatrix<>& m, double scale_x, double scale_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute desired size from scale factor */
   unsigned long size_x = static_cast<unsigned long>(
      math::ceil(scale_x * static_cast<double>(m._dims[0]))
   );
   unsigned long size_y = static_cast<unsigned long>(
      math::ceil(scale_y * static_cast<double>(m._dims[1]))
   );
   /* resample */
   return lib_image::resample_2D(m, size_x, size_y);
}

/*
 * Resample the 2D matrix on a grid of the specified size.
 */
matrix<> lib_image::resample_2D(
   const matrix<>& m, unsigned long size_x, unsigned long size_y)
{  
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* resample */
   matrix<> m_dst(size_x, size_y);
   compute_resample_2D(
      m._data, m_dst._data, m._dims[0], m._dims[1], size_x, size_y
   );
   return m_dst;
}

cmatrix<> lib_image::resample_2D(
   const cmatrix<>& m, unsigned long size_x, unsigned long size_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* resample */
   cmatrix<> m_dst(size_x, size_y);
   compute_resample_2D(
      m._data, m_dst._data, m._dims[0], m._dims[1], size_x, size_y
   );
   return m_dst;
}

/***************************************************************************
 * Matrix rotation (2D).
 ***************************************************************************/

namespace {
/*
 * Compute x/y-support for rotated 2D matrix.
 */
double support_x_rotated(double support_x, double support_y, double ori) {
   const double sx_cos_ori = support_x * math::cos(ori);
   const double sy_sin_ori = support_y * math::sin(ori);
   double x0_mag = math::abs(sx_cos_ori - sy_sin_ori);
   double x1_mag = math::abs(sx_cos_ori + sy_sin_ori);
   return (x0_mag > x1_mag) ? x0_mag : x1_mag;
}

double support_y_rotated(double support_x, double support_y, double ori) {
   const double sx_sin_ori = support_x * math::sin(ori);
   const double sy_cos_ori = support_y * math::cos(ori);
   double y0_mag = math::abs(sx_sin_ori - sy_cos_ori);
   double y1_mag = math::abs(sx_sin_ori + sy_cos_ori);
   return (y0_mag > y1_mag) ? y0_mag : y1_mag;
}

/*
 * Compute integer-valued supports for rotated 2D matrix.
 */
unsigned long support_x_rotated(
   unsigned long support_x, unsigned long support_y, double ori)
{
   return static_cast<unsigned long>(
      math::ceil(support_x_rotated(
         static_cast<double>(support_x), static_cast<double>(support_y), ori
      ))
   );
}

unsigned long support_y_rotated(
   unsigned long support_x, unsigned long support_y, double ori)
{
   return static_cast<unsigned long>(
      math::ceil(support_y_rotated(
         static_cast<double>(support_x), static_cast<double>(support_y), ori
      ))
   );
}

/*
 * Compute x/y-size for rotated 2D matrix.
 */
unsigned long size_x_rotated(
   unsigned long size_x, unsigned long size_y, double ori)
{
   double support_x = static_cast<double>(size_x) / 2;
   double support_y = static_cast<double>(size_y) / 2;
   return static_cast<unsigned long>(
      math::ceil(support_x_rotated(support_x, support_y, ori) * 2)
   );
}

unsigned long size_y_rotated(
   unsigned long size_x, unsigned long size_y, double ori)
{
   double support_x = static_cast<double>(size_x) / 2;
   double support_y = static_cast<double>(size_y) / 2;
   return static_cast<unsigned long>(
      math::ceil(support_y_rotated(support_x, support_y, ori) * 2)
   );
}

/*
 * Compute rotated 2D matrix using bilinear interpolation.
 */
template <typename T>
void compute_rotate_2D(
   const T*      m_src,       /* source matrix */
   T*            m_dst,       /* destination matrix */
   unsigned long size_x_src,  /* size of source */
   unsigned long size_y_src,
   unsigned long size_x_dst,  /* size of destination */
   unsigned long size_y_dst,
   double ori)                /* orientation */
{
   /* check that matrices are nonempty */
   if ((size_x_src > 0) && (size_y_src > 0) &&
       (size_x_dst > 0) && (size_y_dst > 0))
   {
      /* compute sin and cos of rotation angle */
      const double cos_ori = math::cos(ori);
      const double sin_ori = math::sin(ori);
      /* compute location of origin in src */
      const double origin_x_src = static_cast<double>((size_x_src - 1)) / 2;
      const double origin_y_src = static_cast<double>((size_y_src - 1)) / 2;
      /* rotate */
      double u = -(static_cast<double>((size_x_dst - 1)) / 2);
      unsigned long n = 0;
      for (unsigned long dst_x = 0; dst_x < size_x_dst; dst_x++) {
         double v = -(static_cast<double>((size_y_dst - 1)) / 2);
         for (unsigned long dst_y = 0; dst_y < size_y_dst; dst_y++) {
            /* reverse rotate by orientation and shift by origin offset */
            double x = u * cos_ori + v * sin_ori + origin_x_src;
            double y = v * cos_ori - u * sin_ori + origin_y_src;
            /* check that location is in first quadrant */
            if ((x >= 0) && (y >= 0)) {
               /* compute integer bounds on location */
               unsigned long x0 = static_cast<unsigned long>(math::floor(x));
               unsigned long x1 = static_cast<unsigned long>(math::ceil(x));
               unsigned long y0 = static_cast<unsigned long>(math::floor(y));
               unsigned long y1 = static_cast<unsigned long>(math::ceil(y));
               /* check that location is within src matrix */
               /*if ((0 <= x0) && (x1 < size_x_src) &&
                   (0 <= y0) && (y1 < size_y_src))
               {*/
               /* 05/2014 - Modified by Jordi Pont-Tuset <jordi.pont@upc.edu> */
               if ((x1 < size_x_src) && (y1 < size_y_src))
               {
                  /* compute distances to bounds */
                  double dist_x0 = x - x0;
                  double dist_x1 = x1 - x;
                  double dist_y0 = y - y0;
                  double dist_y1 = y1 - y;
                  /* grab matrix elements */
                  const T& m00 = m_src[x0*size_y_src + y0];
                  const T& m01 = m_src[x0*size_y_src + y1];
                  const T& m10 = m_src[x1*size_y_src + y0];
                  const T& m11 = m_src[x1*size_y_src + y1];
                  /* interpolate in x-direction */
                  const T t0 =
                     (x0 != x1) ? (dist_x1 * m00 + dist_x0 * m10) : m00;
                  const T t1 =
                     (x0 != x1) ? (dist_x1 * m01 + dist_x0 * m11) : m01;
                  /* interpolate in y-direction */
                  m_dst[n] = (y0 != y1) ? (dist_y1 * t0 + dist_y0 * t1) : t0;
               }
            }
            /* increment coordinate */
            n++;
            v++;
         }
         u++;
      }
   }
}
} /* namespace */

/*
 * Rotate and pad with zeros.
 */
matrix<> lib_image::rotate_2D(const matrix<>& m, double ori) {
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute size of rotated matrix */
   unsigned long size_x = size_x_rotated(m._dims[0], m._dims[1], ori);
   unsigned long size_y = size_y_rotated(m._dims[0], m._dims[1], ori);
   /* compute rotation */
   return lib_image::rotate_2D_crop(m, ori, size_x, size_y);
}

cmatrix<> lib_image::rotate_2D(const cmatrix<>& m, double ori) {
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute size of rotated matrix */
   unsigned long size_x = size_x_rotated(m._dims[0], m._dims[1], ori);
   unsigned long size_y = size_y_rotated(m._dims[0], m._dims[1], ori);
   /* compute rotation */
   return lib_image::rotate_2D_crop(m, ori, size_x, size_y);
}

/*
 * Rotate and pad with zeros, but crop the result to be the same size as
 * the input matrix.
 */
matrix<> lib_image::rotate_2D_crop(const matrix<>& m, double ori) {
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute rotation */
   return lib_image::rotate_2D_crop(m, ori, m._dims[0], m._dims[1]);
}

cmatrix<> lib_image::rotate_2D_crop(const cmatrix<>& m, double ori) {
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute rotation */
   return lib_image::rotate_2D_crop(m, ori, m._dims[0], m._dims[1]);
}

/*
 * Rotate and pad with zeros, but crop the result to the specified size.
 */
matrix<> lib_image::rotate_2D_crop(
   const matrix<>& m, double ori, unsigned long size_x, unsigned long size_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute rotation */
   matrix<> m_rot(size_x, size_y);
   compute_rotate_2D(
      m._data, m_rot._data, m._dims[0], m._dims[1], size_x, size_y, ori
   );
   return m_rot;
}

cmatrix<> lib_image::rotate_2D_crop(
   const cmatrix<>& m, double ori, unsigned long size_x, unsigned long size_y)
{
   /* check that matrix is 2D */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* compute rotation */
   cmatrix<> m_rot(size_x, size_y);
   compute_rotate_2D(
      m._data, m_rot._data, m._dims[0], m._dims[1], size_x, size_y, ori
   );
   return m_rot;
}

/***************************************************************************
 * Image pyramids.
 ***************************************************************************/

/*
 * Gaussian pyramid.
 *
 * Construct the Gaussian pyramid for the given image.  The pyramid is
 * defined by:
 *
 *                   I[n] = D(G(sigma) * I[n-1])
 *
 * where * denotes convolution, D denotes downsampling by a factor of
 * 1/sigma in the x- and y-dimensions, I[0] is then input image, and 
 * G is a 2D Gaussian kernel.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::gaussian_pyramid(
   const matrix<>& m, double sigma)
{
   /* check arguments */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("image must be 2D");
   if (sigma <= 1)
      throw ex_invalid_argument("sigma must be > 1");
   /* compute gaussian filters */
   matrix<> gx = lib_image::gaussian(sigma);
   matrix<> gy = transpose(gx);
   unsigned long min_size = gx.size(0);
   /* allocate pyramid */
   auto_collection< matrix<>, array_list< matrix<> > > g_pyramid(
      new array_list< matrix<> >()
   );
   /* initialize base level of pyramid */
   auto_ptr< matrix<> > g_base(new matrix<>(m));
   g_pyramid->add(*g_base);
   g_base.release();
   /* compute gaussian pyramid */
   while (true) {
      /* get previous level */
      const matrix<>& g_prev = g_pyramid->tail();
      /* compute dimensions of current level */
      unsigned long size_x = static_cast<unsigned long>(
         math::floor(static_cast<double>(g_prev._dims[0]) / sigma)
      );
      unsigned long size_y = static_cast<unsigned long>(
         math::floor(static_cast<double>(g_prev._dims[1]) / sigma)
      );
      /* check if at top of pyramid */
      if ((size_x < min_size) || (size_y < min_size))
         break;
      /* compute current level */
      auto_ptr< matrix<> > g_curr(new matrix<>(
         lib_image::resample_2D(
            conv_crop(conv_crop(g_prev, gx), gy), size_x, size_y
         )
      ));
      /* store current level */
      g_pyramid->add(*g_curr);
      g_curr.release();
   }
   return g_pyramid;
}

/*
 * Laplacian pyramid.
 *
 * Construct an approximation to the scale-normalized Laplacian of Gaussian
 * pyramid by computing difference of Gaussians.  This pyramid is defined by:
 *
 *             L[n] = I[n] - U(I[n+1])    for 0 <= n < N
 *             L[N] = I[N]                for the coarsest level N
 *
 * where U denotes upsampling by a factor of sigma, and I[n] is a level in
 * the Gaussian pyramid defined above.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::laplacian_pyramid(
   const matrix<>& m, double sigma)
{
   auto_collection< matrix<>, array_list< matrix<> > > g_pyramid =
      lib_image::gaussian_pyramid(m, sigma);
   return lib_image::laplacian_pyramid(*g_pyramid);
}

auto_collection< matrix<>, array_list< matrix<> > > lib_image::laplacian_pyramid(
   const array_list< matrix<> >& g_pyramid)
{
   /* allocate pyramid */
   auto_collection< matrix<>, array_list< matrix<> > > l_pyramid(
      new array_list< matrix<> >()
   );
   /* compute pyramid */
   unsigned long n_levels = g_pyramid.size();
   for (unsigned long n = 0; (n+1) < n_levels; n++) {
      const matrix<>& g_curr = g_pyramid[n];
      const matrix<>& g_next = g_pyramid[n+1];
      auto_ptr< matrix<> > l_curr(new matrix<>(
         g_curr - 
         lib_image::resample_2D(g_next, g_curr._dims[0], g_curr._dims[1])
      ));
      l_pyramid->add(*l_curr);
      l_curr.release();
   }
   /* add coarsest level of gaussian pyramid to top of laplacian pyramid */
   if (n_levels > 0) {
      auto_ptr< matrix<> > l_top(new matrix<>(g_pyramid[n_levels-1]));
      l_pyramid->add(*l_top);
      l_top.release();
   }
   return l_pyramid;
}

/*
 * Resize all images in the pyramid to be the same size.
 * The base (largest) image is scaled by the given factor in both the
 * x- and y-dimensions and all other images are resampled to that size.
 */
void lib_image::pyramid_resize(array_list< matrix<> >& pyramid, double scale) {
   /* check that pyramid is nonempty */
   unsigned long n_levels = pyramid.size();
   if (n_levels == 0)
      return;
   /* resize base level */
   pyramid[0] = lib_image::resample_2D(pyramid[0], scale);
   /* get size of base level */
   unsigned long size_x = pyramid[0]._dims[0];
   unsigned long size_y = pyramid[0]._dims[1];
   /* resize other levels */
   for (unsigned long n = 1; n < n_levels; n++)
      pyramid[n] = lib_image::resample_2D(pyramid[n], size_x, size_y);
}

/*
 * Reconstruct the original image given its laplacian pyramid.
 */
matrix<> lib_image::image_from_laplacian_pyramid(
   const array_list< matrix<> >& l_pyramid)
{
   /* check that pyramid is nonempty */
   unsigned long n_levels = l_pyramid.size();
   if (n_levels == 0)
      return matrix<>();
   /* reconstruct original image */
   matrix<> m(l_pyramid[n_levels-1]);
   for (unsigned long n = (n_levels-1); n > 0; n--) {
      const matrix<>& l_curr = l_pyramid[n-1];
      m = lib_image::resample_2D(m, l_curr._dims[0], l_curr._dims[1]);
      m += l_curr;
   }
   return m;
}

/***************************************************************************
 * Gaussian kernels.
 ***************************************************************************/

/*
 * Gaussian kernel (1D).
 *
 * Specify the standard deviation and (optionally) the support.
 * The length of the returned vector is 2*support + 1.
 * The support defaults to 3*sigma.
 * The kernel is normalized to have unit L1 norm.
 * If returning a 1st or 2nd derivative, the kernel has zero mean.
 */
matrix<> lib_image::gaussian(
   double sigma, unsigned int deriv, bool hlbrt)
{
   unsigned long support = static_cast<unsigned long>(math::ceil(3*sigma));
   return lib_image::gaussian(sigma, deriv, hlbrt, support);
}

matrix<> lib_image::gaussian(
   double sigma, unsigned int deriv, bool hlbrt, unsigned long support)
{
   /* enlarge support so that hilbert transform can be done efficiently */
   unsigned long support_big = support;
   if (hlbrt) {
      support_big = 1;
      unsigned long temp = support;
      while (temp > 0) {
         support_big *= 2;
         temp /= 2;
      }
   }
   /* compute constants */
   const double sigma2_inv = double(1)/(sigma*sigma);
   const double neg_two_sigma2_inv = double(-0.5)*sigma2_inv;
   /* compute gaussian (or gaussian derivative) */
   unsigned long size = 2*support_big + 1;
   matrix<> m(size, 1);
   double x = -(static_cast<double>(support_big)); 
   if (deriv == 0) {
      /* compute guassian */
      for (unsigned long n = 0; n < size; n++, x++)
         m._data[n] = math::exp(x*x*neg_two_sigma2_inv);
   } else if (deriv == 1) {
      /* compute gaussian first derivative */
      for (unsigned long n = 0; n < size; n++, x++)
         m._data[n] = math::exp(x*x*neg_two_sigma2_inv) * (-x);
   } else if (deriv == 2) {
      /* compute gaussian second derivative */
      for (unsigned long n = 0; n < size; n++, x++) {
         double x2 = x * x;
         m._data[n] = math::exp(x2*neg_two_sigma2_inv) * (x2*sigma2_inv - 1);
      }
   } else {
      throw ex_invalid_argument("only derivatives 0, 1, 2 supported");
   }
   /* take hilbert transform (if requested) */
   if (hlbrt) {
      /* grab power of two sized submatrix (ignore last element) */
      m._size--;
      m._dims[0]--;
      /* grab desired submatrix after hilbert transform */
      array<unsigned long> start(2);
      array<unsigned long> end(2);
      start[0] = support_big - support;
      end[0] = start[0] + support + support;
      m = (lib_signal::hilbert(m)).submatrix(start, end);
   }
   /* make zero mean (if returning derivative) */
   if (deriv > 0)
      m -= mean(m);
   /* make unit L1 norm */
   m /= sum(abs(m));
   return m;
}

/*
 * Gaussian kernel (2D).
 *
 * Specify the standard deviation along each axis, and (optionally) the 
 * orientation (in radians).  In addition, optionally specify the support.
 * The support defaults to 3*sigma (automatically adjusted for rotation).
 *
 * The 1st or 2nd derivative and/or the hilbert transform can optionally
 * be taken along the y-axis prior to rotation.
 *
 * The kernel is normalized to have unit L1 norm.
 * If returning a 1st or 2nd derivative, the kernel has zero mean.
 */
matrix<> lib_image::gaussian_2D(
   double sigma_x, double sigma_y, double ori, unsigned int deriv, bool hlbrt)
{   
   /* compute support from rotated corners of original support box */
   unsigned long support_x = static_cast<unsigned long>(math::ceil(3*sigma_x));
   unsigned long support_y = static_cast<unsigned long>(math::ceil(3*sigma_y));
   unsigned long support_x_rot = support_x_rotated(support_x, support_y, ori);
   unsigned long support_y_rot = support_y_rotated(support_x, support_y, ori);
   /* compute gaussian */
   return lib_image::gaussian_2D(
      sigma_x, sigma_y, ori, deriv, hlbrt, support_x_rot, support_y_rot
   );
}

matrix<> lib_image::gaussian_2D(
   double        sigma_x, 
   double        sigma_y,
   double        ori,
   unsigned int  deriv,
   bool          hlbrt,
   unsigned long support_x,
   unsigned long support_y)
{
   /* compute size of larger grid for axis-aligned gaussian   */
   /* (reverse rotate corners of bounding box by orientation) */
   unsigned long support_x_rot = support_x_rotated(support_x, support_y, -ori);
   unsigned long support_y_rot = support_y_rotated(support_x, support_y, -ori);
   /* compute 1D kernels */
   matrix<> mx = lib_image::gaussian(sigma_x, 0,     false, support_x_rot);
   matrix<> my = lib_image::gaussian(sigma_y, deriv, hlbrt, support_y_rot);
   /* compute 2D kernel from product of 1D kernels */
   matrix<> m(mx._size, my._size);
   unsigned long n = 0;
   for (unsigned long n_x = 0; n_x < mx._size; n_x++) {
      for (unsigned long n_y = 0; n_y < my._size; n_y++) {
         m._data[n] = mx._data[n_x] * my._data[n_y];
         n++;
      }
   }
   /* rotate 2D kernel by orientation */
   m = lib_image::rotate_2D_crop(m, ori, 2*support_x + 1, 2*support_y + 1);
   /* make zero mean (if returning derivative) */
   if (deriv > 0)
      m -= mean(m);
   /* make unit L1 norm */
   m /= sum(abs(m));
   return m;
}

/*
 * Gaussian center-surround kernel (2D).
 *
 * Specify the standard deviation along each axis, and (optionally) the 
 * orientation (in radians).  In addition, optionally specify the support.
 * The support defaults to 3*sigma (automatically adjusted for rotation).
 *
 * The center-surround kernel is the difference of a Gaussian with the 
 * specified standard deviation and one with a standard deviation scaled 
 * by the specified factor.
 *
 * The kernel is normalized to have unit L1 norm and zero mean.
 */
matrix<> lib_image::gaussian_cs_2D(
   double sigma_x, double sigma_y, double ori, double scale_factor)
{
   /* compute support from rotated corners of original support box */
   unsigned long support_x = static_cast<unsigned long>(math::ceil(3*sigma_x));
   unsigned long support_y = static_cast<unsigned long>(math::ceil(3*sigma_y));
   unsigned long support_x_rot = support_x_rotated(support_x, support_y, ori);
   unsigned long support_y_rot = support_y_rotated(support_x, support_y, ori);
   /* compute gaussian */
   return lib_image::gaussian_cs_2D(
      sigma_x, sigma_y, ori, scale_factor, support_x_rot, support_y_rot
   );
}

matrix<> lib_image::gaussian_cs_2D(
   double        sigma_x, 
   double        sigma_y,
   double        ori,
   double        scale_factor,
   unsigned long support_x,
   unsigned long support_y)
{
   /* compute standard deviation for center kernel */
   double sigma_x_c = sigma_x / scale_factor;
   double sigma_y_c = sigma_y / scale_factor;
   /* compute center and surround kernels */
   matrix<> m_center = lib_image::gaussian_2D(
      sigma_x_c, sigma_y_c, ori, 0, false, support_x, support_y
   );
   matrix<> m_surround = lib_image::gaussian_2D(
      sigma_x, sigma_y, ori, 0, false, support_x, support_y
   );
   /* compute center-surround kernel */
   matrix<> m = m_surround - m_center;
   /* make zero mean and unit L1 norm */
   m -= mean(m);
   m /= sum(abs(m));
   return m;
}

/***************************************************************************
 * Filter banks.
 ***************************************************************************/

/*
 * Get the standard set of filter orientations.
 * 
 * The standard orientations are k*pi/n for k in [0,n) where n is the 
 * number of orientation requested.
 */
array<double> lib_image::standard_filter_orientations(unsigned long n_ori) {
   array<double> oris(n_ori);
   double ori = 0;
   double ori_step = (n_ori > 0) ? (M_PIl / static_cast<double>(n_ori)) : 0;
   for (unsigned long n = 0; n < n_ori; n++, ori += ori_step)
      oris[n] = ori;
   return oris;
}

namespace {
/*
 * Runnable object for creating a filter.
 */
class filter_creator : public runnable {
public:
   /*
    * Constructor.
    */
   explicit filter_creator(
      double        sigma_x,     /* sigma x */
      double        sigma_y,     /* sigma y */
      double        ori,         /* orientation */
      unsigned int  deriv,       /* derivative in y-direction (0, 1, or 2) */
      bool          hlbrt,       /* take hilbert transform in y-direction? */
      unsigned long support_x,   /* x support */
      unsigned long support_y,   /* y support */
      matrix<>&     f)           /* output filter matrix */
    : _sigma_x(sigma_x),
      _sigma_y(sigma_y),
      _ori(ori),
      _deriv(deriv),
      _hlbrt(hlbrt),
      _support_x(support_x),
      _support_y(support_y),
      _f(f)
   { }

   /*
    * Destructor.
    */
   virtual ~filter_creator() { /* do nothing */ }

   /*
    * Create the filter.
    */
   virtual void run() {
      _f = lib_image::gaussian_2D(
         _sigma_x, _sigma_y, _ori, _deriv, _hlbrt, _support_x, _support_y
      );
   }
   
protected:
   double        _sigma_x;       /* sigma x */
   double        _sigma_y;       /* sigma y */
   double        _ori;           /* orientation */
   unsigned int  _deriv;         /* derivative in y-direction (0, 1, or 2) */
   bool          _hlbrt;         /* take hilbert transform in y-direction? */
   unsigned long _support_x;     /* x support */
   unsigned long _support_y;     /* y support */
   matrix<>&     _f;             /* output filter matrix */
};
} /* namespace */

/*
 * Oriented Gaussian derivative filters.
 *
 * Create a filter set consisting of rotated versions of the Gaussian 
 * derivative filter with the specified parameters.
 *
 * The filters are created at the standard orientations unless custom 
 * orientations are specified.
 *
 * Specify the standard deviation (sigma) along the longer principle axis.
 * The standard deviation along the other principle axis is sigma/r where
 * r is the elongation ratio.
 * 
 * Each returned filter is an (s+1) x (s+1) matrix where s = 3*sigma.
 *
 * Filters are created in parallel when possible.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::gaussian_filters(
   unsigned long        n_ori,
   double               sigma,
   unsigned int         deriv,
   bool                 hlbrt,
   double               elongation)
{
   array<double> oris = lib_image::standard_filter_orientations(n_ori);
   return lib_image::gaussian_filters(oris, sigma, deriv, hlbrt, elongation);
}

auto_collection< matrix<>, array_list< matrix<> > > lib_image::gaussian_filters(
   const array<double>& oris,
   double               sigma,
   unsigned int         deriv,
   bool                 hlbrt,
   double               elongation)
{
   /* compute support from sigma */
   unsigned long support = static_cast<unsigned long>(math::ceil(3*sigma));
   double sigma_x = sigma;
   double sigma_y = sigma_x/elongation;
   /* allocate collection to hold filters */
   auto_collection< matrix<>, array_list< matrix<> > > filters(
      new array_list< matrix<> >()
   );
   /* allocate collection of filter creators */
   auto_collection< runnable, list<runnable> > filter_creators(
      new list<runnable>()
   );
   /* setup filter creators */
   unsigned long n_ori = oris.size();
   for (unsigned long n = 0; n < n_ori; n++) {
      auto_ptr< matrix<> > f(new matrix<>());
      auto_ptr<filter_creator> f_creator(
         new filter_creator(
            sigma_x, sigma_y, oris[n], deriv, hlbrt, support, support, *f
         )
      );
      filters->add(*f);
      f.release();
      filter_creators->add(*f_creator);
      f_creator.release();
   }
   /* create filters */
   child_thread::run(*filter_creators);
   return filters;
}

/*
 * Even and odd-symmetric filters for computing oriented edge energy.
 *
 * The even-symmetric filter is a Gaussian second-derivative and the odd-
 * symmetric filter is its Hilbert transform.  The filters are elongated 
 * by a ratio of 3:1 along their principle axis.  Each filter set consists
 * of rotated versions of the same filter.
 *
 * Each returned filter is an (s+1) x (s+1) matrix where s = 3*sigma and 
 * sigma is the specified standard deviation.
 *
 * Filters are created in parallel when possible.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::oe_filters_even(
   unsigned long n_ori,
   double        sigma)
{
   return lib_image::gaussian_filters(n_ori, sigma, 2, false, 3.0);
}

auto_collection< matrix<>, array_list< matrix<> > > lib_image::oe_filters_odd(
   unsigned long n_ori,
   double        sigma)
{
   return lib_image::gaussian_filters(n_ori, sigma, 2, true, 3.0);
}

void lib_image::oe_filters(
   unsigned long                                        n_ori,
   double                                               sigma,
   auto_collection< matrix<>, array_list< matrix<> > >& filters_even,
   auto_collection< matrix<>, array_list< matrix<> > >& filters_odd)
{
   /* runnable class for creating oe filters */
   class oe_filters_creator : public runnable {
   public:
      /*
       * Constructor.
       */
      explicit oe_filters_creator(
         unsigned long                                        n_ori,
         double                                               sigma,
         bool                                                 even_or_odd,
         auto_collection< matrix<>, array_list< matrix<> > >& filters)
       : _n_ori(n_ori),
         _sigma(sigma),
         _even_or_odd(even_or_odd),
         _filters(filters)
      { }

      /*
       * Destructor.
       */
      virtual ~oe_filters_creator() { /* do nothing */ }

      /*
       * Create the filter set.
       */
      virtual void run() {
         _filters = _even_or_odd ?
            lib_image::oe_filters_odd(_n_ori, _sigma)
          : lib_image::oe_filters_even(_n_ori, _sigma);
      }

   protected:
      unsigned long                                        _n_ori;
      double                                               _sigma;
      bool                                                 _even_or_odd;
      auto_collection< matrix<>, array_list< matrix<> > >& _filters;
   };
   /* create oe filters */
   oe_filters_creator f_even_creator(n_ori, sigma, false, filters_even);
   oe_filters_creator f_odd_creator(n_ori, sigma, true, filters_odd);
   child_thread::run(f_even_creator, f_odd_creator);
}

/*
 * Filters for computing textons.
 *
 * The set of texton filters is the union of the even and odd-symmetric
 * filter sets above in addition to a center-surround difference of 
 * Gaussians filter.
 *
 * Each returned filter is an (s+1) x (s+1) matrix where s = 3*sigma and 
 * sigma is the specified standard deviation.
 *
 * Filters are created in parallel when possible.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::texton_filters(
   unsigned long n_ori,
   double        sigma)
{
   /* allocate collection to hold filters */
   auto_collection< matrix<>, array_list< matrix<> > > filters(
      new array_list< matrix<> >()
   );
   /* get even and odd-symmetric filter sets */
   auto_collection< matrix<>, array_list< matrix<> > > filters_even;
   auto_collection< matrix<>, array_list< matrix<> > > filters_odd;
   lib_image::oe_filters(n_ori, sigma, filters_even, filters_odd);
   /* add even and odd-symmetric filters to collection */
   filters->add(*filters_even); filters_even.release();
   filters->add(*filters_odd);  filters_odd.release();
   /* compute center surround filter */
   unsigned long support = static_cast<unsigned long>(math::ceil(3*sigma));
   matrix<> f_cs = lib_image::gaussian_cs_2D(
      sigma, sigma, 0, M_SQRT2l, support, support
   );
   /* add center surround filter to collection */
   auto_ptr< matrix<> > f_ptr(new matrix<>());
   matrix<>::swap(f_cs, *f_ptr);
   filters->add(*f_ptr);
   f_ptr.release();
   return filters;
}

/***************************************************************************
 * Image filtering.
 ***************************************************************************/

namespace {
/*
 * Abstract base class for runnable filterers.
 */
class filterer_base : public runnable {
public:
   /*
    * Constructor.
    */
   explicit filterer_base(
      const matrix<>& m,      /* image */
      const matrix<>& f,      /* filter */
      matrix<>& m_result)     /* filter output */
    : _m(m), _f(f), _m_result(m_result) { }

   /*
    * Destructor.
    */
   virtual ~filterer_base() { /* do nothing */ }

   /*
    * Run the filter.
    */
   virtual void run() = 0;
   
protected:
   const matrix<>& _m;        /* image */
   const matrix<>& _f;        /* filter */
   matrix<>&       _m_result; /* filter output */
};

/*
 * Runnable object for convolving an image with a filter.
 */
class filterer : public filterer_base {
public:
   /*
    * Constructor.
    */
   explicit filterer(
      const matrix<>& m,      /* image */
      const matrix<>& f,      /* filter */
      matrix<>& m_result)     /* filter output */
    : filterer_base(m, f, m_result) { }

   /*
    * Destructor.
    */
   virtual ~filterer() { /* do nothing */ }

   /*
    * Convolve image with filter.
    */
   virtual void run() {
      _m_result = conv_crop(_m, _f);
   }
};

/*
 * Runnable object for convolving an image with a filter and squaring
 * the result.
 */
class filterer_sq : public filterer_base {
public:
   /*
    * Constructor.
    */
   explicit filterer_sq(
      const matrix<>& m,      /* image */
      const matrix<>& f,      /* filter */
      matrix<>& m_result)     /* filter output */
    : filterer_base(m, f, m_result) { }

   /*
    * Destructor.
    */
   virtual ~filterer_sq() { /* do nothing */ }

   /*
    * Convolve image with filter and square result.
    */
   virtual void run() {
      matrix<> m_conv = conv_crop(_m, _f);
      _m_result = prod(m_conv, m_conv);
   }
};

/*
 * Abstract base class for runnable rectifying filterers.
 */
class filterer_rectified_base : public runnable {
public:
   /*
    * Constructor.
    */
   explicit filterer_rectified_base(
      const matrix<>& m,      /* image */
      const matrix<>& f,      /* filter */
      matrix<>& m_result_neg, /* negative rectified response */
      matrix<>& m_result_pos) /* positive rectified response */
    : _m(m), _f(f), _m_result_neg(m_result_neg), _m_result_pos(m_result_pos)
   { }
   
   /*
    * Destructor.
    */
   virtual ~filterer_rectified_base() { /* do nothing */ }

   /*
    * Run the filter.
    */
   virtual void run() = 0;
   
protected:
   const matrix<>& _m;              /* image */
   const matrix<>& _f;              /* filter */
   matrix<>&       _m_result_neg;   /* negative filter output */
   matrix<>&       _m_result_pos;   /* positive filter output */
};

/*
 * Runnable object for rectified filtering.
 */
class filterer_rectified : public filterer_rectified_base {
public:
   /*
    * Constructor.
    */
   explicit filterer_rectified(
      const matrix<>& m,      /* image */
      const matrix<>& f,      /* filter */
      matrix<>& m_result_neg, /* negative rectified response */
      matrix<>& m_result_pos) /* positive rectified response */
    : filterer_rectified_base(m, f, m_result_neg, m_result_pos) { }

   /*
    * Destructor.
    */
   virtual ~filterer_rectified() { /* do nothing */ }

   /*
    * Convolve image with filter and rectify.
    */
   virtual void run() {
      /* convolve with filter */
      matrix<> m_conv = conv_crop(_m, _f);
      /* resize result matrices */
      _m_result_neg.resize(0);
      _m_result_pos.resize(0);
      array<unsigned long> dims = m_conv.dimensions();
      _m_result_neg.resize(dims);
      _m_result_pos.resize(dims);
      /* compute rectified signals */
      unsigned long n_inds = m_conv.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         double v = m_conv[n];
         if (v < 0)
            _m_result_neg[n] = -v;
         else
            _m_result_pos[n] = v;
      }
   }
};

/*
 * Runnable object for squared rectified filtering.
 */
class filterer_rectified_sq : public filterer_rectified_base {
public:
   /*
    * Constructor.
    */
   explicit filterer_rectified_sq(
      const matrix<>& m,      /* image */
      const matrix<>& f,      /* filter */
      matrix<>& m_result_neg, /* negative rectified response */
      matrix<>& m_result_pos) /* positive rectified response */
    : filterer_rectified_base(m, f, m_result_neg, m_result_pos) { }

   /*
    * Destructor.
    */
   virtual ~filterer_rectified_sq() { /* do nothing */ }

   /*
    * Convolve image with filter, rectify, and square.
    */
   virtual void run() {
      /* convolve with filter */
      matrix<> m_conv = conv_crop(_m, _f);
      /* resize result matrices */
      _m_result_neg.resize(0);
      _m_result_pos.resize(0);
      array<unsigned long> dims = m_conv.dimensions();
      _m_result_neg.resize(dims);
      _m_result_pos.resize(dims);
      /* compute rectified signals */
      unsigned long n_inds = m_conv.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         double v  = m_conv[n];
         double v2 = v*v;
         if (v < 0)
            _m_result_neg[n] = v2;
         else
            _m_result_pos[n] = v2;
      }
   }
};
} /* namespace */

/*
 * Image filtering.
 *
 * Return the convolution of the image with each filter (cropped so that 
 * the result is the same size as the original image).
 *
 * Filter responses are computed in parallel when possible.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::filter(
   const matrix<>&               m,
   const collection< matrix<> >& filters)
{
   /* allocate collection to hold filter responses */
   auto_collection< matrix<>, array_list< matrix<> > > responses(
      new array_list< matrix<> >()
   );
   /* allocate collection of filterers */
   auto_collection< runnable, list<runnable> > filterers(new list<runnable>());
   /* setup filterers */
   for (auto_ptr< iterator< matrix<> > > i = filters.iter_create();
        i->has_next(); )
   {
      matrix<>& f = i->next();
      auto_ptr< matrix<> > m_result(new matrix<>());
      auto_ptr<filterer> filt(new filterer(m, f, *m_result));
      responses->add(*m_result);
      m_result.release();
      filterers->add(*filt);
      filt.release();
   }
   /* run filterers */
   child_thread::run(*filterers);
   return responses;
}

/*
 * Image filtering.
 *
 * Return the square of the convolution of the image with each filter
 * (cropped so that the result is the same size as the original image).
 *
 * Filter responses are computed in parallel when possible.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::filter_sq(
   const matrix<>&               m,
   const collection< matrix<> >& filters)
{
   /* allocate collection to hold filter responses */
   auto_collection< matrix<>, array_list< matrix<> > > responses(
      new array_list< matrix<> >()
   );
   /* allocate collection of filterers */
   auto_collection< runnable, list<runnable> > filterers(new list<runnable>());
   /* setup filterers */
   for (auto_ptr< iterator< matrix<> > > i = filters.iter_create();
        i->has_next(); )
   {
      matrix<>& f = i->next();
      auto_ptr< matrix<> > m_result(new matrix<>());
      auto_ptr<filterer_sq> filt(new filterer_sq(m, f, *m_result));
      responses->add(*m_result);
      m_result.release();
      filterers->add(*filt);
      filt.release();
   }
   /* run filterers */
   child_thread::run(*filterers);
   return responses;
}

/*
 * Image filtering and rectification.
 *
 * Return the rectified components of the convolution of the image with
 * each filter (cropped so that the result is the same size as the original
 * image).
 *
 * Filter responses are computed in parallel when possible.
 */
void lib_image::filter_rectified(
   const matrix<>&                                      m,
   const collection< matrix<> >&                        filters,
   auto_collection< matrix<>, array_list< matrix<> > >& responses_neg,
   auto_collection< matrix<>, array_list< matrix<> > >& responses_pos)
{
   /* allocate collections to hold filter responses */
   responses_neg.reset(new array_list< matrix<> >());
   responses_pos.reset(new array_list< matrix<> >());
   /* allocate collection of filterers */
   auto_collection< runnable, list<runnable> > filterers(new list<runnable>());
   /* setup filterers */
   for (auto_ptr< iterator< matrix<> > > i = filters.iter_create();
        i->has_next(); )
   {
      matrix<>& f = i->next();
      auto_ptr< matrix<> > m_result_neg(new matrix<>());
      auto_ptr< matrix<> > m_result_pos(new matrix<>());
      auto_ptr<filterer_rectified> filt(
         new filterer_rectified(m, f, *m_result_neg, *m_result_pos)
      );
      responses_neg->add(*m_result_neg); m_result_neg.release();
      responses_pos->add(*m_result_pos); m_result_pos.release();
      filterers->add(*filt);             filt.release();
   }
   /* run filterers */
   child_thread::run(*filterers);
}
 
/*
 * Image filtering and rectification.
 *
 * Return the square of the rectified components of the convolution of the
 * image with each filter (cropped so that the result is the same size as
 * the original image).
 *
 * Filter responses are computed in parallel when possible.
 */
void lib_image::filter_rectified_sq(
   const matrix<>&                                      m,
   const collection< matrix<> >&                        filters,
   auto_collection< matrix<>, array_list< matrix<> > >& responses_neg,
   auto_collection< matrix<>, array_list< matrix<> > >& responses_pos)
{
   /* allocate collections to hold filter responses */
   responses_neg.reset(new array_list< matrix<> >());
   responses_pos.reset(new array_list< matrix<> >());
   /* allocate collection of filterers */
   auto_collection< runnable, list<runnable> > filterers(new list<runnable>());
   /* setup filterers */
   for (auto_ptr< iterator< matrix<> > > i = filters.iter_create();
        i->has_next(); )
   {
      matrix<>& f = i->next();
      auto_ptr< matrix<> > m_result_neg(new matrix<>());
      auto_ptr< matrix<> > m_result_pos(new matrix<>());
      auto_ptr<filterer_rectified_sq> filt(
         new filterer_rectified_sq(m, f, *m_result_neg, *m_result_pos)
      );
      responses_neg->add(*m_result_neg); m_result_neg.release();
      responses_pos->add(*m_result_pos); m_result_pos.release();
      filterers->add(*filt);             filt.release();
   }
   /* run filterers */
   child_thread::run(*filterers);
}

namespace {
/*
 * Runnable object for stacking multiple matrices into a matrix of vectors.
 * Automatically split the stacking task into parallel subtasks.
 */
class vectorizer : public runnable {
public:
   /*
    * Constructor.
    * The result matrix must have been initialized to an appropriately sized
    * matrix of length zero vectors.
    */
   explicit vectorizer(
      const array_list< matrix<> >& m_arr,   /* matrices to stack */
      unsigned long                 start,   /* first index to process */
      unsigned long                 end,     /* last index to process */
      matrix< matrix<> >&           m)       /* resulting matrix of vectors */
    : _m_arr(m_arr), _start(start), _end(end), _m(m) { }
   
   /*
    * Destructor.
    */
   virtual ~vectorizer() { /* do nothing */ }

   /*
    * Stack matrices.
    */
   virtual void run() {
      if ((thread::processors() > 1) && (_start < _end)) {
         /* split stacking task */
         unsigned long mid = (_start + _end)/2;
         vectorizer v_left(_m_arr, _start, mid, _m);
         vectorizer v_right(_m_arr, mid+1, _end, _m);
         child_thread::run(v_left, v_right);
      } else {
         /* stack matrices */
         unsigned long length = _m_arr.size();
         for (unsigned long n = _start; n <= _end; n++) {
            matrix<>& v = _m[n];
            v.resize(length);
            for (unsigned long ind = 0; ind < length; ind++)
               v[ind] = (_m_arr[ind])[n];
         }
      }
   }

protected:
   const array_list< matrix<> >& _m_arr;     /* matrices to stack */
   unsigned long                 _start;     /* first index to process */
   unsigned long                 _end;       /* last index to process */
   matrix< matrix<> >&           _m;         /* resulting matrix of vectors */
};
} /* namespace */

/*
 * Stack the output of multiple filters into column vectors.
 */
matrix< matrix<> > lib_image::vectorize_filter_responses(
   const collection< matrix<> >& filter_responses)
{
   /* create array_list of responses */
   array_list< matrix<> > responses(filter_responses);
   /* check number of filters */
   unsigned long n_filt = responses.size();
   if (n_filt == 0)
      return matrix< matrix<> >();
   /* check that all response matrices are the same size */
   for (unsigned long n = 1; n < n_filt; n++)
      matrix<>::assert_dims_equal(responses[0]._dims, responses[n]._dims);
   /* stack matrices */
   matrix< matrix<> > m(responses[0]._dims);
   if (m._size > 0) {
      vectorizer v(responses, 0, (m._size-1), m);
      v.run();
   }
   return m;
}

/***************************************************************************
 * Image quantization.
 ***************************************************************************/

/*
 * Vector quantize filter responses.
 * Return both the cluster assignments and cluster centroids.
 */
matrix<unsigned long> lib_image::cluster_filter_responses(
   const collection< matrix<> >&                        responses,
   const centroid_clusterer< matrix<> >&                clusterer,
   auto_collection< matrix<>, array_list< matrix<> > >& centroids)
{
   /* stack filter responses into vectors */
   matrix< matrix<> > vecs = lib_image::vectorize_filter_responses(responses);
   /* cluster vectors */
   return lib_image::cluster_filter_responses(vecs, clusterer, centroids);
}

matrix<unsigned long> lib_image::cluster_filter_responses(
   const matrix< matrix<> >&                            responses,
   const centroid_clusterer< matrix<> >&                clusterer,
   auto_collection< matrix<>, array_list< matrix<> > >& centroids)
{
   /* create collection of vectors */
   list< matrix<> > vec_list;
   for (unsigned long n = 0; n < responses._size; n++)
      vec_list.add(responses._data[n]);
   /* cluster */
   array<unsigned long> assign = clusterer.cluster(vec_list, centroids);
   /* convert assignment array to assignment matrix */
   matrix<unsigned long> m_assign = matrix<unsigned long>::to_matrix(assign);
   m_assign.reshape(responses._dims);
   return m_assign;
}

/*
 * Cluster image values.
 * Return both the cluster assignments and cluster centroids.
 */
matrix<unsigned long> lib_image::cluster_values(
   const matrix<>&                                m,
   const centroid_clusterer<double>&              clusterer,
   auto_collection< double, array_list<double> >& centroids)
{
   /* create collection of values */
   list<double> value_list;
   for (unsigned long n = 0; n < m._size; n++)
      value_list.add(m._data[n]);
   /* cluster */
   array<unsigned long> assign = clusterer.cluster(value_list, centroids);
   /* convert assignment array to assignment matrix */
   matrix<unsigned long> m_assign = matrix<unsigned long>::to_matrix(assign);
   m_assign.reshape(m._dims);
   return m_assign;
}

/*
 * Quantize image values into uniformly spaced bins in [0,1].
 * Return the assignments and bin centroids.
 */
matrix<unsigned long> lib_image::quantize_values(
   const matrix<>&                                m,
   unsigned long                                  n_bins)
{
   auto_collection< double, array_list<double> > centroids;
   return lib_image::quantize_values(m, n_bins, centroids);
}

matrix<unsigned long> lib_image::quantize_values(
   const matrix<>&                                m,
   unsigned long                                  n_bins,
   auto_collection< double, array_list<double> >& centroids)
{
   /* check arguments */
   if (n_bins == 0)
      throw ex_invalid_argument("n_bins must be > 0");
   /* create centroids */
   centroids.reset(new array_list<double>());
   for (unsigned long n = 0; n < n_bins; n++) {
      auto_ptr<double> val(
         new double(
            (static_cast<double>(n) + 0.5)/static_cast<double>(n_bins)
         )
      );
      centroids->add(*val);
      val.release();
   }
   /* compute assignments */
   matrix<unsigned long> assign(m._dims);
   for (unsigned long n = 0; n < m._size; n++) {
      unsigned long bin = static_cast<unsigned long>(
         math::floor(m._data[n]*static_cast<double>(n_bins))
      );
      if (bin == n_bins) { bin = n_bins-1; }
      assign._data[n] = bin;
   }
   return assign;
}

/*
 * Construct a quantized image by looking up the value of the centroid to
 * which each element is assigned.
 */
matrix<> lib_image::quantized_image(
   const matrix<unsigned long>& assignments,
   const array_list<double>&    centroids)
{
   matrix<> m_quantized(assignments._dims);
   for (unsigned long n = 0; n < assignments._size; n++)
      m_quantized._data[n] = centroids[assignments._data[n]];
   return m_quantized;
}

/***************************************************************************
 * Edge detection.
 ***************************************************************************/

/*
 * Compute oriented edge energy.
 *
 * Given that fe and fo are the even and odd-symmetric filters at a specific
 * orientation, the even and odd components of oriented energy (oe_e, oe_o)
 * at that orientation for an image I are defined as:
 *
 *                oe_e = (I ** fe)
 *                oe_o = (I ** fo)
 * 
 * where ** denotes convolution.
 */
void lib_image::edge_oe(
   const matrix<>& m,
   const matrix<>& fe,
   const matrix<>& fo,
   auto_ptr< matrix<> >& oe_even,
   auto_ptr< matrix<> >& oe_odd)
{
   /* allocate result */
   oe_even.reset(new matrix<>());
   oe_odd.reset(new matrix<>());
   /* setup filterers */
   filterer filt_even(m, fe, *oe_even);
   filterer filt_odd(m, fo, *oe_odd);
   /* compute oe */
   child_thread::run(filt_even, filt_odd);
}

/*
 * Compute oriented edge energy (at multiple orientations/scales).
 *
 * Different orientations/scales are processed in parallel when possible.
 */
void lib_image::edge_oe(
   const matrix<>&                                      m,
   const collection< matrix<> >&                        fe_set,
   const collection< matrix<> >&                        fo_set,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_even,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_odd)
{
   /* check arguments */
   unsigned long n_filt = fe_set.size();
   if (n_filt != fo_set.size())
      throw ex_invalid_argument(
         "even and odd-symmetric filter sets do not match"
      );
   /* place all filters in one collection */
   array_list< matrix<> > filters;
   filters.add(fe_set);
   filters.add(fo_set);
   /* compute filter responses */
   auto_collection< matrix<>, array_list< matrix<> > > responses =
      lib_image::filter(m, filters);
   /* allocate collections to hold rectified oe */
   oes_even.reset(new array_list< matrix<> >());
   oes_odd.reset(new array_list< matrix<> >());
   /* read off even components */
   for (unsigned long n = 0; n < n_filt; n++) {
      auto_ptr< matrix<> > oe_even(&(responses->remove_head()));
      oes_even->add(*oe_even);
      oe_even.release();
   }
   /* read off odd components */
   for (unsigned long n = 0; n < n_filt; n++) {
      auto_ptr< matrix<> > oe_odd(&(responses->remove_head()));
      oes_odd->add(*oe_odd);
      oe_odd.release();
   }
}

/*
 * Compute rectified oriented edge energy.
 *
 * Given that fo is the odd-symmetric filter at a specified orientation, the
 * rectified oriented energy (oe-, oe+) at that orientation for an image I 
 * is defined as:
 *
 *                oe- = [I ** fo]-
 *                oe+ = [I ** fo]+
 *
 * where ** denotes convolution, and -/+ denote negative/positive
 * rectification.
 */
void lib_image::edge_oe_rectified(
   const matrix<>&       m,
   const matrix<>&       fo,
   auto_ptr< matrix<> >& oe_neg,
   auto_ptr< matrix<> >& oe_pos)
{
   /* allocate result */
   oe_neg.reset(new matrix<>());
   oe_pos.reset(new matrix<>());
   /* compute rectified oe */
   filterer_rectified filt_rect(m, fo, *oe_neg, *oe_pos);
   filt_rect.run();
}

/*
 * Compute rectified oriented edge energy (at multiple orientations/scales).
 *
 * Different orientations/scales are processed in parallel when possible.
 */
void lib_image::edge_oe_rectified(
   const matrix<>&                                      m,
   const collection< matrix<> >&                        fo_set,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_neg,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_pos)
{
   lib_image::filter_rectified(m, fo_set, oes_neg, oes_pos);
}

namespace {
/*
 * Runnable class for combining oe components.
 */
class oe_combiner : public runnable {
public:
   /*
    * Constructor.
    */
   explicit oe_combiner(
      const matrix<>& oe_even,      /* oe_e */
      const matrix<>& oe_odd,       /* oe_o */
      matrix<>&       oe_strength,  /* returned oe */
      matrix<>&       oe_phase)     /* returned oe phase */
    : _oe_even(oe_even),
      _oe_odd(oe_odd),
      _oe_strength(oe_strength),
      _oe_phase(oe_phase)
   { }

   /*
    * Destructor.
    */
   virtual ~oe_combiner() { /* do nothing */ }

   /*
    * Compute oe strength and phase.
    */
   virtual void run() {
      _oe_strength =
         sqrt(prod(_oe_even, _oe_even) + prod(_oe_odd, _oe_odd)).real();
      _oe_phase = atan2(_oe_odd, _oe_even);
   }
   
protected:
   const matrix<>& _oe_even;        /* oe_e */
   const matrix<>& _oe_odd;         /* oe_o */
   matrix<>&       _oe_strength;    /* returned oe */
   matrix<>&       _oe_phase;       /* returned oe phase */
};
} /* namespace */

/*
 * Combine components of oriented edge energy into edge maps.
 *
 * For each pair of components, the oriented energy is computed as:
 *
 *             oe_strength = sqrt((oe_e)^2 + (oe_o)^2)
 *             oe_phase    = atan2(oe_o, oe_e)
 *
 * where oe_o and oe_e are the even and odd components of oe, respectively.
 */
void lib_image::combine_edge_oe(
   const collection< matrix<> >&                        oes_even,
   const collection< matrix<> >&                        oes_odd,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_strength,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_phase)
{
   /* check arguments - number of components */
   if (oes_even.size() != oes_odd.size()) {
      throw ex_invalid_argument(
         "number of oe_even matrices different from number of oe_odd matrices"
      );
   }
   /* initialize result */
   oes_strength.reset(new array_list< matrix<> >());
   oes_phase.reset(new array_list< matrix<> >());
   /* allocate collection to hold combiners */
   auto_collection< runnable, list<runnable> > oe_combiners(
      new list<runnable>()
   );
   /* setup oe combiners */
   auto_ptr< iterator< matrix<> > > i_even = oes_even.iter_create();
   auto_ptr< iterator< matrix<> > > i_odd  = oes_odd.iter_create();
   while (i_even->has_next()) {
      auto_ptr< matrix<> > oe_strength(new matrix<>());
      auto_ptr< matrix<> > oe_phase(new matrix<>());
      auto_ptr<oe_combiner> oe_comb(
         new oe_combiner(
            i_even->next(), i_odd->next(), *oe_strength, *oe_phase
         )
      );
      oes_strength->add(*oe_strength); oe_strength.release();
      oes_phase->add(*oe_phase);       oe_phase.release();
      oe_combiners->add(*oe_comb);     oe_comb.release();
   }
   /* run oe combiners */
   child_thread::run(*oe_combiners);
}

namespace {
/*
 * Runnable class for combining rectified oe components.
 */
class oe_rectified_combiner : public runnable {
public:
   /*
    * Constructor.
    */
   explicit oe_rectified_combiner(
      const matrix<>& oe_neg,       /* oe- */
      const matrix<>& oe_pos,       /* oe+ */
      matrix<>&       oe_strength,  /* returned oe */
      matrix<>&       oe_polarity)  /* returned oe polarity */
    : _oe_neg(oe_neg),
      _oe_pos(oe_pos),
      _oe_strength(oe_strength),
      _oe_polarity(oe_polarity)
   { }

   /*
    * Destructor.
    */
   virtual ~oe_rectified_combiner() { /* do nothing */ }

   /*
    * Compute oe strength and phase.
    */
   virtual void run() {
      _oe_strength = _oe_neg + _oe_pos;
      _oe_polarity.resize(0);
      _oe_polarity.resize(_oe_strength.dimensions());
      unsigned long size = _oe_polarity.size();
      for (unsigned long n = 0; n < size; n++) {
         if (_oe_neg[n] > 0)
            _oe_polarity[n] = -1;
         else if (_oe_pos[n] > 0)
            _oe_polarity[n] = 1;
      }
   }
   
protected:
   const matrix<>& _oe_neg;         /* oe- */
   const matrix<>& _oe_pos;         /* oe+ */
   matrix<>&       _oe_strength;    /* returned oe */
   matrix<>&       _oe_polarity;    /* returned oe polarity */
};
} /* namespace */

/*
 * Combine components of rectified oriented edge energy into edge maps.
 *
 * For each pair of components, the oriented energy is computed as:
 *
 *             oe_strength = sqrt((oe-)^2 + (oe+)^2)
 *                         = (oe-) + (oe+)
 *
 *             oe_polarity = sign(-(oe-) + (oe+))
 *
 * where oe- and oe+ are the negative and positive rectified components
 * of oe, respectively.
 *
 * The polarity is -1 at locations for which oe- was the dominant edge, +1
 * at locations for which oe+ was the dominant edge, and 0 at locations 
 * with no edge strength.
 */
void lib_image::combine_edge_oe_rectified(
   const collection< matrix<> >&                        oes_neg,
   const collection< matrix<> >&                        oes_pos,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_strength,
   auto_collection< matrix<>, array_list< matrix<> > >& oes_polarity)
{
   /* check arguments - number of components */
   if (oes_neg.size() != oes_pos.size()) {
      throw ex_invalid_argument(
         "number of oe- matrices different from number of oe+ matrices"
      );
   }
   /* initialize result */
   oes_strength.reset(new array_list< matrix<> >());
   oes_polarity.reset(new array_list< matrix<> >());
   /* allocate collection to hold combiners */
   auto_collection< runnable, list<runnable> > oe_combiners(
      new list<runnable>()
   );
   /* setup oe combiners */
   auto_ptr< iterator< matrix<> > > i_neg = oes_neg.iter_create();
   auto_ptr< iterator< matrix<> > > i_pos = oes_pos.iter_create();
   while (i_neg->has_next()) {
      auto_ptr< matrix<> > oe_strength(new matrix<>());
      auto_ptr< matrix<> > oe_polarity(new matrix<>());
      auto_ptr<oe_rectified_combiner> oe_comb(
         new oe_rectified_combiner(
            i_neg->next(), i_pos->next(), *oe_strength, *oe_polarity
         )
      );
      oes_strength->add(*oe_strength); oe_strength.release();
      oes_polarity->add(*oe_polarity); oe_polarity.release();
      oe_combiners->add(*oe_comb);     oe_comb.release();
   }
   /* run oe combiners */
   child_thread::run(*oe_combiners);
}

namespace {
/*
 * Runnable object for computing the index of the maximum value in a stack of
 * matrices.  Computation is performed in parallel across the result matrix.
 */
class matrix_max_selector : public runnable {
public:
   /*
    * Constructor.
    */
   explicit matrix_max_selector(
      const array_list< matrix<> >& matrices,   /* stack of matrices */
      unsigned long                 start,      /* first index to process */
      unsigned long                 end,        /* last index to process */
      matrix<unsigned long>&        index)      /* returned index matrix */
    : _matrices(matrices), _start(start), _end(end), _index(index)
   { }

   /*
    * Destructor.
    */
   virtual ~matrix_max_selector() { /* do nothing */ }

   /*
    * Find indices of maximum values.
    */
   virtual void run() {
      if ((thread::processors() > 1) && (_start < _end)) {
         /* split task */
         unsigned long mid = (_start + _end)/2;
         matrix_max_selector m_left(_matrices, _start, mid, _index);
         matrix_max_selector m_right(_matrices, mid+1, _end, _index);
         child_thread::run(m_left, m_right);
      } else {
         /* compute indices of maximum values */
         unsigned long n_matrices = _matrices.size();
         if (n_matrices > 0) {
            for (unsigned long n = _start; n <= _end; n++) {
               double        max_val = (_matrices[0])[n];
               unsigned long max_ind = 0;
               for (unsigned long ind = 1; ind < n_matrices; ind++) {
                  double val = (_matrices[ind])[n];
                  if (val > max_val) {
                     max_val = val;
                     max_ind = ind;
                  }
               }
               _index[n] = max_ind;
            }
         }
      }
   }

   /*
    * Select values from a matrix stack given the indices.
    */
   static matrix<> select_values(
      const array_list< matrix<> >& matrices,   /* stack of matrices */
      const matrix<unsigned long>&  index)      /* index */
   {
      matrix<> m(index.dimensions());
      unsigned long size = m.size();
      for (unsigned long n = 0; n < size; n++)
         m[n] = (matrices[index[n]])[n];
      return m;
   }

   /*
    * Select values from an array given the indices.
    */
   static matrix<> select_values(
      const array<double>&          arr,        /* array */
      const matrix<unsigned long>&  index)      /* index */
   {
      matrix<> m(index.dimensions());
      unsigned long size = m.size();
      for (unsigned long n = 0; n < size; n++)
         m[n] = arr[index[n]];
      return m;
   }

protected:
   const array_list< matrix<> >& _matrices;     /* stack of matrices */
   unsigned long                 _start;        /* first index to process */
   unsigned long                 _end;          /* last index to process */
   matrix<unsigned long>&        _index;        /* returned index matrix */
};
} /* namespace */

/*
 * Combine oriented edge energy into a single edge map.
 *
 * Each location is assigned edge strength equal to the maximum oe value
 * over all orientations.  The phase, atan2(oe_o, oe_e), for the maximum
 * oe value at each location is also returned.  In addition, the dominant
 * orientation at each location is returned.
 *
 * Unless otherwise specified, the orientations corresponding to the input
 * collection of n oriented edge maps are assumed to be k*pi/n for integer
 * k in [0,n).
 */
void lib_image::combine_edge_oe(
   const collection< matrix<> >& oes_even,
   const collection< matrix<> >& oes_odd,
   auto_ptr< matrix<> >&         edge,
   auto_ptr< matrix<> >&         edge_phase,
   auto_ptr< matrix<> >&         edge_ori)
{
   unsigned long n_ori = oes_even.size();
   array<double> oris = lib_image::standard_filter_orientations(n_ori);
   lib_image::combine_edge_oe(
      oes_even, oes_odd, oris, edge, edge_phase, edge_ori
   );
}

void lib_image::combine_edge_oe(
   const collection< matrix<> >& oes_even,
   const collection< matrix<> >& oes_odd,
   const array<double>&          oris,
   auto_ptr< matrix<> >&         edge,
   auto_ptr< matrix<> >&         edge_phase,
   auto_ptr< matrix<> >&         edge_ori)
{
   /* check arguments - number of orientations */
   unsigned long n_ori = oris.size();
   if (oes_even.size() != n_ori)
      throw ex_invalid_argument(
         "number of oe_even matrices different from number of orientations"
      );
   if (oes_odd.size() != n_ori)
      throw ex_invalid_argument(
         "number of oe_odd matrices different from number of orientations"
      );
   /* initialize edge map */
   edge.reset(new matrix<>());
   edge_phase.reset(new matrix<>());
   edge_ori.reset(new matrix<>());
   /* check if nontrivial result */
   if (n_ori > 0) {
      /* combine oe components */
      auto_collection< matrix<>, array_list< matrix<> > > oes_strength;
      auto_collection< matrix<>, array_list< matrix<> > > oes_phase;
      lib_image::combine_edge_oe(
         oes_even, oes_odd, oes_strength, oes_phase
      );
      /* check that all combined components have the same size */
      array<unsigned long> dims = (*oes_strength)[0]._dims;
      for (unsigned long n = 0; n < n_ori; n++) {
         matrix<>::assert_dims_equal(dims, (*oes_strength)[n]._dims);
         matrix<>::assert_dims_equal(dims, (*oes_phase)[n]._dims);
      }
      /* compute indices of maximum oe values */
      matrix<unsigned long> index(dims);
      if (index._size > 0) {
         matrix_max_selector m_select(*oes_strength, 0, index._size-1, index);
         m_select.run();
      }
      /* select responses from channel with max oe value */ 
      *edge       = matrix_max_selector::select_values(*oes_strength, index);
      *edge_phase = matrix_max_selector::select_values(*oes_phase, index);
      *edge_ori   = matrix_max_selector::select_values(oris, index);
   }
}

/*
 * Combine rectified oriented edge energy into a single edge map.
 *
 * Each location is assigned edge strength equal to the maximum oe value
 * over all orientations.  The polarity, sign(-(oe-) + (oe+)), for the 
 * maximum oe value at each location is also returned.  In addition, the
 * dominant orientation at each location is returned. 
 *
 * Unless otherwise specified, the orientations corresponding to the input
 * collection of n oriented edge maps are assumed to be k*pi/n for integer
 * k in [0,n).
 */
void lib_image::combine_edge_oe_rectified(
   const collection< matrix<> >& oes_neg,
   const collection< matrix<> >& oes_pos,
   auto_ptr< matrix<> >&         edge,
   auto_ptr< matrix<> >&         edge_polarity,
   auto_ptr< matrix<> >&         edge_ori)
{
   unsigned long n_ori = oes_neg.size();
   array<double> oris = lib_image::standard_filter_orientations(n_ori);
   lib_image::combine_edge_oe_rectified(
      oes_neg, oes_pos, oris, edge, edge_polarity, edge_ori
   );
}

void lib_image::combine_edge_oe_rectified(
   const collection< matrix<> >& oes_neg,
   const collection< matrix<> >& oes_pos,
   const array<double>&          oris,
   auto_ptr< matrix<> >&         edge,
   auto_ptr< matrix<> >&         edge_polarity,
   auto_ptr< matrix<> >&         edge_ori)
{
   /* check arguments - number of orientations */
   unsigned long n_ori = oris.size();
   if (oes_neg.size() != n_ori)
      throw ex_invalid_argument(
         "number of oe- matrices different from number of orientations"
      );
   if (oes_pos.size() != n_ori)
      throw ex_invalid_argument(
         "number of oe+ matrices different from number of orientations"
      );
   /* initialize edge map */
   edge.reset(new matrix<>());
   edge_polarity.reset(new matrix<>());
   edge_ori.reset(new matrix<>());
   /* check if nontrivial result */
   if (n_ori > 0) {
      /* combine oe components */
      auto_collection< matrix<>, array_list< matrix<> > > oes_strength;
      auto_collection< matrix<>, array_list< matrix<> > > oes_polarity;
      lib_image::combine_edge_oe_rectified(
         oes_neg, oes_pos, oes_strength, oes_polarity
      );
      /* check that all combined components have the same size */
      array<unsigned long> dims = (*oes_strength)[0]._dims;
      for (unsigned long n = 0; n < n_ori; n++) {
         matrix<>::assert_dims_equal(dims, (*oes_strength)[n]._dims);
         matrix<>::assert_dims_equal(dims, (*oes_polarity)[n]._dims);
      }
      /* compute indices of maximum oe values */
      matrix<unsigned long> index(dims);
      if (index._size > 0) {
         matrix_max_selector m_select(*oes_strength, 0, index._size-1, index);
         m_select.run();
      }
      /* select responses from channel with max oe value */ 
      *edge          = matrix_max_selector::select_values(*oes_strength, index);
      *edge_polarity = matrix_max_selector::select_values(*oes_polarity, index);
      *edge_ori      = matrix_max_selector::select_values(oris, index);
   }
}

/***************************************************************************
 * Textons.
 ***************************************************************************/

/*
 * Compute textons using the given filter set.
 *
 * Convolve the image with each filter (in parallel) and cluster the 
 * response vectors into textons.  Return the texton assignment for each
 * pixel as well as the textons themselves.
 * 
 * Cluster textons using the parallel K-means clusterer with the L2 distance
 * metric.  The number of textons (K) and maximum number of iterations may
 * be specified below.  Note that setting the number of iterations to zero
 * indicates unlimited iterations (in this case, K-means continues until
 * convergence).
 *
 * In order to speed up cluserting, users may specify that textons should be
 * produced by clustering only a randomly selected subsample of the filter
 * responses.
 */
matrix<unsigned long> lib_image::textons(
   const matrix<>&                                      m,
   const collection< matrix<> >&                        filters,
   auto_collection< matrix<>, array_list< matrix<> > >& textons,
   unsigned long                                        K,
   unsigned long                                        max_iter,
   double                                               subsampling)
{
   return lib_image::textons(
      m,
      filters,
      sample_clusterer< matrix<> >(
         kmeans::matrix_clusterer<>(K, max_iter, matrix_metrics<>::L2_metric()),
         subsampling,
         K
      ),
      textons
   );
}

/*
 * Compute textons as above, however specify a custom clustering routine.
 */
matrix<unsigned long> lib_image::textons(
   const matrix<>&                                      m,
   const collection< matrix<> >&                        filters,
   const centroid_clusterer< matrix<> >&                clusterer,
   auto_collection< matrix<>, array_list< matrix<> > >& textons)
{
   /* convolve image with filters */
   auto_collection< matrix<>, array_list< matrix<> > > responses =
      lib_image::filter(m, filters);
   /* cluster filter responses */
   return lib_image::cluster_filter_responses(*responses, clusterer, textons);
}

/***************************************************************************
 * Histogram gradient helper functions (2D).
 ***************************************************************************/

namespace {
/*
 * Construct region mask for circular disc of the given radius.
 */
matrix<bool> region_mask_disc(unsigned long r) {
   /* initialize matrix */
   unsigned long size = 2*r + 1;
   matrix<bool> mask(size, size);
   /* set values in disc to true */
   long radius = static_cast<long>(r);
   long r_sq = radius * radius;
   unsigned long ind = 0;
   for (long x = -radius; x <= radius; x++) {
      long x_sq = x * x;
      for (long y = -radius; y <= radius; y++) {
         /* check if index is within disc */
         long y_sq = y * y;
         if ((x_sq + y_sq) <= r_sq)
            mask[ind] = true;
         /* increment linear index */
         ind++;
      }
   }
   return mask;
}

/*
 * Construct region mask for annulus with the given inner and outer radii.
 */
matrix<bool> region_mask_annulus(unsigned long r_inner, unsigned long r_outer)
{
   /* initialize matrix */
   unsigned long size = 2*r_outer + 1;
   matrix<bool> mask(size, size);
   /* set values in annulus to true */
   long radius = static_cast<long>(r_outer);
   long r_outer_sq = radius * radius;
   long r_inner_sq = static_cast<long>(r_inner) * static_cast<long>(r_inner);
   unsigned long ind = 0;
   for (long x = -radius; x <= radius; x++) {
      long x_sq = x * x;
      for (long y = -radius; y <= radius; y++) {
         /* check if index is within circular disc */
         long y_sq = y * y;
         long dist_sq = x_sq + y_sq;
         if ((r_inner_sq < dist_sq) && (dist_sq <= r_outer_sq))
            mask[ind] = true;
         /* increment linear index */
         ind++;
      }
   }
   return mask;
}

/*
 * Construct weight matrix for circular disc of the given radius.
 */
matrix<> weight_matrix_disc(unsigned long r) {
   /* initialize matrix */
   unsigned long size = 2*r + 1;
   matrix<> weights(size, size);
   /* set values in disc to 1 */
   long radius = static_cast<long>(r);
   long r_sq = radius * radius;
   unsigned long ind = 0;
   for (long x = -radius; x <= radius; x++) {
      long x_sq = x * x;
      for (long y = -radius; y <= radius; y++) {
         /* check if index is within disc */
         long y_sq = y * y;
         if ((x_sq + y_sq) <= r_sq)
            weights[ind] = 1;
         /* increment linear index */
         ind++;
      }
   }
   return weights;
}

/*
 * Construct orientation slice lookup map.
 */
matrix<unsigned long> orientation_slice_map(
   unsigned long size_x, unsigned long size_y, unsigned long n_ori)
{
   /* initialize map */
   matrix<unsigned long> slice_map(size_x, size_y);
   /* compute orientation of each element from center */
   unsigned long ind = 0;
   double x = -static_cast<double>(size_x)/2;
   for (unsigned long n_x = 0; n_x < size_x; n_x++) {
      double y = -static_cast<double>(size_y)/2;
      for (unsigned long n_y = 0; n_y < size_y; n_y++) {
         /* compute orientation index */
         double ori = math::atan2(y, x) + M_PIl;
         unsigned long idx = static_cast<unsigned long>(
            math::floor(ori / M_PIl * static_cast<double>(n_ori))
         );
         if (idx >= (2*n_ori))
            idx = 2*n_ori - 1;
         slice_map[ind] = idx;
         /* increment position */
         ind++;
         y++;
      }
      /* increment x-coordinate */
      x++;
   }
   return slice_map;
}

/*
 * Compute convolution in place (for 1D matrices).
 */
void conv_in_place_1D(
   const matrix<>& m0,
   const matrix<>& m1,
   matrix<>& m)
{
   /* get size of each matrix */
   unsigned long size0 = m0.size();
   unsigned long size1 = m1.size();
   /* set dimensions for result matrix no larger than left input */
   unsigned long size = ((size0 > 0) && (size1 > 0)) ? (size0) : 0;
   /* set start position for result matrix no larger than left input */
   unsigned long pos_start = size1/2;
   /* initialize position in result */
   unsigned long pos = pos_start;
   for (unsigned long n = 0; n < size; n++) {
      /* compute range of offset */
      unsigned long offset_min = ((pos + 1) > size0) ? (pos + 1 - size0) : 0;
      unsigned long offset_max = (pos < size1) ? pos : (size1 - 1);
      /* multiply and add corresponing elements */
      unsigned long ind0 = pos - offset_min;
      unsigned long ind1 = offset_min;
      while (ind1 <= offset_max) {
         /* update result value */
         m[n] += m0[ind0] * m1[ind1];
         /* update linear positions */
         ind0--;
         ind1++;
      }
      /* update position */
      pos++;
   }
}

/*
 * Compute histograms and histogram differences at each location.
 */
void compute_hist_gradient_2D(
   const matrix<unsigned long>&                 labels,
   const matrix<>&                              weights,
   const matrix<unsigned long>&                 slice_map,
   const matrix<>&                              smoothing_kernel,
   array_list< matrix<> >&                      slice_hist, /* hist per slice */
   array_list< matrix<> >&                      gradients,  /* matrix per ori */
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* get number of orientations */
   unsigned long n_ori = gradients.size();
   /* get label matrix size */ 
   unsigned long size0_x = labels.size(0);
   unsigned long size0_y = labels.size(1);
   /* get window size */
   unsigned long size1_x = weights.size(0);
   unsigned long size1_y = weights.size(1);
   /* set start position for gradient matrices */
   unsigned long pos_start_x = size1_x/2;
   unsigned long pos_start_y = size1_y/2;
   unsigned long pos_bound_y = pos_start_y + size0_y;
   /* initialize position in result */
   unsigned long pos_x = pos_start_x;
   unsigned long pos_y = pos_start_y;
   /* compute initial range of offset_x */
   unsigned long offset_min_x =
      ((pos_x + 1) > size0_x) ? (pos_x + 1 - size0_x) : 0;
   unsigned long offset_max_x =
      (pos_x < size1_x) ? pos_x : (size1_x - 1);
   unsigned long ind0_start_x = (pos_x - offset_min_x) * size0_y;
   unsigned long ind1_start_x = (offset_min_x) * size1_y;
   unsigned long size = labels.size();
   /* determine whether to use smoothing kernel */
   bool use_smoothing = !(smoothing_kernel.is_empty());
   matrix<> temp_conv(slice_hist[0].dimensions());
   unsigned long size_hist = slice_hist[0].size();
   /* allocate half disc histograms */
   matrix<> hist_left(slice_hist[0].dimensions());
   matrix<> hist_right(slice_hist[0].dimensions());
   for (unsigned long n = 0; n < size; n++) {
      /* compute range of offset_y */
      unsigned long offset_min_y =
         ((pos_y + 1) > size0_y) ? (pos_y + 1 - size0_y) : 0;
      unsigned long offset_max_y =
         (pos_y < size1_y) ? pos_y : (size1_y - 1);
      unsigned long offset_range_y = offset_max_y - offset_min_y;
      /* initialize indices */
      unsigned long ind0 = ind0_start_x + (pos_y - offset_min_y);
      unsigned long ind1 = ind1_start_x + offset_min_y;
      /* clear histograms */
      for (unsigned long n_hist = 0; n_hist < 2*n_ori; n_hist++)
         slice_hist[n_hist].fill(0);
      /* update histograms */
      for (unsigned long o_x = offset_min_x; o_x <= offset_max_x; o_x++) {
         for (unsigned long o_y = offset_min_y; o_y < offset_max_y; o_y++) {
            /* update histogram value */
            (slice_hist[slice_map[ind1]])[labels[ind0]] += weights[ind1];
            /* update linear positions */
            ind0--;
            ind1++;
         }
         /* update last histogram value */
         (slice_hist[slice_map[ind1]])[labels[ind0]] += weights[ind1];
         /* update linear positions */
         ind0 = ind0 + offset_range_y - size0_y;
         ind1 = ind1 - offset_range_y + size1_y;
      }
      /* smooth bins */
      if (use_smoothing) {
         for (unsigned long o = 0; o < 2*n_ori; o++) {
            matrix<>& sh = slice_hist[o];
            temp_conv.fill(0);
            conv_in_place_1D(sh, smoothing_kernel, temp_conv);
            for (unsigned long nh = 0; nh < size_hist; nh++)
               sh[nh] = temp_conv[nh];
         }
      }
      /* L1 normalize bins */
      for (unsigned long o = 0; o < 2*n_ori; o++) {
         double sum_slice_hist = sum(slice_hist[o]);
         if (sum_slice_hist != 0)
            slice_hist[o] /= sum_slice_hist;
      }
      /* compute circular gradients - initialize histograms */
      hist_left.fill(0);
      hist_right.fill(0);
      for (unsigned long o = 0; o < n_ori; o++) {
         hist_left  += slice_hist[o];
         hist_right += slice_hist[o+n_ori];
      }
      /* compute circular gradients - spin the disc */
      for (unsigned long o = 0; o < n_ori; o++) {
         (gradients[o])[n] = f_dist(hist_left, hist_right);
         hist_left -= slice_hist[o];
         hist_left += slice_hist[o+n_ori];
         hist_right += slice_hist[o];
         hist_right -= slice_hist[o+n_ori];
      }
      /* update position */
      pos_y++;
      if (pos_y == pos_bound_y) {
         /* reset y position, increment x position */
         pos_y = pos_start_y;
         pos_x++;
         /* update range of offset_x */
         offset_min_x = ((pos_x + 1) > size0_x) ? (pos_x + 1 - size0_x) : 0;
         offset_max_x = (pos_x < size1_x) ? pos_x : (size1_x - 1);
         ind0_start_x = (pos_x - offset_min_x) * size0_y;
         ind1_start_x = (offset_min_x) * size1_y;
      }
   }
}

/*
 * Compute histograms and histogram differences at each location.
 */
void compute_hist_gradient_2D(
   const matrix< matrix<> >&                    histograms,
   const matrix<>&                              weights,
   const matrix<unsigned long>&                 slice_map,
   const matrix<>&                              smoothing_kernel,
   array_list< matrix<> >&                      slice_hist, /* hist per slice */
   array_list< matrix<> >&                      gradients,  /* matrix per ori */
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* get number of orientations */
   unsigned long n_ori = gradients.size();
   /* get label matrix size */ 
   unsigned long size0_x = histograms.size(0);
   unsigned long size0_y = histograms.size(1);
   /* get window size */
   unsigned long size1_x = weights.size(0);
   unsigned long size1_y = weights.size(1);
   /* set start position for gradient matrices */
   unsigned long pos_start_x = size1_x/2;
   unsigned long pos_start_y = size1_y/2;
   unsigned long pos_bound_y = pos_start_y + size0_y;
   /* initialize position in result */
   unsigned long pos_x = pos_start_x;
   unsigned long pos_y = pos_start_y;
   /* compute initial range of offset_x */
   unsigned long offset_min_x =
      ((pos_x + 1) > size0_x) ? (pos_x + 1 - size0_x) : 0;
   unsigned long offset_max_x =
      (pos_x < size1_x) ? pos_x : (size1_x - 1);
   unsigned long ind0_start_x = (pos_x - offset_min_x) * size0_y;
   unsigned long ind1_start_x = (offset_min_x) * size1_y;
   unsigned long size = histograms.size();
   /* determine whether to use smoothing kernel */
   bool use_smoothing = !(smoothing_kernel.is_empty());
   /* allocate half disc histograms */
   matrix<> hist_left(slice_hist[0].dimensions());
   matrix<> hist_right(slice_hist[0].dimensions());
   for (unsigned long n = 0; n < size; n++) {
      /* compute range of offset_y */
      unsigned long offset_min_y =
         ((pos_y + 1) > size0_y) ? (pos_y + 1 - size0_y) : 0;
      unsigned long offset_max_y =
         (pos_y < size1_y) ? pos_y : (size1_y - 1);
      unsigned long offset_range_y = offset_max_y - offset_min_y;
      /* initialize indices */
      unsigned long ind0 = ind0_start_x + (pos_y - offset_min_y);
      unsigned long ind1 = ind1_start_x + offset_min_y;
      /* clear histograms */
      for (unsigned long n_hist = 0; n_hist < 2*n_ori; n_hist++)
         slice_hist[n_hist].fill(0);
      /* update histograms */
      for (unsigned long o_x = offset_min_x; o_x <= offset_max_x; o_x++) {
         for (unsigned long o_y = offset_min_y; o_y < offset_max_y; o_y++) {
            /* update histogram value */
            (slice_hist[slice_map[ind1]]) += weights[ind1] * histograms[ind0];
            /* update linear positions */
            ind0--;
            ind1++;
         }
         /* update last histogram value */
         (slice_hist[slice_map[ind1]]) += weights[ind1] * histograms[ind0];
         /* update linear positions */
         ind0 = ind0 + offset_range_y - size0_y;
         ind1 = ind1 - offset_range_y + size1_y;
      }
      /* smooth bins */
      if (use_smoothing) {
         for (unsigned long o = 0; o < 2*n_ori; o++)
            slice_hist[o] = conv_crop(slice_hist[o], smoothing_kernel);
      }
      /* L1 normalize bins */
      for (unsigned long o = 0; o < 2*n_ori; o++) {
         double sum_slice_hist = sum(slice_hist[o]);
         if (sum_slice_hist != 0)
            slice_hist[o] /= sum_slice_hist;
      }
      /* compute circular gradients - initialize histograms */
      hist_left.fill(0);
      hist_right.fill(0);
      for (unsigned long o = 0; o < n_ori; o++) {
         hist_left  += slice_hist[o];
         hist_right += slice_hist[o+n_ori];
      }
      /* compute circular gradients - spin the disc */
      for (unsigned long o = 0; o < n_ori; o++) {
         (gradients[o])[n] = f_dist(hist_left, hist_right);
         hist_left -= slice_hist[o];
         hist_left += slice_hist[o+n_ori];
         hist_right += slice_hist[o];
         hist_right -= slice_hist[o+n_ori];
      }
      /* update position */
      pos_y++;
      if (pos_y == pos_bound_y) {
         /* reset y position, increment x position */
         pos_y = pos_start_y;
         pos_x++;
         /* update range of offset_x */
         offset_min_x = ((pos_x + 1) > size0_x) ? (pos_x + 1 - size0_x) : 0;
         offset_max_x = (pos_x < size1_x) ? pos_x : (size1_x - 1);
         ind0_start_x = (pos_x - offset_min_x) * size0_y;
         ind1_start_x = (offset_min_x) * size1_y;
      }
   }
}
} /* namespace */

/***************************************************************************
 * Difference of histograms (2D).
 ***************************************************************************/

/*
 * Compute the distance between histograms of label values in the disc
 * and annulus of the specified radii centered at each location in the
 * 2D matrix.  The disc is bounded by the inner radius while the annulus
 * lies between the inner and outer radius.
 *
 * Alternatively, instead of specifying label values at each point, the user
 * may specify a histogram at each point, in which case the histogram for
 * a region is the sum of the histograms at points in that region.
 *
 * The user may optionally specify a nonempty 1D smoothing kernel to 
 * convolve with histograms prior to computing the distance between them.
 *
 * The user may also optionally specify a custom functor for computing the
 * distance between histograms.
 */
matrix<> lib_image::hist_diff_2D(
   const matrix<unsigned long>&                 labels,
   unsigned long                                r_inner,
   unsigned long                                r_outer,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* construct region mask and weight matrix */
   matrix<bool> region_mask = region_mask_annulus(r_inner, r_outer);
   matrix<> weights = weight_matrix_disc(r_outer);
   /* compute histrogram difference */
   return lib_image::hist_diff_2D(
      labels, region_mask, weights, smoothing_kernel, f_dist
   );
}

matrix<> lib_image::hist_diff_2D(
   const matrix< matrix<> >&                    histograms,
   unsigned long                                r_inner,
   unsigned long                                r_outer,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* construct region mask and weight matrix */ 
   matrix<bool> region_mask = region_mask_annulus(r_inner, r_outer);
   matrix<> weights = weight_matrix_disc(r_outer);
   /* compute histrogram difference */
   return lib_image::hist_diff_2D(
      histograms, region_mask, weights, smoothing_kernel, f_dist
   );
}

/*
 * Compute the distance between histograms of label values in regions
 * centered at each location in the 2D matrix.
 *
 * The given 2D binary region mask (which must have odd dimensions) defines
 * the regions.  Each label adds the weight at the corresponding position 
 * in the weight matrix to its histogram bin in the respective region.
 *
 * If the user specifies a histogram at each point instead of a label, then 
 * the histogram for a region is a weighted sum of histograms.
 */
matrix<> lib_image::hist_diff_2D(
   const matrix<unsigned long>&                 labels,
   const matrix<bool>&                          region_mask,
   const matrix<>&                              weights,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* check arguments - labels */
   if (labels._dims.size() != 2)
      throw ex_invalid_argument("label matrix must be 2D");
   /* check arguments - region mask */
   if (region_mask._dims.size() != 2)
      throw ex_invalid_argument("region mask must be 2D");
   unsigned long r_size_x = region_mask._dims[0];
   unsigned long r_size_y = region_mask._dims[1];
   if (((r_size_x/2) == ((r_size_x+1)/2)) ||
       ((r_size_y/2) == ((r_size_y+1)/2)))
      throw ex_invalid_argument("dimensions of region mask must be odd");
   /* check arguments - weights */
   if (weights._dims.size() != 2)
      throw ex_invalid_argument("weight matrix must be 2D");
   unsigned long w_size_x = weights._dims[0];
   unsigned long w_size_y = weights._dims[1];
   if ((w_size_x != r_size_x) || (w_size_y != r_size_y))
      throw ex_invalid_argument(
         "dimensions of weight matrix must match dimensions of region mask"
      );
   /* initialize result histogram difference matrix */
   matrix<> h_diff(labels._dims);
   if (labels._size == 0)
      return h_diff;
   /* initialize collection containing result histogram difference matrix */
   array_list< matrix<> > gradients;
   gradients.add(h_diff);
   /* initialize slice histogram matrices */
   unsigned long hist_length = max(labels) + 1;
   matrix<> slice_hist0(1, hist_length);
   matrix<> slice_hist1(1, hist_length);
   array_list< matrix<> > slice_hist;
   slice_hist.add(slice_hist0);
   slice_hist.add(slice_hist1);
   /* convert region mask into a slice lookup map (with two slices) */
   matrix<unsigned long> slice_map(region_mask);
   /* compute histograms and histogram differences at each location */
   compute_hist_gradient_2D(
      labels, weights, slice_map, vector(smoothing_kernel),
      slice_hist, gradients, f_dist
   );
   return h_diff;
}

matrix<> lib_image::hist_diff_2D(
   const matrix< matrix<> >&                    histograms,
   const matrix<bool>&                          region_mask,
   const matrix<>&                              weights,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* check arguments - labels */
   if (histograms._dims.size() != 2)
      throw ex_invalid_argument("histogram matrix must be 2D");
   /* check arguments - region mask */
   if (region_mask._dims.size() != 2)
      throw ex_invalid_argument("region mask must be 2D");
   unsigned long r_size_x = region_mask._dims[0];
   unsigned long r_size_y = region_mask._dims[1];
   if (((r_size_x/2) == ((r_size_x+1)/2)) ||
       ((r_size_y/2) == ((r_size_y+1)/2)))
      throw ex_invalid_argument("dimensions of region mask must be odd");
   /* check arguments - weights */
   if (weights._dims.size() != 2)
      throw ex_invalid_argument("weight matrix must be 2D");
   unsigned long w_size_x = weights._dims[0];
   unsigned long w_size_y = weights._dims[1];
   if ((w_size_x != r_size_x) || (w_size_y != r_size_y))
      throw ex_invalid_argument(
         "dimensions of weight matrix must match dimensions of region mask"
      );
   /* initialize result histogram difference matrix */
   matrix<> h_diff(histograms._dims);
   if (histograms._size == 0)
      return h_diff;
   /* initialize collection containing result histogram difference matrix */
   array_list< matrix<> > gradients;
   gradients.add(h_diff);
   /* initialize slice histogram matrices */
   matrix<> slice_hist0(histograms[0]._dims);
   matrix<> slice_hist1(histograms[0]._dims);
   array_list< matrix<> > slice_hist;
   slice_hist.add(slice_hist0);
   slice_hist.add(slice_hist1);
   /* convert region mask into a slice lookup map (with two slices) */
   matrix<unsigned long> slice_map(region_mask);
   /* compute histograms and histogram differences at each location */
   compute_hist_gradient_2D(
      histograms, weights, slice_map, vector(smoothing_kernel),
      slice_hist, gradients, f_dist
   );
   return h_diff;
}

/***************************************************************************
 * Oriented gradient of histograms (2D).
 ***************************************************************************/

namespace {
/*
 * Runnable class for resizing matrices in parallel.
 */
class matrix_resizer : public runnable {
public:
   /*
    * Constructor.
    */
   matrix_resizer(
      array_list< matrix<> >&     m_arr,  /* matrices to resize */
      const array<unsigned long>& dims,   /* desired dimensions */
      unsigned long               start,  /* first index to process */
      unsigned long               end)    /* last index to process */
    : _m_arr(m_arr), _dims(dims), _start(start), _end(end) { }
   
   /*
    * Destructor.
    */
   virtual ~matrix_resizer() { /* do nothing */ }

   /*
    * Resize matrices.
    */
   virtual void run() {
      if ((thread::processors() > 1) && (_start < _end)) {
         /* split resizing task */
         unsigned long mid = (_start + _end)/2;
         matrix_resizer r_left(_m_arr, _dims, _start, mid);
         matrix_resizer r_right(_m_arr, _dims, mid+1, _end);
         child_thread::run(r_left, r_right);
      } else {
         /* resize matrices */
         for (unsigned long n = _start; n <= _end; n++)
            _m_arr[n].resize(_dims);
      }
   }
   
protected:
   array_list< matrix<> >&     _m_arr;    /* matrices to resize */
   const array<unsigned long>& _dims;     /* desired dimensions */
   unsigned long               _start;    /* first index to process */
   unsigned long               _end;      /* last index to process */
};
} /* namespace */

/*
 * Compute the distance between histograms of label values in oriented
 * half-dics of the specified radius centered at each location in the 2D
 * matrix.  Return one distance matrix per orientation.
 *
 * Alternatively, instead of specifying label values at each point, the user
 * may specify a histogram at each point, in which case the histogram for
 * a half-disc is the sum of the histograms at points in that half-disc.
 *
 * The half-disc orientations are k*pi/n for k in [0,n) where n is the 
 * number of orientation requested.
 *
 * The user may optionally specify a nonempty 1D smoothing kernel to 
 * convolve with histograms prior to computing the distance between them.
 *
 * The user may also optionally specify a custom functor for computing the
 * distance between histograms.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::hist_gradient_2D(
   const matrix<unsigned long>&                 labels,
   unsigned long                                r,
   unsigned long                                n_ori,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* construct weight matrix for circular disc */
   matrix<> weights = weight_matrix_disc(r);
   /* compute oriented gradient histograms */
   return lib_image::hist_gradient_2D(
      labels, weights, n_ori, smoothing_kernel, f_dist
   );
}

auto_collection< matrix<>, array_list< matrix<> > > lib_image::hist_gradient_2D(
   const matrix< matrix<> >&                    histograms,
   unsigned long                                r,
   unsigned long                                n_ori,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* construct weight matrix for circular disc */
   matrix<> weights = weight_matrix_disc(r);
   /* compute oriented gradient histograms */
   return lib_image::hist_gradient_2D(
      histograms, weights, n_ori, smoothing_kernel, f_dist
   );
}

/*
 * Compute the distance between histograms of label values in oriented
 * half-regions centered at each location in the 2D matrix.  Return one
 * distance matrix per orientation.
 *
 * The given 2D weight matrix (which must have odd dimensions) defines the
 * half-regions.  Each label adds the weight at the corresponding position
 * in this matrix to its histogram bin.
 *
 * If the user specifies a histogram at each point instead of a label, then 
 * the histogram for a half-region is a weighted sum of histograms.
 *
 * The above version of hist_gradient_2D which specifies a radius r is 
 * equivalent to calling this version with a (2*r+1) x (2*r+1) weight matrix
 * in which elements within a distance r from the center have weight one and
 * all other elements have weight zero.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::hist_gradient_2D(
   const matrix<unsigned long>&                 labels,
   const matrix<>&                              weights,
   unsigned long                                n_ori,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* check arguments - labels */
   if (labels._dims.size() != 2)
      throw ex_invalid_argument("label matrix must be 2D");
   /* check arguments - weights */
   if (weights._dims.size() != 2)
      throw ex_invalid_argument("weight matrix must be 2D");
   unsigned long w_size_x = weights._dims[0];
   unsigned long w_size_y = weights._dims[1];
   if (((w_size_x/2) == ((w_size_x+1)/2)) ||
       ((w_size_y/2) == ((w_size_y+1)/2)))
      throw ex_invalid_argument("dimensions of weight matrix must be odd");
   /* allocate result gradient matrices */
   auto_collection< matrix<>, array_list< matrix<> > > gradients(
      new array_list< matrix<> >()
   );
   /* check that result is nontrivial */
   if (n_ori == 0)
      return gradients;
   /* initialize result gradient matrices */
   for (unsigned long n = 0; n < n_ori; n++) {
      auto_ptr< matrix<> > m_ptr(new matrix<>());
      gradients->add(*m_ptr);
      m_ptr.release();
   }
   {
      matrix_resizer m_resizer(*gradients, labels._dims, 0, n_ori-1);
      m_resizer.run();
   }
   /* check that result is nonempty */
   if (labels._size == 0)
      return gradients;
   /* allocate matrices to hold histograms of each slice */
   auto_collection< matrix<>, array_list< matrix<> > > slice_hist(
      new array_list< matrix<> >()
   );
   /* initialize slice histogram matrices */
   for (unsigned long n = 0; n < 2*n_ori; n++) {
      auto_ptr< matrix<> > m_ptr(new matrix<>());
      slice_hist->add(*m_ptr);
      m_ptr.release();
   }
   {
      unsigned long hist_length = max(labels) + 1;
      array<unsigned long> slice_hist_dims(1, hist_length);
      matrix_resizer m_resizer(*slice_hist, slice_hist_dims, 0, 2*n_ori-1);
      m_resizer.run();
   }
   /* build orientation slice lookup map */
   matrix<unsigned long> slice_map = orientation_slice_map(
      w_size_x, w_size_y, n_ori
   );
   /* compute histograms and histogram differences at each location */
   compute_hist_gradient_2D(
      labels, weights, slice_map, vector(smoothing_kernel),
      *slice_hist, *gradients, f_dist
   );
   return gradients;
}

auto_collection< matrix<>, array_list< matrix<> > > lib_image::hist_gradient_2D(
   const matrix< matrix<> >&                    histograms,
   const matrix<>&                              weights,
   unsigned long                                n_ori,
   const matrix<>&                              smoothing_kernel,
   const distanceable_functor<matrix<>,double>& f_dist)
{
   /* check arguments - histograms */
   if (histograms._dims.size() != 2)
      throw ex_invalid_argument("histogram matrix must be 2D");
   /* check arguments - weights */
   if (weights._dims.size() != 2)
      throw ex_invalid_argument("weight matrix must be 2D");
   unsigned long w_size_x = weights._dims[0];
   unsigned long w_size_y = weights._dims[1];
   if (((w_size_x/2) == ((w_size_x+1)/2)) ||
       ((w_size_y/2) == ((w_size_y+1)/2)))
      throw ex_invalid_argument("dimensions of weight matrix must be odd");
   /* allocate result gradient matrices */
   auto_collection< matrix<>, array_list< matrix<> > > gradients(
      new array_list< matrix<> >()
   );
   /* check that result is nontrivial */
   if (n_ori == 0)
      return gradients;
   /* initialize result gradient matrices */
   for (unsigned long n = 0; n < n_ori; n++) {
      auto_ptr< matrix<> > m_ptr(new matrix<>());
      gradients->add(*m_ptr);
      m_ptr.release();
   }
   {
      matrix_resizer m_resizer(*gradients, histograms._dims, 0, n_ori-1);
      m_resizer.run();
   }
   /* check that result is nonempty */
   if (histograms._size == 0)
      return gradients;
   /* allocate matrices to hold histograms of each slice */
   auto_collection< matrix<>, array_list< matrix<> > > slice_hist(
      new array_list< matrix<> >()
   );
   /* initialize slice histogram matrices */
   for (unsigned long n = 0; n < 2*n_ori; n++) {
      auto_ptr< matrix<> > m_ptr(new matrix<>());
      slice_hist->add(*m_ptr);
      m_ptr.release();
   }
   {
      array<unsigned long> slice_hist_dims(histograms[0]._dims);
      matrix_resizer m_resizer(*slice_hist, slice_hist_dims, 0, 2*n_ori-1);
      m_resizer.run();
   }
   /* build orientation slice lookup map */
   matrix<unsigned long> slice_map = orientation_slice_map(
      w_size_x, w_size_y, n_ori
   );
   /* compute histograms and histogram differences at each location */
   compute_hist_gradient_2D(
      histograms, weights, slice_map, vector(smoothing_kernel),
      *slice_hist, *gradients, f_dist
   );
   return gradients;
}

/*
 * Combine oriented gradients of histograms into a single gradient map.
 *
 * Each location is assigned gradient strength eqaul to the maximum gradient
 * value over all orientations.  The dominant orientation at each location
 * is also returned.
 *
 * Unless otherwise specified, the orientations corresponding to the input
 * collection of n oriented gradient maps are assumed to be k*pi/n for
 * integer k in [0,n).
 */
void lib_image::combine_hist_gradients(
   const collection< matrix<> >& gradients,
   auto_ptr< matrix<> >&         gradient,
   auto_ptr< matrix<> >&         gradient_ori)
{
   unsigned long n_ori = gradients.size();
   array<double> oris = lib_image::standard_filter_orientations(n_ori);
   lib_image::combine_hist_gradients(gradients, oris, gradient, gradient_ori);
}

void lib_image::combine_hist_gradients(
   const collection< matrix<> >& gradients,
   const array<double>&          oris,
   auto_ptr< matrix<> >&         gradient,
   auto_ptr< matrix<> >&         gradient_ori)
{
   /* check arguments - number of orientations */
   unsigned long n_ori = oris.size();
   if (gradients.size() != n_ori)
      throw ex_invalid_argument(
         "number of gradient matrices different from number of orientations"
      );
   /* initialize gradient map */
   gradient.reset(new matrix<>());
   gradient_ori.reset(new matrix<>());
   /* check if nontrivial result */
   if (n_ori > 0) {
      /* check that all gradient matrices have the same size */
      array_list< matrix<> > gradients_arr(gradients);
      array<unsigned long> dims = gradients_arr[0]._dims;
      for (unsigned long n = 0; n < n_ori; n++)
         matrix<>::assert_dims_equal(dims, gradients_arr[n]._dims);
      /* compute indices of maximum gradient values */
      matrix<unsigned long> index(dims);
      if (index._size > 0) {
         matrix_max_selector m_select(gradients_arr, 0, index._size-1, index);
         m_select.run();
      }
      *gradient     = matrix_max_selector::select_values(gradients_arr, index);
      *gradient_ori = matrix_max_selector::select_values(oris, index);
   }
}

/***************************************************************************
 * Smoothing.
 ***************************************************************************/

/*
 * Compute "needleness" at each point from derivatives.
 * 
 * The "needleness" for a signal f(x) is defined as 
 *
 *             f(x) * (-f''(x)) / (abs(f'(x)) + epsilon)
 *
 * This function is useful for smoothing and merging double peaks.
 */
matrix<> lib_image::needleness(
   const matrix<>& f0, const matrix<>& f1, const matrix<>& f2, double epsilon)
{
   /* check arguments */
   matrix<>::assert_dims_equal(f0._dims, f1._dims);
   matrix<>::assert_dims_equal(f0._dims, f2._dims);
   /* compute needleness */
   matrix<> f(f0._dims);
   for (unsigned long n = 0; n < f._size; n++) {
      f._data[n] =
         f0._data[n] * (-(f2._data[n])) / (math::abs(f1._data[n]) + epsilon);
   }
   return f;
}

/*
 * Smooth and sharpen 2D oriented gradients using gaussian derivative 
 * filters and the "needleness" function.
 *
 * Smoothing is done along the gradient directions, which are assumed to be
 * k*pi/n for k in [0,n) where n is the number of orientated gradients.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::smooth_grad_2D(
   const collection< matrix<> >& gradients,
   double                        sigma,
   double                        epsilon)
{
   unsigned long n_ori = gradients.size();
   array<double> oris = lib_image::standard_filter_orientations(n_ori);
   return lib_image::smooth_grad_2D(gradients, oris, sigma, epsilon);
}

/*
 * Smooth and sharpen 2D oriented gradients using gaussian derivative 
 * filters and the "needleness" function.
 *
 * Specify the orientations along which to smooth.
 */
auto_collection< matrix<>, array_list< matrix<> > > lib_image::smooth_grad_2D(
   const collection< matrix<> >& gradients,
   const array<double>&          oris,
   double                        sigma,
   double                        epsilon)
{
   /* get number of orientations */
   array_list< matrix<> > gradients_arr(gradients);
   unsigned long n_ori = gradients_arr.size();
   /* smooth each orientation */
   auto_collection< matrix<>, array_list< matrix<> > > f_smoothed(
      new array_list< matrix<> >()
   );
   for (unsigned long n = 0; n < n_ori; n++) {
      /* create oriented gaussian smoothing kernels */
      matrix<> g0 = lib_image::gaussian_2D(sigma, sigma, oris[n], 0);
      matrix<> g1 = lib_image::gaussian_2D(sigma, sigma, oris[n], 1);
      matrix<> g2 = lib_image::gaussian_2D(sigma, sigma, oris[n], 2);
      /* convolve with smoothing kernels */
      matrix<>& f = gradients_arr[n];
      matrix<> f0 = conv_crop(f, g0);
      matrix<> f1 = conv_crop(f, g1);
      matrix<> f2 = conv_crop(f, g2);
      /* take positive component of f0 */
      for (unsigned long ind = 0; ind < f0._size; ind++)
         f0._data[ind] = (f0._data[ind] > 0) ? (f0._data[ind]) : 0;
      /* take negative component of f2 */
      for (unsigned long ind = 0; ind < f2._size; ind++)
         f2._data[ind] = (f2._data[ind] < 0) ? (f2._data[ind]) : 0;
      /* compute needleness */
      auto_ptr< matrix<> > f_result(
         new matrix<>(lib_image::needleness(f0, f1, f2, epsilon))
      );
      f_smoothed->add(*f_result);
      f_result.release();
   }
   return f_smoothed;
}

/***************************************************************************
 * Non-max suppression.
 ***************************************************************************/

/*
 * Non-max suppression.
 *
 * Perform non-max suppression in the local neighborhood of each matrix 
 * element.  Each local neighborhood is a block of 3^d elements where d is 
 * the dimensionality of the matrix.
 *
 * Optionally, a matrix of size 3^d can be used as a mask to determine
 * to which elements of the local neighborhood the center element should 
 * be compared during non-max suppression.
 *
 * Boundary elements are not considered as candidates for a local maximum.
 *
 * Non-max elements are assigned a value of zero.
 *
 * NOTE: The original matrix must be nonnegative.
 */
matrix<> lib_image::nonmax(const matrix<>& m) {
   /* create default mask for full 3^d neighborhood */
   unsigned long n_dims = m._dims.size();
   array<unsigned long> mask_dims(n_dims, 3);
   matrix<bool> mask(mask_dims, true);
   if (n_dims > 0)
      mask._data[mask._size/2] = false;   /* don't compare center to itself */
   /* perform non-max suppression */
   return lib_image::nonmax(m, mask);
}

matrix<> lib_image::nonmax(const matrix<>& m, const matrix<bool>& mask) {
   /* check that mask is valid */
   unsigned long n_dims = m._dims.size();
   bool mask_valid = (n_dims == mask._dims.size());
   for (unsigned long n = 0; ((n < n_dims) && (mask_valid)); n++)
      mask_valid = (mask._dims[n] == 3);
   if (!mask_valid)
      throw ex_invalid_argument(
         "mask for d-dimensional matrix must be size 3^d"
      );
   /* initialize result matrix */
   matrix<> nmax(m._dims);
   /* check that result is nontrivial */
   if (nmax._size > 0) {
      /* compute offsets of neighbors */
      array<long> neighbor_offsets = matrix<>::linear_neighbor_offsets(m._dims);
      /* grab offsets specified by mask */
      array<unsigned long> used_offsets = mask.find();
      array<long> offsets = neighbor_offsets.subarray(used_offsets);
      unsigned long n_offsets = offsets.size();
      /* compute indices of interior elements */
      array<unsigned long> start(n_dims, 1);
      array<unsigned long> end(n_dims);
      for (unsigned long n = 0; n < n_dims; n++)
         end[n] = m._dims[n] - 1;
      array<unsigned long> inds = matrix<>::linear_indices(m._dims, start, end);
      /* perform non-max suppression at each interior element */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         /* get index */
         unsigned long ind_cntr = inds[n];
         /* get element value */
         double v = m._data[ind_cntr];
         if (v < 0)
            throw ex_invalid_argument("matrix must be nonnegative");
         /* compare to neighbors */
         bool is_max = true;
         for (unsigned long n_off = 0;
              ((n_off < n_offsets) && (is_max)); n_off++)
         {
            long ind = static_cast<long>(ind_cntr) + offsets[n_off];
            /* check if non-max */
            if (v <= m._data[ind])
               is_max = false;
         }
         /* suppress non-max */
         if (is_max)
            nmax._data[ind_cntr] = v;
      }
   }
   return nmax;
}

/*
 * Non-max suppression (along specified dimension).
 *
 * Perform non-max suppression along the specified dimension of the matrix.
 *
 * Boundary elements along this dimension are only considered as candidates
 * for a local maximum if the flag to allow boundary candidates is set.
 *
 * Non-max elements are assigned a value of zero.
 *
 * NOTE: The original matrix must be nonnegative.
 */
matrix<> lib_image::nonmax(
   const matrix<>& m, unsigned long d, bool allow_boundary)
{
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* initialize result matrix */
   matrix<> nmax(m._dims);
   /* check that result is nontrivial */
   if (nmax._size > 0) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         unsigned long ind_next = ind + step_size;
         /* check nonnegativitiy */
         if (m._data[ind] < 0)
            throw ex_invalid_argument("matrix must be nonnegative");
         /* check if neighbors exist */
         if (n_steps > 1) {
            /* check first element */
            if (allow_boundary) {
               if (m._data[ind] > m._data[ind_next])
                  nmax._data[ind] = m._data[ind];
            }
            /* check interior elements */
            for (unsigned long n_step = 1; n_step < (n_steps-1); n_step++) {
               unsigned long ind_prev = ind;
               ind = ind_next;
               ind_next += step_size;
               if ((m._data[ind] > m._data[ind_prev]) &&
                   (m._data[ind] > m._data[ind_next]))
                  nmax._data[ind] = m._data[ind];
            }
            /* check last element */
            if (allow_boundary) {
               if (m._data[ind_next] > m._data[ind])
                  nmax._data[ind_next] = m._data[ind_next];
            }
         } else if (allow_boundary) {
            nmax._data[ind] = m._data[ind];
         }
      }
   }
   return nmax;
}

/*
 * Oriented non-max suppression (2D).
 *
 * Perform non-max suppression orthogonal to the specified orientation on
 * the given 2D matrix using linear interpolation in a 3x3 neighborhood.
 *
 * A local maximum must be greater than the interpolated values of its 
 * adjacent elements along the direction orthogonal to this orientation.
 *
 * Elements which have a neighbor on only one side of this direction are 
 * only considered as candidates for a local maximum if the flag to allow 
 * boundary candidates is set.
 *
 * The same orientation angle may be specified for all elements, or the 
 * orientation may be specified for each matrix element.
 *
 * If an orientation is specified per element, then the elements themselves
 * may optionally be treated as oriented vectors by specifying a value less 
 * than pi/2 for the orientation tolerance.  In this case, neighboring 
 * vectors are projected along a line in the local orientation direction and
 * the lengths of the projections are used in determining local maxima.
 * When projecting, the orientation tolerance is subtracted from the true
 * angle between the vector and the line (with a result less than zero
 * causing the length of the projection to be the length of the vector).
 *
 * Non-max elements are assigned a value of zero.
 *
 * NOTE: The original matrix must be nonnegative.
 */
matrix<> lib_image::nonmax_oriented_2D(
   const matrix<>& m, double ori, bool allow_boundary)
{
   /* create per element orientation matrix */
   matrix<> m_ori(m._dims, ori);
   /* perform non-max suppression */
   return lib_image::nonmax_oriented_2D(m, m_ori, M_PI_2l, allow_boundary);
}

matrix<> lib_image::nonmax_oriented_2D(
   const matrix<>& m, const matrix<>& m_ori, double o_tol, bool allow_boundary)
{
   /* check arguments - matrix size */
   if (m._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   matrix<>::assert_dims_equal(m._dims, m_ori._dims);
   /* check arguments - orientation tolerance */
   if (o_tol < 0)
      throw ex_invalid_argument("orientation tolerance must be nonnegative");
   /* get matrix size */
   unsigned long size_x = m._dims[0];
   unsigned long size_y = m._dims[1];
   /* intialize result matrix */
   matrix<> nmax(size_x, size_y);
   /* perform oriented non-max suppression at each element */
   unsigned long n = 0;
   for (unsigned long x = 0; x < size_x; x++) {
      for (unsigned long y = 0; y < size_y; y++) {
         /* compute direction (in [0,pi)) along which to suppress */
         double ori = m_ori._data[n];
         double theta = ori + M_PI_2l;
         theta -= math::floor(theta/M_PIl) * M_PIl;
         /* check nonnegativity */
         double v = m._data[n];
         if (v < 0)
            throw ex_invalid_argument("matrix must be nonnegative");
         /* initialize indices of values in local neighborhood */
         unsigned long ind0a = 0, ind0b = 0, ind1a = 0, ind1b = 0;
         /* initialize distance weighting */
         double d = 0;
         /* initialize boundary flags */
         bool valid0 = false, valid1 = false;
         /* compute interpolation indicies */
         if (theta == 0) {
            valid0 = (x > 0); valid1 = (x < (size_x-1));
            if (valid0) { ind0a = n-size_y; ind0b = ind0a; }
            if (valid1) { ind1a = n+size_y; ind1b = ind1a; }
         } else if (theta < M_PI_4l) {
            d = math::tan(theta);
            valid0 = ((x > 0) && (y > 0));
            valid1 = ((x < (size_x-1)) && (y < (size_y-1)));
            if (valid0) { ind0a = n-size_y; ind0b = ind0a-1; }
            if (valid1) { ind1a = n+size_y; ind1b = ind1a+1; }
         } else if (theta < M_PI_2l) {
            d = math::tan(M_PI_2l - theta);
            valid0 = ((x > 0) && (y > 0));
            valid1 = ((x < (size_x-1)) && (y < (size_y-1)));
            if (valid0) { ind0a = n-1; ind0b = ind0a-size_y; }
            if (valid1) { ind1a = n+1; ind1b = ind1a+size_y; }
         } else if (theta == M_PI_2l) {
            valid0 = (y > 0); valid1 = (y < (size_y-1));
            if (valid0) { ind0a = n-1; ind0b = ind0a; }
            if (valid1) { ind1a = n+1; ind1b = ind1a; }
         } else if (theta < (3.0*M_PI_4l)) {
            d = math::tan(theta - M_PI_2l);
            valid0 = ((x < (size_x-1)) && (y > 0));
            valid1 = ((x > 0) && (y < (size_y-1)));
            if (valid0) { ind0a = n-1; ind0b = ind0a+size_y; }
            if (valid1) { ind1a = n+1; ind1b = ind1a-size_y; }
         } else /* (theta < M_PIl) */ {
            d = math::tan(M_PIl - theta);
            valid0 = ((x < (size_x-1)) && (y > 0));
            valid1 = ((x > 0) && (y < (size_y-1)));
            if (valid0) { ind0a = n+size_y; ind0b = ind0a-1; }
            if (valid1) { ind1a = n-size_y; ind1b = ind1a+1; }
         }
         /* check boundary conditions */
         if (allow_boundary || (valid0 && valid1)) {
            /* initialize values in local neighborhood */
            double v0a = 0,   v0b = 0,   v1a = 0,   v1b = 0;
            /* initialize orientations in local neighborhood */
            double ori0a = 0, ori0b = 0, ori1a = 0, ori1b = 0;
            /* grab values and orientations */
            if (valid0) {
               v0a = m._data[ind0a];
               v0b = m._data[ind0b];
               ori0a = m_ori._data[ind0a] - ori;
               ori0b = m_ori._data[ind0b] - ori;
            }
            if (valid1) {
               v1a = m._data[ind1a];
               v1b = m._data[ind1b];
               ori1a = m_ori._data[ind1a] - ori;
               ori1b = m_ori._data[ind1b] - ori;
            }
            /* place orientation difference in [0,pi/2) range */
            ori0a -= math::floor(ori0a/(2*M_PIl)) * (2*M_PIl);
            ori0b -= math::floor(ori0b/(2*M_PIl)) * (2*M_PIl);
            ori1a -= math::floor(ori1a/(2*M_PIl)) * (2*M_PIl);
            ori1b -= math::floor(ori1b/(2*M_PIl)) * (2*M_PIl);
            if (ori0a >= M_PIl) { ori0a = 2*M_PIl - ori0a; }
            if (ori0b >= M_PIl) { ori0b = 2*M_PIl - ori0b; }
            if (ori1a >= M_PIl) { ori1a = 2*M_PIl - ori1a; }
            if (ori1b >= M_PIl) { ori1b = 2*M_PIl - ori1b; }
            if (ori0a >= M_PI_2l) { ori0a = M_PIl - ori0a; }
            if (ori0b >= M_PI_2l) { ori0b = M_PIl - ori0b; }
            if (ori1a >= M_PI_2l) { ori1a = M_PIl - ori1a; }
            if (ori1b >= M_PI_2l) { ori1b = M_PIl - ori1b; }
            /* correct orientation difference by tolerance */
            ori0a = (ori0a <= o_tol) ? 0 : (ori0a - o_tol);
            ori0b = (ori0b <= o_tol) ? 0 : (ori0b - o_tol);
            ori1a = (ori1a <= o_tol) ? 0 : (ori1a - o_tol);
            ori1b = (ori1b <= o_tol) ? 0 : (ori1b - o_tol);
            /* interpolate */
            double v0 =
               (1.0-d)*v0a*math::cos(ori0a) + d*v0b*math::cos(ori0b);
            double v1 =
               (1.0-d)*v1a*math::cos(ori1a) + d*v1b*math::cos(ori1b);
            /* suppress non-max */
            if ((v > v0) && (v > v1))
               nmax._data[n] = v;
         }
         /* increment linear coordinate */
         n++;
      }
   }
   return nmax;
}

/***************************************************************************
 * Neighbor lookup and comparison (2D).
 ***************************************************************************/

namespace {
/*
 * Compute which neighbors are present in a 2D matrix.
 */
void neighbor_exists_2D(
   unsigned long size_x, unsigned long size_y,  /* matrix size */
   unsigned long x,      unsigned long y,       /* position */
   bool& n,  bool& s,  bool& e,  bool& w,
   bool& ne, bool& nw, bool& se, bool& sw)
{
   n  = ((y+1) < size_y);  s = (y > 0);  e  = ((x+1) < size_x);  w  = (x > 0);
   ne = n && e;            nw = n && w;  se = s && e;            sw = s && w;
}
 
/*
 * Compute indices of neighbors in a 2D matrix.
 */
void neighbor_indices_2D(
   unsigned long size_y,                  /* matrix size in y-direction */
   unsigned long ind,                     /* linear position */
   unsigned long& n,  unsigned long& s,  unsigned long& e,  unsigned long& w,
   unsigned long& ne, unsigned long& nw, unsigned long& se, unsigned long& sw)
{
   n  = ind + 1;  s  = ind - 1;  e  = ind + size_y;  w  = ind - size_y;
   ne = e + 1;    nw = w + 1;    se = e - 1;         sw = w - 1;
}

void neighbor_indices_2D(
   unsigned long size_y,                  /* matrix size in y-direction */
   unsigned long x,   unsigned long y,    /* position */
   unsigned long& n,  unsigned long& s,  unsigned long& e,  unsigned long& w,
   unsigned long& ne, unsigned long& nw, unsigned long& se, unsigned long& sw)
{
   neighbor_indices_2D(
      size_y, (x * size_y + y), n, s, e, w, ne, nw, se, sw
   );
}

/*
 * Lookup values of neighbors in a 2D matrix.
 * Return zero for any neighbors which do not exist.
 */
template <typename T>
void neighbor_values_2D(
   const T* m, unsigned long size_x, unsigned long size_y, /* matrix */
   unsigned long x, unsigned long y,                       /* position */
   unsigned long ind,                                      /* linear position */
   T& n, T& s, T& e, T& w, T& ne, T& nw, T& se, T& sw)
{
   /* compute which neighbors are present */ 
   bool has_n  = false, has_s  = false, has_e  = false, has_w  = false;
   bool has_ne = false, has_nw = false, has_se = false, has_sw = false;
   neighbor_exists_2D(
      size_x, size_y, x, y,
      has_n, has_s, has_e, has_w, has_ne, has_nw, has_se, has_sw
   );
   /* compute indices of neighbors */
   unsigned long ind_n  = 0, ind_s  = 0, ind_e  = 0, ind_w  = 0;
   unsigned long ind_ne = 0, ind_nw = 0, ind_se = 0, ind_sw = 0;
   neighbor_indices_2D(
      size_y, ind,
      ind_n, ind_s, ind_e, ind_w, ind_ne, ind_nw, ind_se, ind_sw
   );
   /* get values of neighbors */
   n  = (has_n  ? m[ind_n]  : T());  s  = (has_s  ? m[ind_s]  : T());
   e  = (has_e  ? m[ind_e]  : T());  w  = (has_w  ? m[ind_w]  : T());
   ne = (has_ne ? m[ind_ne] : T());  nw = (has_nw ? m[ind_nw] : T());
   se = (has_se ? m[ind_se] : T());  sw = (has_sw ? m[ind_sw] : T());
}

template <typename T>
void neighbor_values_2D(
   const T* m, unsigned long size_x, unsigned long size_y, /* matrix */
   unsigned long x, unsigned long y,                       /* position */
   T& n, T& s, T& e, T& w, T& ne, T& nw, T& se, T& sw)
{
   neighbor_values_2D(
      m, size_x, size_y, x, y, (x * size_y + y), n, s, e, w, ne, nw, se, sw
   );
}

/*
 * Compare values of neighbors in a 2D matrix to the given value.
 * Neighbors which do not exist are assumed to have a value of zero.
 */
template <typename T>
void neighbor_values_compare_2D(
   const T* m, unsigned long size_x, unsigned long size_y, /* matrix */
   unsigned long x, unsigned long y,                       /* position */
   unsigned long ind,                                      /* linear position */
   const T& v,                                             /* value */
   bool& n,  bool& s,  bool& e,  bool& w,
   bool& ne, bool& nw, bool& se, bool& sw)
{
   /* get values of neighbors */
   T v_n  = T(), v_s  = T(), v_e  = T(), v_w  = T();
   T v_ne = T(), v_nw = T(), v_se = T(), v_sw = T();
   neighbor_values_2D(
      m, size_x, size_y, x, y, ind, v_n, v_s, v_e, v_w, v_ne, v_nw, v_se, v_sw
   );
   /* compare */
   n  = (v_n  == v);  s  = (v_s  == v);  e  = (v_e  == v);  w  = (v_w  == v);
   ne = (v_ne == v);  nw = (v_nw == v);  se = (v_se == v);  sw = (v_sw == v);
}

template <typename T>
void neighbor_values_compare_2D(
   const T* m, unsigned long size_x, unsigned long size_y, /* matrix */
   unsigned long x, unsigned long y,                       /* position */
   const T& v,                                             /* value */
   bool& n,  bool& s,  bool& e,  bool& w,
   bool& ne, bool& nw, bool& se, bool& sw)
{
   neighbor_values_compare_2D(
      m, size_x, size_y, x, y, (x * size_y + y), v, n, s, e, w, ne, nw, se, sw
   );
}
} /* namespace */

/***************************************************************************
 * Local topology analysis (2D).
 ***************************************************************************/

namespace {
/*
 * Determine whether the center element is part of a thick region (a 2x2 or
 * larger block) given flags indicating which elements are present in the 
 * local 3x3 neighborhood.
 */
bool is_thick_region_2D(
   bool n, bool s, bool e, bool w, bool ne, bool nw, bool se, bool sw)
{
   return (
      (nw && n && w) || (n && ne && e) || (w && sw && s) || (e && s && se)
   );
}

/*
 * Determine whether removal of the center element is a valid erosion
 * operation.  Erosion is only possible if the center element is a local 
 * extrema on the boundary of a thick region.
 */
bool is_erodable_2D(
   bool n, bool s, bool e, bool w, bool ne, bool nw, bool se, bool sw)
{
   if (nw && n && w) {   
      return ((!e) && (!s) && (!se));
   } else if (n && ne && e) {
      return ((!w) && (!s) && (!sw));
   } else if (w && sw && s) {
      return ((!e) && (!n) && (!ne));
   } else if (e && s && se) {
      return ((!w) && (!n) && (!nw));
   } else {
      return false;
   }
}

/*
 * Determine whether removal of the center element preserves structure (does
 * not introduce holes, delete components, disconnect components, or shrink 
 * components) given flags indicating which elements are present in the local
 * 3x3 neighborhood.
 */
bool is_removable_2D(
   bool n, bool s, bool e, bool w, bool ne, bool nw, bool se, bool sw)
{
   /* compute number of adjacent neighbors */
   unsigned long n_adjacent = 
      (n  ? 1 : 0) + (s  ? 1 : 0) +
      (e  ? 1 : 0) + (w  ? 1 : 0);
   /* compute number of diagonal neighbors */
   unsigned long n_diagonal =
      (ne ? 1 : 0) + (nw ? 1 : 0) +
      (se ? 1 : 0) + (sw ? 1 : 0);
   /* compute total number of neighbors */
   unsigned long nn = n_adjacent + n_diagonal;
   /* check removability */
   if (nn <= 2) {
      /* removing center will shrink or disconnect component */
      return false;
   } else if (nn == 3) {
      /* can only remove if neighbors are adjacently connected */
      return (
         (ne && n && nw) || (n && ne && e) ||
         (se && s && sw) || (n && nw && w) ||
         (ne && e && se) || (s && se && e) ||
         (nw && w && sw) || (s && sw && w)
      );
   } else if (n_adjacent == 4) {
      /* removing center will create a hole */
      return false;
   } else if (is_thick_region_2D(n, s, e, w, ne, nw, se, sw)) {
      /* can merge center with neighbors if the neighbors are connected */
      unsigned long n_connect =
         ((n && (ne || e)) ? 1 : 0) +
         ((n && (nw || w)) ? 1 : 0) + 
         ((s && (se || e)) ? 1 : 0) + 
         ((s && (sw || w)) ? 1 : 0) + 
         ((e && ne) ? 1 : 0) + 
         ((e && se) ? 1 : 0) + 
         ((w && nw) ? 1 : 0) + 
         ((w && sw) ? 1 : 0);
      return ((n_connect+1) >= nn);
   } else {
      return false;
   }
}

/*
 * Determine whether the center element must be a vertex when elements are
 * grouped into contours given flags indicating which elements are present
 * in the local 3x3 neighborhood.
 */
bool is_vertex_2D(
   bool n, bool s, bool e, bool w, bool ne, bool nw, bool se, bool sw)
{
   /* compute number of adjacent neighbors */
   unsigned long n_adjacent = 
      (n  ? 1 : 0) + (s  ? 1 : 0) +
      (e  ? 1 : 0) + (w  ? 1 : 0);
   /* compute number of diagonal neighbors */
   unsigned long n_diagonal =
      (ne ? 1 : 0) + (nw ? 1 : 0) +
      (se ? 1 : 0) + (sw ? 1 : 0);
   /* compute total number of neighbors */
   unsigned long nn = n_adjacent + n_diagonal;
   /* check vertex conditions */
   if (nn <= 1) {
      /* start/end vertex - only 0 or 1 neighbors */
      return true;
   } else if (nn == 2) {
      /* start/end vertex if neighbors are adjacent */
      return (
         (n && (ne || nw)) ||
         (s && (se || sw)) ||
         (e && (ne || se)) ||
         (w && (nw || sw))
      );
   } else if ((nn == 3) || (nn == 4)) {
      if ((nw && n && w) || (n && ne && e) || 
          (w && sw && s) || (e && s && se)) {
         /* block of 4 (including center) */
         return true;
      } else if 
         ((nw && n && ne) || (sw && s && se) || 
          (nw && w && sw) || (ne && e && se)) {
         /* check for vertex hanging from three collinear points */
         return (nn == 3);
      } else if (n_adjacent >= 3) {
         /* junction with adjacent neighbors */
         return true;
      } else if (n_diagonal >= 3) {
         /* junction with diagonal neighbors */
         return true;
      } else if
         ((n && se && sw) || (s && ne && nw) ||
          (w && se && ne) || (e && nw && sw)) {
         /* junction with 1 adjacent, 2 diagonal */
         return true;
      } else if
         ((nw && s && e) || (ne && s && w) ||
          (sw && n && e) || (se && n && w)) {
         /* junction with 2 adjacent, 1 diagonal */
         return true;
      } else {
         return false;
      }
   } else {
      /* check not in middle of H or truncated H shape */
      return (!((n_adjacent == 2) && ((n && s) || (e && w))));
   }
}
} /* namespace */

/***************************************************************************
 * Skeletonization (2D).
 ***************************************************************************/

namespace {
/*
 * Skeleton element.
 */
class skeleton_element {
public:
   /*
    * Constructor.
    */
   explicit skeleton_element(
      unsigned long x, unsigned long y, double v, bool is_erode, bool is_rem)
    : x(x), y(y), value(v), is_erodable(is_erode), is_removable(is_rem) { }

   /*
    * Copy constructor.
    */
   explicit skeleton_element(const skeleton_element& el)
    : x(el.x), y(el.y), value(el.value),
      is_erodable(el.is_erodable), is_removable(el.is_removable) { }

   /*
    * Destructor.
    */
   ~skeleton_element() { /* do nothing */ }

   /*
    * Data.
    */
   unsigned long x;              /* x-coordinate */
   unsigned long y;              /* y-coordinate */
   double        value;          /* element value */
   bool          is_erodable;    /* is removal an erosion operation? */
   bool          is_removable;   /* will removal preserve topology? */
};

/*
 * Priority comparison functor for skeleton elements.
 */
class skeleton_element_priority_functor
 : public comparable_functor<skeleton_element> {
public:
   int operator()(
      const skeleton_element& el0, const skeleton_element& el1) const
   {
      static const compare_functor<bool>& f_cmp_bool =
         compare_functors<bool>::f_compare();
      static const compare_functor<double>& f_cmp_double = 
         compare_functors<double>::f_compare();
      int cmp_quality = -(f_cmp_bool(el0.is_erodable, el1.is_erodable));
      if (cmp_quality == 0)
         cmp_quality = -(f_cmp_bool(el0.is_removable, el1.is_removable));
      return (
         (cmp_quality == 0) ?
            f_cmp_double(el0.value, el1.value) : cmp_quality
      );
   }
};

/*
 * Search comparison functor for skeleton elements.
 */
class skeleton_element_compare_functor
 : public comparable_functor<skeleton_element> {
public:
   int operator()(
      const skeleton_element& el0, const skeleton_element& el1) const
   {
      static const compare_functor<unsigned long>& f_cmp =
         compare_functors<unsigned long>::f_compare();
      int cmp_x = f_cmp(el0.x, el1.x);
      return ((cmp_x == 0) ? f_cmp(el0.y, el1.y) : cmp_x);
   }
};

/*
 * Update the status of an element on the skeleton.
 *
 * Return whether immediate removal of the element is a simple erosion
 * (is_erodable), whether removal preserves structure (is_removable), and
 * whether the element may be a candidate for removal after other neighboring
 * elements have been removed (is_candidate).
 */
void skeleton_status(
   const double* m,        /* matrix */
   unsigned long size_x,   /* matrix size in x-direction */
   unsigned long size_y,   /* matrix size in y-direction */
   unsigned long x,        /* x-coordinate of element */
   unsigned long y,        /* y-coordinate of element */
   unsigned long ind,      /* linear index of element */
   bool& is_candidate,     /* is the element a candidate for removal? */
   bool& is_erodable,      /* is the element currently erodable? */
   bool& is_removable)     /* is the element currently removable? */
{
   /* check that element has not already been removed */
   if (m[ind] == 0) {
      is_candidate = false;
      is_erodable  = false;
      is_removable = false;
      return;
   }
   /* check which neighbors are nonzero */
   bool n  = false, s  = false, e  = false, w  = false;
   bool ne = false, nw = false, se = false, sw = false;
   neighbor_values_compare_2D(
      m, size_x, size_y, x, y, ind, 0.0, n, s, e, w, ne, nw, se, sw
   );
   n  = !n;   s  = !s;   e  = !e;   w  = !w;
   ne = !ne;  nw = !nw;  se = !se;  sw = !sw;
   /* check if element should be a candidate if it is not currently */
   is_candidate =
      is_candidate || is_thick_region_2D(n, s, e, w, ne, nw, se, sw);
   /* candidates are erodable if removal is a simple erosion */
   is_erodable = is_candidate && is_erodable_2D(n, s, e, w, ne, nw, se, sw);
   /* candidates are immediately removable if removal preserves structure */
   is_removable =
      (is_erodable) ||
      (is_candidate && is_removable_2D(n, s, e, w, ne, nw, se, sw));
}
} /* namespace */

/*
 * Skeletonize the given 2D matrix, preserving its topology.
 *
 * Nonzero elements are treated as part of an area to be skeletonized and
 * zero elements are treated as background.
 *
 * Skeletonization proceeds by repeatedly removing (setting to zero) the
 * smallest nonzero element whos deletion best preserves the topology.
 */
matrix<> lib_image::skeletonize_2D(const matrix<>& m) {
   /* initialize comparison functors */
   static const skeleton_element_priority_functor* f_prio =
      new skeleton_element_priority_functor();
   static const skeleton_element_compare_functor* f_search =
      new skeleton_element_compare_functor();
   /* initialize skeleton matrix */
   matrix<> skel(m);
   /* check that matrix is 2D */
   if (skel._dims.size() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix size */
   unsigned long size_x = skel._dims[0];
   unsigned long size_y = skel._dims[1];
   /* initialize queue of skeleton elements */
   auto_collection< skeleton_element, queue_set<skeleton_element> > q(
      new queue_set<skeleton_element>(*f_prio, *f_search)
   );
   unsigned long ind = 0;
   for (unsigned long x = 0; x < size_x; x++) {
      for (unsigned long y = 0; y < size_y; y++) {
         /* compute removability of element */
         bool is_candidate = false;
         bool is_erodable  = false;
         bool is_removable = false;
         skeleton_status(
            skel._data, size_x, size_y, x, y, ind,
            is_candidate, is_erodable, is_removable
         );
         /* enqueue candidate elements */
         if (is_candidate) {
            auto_ptr<skeleton_element> el(
               new skeleton_element(
                  x, y, skel._data[ind], is_erodable, is_removable
               )
            );
            q->enqueue(*el);
            el.release();
         }
         /* increment linear index */
         ind++;
      }
   }
   /* repeatedly remove elements and update skeleton */
   while (!(q->is_empty())) {
      /* dequeue element */
      auto_ptr<skeleton_element> el(&(q->dequeue()));
      /* check if skeletonization is complete */
      if (!(el->is_removable))
         break;
      /* update skeleton */
      unsigned long el_x = el->x;
      unsigned long el_y = el->y;
      ind = el_x * size_y + el_y;
      skel._data[ind] = 0;
      /* compute size of local neighborhood */
      unsigned long start_x = (el_x > 0) ? (el_x-1) : (el_x);
      unsigned long start_y = (el_y > 0) ? (el_y-1) : (el_y);
      unsigned long end_x = ((el_x+1) < size_x) ? (el_x+1) : (el_x);
      unsigned long end_y = ((el_y+1) < size_y) ? (el_y+1) : (el_y);
      unsigned long n_y = end_y - start_y + 1;
      if (el_x > 0) { ind -= size_y; }
      if (el_y > 0) { ind--; }
      /* update surrounding elements in queue */
      for (unsigned long x = start_x; x <= end_x; x++) {
         for (unsigned long y = start_y; y <= end_y; y++) {
            /* check if on skeleton */
            if (skel._data[ind] != 0) {
               /* check if element is in queue */
               skeleton_element el_search(x, y, 0, false, false);
               if (q->contains(el_search)) {
                  /* find element in queue */
                  skeleton_element& el_curr = q->find(el_search);
                  /* update element status */
                  bool is_candidate = true;
                  skeleton_status(
                     skel._data, size_x, size_y, x, y, ind,
                     is_candidate, el_curr.is_erodable, el_curr.is_removable
                  );
                  /* update element in queue */
                  q->update(el_curr);
               }
            }
            /* increment linear index */
            ind++;
         }
         /* reset linear index */
         ind -= n_y;
         ind += size_y;
      }
   }
   return skel;
}

/***************************************************************************
 * Connected components.
 ***************************************************************************/

/*
 * Label connected components in a matrix.
 *
 * For real-valued matrices, adjacent nonzero elements are considered
 * connected.  For integer-valued matrices, adjacent elements must also
 * have the same value to be considered connected.  In both cases, zero
 * elements denote empty space.
 *
 * Each element is assigned an integer label indicating the component to 
 * which it belongs.  The label zero is assigned to elements not part of 
 * any component.  Labels range from 0 to n, where n is the number of 
 * components.
 *
 * Optionally, a mask matrix can be specified to determine to which elements
 * of the local neighborhood the center element may possibly be connected.
 * The default mask is the local neighborhood of size 3^d where d is the
 * dimensionality of the matrix.
 *
 * NOTE: The mask must be nonempty and should be symmetric about its center
 *       (connectivity is ill-defined when using an asymmetric mask).
 */
matrix<unsigned long> lib_image::connected_components(
   const matrix<>& m)
{
   /* create integer-valued matrix (map all nonzero values to 1) */
   matrix<unsigned long> m_nonzero(m._dims);
   for (unsigned long n = 0; n < m._size; n++) {
      if (m._data[n] != 0)
         m_nonzero._data[n] = 1;
   }
   /* compute connected components */
   return lib_image::connected_components(m_nonzero);
}

matrix<unsigned long> lib_image::connected_components(
   const matrix<unsigned long>& m)
{
   /* create default mask for full 3^d neighborhood */
   unsigned long n_dims = m._dims.size();
   array<unsigned long> mask_dims(n_dims, 3);
   matrix<bool> mask(mask_dims, true);
   if (n_dims > 0)
      mask._data[mask._size/2] = false;   /* don't compare center to itself */
   /* label connected components */
   return lib_image::connected_components(m, mask);
}

matrix<unsigned long> lib_image::connected_components(
   const matrix<>& m, const matrix<bool>& mask)
{
   /* create integer-valued matrix (map all nonzero values to 1) */
   matrix<unsigned long> m_nonzero(m._dims);
   for (unsigned long n = 0; n < m._size; n++) {
      if (m._data[n] != 0)
         m_nonzero._data[n] = 1;
   }
   /* compute connected components */
   return lib_image::connected_components(m_nonzero, mask);
}

matrix<unsigned long> lib_image::connected_components(
   const matrix<unsigned long>& m, const matrix<bool>& mask)
{
   /* check that mask is nonempty */
   if (mask._size == 0)
      throw ex_invalid_argument("mask must be nonempty");
   /* get number of dimensions */
   unsigned long n_dims_m    = m._dims.size();
   unsigned long n_dims_mask = mask._dims.size();
   unsigned long n_dims = (n_dims_m > n_dims_mask) ? n_dims_m : n_dims_mask;
   /* extend dimension arrays */
   array<unsigned long> dims_m    = matrix<>::extend_dims(m._dims, n_dims);
   array<unsigned long> dims_mask = matrix<>::extend_dims(mask._dims, n_dims);
   /* compute step sizes along each dimension */
   array<unsigned long> sizes_m    = matrix<>::dims_sizes(dims_m);
   array<unsigned long> sizes_mask = matrix<>::dims_sizes(dims_mask);
   /* compute mask size to left/right of center along each dimension */
   array<unsigned long> mask_size_left(n_dims);
   array<unsigned long> mask_size_right(n_dims);
   for (unsigned long d = 0; d < n_dims; d++) {
      unsigned long s = dims_mask[d] - 1;
      mask_size_left[d]  = s / 2;
      mask_size_right[d] = s - mask_size_left[d];
   }  
   /* allocate arrays to store position in local neighborhood */
   array<unsigned long> neighbor_pos(n_dims);
   array<unsigned long> neighbor_pos_start(n_dims);
   array<unsigned long> neighbor_pos_end(n_dims);
   array<unsigned long> neighborhood_diff(n_dims);
   /* initialize queue of elements to process */
   auto_collection< array<unsigned long>, list< array<unsigned long> > > q(
      new list< array<unsigned long> >()
   );
   /* initialize labels */
   matrix<unsigned long> labels(m._dims);
   unsigned long next_label = 1;
   /* initialize position in matrix */
   array<unsigned long> pos(n_dims);
   /* follow connected components from each element */
   for (unsigned long n = 0; n < m._size; n++) {
      /* check if current position begins component  */
      unsigned long val = m._data[n];
      if ((val != 0) && (labels._data[n] == 0)) {
         /* label it */
         labels._data[n] = next_label;
         /* enqueue current position */
         {
            auto_ptr< array<unsigned long> > pos_copy(
               new array<unsigned long>(pos)
            );
            q->add(*pos_copy);
            pos_copy.release();
         }
         /* follow current connected component */
         while (!(q->is_empty())) {
            /* dequeue position */
            auto_ptr< array<unsigned long> > curr_pos(&(q->remove_head()));
            /* initialize indices and position in neighborhood */
            unsigned long ind_m    = 0;
            unsigned long ind_mask = 0;
            for (unsigned long d = 0; d < n_dims; d++) {
               unsigned long p = (*curr_pos)[d];
               unsigned long s_left  = mask_size_left[d];
               unsigned long s_right = mask_size_right[d];
               /* compute start position along dimension d */
               unsigned long pos_start_m = 0;
               if (p >= s_left) {
                  pos_start_m = p - s_left;
               } else {
                  ind_mask += (s_left - p) * sizes_mask[d];
               }
               ind_m += pos_start_m * sizes_m[d];
               /* compute end position along dimension d*/
               unsigned long pos_end_m = 
                  ((p + s_right) < dims_m[d]) ? (p + s_right) : (dims_m[d] - 1);
               /* store position in neighborhood */
               neighbor_pos[d]       = pos_start_m;
               neighbor_pos_start[d] = pos_start_m;
               neighbor_pos_end[d]   = pos_end_m;
               neighborhood_diff[d]  = pos_end_m - pos_start_m;
            }
            /* label and enqueue any unlabeled adjacent positions */
            while (true) {
               /* check if unlabeled and adjacent */
               if ((m._data[ind_m] == val) &&
                   (labels._data[ind_m] == 0) &&
                   (mask._data[ind_mask]))
               {
                  /* label */
                  labels._data[ind_m] = next_label;
                  /* enqueue */
                  {
                     auto_ptr< array<unsigned long> > pos_copy(
                        new array<unsigned long>(neighbor_pos)
                     );
                     q->add(*pos_copy);
                     pos_copy.release();
                  }
               }
               /* update position in local neighborhood */
               unsigned long d = n_dims - 1;
               while ((neighbor_pos[d] == neighbor_pos_end[d]) && (d > 0)) {
                  ind_m    -= neighborhood_diff[d] * sizes_m[d];
                  ind_mask -= neighborhood_diff[d] * sizes_mask[d];
                  neighbor_pos[d] = neighbor_pos_start[d];
                  d--;
               }
               if (neighbor_pos[d] < neighbor_pos_end[d]) {
                  ind_m    += sizes_m[d];
                  ind_mask += sizes_mask[d];
                  neighbor_pos[d]++;
               } else {
                  break;
               }
            }
         }
         /* increment next label */
         next_label++;
      }
      /* update position */
      {
         unsigned long d = n_dims - 1;
         pos[d]++;
         while ((pos[d] == dims_m[d]) && (d > 0)) {
            pos[d] = 0;
            d--;
            pos[d]++;
         }
      }
   }
   return labels;
}

/*
 * Label connected components along the specified dimension of the matrix.
 * 
 * For real-valued matrices, adjacent nonzero elements are considered
 * connected.  For integer-valued matrices, adjacent elements must also
 * have the same value to be considered connected.  In both cases, zero
 * elements denote empty space.
 *
 * As above, labels range from 0 to n, where n is the number of components.
 */
matrix<unsigned long> lib_image::connected_components(
   const matrix<>& m, unsigned long d)
{
   /* create integer-valued matrix (map all nonzero values to 1) */
   matrix<unsigned long> m_nonzero(m._dims);
   for (unsigned long n = 0; n < m._size; n++) {
      if (m._data[n] != 0)
         m_nonzero._data[n] = 1;
   }
   /* compute connected components along dimension */
   return lib_image::connected_components(m_nonzero, d);
}

matrix<unsigned long> lib_image::connected_components(
   const matrix<unsigned long>& m, unsigned long d)
{
   /* get # steps along dimension */
   unsigned long n_dims  = m._dims.size();
   unsigned long n_steps = (d < n_dims) ? m._dims[d] : 1;
   /* initialize labels */
   matrix<unsigned long> labels(m._dims);
   /* check that result is nontrivial */
   if (labels._size > 0) {
      /* compute indices along first slice of dimension */
      array<unsigned long> inds = matrix<>::linear_indices_slice(m._dims, d);
      /* compute step size along dimension */
      unsigned long step_size = matrix<>::dims_size(m._dims, d);
      /* compute result matrix */
      unsigned long label = 0;
      unsigned long n_inds = inds.size();
      for (unsigned long n = 0; n < n_inds; n++) {
         unsigned long ind = inds[n];
         unsigned long val_prev = 0;
         for (unsigned long n_step = 0; n_step < n_steps; n_step++) {
            unsigned long val = m._data[ind];
            if (val != 0) {
               if (val != val_prev)
                  label++;
               labels._data[ind] = label;
            }
            val_prev = val;
            ind += step_size;
         }
      }
   }
   return labels;
}

/*
 * Return the number of connected components (excluding the zero component)
 * in a labeled component matrix generated by one of the above functions.
 *
 * Since components are numbered in order, this is the same as the maximum
 * value in the matrix.
 */
unsigned long lib_image::count_connected_components(
   const matrix<unsigned long>& labels)
{
   return max(labels);
}

/***************************************************************************
 * Contour and region scanning (2D).
 ***************************************************************************/

namespace {
/*
 * Compute and return the indices of the cells in a grid that intersect the
 * specified vertical line segment.  
 *
 * The entire line segment must lie within the grid.
 */
auto_ptr< array<unsigned long> > scan_vertical_line_segment(
   unsigned long x,              /* segment x-coordinate */
   unsigned long y_start,        /* segment starting y-coordinate */
   unsigned long y_end,          /* segment ending y-coordinate */
   unsigned long size_y)         /* grid size in y-direction */
{
   bool direction = (y_end >= y_start);
   unsigned long n_steps = 
      direction ? (y_end - y_start) : (y_start - y_end);
   auto_ptr< array<unsigned long> > inds(new array<unsigned long>(n_steps + 1));
   unsigned long ind = x * size_y + y_start;
   (*inds)[0] = ind;
   for (unsigned long n = 1; n <= n_steps; n++) {
      if (direction) { ind++; } else { ind--; }
      (*inds)[n] = ind;
   }
   return inds;
}

/*
 * Compute and return the indices of the cells in a grid that intersect the
 * specified horizontal line segment.  
 *
 * The entire line segment must lie within the grid.
 */
auto_ptr< array<unsigned long> > scan_horizontal_line_segment(
   unsigned long y,              /* segment y-coordinate */
   unsigned long x_start,        /* segment starting x-coordinate */
   unsigned long x_end,          /* segment ending x-coordinate */
   unsigned long size_y)         /* grid size in y-direction */
{
   bool direction = (x_end >= x_start);
   unsigned long n_steps = 
      direction ? (x_end - x_start) : (x_start - x_end);
   auto_ptr< array<unsigned long> > inds(new array<unsigned long>(n_steps + 1));
   unsigned long ind = x_start * size_y + y;
   (*inds)[0] = ind;
   for (unsigned long n = 1; n <= n_steps; n++) {
      if (direction) { ind += size_y; } else { ind -= size_y; }
      (*inds)[n] = ind;
   }
   return inds;
}

/*
 * Compute and return the indices of the cells in a grid that intersect the
 * straight line segment between with the given start and end points.
 *
 * The entire line segment must lie within the grid.
 */
auto_ptr< array<unsigned long> > scan_line_segment(
   const point_2D& p_start,      /* line segment start */
   const point_2D& p_end,        /* line segment end */
   unsigned long size_y,         /* grid size in y-direction */
   bool allow_diagonal = true)   /* allow diagonal connectivity? */
{
   /* get segment endpoint coordinates */
   double x_start = p_start.x();
   double y_start = p_start.y();
   double x_end = p_end.x();
   double y_end = p_end.y();
   /* compute grid coordinates */ 
   unsigned long x_start_grid = static_cast<unsigned long>(math::round(x_start));
   unsigned long y_start_grid = static_cast<unsigned long>(math::round(y_start));
   unsigned long x_end_grid = static_cast<unsigned long>(math::round(x_end));
   unsigned long y_end_grid = static_cast<unsigned long>(math::round(y_end));
   /* check slope */
   if (x_start_grid == x_end_grid) {
      /* vertical line segment */
      return scan_vertical_line_segment(
         x_start_grid, y_start_grid, y_end_grid, size_y
      );
   } else if (y_start_grid == y_end_grid) {
      /* horizontal line segment */
      return scan_horizontal_line_segment(
         y_start_grid, x_start_grid, x_end_grid, size_y
      );
   } else {
      /* compute direction of motion */
      bool direction_x = (x_end_grid >= x_start_grid);
      bool direction_y = (y_end_grid >= y_start_grid);
      /* compute offsets of next pixel boundaries */
      double delta_x = direction_x ? 0.5 : -0.5;
      double delta_y = direction_y ? 0.5 : -0.5;
      /* compute equation of line */
      double slope = (y_end - y_start) / (x_end - x_start);
      double offset = y_start - slope * x_start;
      /* determine how to break ties when moving along diagonal */
      bool prefer_vertical = (math::abs(slope) >= 1);
      /* initialize index array */
      unsigned long max_steps_x = direction_x ?
         (x_end_grid - x_start_grid) : (x_start_grid - x_end_grid);
      unsigned long max_steps_y = direction_y ?
         (y_end_grid - y_start_grid) : (y_start_grid - y_end_grid);
      unsigned long max_steps = max_steps_x + max_steps_y;
      auto_ptr< array<unsigned long> > inds(
         new array<unsigned long>(max_steps + 1)
      );
      /* initialize indices */
      unsigned long ind = x_start_grid * size_y + y_start_grid;
      unsigned long ind_end = x_end_grid * size_y + y_end_grid;
      /* store starting point */
      (*inds)[0] = ind;
      /* initialize current coordinate */
      double x = static_cast<double>(x_start_grid);
      double y = static_cast<double>(y_start_grid);
      /* initialize segment length */
      unsigned long n_points = 1;
      /* move along line segment */
      while ((ind != ind_end) && (n_points <= max_steps)) {
         /* compute location of next grid cell intersected by segment */
         double y_boundary = y + delta_y;
         double y_at_x_boundary = slope * (x + delta_x) + offset;
         bool move_vertical   = false;
         bool move_horizontal = false;
         if (y_boundary < y_at_x_boundary) {
            move_vertical   = direction_y;
            move_horizontal = !direction_y;
         } else if (y_boundary > y_at_x_boundary) {
            move_vertical   = !direction_y;
            move_horizontal = direction_y;
         } else {
            move_vertical   = allow_diagonal || prefer_vertical;
            move_horizontal = allow_diagonal || !prefer_vertical;
         }
         /* move to next grid cell */
         if (move_vertical) {
            if (direction_y) {
               y++;
               ind++;
            } else {
               y--;
               ind--;
            }
         }
         if (move_horizontal) {
            if (direction_x) {
               x++;
               ind += size_y;
            } else {
               x--;
               ind -= size_y;
            }
         }
         /* store index */
         (*inds)[n_points++] = ind;
      }
      inds->resize(n_points);
      return inds;
   }
}

/*
 * Check whether the specified point lies within the grid of the given size.
 */
bool is_in_grid(
   const point_2D& p,         /* point */
   unsigned long size_x,      /* grid size in x-direction */
   unsigned long size_y)      /* grid size in y-direction */
{
   return ((0 <= p.x()) && (p.x() <= (static_cast<double>(size_x) - 1.0)) &&
           (0 <= p.y()) && (p.y() <= (static_cast<double>(size_y) - 1.0)));
}

/*
 * Clip the line segment so that only the portion lying completely within the
 * grid of the specified size remains.  Update the line segment endpoint.
 *
 * The start point must already lie within the grid.  The end point will be
 * clipped.
 */
void clip_line_segment_to_grid(
   const point_2D& p_start,   /* line segment start (must be within grid) */
   point_2D& p_end,           /* line segment end (will be clipped to grid) */
   unsigned long size_x,      /* grid size in x-direction */
   unsigned long size_y)      /* grid size in y-direction */
{
   /* get segment endpoint coordinates */
   double x_start = p_start.x();
   double y_start = p_start.y();
   double x_end = p_end.x();
   double y_end = p_end.y();
   /* compute bound in each dimension */
   double x_bound = static_cast<double>(size_x) - 1.0;
   double y_bound = static_cast<double>(size_y) - 1.0;
   /* adjust so that x-coordinate is within range */
   if (x_end < 0) {
      y_end = y_start + (y_end-y_start)/(x_end-x_start) * (-x_start);
      x_end = 0;
   } else if (x_end > x_bound) {
      y_end = y_start + (y_end-y_start)/(x_end-x_start) * (x_bound-x_start);
      x_end = x_bound;
   }
   /* adjust so that y-coordinate is within range */
   if (y_end < 0) {
      x_end = x_start + (x_end-x_start)/(y_end-y_start) * (-y_start);
      y_end = 0;
   } else if (y_end > y_bound) {
      x_end = x_start + (x_end-x_start)/(y_end-y_start) * (y_bound-y_start);
      y_end = y_bound;
   }
   /* update endpoint */
   p_end.x() = x_end;
   p_end.y() = y_end;
}

/*
 * Compute and return the indices of the cells in a grid that intersect the
 * circular disc of radius r around the point p.  These indices are computed 
 * by shooting a ray from p to each point with integer coordinates on the 
 * boundary of the disc.  We return a collection of index sets, one for each
 * ray.  The complete set of indices is the union of indices in the returned
 * arrays.
 *
 * Note that p must lie within the grid.
 */
auto_collection< array<unsigned long> > scan_disc(
   const point_2D& p,         /* disc center (must be within grid) */
   double r,                  /* disc radius */
   unsigned long size_x,      /* grid size in x-direction */
   unsigned long size_y)      /* grid size in y-direction */
{
   /* initialize collection of indices */
   auto_collection< array<unsigned long> > inds(
      new list< array<unsigned long> >()
   );
   /* loop over boundary */
   r = math::abs(r);
   unsigned long n_steps = static_cast<unsigned long>(math::ceil(2*M_PIl*r))+1;
   double step_size = 2*M_PIl / static_cast<double>(n_steps);
   double theta = 0;
   for (unsigned long n = 0; n < n_steps; n++, theta += step_size) {
      /* compute point on disc boundary */
      double x_offset = r * math::cos(theta);
      double y_offset = r * math::sin(theta);
      point_2D q(p.x() + x_offset, p.y() + y_offset);
      /* scan line segment from center to boundary */
      clip_line_segment_to_grid(p, q, size_x, size_y);
      auto_ptr< array<unsigned long> > inds_segment = 
         scan_line_segment(p, q, size_y, false);
      /* store indices */
      inds->add(*inds_segment);
      inds_segment.release();
   }
   return inds;
}
} /* namespace */

/***************************************************************************
 * Contour extraction (2D).
 ***************************************************************************/

namespace {
/*
 * Create a contour vertex with the specified coordinates.
 */
auto_ptr<lib_image::contour_vertex> create_contour_vertex(
   unsigned long x, unsigned long y)
{
   auto_ptr<lib_image::contour_vertex> v(new lib_image::contour_vertex());
   /* set vertex identity to defaults */
   v->id = 0;
   v->is_subdivision = false;
   /* set vertex coordinates */
   v->x = x;
   v->y = y;
   return v;
}

/*
 * Create a contour edge between the two specified vertices.
 */
auto_ptr<lib_image::contour_edge> create_contour_edge(
   lib_image::contour_vertex& v_start, lib_image::contour_vertex& v_end)
{
   auto_ptr<lib_image::contour_edge> e(new lib_image::contour_edge());
   /* set edge identity to defaults */
   e->id = 0;
   e->contour_equiv_id = 0;
   e->is_completion = false;
   /* set edge vertices */
   e->vertex_start = &v_start;
   e->vertex_end   = &v_end;
   /* update vertex edge lists */
   e->vertex_start_enum = v_start.edges_start.size();
   e->vertex_end_enum   = v_end.edges_end.size();
   v_start.edges_start.add(*e);
   v_end.edges_end.add(*e);
   return e;
}
} /* namespace */

/*
 * Return the point representing the contour vertex.
 */
point_2D lib_image::contour_vertex::point() const {
   return point_2D(
      static_cast<double>(this->x),
      static_cast<double>(this->y)
   );
}

/*
 * Return the point at the given position along the contour edge.
 */
point_2D lib_image::contour_edge::point(unsigned long n) const {
   return point_2D(
      static_cast<double>(this->x_coords[n]),
      static_cast<double>(this->y_coords[n])
   );
}
   
/*
 * Return the size (number of interior points) of the contour edge.
 */
unsigned long lib_image::contour_edge::size() const {
   return this->x_coords.size();
}

/*
 * Return the length (distance between endpoints) of the contour edge.
 */
double lib_image::contour_edge::length() const {
   return abs((this->vertex_end->point()) - (this->vertex_start->point()));
}

/*
 * Extract contours from labeled connected components.
 *
 * Break each connected component into a set of thin contours that may
 * intersect only at returned contour vertices.
 */
lib_image::contour_set::contour_set(const matrix<unsigned long>& labels)
 : _is_vertex(),
   _is_edge(),
   _assignments(),
   _vertices(new array_list<contour_vertex>()),
   _edges(new array_list<contour_edge>()),
   _edges_equiv(new array_list< list<contour_edge> >())
{
   /* check that label matrix is 2D */
   if (labels.dimensionality() != 2)
      throw ex_invalid_argument("matrix must be 2D");
   /* get matrix size */
   unsigned long size_x = labels.size(0);
   unsigned long size_y = labels.size(1);
   /* initialize pixel assignment map */
   _is_vertex.resize(size_x, size_y);
   _is_edge.resize(size_x, size_y);
   _assignments.resize(size_x, size_y, 0 /* init value */);
   /* initialize map of connected component label -> vertex in it */
   unsigned long n_labels = lib_image::count_connected_components(labels);
   array<bool>          has_vertex(n_labels);
   array<unsigned long> label_x(n_labels);
   array<unsigned long> label_y(n_labels);
   /* find all vertices */
   for (unsigned long ind = 0, x = 0; x < size_x; x++) {
      for (unsigned long y = 0; y < size_y; y++) {
         /* check if part of an edge */
         unsigned long label = labels[ind];
         if (label != 0) {
            /* check if neighbors have the same label */
            bool n  = false, s  = false, e  = false, w  = false;
            bool ne = false, nw = false, se = false, sw = false;
            neighbor_values_compare_2D(
               labels.data(), size_x, size_y, x, y, ind, label,
               n, s, e, w, ne, nw, se, sw
            );
            /* check vertex conditions */
            if (is_vertex_2D(n, s, e, w, ne, nw, se, sw)) {
               /* create vertex and set its id */
               auto_ptr<contour_vertex> v = create_contour_vertex(x, y);
               v->id = _vertices->size();
               /* record vertex assignment */
               _assignments[ind] = v->id;
               /* add vertex to vertex set */
               _vertices->add(*v);
               v.release();
               /* mark vertex */
               _is_vertex[ind] = true;
               /* indicate vertex found for current label */
               has_vertex[label-1] = true;
            }
            /* record position of element with current label */
            label_x[label-1] = x;
            label_y[label-1] = y;
         }
         /* increment linear coordinate */
         ind++;
      }
   }
   /* detect loops and arbitrarily pick a vertex for any loop */
   for (unsigned long n = 0; n < n_labels; n++) {
      if (!has_vertex[n]) {
         /* create vertex and set its id */
         auto_ptr<contour_vertex> v = create_contour_vertex(
            label_x[n], label_y[n]
         );
         v->id = _vertices->size();
         /* record vertex assignment */
         _assignments(label_x[n], label_y[n]) = v->id;
         /* add vertex to vertex set */
         _vertices->add(*v);
         v.release();
         /* mark vertex */
         _is_vertex(label_x[n], label_y[n]) = true;
         /* indicate vertex found for current label */
         has_vertex[n] = true;
      }
   }
   /* process all vertices */
   for (unsigned long v_id = 0; v_id < _vertices->size(); v_id++) {
      /* grab vertex */
      contour_vertex& v = (*_vertices)[v_id];
      unsigned long label = labels(v.x, v.y);
      /* initialize queue of vertex neighbors */
      array<unsigned long> q_x(8 /* max # of neighbors */);
      array<unsigned long> q_y(8 /* max # of neighbors */);
      unsigned long n_neighbors = 0;
      /* enqueue coordinates of adjacent neighbors with the same label */
      unsigned long x_start = (v.x > 0) ? (v.x - 1) : (v.x + 1);
      unsigned long x_end   = (v.x + 1);
      for (unsigned long x = x_start; (x <= x_end) && (x < size_x); x += 2) {
         if (labels(x, v.y) == label) {
            q_x[n_neighbors] = x;
            q_y[n_neighbors] = v.y;
            n_neighbors++;
         }
      }
      unsigned long y_start = (v.y > 0) ? (v.y - 1) : (v.y + 1);
      unsigned long y_end   = (v.y + 1);
      for (unsigned long y = y_start; (y <= y_end) && (y < size_y); y += 2) {
         if (labels(v.x, y) == label) {
            q_x[n_neighbors] = v.x;
            q_y[n_neighbors] = y;
            n_neighbors++;
         }
      }
      /* enqueue coordinates of diagonal neighbors with the same label */
      /* only consider diagonal neighbors unreachable through adjacent ones */
      for (unsigned long x = x_start; (x <= x_end) && (x < size_x); x += 2) {
         for (unsigned long y = y_start; (y <= y_end) && (y < size_y); y += 2) {
            if ((labels(x, y) == label) &&
                (labels(v.x, y) != label) &&
                (labels(x, v.y) != label))
            {
               q_x[n_neighbors] = x;
               q_y[n_neighbors] = y;
               n_neighbors++;
            }
         }
      }
      /* process neighbors */
      for (unsigned long n = 0; n < n_neighbors; n++) {
         /* get neighbor index */
         unsigned long ind = q_x[n] * size_y + q_y[n];
         /* check neighbor status */
         if (_is_vertex[ind]) {
            /* check if vertex with higher id -> make edge */
            if (_assignments[ind] > v_id) {
               /* create edge and set its identity */
               contour_vertex& v_end = (*_vertices)[_assignments[ind]];
               auto_ptr<contour_edge> e = create_contour_edge(v, v_end);
               e->id               = _edges->size();
               e->contour_equiv_id = e->id;
               /* create edge equivalence class */
               auto_ptr< list<contour_edge> > e_equiv(new list<contour_edge>());
               e_equiv->add(*e);
               /* add edge */
               _edges->add(*e);             e.release();
               _edges_equiv->add(*e_equiv); e_equiv.release();
            }
         } else if (!_is_edge[ind]) {
            /* neighbor not assigned to an edge -> trace edge */
            this->trace_edge(labels, v_id, q_x[n], q_y[n]);
         }
      }
   }
   /* find and remove unused vertices (single pixel islands) */
   unsigned long n_vertices = _vertices->size();
   for (unsigned long v_id = 0, n = 0; n < n_vertices; n++) {
      auto_ptr<contour_vertex> v(&(_vertices->remove_head()));
      if ((v->edges_start.is_empty()) && (v->edges_end.is_empty())) {
         /* remove vertex */
         _is_vertex(v->x, v->y) = false;
      } else {
         /* update vertex id */
         v->id = v_id++;
         _assignments(v->x, v->y) = v->id;
         /* add vertex */
         _vertices->add(*v);
         v.release();
      }
   }
}

/*
 * Copy constructor.
 */
lib_image::contour_set::contour_set(const contour_set& contours)
 : _is_vertex(contours._is_vertex),
   _is_edge(contours._is_edge),
   _assignments(contours._assignments),
   _vertices(new array_list<contour_vertex>()),
   _edges(new array_list<contour_edge>()),
   _edges_equiv(new array_list< list<contour_edge> >())
{
   /* copy vertices */
   for (unsigned long n = 0, n_v = contours._vertices->size(); n < n_v; n++) {
      /* create vertex copy */
      const contour_vertex& v = (*(contours._vertices))[n];
      auto_ptr<contour_vertex> v_copy = create_contour_vertex(v.x, v.y);
      v_copy->id             = v.id;
      v_copy->is_subdivision = v.is_subdivision;
      /* store vertex copy */
      _vertices->add(*v_copy);
      v_copy.release();
   }
   /* copy edges */
   for (unsigned long n = 0, n_e = contours._edges->size(); n < n_e; n++) {
      /* create edge copy */
      const contour_edge& e = (*(contours._edges))[n];
      auto_ptr<contour_edge> e_copy = create_contour_edge(
         (*_vertices)[e.vertex_start->id], (*_vertices)[e.vertex_end->id]
      );
      e_copy->id               = e.id;
      e_copy->contour_equiv_id = e.contour_equiv_id;
      e_copy->is_completion    = e.is_completion;
      e_copy->x_coords         = e.x_coords;
      e_copy->y_coords         = e.y_coords;
      /* store edge copy */
      _edges->add(*e_copy);
      e_copy.release();
   }
   /* copy edge equivalence classes */
   for (unsigned long n = 0, n_e = contours._edges_equiv->size(); n < n_e; n++)
   {
      /* create equivalence class copy */
      const list<contour_edge>& edge_equiv = (*(contours._edges_equiv))[n];
      auto_ptr< list<contour_edge> > edge_equiv_copy(new list<contour_edge>());
      for (list<contour_edge>::iterator_t i(edge_equiv); i.has_next(); ) {
         const contour_edge& e = i.next();
         edge_equiv_copy->add((*_edges)[e.id]);
      }
      /* store equivalence class copy */
      _edges_equiv->add(*edge_equiv_copy);
      edge_equiv_copy.release();
   }
}

/*
 * Destructor.
 */
lib_image::contour_set::~contour_set() {
   /* do nothing */
}

/***************************************************************************
 * Contour edge operations.
 ***************************************************************************/

/*
 * Given connected component labels, the starting vertex of an edge, 
 * and the next pixel along the edge, trace the entire contour and 
 * add it to the set.
 */
void lib_image::contour_set::trace_edge(
   const matrix<unsigned long>& labels,
   unsigned long                v_id,
   unsigned long                x,
   unsigned long                y)
{
   /* get matrix size */
   unsigned long size_x = labels.size(0);
   unsigned long size_y = labels.size(1);
   /* grab current vertex and label */
   contour_vertex& v = (*_vertices)[v_id];
   unsigned long ind_v = (v.x) * size_y + (v.y);
   unsigned long label = labels[ind_v];
   /* initialize coordinates along edge */
   auto_collection< unsigned long, array_list<unsigned long> > e_x(
      new array_list<unsigned long>()
   );
   auto_collection< unsigned long, array_list<unsigned long> > e_y(
      new array_list<unsigned long>()
   );
   /* temporarily mark start vertex as inaccessible */
   _is_edge[ind_v] = true;
   /* trace edge */
   unsigned long e_id = _edges->size();
   bool endpoint_is_new_vertex = false;
   do {
      /* add latest coordinate to edge */
      auto_ptr<unsigned long> ptr_x(new unsigned long(x));
      auto_ptr<unsigned long> ptr_y(new unsigned long(y));
      e_x->add(*ptr_x); ptr_x.release();
      e_y->add(*ptr_y); ptr_y.release();
      /* mark as edge */
      unsigned long ind = x * size_y + y;
      _is_edge[ind]     = true;
      _assignments[ind] = e_id;
      /* compute which neighbors are present in matrix */
      bool has_n  = false, has_s  = false, has_e  = false, has_w  = false;
      bool has_ne = false, has_nw = false, has_se = false, has_sw = false;
      neighbor_exists_2D(
         size_x, size_y, x, y,
         has_n, has_s, has_e, has_w, has_ne, has_nw, has_se, has_sw
      );
      /* compute indices of neighbors */
      unsigned long ind_n  = 0, ind_s  = 0, ind_e  = 0, ind_w  = 0;
      unsigned long ind_ne = 0, ind_nw = 0, ind_se = 0, ind_sw = 0;
      neighbor_indices_2D(
         size_y, ind,
         ind_n, ind_s, ind_e, ind_w, ind_ne, ind_nw, ind_se, ind_sw
      );
      /* check which neighbors are available */
      bool av_n = has_n ?
         ((labels[ind_n] == label) && (!_is_edge[ind_n])) : false;
      bool av_s = has_s ?
         ((labels[ind_s] == label) && (!_is_edge[ind_s])) : false;
      bool av_e = has_e ?
         ((labels[ind_e] == label) && (!_is_edge[ind_e])) : false;
      bool av_w = has_w ?
         ((labels[ind_w] == label) && (!_is_edge[ind_w])) : false;
      bool av_ne = has_ne ?
         ((labels[ind_ne] == label) && (!_is_edge[ind_ne])) : false;
      bool av_nw = has_nw ?
         ((labels[ind_nw] == label) && (!_is_edge[ind_nw])) : false;
      bool av_se = has_se ?
         ((labels[ind_se] == label) && (!_is_edge[ind_se])) : false;
      bool av_sw = has_sw ?
         ((labels[ind_sw] == label) && (!_is_edge[ind_sw])) : false;
      /* check for adjacent non-vertex */
      if (av_n) { if (!_is_vertex[ind_n]) { y++; continue; } }
      if (av_s) { if (!_is_vertex[ind_s]) { y--; continue; } }
      if (av_e) { if (!_is_vertex[ind_e]) { x++; continue; } }
      if (av_w) { if (!_is_vertex[ind_w]) { x--; continue; } }
      /* check for adjacent vertex */
      if (av_n) { y++; break; }
      if (av_s) { y--; break; }
      if (av_e) { x++; break; }
      if (av_w) { x--; break; }
      /* check for diagonal non-vertex */
      if (av_ne) { if (!_is_vertex[ind_ne]) { x++; y++; continue; } }
      if (av_nw) { if (!_is_vertex[ind_nw]) { x--; y++; continue; } }
      if (av_se) { if (!_is_vertex[ind_se]) { x++; y--; continue; } }
      if (av_sw) { if (!_is_vertex[ind_sw]) { x--; y--; continue; } }
      /* check for diagonal vertex */
      if (av_ne) { x++; y++; break; }
      if (av_nw) { x--; y++; break; }
      if (av_se) { x++; y--; break; }
      if (av_sw) { x--; y--; break; }
      /* otherwise, we have discovered a new vertex */
      endpoint_is_new_vertex = true;
      _is_edge[ind] = false;  /* last point is a new vertex, not an edge */
   } while (!endpoint_is_new_vertex);
   /* remove inaccessible mark from start vertex */
   _is_edge[ind_v] = false;
   /* add endpoint as vertex (if needed) */
   unsigned long ind_end = x * size_y + y;
   if (endpoint_is_new_vertex) {
      /* create vertex and set its id */
      auto_ptr<contour_vertex> v_end = create_contour_vertex(x, y);
      v_end->id = _vertices->size();
      /* record vertex assignment */
      _assignments[ind_end] = v_end->id;
      /* add vertex to vertex set */
      _vertices->add(*v_end);
      v_end.release();
      /* mark vertex */
      _is_vertex[ind_end] = true;
   }
   /* create edge and set its identity */
   contour_vertex& v_e = (*_vertices)[_assignments[ind_end]];
   auto_ptr<contour_edge> e = create_contour_edge(v, v_e);
   e->id               = e_id;
   e->contour_equiv_id = e_id;
   /* store edge coordinates */
   unsigned long n_edge_points = 
      (endpoint_is_new_vertex ? (e_x->size() - 1) : (e_x->size()));
   e->x_coords.resize(n_edge_points);
   e->y_coords.resize(n_edge_points);
   for (unsigned long n = 0; n < n_edge_points; n++) {
      e->x_coords[n] = (*e_x)[n];
      e->y_coords[n] = (*e_y)[n];
   }
   /* create edge equivalence class */
   auto_ptr< list<contour_edge> > e_equiv(new list<contour_edge>());
   e_equiv->add(*e);
   /* add edge */
   _edges->add(*e);             e.release();
   _edges_equiv->add(*e_equiv); e_equiv.release();
}

/*
 * Split the contour edge with the given id at the specified position.
 */
void lib_image::contour_set::split_edge(unsigned long e_id, unsigned long p) {
   /* get the edge to split */
   contour_edge& e = (*_edges)[e_id];
   /* check that split position is valid */
   unsigned long n_points = e.size();
   if (p >= n_points)
      throw ex_index_out_of_bounds("invalid contour split point", p);
   /* create split vertex */
   auto_ptr<contour_vertex> v_split = create_contour_vertex(
      e.x_coords[p], e.y_coords[p]
   );
   v_split->id             = _vertices->size();
   v_split->is_subdivision = true;
   /* record vertex assignment */
   unsigned long ind_v = (v_split->x) * (this->size_y()) + (v_split->y);
   _assignments[ind_v] = v_split->id;
   /* mark vertex */
   _is_edge[ind_v] = false;
   _is_vertex[ind_v] = true;
   /* create split edge */
   auto_ptr<contour_edge> e_split(new contour_edge());
   e_split->id               = _edges->size();
   e_split->contour_equiv_id = e.contour_equiv_id;
   e_split->is_completion    = e.is_completion;
   /* update start and end vertices */
   e_split->vertex_start = v_split.get();
   e_split->vertex_end   = e.vertex_end;
   e.vertex_end = e_split->vertex_start;
   /* update existing edge linkage */
   e_split->vertex_end->edges_end.replace(e.vertex_end_enum, *e_split);
   e_split->vertex_end_enum = e.vertex_end_enum;
   /* add edge linkage to split vertex */
   v_split->edges_start.add(*e_split);
   v_split->edges_end.add(e);
   e_split->vertex_start_enum = 0;
   e.vertex_end_enum          = 0;
   /* split points along edge */
   unsigned long n_points_left  = p;
   unsigned long n_points_right = (n_points - 1) - p;
   if (n_points_right > 0) {
      e_split->x_coords = e.x_coords.subarray((p + 1), (n_points - 1));
      e_split->y_coords = e.y_coords.subarray((p + 1), (n_points - 1));
   }
   e.x_coords.resize(n_points_left);
   e.y_coords.resize(n_points_left);
   /* update edge assignments */
   for (unsigned long n = 0; n < n_points_right; n++)
      _assignments(e_split->x_coords[n], e_split->y_coords[n]) = e_split->id;
   /* update vertex set */
   _vertices->add(*v_split);
   v_split.release();
   /* update edge set */
   (*_edges_equiv)[e_split->contour_equiv_id].add(*e_split);
   _edges->add(*e_split);
   e_split.release();
}

/*
 * Sort the list of edges in each edge equivalence class so they appear
 * in the order encountered on a traversal of the original contour edge.
 */
void lib_image::contour_set::sort_edges_equiv() {
   /* initialize map of subdivision vertex id -> next edge */
   unsigned long n_vertices = _vertices->size();
   unsigned long n_edges    = _edges->size();
   array<unsigned long> next_edge(n_vertices);
   array<bool>          has_prev(n_vertices);
   /* process each equivalence class */
   unsigned long n_equiv = _edges_equiv->size();
   for (unsigned long equiv_id = 0; equiv_id < n_equiv; equiv_id++) {
      /* get set of edges in equivalence class */
      list<contour_edge>& equiv_set = (*_edges_equiv)[equiv_id];
      /* build vertex map */
      for (list<contour_edge>::iterator_t i(equiv_set); i.has_next(); ) {
         const contour_edge& e = i.next();
         next_edge[e.vertex_start->id] = e.id;
         has_prev[e.vertex_end->id] = true;
      }
      /* find the id of the starting edge */
      unsigned long set_size = equiv_set.size();
      unsigned long e_id = (set_size == 0) ? n_edges : equiv_set.head().id;
      for (unsigned long n = 0; n < set_size; n++) {
         const contour_edge& e = equiv_set.remove_head();
         if (!(has_prev[e.vertex_start->id]))
            e_id = e.id;
      }
      /* reconstruct equivalence class in sorted order */
      for (unsigned long n = 0; n < set_size; n++) {
         contour_edge& e = (*_edges)[e_id];
         equiv_set.add(e);
         e_id = next_edge[e.vertex_end->id];
      }
      /* reset vertex previous edge indicators */
      for (list<contour_edge>::iterator_t i(equiv_set); i.has_next(); ) {
         const contour_edge& e = i.next();
         has_prev[e.vertex_end->id] = false;
      }
   }
}

/***************************************************************************
 * Contour vertex and edge access.
 ***************************************************************************/

/*
 * Get image size in the x-dimension.
 */
unsigned long lib_image::contour_set::size_x() const {
   return _assignments.size(0);
}

/*
 * Get image size in the y-dimension.
 */
unsigned long lib_image::contour_set::size_y() const {
   return _assignments.size(1);
}

/*
 * Get number of vertices in the contour set.
 */
unsigned long lib_image::contour_set::vertices_size() const {
   return _vertices->size();
}

/*
 * Get number of edges in the contour set.
 */
unsigned long lib_image::contour_set::edges_size() const {
   return _edges->size();
}

/*
 * Get the number of edge equivalence classes.
 * Each equivalence class represents a contour before subdivision.
 */
unsigned long lib_image::contour_set::edges_equiv_size() const {
   return _edges_equiv->size();
}

/*
 * Return the specified vertex (by reference).
 */
const lib_image::contour_vertex& lib_image::contour_set::vertex(
   unsigned long v_id) const
{
   if (v_id >= _vertices->size())
      throw ex_index_out_of_bounds("invalid contour vertex id", v_id);
   return (*_vertices)[v_id];
}

/*
 * Return the specified edge (by reference).
 */
const lib_image::contour_edge& lib_image::contour_set::edge(
   unsigned long e_id) const
{
   if (e_id >= _edges->size())
      throw ex_index_out_of_bounds("invalid contour edge id", e_id);
   return (*_edges)[e_id];
}

/*
 * Return the specified edge equivalence class (by reference).
 *
 * The edges are listed in the order encountered on a traversal of the 
 * original contour before subdivision.
 */
const list<lib_image::contour_edge>& lib_image::contour_set::edge_equiv(
   unsigned long equiv_id) const
{
   if (equiv_id >= _edges_equiv->size())
      throw ex_index_out_of_bounds("invalid contour equivalence id", equiv_id);
   return (*_edges_equiv)[equiv_id];
}

/*
 * Check if the specified pixel is a vertex.
 */
bool lib_image::contour_set::is_vertex(unsigned long n) const {
   return _is_vertex[n];
}

/*
 * Check if the specified pixel is a vertex.
 */
bool lib_image::contour_set::is_vertex(unsigned long x, unsigned long y) const {
   return _is_vertex(x,y);
}

/*
 * Check if the specified pixel lies on a (non-completion) edge interior.
 */
bool lib_image::contour_set::is_edge(unsigned long n) const {
   return _is_edge[n];
}

/*
 * Check if the specified pixel lies on a (non-completion) edge interior.
 */
bool lib_image::contour_set::is_edge(unsigned long x, unsigned long y) const {
   return _is_edge(x,y);
}

/*
 * Return the id of the vertex at the given pixel.
 * Throw an exception (ex_not_found) if no such vertex exists.
 */
unsigned long lib_image::contour_set::vertex_id(unsigned long n) const {
   if (!(_is_vertex[n]))
      throw ex_not_found("vertex not found at specified pixel location");
   return _assignments[n];
}

/*
 * Return the id of the vertex at the given pixel.
 * Throw an exception (ex_not_found) if no such vertex exists.
 */
unsigned long lib_image::contour_set::vertex_id(
   unsigned long x, unsigned long y) const
{
   if (!(_is_vertex(x,y)))
      throw ex_not_found("vertex not found at specified pixel location");
   return _assignments(x,y);
}

/*
 * Return the id of the (non-completion) edge whos interior contains the
 * given pixel.  Throw an exception (ex_not_found) if no such edge exists.
 */
unsigned long lib_image::contour_set::edge_id(unsigned long n) const {
   if (!(_is_edge[n]))
      throw ex_not_found("edge not found at specified pixel location");
   return _assignments[n];
}

/*
 * Return the id of the (non-completion) edge whos interior contains the
 * given pixel.  Throw an exception (ex_not_found) if no such edge exists.
 */
unsigned long lib_image::contour_set::edge_id(
   unsigned long x, unsigned long y) const
{
   if (!(_is_edge(x,y)))
      throw ex_not_found("edge not found at specified pixel location");
   return _assignments(x,y);
}

/***************************************************************************
 * Contour subdivision.
 ***************************************************************************/

/*
 * Recursively subdivide the contours until they are approximately linear
 * using the following local criteria:
 *
 * (1) Break a contour if the angle formed by the ray from one contour
 *     endpoint to the other and the ray from that endpoint to another
 *     point on the contour exceeds the maximum deviation angle.
 * 
 * (2) Break a contour if the distance from any point on it to the line
 *     segment connecting its endpoints both exceeds some fraction of the
 *     length of this segment and is greater than the tolerance distance.
 *
 * A contour is always subdivided at the point that is maximally distant
 * from the line segment between its endpoints, regardless of whether
 * criterion (1) or (2) triggered the subdivision.
 */
void lib_image::contour_set::subdivide_local(
   bool   enforce_angle,
   double max_angle,
   bool   enforce_dist,
   double max_dist_dev,
   double tol_dist)
{
   /* intialize queue of edges to check */
   list<contour_edge> q(*_edges);
   /* repeatedly dequeue and process edges */
   while (!(q.is_empty())) {
      /* dequeue edge */
      contour_edge& e = q.remove_head();
      /* get endpoints of line segment approximation */
      point_2D p_start = e.vertex_start->point();
      point_2D p_end   = e.vertex_end->point();
      /* compute distance threshold */
      point_2D r_start_end = p_end - p_start;
      double len_start_end = abs(r_start_end);
      double max_dist = max_dist_dev * len_start_end;
      /* initialize split position and distance */
      unsigned long split_pos = 0;
      double split_dist = 0;
      /* determine whether edge should be split */
      bool split = false;
      unsigned long n_points = e.size();
      for (unsigned long n = 0; n < n_points; n++) {
         /* get current point along contour */
         point_2D p = e.point(n);
         /* compute distance from point to segment */
         double dist = p.distance_to_segment(p_start, p_end);
         /* update optimal split position */
         if (dist > split_dist) {
            split_pos  = n;
            split_dist = dist;
         }
         /* check distance constraint */
         if ((!split) && (enforce_dist))
            split = ((dist > max_dist) && (dist > tol_dist));
         /* check angle constraint */
         if ((!split) && (enforce_angle)) {
            /* compute angle between point and segment */
            point_2D r_start_p = p - p_start;
            point_2D r_end_p   = p - p_end;
            double angle_start = math::acos(
               dot(r_start_p, r_start_end) / (abs(r_start_p) * len_start_end)
            );
            double angle_end = math::acos(
               dot(r_end_p, -r_start_end) / (abs(r_end_p) * len_start_end)
            );
            double angle = (angle_start > angle_end) ? angle_start : angle_end;
            /* update split flag */
            split = (angle > max_angle);
         }
      }
      /* split edge if it violates constraints */
      if (split) {
         /* split edge */
         this->split_edge(e.id, split_pos);
         /* enqueue split edges */
         q.add(e);
         q.add(_edges->tail());
      }
   }
   /* order edges in equivalence classes by traversal of original contour */
   this->sort_edges_equiv();
}

namespace {
/*
 * Find and return the position of the point on the contour closest to the 
 * line passing through the specified points.
 *
 * Throw an exception if the contour has no interior points.
 */
unsigned long find_contour_point_closest_to_line(
   const lib_image::contour_edge& e, const point_2D& p, const point_2D& q)
{
   /* check that contour is nonempty */
   unsigned long n_points = e.size();
   if (n_points == 0)
      throw ex_invalid_argument("contour has no interior points");
   /* initialize position and minimum distance */
   unsigned long pos = 0;
   double min_dist = e.point(pos).distance_to_line(p, q);
   /* find closest point */
   for (unsigned long n = 1; n < n_points; n++) {
      double dist = e.point(n).distance_to_line(p, q);
      if (dist < min_dist) {
         pos = n;
         min_dist = dist;
      }
   }
   return pos;
}

/*
 * Find and return the position of the point on the contour closest to the 
 * line segment passing through the specified points.
 *
 * Throw an exception if the contour has no interior points.
 */
unsigned long find_contour_point_closest_to_segment(
   const lib_image::contour_edge& e, const point_2D& p, const point_2D& q)
{
   /* check that contour is nonempty */
   unsigned long n_points = e.size();
   if (n_points == 0)
      throw ex_invalid_argument("contour has no interior points");
   /* initialize position and minimum distance */
   unsigned long pos = 0;
   double min_dist = e.point(pos).distance_to_segment(p, q);
   /* find closest point */
   for (unsigned long n = 1; n < n_points; n++) {
      double dist = e.point(n).distance_to_segment(p, q);
      if (dist < min_dist) {
         pos = n;
         min_dist = dist;
      }
   }
   return pos;
}

/*
 * Find and return the position of the point on the contour farthest from the
 * line passing through the specified points.
 *
 * Throw an exception if the contour has no interior points.
 */
unsigned long find_contour_point_farthest_from_line(
   const lib_image::contour_edge& e, const point_2D& p, const point_2D& q)
{
   /* check that contour is nonempty */
   unsigned long n_points = e.size();
   if (n_points == 0)
      throw ex_invalid_argument("contour has no interior points");
   /* initialize position and maximum distance */
   unsigned long pos = 0;
   double max_dist = 0;
   /* find closest point */
   for (unsigned long n = 0; n < n_points; n++) {
      double dist = e.point(n).distance_to_line(p, q);
      if (dist > max_dist) {
         pos = n;
         max_dist = dist;
      }
   }
   return pos;
}

/*
 * Find and return the position of the point on the contour farthest from the
 * line segment passing through the specified points.
 *
 * Throw an exception if the contour has no interior points.
 */
unsigned long find_contour_point_farthest_from_segment(
   const lib_image::contour_edge& e, const point_2D& p, const point_2D& q)
{
   /* check that contour is nonempty */
   unsigned long n_points = e.size();
   if (n_points == 0)
      throw ex_invalid_argument("contour has no interior points");
   /* initialize position and maximum distance */
   unsigned long pos = 0;
   double max_dist = 0;
   /* find closest point */
   for (unsigned long n = 0; n < n_points; n++) {
      double dist = e.point(n).distance_to_segment(p, q);
      if (dist > max_dist) {
         pos = n;
         max_dist = dist;
      }
   }
   return pos;
}
} /* namespace */

/*
 * Recursively subdivide the contours until their straight line segment
 * approximation respects the true contour topology.
 *
 * Enforce the global criteria that the line segments between contour
 * endpoints intersect only at contour vertices.
 */
void lib_image::contour_set::subdivide_global() {
   /* create comparison functor for sorting edges by vertex id */
   class edge_compare : public comparable_functor<contour_edge> {
   public:
      edge_compare(){} // EDIT: Jordi Pont-Tuset <jordi.pont@upc.edu> to compile on OSX Mavericks
      int operator()(const contour_edge& e0, const contour_edge& e1) const {
         bool v0_cmp = (e0.vertex_start->id < e0.vertex_end->id);
         bool v1_cmp = (e1.vertex_start->id < e1.vertex_end->id);
         unsigned long min0 = v0_cmp ? e0.vertex_start->id : e0.vertex_end->id;
         unsigned long max0 = v0_cmp ? e0.vertex_end->id : e0.vertex_start->id;
         unsigned long min1 = v1_cmp ? e1.vertex_start->id : e1.vertex_end->id;
         unsigned long max1 = v1_cmp ? e1.vertex_end->id : e1.vertex_start->id;
         int cmp_min = (min0 < min1) ? -1 : ((min0 > min1) ? 1 : 0);
         int cmp_max = (max0 < max1) ? -1 : ((max0 > max1) ? 1 : 0);
         return ((cmp_min == 0) ? cmp_max : cmp_min);
      }
   };
   static const edge_compare e_compare;
   /* subdivide edges which share both vertices */
   {
      /* sort edges by vertex ids */
      array_list<contour_edge> edges_sorted(*_edges);
      edges_sorted.sort(e_compare);
      /* find edges sharing both vertices */
      unsigned long n_edges = edges_sorted.size(); 
      unsigned long prev = 0;
      for (unsigned long curr = 1; curr < n_edges; curr++) {
         /* get edges with possibly identical vertices */
         contour_edge& e_prev = edges_sorted[prev];
         contour_edge& e_curr = edges_sorted[curr];
         /* check if edges overlap at both vertices */
         if (e_compare(e_prev, e_curr) == 0) {
            /* get common endpoints */
            point_2D p = e_prev.vertex_start->point();
            point_2D q = e_prev.vertex_end->point();
            /* determine which edge to split */
            bool can_split_prev = (e_prev.size() > 0);
            bool can_split_curr = (e_curr.size() > 0);
            bool split_prev = false;
            if (can_split_prev && can_split_curr) {
               /* compute split position along each edge */
               unsigned long pos_prev = 
                  find_contour_point_farthest_from_segment(e_prev, p, q);
               unsigned long pos_curr = 
                  find_contour_point_farthest_from_segment(e_curr, p, q);
               /* compute distance from split point to line segment */
               double dist_prev =
                  e_prev.point(pos_prev).distance_to_segment(p, q);
               double dist_curr =
                  e_curr.point(pos_curr).distance_to_segment(p, q);
               /* select the least straight edge */
               split_prev = (dist_prev > dist_curr);
            } else if (can_split_prev) {
               /* set previous edge to be split */
               split_prev = true;
            } else if (!can_split_curr) {
               /* neither edge can be split - cannot make topology consistent */
               throw ex_invalid_argument(
                  "cannot subdivide indentical edges to be nonoverlapping"
               );
            }
            /* split edge */
            contour_edge* e = (split_prev) ? (&e_prev) : (&e_curr);
            unsigned long split_pos =
               find_contour_point_farthest_from_segment(*e, p, q);
            this->split_edge(e->id, split_pos);
            /* update index of previous edge */
            if (split_prev) { prev = curr; }
         } else {
            /* update index of previous edge */
            prev = curr;
         }
      }
   }
   /* repeatedly subdivide intersecting edges until topology is consistent */
   while (true) {
      /* create points corresponding to vertices */
      auto_collection< point_2D, array_list<point_2D> > points(
         new array_list<point_2D>()
      );
      unsigned long n_vertices = _vertices->size();
      for (unsigned long n = 0; n < n_vertices; n++) {
         const contour_vertex& v = (*_vertices)[n];
         auto_ptr<point_2D> p(
            new point_2D(static_cast<double>(v.x), static_cast<double>(v.y))
         );
         points->add(*p);
         p.release();
      }
      /* create segments corresponding to edges */
      unsigned long n_edges = _edges->size();
      array<unsigned long> segment_start_ids(n_edges);
      array<unsigned long> segment_end_ids(n_edges);
      for (unsigned long n = 0; n < n_edges; n++) {
         const contour_edge& e = (*_edges)[n];
         segment_start_ids[n] = e.vertex_start->id;
         segment_end_ids[n]   = e.vertex_end->id;
      }
      /* compute segment intersections */
      seg_intersect seg_int(*points, segment_start_ids, segment_end_ids);
      /* initialize edge modification flags */
      array<bool> do_split(n_edges, false);
      bool found_intersection = false;
      /* mark edges which must be split */
      unsigned long n_v = seg_int.vertices_size();
      for (unsigned long n = 0; n < n_v; n++) {
         /* skip endpoint only intersections */
         if (!(seg_int.is_interior_intersection(n)))
            continue;
         found_intersection = true;
         /* get segments whos interiors/endpoints pass through intersection */
         array<unsigned long> interior_ids = seg_int.interior_intersection(n);
         array<unsigned long> endpoint_ids = seg_int.endpoint_intersection(n);
         /* check intersection type */
         if (endpoint_ids.is_empty()) {
            /* intersection of edge interiors only - mark all but one edge */
            const contour_edge* ea = &((*_edges)[interior_ids[0]]);
            unsigned long n_interior = interior_ids.size();
            for (unsigned long ne = 1; ne < n_interior; ne++) {
               const contour_edge* eb = &((*_edges)[interior_ids[ne]]);
               /* check if contours can be split */
               if (ea->size() == 0) {
                  /* only eb can be split */
                  do_split[eb->id] = true;
                  continue;
               } else if (eb->size() == 0) {
                  /* only ea can be split */
                  do_split[ea->id] = true;
                  ea = eb;
                  continue;
               }
               /* get contour endpoints */
               point_2D pa = ea->vertex_start->point();
               point_2D qa = ea->vertex_end->point();
               point_2D pb = eb->vertex_start->point();
               point_2D qb = eb->vertex_end->point();
               /* compute distance from contour ea to segment defined by eb */
               unsigned long pos_a =
                  find_contour_point_closest_to_segment(*ea, pb, qb);
               double dist_a = ea->point(pos_a).distance_to_segment(pb, qb);
               /* compute distance from contour eb to segment defined by ea */
               unsigned long pos_b =
                  find_contour_point_closest_to_segment(*eb, pa, qa);
               double dist_b = eb->point(pos_b).distance_to_segment(pa, qa);
               /* compare distances */
               if (dist_a < dist_b) {
                  /* contour ea intersects segment eb - split contour eb */
                  do_split[eb->id] = true;
               } else {
                  /* contour eb intersects segment ea - split contour ea */
                  do_split[ea->id] = true;
                  ea = eb;
               }
            }
         } else {
            /* mark all edges with interiors passing through intersection */
            unsigned long n_interior = interior_ids.size();
            for (unsigned long ne = 0; ne < n_interior; ne++)
               do_split[interior_ids[ne]] = true;
         }
      }
      /* check if finished (no interior intersections) */
      if (!found_intersection)
         break;
      /* split marked edges */
      for (unsigned long n = 0; n < n_edges; n++) {
         if (do_split[n]) {
            /* get edge, check that it can be split */
            contour_edge& e = (*_edges)[n];
            if (e.size() == 0) {
               throw ex_invalid_argument(
                  "cannot subdivide edges to make topology consistent"
               );
            }
            /* split edge */
            point_2D p = e.vertex_start->point();
            point_2D q = e.vertex_end->point();
            unsigned long split_pos =
               find_contour_point_farthest_from_segment(e, p, q);
            this->split_edge(e.id, split_pos);
         }
      }
   }
   /* order edges in equivalence classes by traversal of original contour */
   this->sort_edges_equiv();
}

/***************************************************************************
 * Contour completion.
 ***************************************************************************/

/*
 * Add the four bounding corners of the image to the set of contour
 * vertices if they neither currently exist in the vertex set nor lie
 * on an existing contour.
 */
void lib_image::contour_set::add_bounding_vertices() {
   /* get image dimensions */
   unsigned long size_x = _assignments.size(0);
   unsigned long size_y = _assignments.size(1);
   /* check image dimensions */
   if ((size_x == 0) || (size_y == 0))
      return;
   /* assemble corner vertex coordinates */
   unsigned long xs[] = { 0,         0,  (size_x-1), (size_x-1) };
   unsigned long ys[] = { 0, (size_y-1),          0, (size_y-1) };
   /* check if corner vertices are covered */
   for (unsigned long n = 0; n < 4; n++) {
      unsigned long x = xs[n];
      unsigned long y = ys[n];
      if (!(_is_vertex(x,y) || _is_edge(x,y))) {
         /* vertex doesn't exist - create it */
         auto_ptr<contour_vertex> v = create_contour_vertex(x, y);
         v->id = _vertices->size();
         /* record vertex assignment */
         _assignments(x,y) = v->id;
         /* add vertex to vertex set */
         _vertices->add(*v);
         v.release();
         /* mark vertex */
         _is_vertex(x,y) = true;
      }
   }
}

/*
 * Completion edges from constrained Delaunay triangulation (CDT).
 *
 * Compute the CDT of the straight line segment approximation of the
 * contour edges.  Update the contour set to include all new edges
 * appearing in the CDT as completion edges.
 *
 * Optionally, return the CDT vertices (points which correspond to
 * the contour vertices), and the CDT data structure itself (which
 * references the CDT vertices).
 *
 * Note that the line segment approximation must respect the true
 * contour topology in order to guarantee the existence of the CDT.
 * The subdivide_global() method enforces this property and should 
 * be called prior to computing the CDT.
 */
void lib_image::contour_set::completion_cdt() {
   auto_collection< point_2D, array_list<point_2D> > cdt_vertices;
   auto_ptr<triangulation> cdt;
   this->completion_cdt(cdt_vertices, cdt);
}

void lib_image::contour_set::completion_cdt(
   auto_collection< point_2D, array_list<point_2D> >& cdt_vertices,
   auto_ptr<triangulation>&                           cdt)
{
   /* create points corresponding to contour vertices */
   auto_collection< point_2D, array_list<point_2D> > points(
      new array_list<point_2D>()
   );
   unsigned long n_vertices = _vertices->size();
   for (unsigned long n = 0; n < n_vertices; n++) {
      const contour_vertex& v = (*_vertices)[n];
      auto_ptr<point_2D> p(
         new point_2D(static_cast<double>(v.x), static_cast<double>(v.y))
      );
      points->add(*p);
      p.release();
   }
   /* create constraints corresponding to contour edges */
   unsigned long n_edges = _edges->size();
   array<unsigned long> constraint_start_ids(n_edges);
   array<unsigned long> constraint_end_ids(n_edges);
   for (unsigned long n = 0; n < n_edges; n++) {
      const contour_edge& e = (*_edges)[n];
      constraint_start_ids[n] = e.vertex_start->id;
      constraint_end_ids[n]   = e.vertex_end->id;
   }
   /* compute cdt */
   auto_ptr<triangulation> t = triangulation::delaunay(
      *points, constraint_start_ids, constraint_end_ids
   );
   /* create contour edges for cdt completion edges */
   unsigned long n_edges_cdt = t->edges_size();
   for (unsigned long n = n_edges; n < n_edges_cdt; n++) {
      /* create edge and set its identity */
      auto_ptr<contour_edge> e = create_contour_edge(
         (*_vertices)[t->edge_vertex_id(n,0)],
         (*_vertices)[t->edge_vertex_id(n,1)]
      );
      e->id               = n;
      e->contour_equiv_id = _edges_equiv->size();
      e->is_completion    = true;
      /* compute coordinates of points along edge */
      const point_2D& p_start = t->vertex(e->vertex_start->id);
      const point_2D& p_end   = t->vertex(e->vertex_end->id);
      unsigned long size_y = this->size_y();
      auto_ptr< array<unsigned long> > inds = scan_line_segment(
         p_start, p_end, size_y, true
      );
      unsigned long n_inds = inds->size();
      unsigned long n_coords = (n_inds > 2) ? (n_inds - 2) : 0;
      e->x_coords.resize(n_coords);
      e->y_coords.resize(n_coords);
      for (unsigned long n = 0; n < n_coords; n++) {
         unsigned long ind = (*inds)[n+1];
         e->x_coords[n] = ind / size_y;
         e->y_coords[n] = ind % size_y;
      }
      /* create edge equivalence class */
      auto_ptr< list<contour_edge> > e_equiv(new list<contour_edge>());
      e_equiv->add(*e);
      /* add edge */
      _edges->add(*e);             e.release();
      _edges_equiv->add(*e_equiv); e_equiv.release();
   }
   /* return cdt */
   cdt_vertices = points;
   cdt = t;
}

/***************************************************************************
 * Contour traversal state machine.
 ***************************************************************************/

namespace {
/*
 * A contour traversal is a directed traversal of the boundary of the pixels
 * lying on the contour.  The direction of traversal specifies whether the 
 * left or right side (as viewed when moving from the start vertex to the end
 * vertex) of the contour is traced out.
 *
 * At each step of the traversal, the current state is defined by the position
 * along the contour (a pixel), and an oriented boundary of that pixel.  Each
 * pixel has eight oriented boundaries, as shown below in Figure 1.
 *
 * Successive pixels along the contour must be neighbors in a 3x3 grid, and the
 * state machine makes use of the labeling of these neighbors shown in Figure 2.
 *
 *                     0                       4
 *                    --+                     +--
 *               3  +     |  1            5 |     + 7
 *                  |     +                 +     |
 *                    +--                     --+
 *                     2                       6
 *
 *             clockwise states      counterclockwise states
 *
 * Figure 1: Labeling of states representing the directed boundary of a pixel.
 *
 *   
 *                            2   5   8 
 *                       ^ 
 *                       |    1   *   7
 *                       y  
 *                            0   3   6
 * 
 *                              x -->
 *
 * Figure 2: Labeling of neighbors offset from the center pixel (*).
 */

/*
 * Given the coordinates of two adjacent contour points, compute the label 
 * describing the relative position of the second with respect to the first.
 * The labeling scheme is shown above in Figure 2.
 */
unsigned long contour_traversal_compute_offset(
   unsigned long x, unsigned long y, unsigned long x_next, unsigned long y_next)
{
   unsigned long x_offset = 1 + x_next - x;
   unsigned long y_offset = 1 + y_next - y;
   return (x_offset*3 + y_offset);
}

/*
 * Compute the label of the next element along the contour given the index
 * of the current element.  If the current element is the last point on the
 * contour, compute the offset of the end vertex.
 */
unsigned long contour_traversal_compute_offset(
   const lib_image::contour_edge& e,      /* contour edge */
   unsigned long                  ind)    /* current position */
{
   unsigned long n_coords = e.size();
   unsigned long ind_next = ind + 1;
   if (ind_next < n_coords) {
      return contour_traversal_compute_offset(
         e.x_coords[ind],      e.y_coords[ind],
         e.x_coords[ind_next], e.y_coords[ind_next]
      );
   } else {
      return contour_traversal_compute_offset(
         e.x_coords[ind],      e.y_coords[ind],
         e.vertex_end->x,      e.vertex_end->y
      );
   }
}

/*
 * Return the initial state for traversal of the left side of the contour.
 */
int contour_traversal_init_state_left(unsigned long offset) {
   static const int init_state_map[] = { 3, 0, 0, 3, -1, 1, 2, 2, 1 };
   return init_state_map[offset];
}

/*
 * Return the initial state for traversal of the left side of the contour.
 */
int contour_traversal_init_state_left(const lib_image::contour_edge& e) {
   return contour_traversal_init_state_left(
      contour_traversal_compute_offset(
         e.x_coords[0], e.y_coords[0], e.vertex_start->x, e.vertex_start->y
      )
   );
}

/*
 * Return the initial state for traversal of the right side of the contour.
 */
int contour_traversal_init_state_right(unsigned long offset) {
   static const int init_state_map[] = { 6, 6, 5, 7, -1, 5, 7, 4, 4 };
   return init_state_map[offset];
}

/*
 * Return the initial state for traversal of the right side of the contour.
 */
int contour_traversal_init_state_right(const lib_image::contour_edge& e) {
   return contour_traversal_init_state_right(
      contour_traversal_compute_offset(
         e.x_coords[0], e.y_coords[0], e.vertex_start->x, e.vertex_start->y
      )
   );
}
 
/*
 * Move to the next boundary element in a contour traversal.
 */
void contour_traversal_update_state(
   unsigned long& ind,                    /* current position */
   int&           state,                  /* current state */
   unsigned long  offset)                 /* offset of next element */
{
   /* contour traversal state machine */
   static const int state_rotate[] = { 1, 2, 3, 0, 5, 6, 7, 4 };
   static const int state_lift[]   = { 3, 0, 1, 2, 7, 4, 5, 6 };
   /* offset indices which correspond to state changes */
   static const unsigned long offset_continue[]  = { 7, 3, 1, 5, 1, 3, 7, 5 };
   static const unsigned long offset_lift_head[] = { 8, 6, 0, 2, 2, 0, 6, 8 };
   static const unsigned long offset_lift_tail[] = { 5, 7, 3, 1, 5, 1, 3, 7 };
   /* update state */
   if (offset == 4) {
      /* contour contains duplicate pixels */
      ind++;
   } else if (offset_continue[state] == offset) {
      /* next pixel on contour is aligned with current */
      ind++;
   } else if ((offset_lift_head[state] == offset) ||
              (offset_lift_tail[state] == offset)) {
      /* move to next boundary pixel and update boundary direction */
      ind++;
      state = state_lift[state];
   } else {
      /* rotate around boundary of current pixel */
      state = state_rotate[state];
   }
}
} /* namespace */

/***************************************************************************
 * Region extraction (2D).
 ***************************************************************************/

namespace {
/*
 * Scan rays in the direction of the normal from each point along a traversal
 * of the edge.  The returned rays sweep out a region of the specified width.
 */
auto_collection< array<unsigned long> > extract_edge_side(
   const lib_image::contour_edge&   e,       /* contour edge */
   bool                             side,    /* false = left, true = right */
   double                           width,   /* region width */
   unsigned long                    size_x,  /* grid size in x-direction */
   unsigned long                    size_y)  /* grid size in y-direction */
{
   /* allocate collection of rays scanned from the edge */
   auto_collection< array<unsigned long> > rays(
      new list< array<unsigned long> >()
   );
   /* check that contour is nonempty */
   unsigned long n_points = e.size();
   if (n_points > 0) {
      /* compute angle of normal to contour side */
      point_2D p_start = e.vertex_start->point();
      point_2D p_end   = e.vertex_end->point();
      double theta = arg(p_end - p_start) + M_PI_2l;
      if (side) { theta += M_PIl; }
      /* compute angular change for moving ~1/8 of a grid cell along boundary */
      double delta = 0.125 * (2*M_PIl) / (math::ceil(2*M_PIl*width)+1);
      /* compute rays bounding a thin cone in the direction to scan */
      point_2D p_ray_prev(
         width * math::cos(theta - delta), width * math::sin(theta - delta)
      );
      point_2D p_ray_next(
         width * math::cos(theta + delta), width * math::sin(theta + delta)
      );
      /* initialize traversal state */
      int state = side ?
         contour_traversal_init_state_right(e)
       : contour_traversal_init_state_left(e);
      /* traverse edge */
      unsigned long ind = 0;
      while ((ind < n_points) && (state != -1)) {
         /* get current point on edge */
         point_2D p = e.point(ind);
         /* move to adjacent point on correct side of boundary */
         if ((state == 0) || (state == 4)) {
            (p.y())++;
         } else if ((state == 1) || (state == 7)) {
            (p.x())++;
         } else if ((state == 2) || (state == 6)) {
            (p.y())--;
         } else /* ((state == 3) || (state == 5)) */ {
            (p.x())--;
         }
         /* check that point lies within grid */
         if (is_in_grid(p, size_x, size_y)) {
            /* compute segment end points */
            point_2D q_prev = p + p_ray_prev;
            point_2D q_next = p + p_ray_next;
            clip_line_segment_to_grid(p, q_prev, size_x, size_y);
            clip_line_segment_to_grid(p, q_next, size_x, size_y);
            /* scan segments */
            auto_ptr< array<unsigned long> > ray_prev = scan_line_segment(
               p, q_prev, size_y, false
            );
            auto_ptr< array<unsigned long> > ray_next = scan_line_segment(
               p, q_next, size_y, false
            );
            /* store scanned rays */
            rays->add(*ray_prev); ray_prev.release();
            rays->add(*ray_next); ray_next.release();
         }
         /* update traversal state */
         unsigned long offset = contour_traversal_compute_offset(e, ind);
         contour_traversal_update_state(ind, state, offset);
      }
   }
   return rays;
}

/*
 * Scan rays in a full disc at the given subdivision vertex along a contour.
 */
auto_collection< array<unsigned long> > extract_vertex_side(
   const lib_image::contour_vertex&    v,       /* subdivision vertex */
   unsigned long x_prev, unsigned long y_prev,  /* previous point on contour */
   unsigned long x_next, unsigned long y_next,  /* next point on contour */
   bool                                side,    /* false = left, true = right */
   double                              width,   /* region width */
   unsigned long                       size_x,  /* grid size in x-direction */
   unsigned long                       size_y)  /* grid size in y-direction */
{
   /* allocate collection of rays scanned from the vertex */
   auto_collection< array<unsigned long> > rays(
      new list< array<unsigned long> >()
   );
   /* compute offsets of previous/next contour points */
   unsigned long offset_prev = contour_traversal_compute_offset(
      v.x, v.y, x_prev, y_prev
   );
   unsigned long offset_next = contour_traversal_compute_offset(
      v.x, v.y, x_next, y_next
   );
   /* initialize traversal state along boundary of vertex */
   int state = side ?
      contour_traversal_init_state_right(offset_prev)
    : contour_traversal_init_state_left(offset_prev);
   /* traverse vertex boundary */
   unsigned long ind = 0;
   while ((ind == 0) && (state != -1)) {
      /* get point corresponding to vertex */
      point_2D p = v.point();
      /* move to adjacent point on correct side of boundary */
      if ((state == 0) || (state == 4)) {
         (p.y())++;
      } else if ((state == 1) || (state == 7)) {
         (p.x())++;
      } else if ((state == 2) || (state == 6)) {
         (p.y())--;
      } else /* ((state == 3) || (state == 5)) */ {
         (p.x())--;
      }
      /* check that point lies within grid */
      if (is_in_grid(p, size_x, size_y)) {
         /* scan disc */
         auto_collection< array<unsigned long> > disc_rays = scan_disc(
            p, width, size_x, size_y
         );
         rays->add(*disc_rays);
         delete (disc_rays.release());
      }
      /* update traversal state */
      contour_traversal_update_state(ind, state, offset_next);
   }
   return rays;
}

/*
 * Given a collection of scanned rays (each containing indices of grid elements
 * that lie on a line segment), extract the visible region.  Traversing each ray
 * until it is blocked by one of the following conditions carves out the visible
 * region.  A ray is blocked at the first point it hits which (1) has the
 * specified id or (2) has visibility greater than the specified visibility.
 */
void extract_visible_region(
   const collection< array<unsigned long> >& rays,       /* scanned rays */
   unsigned long                             id,         /* current id */
   const matrix<unsigned long>&              id_map,     /* id map */
   double                                    vis,        /* current vis */
   const matrix<>&                           visibility, /* visibility map */
   array<unsigned long>&                     region)     /* returned region */
{
   /* traverse rays, computing where they are blocked */
   unsigned long n_rays = rays.size();
   array<unsigned long> ray_len(n_rays);
   unsigned long total_len = 0;
   {
      auto_ptr< iterator< array<unsigned long> > > i = rays.iter_create();
      for (unsigned long n = 0; n < n_rays; n++) {
         /* get ray */
         const array<unsigned long>& ray = i->next();
         unsigned long len = ray.size();
         /* traverse ray */
         for (unsigned long pos = 0; pos < len; pos++) {
            unsigned long ind = ray[pos];
            if ((id_map[ind] == id) || (visibility[ind] > vis))
               len = pos;
         }
         ray_len[n] = len;
         total_len += len;
      }
   }
   /* assemble region */
   region.resize(total_len);
   {
      unsigned long r_pos = 0;
      auto_ptr< iterator< array<unsigned long> > > i = rays.iter_create();
      for (unsigned long n = 0; n < n_rays; n++) {
         /* get ray */
         const array<unsigned long>& ray = i->next();
         unsigned long len = ray_len[n];
         /* add visible portion to region */
         for (unsigned long pos = 0; pos < len; pos++)
            region[r_pos++] = ray[pos];
      }
   }
   /* compute unique region membership */
   region.unique();
}

/*
 * Edge region extractor.  
 */
class edge_region_extractor : public runnable {
public:
   /*
    * Constructor.
    */
   explicit edge_region_extractor(
      const lib_image::contour_edge&   e,             /* contour edge */
      double                           width,         /* region width */
      const matrix<unsigned long>&     id_map,        /* id map */
      double                           vis,           /* edge visibility */
      const matrix<>&                  visibility,    /* visibility map */
      array<unsigned long>&            region_left,   /* returned left region */
      array<unsigned long>&            region_right)  /* returned right region */
    : _e(e), _width(width), _id_map(id_map), _vis(vis), _visibility(visibility),
      _region_left(region_left), _region_right(region_right)
   { }

   /*
    * Destructor.
    */
   virtual ~edge_region_extractor() { /* do nothing */ }

   /*
    * Extract regions on both sides of the edge.
    */
   virtual void run() {
      /* scan each side of the edge */
      unsigned long size_x = _id_map.size(0);
      unsigned long size_y = _id_map.size(1);
      auto_collection< array<unsigned long> > rays_left = extract_edge_side(
         _e, false, _width, size_x, size_y
      );
      auto_collection< array<unsigned long> > rays_right = extract_edge_side(
         _e, true,  _width, size_x, size_y
      );
      /* extract visible regions */
      unsigned long id = _e.contour_equiv_id;
      extract_visible_region(
         *rays_left,  id, _id_map, _vis, _visibility, _region_left
      );
      extract_visible_region(
         *rays_right, id, _id_map, _vis, _visibility, _region_right
      );
   }

protected:
   const lib_image::contour_edge&   _e;               /* contour edge */
   double                           _width;           /* region width */
   const matrix<unsigned long>&     _id_map;          /* id map */
   double                           _vis;             /* edge visibility */
   const matrix<>&                  _visibility;      /* visibility map */
   array<unsigned long>&            _region_left;     /* returned left region */
   array<unsigned long>&            _region_right;    /* returned right region */
};

/*
 * Vertex region extractor.
 */
class vertex_region_extractor : public runnable {
public:
   /*
    * Constructor.
    */
   explicit vertex_region_extractor(
      const lib_image::contour_vertex& v,             /* contour vertex */
      const lib_image::contour_edge&   e_prev,        /* previous edge */
      const lib_image::contour_edge&   e_next,        /* next edge */
      double                           width,         /* region width */
      const matrix<unsigned long>&     id_map,        /* id map */
      double                           vis,           /* edge visibility */
      const matrix<>&                  visibility,    /* visibility map */
      array<unsigned long>&            region_left,   /* returned left region */
      array<unsigned long>&            region_right)  /* returned right region */
    : _v(v), _e_prev(e_prev), _e_next(e_next),
      _width(width), _id_map(id_map), _vis(vis), _visibility(visibility),
      _region_left(region_left), _region_right(region_right)
   { }

   /*
    * Destructor.
    */
   virtual ~vertex_region_extractor() { /* do nothing */ }

   /*
    * Extract regions on both sides of the subdivision vertex.
    */
   virtual void run() {
      /* get previous/next vertices */
      const lib_image::contour_vertex* v_prev = _e_prev.vertex_start;
      const lib_image::contour_vertex* v_next = _e_next.vertex_end;
      /* compute coordinates of previous point on contour */
      unsigned long n_points_prev = _e_prev.size(); 
      bool use_e_prev = (n_points_prev > 0);
      unsigned long x_prev =
         use_e_prev ? _e_prev.x_coords[n_points_prev-1] : v_prev->x;
      unsigned long y_prev =
         use_e_prev ? _e_prev.y_coords[n_points_prev-1] : v_prev->y;
      /* compute coordinates of next point on contour */
      unsigned long n_points_next = _e_next.size(); 
      bool use_e_next = (n_points_next > 0);
      unsigned long x_next =
         use_e_next ? _e_next.x_coords[0] : v_next->x;
      unsigned long y_next =
         use_e_next ? _e_next.y_coords[0] : v_next->y;
      /* scan each side of the vertex */
      unsigned long size_x = _id_map.size(0);
      unsigned long size_y = _id_map.size(1);
      auto_collection< array<unsigned long> > rays_left = extract_vertex_side(
         _v, x_prev, y_prev, x_next, y_next, false, _width, size_x, size_y
      );
      auto_collection< array<unsigned long> > rays_right = extract_vertex_side(
         _v, x_prev, y_prev, x_next, y_next, true,  _width, size_x, size_y
      );
      /* extract visible regions */
      unsigned long id = _e_prev.contour_equiv_id;
      extract_visible_region(
         *rays_left,  id, _id_map, _vis, _visibility, _region_left
      );
      extract_visible_region(
         *rays_right, id, _id_map, _vis, _visibility, _region_right
      );
   }

protected:
   const lib_image::contour_vertex& _v;               /* contour vertex */
   const lib_image::contour_edge&   _e_prev;          /* previous edge */
   const lib_image::contour_edge&   _e_next;          /* next edge */
   double                           _width;           /* region width */
   const matrix<unsigned long>&     _id_map;          /* id map */
   double                           _vis;             /* edge visibility */
   const matrix<>&                  _visibility;      /* visibility map */
   array<unsigned long>&            _region_left;     /* returned left region */
   array<unsigned long>&            _region_right;    /* returned right region */
};

/*
 * Region combiner.
 */
class region_combiner : public runnable {
public:
   /*
    * Constructor.
    */
   explicit region_combiner(
      const list<lib_image::contour_edge>& equiv_set, /* contour equiv set */
      const array< array<unsigned long> >& v_regions, /* vertex regions */
      const array< array<unsigned long> >& e_regions, /* edge regions */
      array<unsigned long>&                region)    /* combined region */
    : _equiv_set(equiv_set),
      _v_regions(v_regions),
      _e_regions(e_regions),
      _region(region)
   { }

   /*
    * Destructor.
    */
   virtual ~region_combiner() { /* do nothing */ }

   /*
    * Combine regions associated with a subdivided contour to produce the
    * region associated with the original contour before subdivision.
    */
   virtual void run() {
      /* compute total size of region */
      unsigned long r_size = 0;
      for (list<lib_image::contour_edge>::iterator_t i(_equiv_set);
           i.has_next(); )
      {
         const lib_image::contour_edge& e = i.next();
         r_size += _e_regions[e.id].size();
         if (i.has_next())
            r_size += _v_regions[e.vertex_end->id].size();
      }
      /* resize region */
      _region.resize(r_size);
      /* copy membership of component regions */
      unsigned long r_pos = 0; 
      for (list<lib_image::contour_edge>::iterator_t i(_equiv_set);
           i.has_next(); )
      {
         const lib_image::contour_edge& e = i.next();
         const array<unsigned long>& e_region = _e_regions[e.id];
         for (unsigned long n = 0, size = e_region.size(); n < size; n++)
            _region[r_pos++] = e_region[n];
         if (i.has_next()) {
            const array<unsigned long>& v_region = _v_regions[e.vertex_end->id];
            for (unsigned long n = 0, size = v_region.size(); n < size; n++)
               _region[r_pos++] = v_region[n];
         }
      }
      /* compute unique region membership */
      _region.unique();
   }

protected:
   const list<lib_image::contour_edge>& _equiv_set;   /* contour equiv set */
   const array< array<unsigned long> >& _v_regions;   /* vertex regions */
   const array< array<unsigned long> >& _e_regions;   /* edge regions */
   array<unsigned long>&                _region;      /* combined region */
};
} /* namespace */

/*
 * Extract the regions lying on the left and right sides of each of the 
 * given contours.
 * 
 * The width of both the left and right regions is determined as the
 * specified fraction of the length of the corresponding contour edge,
 * with a minimum width used for short edges.  In this manner, region
 * size scales with contour length.
 */
lib_image::region_set::region_set(
   const contour_set& contours,
   double             width,
   double             width_min)
{
   /* set all visibilities to zero (so no boundaries are seen) */
   array<double> edge_vis(contours.edges_size());
   matrix<> visibility(contours.size_x(), contours.size_y());
   /* extract regions */
   this->extract_regions(contours, width, width_min, edge_vis, visibility);
}

/*
 * Extract regions bounded by intervening contours.
 *
 * When computing region membership, take into account boundaries implied
 * by other non-completion contours in the following manner.  Denote the
 * specified relative scale sensitivity by s.  The region surrounding
 * contour C extends to the given width except where it is blocked by
 * a non-completion contour C' such that:
 *
 *                         L(C') > s * L(C)
 *
 * where L(*) denotes contour length.
 * 
 * In other words, a contour sees all boundaries implied by other contours
 * whose length is greater than s times larger than its own length.
 *
 * Setting s = 0 means each contour sees all other contours.  Setting s
 * to infinity is equivalent to the above version of region extraction,
 * in which each contour ignores all others.
 */
lib_image::region_set::region_set(
   const contour_set& contours,
   double             width,
   double             width_min,
   double             scale)
{
   /* compute length of each original contour (equivalence class) */
   array<double> edge_equiv_len(contours.edges_equiv_size());
   unsigned long n_edges = contours.edges_size();
   for (unsigned long e_id = 0; e_id < n_edges; e_id++) {
      const contour_edge& e = (*(contours._edges))[e_id];
      edge_equiv_len[e.contour_equiv_id] += e.length();
   }
   /* set edge and pixel visibilities to contour lengths */
   array<double> edge_vis(n_edges);
   matrix<> visibility(contours.size_x(), contours.size_y());
   for (unsigned long e_id = 0; e_id < n_edges; e_id++) {
      /* set edge visibility */
      const contour_edge& e = (*(contours._edges))[e_id];
      double vis = edge_equiv_len[e.contour_equiv_id];
      edge_vis[e_id] = scale * vis;
      /* skip setting pixel visibilities of completion edges */
      if (e.is_completion)
         continue;
      /* set pixel visibility for points on edge */
      unsigned long n_points = e.size();
      for (unsigned long n = 0; n < n_points; n++)
         visibility(e.x_coords[n], e.y_coords[n]) = vis;
      /* set pixel visibility for edge vertices */
      if (visibility(e.vertex_start->x, e.vertex_start->y) < vis)
         visibility(e.vertex_start->x, e.vertex_start->y) = vis;
      if (visibility(e.vertex_end->x, e.vertex_end->y) < vis)
         visibility(e.vertex_end->x, e.vertex_end->y) = vis;
   }
   /* extract regions */
   this->extract_regions(contours, width, width_min, edge_vis, visibility);
}

/*
 * Extract regions using a custom intervening contour cue.
 *
 * Instead of using contour length to determine visibility, define a per
 * pixel visibility map and a minimum visibility for each contour.  The
 * region surrounding a contour C extends to the given width except where
 * it is blocked by a pixel p such that:
 *
 *                         V(p) > V(C)
 *
 * where V(*) denotes visibility.
 */
lib_image::region_set::region_set(
   const contour_set&   contours,
   double               width,
   double               width_min,
   const array<double>& edge_vis,
   const matrix<>&      visibility)
{
   this->extract_regions(contours, width, width_min, edge_vis, visibility);
}

/*
 * Extract regions using a custom intervening contour cue.
 */
void lib_image::region_set::extract_regions(
   const contour_set&   contours,
   double               width,
   double               width_min,
   const array<double>& edge_vis,
   const matrix<>&      visibility)
{
   /* allocate map of pixel -> id of contour equivalence class */
   unsigned long n_equiv = contours.edges_equiv_size();
   matrix<unsigned long> id_map(contours.size_x(), contours.size_y(), n_equiv);
   /* set id map entries */
   unsigned long n_edges = contours.edges_size();
   for (unsigned long e_id = 0; e_id < n_edges; e_id++) {
      /* get edge, skip if completion */
      const contour_edge& e = (*(contours._edges))[e_id];
      if (e.is_completion)
         continue;
      /* set id for points on edge */
      unsigned long n_points = e.size();
      for (unsigned long n = 0; n < n_points; n++)
         id_map(e.x_coords[n], e.y_coords[n]) = e.contour_equiv_id;
      /* set id for edge vertices (if subdivision) */
      if (e.vertex_start->is_subdivision)
         id_map(e.vertex_start->x, e.vertex_start->y) = e.contour_equiv_id;
      if (e.vertex_end->is_subdivision)
         id_map(e.vertex_end->x, e.vertex_end->y) = e.contour_equiv_id;
   }
   /* compute length of each original contour (equivalence class) */
   array<double> edge_equiv_len(contours.edges_equiv_size());
   for (unsigned long e_id = 0; e_id < n_edges; e_id++) {
      const contour_edge& e = (*(contours._edges))[e_id];
      edge_equiv_len[e.contour_equiv_id] += e.length();
   }
   /* resize arrays to hold returned vertex and edge regions */
   unsigned long n_vertices = contours.vertices_size();
   _vertex_regions_left.resize(n_vertices);
   _vertex_regions_right.resize(n_vertices);
   _edge_regions_left.resize(n_edges);
   _edge_regions_right.resize(n_edges);
   /* allocate collection of region extractors */
   auto_collection< runnable, list<runnable> > region_extractors(
      new list<runnable>()
   );
   /* setup region extractors */
   for (unsigned long equiv_id = 0; equiv_id < n_equiv; equiv_id++) {
      /* get contour equivalence set */
      const list<contour_edge>& equiv_set =
         (*(contours._edges_equiv))[equiv_id];
      /* compute absolute region width to use for contour */
      double w = width * edge_equiv_len[equiv_id];
      if (w < width_min) { w = width_min; }
      /* create region extractors for edges and subdivision vertices */
      const contour_edge* e_prev = NULL;
      for (list<contour_edge>::iterator_t i(equiv_set); i.has_next(); ) {
         /* get current edge */
         const contour_edge& e = i.next();
         /* create extractor for edge region */
         {
            auto_ptr<edge_region_extractor> e_region_extractor(
               new edge_region_extractor(
                  e, w, id_map, edge_vis[e.id], visibility, 
                  _edge_regions_left[e.id], _edge_regions_right[e.id]
               )
            );
            region_extractors->add(*e_region_extractor);
            e_region_extractor.release();
         }
         /* create extractor for vertex region (if subdivision) */
         if ((e.vertex_start->is_subdivision) && (e_prev != NULL)) {
            /* compute vertex visibility */
            double vis_e_prev = edge_vis[e_prev->id];
            double vis_e      = edge_vis[e.id];
            double vis = (vis_e_prev > vis_e) ? vis_e_prev : vis_e;
            /* create region extractor */
            auto_ptr<vertex_region_extractor> v_region_extractor(
               new vertex_region_extractor(
                  *(e.vertex_start), *e_prev, e, w, id_map, vis, visibility, 
                  _vertex_regions_left[e.vertex_start->id],
                  _vertex_regions_right[e.vertex_start->id]
               )
            );
            region_extractors->add(*v_region_extractor);
            v_region_extractor.release();
         }
         /* record previous edge */
         e_prev = &e;
      }
   }
   /* extract regions */
   child_thread::run(*region_extractors);
   /* resize arrays to hold returned contour regions */
   _edge_equiv_regions_left.resize(n_equiv);
   _edge_equiv_regions_right.resize(n_equiv);
   /* allocate collection of region combiners */
   auto_collection< runnable, list<runnable> > region_combiners(
      new list<runnable>()
   );
   /* setup region combiners */
   for (unsigned long equiv_id = 0; equiv_id < n_equiv; equiv_id++) {
      /* create combiner for left side of contour */
      auto_ptr<region_combiner> region_left_combiner(
         new region_combiner(
            (*(contours._edges_equiv))[equiv_id],
            _vertex_regions_left,
            _edge_regions_left,
            _edge_equiv_regions_left[equiv_id]
         )
      );
      region_combiners->add(*region_left_combiner);
      region_left_combiner.release();
      /* create combiner for right side of contour */
      auto_ptr<region_combiner> region_right_combiner(
         new region_combiner(
            (*(contours._edges_equiv))[equiv_id],
            _vertex_regions_right,
            _edge_regions_right,
            _edge_equiv_regions_right[equiv_id]
         )
      );
      region_combiners->add(*region_right_combiner);
      region_right_combiner.release();
   }
   /* combine vertex and edge regions associated with the same contour */
   child_thread::run(*region_combiners);
}

/*
 * Copy constructor.
 */
lib_image::region_set::region_set(const region_set& regions)
 : _vertex_regions_left(regions._vertex_regions_left),
   _vertex_regions_right(regions._vertex_regions_right),
   _edge_regions_left(regions._edge_regions_left),
   _edge_regions_right(regions._edge_regions_right),
   _edge_equiv_regions_left(regions._edge_equiv_regions_left),
   _edge_equiv_regions_right(regions._edge_equiv_regions_right)
{ }

/*
 * Destructor.
 */
lib_image::region_set::~region_set() {
   /* do nothing */
}

/***************************************************************************
 * Region access.
 ***************************************************************************/

/*
 * Get number of vertices in the corresponding contour set.
 */
unsigned long lib_image::region_set::vertices_size() const {
   return _vertex_regions_left.size();
}

/*
 * Get number of edges in the region set.
 */
unsigned long lib_image::region_set::edges_size() const {
   return _edge_regions_left.size();
}

/*
 * Get the number of edge equivalence classes.
 * Each equivalence class represents a contour before subdivision.
 */
unsigned long lib_image::region_set::edges_equiv_size() const {
   return _edge_equiv_regions_left.size();
}

/*
 * Return the region (by reference) on the left side of the vertex. 
 * Only subdivision vertices are associated with nonempty regions.
 */
const array<unsigned long>& lib_image::region_set::vertex_region_left(
   unsigned long v_id) const
{
   if (v_id >= _vertex_regions_left.size())
      throw ex_index_out_of_bounds("invalid contour vertex id", v_id);
   return _vertex_regions_left[v_id];
}

/*
 * Return the region (by reference) on the right side of the vertex. 
 * Only subdivision vertices are associated with nonempty regions.
 */
const array<unsigned long>& lib_image::region_set::vertex_region_right(
   unsigned long v_id) const
{
   if (v_id >= _vertex_regions_right.size())
      throw ex_index_out_of_bounds("invalid contour vertex id", v_id);
   return _vertex_regions_right[v_id];
}

/*
 * Return the region (by reference) on the left side of the edge.
 */
const array<unsigned long>& lib_image::region_set::edge_region_left(
   unsigned long e_id) const
{
   if (e_id >= _edge_regions_left.size())
      throw ex_index_out_of_bounds("invalid contour edge id", e_id);
   return _edge_regions_left[e_id];
}

/*
 * Return the region (by reference) on the right side of the edge.
 */
const array<unsigned long>& lib_image::region_set::edge_region_right(
   unsigned long e_id) const
{
   if (e_id >= _edge_regions_right.size())
      throw ex_index_out_of_bounds("invalid contour edge id", e_id);
   return _edge_regions_right[e_id];
}

/*
 * Return the region (by reference) on the left side of the given
 * edge equivalence class (original contour before subdivision).
 */
const array<unsigned long>& lib_image::region_set::edge_equiv_region_left(
   unsigned long equiv_id) const
{
   if (equiv_id >= _edge_equiv_regions_left.size())
      throw ex_index_out_of_bounds("invalid contour equivalence id", equiv_id);
   return _edge_equiv_regions_left[equiv_id];
}

/*
 * Return the region (by reference) on the right side of the given
 * edge equivalence class (original contour before subdivision).
 */
const array<unsigned long>& lib_image::region_set::edge_equiv_region_right(
   unsigned long equiv_id) const
{
   if (equiv_id >= _edge_equiv_regions_right.size())
      throw ex_index_out_of_bounds("invalid contour equivalence id", equiv_id);
   return _edge_equiv_regions_right[equiv_id];
}

matrix<> lib_image::line_inds(
   double x0, double y0, double x1, double y1, double sx)
{
   point_2D p_start(x0,y0);
   point_2D p_end(x1,y1);
   auto_ptr< array<unsigned long> > a = scan_line_segment(p_start,p_end,static_cast<unsigned long>(sx),false);
   unsigned long s = a->size();
   matrix<> m(s,1);
   for (unsigned long n = 0; n < s; n++)
      m[n] = static_cast<double>((*a)[n]);
   return m;
}

} /* namespace libraries */
} /* namespace math */
