/*
 * Image processing library.
 */
#ifndef MATH__LIBRARIES__LIB_IMAGE_HH
#define MATH__LIBRARIES__LIB_IMAGE_HH

#include "collections/abstract/collection.hh"
#include "collections/array_list.hh"
#include "collections/list.hh"
#include "collections/pointers/auto_collection.hh"
#include "functors/distanceable_functors.hh"
#include "lang/array.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"
#include "math/geometry/triangulation.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "math/matrices/cmatrix.hh"
#include "math/matrices/functors/matrix_distance_functors.hh"
#include "mlearning/clustering/clusterers/abstract/centroid_clusterer.hh"

namespace math {
namespace libraries {
/*
 * Imports.
 */
using collections::abstract::collection;
using collections::array_list;
using collections::list;
using collections::pointers::auto_collection;
using functors::distanceable_functor;
using lang::array;
using lang::pointers::auto_ptr;
using math::geometry::point_2D;
using math::geometry::triangulation;
using math::matrices::matrix;
using math::matrices::cmatrix;
using math::matrices::functors::matrix_distance_functors;
using mlearning::clustering::clusterers::abstract::centroid_clusterer;

/*
 * Image processing functions.
 */
class lib_image {
public:
   /************************************************************************
    * Image color space transforms.
    *
    * Input RGB images should be scaled so that range of possible values for
    * each color channel is [0, 1].
    ************************************************************************/

   /*
    * Compute a grayscale image from an RGB image.
    */
   static matrix<> grayscale(
      const matrix<>&,  /* r */
      const matrix<>&,  /* g */
      const matrix<>&   /* b */
   );

   /*
    * Normalize a grayscale image so that intensity values lie in [0,1].
    */
   static void grayscale_normalize(matrix<>&);
   
   /*
    * Normalize and stretch a grayscale image so that intensity values span 
    * the full [0,1] range.
    */
   static void grayscale_normalize_stretch(matrix<>&);

   /*
    * Gamma correct the RGB image using the given correction value.
    */
   static void rgb_gamma_correct(
      matrix<>&,  /* r */
      matrix<>&,  /* g */
      matrix<>&,  /* b */
      double      /* gamma */
   );
   
   /*
    * Normalize an Lab image so that values for each channel lie in [0,1].
    */
   static void lab_normalize(
      matrix<>&,  /* l */
      matrix<>&,  /* a */
      matrix<>&   /* b */
   );

   /*
    * Convert from RGB color space to XYZ color space.
    */
   static void rgb_to_xyz(
      matrix<>&,  /* r (input) -> x (output) */
      matrix<>&,  /* g (input) -> y (output) */
      matrix<>&   /* b (input) -> z (output) */
   );
   
   /*
    * Convert from RGB color space to Lab color space.
    */
   static void rgb_to_lab(
      matrix<>&,  /* r (input) -> l (output) */
      matrix<>&,  /* g (input) -> a (output) */
      matrix<>&   /* b (input) -> b (output) */
   );

   /*
    * Convert from XYZ color space to RGB color space.
    */
   static void xyz_to_rgb(
      matrix<>&,  /* x (input) -> r (output) */
      matrix<>&,  /* y (input) -> g (output) */
      matrix<>&   /* z (input) -> b (output) */
   );
   
   /*
    * Convert from XYZ color space to Lab color space.
    */
   static void xyz_to_lab(
      matrix<>&,  /* x (input) -> l (output) */
      matrix<>&,  /* y (input) -> a (output) */
      matrix<>&   /* z (input) -> b (output) */
   );

   /*
    * Convert from Lab color space to RGB color space.
    */
   static void lab_to_rgb(
      matrix<>&,  /* l (input) -> r (output) */
      matrix<>&,  /* a (input) -> g (output) */
      matrix<>&   /* b (input) -> b (output) */
   );
   
   /*
    * Convert from Lab color space to XYZ color space.
    */
   static void lab_to_xyz(
      matrix<>&,  /* l (input) -> x (output) */
      matrix<>&,  /* a (input) -> y (output) */
      matrix<>&   /* b (input) -> z (output) */
   );
   
   /************************************************************************
    * Matrix border operations (2D).
    *
    * Expand/remove the border of a 2D matrix.
    ************************************************************************/
  
   /*
    * Expand the border of the 2D matrix by the specified size on all sides.
    * The expanded border is filled with the mirror image of interior data.
    */
   static  matrix<> border_mirror_2D(const  matrix<>&, unsigned long);
   static cmatrix<> border_mirror_2D(const cmatrix<>&, unsigned long);
   
   /*
    * Expand the border of the 2D matrix by the specified sizes in the x- and
    * y-dimensions.  The expanded border is filled with the mirror image of
    * interior data.
    */
   static  matrix<> border_mirror_2D(
      const  matrix<>&, unsigned long, unsigned long
   );
   static cmatrix<> border_mirror_2D(
      const cmatrix<>&, unsigned long, unsigned long
   );
   
   /*
    * Trim the specified border size off all sides of the 2D matrix.
    */
   static  matrix<> border_trim_2D(const  matrix<>&, unsigned long);
   static cmatrix<> border_trim_2D(const cmatrix<>&, unsigned long);

   /*
    * Trim the specified border dimensions off of the 2D matrix.
    */
   static  matrix<> border_trim_2D(
      const  matrix<>&, unsigned long, unsigned long
   );
   static cmatrix<> border_trim_2D(
      const cmatrix<>&, unsigned long, unsigned long
   );

   /************************************************************************
    * Matrix subdivision (2D).
    *
    * Subdivide a 2D matrix into (possibly) overlapping submatrices. 
    ************************************************************************/

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
   static void subdivision_params_2D(
      unsigned long,       /* matrix size in x-dimension */
      unsigned long,       /* matrix size in y-dimension */
      unsigned long,       /* desired number of subdivisions */
      unsigned long&,      /* returned number of subdivisions in x-dimension */
      unsigned long&,      /* returned number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   static void subdivision_params_2D(
      const matrix<>&,
      unsigned long,       /* desired number of subdivisions */
      unsigned long&,      /* returned number of subdivisions in x-dimension */
      unsigned long&,      /* returned number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   static void subdivision_params_2D(
      const cmatrix<>&,
      unsigned long,       /* desired number of subdivisions */
      unsigned long&,      /* returned number of subdivisions in x-dimension */
      unsigned long&,      /* returned number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   /*
    * Subdivide the 2D matrix into a grid of overlapping submatrices.
    * Return the collection of submatrices.
    */
   static auto_collection<  matrix<>, array_list<  matrix<> > > subdivide_2D(
      const matrix<>&,     
      unsigned long,       /* number of subdivisions in x-dimension */
      unsigned long,       /* number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   static auto_collection< cmatrix<>, array_list< cmatrix<> > > subdivide_2D(
      const cmatrix<>&,
      unsigned long,       /* number of subdivisions in x-dimension */
      unsigned long,       /* number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   /*
    * Combine submatrices resulting from matrix subdivision.
    */
   static  matrix<> combine_2D(
      const array_list< matrix<> >&,
      unsigned long,       /* number of subdivisions in x-dimension */
      unsigned long,       /* number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   static cmatrix<> combine_2D(
      const array_list< cmatrix<> >&,
      unsigned long,       /* number of subdivisions in x-dimension */
      unsigned long,       /* number of subdivisions in y-dimension */
      unsigned long = 0,   /* overlap in x-dimension */
      unsigned long = 0    /* overlap in y-dimension */
   );

   /************************************************************************
    * Matrix resampling (2D).
    *
    * Resample the matrix using bilinear interpolation.
    ************************************************************************/

   /*
    * Resample the 2D matrix by the given scale factor in both the x- and y-
    * directions.  A factor less than 1 produces a result smaller than the
    * original matrix, while a factor greater than 1 produces a result larger.
    */
   static  matrix<> resample_2D(const  matrix<>&, double);
   static cmatrix<> resample_2D(const cmatrix<>&, double);
  
   /*
    * Resample the 2D matrix by the given scale factors in the x- and y-
    * directions, respectively.
    */
   static  matrix<> resample_2D(const  matrix<>&, double, double);
   static cmatrix<> resample_2D(const cmatrix<>&, double, double);
    
   /*
    * Resample the 2D matrix on a grid of the specified size.
    */
   static  matrix<> resample_2D(const  matrix<>&, unsigned long, unsigned long);
   static cmatrix<> resample_2D(const cmatrix<>&, unsigned long, unsigned long);

   /************************************************************************
    * Matrix rotation (2D).
    *
    * Rotate the 2D matrix by the given angle (in radians) about its center 
    * point.  Values in the resulting matrix are determined using bilinear 
    * interpolation.
    *
    * The resulting matrix is large enough to contain to the entire rotated 
    * version of the original matrix.  Alternatively, rotate_2D_crop crops 
    * the result matrix to be the same size as the original matrix.  In 
    * either case, extra space is padded with zeros.
    ************************************************************************/
    
   /*
    * Rotate and pad with zeros.
    */
   static  matrix<> rotate_2D(const  matrix<>&, double);
   static cmatrix<> rotate_2D(const cmatrix<>&, double);
 
   /*
    * Rotate and pad with zeros, but crop the result to be the same size as
    * the input matrix.
    */
   static  matrix<> rotate_2D_crop(const  matrix<>&, double);
   static cmatrix<> rotate_2D_crop(const cmatrix<>&, double);

   /*
    * Rotate and pad with zeros, but crop the result to the specified size.
    */
   static  matrix<> rotate_2D_crop(
      const  matrix<>&, double, unsigned long, unsigned long
   );
   static cmatrix<> rotate_2D_crop(
      const cmatrix<>&, double, unsigned long, unsigned long 
   );
 
   /************************************************************************
    * Image pyramids.
    ************************************************************************/

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
   static auto_collection< matrix<>, array_list< matrix<> > > gaussian_pyramid(
      const matrix<>&,                 /* image */
      double = math::sqrt(M_SQRT2l)    /* sigma (must be > 1) */
   );

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
   static auto_collection< matrix<>, array_list< matrix<> > > laplacian_pyramid(
      const matrix<>&,                 /* image */
      double = math::sqrt(M_SQRT2l)    /* sigma (must be > 1) */
   );

   static auto_collection< matrix<>, array_list< matrix<> > > laplacian_pyramid(
      const array_list< matrix<> >&    /* gaussian pyramid */
   );

   /*
    * Resize all images in the pyramid to be the same size.
    * The base (largest) image is scaled by the given factor in both the
    * x- and y-dimensions and all other images are resampled to that size.
    */
   static void pyramid_resize(array_list< matrix<> >&, double = 1 /* scale */);

   /*
    * Reconstruct the original image given its laplacian pyramid.
    */
   static matrix<> image_from_laplacian_pyramid(const array_list< matrix<> >&);
  
   /************************************************************************
    * Gaussian kernels.
    *
    * The kernels are evaluated at integer coordinates in the range [-s, s]
    * (in the 1D case) or [-s_x, s_x] x [-s_y, s_y] (in the 2D case), where
    * s is the specified support.
    ************************************************************************/

   /*
    * Gaussian kernel (1D).
    *
    * Specify the standard deviation and (optionally) the support.
    * The length of the returned vector is 2*support + 1.
    * The support defaults to 3*sigma.
    * The kernel is normalized to have unit L1 norm.
    * If returning a 1st or 2nd derivative, the kernel has zero mean.
    */
   static matrix<> gaussian(
      double = 1,       /* sigma */
      unsigned int = 0, /* derivative (0, 1, or 2) */
      bool = false      /* take hilbert transform? */
   );
   
   static matrix<> gaussian(
      double,           /* sigma */
      unsigned int,     /* derivative (0, 1, or 2) */
      bool,             /* take hilbert transform? */
      unsigned long     /* support */
   );
   
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
   static matrix<> gaussian_2D(
      double = 1,       /* sigma x */
      double = 1,       /* sigma y */
      double = 0,       /* orientation */
      unsigned int = 0, /* derivative in y-direction (0, 1, or 2) */
      bool = false      /* take hilbert transform in y-direction? */
   );

   static matrix<> gaussian_2D(
      double,           /* sigma x */
      double,           /* sigma y */
      double,           /* orientation */
      unsigned int,     /* derivative in y-direction (0, 1, or 2) */
      bool,             /* take hilbert transform in y-direction? */
      unsigned long,    /* x support */
      unsigned long     /* y support */
   );
  
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
   static matrix<> gaussian_cs_2D(
      double = 1,       /* sigma x */
      double = 1,       /* sigma y */
      double = 0,       /* orientation */
      double = M_SQRT2l /* scale factor between center and surround */
   );

   static matrix<> gaussian_cs_2D(
      double,           /* sigma x */
      double,           /* sigma y */
      double,           /* orientation */
      double,           /* scale factor between center and surround */
      unsigned long,    /* x support */
      unsigned long     /* y support */
   );

   /************************************************************************
    * Filter banks.
    ************************************************************************/
    
   /*
    * Get the standard set of filter orientations.
    * 
    * The standard orientations are k*pi/n for k in [0,n) where n is the 
    * number of orientation requested.
    */
   static array<double> standard_filter_orientations(unsigned long /* n */);
   
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
   static auto_collection< matrix<>, array_list< matrix<> > > gaussian_filters(
      unsigned long,          /* number of orientations */
      double,                 /* sigma */
      unsigned int = 0,       /* derivative in y-direction (0, 1, or 2) */
      bool = false,           /* take hilbert transform in y-direction? */
      double = 3.0            /* elongation ratio */
   );
   
   static auto_collection< matrix<>, array_list< matrix<> > > gaussian_filters(
      const array<double>&,   /* orientations */
      double,                 /* sigma */
      unsigned int = 0,       /* derivative in y-direction (0, 1, or 2) */
      bool = false,           /* take hilbert transform in y-direction? */
      double = 3.0            /* elongation ratio */
   );
   
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
   static auto_collection< matrix<>, array_list< matrix<> > > oe_filters_even(
      unsigned long,          /* number of orientations to sample */
      double                  /* sigma */
   );

   static auto_collection< matrix<>, array_list< matrix<> > > oe_filters_odd(
      unsigned long,          /* number of orientations to sample */
      double                  /* sigma */
   );

   static void oe_filters(
      unsigned long,          /* number of orientations to sample */
      double,                 /* sigma */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* even filters */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* odd filters  */
   );
      
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
   static auto_collection< matrix<>, array_list< matrix<> > > texton_filters(
      unsigned long,          /* number of orientations to sample */
      double                  /* sigma */
   );

   /************************************************************************
    * Image filtering.
    ************************************************************************/
    
   /*
    * Image filtering.
    *
    * Return the convolution of the image with each filter (cropped so that 
    * the result is the same size as the original image).
    *
    * Filter responses are computed in parallel when possible.
    */
   static auto_collection< matrix<>, array_list< matrix<> > > filter(
      const matrix<>&,                                      /* image */
      const collection< matrix<> >&                         /* filter set */
   );

   /*
    * Image filtering.
    *
    * Return the square of the convolution of the image with each filter
    * (cropped so that the result is the same size as the original image).
    *
    * Filter responses are computed in parallel when possible.
    */
   static auto_collection< matrix<>, array_list< matrix<> > > filter_sq(
      const matrix<>&,                                      /* image */
      const collection< matrix<> >&                         /* filter set */
   );

   /*
    * Image filtering and rectification.
    *
    * Return the rectified components of the convolution of the image with
    * each filter (cropped so that the result is the same size as the original
    * image).
    *
    * Filter responses are computed in parallel when possible.
    */
   static void filter_rectified(
      const matrix<>&,                                      /* image */
      const collection< matrix<> >&,                        /* filter set */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* responses- */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* responses+ */
   );
    
   /*
    * Image filtering and rectification.
    *
    * Return the square of the rectified components of the convolution of the
    * image with each filter (cropped so that the result is the same size as
    * the original image).
    *
    * Filter responses are computed in parallel when possible.
    */
   static void filter_rectified_sq(
      const matrix<>&,                                      /* image */
      const collection< matrix<> >&,                        /* filter set */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* responses- */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* responses+ */
   );

   /*
    * Stack the output of multiple filters into column vectors.
    */
   static matrix< matrix<> > vectorize_filter_responses(
      const collection< matrix<> >&                         /* responses */
   );

   /************************************************************************
    * Image quantization.
    ************************************************************************/
    
   /*
    * Vector quantize filter responses.
    * Return both the cluster assignments and cluster centroids.
    */
   static matrix<unsigned long> cluster_filter_responses(
      const collection< matrix<> >&,                        /* responses */
      const centroid_clusterer< matrix<> >&,                /* clusterer */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* centroids */
   );
   
   static matrix<unsigned long> cluster_filter_responses(
      const matrix< matrix<> >&,                            /* responses */
      const centroid_clusterer< matrix<> >&,                /* clusterer */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* centroids */
   );

   /*
    * Cluster image values.
    * Return both the cluster assignments and cluster centroids.
    */
   static matrix<unsigned long> cluster_values(
      const matrix<>&,                                      /* image */
      const centroid_clusterer<double>&,                    /* clusterer */
      auto_collection< double, array_list<double> >&        /* centroids */
   );

   /*
    * Quantize image values into uniformly spaced bins in [0,1].
    * Return the assignments and (optionally) bin centroids.
    */
   static matrix<unsigned long> quantize_values(
      const matrix<>&,                                      /* image */
      unsigned long                                         /* number of bins */
   );
   
   static matrix<unsigned long> quantize_values(
      const matrix<>&,                                      /* image */
      unsigned long,                                        /* number of bins */
      auto_collection< double, array_list<double> >&        /* centroids */
   );
   
   /*
    * Construct a quantized image by looking up the value of the centroid to
    * which each element is assigned.
    */
   static matrix<> quantized_image(
      const matrix<unsigned long>&,                         /* assignments */
      const array_list<double>&                             /* centroids */
   ); 

   /************************************************************************
    * Edge detection.
    ************************************************************************/

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
   static void edge_oe(
      const matrix<>&,                 /* image */
      const matrix<>&,                 /* even-symmetric filter */
      const matrix<>&,                 /* odd-symmetric filter */
      auto_ptr< matrix<> >&,           /* returned oe_e */
      auto_ptr< matrix<> >&            /* returned oe_o */
   );
    
   /*
    * Compute oriented edge energy (at multiple orientations/scales).
    *
    * Different orientations/scales are processed in parallel when possible.
    */
   static void edge_oe(
      const matrix<>&,                 /* image */
      const collection< matrix<> >&,   /* even-symmetric filter set */
      const collection< matrix<> >&,   /* odd-symmetric filter set */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* returned oe_e */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* returned oe_o */
   );
 
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
   static void edge_oe_rectified(
      const matrix<>&,                 /* image */
      const matrix<>&,                 /* odd-symmetric filter */
      auto_ptr< matrix<> >&,           /* returned oe- */
      auto_ptr< matrix<> >&            /* returned oe+ */
   );

   /*
    * Compute rectified oriented edge energy (at multiple orientations/scales).
    *
    * Different orientations/scales are processed in parallel when possible.
    */
   static void edge_oe_rectified(
      const matrix<>&,                 /* image */
      const collection< matrix<> >&,   /* odd-symmetric filter set */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* returned oe- */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* returned oe+ */
   );
 
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
   static void combine_edge_oe(
      const collection< matrix<> >&,                        /* edge oe_e */
      const collection< matrix<> >&,                        /* edge oe_o */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* edge strength */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* edge phase */
   );

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
   static void combine_edge_oe_rectified(
      const collection< matrix<> >&,                        /* edge oe- */
      const collection< matrix<> >&,                        /* edge oe+ */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* edge strength */
      auto_collection< matrix<>, array_list< matrix<> > >&  /* edge polarity */
   );
   
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
   static void combine_edge_oe(
      const collection< matrix<> >&,   /* edge oe_e */
      const collection< matrix<> >&,   /* edge oe_o */
      auto_ptr< matrix<> >&,           /* returned edge strength */
      auto_ptr< matrix<> >&,           /* returned edge phase */
      auto_ptr< matrix<> >&            /* returned edge orientation */
   );

   static void combine_edge_oe(
      const collection< matrix<> >&,   /* edge oe_e */
      const collection< matrix<> >&,   /* edge oe_o */
      const array<double>&,            /* orientations */
      auto_ptr< matrix<> >&,           /* returned edge strength */
      auto_ptr< matrix<> >&,           /* returned edge phase */
      auto_ptr< matrix<> >&            /* returned edge orientation */
   );
  
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
   static void combine_edge_oe_rectified(
      const collection< matrix<> >&,   /* edge oe- */
      const collection< matrix<> >&,   /* edge oe+ */
      auto_ptr< matrix<> >&,           /* returned edge strength */
      auto_ptr< matrix<> >&,           /* returned edge polarity */
      auto_ptr< matrix<> >&            /* returned edge orientation */
   );
  
   static void combine_edge_oe_rectified(
      const collection< matrix<> >&,   /* edge oe- */
      const collection< matrix<> >&,   /* edge oe+ */
      const array<double>&,            /* orientations */
      auto_ptr< matrix<> >&,           /* returned edge strength */
      auto_ptr< matrix<> >&,           /* returned edge polarity */
      auto_ptr< matrix<> >&            /* returned edge orientation */
   );
   
   /************************************************************************
    * Textons.
    ************************************************************************/

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
   static matrix<unsigned long> textons(
      const matrix<>&,                                      /* image */
      const collection< matrix<> >&,                        /* texton filters */
      auto_collection< matrix<>, array_list< matrix<> > >&, /* textons */
      unsigned long = 16,                                   /* # textons (K) */
      unsigned long = 10,                                   /* # iterations */
      double        = 0.10                                  /* subsampling */
   );
   
   /*
    * Compute textons as above, however specify a custom clustering routine.
    */
   static matrix<unsigned long> textons(
      const matrix<>&,                                      /* image */
      const collection< matrix<> >&,                        /* texton filters */
      const centroid_clusterer< matrix<> >&,                /* clusterer */ 
      auto_collection< matrix<>, array_list< matrix<> > >&  /* textons */
   );

   /************************************************************************
    * Difference of histograms (2D).
    ************************************************************************/

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
   static matrix<> hist_diff_2D(
      const matrix<unsigned long>&,    /* label matrix */
      unsigned long,                   /* inner radius */
      unsigned long,                   /* outer radius */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );
   
   static matrix<> hist_diff_2D(
      const matrix< matrix<> >&,       /* histogram matrix */
      unsigned long,                   /* inner radius */
      unsigned long,                   /* outer radius */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );

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
   static matrix<> hist_diff_2D(
      const matrix<unsigned long>&,    /* label matrix */
      const matrix<bool>&,             /* region mask */
      const matrix<>&,                 /* weight matrix */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );
   
   static matrix<> hist_diff_2D(
      const matrix< matrix<> >&,       /* histogram matrix */
      const matrix<bool>&,             /* region mask */
      const matrix<>&,                 /* weight matrix */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );
   
   /************************************************************************
    * Oriented gradient of histograms (2D).
    ************************************************************************/

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
   static auto_collection< matrix<>, array_list< matrix<> > > hist_gradient_2D(
      const matrix<unsigned long>&,    /* label matrix */
      unsigned long,                   /* radius */
      unsigned long,                   /* number of orientations to sample */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );
   
   static auto_collection< matrix<>, array_list< matrix<> > > hist_gradient_2D(
      const matrix< matrix<> >&,       /* histogram matrix */
      unsigned long,                   /* radius */
      unsigned long,                   /* number of orientations to sample */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );

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
   static auto_collection< matrix<>, array_list< matrix<> > > hist_gradient_2D(
      const matrix<unsigned long>&,    /* label matrix */
      const matrix<>&,                 /* weight matrix */
      unsigned long,                   /* number of orientations to sample */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );

   static auto_collection< matrix<>, array_list< matrix<> > > hist_gradient_2D(
      const matrix< matrix<> >&,       /* histogram matrix */
      const matrix<>&,                 /* weight matrix */
      unsigned long,                   /* number of orientations to sample */
      const matrix<>& = matrix<>(),    /* smoothing kernel */
      const distanceable_functor<matrix<>,double>& =
         matrix_distance_functors<>::X2_distance() /* histogram distance */
   );

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
   static void combine_hist_gradients(
      const collection< matrix<> >&,   /* oriented gradients */
      auto_ptr< matrix<> >&,           /* returned gradient strength */
      auto_ptr< matrix<> >&            /* returned gradient orientation */
   );

   static void combine_hist_gradients(
      const collection< matrix<> >&,   /* oriented gradients */
      const array<double>&,            /* orientations */
      auto_ptr< matrix<> >&,           /* returned gradient strength */
      auto_ptr< matrix<> >&            /* returned gradient orientation */
   );
   
   /************************************************************************
    * Smoothing.
    ************************************************************************/

   /*
    * Compute "needleness" at each point from derivatives.
    * 
    * The "needleness" for a signal f(x) is defined as 
    *
    *             f(x) * (-f''(x)) / (abs(f'(x)) + epsilon)
    *
    * This function is useful for smoothing and merging double peaks.
    */
   static matrix<> needleness(
      const matrix<>&,                 /* f   */
      const matrix<>&,                 /* f'  */
      const matrix<>&,                 /* f'' */
      double                           /* epsilon */
   );

   /*
    * Smooth and sharpen 2D oriented gradients using gaussian derivative 
    * filters and the "needleness" function.
    *
    * Smoothing is done along the gradient directions, which are assumed to be
    * k*pi/n for k in [0,n) where n is the number of orientated gradients.
    */
   static auto_collection< matrix<>, array_list< matrix<> > > smooth_grad_2D(
      const collection< matrix<> >&,   /* oriented gradients */
      double,                          /* sigma of smoothing filter */
      double                           /* epsilon for "needleness" */
   );

   /*
    * Smooth and sharpen 2D oriented gradients using gaussian derivative 
    * filters and the "needleness" function.
    *
    * Specify the orientations along which to smooth.
    */
   static auto_collection< matrix<>, array_list< matrix<> > > smooth_grad_2D(
      const collection< matrix<> >&,   /* oriented gradients */
      const array<double>&,            /* orientations */
      double,                          /* sigma of smoothing filter */
      double                           /* epsilon for "needleness" */
   );
    
   /************************************************************************
    * Non-max suppression.
    ************************************************************************/

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
   static matrix<> nonmax(const matrix<>&);

   static matrix<> nonmax(const matrix<>&, const matrix<bool>& /* mask */);
   
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
   static matrix<> nonmax(
      const matrix<>&,
      unsigned long,    /* dimension */
      bool = false      /* allow boundary? */
   );

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
   static matrix<> nonmax_oriented_2D(
      const matrix<>&,
      double,           /* orientation */
      bool = false      /* allow boundary? */
   );

   static matrix<> nonmax_oriented_2D(
      const matrix<>&,
      const matrix<>&,  /* orientation */
      double = M_PI_2l, /* orientation tolerance (must be >= 0) */
      bool = false      /* allow boundary? */
   );

   /************************************************************************
    * Skeletonization (2D).
    ************************************************************************/

   /*
    * Skeletonize the given 2D matrix, preserving its topology.
    *
    * Nonzero elements are treated as part of an area to be skeletonized and
    * zero elements are treated as background.
    *
    * Skeletonization proceeds by repeatedly removing (setting to zero) the
    * smallest nonzero element whos deletion best preserves the topology.
    */
   static matrix<> skeletonize_2D(const matrix<>&);

   /************************************************************************
    * Connected components.
    ************************************************************************/

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
   static matrix<unsigned long> connected_components(
      const matrix<>&
   );
   
   static matrix<unsigned long> connected_components(
      const matrix<unsigned long>&
   );

   static matrix<unsigned long> connected_components(
      const matrix<>&,
      const matrix<bool>&  /* mask */
   );

   static matrix<unsigned long> connected_components(
      const matrix<unsigned long>&,
      const matrix<bool>&  /* mask */
   );

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
   static matrix<unsigned long> connected_components(
      const matrix<>&,
      unsigned long        /* dimension */
   );
  
   static matrix<unsigned long> connected_components(
      const matrix<unsigned long>&,
      unsigned long        /* dimension */
   );
   
   /*
    * Return the number of connected components (excluding the zero component)
    * in a labeled component matrix generated by one of the above functions.
    *
    * Since components are numbered in order, this is the same as the maximum
    * value in the matrix.
    */
   static unsigned long count_connected_components(
      const matrix<unsigned long>&  /* labeled components */
   );
   
   /************************************************************************
    * Contour extraction (2D).
    ************************************************************************/
  
   /*
    * Declare contour vertex and edge classes.
    */
   class contour_vertex;
   class contour_edge;
   
   /*
    * Contour vertex.
    *
    * A contour vertex is a location in the pixel grid along with the lists of
    * edges that start or end at that location.
    */
   class contour_vertex {
   public:
      /*********************************************************************
       * Vertex data.
       *********************************************************************/

      /* vertex identity */
      unsigned long id;                      /* id of vertex */
      bool is_subdivision;                   /* is contour subdivision point? */
      /* vertex coordinates */
      unsigned long x;                       /* x-coordinate */
      unsigned long y;                       /* y-coordinate */
      /* vertex connectivity */
      array_list<contour_edge> edges_start;  /* edges that start at vertex */
      array_list<contour_edge> edges_end;    /* edges that end at vertex */

      /*********************************************************************
       * Vertex methods.
       *********************************************************************/

      point_2D point() const; /* return the point representing the vertex */
   };
   
   /*
    * Contour edge.
    *
    * A contour edge is a thin connected path between two vertices in the pixel
    * grid.  Each contour edge stores the coordinates of the pixels along this
    * path (excluding its vertices), the ids of its endpoints, and its own id
    * and contour id.
    *
    * Each contour equivalence class is a longest possible sequence of connected
    * contours that were part of the same original uninterrupted contour.
    */
   class contour_edge {
   public:
      /*********************************************************************
       * Edge data. 
       *********************************************************************/

      /* edge identity */
      unsigned long id;                /* id of edge */
      unsigned long contour_equiv_id;  /* id of edge's equivalence class */
      bool is_completion;              /* is the edge a completion edge? */
      /* edge coordinates */
      array<unsigned long> x_coords;   /* x-coordinates of interior points */
      array<unsigned long> y_coords;   /* y-coordinates of interior points */
      /* edge connectivity */
      contour_vertex* vertex_start;    /* vertex starting the edge */
      contour_vertex* vertex_end;      /* vertex ending the edge */
      unsigned long vertex_start_enum; /* position in start vertex edge list */
      unsigned long vertex_end_enum;   /* position in end vertex edge list */

      /*********************************************************************
       * Edge methods.
       *********************************************************************/
      
      unsigned long size() const;      /* return size (# of interior points) */
      double length() const;           /* return distance between endpoints */
      point_2D point(unsigned long) const; /* return point at given position */
   };

   /*
    * Declare region set class.
    */
   class region_set;

   /*
    * Contour set.
    *
    * A contour set is a collection of contour vertices and contour edges.
    * The edges in the set intersect only at vertices and all intersection
    * vertices are contained in the set.
    *
    * In addition, a contour set maintains a map assigning each pixel to
    * the vertex or (non-completion) edge passing through it (if any).
    *
    * Contours vertices and edges are labeled from 0 to V-1 and from 0 to E-1,
    * where V and E are the number of vertices and edges, respectfully.
    */
   class contour_set {
   public:
      /*
       * Friend classes.
       */
      friend class region_set;

      /*********************************************************************
       * Constructors and destructor.
       *********************************************************************/
      
      /*
       * Extract contours from labeled connected components.
       *
       * Break each connected component into a set of thin contours that may
       * intersect only at returned contour vertices.
       */
      explicit contour_set(
         const matrix<unsigned long>& /* labeled components */
      );
      
      /*
       * Copy constructor.
       */
      explicit contour_set(const contour_set&);
      
      /*
       * Destructor.
       */
      virtual ~contour_set();

      /*********************************************************************
       * Contour vertex and edge access.
       *********************************************************************/

      /*
       * Get image dimensions.
       */
      unsigned long size_x() const;
      unsigned long size_y() const;

      /*
       * Get number of vertices in the contour set.
       */
      unsigned long vertices_size() const;
      
      /*
       * Get number of edges in the contour set.
       */
      unsigned long edges_size() const;

      /*
       * Get the number of edge equivalence classes.
       * Each equivalence class represents a contour before subdivision.
       */
      unsigned long edges_equiv_size() const;

      /*
       * Return the specified vertex (by reference).
       */
      const contour_vertex& vertex(unsigned long /* vertex id */) const;

      /*
       * Return the specified edge (by reference).
       */
      const contour_edge& edge(unsigned long /* edge id */) const;
     
      /*
       * Return the specified edge equivalence class (by reference).
       *
       * The edges are listed in the order encountered on a traversal of the 
       * original contour before subdivision.
       */
      const list<contour_edge>& edge_equiv(unsigned long /* equiv id */) const;
      
      /*
       * Check if the specified pixel is a vertex.
       */
      bool is_vertex(unsigned long /* absolute index */) const;
      bool is_vertex(unsigned long /* x */, unsigned long /* y */) const;
      
      /*
       * Check if the specified pixel lies on a (non-completion) edge interior.
       */
      bool is_edge(unsigned long /* absolute index */) const;
      bool is_edge(unsigned long /* x */, unsigned long /* y */) const;

      /*
       * Return the id of the vertex at the given pixel.
       * Throw an exception (ex_not_found) if no such vertex exists.
       */
      unsigned long vertex_id(
         unsigned long /* absolute index */) const;
      unsigned long vertex_id(
         unsigned long /* x */, unsigned long /* y */) const;
      
      /*
       * Return the id of the (non-completion) edge whos interior contains the
       * given pixel.  Throw an exception (ex_not_found) if no such edge exists.
       */
      unsigned long edge_id(
         unsigned long /* absolute index */) const;
      unsigned long edge_id(
         unsigned long /* x */, unsigned long /* y */) const;
      
      /*********************************************************************
       * Contour subdivision.
       *
       * The following subdivision methods split contours in order to enforce
       * certain properties of the straight line segment approximation to the
       * contours.
       *
       * The "local" property concerns approximation quality.  Contours
       * should locally approximate line segments well.
       *
       * The "global" property concerns topological consistency.  If the 
       * contours were replaced by line segments, all intersection points 
       * should be preserved (and no new intersections introduced).
       *
       * These two properties are independent in the sense that it is possible
       * to satisfy one without satisfying the other.
       *********************************************************************/
     
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
      void subdivide_local(
         bool   = true,       /* enforce angle constraint? */
         double = M_PI_2l,    /* maximum deviation angle (in radians) */
         bool   = true,       /* enforce distance constraint? */
         double = 0.05,       /* maximum fractional deviation */
         double = 2.0         /* tolerance distance (in pixels) */
      );

      /*
       * Recursively subdivide the contours until their straight line segment
       * approximation respects the true contour topology.
       *
       * Enforce the global criteria that the line segments between contour
       * endpoints intersect only at contour vertices.
       */
      void subdivide_global();

      /*********************************************************************
       * Contour completion.
       *
       * Completion edges are contours that are possible candidates for true
       * boundaries, but were not detected through low-level image processing.
       * Such completion edges often fill gaps between existing contours.
       *
       * Edges generated through contour completion are tagged by setting
       * their is_completion flag.
       *
       * Completion edges connect existing contour vertices, however their
       * interiors may overlap with the interiors of previously existing edges.
       * Note that the map from pixels to edges containing them does not
       * include completion edges.
       *********************************************************************/
 
      /*
       * Add the four bounding corners of the image to the set of contour
       * vertices if they neither currently exist in the vertex set nor lie
       * on an existing contour.
       */
      void add_bounding_vertices();

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
      void completion_cdt();

      void completion_cdt(
         auto_collection<
            point_2D, array_list<point_2D> >&,  /* returned CDT vertices */
         auto_ptr<triangulation>&               /* returned CDT */
      );

   protected:
      /*********************************************************************
       * Contour set data.
       *********************************************************************/
       
      /*
       * Pixel assignments.
       */
      matrix<bool> _is_vertex;   /* is pixel a vertex? */
      matrix<bool> _is_edge;     /* is pixel on a (non-completion) edge? */
      matrix<unsigned long> _assignments; /* map of pixel -> vertex/edge id */

      /*
       * Contour vertices and edges.
       */
      auto_collection< contour_vertex, array_list<contour_vertex> > _vertices;
      auto_collection< contour_edge,   array_list<contour_edge>   > _edges;

      /*
       * Contour edge equivalence classes.
       */
      auto_collection< 
         list<contour_edge>, array_list< list<contour_edge> > > _edges_equiv;
         
      /*********************************************************************
       * Contour edge operations.
       *********************************************************************/
       
      /*
       * Given connected component labels, the starting vertex of an edge, 
       * and the next pixel along the edge, trace the entire contour and 
       * add it to the set.
       */
      void trace_edge(
         const matrix<unsigned long>&, /* labeled components */
         unsigned long,                /* starting vertex id */
         unsigned long,                /* x-coordinate of next pixel */
         unsigned long                 /* y-coordinate of next pixel */
      );
      
      /*
       * Split the contour edge with the given id at the specified position.
       */
      void split_edge(unsigned long /* id */, unsigned long /* position */);

      /*
       * Sort the list of edges in each edge equivalence class so they appear
       * in the order encountered on a traversal of the original contour edge.
       */
      void sort_edges_equiv();
   };

   /************************************************************************
    * Region extraction (2D).
    ************************************************************************/

   /*
    * Region set.
    *
    * A region set is a collection of regions bordering contour edges in a
    * corresponding contour set.  Each region is represented by the array of
    * indices (linear matrix offsets) of the pixels contained within it.
    *
    * The set stores separate regions for the left and right sides of each
    * contour.  Traversal of a contour from its start vertex to its end vertex
    * defines the viewpoint from which left/right is determined.
    *
    * During region extraction, we regard each contour edge as an approximately
    * linear element.  Hence, contours should be appropriately subdivided prior
    * to computing regions.  We compute both regions bordering the subdivided
    * edges and regions bordering the original contours before subdivision.
    */
   class region_set {
   public:
      /*********************************************************************
       * Constructors and destructor.
       *********************************************************************/

      /*
       * Extract the regions lying on the left and right sides of each of the 
       * given contours.
       * 
       * The width of both the left and right regions is determined as the
       * specified fraction of the length of the corresponding contour edge,
       * with a minimum width used for short edges.  In this manner, region
       * size scales with contour length.
       */
      explicit region_set(
         const contour_set&,     /* contour set */
         double,                 /* region width (edge fraction) */
         double = 2.0            /* minimum width (in pixels) */
      );

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
      explicit region_set(
         const contour_set&,     /* contour set */
         double,                 /* region width (edge fraction) */
         double,                 /* minimum width (in pixels) */
         double                  /* relative scale sensitivity */
      );

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
      explicit region_set(
         const contour_set&,     /* contour set */
         double,                 /* region width (edge fraction) */
         double,                 /* minimum width (in pixels) */
         const array<double>&,   /* edge visibility */
         const matrix<>&         /* pixel visibility */
      );

      /*
       * Copy constructor.
       */
      explicit region_set(const region_set&);

      /*
       * Destructor.
       */
      virtual ~region_set();

      /*********************************************************************
       * Region access.
       *********************************************************************/

      /*
       * Get number of vertices in the corresponding contour set.
       */
      unsigned long vertices_size() const;
      
      /*
       * Get number of edges in the region set.
       */
      unsigned long edges_size() const;

      /*
       * Get the number of edge equivalence classes.
       * Each equivalence class represents a contour before subdivision.
       */
      unsigned long edges_equiv_size() const;

      /*
       * Return the region (by reference) on the left/right side of the vertex. 
       * Only subdivision vertices are associated with nonempty regions.
       */
      const array<unsigned long>&
         vertex_region_left(unsigned long /* vertex id */) const;
      const array<unsigned long>&
         vertex_region_right(unsigned long /* vertex id */) const;
      
      /*
       * Return the region (by reference) on the left/right side of the edge.
       */
      const array<unsigned long>&
         edge_region_left(unsigned long /* edge id */) const;
      const array<unsigned long>&
         edge_region_right(unsigned long /* edge id */) const;

      /*
       * Return the region (by reference) on the left/right side of the given
       * edge equivalence class (original contour before subdivision).
       */
      const array<unsigned long>&
         edge_equiv_region_left(unsigned long /* equiv id */) const;
      const array<unsigned long>&
         edge_equiv_region_right(unsigned long /* equiv id */) const;

   protected:
      /*********************************************************************
       * Region set data.
       *********************************************************************/

      array< array<unsigned long> > _vertex_regions_left;
      array< array<unsigned long> > _vertex_regions_right;

      array< array<unsigned long> > _edge_regions_left;
      array< array<unsigned long> > _edge_regions_right;

      array< array<unsigned long> > _edge_equiv_regions_left;
      array< array<unsigned long> > _edge_equiv_regions_right;

      /*********************************************************************
       * Region extraction.
       *********************************************************************/
  
      /*
       * Extract regions using a custom intervening contour cue.
       */
      void extract_regions(
         const contour_set&,     /* contour set */
         double,                 /* region width (edge fraction) */
         double,                 /* minimum width (in pixels) */
         const array<double>&,   /* edge visibility */
         const matrix<>&         /* pixel visibility */
      );
   };

   static matrix<> line_inds(double, double, double, double, double);
};

} /* namespace libraries */
} /* namespace math */

#endif
