// ------------------------------------------------------------------------ 
//  Copyright (C)
//  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain
//  University of California Berkeley (UCB) - USA
// 
//  Jordi Pont-Tuset <jordi.pont@upc.edu>
//  Pablo Arbelaez <arbelaez@berkeley.edu>
//  June 2014
// ------------------------------------------------------------------------ 
// This file is part of the MCG package presented in:
//    Arbelaez P, Pont-Tuset J, Barron J, Marques F, Malik J,
//    "Multiscale Combinatorial Grouping,"
//    Computer Vision and Pattern Recognition (CVPR) 2014.
// Please consider citing the paper if you use this code.
// ------------------------------------------------------------------------
/*
 * Contours. 
 */
#include "collections/pointers/auto_collection.hh"
#include "collections/array_list.hh"
#include "io/streams/cout.hh"
#include "io/streams/iomanip.hh"
#include "io/streams/ios.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/libraries/lib_image.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"

//#include "concurrent/threads/thread.hh"

#include <time.h>
#include <mex.h>

using collections::pointers::auto_collection;
using collections::array_list;
using io::streams::cout;
using io::streams::ios;
using io::streams::iomanip::setiosflags;
using io::streams::iomanip::setw;
using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using lang::pointers::auto_ptr;
using math::libraries::lib_image;
using math::matrices::matrix;

//using concurrent::threads::thread;

/********************************** 
 * Matlab matrix conversion routines.
 **********************************/

/*
 * Get a string from an mxArray.
 */
char *mexGetString(const mxArray *arr) {
   char *string;
   int buflen;
   buflen = (mxGetM(arr) * mxGetN(arr) * sizeof(mxChar)) + 1;
   string = new char[buflen];
   mxGetString(arr, string, buflen);
   return string;
}

/*
 * Create a single element Matlab double matrix.
 */
mxArray* mxCreateDoubleScalar(double value) {
   mxArray* pa = mxCreateDoubleMatrix(1, 1, mxREAL);
   *mxGetPr(pa) = value;
   return pa;
}

/* 
 * Convert an mxArray to a matrix.
 */
matrix<> to_matrix(const mxArray *a) {
   unsigned long mrows = static_cast<unsigned long>(mxGetM(a));
   unsigned long ncols = static_cast<unsigned long>(mxGetN(a));
   double *data = mxGetPr(a);
   matrix<> m(mrows, ncols);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         m(r,c) = data[(c*mrows) + r];
      }
   }
   return m;
}
   
/*
 * Convert a 2D matrix to an mxArray.
 */
mxArray* to_mxArray(const matrix<>& m) {
   unsigned long mrows = m.size(0);
   unsigned long ncols = m.size(1);
   mxArray *a = mxCreateDoubleMatrix(
      static_cast<int>(mrows),
      static_cast<int>(ncols),
      mxREAL
   );
   double *data = mxGetPr(a);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         data[(c*mrows) + r] = m(r,c);
      }
   }
   return a;
}

/*
 * Convert an array to an mxArray.
 */
mxArray* to_mxArray(const array<double>& m) {
   unsigned long mrows = m.size();
   unsigned long ncols = 1;
   mxArray *a = mxCreateDoubleMatrix(
      static_cast<int>(mrows),
      static_cast<int>(ncols),
      mxREAL
   );
   double *data = mxGetPr(a);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         data[(c*mrows) + r] = m[r];
      }
   }
   return a;
}



/*
 * Matlab interface.
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
   try {
      /* get nonmax-suppressed image */
      matrix<> im = to_matrix(prhs[0]);
      /* compute discrete contours */
      //cout << setiosflags(ios::left);
      //cout << setw(40) << "compute contours ";
      //cout.flush();
      clock_t start_time = clock();
      clock_t time = start_time;
      matrix<> im_skel = im; //lib_image::skeletonize_2D(im);
      matrix<unsigned long> labels = lib_image::connected_components(im_skel);
      lib_image::contour_set contours(labels);
      contours.subdivide_local();
      //contours.subdivide_global();
      //contours.completion_cdt();
      //cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      //cout.flush();
      /* extract contour sides */
      /*
      //cout << setw(40) << "extracting regions ";
      //cout.flush();
      time = clock();
      lib_image::region_set regions(contours, 0.2, 2.0, 0.5);
      //cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      //cout.flush();
      */
      /* get sizes of regions */
      //cout << setw(40) << "returning results ";
      //cout.flush();
      time = clock();
      /* extract vertex/edge map */
      matrix<bool> is_v(labels.dimensions());
      matrix<bool> is_e(labels.dimensions());
      matrix<unsigned long> assign(labels.dimensions());
      for (unsigned long n = 0; n < labels.size(); n++) {
         is_v[n] = contours.is_vertex(n);
         is_e[n] = contours.is_edge(n);
         if (is_v[n])
            assign[n] = contours.vertex_id(n);
         else if (is_e[n])
            assign[n] = contours.edge_id(n);
      }
      /* extract vertex coordinates */
      unsigned long n_vertices = contours.vertices_size();
      matrix<unsigned long> vertex_coords(n_vertices, 2);
      for (unsigned long n = 0; n < n_vertices; n++) {
         vertex_coords(n,0) = contours.vertex(n).x;
         vertex_coords(n,1) = contours.vertex(n).y;
      }
      /* extract edge endpoints */
      unsigned long n_edges = contours.edges_size();
      matrix<unsigned long> edge_endpoints(n_edges, 2);
      for (unsigned long n = 0; n < n_edges; n++) {
         edge_endpoints(n,0) = contours.edge(n).vertex_start->id;
         edge_endpoints(n,1) = contours.edge(n).vertex_end->id;
      }
      /* extract edge coordinates */
      mxArray* mx_e_x_coords = mxCreateCellMatrix(static_cast<int>(n_edges), 1);
      mxArray* mx_e_y_coords = mxCreateCellMatrix(static_cast<int>(n_edges), 1);
      for (unsigned long n = 0; n < n_edges; n++) {
         array<double> e_x(contours.edge(n).x_coords);
         array<double> e_y(contours.edge(n).y_coords);
         mxSetCell(mx_e_x_coords, static_cast<int>(n), to_mxArray(e_x));
         mxSetCell(mx_e_y_coords, static_cast<int>(n), to_mxArray(e_y));
      }
      /* extract completion flags */
      matrix<bool> is_compl(n_edges, 1);
      for (unsigned long n = 0; n < n_edges; n++)
         is_compl[n] = contours.edge(n).is_completion;
      /* create vertex region cell arrays */
      mxArray* mx_v_region_left  = mxCreateCellMatrix(
         static_cast<int>(n_vertices), 1);
      mxArray* mx_v_region_right = mxCreateCellMatrix(
         static_cast<int>(n_vertices), 1);
      /*
      for (unsigned long n = 0; n < n_vertices; n++) {
         array<double> r_left(regions.vertex_region_left(n));
         array<double> r_right(regions.vertex_region_right(n));
         mxSetCell(mx_v_region_left,  static_cast<int>(n), to_mxArray(r_left));
         mxSetCell(mx_v_region_right, static_cast<int>(n), to_mxArray(r_right));
      }
      */
      /* create edge region cell arrays */
      matrix<unsigned long> edge_equiv_ids(n_edges, 1); 
      mxArray* mx_e_region_left  = mxCreateCellMatrix(
         static_cast<int>(n_edges), 1);
      mxArray* mx_e_region_right = mxCreateCellMatrix(
         static_cast<int>(n_edges), 1);
      /*
      for (unsigned long n = 0; n < n_edges; n++) {
         edge_equiv_ids[n] = contours.edge(n).contour_equiv_id;
         array<double> r_left(regions.edge_region_left(n));
         array<double> r_right(regions.edge_region_right(n));
         mxSetCell(mx_e_region_left,  static_cast<int>(n), to_mxArray(r_left));
         mxSetCell(mx_e_region_right, static_cast<int>(n), to_mxArray(r_right));
      }
      */
      /* create contour region cell arrays */
      unsigned long n_edges_equiv = contours.edges_equiv_size();
      mxArray* mx_c_region_left  = mxCreateCellMatrix(
         static_cast<int>(n_edges_equiv), 1);
      mxArray* mx_c_region_right = mxCreateCellMatrix(
         static_cast<int>(n_edges_equiv), 1);
      /*
      for (unsigned long n = 0; n < n_edges_equiv; n++) {
         array<double> r_left(regions.edge_equiv_region_left(n));
         array<double> r_right(regions.edge_equiv_region_right(n));
         mxSetCell(mx_c_region_left,  static_cast<int>(n), to_mxArray(r_left));
         mxSetCell(mx_c_region_right, static_cast<int>(n), to_mxArray(r_right));
      }
      */
      /* return */
      if (nlhs > 0)
         plhs[0] = to_mxArray(im_skel);
      if (nlhs > 1)
         plhs[1] = to_mxArray(matrix<>(labels));
      if (nlhs > 2)
         plhs[2] = to_mxArray(matrix<>(is_v));
      if (nlhs > 3)
         plhs[3] = to_mxArray(matrix<>(is_e));
      if (nlhs > 4)
         plhs[4] = to_mxArray(matrix<>(assign));
      if (nlhs > 5)
         plhs[5] = to_mxArray(matrix<>(vertex_coords));
      if (nlhs > 6)
         plhs[6] = to_mxArray(matrix<>(edge_endpoints));
      if (nlhs > 7)
         plhs[7] = mx_v_region_left;
      if (nlhs > 8)
         plhs[8] = mx_v_region_right; 
      if (nlhs > 9)
         plhs[9] = mx_e_region_left;
      if (nlhs > 10)
         plhs[10] = mx_e_region_right; 
      if (nlhs > 11)
         plhs[11] = mx_c_region_left;
      if (nlhs > 12)
         plhs[12] = mx_c_region_right; 
      if (nlhs > 13)
         plhs[13] = to_mxArray(matrix<>(edge_equiv_ids)); 
      if (nlhs > 14)
         plhs[14] = to_mxArray(matrix<>(is_compl)); 
      if (nlhs > 15)
         plhs[15] = mx_e_x_coords;
      if (nlhs > 16)
         plhs[16] = mx_e_y_coords;
      //cout << "[" << (double(clock() - time)/CLOCKS_PER_SEC) << " sec]\n";
      //cout.flush();
   } catch (ex_index_out_of_bounds& e) {
      cout << "index: " << e.index() << "\n";
      cout << e << "\n";
   } catch (exception& e) {
      cout << e << "\n";
   }
}
