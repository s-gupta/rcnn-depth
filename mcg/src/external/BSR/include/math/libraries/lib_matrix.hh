/*
 * Matrix library.
 */
#ifndef MATH__LIBRARIES__LIB_MATRIX_HH
#define MATH__LIBRARIES__LIB_MATRIX_HH

#include "lang/pointers/auto_ptr.hh"
#include "math/matrices/cmatrix.hh"
#include "math/matrices/matrix.hh"

namespace math {
namespace libraries {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;
using math::matrices::cmatrix;
using math::matrices::matrix;

/*
 * Matrix functions.
 */
class lib_matrix {
public:
   /* FIXME: add more matrix functions */
   
   /************************************************************************
    * Matrix decomposition.
    ************************************************************************/

   /*
    * Eigenvalues and eigenvectors.
    */
   static void eig(
      const matrix<>&,
      auto_ptr< matrix<> >&,  /* returned eigenvalues  */
      auto_ptr< matrix<> >&   /* returned eigenvectors */
   );
   
   static void eig(
      const cmatrix<>&,
      auto_ptr< cmatrix<> >&, /* returned eigenvalues  */
      auto_ptr< cmatrix<> >&  /* returned eigenvectors */
   );
   
   /*
    * Singular value decomposition.
    * Return a diagonal matrix S of the same dimensions as the input matrix X
    * and unitary matrices U and V, such that X = U*S*V'.  The diagonal 
    * entries of S are nonnegative and in decreasing order.
    */
   static void svd(
      const matrix<>&,
      auto_ptr< matrix<> >&,  /* returned U */
      auto_ptr< matrix<> >&,  /* returned S */
      auto_ptr< matrix<> >&   /* returned V */
   );
   
   static void svd(
      const cmatrix<>&,
      auto_ptr< cmatrix<> >&, /* returned U */
      auto_ptr< cmatrix<> >&, /* returned S */
      auto_ptr< cmatrix<> >&  /* returned V */
   );
};

} /* namespace libraries */
} /* namespace math */

#endif
