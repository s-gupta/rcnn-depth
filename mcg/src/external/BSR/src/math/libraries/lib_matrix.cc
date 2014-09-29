/*
 * Matrix library.
 */
#include "lang/pointers/auto_ptr.hh"
#include "math/complex.hh"
#include "math/libraries/lib_matrix.hh"
#include "math/math.hh"
#include "math/matrices/cmatrix.hh"
#include "math/matrices/matrix.hh"

/* FIXME */
#include "lang/exceptions/ex_not_implemented.hh"

namespace math {
namespace libraries {
/*
 * Imports.
 */
using lang::pointers::auto_ptr;
using math::complex;
using math::matrices::matrix;
using math::matrices::cmatrix;

/* FIXME */
using lang::exceptions::ex_not_implemented;

/*
 * Eigenvalues and eigenvectors.
 */
void lib_matrix::eig(
   const matrix<>& m,
   auto_ptr< matrix<> >& evals,
   auto_ptr< matrix<> >& evecs)
{
   /* FIXME */
   throw ex_not_implemented();
}

void lib_matrix::eig(
   const cmatrix<>& m,
   auto_ptr< cmatrix<> >& evals,
   auto_ptr< cmatrix<> >& evecs)
{
   /* FIXME */
   throw ex_not_implemented();
}

/*
 * Singular value decomposition.
 * Return a diagonal matrix S of the same dimensions as the input matrix X
 * and unitary matrices U and V, such that X = U*S*V'.  The diagonal 
 * entries of S are nonnegative and in decreasing order.
 */
void lib_matrix::svd(
   const matrix<>& m,
   auto_ptr< matrix<> >& u,
   auto_ptr< matrix<> >& s,
   auto_ptr< matrix<> >& v)
{
   /* FIXME */
   throw ex_not_implemented();
}

void lib_matrix::svd(
   const cmatrix<>& m,
   auto_ptr< cmatrix<> >& u,
   auto_ptr< cmatrix<> >& s,
   auto_ptr< cmatrix<> >& v)
{
   /* FIXME */
   throw ex_not_implemented();
}

} /* namespace libraries */
} /* namespace math */
