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


#ifndef __smatrix_h__
#define __smatrix_h__

#include <stdio.h>
#include "array.hh"

class SMatrix
{
  public:
    SMatrix(int n, int* nz, int** col, double** values);
    ~SMatrix();
    SMatrix(FILE* fp);
    void dump(FILE* fp);

    double computeNCut(const double* rowSum, const Util::Array1D<int> membership, const int nsegs) const;

    void symmetrize();

    int getNNZ() const;
    double* getRowSum() const;

    void normalizedLaplacian(const double* rowSum);        //in place, converts into norm laplacian 
    void undoNormalizedLaplacian(const double* rowSum);    //in place, converts back into original matrix

    int n;
    int** col;
    int* nz;
    double** values;

    void mvmul(const double* a, double* b) const;
    void mvmul(const double* a1, const double* a2, double* b1, double* b2) const;
    void mvmul(const double* a1, const double* a2, const double* a3, const double* a4, 
               double* b1, double* b2, double* b3, double* b4) const;

};

#endif

