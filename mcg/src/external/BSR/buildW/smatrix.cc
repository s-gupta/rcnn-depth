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


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "smatrix.hh"

SMatrix::SMatrix (int n, int* nz, int** col, double** values)
{
    this->n = n;
    this->nz = nz;
    this->col = col;
    this->values = values;
    int nnz = 0;
    for (int i = 0; i < n; i++)
      nnz += nz[i];
    //Util::Message::debug(Util::String("creating sparse matrix with %d nonzero entries",nnz));
}

SMatrix::SMatrix(FILE* fp) 
{
  size_t st;
  st = fread(&n,sizeof(int),1,fp);
  int nnz;
  st =fread(&nnz,sizeof(int),1,fp);
  nz = new int[n];
  col = new int*[n];
  values = new double*[n];

  for (int row = 0; row < n; row++) {
    st =fread(&nz[row],sizeof(int),1,fp);
    col[row] = new int[nz[row]]; 
    values[row] = new double[nz[row]]; 
    st =fread(values[row],sizeof(double),nz[row],fp);
    st =fread(col[row],sizeof(int),nz[row],fp);
  }
}

SMatrix::~SMatrix ()
{
  for (int i = 0; i < n; i++)
  {
    delete[] col[i];
    delete[] values[i];
  }
  delete col;
  delete values;
  delete nz;
}

void 
SMatrix::dump(FILE* fp) 
{
  size_t st;
  st = fwrite(&n,sizeof(int),1,fp);
  int nnz = getNNZ();
  st = fwrite(&nnz,sizeof(int),1,fp);
  for (int row = 0; row < n; row++) {
    st = fwrite(&nz[row],sizeof(int),1,fp);
    st = fwrite(values[row],sizeof(double),nz[row],fp);
    st = fwrite(col[row],sizeof(int),nz[row],fp);
  }
}

int SMatrix::getNNZ() const
{
  int nnz = 0;
  for (int i = 0; i < n; i++)
    nnz += nz[i];
  return nnz;
}


void SMatrix::symmetrize()
{
  int* tail = new int[n];  
  memset(tail,0,n*sizeof(int));
  for (int r = 0; r < n; r++) 
  {
    int offset = 0;
    while ((offset < nz[r]) && (col[r][offset] < r+1))
    {
      offset++;
    }
    for (int i = offset; i < nz[r]; i++) 
    {
      int c = col[r][i];
      assert( col[c][tail[c]] == r ); 
      double v_rc = values[r][i];
      double v_cr = values[c][tail[c]];
      values[r][i] = 0.5*(v_rc+v_cr);
      values[c][tail[c]] = 0.5*(v_rc+v_cr);
      tail[c]++;
    }
  }  
}


double*
SMatrix::getRowSum () const
{
    double* rowSum = new double[n];
    memset(rowSum,0,n*sizeof(double));
    //Util::Message::debug("computing row sum");
    for (int row = 0; row < n; row++) {
      double sum = 0;;
      for (int j = 0; j < nz[row]; j++) {
        sum += values[row][j];
      }
      rowSum[row] = sum;
    }
    return rowSum;
}

double
SMatrix::computeNCut(const double* rowSum, const Util::Array1D<int> membership, const int nsegs) const
{
  if (nsegs > 1)
  {
    double* num = new double[nsegs];
    double* den = new double[nsegs];
    memset(num,0,nsegs*sizeof(double));
    memset(den,0,nsegs*sizeof(double));
    for (int row = 0; row < n; row++)
    {
      int segi = membership(row);
      den[segi] += rowSum[row]; 
      for (int j = 0; j < nz[row]; j++)
      {
        int segj = membership(col[row][j]);
        if (segi == segj)
        {
          num[segi] += values[row][j]; 
        }
      }
    }

    double assoc = 0;
    for (int s = 0; s < nsegs; s++)
    {
      if (den[s] > 0)
      {
        assoc += (num[s] / den[s]); 
      }
    }
    delete num;
    delete den;
    return (1 - ((1/(double)nsegs)*assoc));
  }
  else
  {
    //Util::Message::debug("only 1 segment!!");
    return 0;
  }
}


void
SMatrix::normalizedLaplacian(const double* rowSum) 
{
    double* isrd = new double[n];
    for (int i = 0; i < n; i++) {
      isrd[i] = 1.0 / sqrt(rowSum[i]);
    }
    //Util::Message::debug("scaling rows");
    for (int row = 0; row < n; row++) {
      double isrdrow = isrd[row];
      for (int j = 0; j < nz[row]; j++) {
        values[row][j] = isrdrow * values[row][j] * isrd[col[row][j]];
      }
    }
    delete[] isrd;
}

void
SMatrix::undoNormalizedLaplacian(const double* rowSum) 
{
    double* isrd = new double[n];
    for (int i = 0; i < n; i++) {
      isrd[i] = sqrt(rowSum[i]);
    }
    //Util::Message::debug("scaling rows");
    for (int row = 0; row < n; row++) {
      double isrdrow = isrd[row];
      for (int j = 0; j < nz[row]; j++) {
        values[row][j] = isrdrow * values[row][j] * isrd[col[row][j]];
      }
    }
    delete[] isrd;
}

void 
SMatrix::mvmul (const double* a, double* b) const
{
    for (int row = 0; row < n; row++) {
      double bval = 0;
      for (int j = 0; j < nz[row]; j++) {
          bval += a[col[row][j]] * values[row][j];
      }
      b[row] = bval;
    }
}

void 
SMatrix::mvmul (const double* a1, const double* a2, 
			double* b1, double* b2) const
{
    for (int row = 0; row < n; row++) {
      double bval1 = 0;
      double bval2 = 0;
      for (int j = 0; j < nz[row]; j++) {
          double v = values[row][j];
          bval1 += a1[col[row][j]] * v;
          bval2 += a2[col[row][j]] * v;
      }
      b1[row] = bval1;
      b2[row] = bval2;
    }
}

void 
SMatrix::mvmul (const double* a1, const double* a2, 
			const double* a3, const double* a4, 
			double* b1, double* b2,
			double* b3, double* b4) const
{
    for (int row = 0; row < n; row++) {
      double bval1 = 0;
      double bval2 = 0;
      double bval3 = 0;
      double bval4 = 0;
      for (int j = 0; j < nz[row]; j++) {
          double v = values[row][j];
          bval1 += a1[col[row][j]] * v;
          bval2 += a2[col[row][j]] * v;
          bval3 += a3[col[row][j]] * v;
          bval4 += a4[col[row][j]] * v;
      }
      b1[row] = bval1;
      b2[row] = bval2;
      b3[row] = bval3;
      b4[row] = bval4;
    }
}


