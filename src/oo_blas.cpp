//
//  blas.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#include "oo_blas.h"

#include <cassert>

#ifdef NO_BLAS

using namespace std;

double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY){
    assert(false); // not implemented
}

double cblas_dnrm2(const int N, const double *X, const int incX){
    assert(false); // not implemented
}

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){
    assert(false); // not implemented
}

void cblas_dscal(const int N, const double alpha, double *X, const int incX){
    assert(false); // not implemented
}

void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){
    assert(false); // not implemented
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){
    assert(false); // not implemented
}
#endif
