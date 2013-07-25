//
//  blas.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#ifndef __OOCholmod__blas__
#define __OOCholmod__blas__

#include <iostream>

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
	AtlasConj=114};

double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);

double cblas_dnrm2(const int N, const double *X, const int incX);

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);

void cblas_dscal(const int N, const double alpha, double *X, const int incX);

void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY);

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);



#endif /* defined(__OOCholmod__blas__) */
