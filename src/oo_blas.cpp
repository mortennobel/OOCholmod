//
//  blas.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#include "oo_blas.h"

#include <cassert>
#include <cmath>
#include <iostream>

#ifdef NO_BLAS

using namespace std;

double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY){
    assert(incX == 1);
    assert(incY == 1);
    double dotProduct = 0;
    for (int i=0;i<N;i++){
        dotProduct += X[i]*Y[i];
    }
    return dotProduct;
}

double cblas_dnrm2(const int N, const double *X, const int incX){
    assert(incX == 1);
    double sqSum = 0;
    for (int i=0;i<N;i++){
        sqSum += X[i]*X[i];
    }
    return sqrt(sqSum);
}

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){
    assert(alpha == 1);
    assert(incX == 1);
    assert(incY == 1);
    for (int i=0;i<N;i++){
        Y[i] += X[i];
    }
}

void cblas_dscal(const int N, const double alpha, double *X, const int incX){
    assert(incX == 1);
    for (int i=0;i<N;i++){
        X[i] *= alpha;
    }
}

void cblas_dgemv(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){
    assert(Order == CblasColMajor);
    assert(TransA == CblasNoTrans);
    assert(alpha == 1);
    assert(beta == 0);
    for (int row=0;row<M;row++){
        double sum = 0;
        for (int column=0;column<N;column++){
            cout << A[row*M + column]<<" x "<<X[column]<<endl;
            sum += A[row*M + column] * X[column];
        }
        cout << " set row " << row << " to " << sum << endl;
        Y[row] = sum;
    }
}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){
    assert(Order == CblasColMajor);
    assert(TransA == CblasNoTrans);
    assert(alpha == 1);
    assert(beta == 0);
    for (int row=0;row<M;row++){
        for (int column=0;column<N;column++){
            double sum = 0;
            for (int j=0;j<K;j++){
                sum += A[row+M*j]*B[j+K*column];
            }
            C[row+column*M] = sum;
        }
    }
}
#endif
