//
//  CholmodDenseVector.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//  License: LGPL 3.0 

#include "CholmodDenseVector.h"
#include <cassert>
#include <vecLib/cblas.h>


CholmodDenseVector::CholmodDenseVector(int size, cholmod_common *Common)
:Common(Common), size(size)
{
    x = cholmod_allocate_dense(size, 1, size /* leading dimension (equal rows) */ , CHOLMOD_REAL, Common);
}

CholmodDenseVector::CholmodDenseVector(cholmod_dense *x, cholmod_common *Common, int size)
:x(x), Common(Common), size(size)
{
}

CholmodDenseVector::~CholmodDenseVector(){
    cholmod_free_dense(&x, Common);
    x = NULL;
}

void CholmodDenseVector::zero(){
#ifdef DEBUG
    assert(x != NULL);
#endif
    memset(x->x, 0, size * sizeof(double));
}

double CholmodDenseVector::dot(CholmodDenseVector *b){
    return cblas_ddot(getSize(), getData(), 1, b->getData(), 1);
}

void CholmodDenseVector::fill(double value){
    double *data = getData();
    for (int i=0;i<size;i++){
        data[i] = value;
    }
}

double CholmodDenseVector::length(){
    return cblas_dnrm2(getSize(), getData(), 1);
}

void CholmodDenseVector::scale(double alpha){
    cblas_dscal (getSize(), alpha, getData(), 1);
}

void CholmodDenseVector::divideBy(CholmodDenseVector *b){
    assert(b->getSize() >= getSize());
    double *thisData = getData();
    double *bData = b->getData();
    for (int i=0;i<getSize();i++){
        thisData[i] /= bData[i];
    }
}

void CholmodDenseVector::multiplyWith(CholmodDenseVector *b){
    assert(b->getSize() >= getSize());
    double *thisData = getData();
    double *bData = b->getData();
    for (int i=0;i<getSize();i++){
        thisData[i] *= bData[i];
    }
}

void CholmodDenseVector::copyTo(CholmodDenseVector *dest){
    assert(dest->getSize() >= getSize());
    const double *srcPtr = getData();
    double *destPtr = dest->getData();
    memcpy(destPtr, srcPtr, sizeof(double) * getSize());
}

void CholmodDenseVector::set(float *inData){
#ifdef DEBUG
    assert(x != NULL);
#endif
    double *data = getData();
    for (int i=0;i<size;i++){
        data[i] = inData[i];
    }
}

void CholmodDenseVector::set(double *data){
#ifdef DEBUG
    assert(x != NULL);
#endif
    memcpy(x->x, data, sizeof(double)*size);
}

void CholmodDenseVector::get(double *outData){
#ifdef DEBUG
    assert(x != NULL);
#endif
    memcpy(outData, x->x, sizeof(double)*size);
}

void CholmodDenseVector::get(float *outData){
    double *data = getData();
    for (int i=0;i<size;i++){
        outData[i] = (float)data[i];
    }
}

void CholmodDenseVector::print(const char* name){
    cholmod_print_dense(x, name, Common);
    int n_rows = (int)x->nrow;
    int n_cols = (int)x->ncol;
    for (int r = 0; r  < n_rows; r++)
    {
        for (int c = 0; c  < n_cols; c++)
        {
            std::cout << ((double*)x->x)[c*n_rows + r] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}