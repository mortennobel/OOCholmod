//
//  DenseMatrix.cpp
//  OOCholmod
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include "dense_matrix.h"
#include <cassert>
#include <vecLib/cblas.h>
#include "config_singleton.h"

namespace oocholmod {
    
    // bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL
    
    DenseMatrix::DenseMatrix(unsigned int rows, unsigned int cols)
    :nrow(rows), ncol(cols)
#ifdef DEBUG
    ,magicNumber(MAGIC_NUMBER)
#endif
    {
        x = cholmod_allocate_dense(rows, cols, rows /* leading dimension (equal rows) */ , CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
    }
    
    DenseMatrix::DenseMatrix(cholmod_dense *x)
    :x(x), nrow(static_cast<unsigned int>(x->nrow)), ncol(static_cast<unsigned int>(x->ncol))
#ifdef DEBUG
    ,magicNumber(MAGIC_NUMBER)
#endif
    {
    }
    
    DenseMatrix::DenseMatrix(DenseMatrix&& move)
    :x(move.x), nrow(move.nrow), ncol(move.ncol)
#ifdef DEBUG
    ,magicNumber(move.magicNumber)
#endif
    {
        move.x = nullptr;
        move.nrow = 0;
        move.ncol = 0;
#ifdef DEBUG
        magicNumber = 0;
#endif
    }
    
    DenseMatrix& DenseMatrix::operator=(DenseMatrix&& other)
    {
        if (this != &other)
        {
            if (x != nullptr){
                cholmod_free_dense(&x, ConfigSingleton::getCommonPtr());
            }
            x = other.x;
            nrow = other.nrow;
            ncol = other.ncol;
#ifdef DEBUG
            magicNumber = other.magicNumber;
#endif
            
            other.x = nullptr;
            other.nrow = 0;
            other.ncol = 0;
#ifdef DEBUG
            other.magicNumber = 0;
#endif
            
        }
        return *this;
    }
    
    DenseMatrix::~DenseMatrix(){
        if (x != nullptr){
#ifdef DEBUG
            assert(magicNumber == MAGIC_NUMBER);
            magicNumber = 0;
#endif
            cholmod_free_dense(&x, ConfigSingleton::getCommonPtr());
            x = NULL;
        }
    }
    
    void DenseMatrix::zero(){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(x);
#endif
        memset(x->x, 0, nrow * ncol * sizeof(double));
    }
    
    double DenseMatrix::dot(const DenseMatrix& b){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        return cblas_ddot(nrow*ncol, getData(), 1, b.getData(), 1);
    }
    
    void DenseMatrix::fill(double value){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                data[c*nrow + r] = value;
            }
        }
    }
    
    double DenseMatrix::length(){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(ncol == 1 || nrow == 1);
#endif
        return cblas_dnrm2(ncol*nrow, getData(), 1);
    }
    
    void DenseMatrix::scale(double alpha){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        cblas_dscal (nrow*ncol, alpha, getData(), 1);
    }
    
    void DenseMatrix::elem_divide(const DenseMatrix& b){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(nrow == b.getRows() && ncol == b.getColumns());
#endif
        double *thisData = getData();
        double *bData = b.getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                thisData[c*nrow + r] /= bData[c*nrow + r];
            }
        }
    }
    
    void DenseMatrix::elem_multiply(const DenseMatrix& b){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(nrow == b.getRows() && ncol == b.getColumns());
#endif
        double *thisData = getData();
        double *bData = b.getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                thisData[c*nrow + r] *= bData[c*nrow + r];
            }
        }
    }
    
    DenseMatrix DenseMatrix::copy(){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        DenseMatrix dest(nrow,ncol);
        const double *srcPtr = getData();
        double *destPtr = dest.getData();
        memcpy(destPtr, srcPtr, nrow*ncol*sizeof(double));
        return dest;
    }
    
    
    void DenseMatrix::set(float *inData){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(x);
#endif
        double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                data[c*nrow + r] = inData[c*nrow + r];
            }
        }
    }
    
    void DenseMatrix::set(double *data){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(x != NULL);
#endif
        memcpy(x->x, data, nrow*ncol*sizeof(double));
    }
    
    void DenseMatrix::get(double *outData){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        assert(x != NULL);
#endif
        memcpy(outData, x->x, nrow*ncol*sizeof(double));
    }
    
    void DenseMatrix::get(float *outData){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                outData[c*nrow + r] = (float)data[c*nrow + r];
            }
        }
    }
    
    void DenseMatrix::print(const char* name) const{
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        cholmod_print_dense(x, name, ConfigSingleton::getCommonPtr());
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
    
    // Addition
    DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS)
    {
        assert(LHS.x && RHS.x);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
    }
    
    DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS)
    {
        assert(LHS.x && RHS.x);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
    }
    
    DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS);
    DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS);
    
    // Multiplication
    DenseMatrix operator*(const DenseMatrix& LHS, const double& RHS);
    DenseMatrix&& operator*(DenseMatrix&& LHS, const double& RHS);
    DenseMatrix operator*(const double& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator*(const double& LHS, DenseMatrix&& RHS);
    
    DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator*(DenseMatrix&& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator*(const DenseMatrix& LHS, DenseMatrix&& RHS);
    DenseMatrix&& operator*(DenseMatrix&& LHS, DenseMatrix&& RHS);
    
    // Transpose
    DenseMatrix transposed(const DenseMatrix& M);
    DenseMatrix&& transposed(DenseMatrix&& M);
    
}
