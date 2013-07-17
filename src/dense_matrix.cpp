//
//  DenseMatrix.cpp
//  OOCholmod
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include <vecLib/cblas.h>
#include <vecLib/clapack.h>

#include "dense_matrix.h"
#include "config_singleton.h"

namespace oocholmod {
    
    // bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL
    
    DenseMatrix::DenseMatrix(unsigned int rows, unsigned int cols, double value)
    :nrow{rows}, ncol{cols}
    {
        dense = cholmod_allocate_dense(rows, cols, rows /* leading dimension (equal rows) */ , CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        if (!isnan(value)) {
            fill(value);
        }
    }
    
    DenseMatrix::DenseMatrix(cholmod_dense *dense_)
    :dense{dense_}, nrow{static_cast<unsigned int>(dense_->nrow)}, ncol{static_cast<unsigned int>(dense_->ncol)}
    {
    }
    
    DenseMatrix::DenseMatrix(DenseMatrix&& move)
    :dense{move.dense}, nrow{move.nrow}, ncol{move.ncol}
    {
        move.dense = nullptr;
        move.nrow = 0;
        move.ncol = 0;
    }
    
    DenseMatrix& DenseMatrix::operator=(DenseMatrix&& other)
    {
        if (this != &other)
        {
            if (dense){
                cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
            }
            dense = other.dense;
            nrow = other.nrow;
            ncol = other.ncol;
            
            other.dense = nullptr;
            other.nrow = 0;
            other.ncol = 0;
        }
        return *this;
    }
    
    DenseMatrix::~DenseMatrix()
    {
        if (dense){
            cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
            dense = nullptr;
        }
    }
    
    void DenseMatrix::zero(){
        memset(dense->x, 0, nrow * ncol * sizeof(double));
    }
    
    void DenseMatrix::fill(double value)
    {
        double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                data[c*nrow + r] = value;
            }
        }
    }
    
    double DenseMatrix::dot(const DenseMatrix& b) const {
#ifdef DEBUG
        assert(ncol == 1 || nrow == 1);
        assert(b.ncol == ncol && b.nrow == nrow);
#endif
        return cblas_ddot(nrow*ncol, getData(), 1, b.getData(), 1);
    }
    
    double DenseMatrix::length() const {
#ifdef DEBUG
        assert(ncol == 1 || nrow == 1);
#endif
        return cblas_dnrm2(ncol*nrow, getData(), 1);
    }
    
    void DenseMatrix::elemDivide(const DenseMatrix& b, DenseMatrix& dest) const {
#ifdef DEBUG
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
    
    void DenseMatrix::elemDivide(const DenseMatrix& b){
        elemDivide(b, *this);
    }
    
    void DenseMatrix::elemMultiply(const DenseMatrix& b){
        elemMultiply(b, *this);
    }

    void DenseMatrix::elemMultiply(const DenseMatrix& b, DenseMatrix& dest) const {
#ifdef DEBUG
        assert(nrow == b.getRows() && ncol == b.getColumns());
#endif
        double *thisData = dest.getData();
        double *bData = b.getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                thisData[c*nrow + r] *= bData[c*nrow + r];
            }
        }
    }
    
    DenseMatrix DenseMatrix::copy() const{
        DenseMatrix dest(nrow,ncol);
        const double *srcPtr = getData();
        double *destPtr = dest.getData();
        memcpy(destPtr, srcPtr, nrow*ncol*sizeof(double));
        return dest;
    }
    
    void DenseMatrix::set(float *inData){
#ifdef DEBUG
        assert(dense);
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
        assert(dense);
#endif
        memcpy(dense->x, data, nrow*ncol*sizeof(double));
    }
    
    void DenseMatrix::get(double *outData) const {
#ifdef DEBUG
        assert(dense);
#endif
        memcpy(outData, dense->x, nrow*ncol*sizeof(double));
    }
    
    void DenseMatrix::get(float *outData) const {
        double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                outData[c*nrow + r] = (float)data[c*nrow + r];
            }
        }
    }
    
    void DenseMatrix::print(const char* name) const{
        cholmod_print_dense(dense, name, ConfigSingleton::getCommonPtr());
        int n_rows = (int)dense->nrow;
        int n_cols = (int)dense->ncol;
        for (int r = 0; r  < n_rows; r++)
        {
            for (int c = 0; c  < n_cols; c++)
            {
                std::cout << ((double*)dense->x)[c*n_rows + r] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    void DenseMatrix::swap(DenseMatrix& other){
        std::swap(dense, other.dense);
        std::swap(nrow, other.nrow);
        std::swap(ncol, other.ncol);
    }
    
    void swap(DenseMatrix& v1, DenseMatrix& v2) {
        v1.swap(v2);
    }
    
    // Addition
    DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        DenseMatrix res = RHS.copy();
        cblas_daxpy(LHS.nrow*LHS.ncol, 1., LHS.getData(), 1, res.getData(), 1);
        return res;
    }
    
    DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS)
    {
        return RHS+std::move(LHS);
    }
    
    DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        cblas_daxpy(LHS.nrow*LHS.ncol, 1., LHS.getData(), 1, RHS.getData(), 1);
        return std::move(RHS);
    }
    
    DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS)
    {
        return LHS+std::move(RHS);
    }
    
    DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& RHS)
    {
        RHS + std::move(*this);
        return *this;
    }
    
    // Multiplication
    DenseMatrix operator*(const DenseMatrix& LHS, const double& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense);
#endif
        cholmod_dense *dense = cholmod_copy_dense(LHS.dense, ConfigSingleton::getCommonPtr());
        DenseMatrix res(dense);
        cblas_dscal (LHS.nrow*LHS.ncol, RHS, res.getData(), 1);
        return res;
    }
    
    DenseMatrix&& operator*(DenseMatrix&& LHS, const double& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense);
#endif
        cblas_dscal (LHS.nrow*LHS.ncol, RHS, LHS.getData(), 1);
        return std::move(LHS);
    }
    
    DenseMatrix operator*(const double& LHS, const DenseMatrix& RHS)
    {
        return RHS*LHS;
    }
    
    DenseMatrix&& operator*(const double& LHS, DenseMatrix&& RHS)
    {
        return std::move(RHS)*LHS;
    }
    
    DenseMatrix& DenseMatrix::operator*=(const double& RHS)
    {
        std::move(*this) * RHS;
        return *this;
    }
    
    DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.ncol == RHS.nrow);
#endif
        DenseMatrix res(LHS.nrow, RHS.ncol);
        
        if(RHS.ncol == 1)
        {
            cblas_dgemv(CblasColMajor, CblasNoTrans, LHS.nrow, LHS.ncol, 1., LHS.getData(), LHS.nrow, RHS.getData(), 1, 0., res.getData(), 1);
        }
        else {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, LHS.nrow, RHS.ncol, LHS.ncol, 1., LHS.getData(), LHS.nrow, RHS.getData(), RHS.nrow, 0., res.getData(), res.nrow);
        }
        return res;
    }
    
    // Transpose
    void DenseMatrix::transpose()
    {
        double *data = getData();
        cholmod_dense *d = cholmod_allocate_dense(ncol, nrow, ncol, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        double *outData = (double*)d->x;
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                outData[r*ncol + c] = data[c*nrow + r];
            }
        }
        cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
        dense = d;
        
        int temp = nrow;
        nrow = ncol;
        ncol = temp;
    }
    
    DenseMatrix transposed(const DenseMatrix& M)
    {
        DenseMatrix res(M.ncol, M.nrow);
        double *data = M.getData();
        double *outData = res.getData();
        for (int c = 0; c < M.ncol; c++){
            for (int r = 0; r < M.nrow; r++){
                outData[r*M.ncol + c] = data[c*M.nrow + r];
            }
        }
        return res;
    }
    
    DenseMatrix&& transposed(DenseMatrix&& M)
    {
        M.transpose();
        return std::move(M);
    }
    
    DenseMatrix solve(const DenseMatrix& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        int N = b.nrow;
        int nrhs = b.ncol;
        int lda = A.nrow;
        int ldb = b.nrow;
        int ipiv[N];
        int info;
        
        cholmod_dense *a = cholmod_copy_dense(A.dense, ConfigSingleton::getCommonPtr());
        cholmod_dense *res = cholmod_copy_dense(b.dense, ConfigSingleton::getCommonPtr());
        
        dgesv_(&N, &nrhs, (double*)a->x, &lda, ipiv, (double*)res->x, &ldb, &info);
#ifdef DEBUG
        assert(info == 0);
#endif
        return DenseMatrix(res);
    }
    
    DenseMatrix solve(DenseMatrix&& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        int N = b.nrow;
        int nrhs = b.ncol;
        int lda = A.nrow;
        int ldb = b.nrow;
        int ipiv[N];
        int info;
        
        cholmod_dense *res = cholmod_copy_dense(b.dense, ConfigSingleton::getCommonPtr());
        dgesv_(&N, &nrhs, A.getData(), &lda, ipiv, (double*)res->x, &ldb, &info);
#ifdef DEBUG
        assert(info == 0);
#endif
        return DenseMatrix(res);
    }
    
    DenseMatrix&& solve(const DenseMatrix& A, DenseMatrix&& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        int N = b.nrow;
        int nrhs = b.ncol;
        int lda = A.nrow;
        int ldb = b.nrow;
        int ipiv[N];
        int info;
        
        cholmod_dense *a = cholmod_copy_dense(A.dense, ConfigSingleton::getCommonPtr());
        dgesv_(&N, &nrhs, (double*)a->x, &lda, ipiv, b.getData(), &ldb, &info);
#ifdef DEBUG
        assert(info == 0);
#endif
        return std::move(b);
    }
    
    DenseMatrix&& solve(DenseMatrix&& A, DenseMatrix&& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        int N = b.nrow;
        int nrhs = b.ncol;
        int lda = A.nrow;
        int ldb = b.nrow;
        int ipiv[N];
        int info;
        
        dgesv_(&N, &nrhs, A.getData(), &lda, ipiv, b.getData(), &ldb, &info);
#ifdef DEBUG
        assert(info == 0);
#endif
        return std::move(b);
    }
    
}
