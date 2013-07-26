//
//  DenseMatrix.cpp
//  OOCholmod
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#ifdef NO_BLAS
#include "oo_blas.h"
#else
#include <vecLib/cblas.h>
#endif
#ifdef NO_LAPACK
#include "oo_lapack.h"
#else
#include <vecLib/clapack.h>
#endif

#include "dense_matrix.h"
#include "config_singleton.h"
#include "sparse_matrix.h"

using namespace std;

namespace oocholmod {
   
    DenseMatrix::DenseMatrix(unsigned int rows, unsigned int cols, double value)
    :nrow{rows}, ncol{cols}
    {
        dense = cholmod_allocate_dense(rows, cols, rows /* leading dimension (equal rows) */ , CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        if (!std::isnan(value)) {
            fill(value);
        } else {
#ifdef DEBUG
            fill(NAN);
#endif
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
    
    SparseMatrix DenseMatrix::toSparse() const {
        return SparseMatrix(cholmod_dense_to_sparse(dense, true, ConfigSingleton::getCommonPtr()));
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
    
    double DenseMatrix::norm(int norm) const{
        return cholmod_norm_dense(dense, norm, ConfigSingleton::getCommonPtr());
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
        return RHS+move(LHS);
    }
    
    DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        cblas_daxpy(LHS.nrow*LHS.ncol, 1., LHS.getData(), 1, RHS.getData(), 1);
        return move(RHS);
    }
    
    DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS)
    {
        return LHS+move(RHS);
    }
    
    DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& RHS)
    {
        RHS + move(*this);
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
        return move(LHS);
    }
    
    DenseMatrix operator*(const double& LHS, const DenseMatrix& RHS)
    {
        return RHS*LHS;
    }
    
    DenseMatrix&& operator*(const double& LHS, DenseMatrix&& RHS)
    {
        return move(RHS)*LHS;
    }
    
    DenseMatrix& DenseMatrix::operator*=(const double& RHS)
    {
        move(*this) * RHS;
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
    
    bool DenseMatrix::operator==(const DenseMatrix& RHS) const {
        if (nrow != RHS.nrow || ncol != RHS.ncol){
            return false;
        }
        double *data = getData();
        double *rhsData = RHS.getData();
        for (int i=0;i<nrow*ncol;i++){
            if (data[i] != rhsData[i]){
                return false;
            }
        }
        return true;
    }
    
    bool DenseMatrix::operator!=(const DenseMatrix& RHS) const{
        return !(*this == RHS);
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
        return move(M);
    }
    
    DenseMatrix solve(const DenseMatrix& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer ipiv[N];
        __CLPK_integer info;
        
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
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer ipiv[N];
        __CLPK_integer info;
        
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
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer ipiv[N];
        __CLPK_integer info;
        
        cholmod_dense *a = cholmod_copy_dense(A.dense, ConfigSingleton::getCommonPtr());
        dgesv_(&N, &nrhs, (double*)a->x, &lda, ipiv, b.getData(), &ldb, &info);
#ifdef DEBUG
        assert(info == 0);
#endif
        return move(b);
    }
    
    DenseMatrix&& solve(DenseMatrix&& A, DenseMatrix&& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer ipiv[N];
        __CLPK_integer info;
        
        dgesv_(&N, &nrhs, A.getData(), &lda, ipiv, b.getData(), &ldb, &info);
#ifdef DEBUG
        assert(info == 0);
#endif
        return move(b);
    }
        
    ostream& operator<<(ostream& os, const DenseMatrix& A)
    {
        os << endl;
        if (A.dense)
        {
            cholmod_print_dense(A.dense, "", ConfigSingleton::getCommonPtr());
            os << "[";
            for (int r = 0; r < A.getRows(); r++)
            {
                for (int c = 0; c < A.getColumns(); c++)
                {
                    os << A(r,c);
                    if(c < A.getColumns()-1) {
                        os << ", ";
                    }
                }
                if(r < A.getRows()-1)
                {
                    os << ";" << endl << " ";
                }
            }
            os << "];" << endl;
        }
        else {
            os << "[Empty dense matrix]" << endl;
        }
        return os;
    }
}
