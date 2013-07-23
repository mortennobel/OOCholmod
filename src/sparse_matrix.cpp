//
//  SparseMatrix.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include <cassert>

#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "factor.h"

using namespace std;

namespace oocholmod {
    
    SparseMatrix::SparseMatrix(unsigned int nrow, unsigned int ncol, bool symmetric, int maxSize)
    :sparse{nullptr}, triplet{nullptr}, nrow{nrow}, ncol{ncol}
    {
        if (symmetric && nrow == ncol) {
            symmetry = SYMMETRIC_UPPER;
        }
        else {
            symmetry = ASYMMETRIC;
        }
        
        if (maxSize == 0 && symmetric) {
            maxTripletElements = (nrow*(ncol+1))/2; // triangular number
        }
        else if (maxSize == 0 && !symmetric) {
            maxTripletElements = nrow*ncol; // triangular number
        }
        else {
            maxTripletElements = maxSize;
        }
    }
    
    SparseMatrix::SparseMatrix(cholmod_sparse *sparse)
    :sparse{sparse}, triplet{nullptr}, nrow{static_cast<unsigned int>(sparse->nrow)},
    ncol{static_cast<unsigned int>(sparse->ncol)},
    values{(double*)sparse->x}, iRow{(int*)sparse->i}, jColumn{(int*)sparse->p}, symmetry{static_cast<Symmetry>(sparse->stype)}, maxTripletElements{0}
    {
#ifdef DEBUG
        assert(sparse->itype == CHOLMOD_INT);
#endif
    }
    
    SparseMatrix::SparseMatrix(SparseMatrix&& other)
    :sparse{other.sparse}, triplet{other.triplet}, nrow{other.nrow}, ncol{other.ncol}, values{other.values}, iRow{other.iRow}, jColumn{other.jColumn}, symmetry{other.symmetry}, maxTripletElements{other.maxTripletElements}
    {
        other.sparse = nullptr;
        other.triplet = nullptr;
        other.values = nullptr;
        other.iRow = nullptr;
        other.jColumn = nullptr;
        other.nrow = 0;
        other.ncol = 0;
    }
    
    SparseMatrix& SparseMatrix::operator=(SparseMatrix&& other){
        if (this != &other){
            if (sparse != nullptr){
                cholmod_free_sparse(&sparse, ConfigSingleton::getCommonPtr());
            }
            if (triplet != nullptr){
                cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
            }
            
            sparse = other.sparse;
            triplet = other.triplet;
            nrow = other.nrow;
            ncol = other.ncol;
            values = other.values;
            iRow = other.iRow;
            jColumn = other.jColumn;
            symmetry = other.symmetry;
            maxTripletElements = other.maxTripletElements;

            other.sparse = nullptr;
            other.triplet = nullptr;
            other.values = nullptr;
            other.iRow = nullptr;
            other.jColumn = nullptr;
            other.nrow = 0;
            other.ncol = 0;
        }
        return *this;
    }
    
    
    SparseMatrix::~SparseMatrix(){
        if (sparse != nullptr || triplet != nullptr){

            if (sparse != nullptr){
                cholmod_free_sparse(&sparse, ConfigSingleton::getCommonPtr());
                sparse = nullptr;
            }
            if (triplet != nullptr){
                cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
                triplet = nullptr;
            }
        }
    }
    
    DenseMatrix SparseMatrix::toDense()
    {
#ifdef DEBUG
        assert(sparse);
#endif
        cholmod_dense *dense = cholmod_sparse_to_dense(sparse, ConfigSingleton::getCommonPtr());
        return DenseMatrix(dense);
    }

    bool SparseMatrix::hasElement(unsigned int row, unsigned int column) const {
        return getIndex(row, column) != -1;
    }
    
    void SparseMatrix::dropSmallEntries(double tol){
#ifdef DEBUG
        assert(getMatrixState() == BUILT);
#endif
        cholmod_drop(tol, sparse, ConfigSingleton::getCommonPtr());
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);
    }
    
    MatrixState SparseMatrix::getMatrixState() const {
        if (nrow == 0 && ncol == 0){
            return DESTROYED;
        }
        if (sparse == nullptr && triplet == nullptr){
            return UNINITIALIZED;
        }
        if (triplet != nullptr){
            return INIT;
        }
        return BUILT;
    }
    
    void SparseMatrix::setSymmetry(Symmetry symmetry){
#ifdef DEBUG
        assert(sparse == nullptr);
        assert(triplet == nullptr);
#endif
        this->symmetry = symmetry;
    }
    
    void SparseMatrix::build(){
#ifdef DEBUG
        assert(sparse == nullptr);
#endif
        sparse = cholmod_triplet_to_sparse(triplet, triplet->nnz, ConfigSingleton::getCommonPtr());
        cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
        triplet = nullptr;
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);

        
#ifdef DEBUG
        assert(sparse->itype == CHOLMOD_INT);
        assert(sparse->stype == symmetry);
        assert(sparse->packed);
#endif
    }
    
    bool SparseMatrix::operator==(const SparseMatrix& RHS)
    {
        for(int r = 0; r < nrow; r++)
        {
            for(int c = 0; c < ncol; c++)
            {
                if((*this)(r,c) != RHS(r,c))
                {
                    return false;
                }
            }
        }
        return true;
    }
    
    void SparseMatrix::swap(SparseMatrix& other){
        std::swap(sparse, other.sparse);
        std::swap(triplet, other.triplet);
        std::swap(nrow, other.nrow);
        std::swap(ncol, other.ncol);
        std::swap(values, other.values);
        std::swap(iRow, other.iRow);
        std::swap(jColumn, other.jColumn);
        std::swap(symmetry, other.symmetry);
        std::swap(maxTripletElements, other.maxTripletElements);
    }
    
    void swap(SparseMatrix& v1, SparseMatrix& v2) {
        v1.swap(v2);
    }
   
    
    Factor SparseMatrix::analyze() const
    {
#ifdef DEBUG
        assertHasSparse();
#endif
        cholmod_factor *L = cholmod_analyze(sparse, ConfigSingleton::getCommonPtr());
        return Factor(L);
    }
    
    void SparseMatrix::zero(){
#ifdef DEBUG
        assertHasSparse();
#endif
        memset(values, 0, sparse->nzmax * sizeof(double));
    }
    
    void SparseMatrix::print(const char* name) const {
        if (sparse){
            cholmod_print_sparse(sparse, name, ConfigSingleton::getCommonPtr());
            cholmod_dense *dense = cholmod_sparse_to_dense(sparse, ConfigSingleton::getCommonPtr());
            int n_rows = (int)dense->nrow;
            int n_cols = (int)dense->ncol;
            for (int r = 0; r  < n_rows; r++)
            {
                for (int c = 0; c  < n_cols; c++)
                {
                    cout << ((double*)dense->x)[c*n_rows + r] << " ";
                }
                cout << endl;
            }
            cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
            cout << endl;
            cout << "Packed "<<sparse->packed<< endl;
            cout << "p: ";
            for (int i=0;i<=sparse->ncol;i++){
                printf("%4i ", ((int*)sparse->p)[i]);
            }
            cout << endl;
            cout << "i: ";
            for (int i=0;i<sparse->nzmax;i++){
                printf("%4i ", ((int*)sparse->i)[i]);
            }
            cout << endl;
            cout << "x: ";
            for (int i=0;i<sparse->nzmax;i++){
                printf("%3.3f ", ((double*)sparse->x)[i]);
            }
            cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
            cout << endl;
        }
        else if (triplet){
            cholmod_print_triplet(triplet, name, ConfigSingleton::getCommonPtr());
            for (int i=0;i<triplet->nnz;i++){
                cout << ((int*)triplet->i)[i]<<" "<< ((int*)triplet->j)[i]<<" "<< ((double*)triplet->x)[i]<<" "<<endl;
            }
        } else {
            cout << "[Empty sparse matrix]"<<endl;
        }
    }
    
    void SparseMatrix::assertValidIndex(unsigned int row, unsigned int column) const
    {
#ifdef DEBUG
        assert(row < nrow);
        assert(column < ncol);
#endif
    }
    
    double SparseMatrix::norm(int norm) const{
#ifdef DEBUG
        assertHasSparse();
#endif
       return cholmod_norm_sparse(sparse, norm, ConfigSingleton::getCommonPtr());
    }
    
    SparseMatrix SparseMatrix::copy() const{
        SparseMatrix res;
        if (sparse){
            res.sparse = cholmod_copy_sparse(sparse, ConfigSingleton::getCommonPtr());
            res.values = ((double*)res.sparse->x);
            res.iRow = ((int*)res.sparse->i);
            res.jColumn = ((int*)res.sparse->p);
        }
        if (triplet){
            res.triplet = cholmod_copy_triplet(triplet, ConfigSingleton::getCommonPtr());
            res.values = (double *)res.triplet->x;
            res.iRow = (int *)res.triplet->i;
            res.jColumn = (int *)res.triplet->j;
        }
        res.nrow = nrow;
        res.ncol = ncol;
        res.symmetry = symmetry;
        res.maxTripletElements = maxTripletElements;
        return move(res);
    }
    
    void SparseMatrix::assertHasSparse() const
    {
#ifdef DEBUG
        assert(sparse != nullptr); // matrix must be build
#endif
    }
    
    void SparseMatrix::increaseTripletCapacity(){
        // grow size with factor 1.5
        maxTripletElements = (int)ceil(1.5f*maxTripletElements);
        auto newTriplet = cholmod_allocate_triplet(nrow, ncol, maxTripletElements, symmetry, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        memcpy(newTriplet->x, triplet->x, triplet->nzmax*sizeof(double));
        memcpy(newTriplet->i, triplet->i, triplet->nzmax*sizeof(int));
        memcpy(newTriplet->j, triplet->j, triplet->nzmax*sizeof(int));
        newTriplet->nnz = triplet->nnz;
        cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
        triplet = newTriplet;
        values = (double *)triplet->x;
        iRow = (int *)triplet->i;
        jColumn = (int *)triplet->j;
    }
    
    void SparseMatrix::assertValidInitAddValue(unsigned int row, unsigned int column) const {
#ifdef DEBUG
        assert(sparse == nullptr); // must be called before matrix build
        assert(row < nrow);
        assert(column < ncol);
#endif
    }
    
    
    SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double scale[2] = {1.,1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, scale, scale, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double scale[2] = {1.,1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, scale, scale, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&LHS.sparse, ConfigSingleton::getCommonPtr());
        LHS.sparse = sparse;
        LHS.values = ((double*)sparse->x);
        LHS.iRow = ((int*)sparse->i);
        LHS.jColumn = ((int*)sparse->p);
        return move(LHS);
    }
    
    SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS)
    {
        return move(RHS) + LHS;
    }
    
    SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS)
    {
        return move(LHS) + RHS;
    }
    
    SparseMatrix operator-(const SparseMatrix& M)
    {
        return -1.*M;
    }
    
    SparseMatrix&& operator-(SparseMatrix&& M)
    {
        return -1.*move(M);
    }
    
    SparseMatrix operator-(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double alpha[2] = {1.,1.};
        double beta[2] = {-1.,-1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, alpha, beta, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator-(SparseMatrix&& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double alpha[2] = {1.,1.};
        double beta[2] = {-1.,-1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, alpha, beta, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&LHS.sparse, ConfigSingleton::getCommonPtr());
        LHS.sparse = sparse;
        LHS.values = ((double*)sparse->x);
        LHS.iRow = ((int*)sparse->i);
        LHS.jColumn = ((int*)sparse->p);
        return move(LHS);
    }
    
    SparseMatrix&& operator-(const SparseMatrix& LHS, SparseMatrix&& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double alpha[2] = {1.,1.};
        double beta[2] = {-1.,-1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, alpha, beta, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&RHS.sparse, ConfigSingleton::getCommonPtr());
        RHS.sparse = sparse;
        RHS.values = ((double*)sparse->x);
        RHS.iRow = ((int*)sparse->i);
        RHS.jColumn = ((int*)sparse->p);
        return move(RHS);
    }
    
    SparseMatrix&& operator-(SparseMatrix&& LHS, SparseMatrix&& RHS)
    {
        return move(LHS) - RHS;
    }
    
    SparseMatrix operator*(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
        cholmod_sparse *sparse = cholmod_ssmult(LHS.sparse, RHS.sparse, ASYMMETRIC, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix operator*(const SparseMatrix& LHS, const double& RHS)
    {
        assert(LHS.sparse);
        cholmod_dense *dense = cholmod_zeros(1, 1, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        ((double*)dense->x)[0] = RHS;
        cholmod_sparse *sparse = cholmod_copy_sparse(LHS.sparse, ConfigSingleton::getCommonPtr());
        cholmod_scale(dense, CHOLMOD_SCALAR, sparse, ConfigSingleton::getCommonPtr());
        cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator*(SparseMatrix&& LHS, const double& RHS)
    {
        assert(LHS.sparse);
        cholmod_dense *dense = cholmod_zeros(1, 1, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        ((double*)dense->x)[0] = RHS;
        cholmod_scale(dense, CHOLMOD_SCALAR, LHS.sparse, ConfigSingleton::getCommonPtr());
        cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
        return move(LHS);
    }
    
    SparseMatrix operator*(const double& LHS, const SparseMatrix& RHS)
    {
        return RHS * LHS;
    }
    
    SparseMatrix&& operator*(const double& LHS, SparseMatrix&& RHS)
    {
        return move(RHS) * LHS;
    }
    
    
    DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
#endif
        double alpha[2] = {1.,1.};
        double beta[2] = {0.,0.};
        DenseMatrix res(RHS.ncol, LHS.nrow);
        cholmod_sdmult(RHS.sparse, true, alpha, beta, transposed(LHS).dense, res.dense, ConfigSingleton::getCommonPtr());
        return transposed(res);
    }
    
    DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.sparse && RHS.dense);
        assert(LHS.ncol == RHS.nrow);
#endif
        double alpha[2] = {1.,1.};
        double beta[2] = {0.,0.};
        DenseMatrix res(LHS.nrow, RHS.ncol);
        cholmod_sdmult(LHS.sparse, false, alpha, beta, RHS.dense, res.dense, ConfigSingleton::getCommonPtr());
        return res;
    }
    
    void SparseMatrix::transpose()
    {
        assert(symmetry == ASYMMETRIC);
        sparse = cholmod_transpose(sparse, 1, ConfigSingleton::getCommonPtr());
        nrow = static_cast<int>(sparse->nrow);
        ncol = static_cast<int>(sparse->ncol);
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);
    }
    
    SparseMatrix transposed(const SparseMatrix& M)
    {
        assert(M.symmetry == ASYMMETRIC);
        cholmod_sparse *sparse = cholmod_transpose(M.sparse, 1, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& transposed(SparseMatrix&& M)
    {
        assert(M.symmetry == ASYMMETRIC);
        cholmod_sparse *sparse = cholmod_transpose(M.sparse, 1, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&M.sparse, ConfigSingleton::getCommonPtr());
        M.sparse = sparse;
        M.nrow = static_cast<int>(sparse->nrow);
        M.ncol = static_cast<int>(sparse->ncol);
        M.values = ((double*)sparse->x);
        M.iRow = ((int*)sparse->i);
        M.jColumn = ((int*)sparse->p);
        return move(M);
    }
    
    DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.sparse && b.dense);
        assert(A.nrow == b.nrow);
#endif
        Factor F = A.analyze();
        F.factorize(A);
        return solve(F, b);
    }
    
    SparseMatrix solve(const SparseMatrix& A, const SparseMatrix& b)
    {
#ifdef DEBUG
        assert(A.sparse && b.sparse);
        assert(A.nrow == b.nrow);
#endif
        Factor F = A.analyze();
        F.factorize(A);
        return solve(F, b);
    }
}
