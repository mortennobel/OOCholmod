//
//  SparseMatrix.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include "sparse_matrix.h"
#include <cassert>

#include "dense_matrix.h"
#include "config_singleton.h"

using namespace std;

// bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL

namespace oocholmod {
    
    SparseMatrix::SparseMatrix(unsigned int nrow, unsigned int ncol, int maxSize)
    :sparse(nullptr), triplet(nullptr), nrow(nrow), ncol(ncol)
#ifdef DEBUG
    ,magicNumber(MAGIC_NUMBER)
#endif
    {
        if (maxSize == 0) {
            maxTripletElements = (nrow*(ncol+1))/2; // triangular number
        } else {
            maxTripletElements = maxSize;
        }
        this->symmetry = SYMMETRIC_UPPER;
#ifdef DEBUG
        assert(nrow == ncol); // must be square
#endif
    }
    
    SparseMatrix::SparseMatrix(cholmod_sparse *sparse)
    :sparse(sparse), triplet(nullptr), nrow((unsigned int)sparse->nrow), ncol((unsigned int)sparse->ncol), maxTripletElements(0)
#ifdef DEBUG
    ,magicNumber(MAGIC_NUMBER)
#endif
    {
        buildLookupIndexFromSparse();
    }
    
    SparseMatrix::SparseMatrix(SparseMatrix&& other)
    :sparse(other.sparse), triplet(other.triplet), nrow(other.nrow), ncol(other.ncol), lookupIndex(std::move(other.lookupIndex)), iRow(other.iRow), jColumn(other.jColumn), symmetry(other.symmetry), maxTripletElements(other.maxTripletElements)
#ifdef DEBUG
    ,magicNumber(other.magicNumber)
#endif
    {
        other.sparse = nullptr;
        other.triplet = nullptr;
        other.values = nullptr;
        other.iRow = nullptr;
        other.jColumn = nullptr;
        other.maxTripletElements = 0;
#ifdef DEBUG
        other.magicNumber = 0L;
#endif
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
            lookupIndex = std::move(other.lookupIndex);
            iRow = other.iRow;
            jColumn = other.jColumn;
            symmetry = other.symmetry;
            maxTripletElements = other.maxTripletElements;
#ifdef DEBUG
            magicNumber = other.magicNumber;
#endif
            other.sparse = nullptr;
            other.triplet = nullptr;
            other.values = nullptr;
            other.iRow = nullptr;
            other.jColumn = nullptr;
            other.maxTripletElements = 0;
#ifdef DEBUG
            other.magicNumber = 0L;
#endif
        }
        return *this;
    }
    
    
    SparseMatrix::~SparseMatrix(){
        if (sparse != nullptr || triplet != nullptr){
#ifdef DEBUG
            assert(magicNumber == MAGIC_NUMBER);
            magicNumber = 0;
#endif
            if (sparse != NULL){
                cholmod_free_sparse(&sparse, ConfigSingleton::getCommonPtr());
                sparse = NULL;
            }
            if (triplet != NULL){
                cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
                triplet = NULL;
            }
        }
    }
    
    void SparseMatrix::setSymmetry(Symmetry symmetry){
#ifdef DEBUG
        assert(sparse != NULL);
        assert(triplet != NULL);
#endif
        this->symmetry = symmetry;
    }
    
    void SparseMatrix::setNullSpace(const DenseMatrix& v){
#ifdef DEBUG
        assert(sparse != NULL);
        assert(magicNumber == MAGIC_NUMBER);
#endif
        // naive implementation: Todo run fast
        int idx = 0;
        for (int j=0;j<ncol;j++){
            int iFrom = ((int*)sparse->p)[j];
            int iTo = ((int*)sparse->p)[j+1]-1;
            for (int i=iFrom;i<=iTo;i++){
                int row = ((int*)sparse->i)[i];
                ((double*)sparse->x)[idx] *= v[row]*v[j];
                idx++;
            }
        }
        for (int i=0;i<v.getSize();i++){
            if (v[i] == 0){
                (*this)(i,i) = 1;
            }
        }
    }
    
    void SparseMatrix::build(bool readOnly){
#ifdef DEBUG
        assert(sparse == NULL);
        assert(magicNumber == MAGIC_NUMBER);
#endif
        sparse = cholmod_triplet_to_sparse(triplet, triplet->nnz, ConfigSingleton::getCommonPtr());
        cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
        triplet = NULL;
        values = NULL;
        iRow = NULL;
        jColumn = NULL;
        
        // build lookup index
#ifdef DEBUG
        assert(sparse->stype == symmetry);
        assert(sparse->packed);
#endif
        if (!readOnly){
            buildLookupIndexFromSparse();
        }
    }
    
    void SparseMatrix::buildLookupIndexFromSparse(){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
        lookupIndex.clear();
        // In packed form, the nonzero pattern of column j is in A->i [A->p [j] ... A->p [j+1]-1]
        int idx = 0;
        for (int j=0;j<ncol;j++){
            int iFrom = ((int*)sparse->p)[j];
            int iTo = ((int*)sparse->p)[j+1]-1;
            for (int i=iFrom;i<=iTo;i++){
                int row = ((int*)sparse->i)[i];
                lookupIndex[key(row, j)] = idx;
                idx++;
            }
        }
    }
    
    Factor SparseMatrix::analyze(){
#ifdef DEBUG
        assert(sparse != NULL);
        assert(magicNumber == MAGIC_NUMBER);
#endif
        cholmod_factor *L = cholmod_analyze(sparse, ConfigSingleton::getCommonPtr());
        return Factor(L);
    }
    
    
    void SparseMatrix::zero(){
#ifdef DEBUG
        assert(sparse != NULL);
        assert(magicNumber == MAGIC_NUMBER);
#endif
        memset(sparse->x, 0, sparse->nzmax * sizeof(double));
    }
    
    void SparseMatrix::print(const char* name){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
#endif
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
    }
    
    void SparseMatrix::assertValidIndex(unsigned int row, unsigned int column){
#ifdef DEBUG
        assert(symmetry == SYMMETRIC_UPPER);
        assert(row < nrow);
        assert(column < ncol);
        if (symmetry == SYMMETRIC_UPPER) {
            assert(row <= column);
        } else if (symmetry == SYMMETRIC_LOWER) {
            assert(row >= column);
        }
        
#endif
    }
    
    void SparseMatrix::assertHasSparse(){
#ifdef DEBUG
        assert(sparse != nullptr); // matrix must be build
#endif
    }
    
    void SparseMatrix::assertValidInitAddValue(unsigned int row, unsigned int column, double value){
#ifdef DEBUG
        assert(sparse == nullptr); // must be called before matrix build
        assert(triplet->nnz < maxTripletElements);
        assert(row < nrow);
        assert(column < ncol);
        if (symmetry == SYMMETRIC_UPPER) {
            assert(row <= column);
        } else if (symmetry == SYMMETRIC_LOWER) {
            assert(row >= column);
        }
        int shiftBits = sizeof(long)*8/2; // shift half of the bits of a long
        long maxId = (long)pow(2, shiftBits);
        assert (row < maxId);
        assert (column < maxId);
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
        return std::move(LHS);
    }
    
    SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS)
    {
        return std::move(RHS) + LHS;
    }
    
    SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS)
    {
        return std::move(LHS) + RHS;
    }
    
    SparseMatrix operator*(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
        cholmod_sparse *sparse = cholmod_ssmult(LHS.sparse, RHS.sparse, ASYMMETRIC, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator*(SparseMatrix&& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
        cholmod_sparse *sparse = cholmod_ssmult(LHS.sparse, RHS.sparse, ASYMMETRIC, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&LHS.sparse, ConfigSingleton::getCommonPtr());
        LHS.sparse = sparse;
        LHS.ncol = RHS.ncol;
        LHS.symmetry = ASYMMETRIC;
        return std::move(LHS);
    }
    
    SparseMatrix&& operator*(const SparseMatrix& LHS, SparseMatrix&& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
        cholmod_sparse *sparse = cholmod_ssmult(LHS.sparse, RHS.sparse, ASYMMETRIC, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&RHS.sparse, ConfigSingleton::getCommonPtr());
        RHS.sparse = sparse;
        RHS.ncol = RHS.ncol;
        RHS.symmetry = ASYMMETRIC;
        return std::move(RHS);
    }
    
    SparseMatrix&& operator*(SparseMatrix&& LHS, SparseMatrix&& RHS)
    {
        return std::move(LHS) * RHS;
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
        return std::move(LHS);
    }
    
    SparseMatrix operator*(const double& LHS, const SparseMatrix& RHS)
    {
        return RHS * LHS;
    }
    
    SparseMatrix&& operator*(const double& LHS, SparseMatrix&& RHS)
    {
        return std::move(RHS) * LHS;
    }
    
    
    DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS)
    {
        
    }
    
    DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.sparse && RHS.x);
        assert(LHS.magicNumber == MAGIC_NUMBER);
        assert(RHS.magicNumber == MAGIC_NUMBER);
        assert(LHS.ncol == RHS.nrow);
#endif
        double alpha[2] = {1.,1.};
        double beta[2] = {0.,0.};
        cholmod_dense *dense = cholmod_copy_dense(RHS.x, ConfigSingleton::getCommonPtr());
        cholmod_sdmult(LHS.sparse, false, alpha, beta, RHS.x, dense, ConfigSingleton::getCommonPtr());
        return DenseMatrix(dense);
    }
    
    DenseMatrix&& operator*(DenseMatrix&& LHS, const SparseMatrix& RHS);
    DenseMatrix operator*(SparseMatrix&& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator*(const SparseMatrix& LHS, DenseMatrix&& RHS);
    DenseMatrix operator*(const DenseMatrix& LHS, SparseMatrix&& RHS);
    DenseMatrix&& operator*(SparseMatrix&& LHS, DenseMatrix&& RHS);
    DenseMatrix&& operator*(DenseMatrix&& LHS, SparseMatrix&& RHS);
    
    
    void SparseMatrix::transpose()
    {
        assert(symmetry == ASYMMETRIC);
        sparse = cholmod_transpose(sparse, 1, ConfigSingleton::getCommonPtr());
        nrow = static_cast<int>(sparse->nrow);
        ncol = static_cast<int>(sparse->ncol);
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
        return std::move(M);
    }
}
