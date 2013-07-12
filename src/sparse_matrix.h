//
//  sparse_matrix.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0 

#pragma once

#include <iostream>
#include <set>
#include <map>
#include <cmath>

#include <cholmod.h>
#include "factor.h"

namespace oocholmod {
    
// forward declaration
class CholmodDenseVector;

enum Symmetry {
    SYMMETRIC_LOWER = -1, // Lower triangular part stored
    ASYMMETRIC = 0,
    SYMMETRIC_UPPER = 1, // Upper triangular part stored
};

///
/// Currently only real (double), upper symmetric matrices are supported.
///
/// The sparse matrix must be used in the following way:
/// 1. Fill the matrix elements using the initAddValue() method
/// 2. Call build()
/// 3. Fill matrix with elements using setValue or addValue
///
/// At any point after the matrix has been build, you can call 
class SparseMatrix {
public:
    /// nrow # of rows of A
    /// ncol # of columns of A
    /// maxSize (size allocated before build). 0 means triangular
    SparseMatrix(cholmod_sparse *sparse);
    SparseMatrix(SparseMatrix&& move);
    SparseMatrix& operator=(SparseMatrix&& other);
    
    virtual ~SparseMatrix();
    
    friend SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS);
    friend SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS);
    friend SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS);
    friend SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS);
    
    void build(bool readOnly = false);
    
    CholmodFactor *analyzePtr();
    CholmodFactor analyze();
    
    void zero();
    
    Symmetry getSymmetry() { return symmetry; }
    
    int getRows(){ return nrow; }
    
    int getColumns(){ return ncol; }
    
    /// Set the nullspace
    /// Null = sparseDiagonal(N)
    /// K = N^T * K * N - (Null-I)
    ///
    /// In other words it gives matrixes with the following patterns
    /// where 0 elements in N marks the columns and rows to 'null'
    ///
    ///     0
    ///     0
    /// 00001000
    ///     0
    ///
    void setNullSpace(CholmodDenseVector *N);
    void setNullSpace(CholmodDenseVector& N);
    
    inline void initAddValue(unsigned int row, unsigned int column, double value=0) {
        assertValidInitAddValue(row, column, value);
        long k = key(row, column);
        auto res = lookupIndex.find(k);
        if (res != lookupIndex.end()) {
            values[(*res).second] += value;
            return;
        }
        lookupIndex[k] = (int)triplet->nnz;
        iRow[triplet->nnz] = row;
        jColumn[triplet->nnz] = column;
        values[triplet->nnz] = value;
        
        triplet->nnz++;
    }
    
    // computes alpha*(A*X) + beta*Y
    // res is result
    // alpha is optional (default 1)
    // beta is optional (default 0)
    void multiply(CholmodDenseVector *X, CholmodDenseVector *res, double alpha = 1, double beta = 0);
    void multiply(CholmodDenseVector& X, CholmodDenseVector& res, double alpha = 1, double beta = 0);
    CholmodDenseVector multiply(CholmodDenseVector& X, double alpha = 1, double beta = 0);
    
    inline double getValue(unsigned int row, unsigned int column){
        assertHasSparse();
        int index = getIndex(row, column);
        return ((double*)sparse->x)[index];
    }

    inline void addValue(unsigned int row, unsigned int column, double value){
        assertHasSparse();
        int index = getIndex(row, column);
        ((double*)sparse->x)[index] += value;
    }
    
    
    /// Set value of sparse matrix
    /// Must be invoked after build
    inline void setValue(unsigned int row, unsigned int column, double value){
        assertHasSparse();
        int index = getIndex(row, column);
        ((double*)sparse->x)[index] += value;
    }
    
    
    /// Get cholmod_sparse pointer
    inline cholmod_sparse *getHandle() { return sparse; }
    
    /// Print debugging information
    void print(const char* name = "");

private:
    SparseMatrix(const SparseMatrix& that) = delete; // prevent copy constructor
    void buildLookupIndexFromSparse();
    cholmod_sparse *sparse;
    cholmod_triplet *triplet;
    unsigned int nrow;
    unsigned int ncol;
    double *values;
    std::map<unsigned long, unsigned int> lookupIndex;
    int *iRow;
	int *jColumn;
    inline long key(unsigned int row, unsigned int column){
        int shiftBits = sizeof(long)*8/2; // shift half of the bits of a long
        return (((long)row)<<shiftBits)+column;
    }
    void assertValidIndex(unsigned int row, unsigned int column);
    void assertHasSparse();
    void assertValidInitAddValue(unsigned int row, unsigned int column, double value);
    inline int getIndex(unsigned int row, unsigned int column) {
        assertValidIndex(row, column);
        return lookupIndex[key(row, column)];
    }
    Symmetry symmetry;
#ifdef DEBUG
    unsigned long magicNumber;
    int maxElements;
#endif
};
    
    
    SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS);
    SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS);
    SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS);
    SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS);
    
}

