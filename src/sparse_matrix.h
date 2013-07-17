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
#include "config_singleton.h"

namespace oocholmod {
    
    // forward declaration
    class DenseMatrix;
    
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
        friend class Factor;
    public:
        /// nrow # of rows of A
        /// ncol # of columns of A
        /// maxSize (size allocated before build). 0 means triangular
        SparseMatrix(unsigned int nrow = 0, unsigned int ncol = 1, bool symmetric = false, int maxSize = 0);
        SparseMatrix(cholmod_sparse *sparse);
        SparseMatrix(SparseMatrix&& move);
        SparseMatrix& operator=(SparseMatrix&& other);
        
        virtual ~SparseMatrix();
        
        // Addition
        friend SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS);
        friend SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS);
        
        // Multiplication
        friend SparseMatrix operator*(const SparseMatrix& LHS, const double& RHS);
        friend SparseMatrix&& operator*(SparseMatrix&& LHS, const double& RHS);
        friend SparseMatrix operator*(const double& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator*(const double& LHS, SparseMatrix&& RHS);
        
        friend SparseMatrix operator*(const SparseMatrix& LHS, const SparseMatrix& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS);
        friend DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS);
        
        // Transpose
        void transpose();
        friend SparseMatrix transposed(const SparseMatrix& M);
        friend SparseMatrix&& transposed(SparseMatrix&& M);
        
        // Solve
        friend DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b);
        friend SparseMatrix solve(const SparseMatrix& A, const SparseMatrix& b);
        friend SparseMatrix solve(const Factor& F, const SparseMatrix& b);
        
        void build(bool readOnly = false);
        
        Factor analyze() const;
        
        void zero();
        
        void setSymmetry(Symmetry symmetry);
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
        void setNullSpace(const DenseMatrix& N);
        
        /// Print debugging information
        void print(const char* name = "");
        
        void swap(SparseMatrix& other);
        
        double& operator()(unsigned int row, unsigned int column = 0) {
            if (sparse != nullptr){
                return getValue(row, column);
            } else {
                return initAddValue(row, column);
            }
            
        }
    private:
        SparseMatrix(const SparseMatrix& that) = delete; // prevent copy constructor
        SparseMatrix operator=(const SparseMatrix& other) = delete; // prevent copy assignment operator
        void buildLookupIndexFromSparse();
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
        
        double& initAddValue(unsigned int row, unsigned int column, double value=0)
        {
            if (!triplet){
                triplet = cholmod_allocate_triplet(nrow, ncol, maxTripletElements, symmetry, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
                values = (double *)triplet->x;
                iRow = (int *)triplet->i;
                jColumn = (int *)triplet->j;
            }
            assertValidInitAddValue(row, column, value);
            long k = key(row, column);
            auto res = lookupIndex.find(k);
            if (res != lookupIndex.end()) {
                auto& pos = (*res).second;
                values[pos] += value;
                return values[pos];
            }
            lookupIndex[k] = (int)triplet->nnz;
            iRow[triplet->nnz] = row;
            jColumn[triplet->nnz] = column;
            values[triplet->nnz] = value;
            
            triplet->nnz++;
            return values[triplet->nnz - 1];
        }
        
        double& getValue(unsigned int row, unsigned int column){
            assertHasSparse();
            int index = getIndex(row, column);
            return ((double*)sparse->x)[index];
        }
        
        cholmod_sparse *sparse;
        cholmod_triplet *triplet;
        unsigned int nrow;
        unsigned int ncol;
        double *values;
        std::map<unsigned long, unsigned int> lookupIndex;
        int *iRow;
        int *jColumn;
        Symmetry symmetry;
        int maxTripletElements;
    };
    
    // Addition
    SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS);
    SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS);
    SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS);
    SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS);
    
    // Multiplication
    SparseMatrix operator*(const SparseMatrix& LHS, const double& RHS);
    SparseMatrix&& operator*(SparseMatrix&& LHS, const double& RHS);
    SparseMatrix operator*(const double& LHS, const SparseMatrix& RHS);
    SparseMatrix&& operator*(const double& LHS, SparseMatrix&& RHS);
    
    SparseMatrix operator*(const SparseMatrix& LHS, const SparseMatrix& RHS);
    
    // DenseMatrix times SparseMatrix (Note that SparseMatrix times DenseMatrix may be faster).
    DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS);
    DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS);
    
    // Transpose
    SparseMatrix transposed(const SparseMatrix& M);
    SparseMatrix&& transposed(SparseMatrix&& M);
    
    // Swap
    void swap(SparseMatrix& v1, SparseMatrix& v2);
    
    // Solve
    DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b);
    SparseMatrix solve(const SparseMatrix& A, const SparseMatrix& b);
}

