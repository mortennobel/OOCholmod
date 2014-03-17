//
//  sparse_matrix.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#pragma once

#include <map>
#include <string>
#include <cholmod.h>

#include "config_singleton.h"

namespace oocholmod {
    
    // forward declaration
    class DenseMatrix;
    class Factor;
    
    enum Symmetry {
        SYMMETRIC_LOWER = -1, // Lower triangular part stored
        ASYMMETRIC = 0,
        SYMMETRIC_UPPER = 1, // Upper triangular part stored
    };
    
    enum MatrixState {
        UNINITIALIZED,
        INIT,
        BUILT,
        DESTROYED
    };
    
    class SparseMatrix;
    
    class SparseMatrixIter {
    public:
        SparseMatrixIter(const SparseMatrix* sparseMatrix, int pos);
        
        bool operator!= (const SparseMatrixIter& other) const{
            return pos != other.pos;
        }
        
        inline double &operator* () const;
        
        const SparseMatrixIter& operator++() {
            ++pos;
            return *this;
        }
    private:
        const SparseMatrix* sparseMatrix;
        int pos;
    };
    
    /// The sparse matrix must be used in the following way:
    /// 1. Fill the matrix elements using the (unsigned int row, unsigned int column) function operator
    /// 2. Call build()
    /// 3. Update matrix with elements using the (unsigned int row, unsigned int column) function operator
    ///
	class SparseMatrix {
        friend class Factor;
    public:
        /// nrow # of rows of A
        /// ncol # of columns of A
        /// initialNumberOfElements. If exceeded (during initialization of the matrix) the number of elements will automatically grow with a factor of 1.5
        SparseMatrix(unsigned int nrow = 0, unsigned int ncol = 1, bool symmetric = false, int initialNumberOfElements = 200);
        SparseMatrix(cholmod_sparse *sparse);
        SparseMatrix(SparseMatrix&& move);
        SparseMatrix& operator=(SparseMatrix&& other);
        
        ~SparseMatrix();
        
        SparseMatrixIter begin();
        SparseMatrixIter end();
        
        void symmetrize();
        
        DenseMatrix toDense() const;
        
        MatrixState getMatrixState() const;
        
        // Addition
        friend SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS);
        friend SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS);
        
        // Subtraction
        friend SparseMatrix operator-(const SparseMatrix& M);
        friend SparseMatrix&& operator-(SparseMatrix&& M);
        
        friend SparseMatrix operator-(const SparseMatrix& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator-(SparseMatrix&& LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator-(const SparseMatrix& LHS, SparseMatrix&& RHS);
        friend SparseMatrix&& operator-(SparseMatrix&& LHS, SparseMatrix&& RHS);
        
        // Multiplication
        friend SparseMatrix operator*(const SparseMatrix& LHS, double RHS);
        friend SparseMatrix&& operator*(SparseMatrix&& LHS, double RHS);
        friend SparseMatrix operator*(double LHS, const SparseMatrix& RHS);
        friend SparseMatrix&& operator*(double LHS, SparseMatrix&& RHS);
        
        friend SparseMatrix operator*(const SparseMatrix& LHS, const SparseMatrix& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS);
        friend DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS);
 
        // Print
        friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);
        
        // hasElement only valid on built matrix
        bool hasElement(unsigned int row, unsigned int column) const;

        /// append an another matrix in the init state
        /// useful for parallel matrix assembly
        void append(const SparseMatrix& m);
        
        /// Returns the infinity-norm or 1-norm of a sparse matrix. All xtypes are supported.
        ///  type of norm: 0: inf. norm, 1: 1-norm
        double norm(int norm) const;
        
        ///
        /// Drop small entries from A, and entries in the ignored part of A if A is symmetric.
        /// keep entries with absolute values > tol
        void dropSmallEntries(double tol = 1e-7f);
        
        // in init state return the number of triplets
        // in built state returns the number of elements
        size_t getNumberOfElements();
        
        // Computes Y = alpha*(A*X) + beta*Y (where this == A) or Y = alpha*(transposed(A)*X) + beta*Y when transpose is true
        // transpose is ignores if matrix is symmetric of Hermitian
        void multiply(bool transpose, double alpha, double beta, const DenseMatrix& X, DenseMatrix& Y);
        // Computes result = alpha*(A*X) (where this == A) or Y = alpha*(transposed(A)*X) when transpose is true
        // transpose is ignores if matrix is symmetric of Hermitian
        DenseMatrix multiply(bool transpose, double alpha, const DenseMatrix& X);
        
        // Transpose
        void transpose();
        friend SparseMatrix transposed(const SparseMatrix& M);
        friend SparseMatrix&& transposed(SparseMatrix&& M);
        
        // Solve
        friend DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b);
        friend SparseMatrix solve(const SparseMatrix& A, const SparseMatrix& b);
        friend SparseMatrix solve(const Factor& F, const SparseMatrix& b);

        // Sum the rows and return a vector
        void sumRows(DenseMatrix& outVector);

        // Hard coded method to perform: sparse = spdiags(N)^T * sparse * spdiags(N) - (spdiags(N) - speye())
        void setNullSpace( DenseMatrix& N);
 
        void build();
        
        SparseMatrix copy() const;
        
        Factor analyze() const;
        
        void zero();
        
        void setSymmetry(Symmetry symmetry);
        Symmetry getSymmetry() const { return symmetry; }
        
        unsigned int getRows() const { return nrow; }
        
        unsigned int getColumns() const { return ncol; }
        
        void swap(SparseMatrix& other);
        
        double operator()(unsigned int row, unsigned int column = 0) const;
        
        double& operator()(unsigned int row, unsigned int column = 0);
        
        bool operator==(const SparseMatrix& RHS) const;
        bool operator!=(const SparseMatrix& RHS) const;
        
        friend class SparseMatrixIter;
    private:
        SparseMatrix(const SparseMatrix& that); // prevent copy constructor (no implementation)
        SparseMatrix operator=(const SparseMatrix& other); // prevent copy assignment operator (no implementation)
        long key(unsigned int row, unsigned int column) const;
        void assertValidIndex(unsigned int row, unsigned int column) const;
        void assertHasSparse() const;
        void increaseTripletCapacity();
        void assertValidInitAddValue(unsigned int row, unsigned int column) const;
        
        int binarySearch(int *array, int low, int high, unsigned int value) const;
        
        int getIndex(unsigned int row, unsigned int column) const;
        
        void createTriplet();
        
        double& initAddValue(unsigned int row, unsigned int column);
        
        double& getValue(unsigned int row, unsigned int column);
        
        double getValue(unsigned int row, unsigned int column) const;
        
        cholmod_sparse *sparse;
        cholmod_triplet *triplet;
        unsigned int nrow;
        unsigned int ncol;
        double *values;
        int *iRow;
        int *jColumn;
        Symmetry symmetry;
        int maxTripletElements;
    };
    
    // ----------- inline functions ------------
    
    double& SparseMatrixIter::operator* () const{
        if (sparseMatrix->getMatrixState() == INIT){
            double * ptr = static_cast<double*>(sparseMatrix->triplet->x);
            return *(ptr + pos);
        } else if (sparseMatrix->getMatrixState() == BUILT){
            double * ptr = static_cast<double*>(sparseMatrix->sparse->x);
            return *(ptr + pos);
        } else {
            throw std::runtime_error("Invalid matrix state");
        }
    }
    
    
    
    inline double SparseMatrix::operator()(unsigned int row, unsigned int column) const
    {
        return getValue(row, column);
    }
    
    inline double& SparseMatrix::operator()(unsigned int row, unsigned int column)
    {
        if (sparse != nullptr){
            return getValue(row, column);
        } else {
            return initAddValue(row, column);
        }
    }
    
    inline int SparseMatrix::binarySearch(int *array, int low, int high, unsigned int value) const {
        while (low <= high)
        {
            // http://googleresearch.blogspot.dk/2006/06/extra-extra-read-all-about-it-nearly.html
            int midpoint = (((unsigned int)high + (unsigned int)low) >> 1);
            int midpointValue = array[midpoint];
            if (value == midpointValue) {
                return midpoint;
            } else if (value < midpointValue) {
                high = midpoint - 1;
            } else {
                low = midpoint + 1;
            }
        }
        return -1;
    }
    
    inline int SparseMatrix::getIndex(unsigned int row, unsigned int column) const
    {
#if DEBUG
        assertValidIndex(row, column);
#endif
        if ((symmetry == SYMMETRIC_UPPER && row > column) || (symmetry == SYMMETRIC_LOWER && row < column)) {
            std::swap(row, column);
        }
        
        int iFrom = jColumn[column];
        int iTo = jColumn[column+1]-1;
        
        return binarySearch(iRow, iFrom, iTo, row);
    }
    
    inline long SparseMatrix::key(unsigned int row, unsigned int column) const {
        int shiftBits = sizeof(long)*8/2; // shift half of the bits of a long
        return (((long)row)<<shiftBits)+column;
    }
    
    inline void SparseMatrix::createTriplet(){
        triplet = cholmod_allocate_triplet(nrow, ncol, maxTripletElements, symmetry, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        values = (double *)triplet->x;
        iRow = (int *)triplet->i;
        jColumn = (int *)triplet->j;
    }
    
    inline double& SparseMatrix::initAddValue(unsigned int row, unsigned int column)
    {
        if (!triplet){
            createTriplet();
        } else if (triplet->nnz == maxTripletElements ){
            increaseTripletCapacity();
        }
        assertValidInitAddValue(row, column);
        iRow[triplet->nnz] = row;
        jColumn[triplet->nnz] = column;
        values[triplet->nnz] = 0;
        
        triplet->nnz++;
        return values[triplet->nnz - 1];
    }
    
    inline double& SparseMatrix::getValue(unsigned int row, unsigned int column)
    {
#ifdef DEBUG
        assertHasSparse();
#endif
        int index = getIndex(row, column);
        if (index == -1){
            static double zero = 0;
            zero = 0;
            return zero;
        }
        return values[index];
    }
    
    inline double SparseMatrix::getValue(unsigned int row, unsigned int column) const
    {
#ifdef DEBUG
        assertHasSparse();
#endif
        int index = getIndex(row, column);
        if (index == -1){
            return 0;
        }
        return values[index];
    }
}

