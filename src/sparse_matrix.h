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

#ifndef ndebug_noexcept
#   ifdef DEBUG
#       define ndebug_noexcept
#   else
#       define ndebug_noexcept noexcept
#   endif
#endif

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
        SparseMatrix(SparseMatrix&& move) ndebug_noexcept;
        SparseMatrix& operator=(SparseMatrix&& other) ndebug_noexcept;
        
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

        void absSumRows(DenseMatrix& outVector);

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
        
        double operator()(unsigned int row, unsigned int column = 0) const ndebug_noexcept;
        
        double& operator()(unsigned int row, unsigned int column = 0) ndebug_noexcept;
        
        bool operator==(const SparseMatrix& RHS) const;
        bool operator!=(const SparseMatrix& RHS) const;
        
        friend class SparseMatrixIter;
    private:
        SparseMatrix(const SparseMatrix& that); // prevent copy constructor (no implementation)
        SparseMatrix operator=(const SparseMatrix& other); // prevent copy assignment operator (no implementation)
        long key(unsigned int row, unsigned int column) const ndebug_noexcept;
        void assertValidIndex(unsigned int row, unsigned int column) const;
        void assertHasSparse() const;
        void increaseTripletCapacity() ndebug_noexcept;
        void assertValidInitAddValue(unsigned int row, unsigned int column) const ndebug_noexcept;
        
        int binarySearch(int *array, int low, int high, unsigned int value) const ndebug_noexcept;
        
        int getIndex(unsigned int row, unsigned int column) const ndebug_noexcept;
        
        void createTriplet() ndebug_noexcept;
        
        double& initAddValue(unsigned int row, unsigned int column) ndebug_noexcept;
        
        double& getValue(unsigned int row, unsigned int column) ndebug_noexcept;
        
        double getValue(unsigned int row, unsigned int column) const ndebug_noexcept;
        
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
}

#include "sparse_matrix.inc"