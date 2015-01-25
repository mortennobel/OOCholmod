//
//  dense_vector.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#pragma once

#include <iostream>
#include <cassert>
#include <cmath>
#include <string> 
#include <cholmod.h>

#ifdef WIN32
#ifndef NAN
	static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
    #define NAN (*(const float *) __nan)
#endif
#endif

#ifndef ndebug_noexcept
#   ifdef DEBUG
#       define ndebug_noexcept
#   else
#       define ndebug_noexcept noexcept
#   endif
#endif


namespace oocholmod {
    
    // forward declaration
    class SparseMatrix;
    class Factor;
    
	class DenseMatrix {
    public:
        // In debug the matrix will be initialized to NAN
        // In release mode, NAN will leave the matrix uninitialized
        DenseMatrix(unsigned int rows = 0, unsigned int cols = 1, double value = NAN);
        
        DenseMatrix(cholmod_dense *x);
        
        DenseMatrix(DenseMatrix&& move) ndebug_noexcept;
        
        DenseMatrix& operator=(DenseMatrix&& other) ndebug_noexcept;
        
        ~DenseMatrix();
        
        
        double& operator()(unsigned int row, unsigned int col = 0) ndebug_noexcept;
        
        double operator()(unsigned int row, unsigned int col = 0) const ndebug_noexcept;
        
        double *begin();
        double *end();
        
        // OPERATORS
        
        
        // Addition
        friend DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS);
        friend DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS);
        
        DenseMatrix& operator+=(const DenseMatrix& RHS);
       
	// Subtraction 
        friend DenseMatrix operator-(const DenseMatrix& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator-(DenseMatrix&& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator-(const DenseMatrix& LHS, DenseMatrix&& RHS);
        friend DenseMatrix&& operator-(DenseMatrix&& LHS, DenseMatrix&& RHS);

        DenseMatrix& operator-=(const DenseMatrix& RHS);
        
        friend DenseMatrix operator-(const DenseMatrix& M);
        friend DenseMatrix&& operator-(DenseMatrix&& M);
 
        // Multiplication
        friend DenseMatrix operator*(const DenseMatrix& LHS, double RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, double RHS);
        friend DenseMatrix operator*(double LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(double LHS, DenseMatrix&& RHS);
        
        DenseMatrix& operator*=(double RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS);
        friend DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS);
        
        /// Returns the infinity-norm, 1-norm, or 2-norm of a dense matrix. Can compute the 2-norm only for a dense column vector. 
        ///  type of norm: 0: inf. norm, 1: 1-norm, 2: 2-norm 
        double norm(int norm) const;
        
        // Returns the determinant. For square matrices ONLY.
        double determinant() const;
        
        // Transpose
        void transpose();
        friend DenseMatrix transposed(const DenseMatrix& M);
        friend DenseMatrix&& transposed(DenseMatrix&& M);
        
        // Inverse
        void inverse();
        friend DenseMatrix inversed(const DenseMatrix& M);
        friend DenseMatrix&& inversed(DenseMatrix&& M);
        
        // Solve
        friend DenseMatrix solve(const DenseMatrix& A, const DenseMatrix& b);
        friend DenseMatrix solve(DenseMatrix&& A, const DenseMatrix& b);
        friend DenseMatrix&& solve(const DenseMatrix& A, DenseMatrix&& b);
        friend DenseMatrix&& solve(DenseMatrix&& A, DenseMatrix&& b);
        
        friend DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b);
        friend DenseMatrix solve(const Factor& F, const DenseMatrix& b);
        
        friend class SparseMatrix;
        
        // Print
        friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& A);
        
        double *getData();
        const double *getData() const;
        
        unsigned int getRows() const;
        
        unsigned int getColumns() const;
        
        DenseMatrix copy() const;
        void zero();
        
        // C = alpha * A + C
        void multiplyAddTo(double alpha, DenseMatrix& C);
        
        // Computes  perform one of the matrix-matrix operations   C := alpha*op( A )*op( B ) + beta*C,
        // this is matrix A
        // op(N) is transposed if transposeN is true
        void multiply(bool transposeA, bool transposeB, double alpha, double beta, const DenseMatrix& B, DenseMatrix& C);
        DenseMatrix multiply(bool transposeA, bool transposeB, double alpha, const DenseMatrix& B);
        
        // assumes that matrix is 1xN or Mx1
        // computes the L^2 norm of the vector
        double length() const;
        
        // elementwise division
        void elemDivide(const DenseMatrix& b);
        void elemDivide(const DenseMatrix& b, DenseMatrix& dest) const;
        
        // elementwise multiplication
        void elemMultiply(const DenseMatrix& b);
        void elemMultiply(const DenseMatrix& b , DenseMatrix& dest) const;
        
        SparseMatrix toSparse() const;
        
        double dot(const DenseMatrix& b) const;
        void fill(double value);
        void set(float *data);
        void set(double *data);
        void get(double *outData) const;
        void get(float *outData) const;
        
        void swap(DenseMatrix& other);
        
        bool operator==(const DenseMatrix& RHS) const;
        bool operator!=(const DenseMatrix& RHS) const;
    private:
        DenseMatrix(const DenseMatrix& that); // prevent copy constructor (no implementation)
        DenseMatrix operator=(const DenseMatrix& other);  // prevent copy assignment operator (no implementation)
        cholmod_dense *dense;
        unsigned int nrow;
        unsigned int ncol;
    };
    
    // ---------- inline functions -----------
    
    inline double& DenseMatrix::operator()(unsigned int row, unsigned int col) ndebug_noexcept
    {
#ifdef DEBUG
        assert(dense);
        assert(row < nrow && col < ncol);
#endif
        return ((double*)dense->x)[col*nrow + row];
    }
    
    inline double DenseMatrix::operator()(unsigned int row, unsigned int col) const ndebug_noexcept
    {
#ifdef DEBUG
        assert(dense);
        assert(row < nrow && col < ncol);
#endif
        return ((double*)dense->x)[col*nrow + row];
    }
    
    inline double *DenseMatrix::getData(){ return (double *)(dense->x); };
    inline const double *DenseMatrix::getData() const { return (double *)(dense->x); };
    
    inline unsigned int DenseMatrix::getRows() const{ return nrow; }
    
    inline unsigned int DenseMatrix::getColumns() const{ return ncol; }
    
}




