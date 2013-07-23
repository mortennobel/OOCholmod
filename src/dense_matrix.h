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
#include <cholmod.h>

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
        
        DenseMatrix(DenseMatrix&& move);
        
        DenseMatrix& operator=(DenseMatrix&& other);
        
        ~DenseMatrix();
        
        
        inline double& operator()(unsigned int row, unsigned int col = 0)
        {
#ifdef DEBUG
            assert(dense);
            assert(row < nrow && col < ncol);
#endif
            return ((double*)dense->x)[col*nrow + row];
        }
        
        inline double operator()(unsigned int row, unsigned int col = 0) const
        {
#ifdef DEBUG
            assert(dense);
            assert(row < nrow && col < ncol);
#endif
            return ((double*)dense->x)[col*nrow + row];
        }
        
        // OPERATORS
        
        // Addition
        friend DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS);
        friend DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS);
        
        DenseMatrix& operator+=(const DenseMatrix& RHS);
        
        // Multiplication
        friend DenseMatrix operator*(const DenseMatrix& LHS, const double& RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, const double& RHS);
        friend DenseMatrix operator*(const double& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(const double& LHS, DenseMatrix&& RHS);
        
        DenseMatrix& operator*=(const double& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS);
        friend DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS);
        
        /// Returns the infinity-norm, 1-norm, or 2-norm of a dense matrix. Can compute the 2-norm only for a dense column vector. 
        ///  type of norm: 0: inf. norm, 1: 1-norm, 2: 2-norm 
        double norm(int norm) const;
        
        // Transpose
        void transpose();
        friend DenseMatrix transposed(const DenseMatrix& M);
        friend DenseMatrix&& transposed(DenseMatrix&& M);
        
        // Solve
        friend DenseMatrix solve(const DenseMatrix& A, const DenseMatrix& b);
        friend DenseMatrix solve(DenseMatrix&& A, const DenseMatrix& b);
        friend DenseMatrix&& solve(const DenseMatrix& A, DenseMatrix&& b);
        friend DenseMatrix&& solve(DenseMatrix&& A, DenseMatrix&& b);
        
        friend DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b);
        friend DenseMatrix solve(const Factor& F, const DenseMatrix& b);
        
        // Print
        friend std::ostream& operator<<(std::ostream& os, const DenseMatrix& A);
        
        inline double *getData(){ return (double *)(dense->x); };
        inline double *getData() const { return (double *)(dense->x); };
        
        int getRows() const{ return nrow; }
        
        int getColumns() const{ return ncol; }
        
        DenseMatrix copy() const;
        void zero();
        // computes the L^2 norm of the vector
        double length() const;
        
        // elementwise division
        void elemDivide(const DenseMatrix& b);
        void elemDivide(const DenseMatrix& b, DenseMatrix& dest) const;
        
        // elementwise multiplication
        void elemMultiply(const DenseMatrix& b);
        void elemMultiply(const DenseMatrix& b , DenseMatrix& dest) const;
        
        
        double dot(const DenseMatrix& b) const;
        void fill(double value);
        void set(float *data);
        void set(double *data);
        void get(double *outData) const;
        void get(float *outData) const;
        
        void swap(DenseMatrix& other);
    private:
        DenseMatrix(const DenseMatrix& that) = delete; // prevent copy constructor
        DenseMatrix operator=(const DenseMatrix& other) = delete; // prevent copy assignment operator
        cholmod_dense *dense;
        unsigned int nrow;
        unsigned int ncol;
    };
    
    // Addition
    DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS);
    DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS);
    
    // Multiplication
    DenseMatrix operator*(const DenseMatrix& LHS, const double& RHS);
    DenseMatrix&& operator*(DenseMatrix&& LHS, const double& RHS);
    DenseMatrix operator*(const double& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator*(const double& LHS, DenseMatrix&& RHS);
    
    DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS);
    
    // Transpose
    DenseMatrix transposed(const DenseMatrix& M);
    DenseMatrix&& transposed(DenseMatrix&& M);
    
    // Swap
    void swap(DenseMatrix& v1, DenseMatrix& v2);
    
    // Solve
    DenseMatrix solve(const DenseMatrix& A, const DenseMatrix& b);
    DenseMatrix solve(DenseMatrix&& A, const DenseMatrix& b);
    DenseMatrix&& solve(const DenseMatrix& A, DenseMatrix&& b);
    DenseMatrix&& solve(DenseMatrix&& A, DenseMatrix&& b);
    
    // Print
    std::ostream& operator<<(std::ostream& os, const DenseMatrix& A);
}




