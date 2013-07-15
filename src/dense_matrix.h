//
//  dense_vector.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#pragma once

#include <iostream>
#include <cassert>


#include <cholmod.h>

namespace oocholmod {
    
    // forward declaration
    class SparseMatrix;
    
    class DenseMatrix {
    public:
        DenseMatrix(unsigned int rows, unsigned int cols = 1);
        
        DenseMatrix(cholmod_dense *x);
        
        DenseMatrix(DenseMatrix&& move);
        
        DenseMatrix& operator=(DenseMatrix&& other);
        
        ~DenseMatrix();
        
        // OPERATORS
        
        // Addition
        friend DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS);
        friend DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS);
        
        // Multiplication
        friend DenseMatrix operator*(const DenseMatrix& LHS, const double& RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, const double& RHS);
        friend DenseMatrix operator*(const double& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(const double& LHS, DenseMatrix&& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(const DenseMatrix& LHS, DenseMatrix&& RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, DenseMatrix&& RHS);
        
        friend DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS);
        friend DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, const SparseMatrix& RHS);
        friend DenseMatrix operator*(SparseMatrix&& LHS, const DenseMatrix& RHS);
        friend DenseMatrix&& operator*(const SparseMatrix& LHS, DenseMatrix&& RHS);
        friend DenseMatrix operator*(const DenseMatrix& LHS, SparseMatrix&& RHS);
        friend DenseMatrix&& operator*(SparseMatrix&& LHS, DenseMatrix&& RHS);
        friend DenseMatrix&& operator*(DenseMatrix&& LHS, SparseMatrix&& RHS);
        
        // Transpose
        void transpose();
        friend DenseMatrix transposed(const DenseMatrix& M);
        friend DenseMatrix&& transposed(DenseMatrix&& M);
        
        
        inline double *getData(){ return (double *)(dense->x); };
        inline double *getData() const { return (double *)(dense->x); };
        
        int getRows() const{ return nrow; }
        
        int getColumns() const{ return ncol; }
        
        DenseMatrix copy() const;
        void zero();
        // computes the L^2 norm of the vector
        double length();
        
        void scale(double alpha);
        
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
        inline void add(unsigned int row, unsigned int col, double value) {
#ifdef DEBUG
            assert(row < nrow && col < ncol);
#endif
            getData()[col*nrow + row] += value;
        };
        
        inline double operator [](int i) const    {return getData()[i];}
        inline double & operator [](int i) {return getData()[i];}
        
        void swap(DenseMatrix& other);
        inline cholmod_dense *getHandle() const { return dense; }
        void print(const char* name = "") const;
    private:
        DenseMatrix(const DenseMatrix& that) = delete; // prevent copy constructor
        DenseMatrix operator=(const DenseMatrix& other) = delete; // prevent copy assignment operator
        cholmod_dense *dense;
        unsigned int nrow;
        unsigned int ncol;
#ifdef DEBUG
        unsigned long magicNumber;
#endif
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
    DenseMatrix&& operator*(DenseMatrix&& LHS, const DenseMatrix& RHS);
    DenseMatrix&& operator*(const DenseMatrix& LHS, DenseMatrix&& RHS);
    DenseMatrix&& operator*(DenseMatrix&& LHS, DenseMatrix&& RHS);
    
    // Transpose
    DenseMatrix transposed(const DenseMatrix& M);
    DenseMatrix&& transposed(DenseMatrix&& M);
    
    void swap(DenseMatrix& v1, DenseMatrix& v2);
}




