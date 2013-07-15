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
    
    class DenseMatrix {
    public:
        DenseMatrix(unsigned int size);
        
        DenseMatrix(cholmod_dense *x, unsigned int size);
        
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
        
        // Transpose
        friend DenseMatrix transposed(const DenseMatrix& M);
        friend DenseMatrix&& transposed(DenseMatrix&& M);
        
        
        inline double *getData(){ return (double *)(x->x); };
        inline double *getData() const { return (double *)(x->x); };
        
        int getRows() const{ return nrow; }
        
        int getColumns() const{ return ncol; }
        
        int getSize() const { return size; }
        
        void copyTo(DenseMatrix *dest);
        void copyTo(DenseMatrix& dest);
        void zero();
        // computes the L^2 norm of the vector
        double length();
        void scale(double alpha);
        // elementwise division
        void divideBy(const DenseMatrix& b);
        // elementwise multiplication
        void multiplyWith(const DenseMatrix& b);
        double dot(const DenseMatrix& b);
        void fill(double value);
        void set(float *data);
        void set(double *data);
        void get(double *outData);
        void get(float *outData);
        inline void add(unsigned int index, double value) {
#ifdef DEBUG
            assert(index < size);
#endif
            getData()[index] += value;
        };
        
        inline double operator [](int i) const    {return getData()[i];}
        inline double & operator [](int i) {return getData()[i];}
        inline cholmod_dense *getHandle() { return x; }
        void print(const char* name = "");
    private:
        DenseMatrix(const DenseMatrix& that) = delete; // prevent copy constructor
        DenseMatrix operator=(const DenseMatrix& other) = delete; // prevent copy assignment operator
        cholmod_dense *x;
        unsigned int size;
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
    
}

