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
    
    class DenseVector {
    public:
        DenseVector(unsigned int size);
        DenseVector(cholmod_dense *x, unsigned int size);
        DenseVector(DenseVector&& move);
        DenseVector& operator=(DenseVector&& other);
        virtual ~DenseVector();
        inline double *getData(){ return (double *)(x->x); };
        inline double *getData() const { return (double *)(x->x); };
        inline int getSize() const { return size; }
        void copyTo(DenseVector *dest);
        void copyTo(DenseVector& dest);
        void zero();
        // computes the L^2 norm of the vector
        double length();
        void scale(double alpha);
        // elementwise division
        void divideBy(const DenseVector& b);
        // elementwise multiplication
        void multiplyWith(const DenseVector& b);
        double dot(const DenseVector& b);
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
        void print(const char* name);
    private:
        DenseVector(const DenseVector& that) = delete; // prevent copy constructor
        cholmod_dense *x;
        unsigned int size;
#ifdef DEBUG
        unsigned long magicNumber;
#endif
    };
    
}

