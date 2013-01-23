//
//  CholmodSparse.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//  License: LGPL 3.0 

#ifndef __CholmodTest__CholmodSparse__
#define __CholmodTest__CholmodSparse__

#include <iostream>
#include <cassert>

#include "cholmod.h"
#include "CholmodFactor.h"

enum Symmetry {
    SYMMETRIC_LOWER = -1, // Lower triangular part stored
    ASYMMETRIC = 0,
    SYMMETRIC_UPPER = 1, // Upper triangular part stored
};


///
/// Currently only real (double), upper symmetric matrices are supported.
///
/// The sparse matrix must be used in the following way:
/// 1. Fill the matrix elements using the mark() method
/// 2. Call build()
/// 3. Fill matrix with elements using setValue or addValue
///
/// At any point after the matrix has been build, you can call 
class CholmodSparse {
public:
    /// nrow # of rows of A
    /// ncol # of columns of A
    CholmodSparse(int nrow, int ncol, cholmod_common *Common);
    ~CholmodSparse();
    
    void build();
    CholmodFactor *analyze();
    
    void zero();
    
    Symmetry getSymmetry() { return symmetry; }
    
    inline void mark(int row, int column) {
#ifdef DEBUG        
        assert(sparse == NULL); // must be called before matrix build
        assert(triplet->nnz < maxElements);
        assert(row < nrow);
        assert(column < ncol);
        if (symmetry == SYMMETRIC_UPPER) {
            assert(row <= column);
        } else if (symmetry == SYMMETRIC_LOWER) {
            assert(row >= column);
        }
#endif
        for (int i = 0; i < triplet->nnz; i++) {
            if (iRow[i] == row && jColumn[i] == column){
                return; // already marked skip
            }
        }
        iRow[triplet->nnz] = row;
        jColumn[triplet->nnz] = column;
#ifdef DEBUG
        values[triplet->nnz] = triplet->nnz + 1; // easy to tract value
#else
        values[triplet->nnz] = 1.0;
#endif
        
        triplet->nnz++;
    }
    
    inline double getValue(int row, int column){
#ifdef DEBUG
        assert(sparse != NULL); // matrix must be build
#endif
        int index = getIndex(row, column);
        std::cout << "index "<<index<< std::endl;
        return ((double*)sparse->x)[index];
    }

    inline void addValue(int row, int column, double value){
#ifdef DEBUG
        assert(sparse != NULL); // matrix must be build
#endif
        int index = getIndex(row, column);
        ((double*)sparse->x)[index] += value;
    }
    
    
    /// Set value of sparse matrix
    /// Must be invoked after build
    inline void setValue(int row, int column, double value){
#ifdef DEBUG
        assert(sparse != NULL); // matrix must be build
#endif
        int index = getIndex(row, column);
        ((double*)sparse->x)[index] += value;
    }
    
    
    /// Get cholmod_sparse pointer
    inline cholmod_sparse *getHandle() { return sparse; }
    
    /// Print debugging information
    void print(const char* name);

private:
    CholmodSparse(const CholmodSparse& that); // prevent copy constructor
    cholmod_common *Common;
    cholmod_sparse *sparse;
    cholmod_triplet *triplet;
    int nrow;
    int ncol;
    double *values;
    int *iRow;
	int *jColumn;
    int *lookupPosition; // provides fast access to matrix after build (lookupPosition[((column*(column+1))/2 + row])
    inline int getIndex(int row, int column) {
#ifdef DEBUG
        assert(symmetry == SYMMETRIC_UPPER);
        assert(row < nrow);
        assert(column < ncol);
        if (symmetry == SYMMETRIC_UPPER) {
            assert(row <= column);
        } else if (symmetry == SYMMETRIC_LOWER) {
            assert(row >= column);
        }
        assert(lookupPosition != NULL);
#endif
        // Fast way to compute the column offset http://en.wikipedia.org/wiki/Triangular_number
        int columnOffset = getTriangularNumber(column);
        int position = lookupPosition[columnOffset + row];
        assert(position>=0);
        return position;
    }
    
    inline int getTriangularNumber(int i) {
        return (i*(i+1))/2;
    }
    Symmetry symmetry;
#ifdef DEBUG
    
    int maxElements;
#endif
};

#endif /* defined(__CholmodTest__CholmodSparse__) */
