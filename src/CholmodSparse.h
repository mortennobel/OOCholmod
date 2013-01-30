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
#include <set>
#include <map>
#include <cmath>

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
    /// maxSize (size allocated before build). 0 means triangular
    CholmodSparse(int nrow, int ncol, cholmod_common *Common, int maxSize = 0);
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
        int shiftBits = sizeof(long)*4;
        long maxId = (long)pow(2, shiftBits);
        assert (row < maxId);
        assert (column < maxId);
#endif
        long k = key(row, column);
        if (usedMap.find(k) != usedMap.end()) return;
        usedMap.insert(k);
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
    std::set<long> usedMap;
    std::map<long, int> lookupIndex;
    int *iRow;
	int *jColumn;
//    int *lookupPosition; // provides fast access to matrix after build (lookupPosition[((column*(column+1))/2 + row])
    inline long key(int row, int column){
        int shiftBits = sizeof(long)*4;
        return (row<<shiftBits)+column;
    }
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
        
#endif
        return lookupIndex[key(row, column)];
    }
    Symmetry symmetry;
#ifdef DEBUG
    
    int maxElements;
#endif
};

#endif /* defined(__CholmodTest__CholmodSparse__) */
