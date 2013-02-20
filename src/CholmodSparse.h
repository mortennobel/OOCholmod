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

// forward declaration
class CholmodDenseVector;

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
class CholmodSparse {
public:
    /// nrow # of rows of A
    /// ncol # of columns of A
    /// maxSize (size allocated before build). 0 means triangular
    CholmodSparse(unsigned int nrow, unsigned int ncol, cholmod_common *Common, int maxSize = 0);
    CholmodSparse(cholmod_sparse *sparse, cholmod_common *Common);
    virtual ~CholmodSparse();
    
    void build(bool readOnly = false);
    
    CholmodFactor *analyze();
    
    void zero();
    
    Symmetry getSymmetry() { return symmetry; }
    
    int getRows(){ return nrow; }
    
    int getColumns(){ return ncol; }
    
    inline void initAddValue(unsigned int row, unsigned int column, double value=0) {
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
        int shiftBits = sizeof(long)*8/2; // shift half of the bits of a long
        long maxId = (long)pow(2, shiftBits);
        assert (row < maxId);
        assert (column < maxId);
#endif
        long k = key(row, column);
        std::map<unsigned long, unsigned int>::iterator res = lookupIndex.find(k);
        if (res != lookupIndex.end()) {
            values[(*res).second] += value;
            return;
        }
        lookupIndex[k] = (int)triplet->nnz;
        iRow[triplet->nnz] = row;
        jColumn[triplet->nnz] = column;
        values[triplet->nnz] = value;
        
        triplet->nnz++;
    }
    
    // computes alpha*(A*X) + beta*Y
    // res is result
    // alpha is optional (default 1)
    // beta is optional (default 0)
    void multiply(CholmodDenseVector *X, CholmodDenseVector *res, double alpha = 1, double beta = 0);
    
    inline double getValue(unsigned int row, unsigned int column){
#ifdef DEBUG
        assert(sparse != NULL); // matrix must be build
#endif
        int index = getIndex(row, column);
        return ((double*)sparse->x)[index];
    }

    inline void addValue(unsigned int row, unsigned int column, double value){
#ifdef DEBUG
        assert(sparse != NULL); // matrix must be build
#endif
        int index = getIndex(row, column);
        ((double*)sparse->x)[index] += value;
    }
    
    
    /// Set value of sparse matrix
    /// Must be invoked after build
    inline void setValue(unsigned int row, unsigned int column, double value){
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
    void buildLookupIndexFromSparse();
    cholmod_common *Common;
    cholmod_sparse *sparse;
    cholmod_triplet *triplet;
    unsigned int nrow;
    unsigned int ncol;
    double *values;
    std::map<unsigned long, unsigned int> lookupIndex;
    int *iRow;
	int *jColumn;
    inline long key(unsigned int row, unsigned int column){
        int shiftBits = sizeof(long)*8/2; // shift half of the bits of a long
        return (((long)row)<<shiftBits)+column;
    }
    inline int getIndex(unsigned int row, unsigned int column) {
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
    unsigned long magicNumber;
    int maxElements;
#endif
};

#endif /* defined(__CholmodTest__CholmodSparse__) */
