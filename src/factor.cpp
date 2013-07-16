//
//  Factor.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include "factor.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "cpp14.h"
#include "config_singleton.h"

using namespace std;

// bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL

namespace oocholmod {
    
    Factor::Factor()
    :factor(nullptr)
#ifdef DEBUG
    ,magicNumber(MAGIC_NUMBER)
#endif
    {
    }
    
    Factor::Factor(cholmod_factor *factor)
    :factor(factor)
#ifdef DEBUG
    ,magicNumber(MAGIC_NUMBER)
#endif
    {
    }
    
    Factor::Factor(Factor&& move)
    :factor(move.factor)
#ifdef DEBUG
    , magicNumber(move.magicNumber)
#endif
    {
        move.factor = nullptr;
#ifdef DEBUG
        move.magicNumber = 0;
#endif
    }
    
    bool Factor::isInitialized(){
        return factor != nullptr;
    }
    
    Factor& Factor::operator=(Factor&& other){
        if (this != &other){
            if (factor != nullptr){
                cholmod_free_factor(&factor, ConfigSingleton::getCommonPtr()) ;
            }
            // copy
            factor = other.factor;
#ifdef DEBUG
            magicNumber = other.magicNumber;
#endif
            // clean up
            other.factor = nullptr;
#ifdef DEBUG
            other.magicNumber = 0;
#endif
        }
        
        return *this;
    }
    
    Factor::~Factor(){
        if (factor){
#ifdef DEBUG
            assert(magicNumber == MAGIC_NUMBER);
            magicNumber = 0;
#endif
            cholmod_free_factor(&factor, ConfigSingleton::getCommonPtr()) ;
        }
    }
    
    bool Factor::factorize(const SparseMatrix& A){
#ifdef DEBUG
        assert(A.symmetry != ASYMMETRIC);
        assert(A.magicNumber == MAGIC_NUMBER);
        assert(A.sparse);
        assert(magicNumber == MAGIC_NUMBER);
        assert(factor);
#endif
        auto Common = ConfigSingleton::getCommonPtr();
        cholmod_factorize(A.sparse, factor, Common) ; /* factorize */
        if (Common->status == CHOLMOD_OK){
            return true;
        }
        Common->status = 0;
        return false;
    }
    
    DenseMatrix solve(const Factor& F, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(F.magicNumber == MAGIC_NUMBER);
        assert(F.factor);
        assert(b.magicNumber == MAGIC_NUMBER);
        assert(b.dense);
#endif
        cholmod_dense *x = cholmod_solve(CHOLMOD_A, F.factor, b.dense, ConfigSingleton::getCommonPtr());
        return DenseMatrix(x);
    }
    
    SparseMatrix solve(const Factor& F, const SparseMatrix& b)
    {
#ifdef DEBUG
        assert(F.magicNumber == MAGIC_NUMBER);
        assert(F.factor);
        assert(b.magicNumber == MAGIC_NUMBER);
        assert(b.sparse);
#endif
        cholmod_sparse *x = cholmod_spsolve(CHOLMOD_A, F.factor, b.sparse, ConfigSingleton::getCommonPtr());
        return SparseMatrix(x);
    }
    
}
