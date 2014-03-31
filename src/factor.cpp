//
//  Factor.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include "factor.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "config_singleton.h"
#include "ooc_exception.h"

using namespace std;

namespace oocholmod {
    
    Factor::Factor()
    :factor(nullptr)
    {
    }
    
    Factor::Factor(cholmod_factor *factor)
    :factor(factor)
    {
    }
    
    Factor::Factor(Factor&& move)
    :factor(move.factor)
    {
        move.factor = nullptr;
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

            // clean up
            other.factor = nullptr;
        }
        
        return *this;
    }
    
    Factor::~Factor(){
        if (factor){
            cholmod_free_factor(&factor, ConfigSingleton::getCommonPtr()) ;
        }
    }
    
    void Factor::factorize(const SparseMatrix& A){
#ifdef DEBUG
        assert(A.symmetry != ASYMMETRIC);
        assert(A.sparse);
        assert(factor);
#endif
        auto Common = ConfigSingleton::getCommonPtr();
        cholmod_factorize(A.sparse, factor, Common) ; /* factorize */
        auto status = Common->status;
        Common->status = 0;
        if (status != CHOLMOD_OK){
            OOCException::createOOCException(ConfigSingleton::getLastError());
        }
    }
    
    DenseMatrix solve(const Factor& F, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(F.factor);
        assert(b.dense);
#endif
        cholmod_dense *x = cholmod_solve(CHOLMOD_A, F.factor, b.dense, ConfigSingleton::getCommonPtr());
        return DenseMatrix(x);
    }
    
    SparseMatrix solve(const Factor& F, const SparseMatrix& b)
    {
#ifdef DEBUG
        assert(F.factor);
        assert(b.sparse);
#endif
        cholmod_sparse *x = cholmod_spsolve(CHOLMOD_A, F.factor, b.sparse, ConfigSingleton::getCommonPtr());
        return SparseMatrix(x);
    }
    
}
