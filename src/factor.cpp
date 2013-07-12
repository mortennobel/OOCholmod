//
//  CholmodFactor.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0 

#include "factor.h"
#include "sparse_matrix.h"
#include "dense_vector.h"
#include "cpp14.h"
#include "config_singleton.h"

using namespace std;

// bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL

namespace oocholmod {

CholmodFactor::CholmodFactor(cholmod_factor *factor)
:factor(factor)
#ifdef DEBUG
,magicNumber(MAGIC_NUMBER)
#endif
{
}

CholmodFactor::CholmodFactor(CholmodFactor&& move)
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

CholmodFactor& CholmodFactor::operator=(CholmodFactor&& other){
    if (this != &other){
        if (factor != nullptr){
            cholmod_free_factor(&factor, ConfigSingleton::getCommonPtr()) ;
        }
#ifdef DEBUG
        magicNumber = other.magicNumber;
#endif
        
        other.factor = nullptr;
#ifdef DEBUG
        magicNumber = 0;
#endif
    }
    
    return *this;
}

CholmodFactor::~CholmodFactor(){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
    magicNumber = 0;
#endif
    cholmod_free_factor(&factor, ConfigSingleton::getCommonPtr()) ;
}

bool CholmodFactor::factorize(SparseMatrix& sparse){
    return factorize(&sparse);
}

bool CholmodFactor::factorize(SparseMatrix *sparse){
#ifdef DEBUG
    assert(sparse->getSymmetry() != ASYMMETRIC);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    auto Common = ConfigSingleton::getCommonPtr();
    cholmod_factorize(sparse->getHandle(), factor, Common) ; /* factorize */
    if (Common->status == CHOLMOD_OK){
        return true;
    }
    Common->status = 0;
    return false;
}

void CholmodFactor::solve(CholmodDenseVector* b, CholmodDenseVector **res){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    delete *res;
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b->getHandle(), ConfigSingleton::getCommonPtr());
    *res = new CholmodDenseVector(x, b->getSize());
}

void CholmodFactor::solve(CholmodDenseVector* b, std::unique_ptr<CholmodDenseVector> &res){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b->getHandle(), ConfigSingleton::getCommonPtr());
    res = make_unique<CholmodDenseVector>(x, b->getSize());
}

CholmodDenseVector CholmodFactor::solve(CholmodDenseVector& b){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b.getHandle(), ConfigSingleton::getCommonPtr());
    return CholmodDenseVector(x, b.getSize());
}
    
}
