//
//  CholmodFactor.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0 

#include "CholmodFactor.h"
#include "CholmodSparse.h"
#include "CholmodDenseVector.h"
#include "cpp14.h"

using namespace std;

// bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL

CholmodFactor::CholmodFactor(cholmod_factor *factor, cholmod_common *Common)
:factor(factor), Common(Common)
#ifdef DEBUG
,magicNumber(MAGIC_NUMBER)
#endif
{
}

CholmodFactor::CholmodFactor(CholmodFactor&& move)
:factor(move.factor), Common(move.Common)
#ifdef DEBUG
, magicNumber(move.magicNumber)
#endif
{
    move.factor = nullptr;
    move.Common = nullptr;
#ifdef DEBUG
    move.magicNumber = 0;
#endif
}

CholmodFactor& CholmodFactor::operator=(CholmodFactor&& other){
    if (this != &other){
        if (factor != nullptr){
            cholmod_free_factor(&factor, Common) ;
        }
        Common = other.Common;
#ifdef DEBUG
        magicNumber = other.magicNumber;
#endif
        
        other.factor = nullptr;
        other.Common = nullptr;
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
    cholmod_free_factor(&factor, Common) ;
}

bool CholmodFactor::factorize(CholmodSparse& sparse){
    return factorize(&sparse);
}

bool CholmodFactor::factorize(CholmodSparse *sparse){
#ifdef DEBUG
    assert(sparse->getSymmetry() != ASYMMETRIC);
    assert(magicNumber == MAGIC_NUMBER);
#endif
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
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b->getHandle(), Common);
    *res = new CholmodDenseVector(x, Common, b->getSize());
}

void CholmodFactor::solve(CholmodDenseVector* b, std::unique_ptr<CholmodDenseVector> &res){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b->getHandle(), Common);
    res = make_unique<CholmodDenseVector>(x, Common, b->getSize());
}

CholmodDenseVector CholmodFactor::solve(CholmodDenseVector& b){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b.getHandle(), Common);
    return CholmodDenseVector(x, Common, b.getSize());
}