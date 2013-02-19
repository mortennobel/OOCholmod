//
//  CholmodFactor.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//  License: LGPL 3.0 

#include "CholmodFactor.h"
#include "CholmodSparse.h"
#include "CholmodDenseVector.h"

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

CholmodFactor::~CholmodFactor(){
#ifdef DEBUG
    magicNumber = 0;
#endif
    cholmod_free_factor(&factor, Common) ;
}

void CholmodFactor::factorize(CholmodSparse *sparse){
#ifdef DEBUG
    assert(sparse->getSymmetry() != ASYMMETRIC);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_factorize(sparse->getHandle(), factor, Common) ; /* factorize */
}

void CholmodFactor::solve(CholmodDenseVector* b, CholmodDenseVector **res){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    delete *res;
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b->getHandle(), Common);
    *res = new CholmodDenseVector(x, Common, b->getSize());
}