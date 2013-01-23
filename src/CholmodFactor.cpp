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

CholmodFactor::CholmodFactor(cholmod_factor *factor, cholmod_common *Common)
:factor(factor), Common(Common) {
}

CholmodFactor::~CholmodFactor(){
    cholmod_free_factor(&factor, Common) ;
}

void CholmodFactor::factorize(CholmodSparse *sparse){
#ifdef DEBUG
    assert(sparse->getSymmetry() != ASYMMETRIC);
#endif
    cholmod_factorize(sparse->getHandle(), factor, Common) ; /* factorize */
}

CholmodDenseVector *CholmodFactor::solve(CholmodDenseVector* b){
    cholmod_dense *x = cholmod_solve(CHOLMOD_A, factor, b->getHandle(), Common);
    return new CholmodDenseVector(x, Common, b->getSize());
}