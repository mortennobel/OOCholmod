//
//  CholmodFactor.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//  License: LGPL 3.0 

#ifndef __CholmodTest__CholmodFactor__
#define __CholmodTest__CholmodFactor__

#include <iostream>

#include "cholmod.h"

class CholmodSparse; // forward declaration
class CholmodDenseVector;

class CholmodFactor {
public:
    CholmodFactor(cholmod_factor *factor, cholmod_common *Common);
    virtual ~CholmodFactor();
    void factorize(CholmodSparse *sparse);

    cholmod_factor *getFactorHandle() { return factor; };
    
    // solves Ax=b
    void solve(CholmodDenseVector* b, CholmodDenseVector** res);
private:
    CholmodFactor(const CholmodFactor& that); // prevent copy constructor
    cholmod_factor *factor;
    cholmod_common *Common;
#ifdef DEBUG
    long magicNumber;
#endif
};

#endif /* defined(__CholmodTest__CholmodFactor__) */
