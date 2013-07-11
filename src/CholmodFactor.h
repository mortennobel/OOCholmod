//
//  CholmodFactor.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
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
    CholmodFactor(CholmodFactor&& move);
    CholmodFactor& operator=(CholmodFactor&& other);
    virtual ~CholmodFactor();
    
    // returns true if factorization is done
    // Return false if matrix is not positive definite
    bool factorize(CholmodSparse *sparse);
    bool factorize(CholmodSparse& sparse);

    cholmod_factor *getFactorHandle() { return factor; };
    
    // solves Ax=b
    void solve(CholmodDenseVector* b, CholmodDenseVector** res);
    void solve(CholmodDenseVector* b, std::unique_ptr<CholmodDenseVector> &res);
    CholmodDenseVector solve(CholmodDenseVector& b);
private:
    CholmodFactor(const CholmodFactor& that); // prevent copy constructor
    cholmod_factor *factor;
    cholmod_common *Common;
#ifdef DEBUG
    unsigned long magicNumber;
#endif
};

#endif /* defined(__CholmodTest__CholmodFactor__) */
