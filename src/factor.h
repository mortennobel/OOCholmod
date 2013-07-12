//
//  factor.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#pragma once

#include <iostream>

#include "cholmod.h"

namespace oocholmod {
    
    
    class SparseMatrix; // forward declaration
    class CholmodDenseVector;
    
    class CholmodFactor {
    public:
        CholmodFactor(cholmod_factor *factor);
        CholmodFactor(CholmodFactor&& move);
        CholmodFactor& operator=(CholmodFactor&& other);
        virtual ~CholmodFactor();
        
        // returns true if factorization is done
        // Return false if matrix is not positive definite
        bool factorize(SparseMatrix *sparse);
        bool factorize(SparseMatrix& sparse);
        
        cholmod_factor *getFactorHandle() { return factor; };
        
        // solves Ax=b
        void solve(CholmodDenseVector* b, CholmodDenseVector** res);
        void solve(CholmodDenseVector* b, std::unique_ptr<CholmodDenseVector> &res);
        CholmodDenseVector solve(CholmodDenseVector& b);
    private:
        CholmodFactor(const CholmodFactor& that) = delete; // prevent copy constructor
        cholmod_factor *factor;
#ifdef DEBUG
        unsigned long magicNumber;
#endif
    };
    
}

