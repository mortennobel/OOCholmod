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
    class DenseMatrix;
    
    class Factor {
        friend class SparseMatrix;
        Factor(cholmod_factor *factor);
    public:
        Factor();
        Factor(Factor&& move);
        Factor& operator=(Factor&& other);
        virtual ~Factor();
        
        // returns true if factorization is done
        // Return false if matrix is not positive definite
        bool factorize(SparseMatrix& sparse);
        
        friend DenseMatrix solve(Factor& F, DenseMatrix& b);
        
        bool isInitialized();
    private:
        Factor(const Factor& that) = delete; // prevent copy constructor
        cholmod_factor *factor;
#ifdef DEBUG
        unsigned long magicNumber;
#endif
    };
    
    DenseMatrix solve(Factor& F, DenseMatrix& b);
}

