//
//  factor.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#pragma once

#include <iostream>

#include <cholmod.h>


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
        // Throws a OOCException if matrix cannot be factorized
        void factorize(const SparseMatrix& sparse);
        
        friend DenseMatrix solve(const Factor& F, const DenseMatrix& b);
        friend SparseMatrix solve(const Factor& F, const SparseMatrix& b);
        
        bool isInitialized();
    private:
        Factor(const Factor& that) {} // prevent copy constructor
        cholmod_factor *factor;
    };
    
	DenseMatrix solve(const Factor& F, const DenseMatrix& b);
	SparseMatrix solve(const Factor& F, const SparseMatrix& b);
}

#include "factor.inc"