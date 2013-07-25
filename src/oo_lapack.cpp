//
//  lapack.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#include "oo_lapack.h"
#include <cassert>

#ifdef NO_LAPACK

int dgesv_(__CLPK_integer *n, __CLPK_integer *nrhs, __CLPK_doublereal *a, __CLPK_integer
           *lda, __CLPK_integer *ipiv, __CLPK_doublereal *b, __CLPK_integer *ldb, __CLPK_integer *info){
    assert(false); // not implemented
}

#endif