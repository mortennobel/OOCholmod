//
//  lapack.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#ifndef __OOCholmod__lapack__
#define __OOCholmod__lapack__

#include <iostream>

#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int 		__CLPK_integer;
typedef int 		__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef int 		__CLPK_ftnlen;
#else
typedef long int 	__CLPK_integer;
typedef long int 	__CLPK_logical;
typedef float 		__CLPK_real;
typedef double 		__CLPK_doublereal;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef long int 	__CLPK_ftnlen;
#endif


int dgesv_(__CLPK_integer *n, __CLPK_integer *nrhs, __CLPK_doublereal *a, __CLPK_integer
           *lda, __CLPK_integer *ipiv, __CLPK_doublereal *b, __CLPK_integer *ldb, __CLPK_integer *info);

#endif /* defined(__OOCholmod__lapack__) */
