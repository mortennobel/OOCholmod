//
//  objects.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 7/11/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#pragma once

#include <iostream>
#include <vector>
#include <cmath>

#include "sparse_matrix.h"
#include "factor.h"
#include "dense_matrix.h"
#include "cpp14.h"


using namespace std;
using namespace oocholmod;

int TestCaseObj(){
    
    SparseMatrix A{3,3};
    
    A(0, 0) = 0;
    A(0, 1) = 0;
    A(0, 2) = 0;
    A(1, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    
    A.build();
    
    // Ax = b
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    
    DenseMatrix b{3};
    b[0] = 6;
    b[1] = -4;
    b[2] = 27;
    //b->print("b");
    
    //A->print("A");
    Factor factor = A.analyze();
    bool res = factor.factorize(A);
    //cout << "factor->factorize(A) "<<res<<endl;
    DenseMatrix x = factor.solve(b);
    //x->print("x");
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, x.getData(), 3);
    
    // update values
    A.zero();
    A(0, 0) = 2;
    A(0, 1) = 9;
    A(0, 2) = 7;
    A(1, 2) = 8;
    A(2, 2) = -3;
    
    factor.factorize(A);
    //A->print("A");
    
    x = factor.solve(b);
    //x->print("x");
    
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, x.getData(), 3);
    
    return 1;
}

int AddTestObj()
{
    SparseMatrix A(3,3);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    SparseMatrix B(3,3);
    B(0, 0) = 2;
    B(0, 1) = -1;
    B(0, 2) = 3;
    B(1, 2) = -5;
    B(2, 2) = -1;
    B.build();
    
    SparseMatrix C = A + B;
    
    SparseMatrix D;
    D = (A+B) + B;
    
    SparseMatrix E = std::move(A) + (B+A);
    
    assert(E(2,2) == -3);
    assert(E(0,2) == 5);
    
    return 1;
}

int TransposeTestObj()
{
    DenseMatrix b{3};
    b[0] = 6;
    b[1] = -4;
    b[2] = 27;
    
    b.transpose();
    
    DenseMatrix c = transposed(b);
    
    DenseMatrix d = transposed(std::move(c));
    
    assertEqual(b.getData(), d.getData(), 3)
    
    // Cannot do before asymmetric matrices are supported.
    return 1;
}

int TestCaseFunctionOperatorObj(){
    
    SparseMatrix A{3,3};
    
    A(0, 0) = 0;
    A(0, 1) = 0;
    A(0, 2) = 0;
    A(1, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    A(2, 2) = 0;
    
    A.build();
    
    // Ax = b
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    
    DenseMatrix b{3};
    b[0] = 6;
    b[1] = -4;
    b[2] = 27;
    //b->print("b");
    
    //A->print("A");
    Factor factor = A.analyze();
    bool res = factor.factorize(A);
    TINYTEST_ASSERT(res);
    //cout << "factor->factorize(A) "<<res<<endl;
    DenseMatrix x = factor.solve(b);
    //x->print("x");
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, x.getData(), 3);
    
    // update values
    A.zero();
    A(0, 0) = 2;
    A(0, 1) = 9;
    A(0, 2) = 7;
    A(1, 2) = 8;
    A(2, 2) = -3;
    
    factor.factorize(A);
    //A->print("A");
    
    x = factor.solve(b);
    //x->print("x");
    
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, x.getData(), 3);
    
    return 1;
}

int MultiplyMatrixMatrixTestObj(){
    SparseMatrix A(3,3);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    SparseMatrix B(3,3);
    B(0, 0) = 2;
    B(0, 1) = -1;
    B(0, 2) = 3;
    B(1, 1) = 6;
    B(1, 2) = -5;
    B(2, 2) = -1;
    B.build();
    
    SparseMatrix C = A*B;
    
    SparseMatrix D = B*(A*B);
//    assert(D(2,2) == 22);
//    assert(D(0,2) == -67);
    
    SparseMatrix E;
    E = (A*B)*C;
    
    return 1;
}

int MultiplyScalarMatrixTestObj(){
    SparseMatrix A(3,3);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    SparseMatrix B = A*0.1;
    SparseMatrix C = -13.4*(A*B);
    
    A = 0.3 * A;
//    A.print();
//    assert(C(2,2) == 281.4);
//    assert(C(0,2) == 40.2);
    
    return 1;
}


int MultiplyTestObj(){
    SparseMatrix A{3,3};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    DenseMatrix x{3};
    x[0] = 3;
    x[1] = 7;
    x[2] = 9;
    
    DenseMatrix res = A*x;
    double expected[3] = {19, 48, 29};
    assertEqual(expected, res.getData(), 3);
    
    x.transpose();
    DenseMatrix res2 = x*A;
    assertEqual(expected, res2.getData(), 3);
    return 1;
}

int FillTestObj(){
    DenseMatrix res(3);
    res.fill(123);
    //res->print("Fill 123 test");
    double expected[3] = {123, 123, 123};
    assertEqual(expected, res.getData(), 3);
    return 1;
}

int DotTestObj(){
    DenseMatrix a{3};
    DenseMatrix b{3};
    a[0] = 1;
    a[1] = 2;
    a[2] = 3;
    
    b[0] = 4;
    b[1] = 5;
    b[2] = 6;
    
    double res = a.dot(b);
    double expected = 32;
    assertEqual(&expected, &res, 1);
    return 1;
}

int LengthTestObj(){
    DenseMatrix a{3};
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    double res = a.length();
    double expected = sqrt(4*4+5*5+6*6);
    assertEqual(&expected, &res, 1);
    return 1;
}

int ScaleTestObj(){
    DenseMatrix a{3};
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    DenseMatrix b{3};
    b[0] = -8;
    b[1] = -10;
    b[2] = -12;
    a.scale(-2);
    assertEqual(a.getData(), b.getData(), 3);
    return 1;
}


int DivideTestObj(){
    DenseMatrix a{3};
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    DenseMatrix b{3};
    b[0] = -8;
    b[1] = -10;
    b[2] = -12;
    a.elemDivide(b);
    
    double expected[] = {-0.5,-0.5,-0.5};
    
    assertEqual(expected, a.getData(), 3);
    return 1;
}

int MultiplyVectorTestObj(){
    DenseMatrix a{3};
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    DenseMatrix b{3};
    b[0] = -8;
    b[1] = -10;
    b[2] = -12;
    a.elemMultiply(b);
    double expected[] = {4 * -8,5 * -10, 6 * -12};
    
    assertEqual(expected, a.getData(), 3);
    return 1;
}

int SingularTestObj(){
    SparseMatrix A{3,3};
    
    A(2, 2) = 1;
    
    A.build();
    
    DenseMatrix b{3};
    b[0] = 0;
    b[1] = 1;
    b[2] = 0;
    //b->print("b");
    
    //A->print("A");
    Factor factor = A.analyze();
    bool res = factor.factorize(A);
    TINYTEST_ASSERT(!res);
    //cout << "Factorize ok "<< res << endl;
    DenseMatrix x = factor.solve(b);
    return 1;
}