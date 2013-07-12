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
#include "dense_vector.h"
#include "cpp14.h"


using namespace std;
using namespace oocholmod;

int TestCaseObj(){
    
    SparseMatrix A(3,3);
    
    A.initAddValue(0, 0);
    A.initAddValue(0, 1);
    A.initAddValue(0, 2);
    A.initAddValue(1, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    A.initAddValue(2, 2);
    
    A.build();
    
    // Ax = b
    A.setValue(0, 0, 1);
    A.setValue(0, 1, 1);
    A.setValue(0, 2, 1);
    A.setValue(1, 2, 5);
    A.setValue(2, 2, -1);
    
    CholmodDenseVector b(3);
    b[0] = 6;
    b[1] = -4;
    b[2] = 27;
    //b->print("b");
    
    //A->print("A");
    CholmodFactor factor = A.analyze();
    bool res = factor.factorize(A);
    //cout << "factor->factorize(A) "<<res<<endl;
    CholmodDenseVector x = factor.solve(b);
    //x->print("x");
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, x.getData(), 3);
    
    // update values
    A.zero();
    A.setValue(0, 0, 2);
    A.setValue(0, 1, 9);
    A.setValue(0, 2, 7);
    A.setValue(1, 2, 8);
    A.setValue(2, 2, -3);
    
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
    A.initAddValue(0, 0, 1);
    A.initAddValue(0, 1, 1);
    A.initAddValue(0, 2, 1);
    A.initAddValue(1, 2, 5);
    A.initAddValue(2, 2, -1);
    A.build();
    
    SparseMatrix B(3,3);
    B.initAddValue(0, 0, 2);
    B.initAddValue(0, 1, -1);
    B.initAddValue(0, 2, 3);
    B.initAddValue(1, 2, -5);
    B.initAddValue(2, 2, -1);
    B.build();
    
    SparseMatrix res = A + B;
    res.print("res");
    return 1;
}

int MultiplyTestObj(){
    SparseMatrix A(3,3);
    A.initAddValue(0, 0, 1);
    A.initAddValue(0, 1, 1);
    A.initAddValue(0, 2, 1);
    A.initAddValue(1, 2, 5);
    A.initAddValue(2, 2, -1);
    A.build();
    //A->print("A");
    
    CholmodDenseVector x(3);
    x[0] = 3;
    x[1] = 7;
    x[2] = 9;
    
    CholmodDenseVector res = A.multiply(x);
    //res->print("b");
    double expected[3] = {19, 48, 29};
    assertEqual(expected, res.getData(), 3);
    return 1;
}

int FillTestObj(){
    CholmodDenseVector res(3);
    res.fill(123);
    //res->print("Fill 123 test");
    double expected[3] = {123, 123, 123};
    assertEqual(expected, res.getData(), 3);
    return 1;
}

int DotTestObj(){
    CholmodDenseVector a(3);
    CholmodDenseVector b(3);
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
    CholmodDenseVector a(3);
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    double res = a.length();
    double expected = sqrt(4*4+5*5+6*6);
    assertEqual(&expected, &res, 1);
    return 1;
}

int ScaleTestObj(){
    CholmodDenseVector a(3);
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    CholmodDenseVector b(3);
    b[0] = -8;
    b[1] = -10;
    b[2] = -12;
    a.scale(-2);
    assertEqual(a.getData(), b.getData(), 3);
    return 1;
}


int DivideTestObj(){
    CholmodDenseVector a(3);
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    CholmodDenseVector b(3);
    b[0] = -8;
    b[1] = -10;
    b[2] = -12;
    a.divideBy(b);
    
    double expected[] = {-0.5,-0.5,-0.5};
    
    assertEqual(expected, a.getData(), 3);
    return 1;
}

int MultiplyVectorTestObj(){
    CholmodDenseVector a(3);
    a[0] = 4;
    a[1] = 5;
    a[2] = 6;
    
    CholmodDenseVector b(3);
    b[0] = -8;
    b[1] = -10;
    b[2] = -12;
    a.multiplyWith(b);
    
    double expected[] = {4 * -8,5 * -10, 6 * -12};
    
    assertEqual(expected, a.getData(), 3);
    return 1;
}

int SingularTestObj(){
    SparseMatrix A(3,3);
    
    A.initAddValue(2, 2, 1);
    
    A.build();
    
    CholmodDenseVector b(3);
    b[0] = 0;
    b[1] = 1;
    b[2] = 0;
    //b->print("b");
    
    //A->print("A");
    CholmodFactor factor = A.analyze();
    bool res = factor.factorize(A);
    TINYTEST_ASSERT(!res);
    //cout << "Factorize ok "<< res << endl;
    CholmodDenseVector x = factor.solve(b);
    return 1;
}