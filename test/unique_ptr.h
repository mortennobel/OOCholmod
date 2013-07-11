//
//  unique_ptr.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 7/11/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//
#pragma once
#include <iostream>
#include <vector>
#include <cmath>

#include "CholmodSparse.h"
#include "CholmodFactor.h"
#include "CholmodDenseVector.h"
#include "cpp14.h"

using namespace std;


int TestCaseUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto A = make_unique<CholmodSparse>(3,3,&com);
    
    A->initAddValue(0, 0);
    A->initAddValue(0, 1);
    A->initAddValue(0, 2);
    A->initAddValue(1, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    A->initAddValue(2, 2);
    
    A->build();
    
    // Ax = b
    A->setValue(0, 0, 1);
    A->setValue(0, 1, 1);
    A->setValue(0, 2, 1);
    A->setValue(1, 2, 5);
    A->setValue(2, 2, -1);
    
    auto b = make_unique<CholmodDenseVector>(3, &com);
    (*b)[0] = 6;
    (*b)[1] = -4;
    (*b)[2] = 27;
    b->print("b");
    
    A->print("A");
    auto factor = unique_ptr<CholmodFactor>( A->analyzePtr());
    bool res = factor->factorize(A.get());
    cout << "factor->factorize(A) "<<res<<endl;
    auto x = unique_ptr<CholmodDenseVector>();
    factor->solve(b.get(), x);
    x->print("x");
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, &((*x)[0]), 3);
    
    // update values
    A->zero();
    A->setValue(0, 0, 2);
    A->setValue(0, 1, 9);
    A->setValue(0, 2, 7);
    A->setValue(1, 2, 8);
    A->setValue(2, 2, -3);
    
    factor->factorize(A.get());
    A->print("A");
    
    factor->solve(b.get(), x);
    x->print("x");
    
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, &((*x)[0]), 3);
    return 1;
}

int MultiplyTestUniquePtr(){
    cout << "Multiply test"<<endl;
    cholmod_common com;
    cholmod_start(&com);
    auto A = make_unique<CholmodSparse>(3,3,&com);
    A->initAddValue(0, 0, 1);
    A->initAddValue(0, 1, 1);
    A->initAddValue(0, 2, 1);
    A->initAddValue(1, 2, 5);
    A->initAddValue(2, 2, -1);
    A->build();
    A->print("A");
    
    auto x = make_unique<CholmodDenseVector>(3, &com);
    (*x)[0] = 3;
    (*x)[1] = 7;
    (*x)[2] = 9;
    
    auto res = make_unique<CholmodDenseVector>(3, &com);
    
    A->multiply(x.get(), res.get());
    res->print("b");
    double expected[3] = {19, 48, 29};
    assertEqual(expected, &((*res)[0]), 3);
    return 1;
}

int FillTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto res = make_unique<CholmodDenseVector>(3, &com);
    res->fill(123);
    res->print("Fill 123 test");
    return 1;
}

int DotTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto a = make_unique<CholmodDenseVector>(3, &com);
    auto b = make_unique<CholmodDenseVector>(3, &com);
    (*a)[0] = 1;
    (*a)[1] = 2;
    (*a)[2] = 3;
    
    (*b)[0] = 4;
    (*b)[1] = 5;
    (*b)[2] = 6;
    
    double res = a->dot(b.get());
    double expected = 32;
    cout << "Dot test" << endl;
    assertEqual(&expected, &res, 1);
    return 1;
}

int LengthTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto a = make_unique<CholmodDenseVector>(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    double res = a->length();
    double expected = sqrt(4*4+5*5+6*6);
    cout << "Length test" << endl;
    assertEqual(&expected, &res, 1);
    return 1;
}

int ScaleTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto a = make_unique<CholmodDenseVector>(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    auto b = make_unique<CholmodDenseVector>(3, &com);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->scale(-2);
    cout << "scale test"<<endl;
    assertEqual(&((*a)[0]), &((*b)[0]), 3);
    return 1;
}


int DivideTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto a = make_unique<CholmodDenseVector>(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    auto b = make_unique<CholmodDenseVector>(3, &com);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->divideBy(b.get());
    cout << "divide test"<<endl;
    
    double expected[] = {-0.5,-0.5,-0.5};
    
    assertEqual(expected, &((*a)[0]), 3);
    return 1;
}

int MultiplyVectorTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto a = make_unique<CholmodDenseVector>(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    auto b = make_unique<CholmodDenseVector>(3, &com);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->multiplyWith(b.get());
    cout << "multiply test"<<endl;
    
    double expected[] = {4 * -8,5 * -10, 6 * -12};
    
    assertEqual(expected, &((*a)[0]), 3);
    return 1;
}

int SingularTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto A = make_unique<CholmodSparse>(3,3,&com);
    
    A->initAddValue(2, 2, 1);
    
    A->build();
    
    auto b = make_unique<CholmodDenseVector>(3, &com);
    (*b)[0] = 0;
    (*b)[1] = 1;
    (*b)[2] = 0;
    b->print("b");
    
    A->print("A");
    auto factor = unique_ptr<CholmodFactor>(A->analyzePtr());
    bool res = factor->factorize(A.get());
    cout << "Factorize ok "<< res << endl;
    unique_ptr<CholmodDenseVector> x;
    factor->solve(b.get(), x);
    return 1;
}