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

#include "sparse_matrix.h"
#include "factor.h"
#include "dense_vector.h"
#include "cpp14.h"

using namespace std;
using namespace oocholmod;


int TestCaseUniquePtr(){
    auto A = make_unique<SparseMatrix>(3,3);
    
    (*A)(0, 0);
    (*A)(0, 1);
    (*A)(0, 2);
    (*A)(1, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    (*A)(2, 2);
    
    A->build();
    
    // Ax = b
    (*A)(0, 0) = 1;
    (*A)(0, 1) = 1;
    (*A)(0, 2) = 1;
    (*A)(1, 2) = 5;
    (*A)(2, 2) = -1;
    
    auto b = make_unique<CholmodDenseVector>(3);
    (*b)[0] = 6;
    (*b)[1] = -4;
    (*b)[2] = 27;
    //b->print("b");
    
    //A->print("A");
    auto factor = unique_ptr<CholmodFactor>( A->analyzePtr());
    bool res = factor->factorize(A.get());
    //cout << "factor->factorize(A) "<<res<<endl;
    auto x = unique_ptr<CholmodDenseVector>();
    factor->solve(b.get(), x);
    //x->print("x");
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, &((*x)[0]), 3);
    
    // update values
    A->zero();
    (*A)(0, 0) = 2;
    (*A)(0, 1) = 9;
    (*A)(0, 2) = 7;
    (*A)(1, 2) = 8;
    (*A)(2, 2) = -3;
    
    factor->factorize(A.get());
    //A->print("A");
    
    factor->solve(b.get(), x);
    //x->print("x");
    
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, &((*x)[0]), 3);
    return 1;
}

int MultiplyTestUniquePtr(){
    auto A = make_unique<SparseMatrix>(3,3);
    (*A)(0, 0) = 1;
    (*A)(0, 1) = 1;
    (*A)(0, 2) = 1;
    (*A)(1, 2) = 5;
    (*A)(2, 2) = -1;
    A->build();
    //A->print("A");
    
    auto x = make_unique<CholmodDenseVector>(3);
    (*x)[0] = 3;
    (*x)[1] = 7;
    (*x)[2] = 9;
    
    auto res = make_unique<CholmodDenseVector>(3);
    
    A->multiply(x.get(), res.get());
    //res->print("b");
    double expected[3] = {19, 48, 29};
    assertEqual(expected, &((*res)[0]), 3);
    return 1;
}

int FillTestUniquePtr(){
    auto res = make_unique<CholmodDenseVector>(3);
    res->fill(123);
    double expected[3] = {123, 123, 123};
    assertEqual(expected, res->getData(), 3);
    return 1;
}

int DotTestUniquePtr(){
    auto a = make_unique<CholmodDenseVector>(3);
    auto b = make_unique<CholmodDenseVector>(3);
    (*a)[0] = 1;
    (*a)[1] = 2;
    (*a)[2] = 3;
    
    (*b)[0] = 4;
    (*b)[1] = 5;
    (*b)[2] = 6;
    
    double res = a->dot(b.get());
    double expected = 32;
    //cout << "Dot test" << endl;
    assertEqual(&expected, &res, 1);
    return 1;
}

int LengthTestUniquePtr(){
    auto a = make_unique<CholmodDenseVector>(3);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    double res = a->length();
    double expected = sqrt(4*4+5*5+6*6);
    //cout << "Length test" << endl;
    assertEqual(&expected, &res, 1);
    return 1;
}

int ScaleTestUniquePtr(){
    auto a = make_unique<CholmodDenseVector>(3);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    auto b = make_unique<CholmodDenseVector>(3);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->scale(-2);
    //cout << "scale test"<<endl;
    assertEqual(&((*a)[0]), &((*b)[0]), 3);
    return 1;
}


int DivideTestUniquePtr(){
    auto a = make_unique<CholmodDenseVector>(3);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    auto b = make_unique<CholmodDenseVector>(3);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->divideBy(b.get());
    //cout << "divide test"<<endl;
    
    double expected[] = {-0.5,-0.5,-0.5};
    
    assertEqual(expected, &((*a)[0]), 3);
    return 1;
}

int MultiplyVectorTestUniquePtr(){
    auto a = make_unique<CholmodDenseVector>(3);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    auto b = make_unique<CholmodDenseVector>(3);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->multiplyWith(b.get());
    //cout << "multiply test"<<endl;
    
    double expected[] = {4 * -8,5 * -10, 6 * -12};
    
    assertEqual(expected, &((*a)[0]), 3);
    return 1;
}

int SingularTestUniquePtr(){
    auto A = make_unique<SparseMatrix>(3,3);
    
    (*A)(2, 2) = 1;
    
    A->build();
    
    auto b = make_unique<CholmodDenseVector>(3);
    (*b)[0] = 0;
    (*b)[1] = 1;
    (*b)[2] = 0;
    //b->print("b");
    
    //A->print("A");
    auto factor = unique_ptr<CholmodFactor>(A->analyzePtr());
    bool res = factor->factorize(A.get());
    TINYTEST_ASSERT(!res);
    unique_ptr<CholmodDenseVector> x;
    factor->solve(b.get(), x);
    return 1;
}