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
#include <algorithm>
#include <climits>

#include "sparse_matrix.h"
#include "factor.h"
#include "dense_matrix.h"


using namespace std;
using namespace oocholmod;

int BuildSparseTestObj()
{
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 1) = -1;
    A(2, 0) = -3;
    A(0, 0) = 3;
    A.build();
    
    assert(A(2,1) == 4);
    assert(A(1,2) == 4);
    assert(A(0,2) == -2);
    assert(A(2,0) == -2);
    assert(A(0,0) == 4);
    
    SparseMatrix B{3,3};
    B(0, 0) = 1;
    B(0, 1) = 1;
    B(0, 2) = 1;
    B(1, 2) = 5;
    B(2, 1) = -1;
    B(2, 0) = -3;
    B(0, 0) = 3;
    B.build();
    
    assert(B(2,1) == -1);
    assert(B(1,2) == 5);
    assert(B(0,2) == 1);
    assert(B(2,0) == -3);
    assert(B(0,0) == 4);
    return 1;
}

int EqualSparseTestObj()
{
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1.356;
    A(1, 2) = 5;
    A(2, 1) = -1.12;
    A(2, 0) = -3;
    A.build();
    
    SparseMatrix B{3,3, true};
    B(0, 0) = 1;
    B(0, 1) = 1;
    B(0, 2) = 1.356;
    B(1, 2) = 5;
    B(2, 1) = -1.12;
    B(2, 0) = -3;
    B.build();
    assert(A == B);
    return 1;
}

int SolveDenseDenseTestObj()
{
    DenseMatrix A{3,3, 1.};
    A(0, 0) = -1;
    A(0, 1) = 5;
    A(0, 2) = -3;
    A(1, 2) = 5;
    A(2, 0) = -1;
    
    DenseMatrix b{3, 2, 1.};
    b(0) = 6;
    b(1) = -4;
    b(2) = 27;
    
    // Ax = b
    DenseMatrix x = solve(A, b);
    double expected[6] = {-23.8750, -1.0625, 4.1875, -0.5, 0.25, 0.25};
    assertEqual(expected, x.getData(), 6);
    
    DenseMatrix x2 = solve(std::move(A), std::move(b));
    assertEqual(expected, x2.getData(), 6);
    return 1;
}

int SolveSparseDenseTestObj()
{    
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    DenseMatrix b{3};
    b(0) = 6;
    b(1) = -4;
    b(2) = 27;
    
    // Ax = b
    DenseMatrix x = solve(A, b);
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, x.getData(), 3);
    return 1;
}

int SolveSparseSparseTestObj()
{
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    SparseMatrix b{3,2};
    b(0,0) = 6;
    b(0,1) = -4;
    b(2,1) = 27;
    b(1,1) = -2;
    b.build();
    
    // Ax = b
    SparseMatrix x = solve(A, b);
    double expected[6] = {10.7143, -2.5714, -2.1429, -15.9286, 9.1429, 2.7857};
    assert(std::abs(expected[0] - x(0,0)) < 0.0001);
    assert(std::abs(expected[2] - x(2,0)) < 0.0001);
    return 1;
}

int SolveSparseDenseFactorTestObj()
{    
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    DenseMatrix b{3};
    b(0) = 6;
    b(1) = -4;
    b(2) = 27;
    
    Factor F;
    F = A.analyze();
    bool res = F.factorize(A);
    DenseMatrix x = solve(F, b);
    double expected[] = {2.78571f,4.57143f,-1.35714f};
    assertEqual(expected, x.getData(), 3);
    
    // update values
    A.zero();
    A(0, 0) = 2;
    A(0, 1) = 9;
    A(0, 2) = 7;
    A(1, 2) = 8;
    A(2, 2) = -3;
    
    F.factorize(A);
    x = solve(F, b);
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, x.getData(), 3);
    
    return 1;
}

int AddSparseSparseTestObj()
{
    SparseMatrix A(3,3, true);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    SparseMatrix B(3,3, true);
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

int AddDenseDenseTestObj()
{
    DenseMatrix A(3, 3, 0.);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    
    DenseMatrix B(3, 3, 1.);
    B(0, 0) = 2;
    B(0, 1) = -1;
    B(0, 2) = 3;
    B(1, 2) = -5;
    B(2, 2) = -1;
    
    DenseMatrix C = A + B;
    DenseMatrix D = (A+B) + B;
    DenseMatrix E = std::move(A) + (B+A);
    
    double expected[9] = {4, 1, 1, 1, 1, 1, 5, 5, -3};
    assertEqual(expected, E.getData(), 9);
    
    return 1;
}

int AddEqualDenseDenseTestObj()
{
    DenseMatrix A(3, 3, 0.);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    
    DenseMatrix B(3, 3, 1.);
    B(0, 0) = 2;
    B(0, 1) = -1;
    B(0, 2) = 3;
    B(1, 2) = -5;
    B(2, 2) = -1;
    
    A += B+A;
    
    double expected[9] = {4, 1, 1, 1, 1, 1, 5, 5, -3};
    assertEqual(expected, A.getData(), 9);
    
    return 1;
}

int TransposeDenseTestObj()
{
    DenseMatrix b{3};
    b(0) = 6;
    b(1) = -4;
    b(2) = 27;
    
    b.transpose();
    
    DenseMatrix c = transposed(b);
    
    DenseMatrix d = transposed(std::move(c));
    
    assertEqual(b.getData(), d.getData(), 3)
    return 1;
}

int TestCaseFunctionOperatorObj(){
    
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    DenseMatrix b{3};
    b(0) = 6;
    b(1) = -4;
    b(2) = 27;
    
    // Ax = b
    Factor factor = A.analyze();
    bool res = factor.factorize(A);
    TINYTEST_ASSERT(res);
    
    DenseMatrix x = solve(factor, b);
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
    x = solve(factor, b);
    
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, x.getData(), 3);
    
    return 1;
}

int MultiplySparseSparseTestObj(){
    SparseMatrix A(2,3);
    A(0, 0) = 1;
    A(0, 1) = 2;
    A(0, 2) = 1;
    A(1, 1) = -2;
    A(1, 2) = 5;
    A.build();
    
    SparseMatrix B(3,3, true);
    B(0, 0) = 2;
    B(0, 1) = -1;
    B(0, 2) = 3;
    B(1, 1) = 6;
    B(1, 2) = -5;
    B(2, 2) = -1;
    B.build();
    
    SparseMatrix C = A*B;
    assert(C(1,1) == -37);
    assert(C(0,2) == -8);
    
    SparseMatrix D;
    D = (C*B)*transposed(A);

    return 1;
}

int MultiplyScalarSparseTestObj(){
    SparseMatrix A(3,3, true);
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();

    SparseMatrix B = A*0.1;
    SparseMatrix C = -10.*(A*B);
    A = 0.3 * A;
    
    assert(C(2,2) == -27);
    assert(C(0,2) == -5);
    
    return 1;
}


int MultiplySparseDenseTestObj(){
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 5;
    A(2, 2) = -1;
    A.build();
    
    DenseMatrix x{3,2, 1.};
    x(0,1) = 3;
    x(1,1) = 7;
    x(2,0) = 9;
    
    DenseMatrix y{1,3};
    y(0,0) = 3;
    y(0,1) = 7;
    y(0,2) = 9;
    
    DenseMatrix res = A*x;
    double expected[6] = {11, 46, -3, 11, 8, 37};
    assertEqual(expected, res.getData(), 6);
    
    DenseMatrix res2 = y*A;
    double expected2[3] = {19, 48, 29};
    assertEqual(expected2, res2.getData(), 3);
    return 1;
}

int MultiplyScalarDenseTestObj()
{
    DenseMatrix x{3,2,1.};
    x(0,1) = -3;
    x(1,1) = 7;
    x(2,0) = 9;

    DenseMatrix res1 = -3.*x;
    
    DenseMatrix res2 = 0.1*std::move(res1);
    
    double expected2[6] = {-0.3, -0.3, -2.7, 0.9, -2.1, -0.3};
    assertEqual(expected2, res2.getData(), 6);
    
    return 1;
}

int MultiplyEqualScalarDenseTestObj()
{
    DenseMatrix x{3,2,1.};
    x(0,1) = -3;
    x(1,1) = 7;
    x(2,0) = 9;
    
    x *= -3.;
    x *= 0.1;
    
    double expected2[6] = {-0.3, -0.3, -2.7, 0.9, -2.1, -0.3};
    assertEqual(expected2, x.getData(), 6);
    
    return 1;
}

int MultiplyDenseDenseTestObj(){
    
    DenseMatrix x{3, 2, 1.};
    x(0,1) = 3;
    x(1,1) = 7;
    x(2,0) = 9;
    
    DenseMatrix y{2, 3, 1.};
    y(0,0) = -3;
    y(0,1) = 7;
    y(1,1) = 2;
    
    DenseMatrix res = y*x;
    double expected[4] = {13, 12, 41, 18};
    assertEqual(expected, res.getData(), 3);
    
    DenseMatrix v{3,1};
    v(0) = 3;
    v(1) = 7;
    v(2) = 9;
    
    DenseMatrix w{1,3};
    w(0,0) = -3;
    w(0,1) = 7;
    w(0,2) = 2;
    
    DenseMatrix res2 = w*v;
    assert(58 == res2(0,0));
    
    DenseMatrix res3 = v*w;
    assert(-9 == res3(0,0));
    
    DenseMatrix res4 = w*x;
    double expected2[2] = {22, 42};
    assertEqual(expected2, res4.getData(), 2);
    return 1;
}

int FillTestObj(){
    DenseMatrix res(3);
    res.fill(123);
    double expected[3] = {123, 123, 123};
    assertEqual(expected, res.getData(), 3);
    
    DenseMatrix res2(3, 1, 123);
    assertEqual(expected, res2.getData(), 3);
    return 1;
}

int DotTestObj(){
    DenseMatrix a{3};
    DenseMatrix b{3};
    a(0) = 1;
    a(1) = 2;
    a(2) = 3;
    
    b(0) = 4;
    b(1) = 5;
    b(2) = 6;
    
    double res = a.dot(b);
    double expected = 32;
    assertEqual(&expected, &res, 1);
    return 1;
}

int LengthTestObj(){
    DenseMatrix a{3};
    a(0) = 4;
    a(1) = 5;
    a(2) = 6;
    
    double res = a.length();
    double expected = sqrt(4*4+5*5+6*6);
    assertEqual(&expected, &res, 1);
    return 1;
}

int ElemDivideTestObj(){
    DenseMatrix a{3};
    a(0) = 4;
    a(1) = 5;
    a(2) = 6;
    
    DenseMatrix b{3};
    b(0) = -8;
    b(1) = -10;
    b(2) = -12;
    a.elemDivide(b);
    
    double expected[] = {-0.5,-0.5,-0.5};
    
    assertEqual(expected, a.getData(), 3);
    return 1;
}

int ElemMultiplyTestObj(){
    DenseMatrix a{3};
    a(0) = 4;
    a(1) = 5;
    a(2) = 6;
    
    DenseMatrix b{3};
    b(0) = -8;
    b(1) = -10;
    b(2) = -12;
    a.elemMultiply(b);
    double expected[] = {4 * -8,5 * -10, 6 * -12};
    
    assertEqual(expected, a.getData(), 3);
    return 1;
}

int SingularTestObj(){
    SparseMatrix A{3,3, true};
    
    A(2, 2) = 1;
    
    A.build();
    
    DenseMatrix b{3};
    b(0) = 0;
    b(1) = 1;
    b(2) = 0;
    
    Factor factor = A.analyze();
    bool res = factor.factorize(A);
    TINYTEST_ASSERT(!res);
    DenseMatrix x = solve(factor, b);
    return 1;
}

int SwapTest(){
    DenseMatrix AD{1,1};
    AD(0) = 1;
    DenseMatrix BD{1,1};
    BD(0) = 2;
    swap(AD, BD);
    TINYTEST_ASSERT(AD(0) == 2);
    TINYTEST_ASSERT(BD(0) == 1);
    
    SparseMatrix A{1,1};
    A(0,0) = 1;
    A.build();
    SparseMatrix B{1,1};
    B(0,0) = 2;
    B.build();
    std::swap(A, B);
    TINYTEST_ASSERT(A(0,0) == 2);
    TINYTEST_ASSERT(B(0,0) == 1);
    
    return 1;
}

int SparseStateTest(){
    SparseMatrix A{1,1};
    TINYTEST_ASSERT(A.getMatrixState() == UNINITIALIZED);
    A(0,0) = 1;
    TINYTEST_ASSERT(A.getMatrixState() == INIT);
    A.build();
    TINYTEST_ASSERT(A.getMatrixState() == BUILT);
    SparseMatrix B = move(A);
    TINYTEST_ASSERT(A.getMatrixState() == DESTROYED);
    return 1;
}

int DynamicTripletGrow(){
    int size = 1000;
    SparseMatrix A{size,size,true, 1};
    for (int i=0;i<size;i++){
        A(i,i) = 123;
    }
    A.build();
    for (int x=0;x<size;x++) for (int y=x;y<size;y++){
        double expected = x==y?123:0;
        double val = A(x,y);
        TINYTEST_ASSERT(expected == val);
    }
    return 1;
}

int LargeSparseMatrix(){
    int size = UINT_MAX/100;
    SparseMatrix A{size,size,true, 1};
    A(0,0) = 1;
    A(size-1,size-1) = 2;
    A.build();

    TINYTEST_ASSERT(A(0,0) == 1);
    TINYTEST_ASSERT(A(size-1,size-1) == 2);
    return 1;
}

int IndexTest2(){
    int size = 3;
    SparseMatrix C{size,size,false};
    int value = 1;
    for (int x=0;x<size;x++) for (int y=0;y<size;y++){
        C(x,y) = value;
        value++;
    }
    C.build();
    value = 1;
    for (int x=0;x<size;x++) for (int y=0;y<size;y++){
        TINYTEST_ASSERT(C(x,y) == value);
        value++;
    }
    return 1;
}

int IndexTest(){
    int size = 3;
    SparseMatrix A{size,size,true};
    DenseMatrix B{size,size};
    SparseMatrix C{size,size,false};
    int value = 1;
    for (int x=0;x<size;x++) for (int y=x;y<size;y++){
        A(x,y) = value;
        B(x,y) = value;
        B(y,x) = value;
        C(x,y) = value;
        if(x != y) {
            C(y,x) = value;
        }
        value++;
    }
    A.build();
    C.build();
    for (int x=0;x<size;x++) for (int y=x;y<size;y++){
        TINYTEST_ASSERT(A(x,y) == A(y,x));
        TINYTEST_ASSERT(B(x,y) == B(y,x));
        TINYTEST_ASSERT(C(x,y) == C(y,x));
        TINYTEST_ASSERT(A(x,y) == B(y,x));
        TINYTEST_ASSERT(A(x,y) == C(y,x));
    }
    
    return 1;
}

