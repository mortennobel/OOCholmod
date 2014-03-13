//
//  oocholmodtest.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 28/10/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#include "oocholmodtest.h"
#include "tinytest.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <climits>
#include <sstream>
#include "sparse_matrix.h"
#include "factor.h"
#include "dense_matrix.h"
#include "timer.h"

using namespace std;
using namespace oocholmod;

#define assertEqual(E, A, L) { \
const double *expected_ = (E); \
const double *actual_ = (A); \
int length_ = (L); \
for (int i=0;i<length_;i++){ \
if (fabs(expected_[i] - actual_[i]) > 0.01){ \
stringstream str_;\
str_ << "Vectors are not equal. Element "; \
str_ << i; \
str_ << " was "; \
str_ << actual_[i]; \
str_ << " but expected "; \
str_ << expected_[i];\
TINYTEST_ASSERT_MSG(false, str_.str().c_str()); \
} \
if (std::isnan(expected_[i]) ^ std::isnan(actual_[i])){ \
TINYTEST_ASSERT_MSG(false, "NaN"); \
} \
} \
}


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
    
    TINYTEST_ASSERT(A(2,1) == 4);
    TINYTEST_ASSERT(A(1,2) == 4);
    TINYTEST_ASSERT(A(0,2) == -2);
    TINYTEST_ASSERT(A(2,0) == -2);
    TINYTEST_ASSERT(A(0,0) == 4);
    
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
    
    SparseMatrix C{3,3};
    C.build();
    
    assert(C(0,0) == 0);
    assert(C(0,1) == 0);
    
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
    
    DenseMatrix x2 = solve(move(A), move(b));
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
    
    SparseMatrix E = move(A) + (B+A);
    
    assert(E(2,2) == -3);
    assert(E(0,2) == 5);
    
    return 1;
}

int SubtractSparseSparseTestObj()
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
    
    SparseMatrix C = -A - B;
    
    SparseMatrix D;
    D = (A-B) - B;
    
    SparseMatrix E = -A - (A - move(B));
    assert(E(2,2) == 1);
    assert(E(0,2) == 1);
    
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
    DenseMatrix E = move(A) + (B+A);
    
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
    
    DenseMatrix d = transposed(move(c));
    
    assertEqual(b.getData(), d.getData(), 3)
    return 1;
}

int InverseDenseTestObj()
{
    DenseMatrix A{3,3};
    A(0,0) = 10;
    A(0,1) = 3;
    A(0,2) = 1;
    A(1,0) = -1;
    A(1,1) = 4;
    A(1,2) = 1;
    A(2,0) = 5.1;
    A(2,1) = -2.1;
    A(2,2) = 4.5;
    
    DenseMatrix B = inversed(A);
    
    DenseMatrix Ai(3,3);
    Ai(0,0) = 0.0950;
    Ai(0,1) = -0.0738;
    Ai(0,2) = -0.0047;
    Ai(1,0) = 0.0454;
    Ai(1,1) = 0.1887;
    Ai(1,2) = -0.0520;
    Ai(2,0) = -0.0865;
    Ai(2,1) = 0.1716;
    Ai(2,2) = 0.2033;
    
    assertEqual(B.getData(), Ai.getData(), 9)
    
    B.inverse();
    
    assertEqual(B.getData(), A.getData(), 9)
    
    DenseMatrix C = inversed(move(B));
    
    assertEqual(C.getData(), Ai.getData(), 9)
    
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
    
    DenseMatrix res2 = 0.1*move(res1);
    
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
    
    TINYTEST_ASSERT(58 == res2(0,0));
    
    DenseMatrix res3 = v*w;
    TINYTEST_ASSERT(-9 == res3(0,0));
    
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
    swap(A, B);
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
    unsigned int size = 1000;
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
    unsigned int size = UINT_MAX/100;
    SparseMatrix A{size,size,true, 1};
    A(0,0) = 1;
    A(size-1,size-1) = 2;
    A.build();
    
    TINYTEST_ASSERT(A(0,0) == 1);
    TINYTEST_ASSERT(A(size-1,size-1) == 2);
    return 1;
}

int IndexTest2(){
    unsigned int size = 3;
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
    unsigned int size = 3;
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

int LargeMatrixPerformance(){
    Timer timer;
    
    unsigned int size = 30000;
    int testSize = 5000;
    SparseMatrix A{size,size};
    timer.start();
    int offsetX = 0;
    int offsetY = 0;
    for (int x=0;x<testSize;x++) for (int y=0;y<testSize;y++){
        // using non linear access to simulate real usage
        offsetX+=23;
        offsetX+=47;
        A((x+offsetX)%size, (y+offsetY)%size) = x+y;
    }
    timer.stop();
    cout <<"Initializing "<<timer.getElapsedTimeInMilliSec()<<endl;
    timer.start();
    A.build();
    timer.stop();
    cout <<"Building "<<timer.getElapsedTimeInMilliSec()<<endl;
    
    timer.start();
    // build index
    double x = A(0,0);
    timer.stop();
    cout <<"Indexing "<<timer.getElapsedTimeInMilliSec()<<endl;
    
    // set new values
    timer.start();
    offsetX = 0;
    offsetY = 0;
    for (int x=0;x<testSize;x++) for (int y=0;y<testSize;y++){
        // using non linear access to simulate real usage
        offsetX+=23;
        offsetX+=47;
        A((x+offsetX)%size, (y+offsetY)%size) = x+y;
    }
    timer.stop();
    cout <<"Updating "<<timer.getElapsedTimeInMilliSec()<<endl;
    
    return 1;
}

int DropSmallEntriesTest(){
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 0.5;
    A(2, 2) = -0.5;
    A.build();
    A.dropSmallEntries(0.5);
    TINYTEST_ASSERT(A.hasElement(0,0));
    TINYTEST_ASSERT(A.hasElement(0,1));
    TINYTEST_ASSERT(A.hasElement(0,2));
    TINYTEST_ASSERT(!A.hasElement(1,2));
    TINYTEST_ASSERT(!A.hasElement(2,2));
    
    return 1;
}

int CopyTest(){
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 0.5;
    A(2, 2) = -0.5;
    A.build();
    SparseMatrix B = A.copy();
    
    TINYTEST_ASSERT(A == B);
    return 1;
}

int NormTest(){
    SparseMatrix A{3,3, true};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 0.5;
    A(2, 2) = -0.5;
    A.build();
    double norm = A.norm(1);
    TINYTEST_ASSERT(norm == 3);
    
    
    DenseMatrix B{3,3,0};
    B(0, 0) = 1;
    B(2, 2) = 1;
    B(0, 1) = 1;
    B(1, 0) = 1;
    B(0, 2) = 1;
    B(2, 0) = 1;
    B(1, 2) = 0.5;
    B(2, 1) = 0.5;
    B(2, 2) = -0.5;
    norm = B.norm(1);
    TINYTEST_ASSERT(norm == 3);
    return 1;
}

int DeterminantTest()
{
    DenseMatrix A{2,2,0};
    A(0, 0) = -7;
    A(1, 1) = 11;
    A(0, 1) = -1;
    A(1, 0) = 1.6;
    TINYTEST_ASSERT(A.determinant() == A(0,0)*A(1,1) - A(1,0)*A(0,1));
    
    DenseMatrix B{3,3,0};
    B(0, 0) = 1;
    B(2, 2) = -7.7;
    B(0, 1) = -1;
    B(1, 0) = 1;
    B(0, 2) = 10.5;
    B(2, 0) = 5;
    B(1, 2) = 0.5;
    B(2, 1) = 0.5;
    B(2, 2) = -0.5;
    TINYTEST_ASSERT(B.determinant() == B(0,0)*B(1,1)*B(2,2) + B(0,1)*B(1,2)*B(2,0) + B(0,2)*B(1,0)*B(2,1) - B(0,0)*B(1,2)*B(2,1) - B(0,1)*B(1,0)*B(2,2) - B(0,2)*B(1,1)*B(2,0));
    return 1;
}

int AppendTest(){
    vector<SparseMatrix> AS;
    for(int i=0;i<5;i++){
        AS.push_back(SparseMatrix{3,3});
    }
    AS[0](0, 0) = 1;
    AS[1](0, 1) = 1;
    AS[2](0, 2) = 1;
    AS[3](1, 2) = 0.5;
    AS[4](2, 2) = -0.5;
    SparseMatrix A{3,3};
    for (SparseMatrix &a : AS){
        A.append(a);
    }
    A.build();
    TINYTEST_ASSERT(A.hasElement(0,0));
    TINYTEST_ASSERT(A.hasElement(0,1));
    TINYTEST_ASSERT(!A.hasElement(1,0));
    TINYTEST_ASSERT(A.hasElement(0,2));
    TINYTEST_ASSERT(!A.hasElement(2,0));
    TINYTEST_ASSERT(A.hasElement(1,2));
    TINYTEST_ASSERT(!A.hasElement(2,1));
    TINYTEST_ASSERT(A.hasElement(2,2));
    return 1;
}

int NumberOfElementsTest(){
    SparseMatrix A{3,3};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 0.5;
    A(2, 2) = -0.5;
    A(2, 2) = 0; // add extra element
    TINYTEST_ASSERT(A.getNumberOfElements()==6);
    A.build();
    TINYTEST_ASSERT(A.getNumberOfElements()==5);
    return 1;
}

int ZeroTest(){
    SparseMatrix A{3,3};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 0.5;
    A(2, 2) = -0.5;
    A.build();
    A.zero();
    TINYTEST_ASSERT(A(0, 0) == 0);
    TINYTEST_ASSERT(A(0, 1) == 0);
    TINYTEST_ASSERT(A(0, 2) == 0);
    TINYTEST_ASSERT(A(1, 2) == 0);
    TINYTEST_ASSERT(A(2, 2) == 0);
    TINYTEST_ASSERT(A(2, 2) == 0);
    
    DenseMatrix B{3,3};
    B(0, 0) = 1;
    B(0, 1) = 1;
    B(0, 2) = 1;
    B(1, 2) = 0.5;
    B(2, 2) = -0.5;
    B.zero();
    TINYTEST_ASSERT(B(0, 0) == 0);
    TINYTEST_ASSERT(B(0, 1) == 0);
    TINYTEST_ASSERT(B(0, 2) == 0);
    TINYTEST_ASSERT(B(1, 2) == 0);
    TINYTEST_ASSERT(B(2, 2) == 0);
    TINYTEST_ASSERT(B(2, 2) == 0);
    
    return 1;
}

int DenseSetGetTest(){
    DenseMatrix A{3,3};
    float floats[9] = {1,2,3,4,5,6,7,8,9};
    double doubles[9] = {2,3,4,5,6,7,8,9,1};
    float floatRes[9];
    double doubleRes[9];
    A.set(floats);
    A.get(floatRes);
    A.set(doubles);
    A.get(doubleRes);
    for (int i=0;i<9;i++){
        TINYTEST_ASSERT(floats[i] == floatRes[i]);
        TINYTEST_ASSERT(doubles[i] == doubleRes[i]);
    }
    return 1;
}

int SparseToDense(){
    SparseMatrix A{3,3};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(0, 2) = 1;
    A(1, 2) = 0.5;
    A(2, 2) = -0.5;
    A.build();
    
    DenseMatrix B{3,3,0};
    B(0, 0) = 1;
    B(0, 1) = 1;
    B(0, 2) = 1;
    B(1, 2) = 0.5;
    B(2, 2) = -0.5;
    TINYTEST_ASSERT(A.toDense() == B);
    
    TINYTEST_ASSERT(A.toDense().toSparse() == A);
    
    B(2, 2) = 999;
    TINYTEST_ASSERT(A.toDense() != B);
    return 1;
}

int DenseSubstractionTest(){
    DenseMatrix B{2,2,0};
    B(0, 0) = 1;
    B(1, 0) = 2;
    B(0, 1) = 3;
    B(1, 1) = 4;
    
    DenseMatrix B2{2,2,0};
    B2(0, 0) = 5;
    B2(1, 0) = 6;
    B2(0, 1) = 7;
    B2(1, 1) = 8;
    
    DenseMatrix R = B2-B;
    for (int i=0;i<4;i++){
        TINYTEST_ASSERT(R.getData()[i] == 4);
    }
    
    R = (-B) + B2;
    for (int i=0;i<4;i++){
        TINYTEST_ASSERT(R.getData()[i] == 4);
    }
    
    R = B2.copy();
    
    R-=B;
    for (int i=0;i<4;i++){
        TINYTEST_ASSERT(R.getData()[i] == 4);
    }
    
    R = move(B2)-B;
    for (int i=0;i<4;i++){
        TINYTEST_ASSERT(R.getData()[i] == 4);
    }
    
    DenseMatrix B3{2,2,0};
    B3(0, 0) = 5;
    B3(1, 0) = 6;
    B3(0, 1) = 7;
    B3(1, 1) = 8;
    
    R = B3-move(B);
    for (int i=0;i<4;i++){
        TINYTEST_ASSERT(R.getData()[i] == 4);
    }
    
    DenseMatrix B4{2,2,0};
    B4(0, 0) = 1;
    B4(1, 0) = 2;
    B4(0, 1) = 3;
    B4(1, 1) = 4;
    
    DenseMatrix B5{2,2,0};
    B5(0, 0) = 5;
    B5(1, 0) = 6;
    B5(0, 1) = 7;
    B5(1, 1) = 8;
    
    R = move(B5)-move(B4);
    for (int i=0;i<4;i++){
        TINYTEST_ASSERT(R.getData()[i] == 4);
    }
    
    return 1;
}

int SparseSymmetrize(){
    SparseMatrix A{3,3};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(1, 0) = 1;
    A(0, 2) = 1;
    A(2, 0) = 1;
    A(1, 2) = 0.5;
    A(2, 1) = 0.5;
    A(2, 2) = -0.5;
    A.build();
    
    SparseMatrix B = A.copy();
    B.symmetrize();
    
    TINYTEST_ASSERT(A.toDense() == B.toDense());
    
    return 1;
}

int SparseBeginEnd(){
    SparseMatrix A{3, 3};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(1, 0) = 1;
    A(0, 2) = 1;
    A(2, 0) = 1;
    
    double sum = 0;
    for (auto d : A){
        sum += d;
    }
    TINYTEST_ASSERT(fabs(sum-5)<0.001);
    A.build();
    sum = 0;
    for (auto d : A){
        sum += d;
    }
    TINYTEST_ASSERT(fabs(sum-5)<0.001);
    
    return 1;
}

int DenseBeginEnd(){
    
    DenseMatrix A{3,3,0};
    A(0, 0) = 1;
    A(0, 1) = 1;
    A(1, 0) = 1;
    A(0, 2) = 1;
    A(2, 2) = 1;
    
    double sum = 0;
    for (auto d : A){
        sum += d;
    }
    TINYTEST_ASSERT(fabs(sum-5)<0.001);
    
    return 1;
}

int SparseMultiply(){
    SparseMatrix S{3,3};
    S(0, 0) = 1;
    S(1, 0) = 2;
    S(2, 0) = 3;
    S(0, 1) = 4;
    S(1, 1) = 5;
    S(2, 1) = 6;
    S(0, 2) = 7;
    S(1, 2) = 8;
    S(2, 2) = 9;
    S.build();
    
    DenseMatrix D{3,3,0};
    D(0, 0) = 10;
    D(1, 0) = 11;
    D(2, 0) = 12;
    D(0, 1) = 13;
    D(1, 1) = 14;
    D(2, 1) = 15;
    D(0, 2) = 16;
    D(1, 2) = 17;
    D(2, 2) = 18;

    DenseMatrix X{3,3,0};
    X(0, 0) = 19;
    X(1, 0) = 20;
    X(2, 0) = 21;
    X(0, 1) = 22;
    X(1, 1) = 23;
    X(2, 1) = 24;
    X(0, 2) = 25;
    X(1, 2) = 26;
    X(2, 2) = 27;

    // alpha*(A*X) + beta*Y
    S.multiply(false, 2, 3, D, X);
    
    DenseMatrix exp{3,3,0};
    exp(0, 0) = 333;
    exp(1, 0) = 402;
    exp(2, 0) = 471;
    exp(0, 1) = 414;
    exp(1, 1) = 501;
    exp(2, 1) = 588;
    exp(0, 2) = 495;
    exp(1, 2) = 600;
    exp(2, 2) = 705;

    TINYTEST_ASSERT(exp == X);
    
    X(0, 0) = 19;
    X(1, 0) = 20;
    X(2, 0) = 21;
    X(0, 1) = 22;
    X(1, 1) = 23;
    X(2, 1) = 24;
    X(0, 2) = 25;
    X(1, 2) = 26;
    X(2, 2) = 27;
    
    // Y = alpha*(A'*X) + beta*Y
    S.multiply(true, 2, 3, D, X);
    
    exp(0, 0) = 193;
    exp(1, 0) = 394;
    exp(2, 0) = 595;
    exp(0, 1) = 238;
    exp(1, 1) = 493;
    exp(2, 1) = 748;
    exp(0, 2) = 283;
    exp(1, 2) = 592;
    exp(2, 2) = 901;
    
    TINYTEST_ASSERT(exp == X);
    
    return 1;
}

int SparseMultiplyNoRes(){
    SparseMatrix S{3,3};
    S(0, 0) = 1;
    S(1, 0) = 2;
    S(2, 0) = 3;
    S(0, 1) = 4;
    S(1, 1) = 5;
    S(2, 1) = 6;
    S(0, 2) = 7;
    S(1, 2) = 8;
    S(2, 2) = 9;
    S.build();
    
    DenseMatrix D{3,3,0};
    D(0, 0) = 10;
    D(1, 0) = 11;
    D(2, 0) = 12;
    D(0, 1) = 13;
    D(1, 1) = 14;
    D(2, 1) = 15;
    D(0, 2) = 16;
    D(1, 2) = 17;
    D(2, 2) = 18;
    
    // alpha*(A*X)
    DenseMatrix X = S.multiply(false, 2, D);
    
    DenseMatrix exp{3,3,0};
    exp(0, 0) = 276;
    exp(1, 0) = 342;
    exp(2, 0) = 408;
    exp(0, 1) = 348;
    exp(1, 1) = 432;
    exp(2, 1) = 516;
    exp(0, 2) = 420;
    exp(1, 2) = 522;
    exp(2, 2) = 624;
    
    TINYTEST_ASSERT(exp == X);
    
    X(0, 0) = 19;
    X(1, 0) = 20;
    X(2, 0) = 21;
    X(0, 1) = 22;
    X(1, 1) = 23;
    X(2, 1) = 24;
    X(0, 2) = 25;
    X(1, 2) = 26;
    X(2, 2) = 27;
    
    // Y = alpha*(A'*X) + beta*Y
    X = S.multiply(true, 2, D);
    
    exp(0, 0) = 136;
    exp(1, 0) = 334;
    exp(2, 0) = 532;
    exp(0, 1) = 172;
    exp(1, 1) = 424;
    exp(2, 1) = 676;
    exp(0, 2) = 208;
    exp(1, 2) = 514;
    exp(2, 2) = 820;
    
    TINYTEST_ASSERT(exp == X);
    
    return 1;
}

int SparseMultiplyTranspose(){
    SparseMatrix S{1,3};
    S(0, 0) = 1;
    S(0, 1) = 4;
    S(0, 2) = 7;
    S.build();
    
    DenseMatrix D{3,1,0};
    D(0, 0) = 10;
    D(1, 0) = 11;
    D(2, 0) = 12;
    
    DenseMatrix X{1,1,0};
    
    // alpha*(A*X) + beta*Y
    S.multiply(false, 1, 0, D, X);
    
    DenseMatrix exp{1,1,0};
    exp(0, 0) = 138;
    TINYTEST_ASSERT(exp == X);
    
    X = S.multiply(false, 1, D);
    TINYTEST_ASSERT(exp == X);
    
    D.transpose();
    X = DenseMatrix{3,3,0};
    S.multiply(true, 1, 0, D, X);
    
    exp = DenseMatrix{3,3,0};
    exp(0, 0) = 10;
    exp(1, 0) = 40;
    exp(2, 0) = 70;
    exp(0, 1) = 11;
    exp(1, 1) = 44;
    exp(2, 1) = 77;
    exp(0, 2) = 12;
    exp(1, 2) = 48;
    exp(2, 2) = 84;
    TINYTEST_ASSERT(exp == X);
    
    X = S.multiply(true, 1, D);
    TINYTEST_ASSERT(exp == X);
    return 1;
}

