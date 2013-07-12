//
//  main.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0 

#include <iostream>
#include <vector>
#include <cmath>

#include "CholmodSparse.h"
#include "CholmodFactor.h"
#include "CholmodDenseVector.h"
#include "cpp14.h"

using namespace std;

void assertEqual(double *expected, double *actual, int length){
    for (int i=0;i<length;i++){
        if (fabs(expected[i] - actual[i]) > 0.01){
            cout << "Vectors are not equal. Element "<<i<< " was " << actual[i]<< " but expected "<< expected[i] << endl;
        }
    }
    cout << "Equal vectors"<< endl;
}

void TestCase(){
    cholmod_common com;
    cholmod_start(&com);
    //auto A = make_unique<CholmodSparse>(3,3,&com);
    auto A = new CholmodSparse(3,3,&com);
    
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
    
    CholmodDenseVector * b = new CholmodDenseVector(3, &com);
    (*b)[0] = 6;
    (*b)[1] = -4;
    (*b)[2] = 27;
    b->print("b");
    
    A->print("A");
    CholmodFactor *factor = A->analyzePtr();
    bool res = factor->factorize(A);
    cout << "factor->factorize(A) "<<res<<endl;
    CholmodDenseVector * x = NULL;
    factor->solve(b, &x);
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
    
    factor->factorize(A);
    A->print("A");

    factor->solve(b,&x);
    x->print("x");
    
    double expected2[] = {1.0935,1.76937,-1.73019};
    assertEqual(expected2, &((*x)[0]), 3);
    // note leaking memory here
}

void MultiplyTest(){
    cout << "Multiply test"<<endl;
    cholmod_common com;
    cholmod_start(&com);
    CholmodSparse *A = new CholmodSparse(3,3,&com);
    A->initAddValue(0, 0, 1);
    A->initAddValue(0, 1, 1);
    A->initAddValue(0, 2, 1);
    A->initAddValue(1, 2, 5);
    A->initAddValue(2, 2, -1);
    A->build();
    A->print("A");
    
    CholmodDenseVector *x = new CholmodDenseVector(3, &com);
    (*x)[0] = 3;
    (*x)[1] = 7;
    (*x)[2] = 9;

    CholmodDenseVector *res = new CholmodDenseVector(3, &com);

    A->multiply(x, res);
    res->print("b");
    double expected[3] = {19, 48, 29};
    assertEqual(expected, &((*res)[0]), 3);
}

void FillTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodDenseVector *res = new CholmodDenseVector(3, &com);
    res->fill(123);
    res->print("Fill 123 test");
}

void DotTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodDenseVector *a = new CholmodDenseVector(3, &com);
    CholmodDenseVector *b = new CholmodDenseVector(3, &com);
    (*a)[0] = 1;
    (*a)[1] = 2;
    (*a)[2] = 3;
    
    (*b)[0] = 4;
    (*b)[1] = 5;
    (*b)[2] = 6;
    
    double res = a->dot(b);
    double expected = 32;
    cout << "Dot test" << endl;
    assertEqual(&expected, &res, 1);
}

void LengthTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodDenseVector *a = new CholmodDenseVector(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    double res = a->length();
    double expected = sqrt(4*4+5*5+6*6);
    cout << "Length test" << endl;
    assertEqual(&expected, &res, 1);
}

void ScaleTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodDenseVector *a = new CholmodDenseVector(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;

    CholmodDenseVector *b = new CholmodDenseVector(3, &com);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->scale(-2);
    cout << "scale test"<<endl;
    assertEqual(&((*a)[0]), &((*b)[0]), 3);
}


void DivideTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodDenseVector *a = new CholmodDenseVector(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    CholmodDenseVector *b = new CholmodDenseVector(3, &com);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->divideBy(b);
    cout << "divide test"<<endl;
    
    double expected[] = {-0.5,-0.5,-0.5};
    
    assertEqual(expected, &((*a)[0]), 3);
}

void MultiplyVectorTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodDenseVector *a = new CholmodDenseVector(3, &com);
    (*a)[0] = 4;
    (*a)[1] = 5;
    (*a)[2] = 6;
    
    CholmodDenseVector *b = new CholmodDenseVector(3, &com);
    (*b)[0] = -8;
    (*b)[1] = -10;
    (*b)[2] = -12;
    a->multiplyWith(b);
    cout << "multiply test"<<endl;
    
    double expected[] = {4 * -8,5 * -10, 6 * -12};
    
    assertEqual(expected, &((*a)[0]), 3);
}

void SingularTest(){
    cholmod_common com;
    cholmod_start(&com);
    CholmodSparse *A = new CholmodSparse(3,3,&com);

    A->initAddValue(2, 2, 1);
    
    A->build();
    
    CholmodDenseVector * b = new CholmodDenseVector(3, &com);
    (*b)[0] = 0;
    (*b)[1] = 1;
    (*b)[2] = 0;
    b->print("b");
    
    A->print("A");
    CholmodFactor *factor = A->analyzePtr();
    bool res = factor->factorize(A);
    cout << "Factorize ok "<< res << endl;
    CholmodDenseVector * x = NULL;
    factor->solve(b, &x);

}

void TestCaseUniquePtr(){
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
}

void MultiplyTestUniquePtr(){
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
}

void FillTestUniquePtr(){
    cholmod_common com;
    cholmod_start(&com);
    auto res = make_unique<CholmodDenseVector>(3, &com);
    res->fill(123);
    res->print("Fill 123 test");
}

void DotTestUniquePtr(){
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
}

void LengthTestUniquePtr(){
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
}

void ScaleTestUniquePtr(){
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
}


void DivideTestUniquePtr(){
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
}

void MultiplyVectorTestUniquePtr(){
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
}

void SingularTestUniquePtr(){
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
}

int main(int argc, const char * argv[])
{
    TestCase();
    MultiplyTest();
    FillTest();
    DotTest();
    LengthTest();
    ScaleTest();
    DivideTest();
    MultiplyVectorTest();
    SingularTest();
    
    TestCaseUniquePtr();
    MultiplyTestUniquePtr();
    FillTestUniquePtr();
    DotTestUniquePtr();
    LengthTestUniquePtr();
    ScaleTestUniquePtr();
    DivideTestUniquePtr();
    MultiplyVectorTestUniquePtr();
    SingularTestUniquePtr();
    
    
    
    cout << flush;
    return 0;
}
