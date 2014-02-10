//
//  oocholmodtest.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 28/10/13.
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#ifndef __OOCholmod__oocholmodtest__
#define __OOCholmod__oocholmodtest__


int BuildSparseTestObj();

int EqualSparseTestObj();

int SolveDenseDenseTestObj();

int SolveSparseDenseTestObj();

int SolveSparseSparseTestObj();

int SolveSparseDenseFactorTestObj();

int AddSparseSparseTestObj();

int SubtractSparseSparseTestObj();

int AddDenseDenseTestObj();

int AddEqualDenseDenseTestObj();

int TransposeDenseTestObj();

int InverseDenseTestObj();

int TestCaseFunctionOperatorObj();

int MultiplySparseSparseTestObj();

int MultiplyScalarSparseTestObj();

int MultiplySparseDenseTestObj();

int MultiplyScalarDenseTestObj();

int MultiplyEqualScalarDenseTestObj();

int MultiplyDenseDenseTestObj();

int FillTestObj();

int DotTestObj();

int LengthTestObj();

int ElemDivideTestObj();

int ElemMultiplyTestObj();

int SingularTestObj();

int SwapTest();
int SparseStateTest();

int DynamicTripletGrow();

int LargeSparseMatrix();

int IndexTest2();

int IndexTest();

int LargeMatrixPerformance();

int DropSmallEntriesTest();

int CopyTest();

int NormTest();

int DeterminantTest();

int AppendTest();

int NumberOfElementsTest();

int ZeroTest();

int DenseSetGetTest();

int SparseToDense();

int DenseSubstractionTest();

int SparseSymmetrize();

int SparseBeginEnd();

int DenseBeginEnd();

#endif /* defined(__OOCholmod__oocholmodtest__) */
