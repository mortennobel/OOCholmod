#include <sstream>
#include <cmath>
#include "tinytest.h"
#include "oocholmodtest.h"

TINYTEST_START_SUITE(ObjSuite);
TINYTEST_ADD_TEST(BuildSparseTestObj);
TINYTEST_ADD_TEST(EqualSparseTestObj);
TINYTEST_ADD_TEST(TestCaseFunctionOperatorObj);
TINYTEST_ADD_TEST(SolveDenseDenseTestObj);
TINYTEST_ADD_TEST(SolveSparseDenseTestObj);
TINYTEST_ADD_TEST(SolveSparseSparseTestObj);
TINYTEST_ADD_TEST(SolveSparseDenseFactorTestObj);
TINYTEST_ADD_TEST(AddSparseSparseTestObj);
TINYTEST_ADD_TEST(AddDenseDenseTestObj);
TINYTEST_ADD_TEST(AddEqualDenseDenseTestObj);
TINYTEST_ADD_TEST(SubtractSparseSparseTestObj);
TINYTEST_ADD_TEST(TransposeDenseTestObj);
TINYTEST_ADD_TEST(InverseDenseTestObj);
TINYTEST_ADD_TEST(MultiplySparseSparseTestObj);
TINYTEST_ADD_TEST(MultiplyScalarSparseTestObj);
TINYTEST_ADD_TEST(MultiplyScalarDenseTestObj);
TINYTEST_ADD_TEST(MultiplyEqualScalarDenseTestObj);
TINYTEST_ADD_TEST(MultiplyDenseDenseTestObj);
TINYTEST_ADD_TEST(MultiplySparseDenseTestObj);
TINYTEST_ADD_TEST(FillTestObj);
TINYTEST_ADD_TEST(DotTestObj);
TINYTEST_ADD_TEST(LengthTestObj);
TINYTEST_ADD_TEST(ElemDivideTestObj);
TINYTEST_ADD_TEST(ElemMultiplyTestObj);
TINYTEST_ADD_TEST(SingularTestObj);
TINYTEST_ADD_TEST(SwapTest);
TINYTEST_ADD_TEST(SparseStateTest);
TINYTEST_ADD_TEST(DynamicTripletGrow);
TINYTEST_ADD_TEST(IndexTest2);
TINYTEST_ADD_TEST(IndexTest);
TINYTEST_ADD_TEST(LargeSparseMatrix);
TINYTEST_ADD_TEST(LargeMatrixPerformance);
TINYTEST_ADD_TEST(DropSmallEntriesTest);
TINYTEST_ADD_TEST(CopyTest);
TINYTEST_ADD_TEST(NormTest);
TINYTEST_ADD_TEST(DeterminantTest);
TINYTEST_ADD_TEST(AppendTest);
TINYTEST_ADD_TEST(NumberOfElementsTest);
TINYTEST_ADD_TEST(ZeroTest);
TINYTEST_ADD_TEST(DenseSetGetTest);
TINYTEST_ADD_TEST(SparseToDense);
TINYTEST_ADD_TEST(DenseSubstractionTest);
TINYTEST_ADD_TEST(SparseSymmetrize);
TINYTEST_ADD_TEST(SparseBeginEnd);
TINYTEST_ADD_TEST(DenseBeginEnd);

TINYTEST_ADD_TEST(SparseMultiply);
TINYTEST_ADD_TEST(SparseMultiplyNoRes);
TINYTEST_ADD_TEST(SparseMultiplyTranspose);

TINYTEST_ADD_TEST(DenseMultiply);
TINYTEST_ADD_TEST(DenseMultiplyNoRes);
TINYTEST_ADD_TEST(DenseMultiplyTranspose);

TINYTEST_END_SUITE();

TINYTEST_START_MAIN();
TINYTEST_RUN_SUITE(ObjSuite);
TINYTEST_END_MAIN();