#include <sstream>
#include <cmath>
#include "tinytest.h"

using namespace std;

#define assertEqual(E, A, L) { \
    double *expected_ = (E); \
    double *actual_ = (A); \
    int length_ = (L); \
    for (int i=0;i<length_;i++){ \
        if (fabs(expected_[i] - actual_[i]) > 0.01){ \
            stringstream str_;\
            str_ << "Vectors are not equal. Element "; \
            str_ << i; \
            str_ << " was "; \
            str_ << actual_[i]; \
            str_ << " but expected "; \
            str_ << expected_[i]; \
            TINYTEST_ASSERT_MSG(false, str_.str().c_str()); \
        } \
    } \
}

#include "objects.h"

TINYTEST_START_SUITE(ObjSuite);
TINYTEST_ADD_TEST(TestCaseFunctionOperatorObj);
TINYTEST_ADD_TEST(SolveDenseDenseTestObj);
TINYTEST_ADD_TEST(SolveSparseDenseTestObj);
TINYTEST_ADD_TEST(SolveSparseSparseTestObj);
TINYTEST_ADD_TEST(SolveSparseDenseFactorTestObj);
TINYTEST_ADD_TEST(AddSparseSparseTestObj);
TINYTEST_ADD_TEST(AddDenseDenseTestObj);
TINYTEST_ADD_TEST(AddEqualDenseDenseTestObj);
TINYTEST_ADD_TEST(TransposeDenseTestObj);
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
TINYTEST_END_SUITE();

TINYTEST_START_MAIN();
TINYTEST_RUN_SUITE(ObjSuite);
TINYTEST_END_MAIN();