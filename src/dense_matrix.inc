
#include <float.h>
#ifdef WIN32
#define isnan _isnan
#endif

#ifdef USE_ACML
#include "acml.h"
typedef int __CLPK_integer;
typedef double __CLPK_doublereal;
#define CblasColMajor
#define CblasTrans 'T'
#define CblasNoTrans 'N'
#define CblasRight 'R'
#define CblasLeft 'L'
#define CblasUpper 'U'
#else
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif
#ifdef _WIN32
#include <clapack.h>
#endif
#endif
#include "dense_matrix.h"
#include "config_singleton.h"
#include "sparse_matrix.h"
#include "ooc_exception.h"



namespace oocholmod {

    using namespace std;

    inline DenseMatrix::DenseMatrix(unsigned int rows, unsigned int cols, double value)
            :nrow(rows), ncol(cols)
    {
        dense = cholmod_allocate_dense(rows, cols, rows /* leading dimension (equal rows) */ , CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        if (!isnan(value)) {
            fill(value);
        } else {
#ifdef DEBUG
            fill(NAN);
#endif
        }
    }

    inline DenseMatrix::DenseMatrix(cholmod_dense *dense_)
            :dense(dense_), nrow(static_cast<unsigned int>(dense_->nrow)), ncol(static_cast<unsigned int>(dense_->ncol))
    {
    }

    inline DenseMatrix::DenseMatrix(DenseMatrix&& move) ndebug_noexcept
            :dense(move.dense), nrow(move.nrow), ncol(move.ncol)
    {
        move.dense = nullptr;
        move.nrow = 0;
        move.ncol = 0;
    }

    inline DenseMatrix& DenseMatrix::operator=(DenseMatrix&& other) ndebug_noexcept
    {
        if (this != &other)
        {
            if (dense){
                cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
            }
            dense = other.dense;
            nrow = other.nrow;
            ncol = other.ncol;

            other.dense = nullptr;
            other.nrow = 0;
            other.ncol = 0;
        }
        return *this;
    }

    inline DenseMatrix::~DenseMatrix()
    {
        if (dense){
            cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
            dense = nullptr;
        }
    }

    inline double *DenseMatrix::begin(){
        return getData();
    }

    inline double *DenseMatrix::end(){
        return getData()+(nrow * ncol);
    }


    inline void DenseMatrix::zero(){
        memset(dense->x, 0, nrow * ncol * sizeof(double));
    }

    inline SparseMatrix DenseMatrix::toSparse() const {
        return SparseMatrix(cholmod_dense_to_sparse(dense, true, ConfigSingleton::getCommonPtr()));
    }

    inline void DenseMatrix::fill(double value)
    {
        double *data = getData();
        for (int i = 0; i < ncol * nrow; i++){
            data[i] = value;
        }
    }

    inline double DenseMatrix::norm(int norm) const{
        return cholmod_norm_dense(dense, norm, ConfigSingleton::getCommonPtr());
    }

    inline double DenseMatrix::determinant() const
    {
#ifdef DEBUG
        assert(ncol == nrow);
#endif
        const double *data = getData();
        if(ncol == 1)
        {
            return data[0];
        }
        if(ncol == 2)
        {
            return data[0]*data[3] - data[1]*data[2];
        }
        if(ncol == 3)
        {
            return data[0]*data[4]*data[8] + data[3]*data[7]*data[2] + data[6]*data[1]*data[5] - data[0]*data[7]*data[5] - data[3]*data[1]*data[8] - data[6]*data[4]*data[2];
        }

        double det = 0.0;
        for (int i = 0; i < ncol; i++)
        {
            double a = 1.0, b = 1.0;
            for (int row = 0; row < ncol; row++)
            {
                a *= data[row + ncol * ((i+row)%ncol)];
                b *= data[row + ncol * ((ncol-1) - (i+row)%ncol)];
            }
            det += a - b;
        }
        return det;
    }

    inline double DenseMatrix::dot(const DenseMatrix& b) const {
#ifdef DEBUG
        assert(ncol == 1 || nrow == 1);
        assert(b.ncol == ncol && b.nrow == nrow);
#endif
#ifdef USE_ACML
        return ddot(nrow*ncol, const_cast<double*>(getData()), 1, const_cast<double*>(b.getData()), 1);
#else
        return cblas_ddot(nrow*ncol, getData(), 1, b.getData(), 1);
#endif
    }

    inline double DenseMatrix::length() const {
#ifdef DEBUG
        assert(ncol == 1 || nrow == 1);
#endif
#ifdef USE_ACML
        return dnrm2(ncol*nrow, const_cast<double*>(getData()), 1);
#else
        return cblas_dnrm2(ncol*nrow, getData(), 1);
#endif
    }

    inline void DenseMatrix::elemDivide(const DenseMatrix& b, DenseMatrix& dest) const {
#ifdef DEBUG
        assert(nrow == b.getRows() && ncol == b.getColumns());
#endif
        const double *thisData = getData();
        double *destData = dest.getData();
        const double *bData = b.getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                destData[c*nrow + r] = thisData[c*nrow + r] / bData[c*nrow + r];
            }
        }
    }

    inline void DenseMatrix::elemDivide(const DenseMatrix& b){
        elemDivide(b, *this);
    }

    inline void DenseMatrix::elemMultiply(const DenseMatrix& b){
        elemMultiply(b, *this);
    }

    inline void DenseMatrix::elemMultiply(const DenseMatrix& b, DenseMatrix& dest) const {
#ifdef DEBUG
        assert(nrow == b.getRows() && ncol == b.getColumns());
#endif
        const double *thisData = getData();
        double *destData = dest.getData();
        const double *bData = b.getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                destData[c*nrow + r] = thisData[c*nrow + r]*bData[c*nrow + r];
            }
        }
    }

    inline DenseMatrix DenseMatrix::copy() const{
        DenseMatrix dest(nrow,ncol);
        const double *srcPtr = getData();
        double *destPtr = dest.getData();
        memcpy(destPtr, srcPtr, nrow*ncol*sizeof(double));
        return dest;
    }

    inline void DenseMatrix::set(float *inData){
#ifdef DEBUG
        assert(dense);
#endif
        double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                data[c*nrow + r] = inData[c*nrow + r];
            }
        }
    }

    inline void DenseMatrix::set(double *data){
#ifdef DEBUG
        assert(dense);
#endif
        memcpy(dense->x, data, nrow*ncol*sizeof(double));
    }

    inline void DenseMatrix::get(double *outData) const {
#ifdef DEBUG
        assert(dense);
#endif
        memcpy(outData, dense->x, nrow*ncol*sizeof(double));
    }

    inline void DenseMatrix::get(float *outData) const {
        const double *data = getData();
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                outData[c*nrow + r] = (float)data[c*nrow + r];
            }
        }
    }

    inline void DenseMatrix::swap(DenseMatrix& other){
        std::swap(dense, other.dense);
        std::swap(nrow, other.nrow);
        std::swap(ncol, other.ncol);
    }

    inline void swap(DenseMatrix& v1, DenseMatrix& v2) {
        v1.swap(v2);
    }

    // Addition
    inline DenseMatrix operator+(const DenseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        DenseMatrix res = RHS.copy();
#ifdef USE_ACML
        daxpy(LHS.nrow*LHS.ncol, 1., const_cast<double*>(LHS.getData()), 1, res.getData(), 1);
#else
        cblas_daxpy(LHS.nrow*LHS.ncol, 1., LHS.getData(), 1, res.getData(), 1);
#endif
        return res;
    }

    inline DenseMatrix&& operator+(DenseMatrix&& LHS, const DenseMatrix& RHS)
    {
        return RHS+move(LHS);
    }

    inline DenseMatrix&& operator+(const DenseMatrix& LHS, DenseMatrix&& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
#ifdef USE_ACML
        daxpy(LHS.nrow*LHS.ncol, 1., const_cast<double*>(LHS.getData()), 1, RHS.getData(), 1);
#else
        cblas_daxpy(LHS.nrow*LHS.ncol, 1., LHS.getData(), 1, RHS.getData(), 1);
#endif
        return move(RHS);
    }

    inline DenseMatrix&& operator+(DenseMatrix&& LHS, DenseMatrix&& RHS)
    {
        return move(LHS)+RHS;
    }

    inline DenseMatrix& DenseMatrix::operator+=(const DenseMatrix& RHS)
    {
        RHS + move(*this);
        return *this;
    }

    // Subtraction
    inline DenseMatrix operator-(const DenseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        return LHS + (-RHS);
    }

    inline DenseMatrix&& operator-(DenseMatrix&& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        return move(LHS) + (-RHS);
    }

    inline DenseMatrix&& operator-(const DenseMatrix& LHS, DenseMatrix&& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
#endif
        return LHS + (-move(RHS));
    }

    inline DenseMatrix&& operator-(DenseMatrix&& LHS, DenseMatrix&& RHS)
    {
        return move(LHS)-RHS;
    }

    inline DenseMatrix& DenseMatrix::operator-=(const DenseMatrix& RHS)
    {
        move(*this) - RHS;
        return *this;
    }

    inline DenseMatrix operator-(const DenseMatrix& M)
    {
        return -move(M.copy());
    }

    inline DenseMatrix&& operator-(DenseMatrix&& M)
    {
        for (int r = 0; r < M.nrow; r++) {
            for (int c = 0; c < M.ncol; c++) {
                M(r,c) = -M(r,c);
            }
        }
        return move(M);
    }

    // Multiplication
    inline DenseMatrix operator*(const DenseMatrix& LHS, double RHS)
    {
#ifdef DEBUG
        assert(LHS.dense);
#endif
        cholmod_dense *dense = cholmod_copy_dense(LHS.dense, ConfigSingleton::getCommonPtr());
        DenseMatrix res(dense);
#ifdef USE_ACML
        dscal (LHS.nrow*LHS.ncol, RHS, res.getData(), 1);
#else
        cblas_dscal (LHS.nrow*LHS.ncol, RHS, res.getData(), 1);
#endif
        return res;
    }

    inline DenseMatrix&& operator*(DenseMatrix&& LHS, double RHS)
    {
#ifdef DEBUG
        assert(LHS.dense);
#endif
#ifdef USE_ACML
        dscal(LHS.nrow*LHS.ncol, RHS, LHS.getData(), 1);
#else
        cblas_dscal (LHS.nrow*LHS.ncol, RHS, LHS.getData(), 1);
#endif
        return move(LHS);
    }

    inline DenseMatrix operator*(double LHS, const DenseMatrix& RHS)
    {
        return RHS*LHS;
    }

    inline DenseMatrix&& operator*(double LHS, DenseMatrix&& RHS)
    {
        return move(RHS)*LHS;
    }

    inline DenseMatrix& DenseMatrix::operator*=(double RHS)
    {
        move(*this) * RHS;
        return *this;
    }

    inline void DenseMatrix::multiply(bool transposeA, bool transposeB, double alpha, double beta, const DenseMatrix& B,
            DenseMatrix& C){
        auto transA = transposeA ? CblasTrans : CblasNoTrans;
        auto transB = transposeB ? CblasTrans : CblasNoTrans;
        int M = transposeA ? ncol : nrow;
        int N = transposeB ? B.nrow : B.ncol;
        int K = transposeA ? nrow : ncol;
        int LDA = nrow; // stride
        int LDB = B.nrow; // stride
        int LDC = C.nrow; // stride
#ifdef DEBUG
        assert(M == C.nrow);
        assert(N == C.ncol);
#endif
#ifdef USE_ACML
        dgemm(transA, transB, M, N, K, alpha, getData(), LDA, const_cast<double*>(B.getData()), LDB, beta, C.getData(), LDC);
#else
        cblas_dgemm(CblasColMajor, transA, transB, M, N, K, alpha, getData(), LDA, B.getData(), LDB, beta, C.getData(), LDC);
#endif
    }

    inline DenseMatrix DenseMatrix::multiply(bool transposeA, bool transposeB, double alpha, const DenseMatrix& B){
        unsigned int M = transposeA ? ncol : nrow;
        unsigned int N = transposeB ? B.nrow : B.ncol;
        DenseMatrix C{M, N, 0};
        multiply(transposeA, transposeB, alpha, 0, B, C);
        return C;
    }

    inline DenseMatrix operator*(const DenseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.ncol == RHS.nrow);
#endif
        DenseMatrix res(LHS.nrow, RHS.ncol);

        if(RHS.ncol == 1)
        {
#ifdef USE_ACML
            dgemv(CblasNoTrans, LHS.nrow, LHS.ncol, 1., const_cast<double*>(LHS.getData()), LHS.nrow, const_cast<double*>(RHS.getData()), 1, 0., res.getData(), 1);
#else
            cblas_dgemv(CblasColMajor, CblasNoTrans, LHS.nrow, LHS.ncol, 1., LHS.getData(), LHS.nrow, RHS.getData(), 1, 0., res.getData(), 1);
#endif
        }
        else {
#ifdef USE_ACML
            dgemm(CblasNoTrans, CblasNoTrans, LHS.nrow, RHS.ncol, LHS.ncol, 1., const_cast<double*>(LHS.getData()), LHS.nrow, const_cast<double*>(RHS.getData()), RHS.nrow, 0., res.getData(), res.nrow);
#else
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, LHS.nrow, RHS.ncol, LHS.ncol, 1., LHS.getData(), LHS.nrow, RHS.getData(), RHS.nrow, 0., res.getData(), res.nrow);
#endif
        }
        return res;
    }

    inline void DenseMatrix::multiplyAddTo(double alpha, DenseMatrix& C){
#ifdef DEBUG
        assert(getRows() == C.getRows());
        assert(getColumns() == C.getColumns());
#endif
        int length = nrow * ncol;
        double *data = getData();
        double *cData = C.getData();
#ifdef USE_ACML
        daxpy(length, alpha, data, 1, cData, 1);
#else
        cblas_daxpy(length, alpha, data, 1, cData, 1);
#endif
    }


    inline bool DenseMatrix::operator==(const DenseMatrix& RHS) const {
        if (nrow != RHS.nrow || ncol != RHS.ncol){
            return false;
        }
        const double *data = getData();
        const double *rhsData = RHS.getData();
        for (int i=0;i<nrow*ncol;i++){
            if (data[i] != rhsData[i]){
                return false;
            }
        }
        return true;
    }

    inline bool DenseMatrix::operator!=(const DenseMatrix& RHS) const{
        return !(*this == RHS);
    }

    // Transpose
    inline void DenseMatrix::transpose()
    {
        double *data = getData();
        cholmod_dense *d = cholmod_allocate_dense(ncol, nrow, ncol, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        double *outData = (double*)d->x;
        for (int c = 0; c < ncol; c++){
            for (int r = 0; r < nrow; r++){
                outData[r*ncol + c] = data[c*nrow + r];
            }
        }
        cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
        dense = d;

        std::swap(nrow, ncol);
    }

    inline DenseMatrix transposed(const DenseMatrix& M)
    {
        DenseMatrix res(M.ncol, M.nrow);
        const double *data = M.getData();
        double *outData = res.getData();
        for (int c = 0; c < M.ncol; c++){
            for (int r = 0; r < M.nrow; r++){
                outData[r*M.ncol + c] = data[c*M.nrow + r];
            }
        }
        return res;
    }

    inline DenseMatrix&& transposed(DenseMatrix&& M)
    {
        M.transpose();
        return move(M);
    }

    inline void DenseMatrix::inverse()
    {
#ifdef DEBUG
        assert(dense);
        assert(nrow == ncol);
#endif
        __CLPK_integer N = nrow;
        __CLPK_integer lwork = N*N;
        __CLPK_integer *ipiv = new __CLPK_integer[N+1];
        __CLPK_doublereal *work = new __CLPK_doublereal[lwork];
        __CLPK_integer info;
#ifdef USE_ACML
        dgetrf(N, N, (double*)dense->x, N, ipiv, &info);
#else
        dgetrf_(&N, &N, (double*)dense->x, &N, ipiv, &info);
#endif
        if (info != 0){
            if (info < 0){
                OOCException::createOOCException(string{"The "}+to_string(-info)+" -th argument had an illegal value");
            }
            if (info > 0){
                OOCException::createOOCException(string{"U("}+to_string(-info)+","+to_string(-info)+") is exactly zero. The factor U is exactly singular.");
            }
        }
#ifdef USE_ACML
        dgetri(N, (double*)dense->x, N, ipiv, &info);
#else
        dgetri_(&N, (double*)dense->x, &N, ipiv, work, &lwork, &info);
#endif
        delete[] ipiv;
        delete[] work;
    }

    inline DenseMatrix inversed(const DenseMatrix& M)
    {
#ifdef DEBUG
        assert(M.dense);
        assert(M.nrow == M.ncol);
#endif


        __CLPK_integer N = M.nrow;
        __CLPK_integer lwork = N*N;
        __CLPK_integer *ipiv = new __CLPK_integer[N+1];
        __CLPK_doublereal *work = new __CLPK_doublereal[lwork];
        __CLPK_integer info;

        cholmod_dense *res = cholmod_copy_dense(M.dense, ConfigSingleton::getCommonPtr());
#ifdef USE_ACML
        dgetrf(N, N, (double*)res->x, N, ipiv, &info);
#else
        dgetrf_(&N, &N, (double*)res->x, &N, ipiv, &info);
#endif
        if (info != 0){
            if (info < 0){
                OOCException::createOOCException(string{"The "}+to_string(-info)+" -th argument had an illegal value");
            }
            if (info > 0){
                OOCException::createOOCException(string{"U("}+to_string(-info)+","+to_string(-info)+") is exactly zero. The factor U is exactly singular.");
            }
        }
#ifdef USE_ACML
        dgetri(N, (double*)res->x, N, ipiv, &info);
#else
        dgetri_(&N, (double*)res->x, &N, ipiv, work, &lwork, &info);
#endif
        delete[] ipiv;
        delete[] work;
        return DenseMatrix(res);
    }

    inline DenseMatrix&& inversed(DenseMatrix&& M)
    {
        M.inverse();
        return move(M);
    }

    inline DenseMatrix solve(const DenseMatrix& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer* ipiv = new __CLPK_integer[N];
        __CLPK_integer info;

        cholmod_dense *a = cholmod_copy_dense(A.dense, ConfigSingleton::getCommonPtr());
        cholmod_dense *res = cholmod_copy_dense(b.dense, ConfigSingleton::getCommonPtr());
#ifdef USE_ACML
        dgesv(N, nrhs, (double*)a->x, lda, ipiv, (double*)res->x, ldb, &info);
#else
        dgesv_(&N, &nrhs, (double*)a->x, &lda, ipiv, (double*)res->x, &ldb, &info);
#endif
        if (info != 0){
            if (info < 0){
                OOCException::createOOCException(string{"The "}+to_string(-info)+" -th argument had an illegal value");
            }
            if (info > 0){
                OOCException::createOOCException(string{"U("}+to_string(-info)+","+to_string(-info)+") is exactly zero. The factor U is exactly singular.");
            }
        }
        delete [] ipiv;
#ifdef DEBUG
        assert(info == 0);
#endif
        return DenseMatrix(res);
    }

    inline DenseMatrix solve(DenseMatrix&& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer* ipiv = new __CLPK_integer[N];
        __CLPK_integer info;

        cholmod_dense *res = cholmod_copy_dense(b.dense, ConfigSingleton::getCommonPtr());
#ifdef USE_ACML
        dgesv(N, nrhs, A.getData(), lda, ipiv, (double*)res->x, ldb, &info);
#else
        dgesv_(&N, &nrhs, A.getData(), &lda, ipiv, (double*)res->x, &ldb, &info);
#endif
        if (info != 0){
            if (info < 0){
                OOCException::createOOCException(string{"The "}+to_string(-info)+" -th argument had an illegal value");
            }
            if (info > 0){
                OOCException::createOOCException(string{"U("}+to_string(-info)+","+to_string(-info)+") is exactly zero. The factor U is exactly singular.");
            }
        }
        delete [] ipiv;
#ifdef DEBUG
        assert(info == 0);
#endif
        return DenseMatrix(res);
    }

    inline DenseMatrix&& solve(const DenseMatrix& A, DenseMatrix&& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer* ipiv = new __CLPK_integer[N];
        __CLPK_integer info;

        cholmod_dense *a = cholmod_copy_dense(A.dense, ConfigSingleton::getCommonPtr());
#ifdef USE_ACML
        dgesv(N, nrhs, (double*)a->x, lda, ipiv, b.getData(), ldb, &info);
#else
        dgesv_(&N, &nrhs, (double*)a->x, &lda, ipiv, b.getData(), &ldb, &info);
#endif
        if (info != 0){
            if (info < 0){
                OOCException::createOOCException(string{"The "}+to_string(-info)+" -th argument had an illegal value");
            }
            if (info > 0){
                OOCException::createOOCException(string{"U("}+to_string(-info)+","+to_string(-info)+") is exactly zero. The factor U is exactly singular.");
            }
        }
        delete [] ipiv;
#ifdef DEBUG
        assert(info == 0);
#endif
        return move(b);
    }

    inline DenseMatrix&& solve(DenseMatrix&& A, DenseMatrix&& b)
    {
#ifdef DEBUG
        assert(A.dense && b.dense);
        assert(A.nrow == A.ncol);
        assert(A.nrow == b.nrow);
#endif
        __CLPK_integer N = b.nrow;
        __CLPK_integer nrhs = b.ncol;
        __CLPK_integer lda = A.nrow;
        __CLPK_integer ldb = b.nrow;
        __CLPK_integer* ipiv = new __CLPK_integer[N];
        __CLPK_integer info;
#ifdef USE_ACML
        dgesv(N, nrhs, A.getData(), lda, ipiv, b.getData(), ldb, &info);
#else
        dgesv_(&N, &nrhs, A.getData(), &lda, ipiv, b.getData(), &ldb, &info);
#endif
        if (info != 0){
            if (info < 0){
                OOCException::createOOCException(string{"The "}+to_string(-info)+" -th argument had an illegal value");
            }
            if (info > 0){
                OOCException::createOOCException(string{"U("}+to_string(-info)+","+to_string(-info)+") is exactly zero. The factor U is exactly singular.");
            }
        }
        delete [] ipiv;
#ifdef DEBUG
        assert(info == 0);
#endif
        return move(b);
    }

    inline ostream& operator<<(ostream& os, const DenseMatrix& A)
    {
        os << endl;
        if (A.dense)
        {
            cholmod_print_dense(A.dense, "", ConfigSingleton::getCommonPtr());
            os << "[";
            for (int r = 0; r < A.getRows(); r++)
            {
                for (int c = 0; c < A.getColumns(); c++)
                {
                    os << A(r,c);
                    if(c < A.getColumns()-1) {
                        os << ", ";
                    }
                }
                if(r < A.getRows()-1)
                {
                    os << ";" << endl << " ";
                }
            }
            os << "];" << endl;
        }
        else {
            os << "[Empty dense matrix]" << endl;
        }
        return os;
    }
}
