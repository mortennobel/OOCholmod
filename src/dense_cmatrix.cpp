//
// Created by Morten Nobel-Jørgensen on 16/02/15.
//

#include "dense_cmatrix.h"
#include "config_singleton.h"
#include <float.h>
#include <cmath>
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



namespace oocholmod {

    DenseCMatrix::DenseCMatrix(unsigned int rows, unsigned int cols)
    : nrow(rows), ncol(cols) {
        dense = cholmod_allocate_dense(rows, cols, rows /* leading dimension (equal rows) */ , CHOLMOD_COMPLEX, ConfigSingleton::getCommonPtr());
        std::complex<double> value(0,0);
        fill(value);
    }
    
    DenseCMatrix::DenseCMatrix(DenseCMatrix&& other) ndebug_noexcept {
        dense = other.dense;
        other.dense = nullptr;
        nrow = other.nrow;
        ncol = other.ncol;
    }
    
    DenseCMatrix& DenseCMatrix::operator=(DenseCMatrix&& other) ndebug_noexcept {
        if (&other != this){
            dense = other.dense;
            other.dense = nullptr;
            nrow = other.nrow;
            ncol = other.ncol;
        }
        return *this;
    }


    void DenseCMatrix::fill(std::complex<double> value){
        auto data = getData();
        for (int i = 0; i < ncol * nrow; i++){
            data[i] = value;
        }
    }

    DenseCMatrix operator*(const DenseCMatrix &LHS, const DenseCMatrix &RHS) {
#ifdef DEBUG
        assert(LHS.dense && RHS.dense);
        assert(LHS.ncol == RHS.nrow);
#endif
        DenseCMatrix res(LHS.nrow, RHS.ncol);

        double scaling[2] ={1,0};
        double scalingC[2] ={0,0};
#ifdef USE_ACML
        zgemm(CblasNoTrans, CblasNoTrans, LHS.nrow, RHS.ncol, LHS.ncol, scaling, const_cast<double*>(LHS.dense->x), LHS.nrow, const_cast<double*>(RHS.dense->x), RHS.nrow, scalingC, const_cast<double*>(res.dense->x), res.nrow);
#else
        cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, LHS.nrow, RHS.ncol, LHS.ncol, scaling, LHS.dense->x, LHS.nrow, RHS.dense->x, RHS.nrow, scalingC, res.dense->x, res.nrow);
#endif

        return res;
    }
}
