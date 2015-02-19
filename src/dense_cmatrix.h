//
// Created by Morten Nobel-Jørgensen on 16/02/15.
//


#pragma once

#include <complex>
#include <cassert>
#include <cholmod.h>

#ifdef WIN32
#ifndef NAN
	static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
    #define NAN (*(const float *) __nan)
#endif
#endif

#ifndef ndebug_noexcept
#   ifdef DEBUG
#       define ndebug_noexcept
#   else
#       define ndebug_noexcept noexcept
#   endif
#endif

namespace oocholmod {
    class DenseCMatrix {
    public:
        DenseCMatrix(unsigned int rows = 0, unsigned int cols = 1);

        DenseCMatrix(DenseCMatrix&& move) ndebug_noexcept;
        DenseCMatrix& operator=(DenseCMatrix&& other) ndebug_noexcept;
        
        std::complex<double>& operator()(unsigned int row, unsigned int col = 0) ndebug_noexcept;

        std::complex<double> operator()(unsigned int row, unsigned int col = 0) const ndebug_noexcept;

        void fill(std::complex<double> value);

        std::complex<double> *getData();

        std::complex<double> const *getData() const;

        friend DenseCMatrix operator*(const DenseCMatrix& LHS, const DenseCMatrix& RHS);
    private:
        cholmod_dense *dense;
        unsigned int nrow;
        unsigned int ncol;


    };

    inline std::complex<double> &DenseCMatrix::operator()(unsigned int row, unsigned int col) {
#ifdef DEBUG
        assert(dense);
        assert(row < nrow && col < ncol);
#endif
        return ((std::complex<double>*)dense->x)[col*nrow + row];
    }

    inline std::complex<double> DenseCMatrix::operator()(unsigned int row, unsigned int col) const {
#ifdef DEBUG
        assert(dense);
        assert(row < nrow && col < ncol);
#endif
        return ((std::complex<double>*)dense->x)[col*nrow + row];
    }

    inline std::complex<double> *DenseCMatrix::getData(){ return (std::complex<double> *)(dense->x); };
    inline const std::complex<double> *DenseCMatrix::getData() const { return (std::complex<double> *)(dense->x); };

    DenseCMatrix operator*(const DenseCMatrix &LHS, const DenseCMatrix &RHS);
}


