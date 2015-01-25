//
//  SparseMatrix.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0

#include <cassert>
#include <algorithm>
#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "factor.h"

using namespace std;

namespace oocholmod {
    
    SparseMatrixIter::SparseMatrixIter(const SparseMatrix* sparseMatrix, int pos)
    :sparseMatrix(sparseMatrix), pos(pos)
    {}
    
    SparseMatrix::SparseMatrix(unsigned int nrow, unsigned int ncol, bool symmetric, int maxSize)
    :sparse(nullptr), triplet(nullptr), nrow(nrow), ncol(ncol)
    {
        if (symmetric && nrow == ncol) {
            symmetry = SYMMETRIC_UPPER;
        }
        else {
            symmetry = ASYMMETRIC;
        }
        
        if (maxSize == 0 && symmetric) {
            maxTripletElements = min<int>(100000,(nrow*(ncol+1))/2); // triangular number
        }
        else if (maxSize == 0 && !symmetric) {
            maxTripletElements = min<int>(100000,nrow*ncol); // triangular number
        }
        else {
            maxTripletElements = maxSize;
        }
    }
    
    SparseMatrix::SparseMatrix(cholmod_sparse *sparse)
    :sparse(sparse), triplet(nullptr), nrow(static_cast<unsigned int>(sparse->nrow)),
    ncol(static_cast<unsigned int>(sparse->ncol)),
    values((double*)sparse->x), iRow((int*)sparse->i), jColumn((int*)sparse->p), symmetry(static_cast<Symmetry>(sparse->stype)), maxTripletElements(0)
    {
#ifdef DEBUG
        assert(sparse->itype == CHOLMOD_INT);
#endif
    }
    
    SparseMatrix::SparseMatrix(SparseMatrix&& other) ndebug_noexcept
    :sparse(other.sparse), triplet(other.triplet), nrow(other.nrow), ncol(other.ncol), values(other.values), iRow(other.iRow), jColumn(other.jColumn), symmetry(other.symmetry), maxTripletElements(other.maxTripletElements)
    {
        other.sparse = nullptr;
        other.triplet = nullptr;
        other.values = nullptr;
        other.iRow = nullptr;
        other.jColumn = nullptr;
        other.nrow = 0;
        other.ncol = 0;
    }
    
    SparseMatrix& SparseMatrix::operator=(SparseMatrix&& other) ndebug_noexcept{
        if (this != &other){
            if (sparse != nullptr){
                cholmod_free_sparse(&sparse, ConfigSingleton::getCommonPtr());
            }
            if (triplet != nullptr){
                cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
            }
            
            sparse = other.sparse;
            triplet = other.triplet;
            nrow = other.nrow;
            ncol = other.ncol;
            values = other.values;
            iRow = other.iRow;
            jColumn = other.jColumn;
            symmetry = other.symmetry;
            maxTripletElements = other.maxTripletElements;

            other.sparse = nullptr;
            other.triplet = nullptr;
            other.values = nullptr;
            other.iRow = nullptr;
            other.jColumn = nullptr;
            other.nrow = 0;
            other.ncol = 0;
        }
        return *this;
    }
    
    
    SparseMatrix::~SparseMatrix(){
        if (sparse != nullptr || triplet != nullptr){
            if (sparse != nullptr){
                cholmod_free_sparse(&sparse, ConfigSingleton::getCommonPtr());
                sparse = nullptr;
            }
            if (triplet != nullptr){
                cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
                triplet = nullptr;
            }
        }
    }
    
    SparseMatrixIter SparseMatrix::begin(){
        return SparseMatrixIter(this, 0);
    }
    
    SparseMatrixIter SparseMatrix::end(){
        return SparseMatrixIter(this, getNumberOfElements());
    }
    
    size_t SparseMatrix::getNumberOfElements(){
        switch (getMatrixState()) {
            case INIT:
                return triplet->nnz;
            case BUILT:
                return jColumn[getColumns()];
            default:
                return 0;
        }
    }
    
    DenseMatrix SparseMatrix::toDense() const
    {
#ifdef DEBUG
        assert(sparse);
#endif
        cholmod_dense *dense = cholmod_sparse_to_dense(sparse, ConfigSingleton::getCommonPtr());
        return DenseMatrix(dense);
    }

    void SparseMatrix::append(const SparseMatrix& m) {
#ifdef DEBUG
        assert(sparse == nullptr);
        assert(m.triplet);
#endif
        if (!triplet){
            createTriplet();
        }
        size_t newSize = triplet->nnz + m.triplet->nnz;
        while (triplet->nzmax < newSize){
            increaseTripletCapacity();
        }
        memcpy(values+triplet->nnz, m.values, sizeof(double)*m.triplet->nnz);
        memcpy(jColumn+triplet->nnz, m.jColumn, sizeof(int)*m.triplet->nnz);
        memcpy(iRow+triplet->nnz, m.iRow, sizeof(int)*m.triplet->nnz);
        triplet->nnz += m.triplet->nnz;
    }
    
    bool SparseMatrix::hasElement(unsigned int row, unsigned int column) const {
#ifdef DEBUG
        assert(sparse);
#endif
        return getIndex(row, column) != -1;
    }
    
    void SparseMatrix::dropSmallEntries(double tol){
#ifdef DEBUG
        assert(getMatrixState() == BUILT);
#endif
        cholmod_drop(tol, sparse, ConfigSingleton::getCommonPtr());
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);
    }
    
    MatrixState SparseMatrix::getMatrixState() const {
        if (nrow == 0 && ncol == 0){
            return DESTROYED;
        }
        if (sparse == nullptr && triplet == nullptr){
            return UNINITIALIZED;
        }
        if (triplet != nullptr){
            return INIT;
        }
        return BUILT;
    }
    
    void SparseMatrix::setSymmetry(Symmetry symmetry){
#ifdef DEBUG
        assert(sparse == nullptr);
        assert(triplet == nullptr);
#endif
        this->symmetry = symmetry;
    }
    
    void SparseMatrix::build(){
#ifdef DEBUG
        assert(sparse == nullptr);
#endif
        if (triplet){
            sparse = cholmod_triplet_to_sparse(triplet, triplet->nnz, ConfigSingleton::getCommonPtr());
            cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
            triplet = nullptr;
        }
        else {
            sparse = cholmod_allocate_sparse(nrow, ncol, 0, false, true, symmetry, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        }
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);

        
#ifdef DEBUG
        assert(sparse->itype == CHOLMOD_INT);
        assert(sparse->stype == symmetry);
        assert(sparse->packed);
#endif
    }
   
    void SparseMatrix::sumRows(DenseMatrix& x){
        int idx = 0;
        for (int j=0;j<ncol;j++){
                int iFrom = ((int*)sparse->p)[j];
                int iTo = ((int*)sparse->p)[j+1]-1;
                for (int i=iFrom;i<=iTo;i++){
                        int row = ((int*)sparse->i)[i];
                        x(row)+=((double*)sparse->x)[idx];
                        if (symmetry == SYMMETRIC_UPPER){
                                if (row!=j){
                                        x(j)+=((double*)sparse->x)[idx];
                                }
                        }
                        idx++;
                }
        }

    }


    void SparseMatrix::setNullSpace(DenseMatrix& v){
        int idx = 0;
        for (int j=0;j<ncol;j++){
                int iFrom = ((int*)sparse->p)[j];
                int iTo = ((int*)sparse->p)[j+1]-1;
                for (int i=iFrom;i<=iTo;i++){
                        int row = ((int*)sparse->i)[i];
                        ((double*)sparse->x)[idx] *= v(row)*v(j);
                        idx++;
                }
        }
        for (int i=0;i<ncol;i++){
                if (v(i) == 0){
                        (*this)(i,i)=1.0;
                }
    }
    }
 
    bool SparseMatrix::operator==(const SparseMatrix& RHS) const
    {
        for(int r = 0; r < nrow; r++)
        {
            for(int c = 0; c < ncol; c++)
            {
                if((*this)(r,c) != RHS(r,c))
                {
                    return false;
                }
            }
        }
        return true;
    }
    
    bool SparseMatrix::operator!=(const SparseMatrix& RHS) const {
        return !(*this == RHS);
    }
    
    void SparseMatrix::swap(SparseMatrix& other){
        std::swap(sparse, other.sparse);
        std::swap(triplet, other.triplet);
        std::swap(nrow, other.nrow);
        std::swap(ncol, other.ncol);
        std::swap(values, other.values);
        std::swap(iRow, other.iRow);
        std::swap(jColumn, other.jColumn);
        std::swap(symmetry, other.symmetry);
        std::swap(maxTripletElements, other.maxTripletElements);
    }
    
    void swap(SparseMatrix& v1, SparseMatrix& v2) {
        v1.swap(v2);
    }
    
    void SparseMatrix::symmetrize(){        
#ifdef DEBUG
        // The original must be asymmetric - otherwise return unchanged
        assert(symmetry == ASYMMETRIC);
        assertHasSparse();
#endif
        
        // allocate a new SYMMETRIX UPPER triplet matrix !
        int numberOfNonZeros = (int)cholmod_nnz(sparse,ConfigSingleton::getCommonPtr());
        cholmod_triplet *triplet_symm = cholmod_allocate_triplet(nrow,ncol,
                                                                 numberOfNonZeros,SYMMETRIC_UPPER,CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        
        // Insert the upper part of sparse into the new triplet
        int idx=0;
        for (int col=0;col<ncol;col++){
            int iFrom = ((int*)sparse->p)[col];
            int iTo = ((int*)sparse->p)[col+1]-1;
            for (int i=iFrom;i<=iTo;i++){
                int row = ((int*)sparse->i)[i];
                if (col >= row){ // Upper half
                    ((int*)triplet_symm->i)[triplet_symm->nnz] = row;
                    ((int*)triplet_symm->j)[triplet_symm->nnz] = col;
                    ((double*)triplet_symm->x)[triplet_symm->nnz] = ((double*)sparse->x)[idx];
                    triplet_symm->nnz++;
                }
                idx++;
            }
        }
        
        // Make sure not to leak
        cholmod_free_sparse(&sparse, ConfigSingleton::getCommonPtr());
        
        // build the new sparse matrix
        sparse = cholmod_triplet_to_sparse(triplet_symm, triplet_symm->nnz, ConfigSingleton::getCommonPtr());
        // deallocate the triplet
        cholmod_free_triplet(&triplet_symm, ConfigSingleton::getCommonPtr());
        triplet_symm = nullptr;
        // Set pointers to the sparse data
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);
        
        symmetry = SYMMETRIC_UPPER;
    }
    
    Factor SparseMatrix::analyze() const
    {
#ifdef DEBUG
        assertHasSparse();
#endif
        cholmod_factor *L = cholmod_analyze(sparse, ConfigSingleton::getCommonPtr());
        return Factor(L);
    }
    
    void SparseMatrix::zero(){
#ifdef DEBUG
        assertHasSparse();
#endif
        memset(values, 0, sparse->nzmax * sizeof(double));
    }
    
    void SparseMatrix::assertValidIndex(unsigned int row, unsigned int column) const
    {
#ifdef DEBUG
        assert(row < nrow);
        assert(column < ncol);
#endif
    }
    
    double SparseMatrix::norm(int norm) const{
#ifdef DEBUG
        assertHasSparse();
#endif
       return cholmod_norm_sparse(sparse, norm, ConfigSingleton::getCommonPtr());
    }
    
    SparseMatrix SparseMatrix::copy() const{
        SparseMatrix res;
        if (sparse){
            res.sparse = cholmod_copy_sparse(sparse, ConfigSingleton::getCommonPtr());
            res.values = ((double*)res.sparse->x);
            res.iRow = ((int*)res.sparse->i);
            res.jColumn = ((int*)res.sparse->p);
        }
        if (triplet){
            res.triplet = cholmod_copy_triplet(triplet, ConfigSingleton::getCommonPtr());
            res.values = (double *)res.triplet->x;
            res.iRow = (int *)res.triplet->i;
            res.jColumn = (int *)res.triplet->j;
        }
        res.nrow = nrow;
        res.ncol = ncol;
        res.symmetry = symmetry;
        res.maxTripletElements = maxTripletElements;
        return move(res);
    }
    
    void SparseMatrix::assertHasSparse() const
    {
#ifdef DEBUG
        assert(sparse != nullptr); // matrix must be build
#endif
    }
    
    void SparseMatrix::increaseTripletCapacity(){
        // grow size with factor 1.5
        maxTripletElements = (int)ceil(1.5f*maxTripletElements);
        auto newTriplet = cholmod_allocate_triplet(nrow, ncol, maxTripletElements, symmetry, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        memcpy(newTriplet->x, triplet->x, triplet->nzmax*sizeof(double));
        memcpy(newTriplet->i, triplet->i, triplet->nzmax*sizeof(int));
        memcpy(newTriplet->j, triplet->j, triplet->nzmax*sizeof(int));
        newTriplet->nnz = triplet->nnz;
        cholmod_free_triplet(&triplet, ConfigSingleton::getCommonPtr());
        triplet = newTriplet;
        values = (double *)triplet->x;
        iRow = (int *)triplet->i;
        jColumn = (int *)triplet->j;
    }
    
    void SparseMatrix::assertValidInitAddValue(unsigned int row, unsigned int column) const {
#ifdef DEBUG
        assert(sparse == nullptr); // must be called before matrix build
        assert(row < nrow);
        assert(column < ncol);
#endif
    }
    
    
    SparseMatrix operator+(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double scale[2] = {1.,1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, scale, scale, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator+(SparseMatrix&& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double scale[2] = {1.,1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, scale, scale, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&LHS.sparse, ConfigSingleton::getCommonPtr());
        LHS.sparse = sparse;
        LHS.values = ((double*)sparse->x);
        LHS.iRow = ((int*)sparse->i);
        LHS.jColumn = ((int*)sparse->p);
        return move(LHS);
    }
    
    SparseMatrix&& operator+(const SparseMatrix& LHS, SparseMatrix&& RHS)
    {
        return move(RHS) + LHS;
    }
    
    SparseMatrix&& operator+(SparseMatrix&& LHS, SparseMatrix&& RHS)
    {
        return move(LHS) + RHS;
    }
    
    SparseMatrix operator-(const SparseMatrix& M)
    {
        return -1.*M;
    }
    
    SparseMatrix&& operator-(SparseMatrix&& M)
    {
        return -1.*move(M);
    }
    
    SparseMatrix operator-(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double alpha[2] = {1.,1.};
        double beta[2] = {-1.,-1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, alpha, beta, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator-(SparseMatrix&& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double alpha[2] = {1.,1.};
        double beta[2] = {-1.,-1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, alpha, beta, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&LHS.sparse, ConfigSingleton::getCommonPtr());
        LHS.sparse = sparse;
        LHS.values = ((double*)sparse->x);
        LHS.iRow = ((int*)sparse->i);
        LHS.jColumn = ((int*)sparse->p);
        return move(LHS);
    }
    
    SparseMatrix&& operator-(const SparseMatrix& LHS, SparseMatrix&& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.nrow == RHS.nrow && LHS.ncol == RHS.ncol);
        double alpha[2] = {1.,1.};
        double beta[2] = {-1.,-1.};
        cholmod_sparse *sparse = cholmod_add(LHS.sparse, RHS.sparse, alpha, beta, true, true, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&RHS.sparse, ConfigSingleton::getCommonPtr());
        RHS.sparse = sparse;
        RHS.values = ((double*)sparse->x);
        RHS.iRow = ((int*)sparse->i);
        RHS.jColumn = ((int*)sparse->p);
        return move(RHS);
    }
    
    SparseMatrix&& operator-(SparseMatrix&& LHS, SparseMatrix&& RHS)
    {
        return move(LHS) - RHS;
    }
    
    SparseMatrix operator*(const SparseMatrix& LHS, const SparseMatrix& RHS)
    {
        assert(LHS.sparse && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
        cholmod_sparse *sparse = cholmod_ssmult(LHS.sparse, RHS.sparse, ASYMMETRIC, true, true, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix operator*(const SparseMatrix& LHS, double RHS)
    {
        assert(LHS.sparse);
        cholmod_dense *dense = cholmod_zeros(1, 1, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        ((double*)dense->x)[0] = RHS;
        cholmod_sparse *sparse = cholmod_copy_sparse(LHS.sparse, ConfigSingleton::getCommonPtr());
        cholmod_scale(dense, CHOLMOD_SCALAR, sparse, ConfigSingleton::getCommonPtr());
        cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& operator*(SparseMatrix&& LHS, double RHS)
    {
        assert(LHS.sparse);
        cholmod_dense *dense = cholmod_zeros(1, 1, CHOLMOD_REAL, ConfigSingleton::getCommonPtr());
        ((double*)dense->x)[0] = RHS;
        cholmod_scale(dense, CHOLMOD_SCALAR, LHS.sparse, ConfigSingleton::getCommonPtr());
        cholmod_free_dense(&dense, ConfigSingleton::getCommonPtr());
        return move(LHS);
    }
    
    SparseMatrix operator*(double LHS, const SparseMatrix& RHS)
    {
        return RHS * LHS;
    }
    
    SparseMatrix&& operator*(double LHS, SparseMatrix&& RHS)
    {
        return move(RHS) * LHS;
    }
    
    // Computes result = alpha*(A*X) + beta*Y (where this == A) or Y = alpha*(A'*X) + beta*Y when transpose is true
    // transpose is ignores if matrix is symmetric of Hermitian
    void SparseMatrix::multiply(bool transpose, double alpha, double beta, const DenseMatrix& X, DenseMatrix& result){
#ifdef DEBUG
        assert(this->symmetry == Symmetry::ASYMMETRIC || !transpose);
        if (!transpose){
            assert(getRows() == result.getRows());
            assert(X.getColumns() == result.getColumns());
            assert(getColumns() == X.getRows());
        } else {
            assert(getColumns() == result.getRows());
            assert(X.getColumns() == result.getColumns());
            assert(getRows() == X.getRows());
        }
#endif
        cholmod_sdmult(sparse, transpose, &alpha, &beta, X.dense, result.dense, ConfigSingleton::getCommonPtr());
    }
    
    DenseMatrix SparseMatrix::multiply(bool transpose, double alpha, const DenseMatrix& X){
        unsigned int rows = transpose?getColumns():getRows();
        DenseMatrix res{rows, X.getColumns(), 0.0};
        multiply(transpose, alpha, 0, X, res);
        return res;
    }
    
    
    DenseMatrix operator*(const DenseMatrix& LHS, const SparseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.dense && RHS.sparse);
        assert(LHS.ncol == RHS.nrow);
#endif
        double alpha[2] = {1.,1.};
        double beta[2] = {0.,0.};
        DenseMatrix res(RHS.ncol, LHS.nrow);
        cholmod_sdmult(RHS.sparse, true, alpha, beta, transposed(LHS).dense, res.dense, ConfigSingleton::getCommonPtr());
        return transposed(res);
    }
    
    DenseMatrix operator*(const SparseMatrix& LHS, const DenseMatrix& RHS)
    {
#ifdef DEBUG
        assert(LHS.sparse && RHS.dense);
        assert(LHS.ncol == RHS.nrow);
#endif
        double alpha[2] = {1.,1.};
        double beta[2] = {0.,0.};
        DenseMatrix res(LHS.nrow, RHS.ncol);
        cholmod_sdmult(LHS.sparse, false, alpha, beta, RHS.dense, res.dense, ConfigSingleton::getCommonPtr());
        return res;
    }
    
    void SparseMatrix::transpose()
    {
        assert(symmetry == ASYMMETRIC);
        sparse = cholmod_transpose(sparse, 1, ConfigSingleton::getCommonPtr());
        nrow = static_cast<int>(sparse->nrow);
        ncol = static_cast<int>(sparse->ncol);
        values = ((double*)sparse->x);
        iRow = ((int*)sparse->i);
        jColumn = ((int*)sparse->p);
    }
    
    SparseMatrix transposed(const SparseMatrix& M)
    {
        assert(M.symmetry == ASYMMETRIC);
        cholmod_sparse *sparse = cholmod_transpose(M.sparse, 1, ConfigSingleton::getCommonPtr());
        return SparseMatrix(sparse);
    }
    
    SparseMatrix&& transposed(SparseMatrix&& M)
    {
        assert(M.symmetry == ASYMMETRIC);
        cholmod_sparse *sparse = cholmod_transpose(M.sparse, 1, ConfigSingleton::getCommonPtr());
        cholmod_free_sparse(&M.sparse, ConfigSingleton::getCommonPtr());
        M.sparse = sparse;
        M.nrow = static_cast<int>(sparse->nrow);
        M.ncol = static_cast<int>(sparse->ncol);
        M.values = ((double*)sparse->x);
        M.iRow = ((int*)sparse->i);
        M.jColumn = ((int*)sparse->p);
        return move(M);
    }
    
    DenseMatrix solve(const SparseMatrix& A, const DenseMatrix& b)
    {
#ifdef DEBUG
        assert(A.sparse && b.dense);
        assert(A.nrow == b.nrow);
#endif
        Factor F = A.analyze();
        F.factorize(A);
        return solve(F, b);
    }
    
    SparseMatrix solve(const SparseMatrix& A, const SparseMatrix& b)
    {
#ifdef DEBUG
        assert(A.sparse && b.sparse);
        assert(A.nrow == b.nrow);
#endif
        Factor F = A.analyze();
        F.factorize(A);
        return solve(F, b);
    }
    
    ostream& operator<<(ostream& os, const SparseMatrix& A)
    {
        os << endl;
        if (A.sparse){
            cholmod_print_sparse(A.sparse, "", ConfigSingleton::getCommonPtr());
            DenseMatrix Dense = A.toDense();
            os << "[";
            for (int r = 0; r < Dense.getRows(); r++)
            {
                for (int c = 0; c < Dense.getColumns(); c++)
                {
                    os << Dense(r,c);
                    if(c < Dense.getColumns()-1) {
                        os << ", ";
                    }
                }
                if(r < Dense.getRows()-1)
                {
                    os << ";" << endl << " ";
                }
            }
            os << "];" << endl;
        }
        else if (A.triplet){
            cholmod_print_triplet(A.triplet, "", ConfigSingleton::getCommonPtr());
            for (int i = 0; i < A.triplet->nnz; i++){
                os << "( " << ((int*)A.triplet->i)[i] << ", " << ((int*)A.triplet->j)[i]<<" ) =\t"<< ((double*)A.triplet->x)[i] << endl;
            }
        } else {
            os << "[Empty sparse matrix]" << endl;
        }
        return os;
    }
}
