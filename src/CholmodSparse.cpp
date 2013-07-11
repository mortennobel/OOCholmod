//
//  CholmodSparse.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 1/21/13.
//  Copyright (c) 2013 DTU Compute. All rights reserved.
//  License: LGPL 3.0 

#include "CholmodSparse.h"
#include <cassert>

#include "CholmodDenseVector.h"


using namespace std;

// bad coffee odd food
#define MAGIC_NUMBER (unsigned long)0xBADC0FFEE0DDF00DL

CholmodSparse::CholmodSparse(unsigned int nrow, unsigned int ncol, cholmod_common *Common, int maxSize)
:Common(Common), sparse(NULL), nrow(nrow), ncol(ncol)
#ifdef DEBUG
,magicNumber(MAGIC_NUMBER)
#endif
{
    int elements;
    if (maxSize == 0) {
        elements = (nrow*(ncol+1))/2; // triangular number
    } else {
        elements = maxSize;
    }
    triplet = cholmod_allocate_triplet(nrow, ncol, elements, (int)SYMMETRIC_UPPER, CHOLMOD_REAL, Common);
    values = (double *)triplet->x;
    iRow = (int *)triplet->i;
	jColumn = (int *)triplet->j;
    this->symmetry = SYMMETRIC_UPPER;
#ifdef DEBUG
    assert(nrow == ncol); // must be square
    maxElements = elements;
#endif
}

CholmodSparse::CholmodSparse(cholmod_sparse *sparse, cholmod_common *Common)
:sparse(sparse), Common(Common),  nrow((unsigned int)sparse->nrow), ncol((unsigned int)sparse->ncol)
#ifdef DEBUG
,magicNumber(MAGIC_NUMBER)
#endif
{
    buildLookupIndexFromSparse();
}

CholmodSparse::CholmodSparse(CholmodSparse&& move)
:Common(move.Common), sparse(move.sparse), triplet(move.triplet), nrow(move.nrow), ncol(move.ncol), lookupIndex(std::move(move.lookupIndex)), iRow(move.iRow), jColumn(move.jColumn), symmetry(move.symmetry)
#ifdef DEBUG
    ,magicNumber(move.magicNumber), maxElements(move.maxElements)
#endif
{
    move.Common = nullptr;
    move.sparse = nullptr;
    move.triplet = nullptr;
    move.values = nullptr;
    move.iRow = nullptr;
    move.jColumn = nullptr;
#ifdef DEBUG
    move.magicNumber = 0L;
    move.maxElements = 0;
#endif
}

CholmodSparse& CholmodSparse::operator=(CholmodSparse&& other){
    if (this != &other){
        if (sparse != NULL){
            cholmod_free_sparse(&sparse, Common);
        }
        if (triplet != NULL){
            cholmod_free_triplet(&triplet, Common);
        }
        
        Common = other.Common;
        sparse = other.sparse;
        triplet = other.triplet;
        nrow = other.nrow;
        ncol = other.ncol;
        lookupIndex = std::move(other.lookupIndex);
        iRow = other.iRow;
        jColumn = other.jColumn;
        symmetry = other.symmetry;
#ifdef DEBUG
        magicNumber = other.magicNumber;
        maxElements = other.maxElements;
#endif
    }
    return *this;
}


CholmodSparse::~CholmodSparse(){
    if (Common != nullptr){
#ifdef DEBUG
        assert(magicNumber == MAGIC_NUMBER);
        magicNumber = 0;
#endif
        if (sparse != NULL){
            cholmod_free_sparse(&sparse, Common);
            sparse = NULL;
        }
        if (triplet != NULL){
            cholmod_free_triplet(&triplet, Common);
            triplet = NULL;
        }
    }
}

void CholmodSparse::setNullSpace(CholmodDenseVector *N){
#ifdef DEBUG
    assert(sparse != NULL);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    // naive implementation: Todo run fast
    CholmodDenseVector &v = *N;
    int idx = 0;
    for (int j=0;j<ncol;j++){
        int iFrom = ((int*)sparse->p)[j];
        int iTo = ((int*)sparse->p)[j+1]-1;
        for (int i=iFrom;i<=iTo;i++){
            int row = ((int*)sparse->i)[i];
            ((double*)sparse->x)[idx] *= v[row]*v[j];
            idx++;
        }
    }
    for (int i=0;i<v.getSize();i++){
        if (v[i] == 0){
            this->setValue(i,i,1);
        }
    }
}

void CholmodSparse::setNullSpace(CholmodDenseVector& N){
    setNullSpace(&N);
}

// computes this * X and store the result in res
void CholmodSparse::multiply(CholmodDenseVector *X, CholmodDenseVector *res, double alpha, double beta){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
    assert(res != NULL);
    assert(X->getSize() == res->getSize());
    assert(X != res);
    assert(nrow == X->getSize());
#endif    
    
    // alpha*(A*X) + beta*Y
    double _alpha[2] = {alpha,alpha};
    double _beta[2] = {beta,beta};
    // int cholmod_sdmult(cholmod_sparse *A, ￼￼int transpose, double alpha [2], double beta [2], cholmod_dense *X, cholmod_dense *Y, cholmod_common *Common );
    cholmod_sdmult(sparse, false, _alpha, _beta, X->getHandle(), res->getHandle(), Common);
}

CholmodDenseVector CholmodSparse::multiply(CholmodDenseVector& X, double alpha, double beta){
    CholmodDenseVector res(X.getSize(), Common);
    multiply(X, res, alpha, beta);
    return res;
}

void CholmodSparse::multiply(CholmodDenseVector& X, CholmodDenseVector& res, double alpha, double beta){
    multiply(&X, &res, alpha, beta);
}

void CholmodSparse::build(bool readOnly){
#ifdef DEBUG
    assert(sparse == NULL);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    sparse = cholmod_triplet_to_sparse(triplet, triplet->nnz, Common);
    cholmod_free_triplet(&triplet, Common);
    triplet = NULL;
    values = NULL;
    iRow = NULL;
	jColumn = NULL;
    
    // build lookup index
#ifdef DEBUG
    assert(sparse->stype == symmetry);
    assert(sparse->packed);
#endif
    if (!readOnly){
        buildLookupIndexFromSparse();
    }
}

void CholmodSparse::buildLookupIndexFromSparse(){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    lookupIndex.clear();
    // In packed form, the nonzero pattern of column j is in A->i [A->p [j] ... A->p [j+1]-1]
    int idx = 0;
    for (int j=0;j<ncol;j++){
        int iFrom = ((int*)sparse->p)[j];
        int iTo = ((int*)sparse->p)[j+1]-1;
        for (int i=iFrom;i<=iTo;i++){
            int row = ((int*)sparse->i)[i];
            lookupIndex[key(row, j)] = idx;
            idx++;
        }
    }
}

CholmodFactor *CholmodSparse::analyzePtr(){
#ifdef DEBUG
    assert(sparse != NULL);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_factor *L = cholmod_analyze(sparse, Common);
    return new CholmodFactor(L, Common);
}

CholmodFactor CholmodSparse::analyze(){
#ifdef DEBUG
    assert(sparse != NULL);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    cholmod_factor *L = cholmod_analyze(sparse, Common);
    return CholmodFactor(L, Common);
}


void CholmodSparse::zero(){
#ifdef DEBUG
    assert(sparse != NULL);
    assert(magicNumber == MAGIC_NUMBER);
#endif
    memset(sparse->x, 0, sparse->nzmax * sizeof(double));
}

void CholmodSparse::print(const char* name){
#ifdef DEBUG
    assert(magicNumber == MAGIC_NUMBER);
#endif
    if (sparse){
        cholmod_print_sparse(sparse, name, Common);
        cholmod_dense *dense = cholmod_sparse_to_dense(sparse, Common);
        int n_rows = (int)dense->nrow;
        int n_cols = (int)dense->ncol;
        for (int r = 0; r  < n_rows; r++)
        {
            for (int c = 0; c  < n_cols; c++)
            {
                cout << ((double*)dense->x)[c*n_rows + r] << " ";
            }
            cout << endl;
        }
        cholmod_free_dense(&dense, Common);
        cout << endl;
        cout << "Packed "<<sparse->packed<< endl;
        cout << "p: ";
        for (int i=0;i<=sparse->ncol;i++){
            printf("%4i ", ((int*)sparse->p)[i]);
        }
        cout << endl;
        cout << "i: ";
        for (int i=0;i<sparse->nzmax;i++){
            printf("%4i ", ((int*)sparse->i)[i]);
        }
        cout << endl;
        cout << "x: ";
        for (int i=0;i<sparse->nzmax;i++){
            printf("%3.3f ", ((double*)sparse->x)[i]);
        }
        cholmod_free_dense(&dense, Common);
        cout << endl;
    }
}

void CholmodSparse::assertValidIndex(unsigned int row, unsigned int column){
#ifdef DEBUG
    assert(symmetry == SYMMETRIC_UPPER);
    assert(row < nrow);
    assert(column < ncol);
    if (symmetry == SYMMETRIC_UPPER) {
        assert(row <= column);
    } else if (symmetry == SYMMETRIC_LOWER) {
        assert(row >= column);
    }
    
#endif
}

void CholmodSparse::assertHasSparse(){
#ifdef DEBUG
    assert(sparse != NULL); // matrix must be build
#endif
}

void CholmodSparse::assertValidInitAddValue(unsigned int row, unsigned int column, double value){
#ifdef DEBUG
    assert(sparse == NULL); // must be called before matrix build
    assert(triplet->nnz < maxElements);
    assert(row < nrow);
    assert(column < ncol);
    if (symmetry == SYMMETRIC_UPPER) {
        assert(row <= column);
    } else if (symmetry == SYMMETRIC_LOWER) {
        assert(row >= column);
    }
    int shiftBits = sizeof(long)*8/2; // shift half of the bits of a long
    long maxId = (long)pow(2, shiftBits);
    assert (row < maxId);
    assert (column < maxId);
#endif
}