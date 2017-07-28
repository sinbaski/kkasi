#ifndef XMATRIX_HPP
#define XMATRIX_HPP

#include <armadillo>
#include <vector>
#include "ExtremeNumber2.hpp"

using namespace std;
using namespace arma;

class XMatrix
{
public:
   
    vector< vector<ExtremeNumber> > entry;

public:

    XMatrix(void);
    XMatrix(const XMatrix& M);
    XMatrix(const mat& M);
    XMatrix(unsigned int m, unsigned int n);
    XMatrix(double** data, unsigned int m, unsigned int n);
    ExtremeNumber& operator() (unsigned i, unsigned j);
    const ExtremeNumber& operator() (unsigned i, unsigned j) const;
    const XMatrix& operator *= (const XMatrix& Y);
    const XMatrix& operator += (const XMatrix& Y);
    const XMatrix& operator = (const XMatrix& Y);
    const XMatrix& operator ^ (unsigned n);
    void set_size(unsigned int m, unsigned int n);
    XMatrix operator ~ (void);
    mat comptify(long *power);
};

XMatrix operator * (const XMatrix& X, const XMatrix& Y);
XMatrix operator + (const XMatrix& X, const XMatrix& Y);

#endif

