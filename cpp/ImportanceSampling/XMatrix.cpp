#include <cassert>
#include <vector>
#include "ExtremeNumber2.hpp"
#include "XMatrix.hpp"

XMatrix::XMatrix(void)
    :entry()
{
}

XMatrix::XMatrix(unsigned m, unsigned n)
    :entry(m)
{
    for_each(entry.begin(), entry.end(),
	     [n](vector<ExtremeNumber>& row) {
		 row.resize(n);
	     });
}

XMatrix::XMatrix(double** data, unsigned m, unsigned n)
    :entry(m)
{
    for (unsigned i = 0; i < m; i++) {
	entry[i].assign(data[i], data[i] + n);
    }
}


XMatrix::XMatrix(const XMatrix& M)
{
    // entry.resize(M.entry.size());
    // copy(M.entry.begin(), M.entry.end(), entry.begin());
    operator= (M);
}

ExtremeNumber& XMatrix::operator() (unsigned i, unsigned j)
{
    assert(entry.size() > i);
    for_each(entry.begin(), entry.end(),
	     [j](vector<ExtremeNumber> &row) {
		 assert(row.size() > j);
	     });
    return entry[i][j];
}

const ExtremeNumber& XMatrix::operator() (unsigned i, unsigned j) const
{
    return entry[i][j];
}

const XMatrix& XMatrix::operator = (const XMatrix& Y)
{
    entry.resize(Y.entry.size());
    copy(Y.entry.begin(), Y.entry.end(), entry.begin());
    return *this;
}

XMatrix operator * (const XMatrix& X, const XMatrix& Y)
{
    unsigned n = Y.entry.size();
    assert(X.entry[0].size() == n);

    XMatrix Z(X.entry.size(), Y.entry[0].size());

#ifndef DEBUG
#pragma omp parallel for
#endif
    for (unsigned i = 0; i < Z.entry.size(); i++) {
	for (unsigned j = 0; j < Z.entry[0].size(); j++) {
	    Z.entry[i][j] = 0;
	    for (unsigned k = 0; k < n; k++) {
		Z.entry[i][j] += X(i, k) * Y(k, j);
	    }
	}
    }
    return Z;
}

const XMatrix& XMatrix::operator *= (const XMatrix& Y)
{
    unsigned n = Y.entry.size();
    assert(entry[0].size() == n);
    operator= (*this * Y);
    return *this;
}

const XMatrix& XMatrix::operator += (const XMatrix& Y)
{
    assert(entry.size() == Y.entry.size());
    assert(entry[0].size() == Y.entry[0].size());

#pragma omp parallel for
    for (unsigned i = 0; i < entry.size(); i++) {
	for (unsigned j = 0; j < entry[0].size(); j++) {
	    entry[i][j] += Y.entry[i][j];
	}
    }
    return *this;
}

XMatrix operator + (const XMatrix& X, const XMatrix& Y)
{
    XMatrix Z(X);
    Z += Y;
    return Z;
}

const XMatrix& XMatrix::operator ^ (unsigned n)
{
    if (entry.empty() || n == 1) return *this;
    if (n == 0) {
	for (unsigned i = 0; i < entry.size(); i++) {
	    entry[i][i] = 1;
	}
    }
    XMatrix X(*this);
    for (unsigned i = 1; i < n; i++) {
//	int p;
	operator *= (X);
	// mat V = comptify(&p);
	// printf("At %d power: Order: %d\n", i + 1, p);
	// V.print();
	// printf("\n");
    }
    return *this;
}

XMatrix XMatrix::operator ~ (void)
{
    XMatrix X(*this);
    for (unsigned i = 0; i < X.entry.size(); i++) {
	for (unsigned j = i + 1; j < X.entry[0].size(); j++) {
	    swap(X(i, j), X(j, i));
	}
    }
    return X;
}

mat XMatrix::comptify(long *power)
{
    *power = LONG_MAX;
    for (unsigned i = 0; i < entry.size(); i++) {
	for (unsigned j = 0; j < entry[0].size(); j++) {
	    *power = min(*power, entry[i][j].power());
	}
    }
    mat X(entry.size(), entry[0].size());
    for (unsigned i = 0; i < entry.size(); i++) {
	for (unsigned j = 0; j < entry[0].size(); j++) {
	    X(i, j) = entry[i][j](*power);
	}
    }
    return X;
}
