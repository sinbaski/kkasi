#include <assert.h>
#include <vector>
#include "SmallNumber.hpp"

using namespace std;
using namespace arma;

inline int bottom (double x)
{
    return x >= 0 ? (int)x : (int)floor(x);
}

SmallNumber::SmallNumber(void)
    :logs()
{
}

SmallNumber::SmallNumber(double x)
    :logs()
{
    assert(x >= 0);
    if (x > 0)
	logs.push_back(log10(x));
}

SmallNumber::SmallNumber(const vector<double>& logs)
    : logs(logs)
{
}

SmallNumber::SmallNumber(const SmallNumber& x)
    : logs(x.logs)
{
}

void carry(vector<double>::iterator o1,
	   SmallNumber &V)
{
    int f = bottom(*o1);
    double r = log10(pow(10, *o1 - f) + 1);
    if (r < 1) {
	*o1 = f + r;
    } else {
	*o1 = f + r - 1;
	vector<double>::iterator i = 
	    upper_bound(o1 + 1, V.logs.end(), f + 1,
			[](int y, double x) {
			    return y < bottom(x);
			});
	if (bottom(*(i - 1)) < f + 1) {
	    V.logs.insert(i - 1, f + 1);
	} else {
	    carry(i - 1, V);
	}
    }
}

void add_to (vector<double>::const_iterator i1,
	     SmallNumber& V)
{
    int f = bottom(*i1);
    vector<double>::iterator k =
	upper_bound(V.logs.begin(), V.logs.end(), f,
		    [](int y, double x){
			return y < bottom(x);
		    });
    if (k == V.logs.begin()) {
	V.logs.insert(k, *i1);
	return;
    }
    k--;
    if (bottom(*k) < f) {
	V.logs.insert(k + 1, *i1);
	return;
    }
    double r = log10(pow(10, *k - f) + pow(10, *i1 - f));
    if (r >= 1) {
	*k = f + r - 1;
	if (k + 1 == V.logs.end() || bottom(*(k + 1)) - f > 1)
	    V.logs.insert(k + 1, f + 1);
	else
	    carry(k + 1, V);
    } else {
	*k = f + r;
    }
}

void adjust_numbers(SmallNumber &V)
{
    if (V.logs.size() < 2) return;
    for (vector<double>::iterator i = V.logs.begin();
	 i < V.logs.end() - 1; ) {
	int f = bottom(*i);
	if (bottom(*(i + 1)) > f) {
	    i++;
	    continue;
	}
	double r = pow(10, *i - f) + pow(10, *(i + 1) - f);
	r = log10(r);
	if (r < 1) {
	    *i = f + r;
	    i = V.logs.erase(i + 1);
	} else {
	    *i = f + r - 1;
	    *(i + 1) = f + 1;
	    i++;
	}
    }
}

SmallNumber multiply_to (vector<double>::const_iterator i,
			 const SmallNumber& X)
{
    SmallNumber V(X);
    for (vector<double>::iterator j = V.logs.begin();
	   j < V.logs.end(); j++) {
	*j += *i;
    }
    adjust_numbers(V);
    return V;
}

const SmallNumber& SmallNumber::erase_small(void)
{
    int f = bottom(*logs.rbegin()) - 20;
    vector<double>::const_iterator k =
	upper_bound(logs.begin(), logs.end(), f);
    if (k == logs.begin()) return *this;
    k--;
    if (bottom(*k) == f) {
	if (k != logs.begin())
	    logs.erase(logs.begin(), k - 1);
    } else {
	logs.erase(logs.begin(), k);
    }
    return *this;
}

const SmallNumber& SmallNumber::operator = (double y)
{
    assert(y >= 0);
    if (y == 0)
	logs.clear();
    else 
	logs.push_back(log10(y));
    return *this;
}

const SmallNumber& SmallNumber::operator = (const SmallNumber &y) {
    logs.assign(y.logs.begin(), y.logs.end());
    return *this;
}

const SmallNumber& SmallNumber::operator *= (double y)
{
    assert(y >= 0);

    if (y == 0) {
	logs.clear();
	return *this;
    }
    
    SmallNumber Y(y);
    operator= (multiply_to(Y.logs.begin(), *this));
    erase_small();
    return *this;
}

SmallNumber operator * (const SmallNumber& x, double y)
{
    assert(y >= 0);

    SmallNumber z(x);
    z *= y;
    return z;
}
    
SmallNumber operator * (const SmallNumber& x, const SmallNumber& y)
{
    if (x.logs.empty() || y.logs.empty()) return SmallNumber();
    SmallNumber z;
    for (vector<double>::const_iterator i = y.logs.begin();
	 i < y.logs.end(); i++) {
	z += multiply_to(i, x);
    }
    z.erase_small();
    return z;
}

const SmallNumber& SmallNumber::operator *= (const SmallNumber& y)
{
    if (logs.empty() || y.logs.empty()) {
	logs.clear();
	return *this;
    }

    operator= (*this * y);
    return *this;
}

const SmallNumber& SmallNumber::operator += (const SmallNumber& y)
{
    for (vector<double>::const_iterator i = y.logs.begin();
	 i < y.logs.end(); i++) {
	add_to(i, *this);
    }
    erase_small();
    return *this;
}
    
SmallNumber operator + (const SmallNumber& x, const SmallNumber& y)
{
    SmallNumber z(x);
    z += y;
    return z;
}

const SmallNumber& SmallNumber::operator += (double y)
{
    assert(y >= 0);
    
    if (y == 0) return *this;
    SmallNumber Y(y);
    return *this += Y;
}
    
SmallNumber operator + (const SmallNumber& x, double y)
{
    SmallNumber z(x);
    z += y;
    return z;
}
    
XMatrix::XMatrix(void)
    :entry()
{
}

XMatrix::XMatrix(unsigned m, unsigned n)
    :entry(m)
{
    for_each(entry.begin(), entry.end(),
	     [n](vector<SmallNumber>& row) {
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

SmallNumber& XMatrix::operator() (unsigned i, unsigned j)
{
    assert(entry.size() > i);
    for_each(entry.begin(), entry.end(),
	     [j](vector<SmallNumber> &row) {
		 assert(row.size() > j);
	     });
    return entry[i][j];
}

const SmallNumber& XMatrix::operator() (unsigned i, unsigned j) const
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

// #pragma omp parallel for
    for (unsigned i = 0; i < Z.entry.size(); i++) {
	for (unsigned j = 0; j < Z.entry[0].size(); j++) {
	    Z.entry[i][j] = 0;
	    for (unsigned k = 0; k < n; k++) {
		Z.entry[i][j] += X(i, k) * Y(k, j);
		double u = *Z.entry[i][j].logs.rbegin();
		;
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

// #pragma omp parallel for
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
	operator *= (X);
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
