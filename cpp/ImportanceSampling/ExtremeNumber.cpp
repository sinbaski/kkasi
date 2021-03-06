#include <assert.h>
#include <climits>
#include <stdio.h>
#include "ExtremeNumber.hpp"
#include <iterator>

using namespace std;
using namespace arma;

inline long bottom (double x)
{
    return x >= 0 ? (long)x : (long)floor(x);
}

ExtremeNumber::ExtremeNumber(void)
    :logs()
{
}

ExtremeNumber::ExtremeNumber(double x)
    :logs()
{
    assert(x >= 0);
    if (x > 0)
	logs.push_back(log10(x));
}

ExtremeNumber::ExtremeNumber(const container_t<double>& logs)
    : logs(logs)
{
}

ExtremeNumber::ExtremeNumber(const ExtremeNumber& x)
    : logs(x.logs)
{
}

void adjust_numbers(container_t<double>::iterator p, bool checkall, ExtremeNumber &V)
{
    if (V.logs.size() < 2) return;
    auto i = p;
    while (distance(i, V.logs.end()) >= 2) {
	long f = bottom(*i);
	if (bottom(*next(i)) > f) {
	    if (!checkall) return;
	    i++;
	    continue;
	}
	double r = pow(10, *i - f) + pow(10, *next(i) - f);
	r = log10(r);
	*i = r + f;
	auto j = V.logs.erase(next(i));
	if (r < 1) {
	    if (checkall) i = j;
	    else return;
	}
    }
}

void add_to (container_t<double>::const_iterator i1,
	     ExtremeNumber& V)
{
    long f = bottom(*i1);
    container_t<double>::iterator k =
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
	V.logs.insert(next(k), *i1);
	return;
    }
    double r = log10(pow(10, *k - f) + pow(10, *i1 - f));
    *k = r + f; // bottom(*k) in [f, f + 1]
    adjust_numbers(k, false, V);
}


ExtremeNumber multiply_to (container_t<double>::const_iterator i,
			 const ExtremeNumber& X)
{
    ExtremeNumber V(X);
    for (container_t<double>::iterator j = V.logs.begin();
	 distance(j, V.logs.end()); j++) {
	*j += *i;
    }
    adjust_numbers(V.logs.begin(), true, V);
    return V;
}

const ExtremeNumber& ExtremeNumber::erase_small(void)
{
    long f = bottom(*logs.rbegin()) - 20;
    container_t<double>::iterator k =
	upper_bound(logs.begin(), logs.end(), f,
		    [](int b, double x) {
			return b < bottom(x);
		    });
    if (k == logs.begin()) return *this;
    k--;
    if (bottom(*k) == f) {
	if (k != logs.begin())
	    logs.erase(logs.begin(), k);
    } else {
	logs.erase(logs.begin(), next(k));
    }
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator = (double y)
{
    assert(y >= 0);
    if (y == 0)
	logs.clear();
    else 
	logs.push_back(log10(y));
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator = (const ExtremeNumber &y) {
    logs.assign(y.logs.begin(), y.logs.end());
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator *= (double y)
{
    assert(y >= 0);

    if (y == 0) {
	logs.clear();
	return *this;
    }
    
    ExtremeNumber Y(y);
    operator= (multiply_to(Y.logs.begin(), *this));
    erase_small();
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator /= (double y)
{
    return operator *= (1/y);
}

ExtremeNumber operator * (const ExtremeNumber& x, double y)
{
    assert(y >= 0);

    ExtremeNumber z(x);
    z *= y;
    return z;
}

ExtremeNumber operator / (const ExtremeNumber& x, double y)
{
    return x * (1/y);
}
    
ExtremeNumber operator * (const ExtremeNumber& x, const ExtremeNumber& y)
{
    if (x.logs.empty() || y.logs.empty()) return ExtremeNumber();
    ExtremeNumber z;
    for (container_t<double>::const_iterator i = y.logs.begin();
	 distance(i, y.logs.end()); i++) {
	z += multiply_to(i, x);
    }
    z.erase_small();
    return z;
}

const ExtremeNumber& ExtremeNumber::operator *= (const ExtremeNumber& y)
{
    if (logs.empty() || y.logs.empty()) {
	logs.clear();
	return *this;
    }

    operator= (*this * y);
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator += (const ExtremeNumber& y)
{
    for (container_t<double>::const_iterator i = y.logs.begin();
	 distance(i, y.logs.end()); i++) {
	add_to(i, *this);
    }
    erase_small();
    return *this;
}
    
ExtremeNumber operator + (const ExtremeNumber& x, const ExtremeNumber& y)
{
    ExtremeNumber z(x);
    z += y;
    return z;
}

const ExtremeNumber& ExtremeNumber::operator += (double y)
{
    assert(y >= 0);
    
    if (y == 0) return *this;
    ExtremeNumber Y(y);
    return *this += Y;
}
    
ExtremeNumber operator + (const ExtremeNumber& x, double y)
{
    ExtremeNumber z(x);
    z += y;
    return z;
}

bool ExtremeNumber::operator < (const ExtremeNumber& x) const
{
    auto i = logs.crbegin(), j = x.logs.crbegin();
    while (distance(i, logs.crend()) && distance(j, x.logs.crend())) {
	if (*i < *j) return true;
	if (*i > *j) return false;
	i++, j++;
    }
    return distance(j, x.logs.crend());
}

long ExtremeNumber::leading_power(void) const
{
    return bottom(*logs.crbegin());
}
    
double ExtremeNumber::comptify(long power) const
{
     return
	accumulate(logs.begin(), logs.end(), 0.0,
		   [power](double s, double y) {
		       return s + pow(10, y - power);
		   });
}

double ExtremeNumber::comptify(void) const
{
    return comptify(0);
}

double log10(const ExtremeNumber& x)
{
    long p = x.leading_power();
    double y = x.comptify(p);
    return p + log10(y);
}

double log(const ExtremeNumber& x)
{
    return log(10) * log10(x);
}

ExtremeNumber ExtremeNumber::operator ^ (unsigned long n) const
{
    ExtremeNumber X(*this);
    for (unsigned i = 0; i < n; i++) {
	X *= *this;
    }
    return X;
}

double ExtremeNumber::nth_root (unsigned long n) const
{
    long p = leading_power();
    long q = (double)n;
    double x = comptify(p);
    return pow(x, 1/q) * pow(10, p/q);
}

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

ExtremeNumber& XMatrix::operator () (unsigned i, unsigned j)
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

#if !defined DEBUG
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

#if !defined DEBUG
#pragma omp parallel for
#endif
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
    *power = INT_MAX;
    for (unsigned i = 0; i < entry.size(); i++) {
	for (unsigned j = 0; j < entry[0].size(); j++) {
	    *power = min(*power, entry[i][j].leading_power());
	}
    }
    mat X(entry.size(), entry[0].size());
    double p = *power;
    for (unsigned i = 0; i < entry.size(); i++) {
	for (unsigned j = 0; j < entry[0].size(); j++) {
	    double a = 0;
	    for_each(entry[i][j].logs.begin(),
		     entry[i][j].logs.end(),
		     [p, &a](double d){
			 a += pow(10, d - p);
		     });
	    X(i, j) = a;
	}
    }
    return X;
}
