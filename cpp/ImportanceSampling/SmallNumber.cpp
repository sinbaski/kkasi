#include <assert.h>
#include <vector>
#include "SmallNumber.hpp"

using namespace std;
using namespace arma;

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
    : logs(logs) {
}

SmallNumber::SmallNumber(const SmallNumber& x)
    : logs(x.logs) {
}

void carry(vector<double>::iterator o1,
	   vector<double>::iterator o2,
	   SmallNumber &V)
{
    int f = (int)(*o1);
    double r = log10(pow(10, *o1 - f) + 1);
    if (r < 1) {
	*o1 = f + r;
    } else {
	*o1 = f + r - 1;
	vector<double>::iterator i = 
	    upper_bound(o1 + 1, o2, f + 1,
			[](double x, int y) {
			    return (int)x < y;
			});
	if ((int)*(i - 1) < f + 1) {
	    V.logs.insert(i - 1, f + 1);
	} else {
	    carry(i - 1, o2, V);
	}
    }
}

void add_to (vector<double>::const_iterator i1,
	     SmallNumber& V)
{
    int f = (int)(*i1);
    vector<double>::iterator k =
	upper_bound(V.logs.begin(), V.logs.end(), f,
		    [](double x, int y){
			return (int)x < y;
		    });
    if (k == V.logs.end()) {
	V.logs.push_back(*i1);
	return;
    }
    if (k == V.logs.begin()) {
	V.logs.insert(k, *i1);
	return;
    }
    k--;
    if ((int)*k < f) {
	V.logs.insert(k + 1, *i1);
    } else {
	double r = log10(pow(10, *k - f) + pow(10, *i1 - f));
	if (r >= 1) {
	    *k = f + r - 1;
	    if ((int)*(k + 1) - f == 1)
		carry(k + 1, V.logs.end(), V);
	    else
		V.logs.insert(k + 1, f + 1);
	} else {
	    *k = f + r;
	}
    }
}

void adjust_numbers(SmallNumber &V)
{
    if (V.logs.size() < 2) return;
    for (vector<double>::iterator i = V.logs.begin();
	 i < V.logs.end() - 1; ) {
	int f = (int)(*i);
	if ((int)*(i + 1) > f) {
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

void multiply_to (vector<double>::const_iterator i,
		  SmallNumber& V)
{
    for (vector<double>::iterator j = V.logs.begin();
	   j < V.logs.end(); j++) {
	*j++ += *i;
    }
    adjust_numbers(V);
}

const SmallNumber& SmallNumber::operator= (double y)
{
    assert(y >= 0);
    if (y == 0)
	logs.clear();
    else 
	logs.push_back(log10(y));
    return *this;
}

const SmallNumber& SmallNumber::operator= (const SmallNumber &y) {
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
    multiply_to(Y.logs.begin(), *this);
    return *this;
}

SmallNumber operator* (const SmallNumber& x, double y)
{
    assert(y >= 0);

    SmallNumber z(x);
    z *= y;
    return z;
}
    
const SmallNumber& SmallNumber::operator *= (const SmallNumber& y)
{
    if (logs.empty() || y.logs.empty()) {
	logs.clear();
	return *this;
    }

    for (vector<double>::const_iterator i = y.logs.begin();
	 i < y.logs.end(); i++) {
	multiply_to(i, *this);
    }
    return *this;
}

SmallNumber operator* (const SmallNumber& x, const SmallNumber& y)
{
    if (x.logs.empty() || y.logs.empty())
	return SmallNumber();
    SmallNumber z(x);
    for (vector<double>::const_iterator i = y.logs.begin();
	 i < y.logs.end(); i++) {
	multiply_to(i, z);
    }
    return z;
}

const SmallNumber& SmallNumber::operator += (const SmallNumber& y)
{
    for (vector<double>::const_iterator i = y.logs.begin();
	 i < y.logs.end(); i++) {
	add_to(i, *this);
    }
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

XMatrix::XMatrix(int m, int n)
    :entry(m)
{
    for_each(entry.begin(), entry.end(),
	     [n](vector<SmallNumber>& row) {
		 row.resize(n);
	     });
}

XMatrix::XMatrix(const XMatrix& M)
{
    entry.resize(M.entry.size());
    copy(M.entry.begin(), M.entry.end(), entry.begin());
}

SmallNumber& XMatrix::operator() (unsigned i, unsigned j)
{
    return entry[i][j];
}

const SmallNumber& XMatrix::operator() (unsigned i, unsigned j) const
{
    return entry[i][j];
}

XMatrix operator * (const XMatrix& X, const XMatrix& Y)
{
    XMatrix Z(X.entry.size(), Y.entry[0].size());
    for (unsigned i = 0; i < Z.entry.size(); i++) {
	for (unsigned j = 0; j < Z.entry[0].size(); j++) {
	    Z(i, j) = 0;
	    for (unsigned k = 0; k < Z.entry[0].size(); k++) {
		Z(i ,j) += X(i, k) * Y(k, i);
	    }
	}
    }
    return Z;
}

const XMatrix& XMatrix::operator *= (const XMatrix& Y)
{
    XMatrix X(*this);
    for (unsigned i = 0; i < entry.size(); i++) {
	for (unsigned j = 0; j < entry[0].size(); j++) {
	    entry[i][j] = 0;
	    for (unsigned k = 0; k < entry[0].size(); k++) {
		SmallNumber z = X(i, k) * Y(k, i);;
		entry[i][j] += z;
	    }
	}
    }
    return *this;
}
