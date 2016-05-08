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
	   vector<double> &V)
{
    double f = floor(*o1);
    double r = log10(pow(10, *o1 - f) + 1);
    if (r < 1) {
	*o1 = f + r;
	return false;
    } else {
	*o1 = f + r - 1;
	vector<double>::iterator i = 
	    upper_bound(o1 + 1, o2, f + 1,
			[](double x, double y) {
			    return floor(x) < y;
			});
	if (floor(i - 1) < f + 1) {
	    V.insert(i - 1, f + 1);
	} else {
	    carry(i - 1, o2, V);
	}
    }
}

void SmallNumber::add_to (vector<double>::const_iterator i1,
			  vector<double>::iterator o1,
			  vector<double>::iterator o2,
			  vector<double>& V)
{
    double f = floor(*i1);
    vector<double>::iterator k =
	lower_bound(o1, o2, f,
		    [](double x, double y){
			return floor(x) < y;
		    });
    if (floor(*k) < f) {
	V.insert(k + 1, *i1);
    } else {
	double r = log10(pow(10, *k - f) + pow(10, *i1 - f));
	if (r >= 1) {
	    *k = f + r - 1;
	    if ((int)(floor(*(k + 1)) - f) == 1)
		carry(k + 1, o2, V);
	    else
		V.insert(k + 1, f + 1);
	} else {
	    *k = f + r;
	}
    }
}

void add_to (vector<double>::const_iterator i1,
	     vector<double>::const_iterator i2,
	     vector<double>::iterator o1,
	     vector<double>::iterator o2,
	     vector<double>& V)
{
    while (i1 < i2) add_to(i1++, o1, o2, V);
}

void adjust_numbers(vector<double>::iterator o1,
		    vector<double>& V)
{
    vector<double>::iterator k;
    if (o1 >= V.end() - 1) return;
    while (o1 < V.end()) {
	double f = floor(*o1);
	k = upper_bound(o1 + 1, V.end(), f
			[](double x, double y) {
			    return floor(x) < y;
			});
	double x = pow(10, *o1 - f);
	double r = accumulate(o1+1, k, x,
			      [](double y) {
				  return x + pow(10, y - f);
			      });
	r = log10(r);
	if (r < 1) {
	    *o1 = f + r;
	    V.erase(o1 + 1, k);
	} else { // r >= 1
	    double d = floor(r);
	    *o1 = f + r - d;
	    vector<double>::const_iterator i;
	    i = upper_bound(k, V.end(), f + d,
			    [](double x, double y) {
				return floor(x) < y;
			    });
	    if (i == k)
	    i--;
	    if (floor(*i) < f + d) {
		V.insert(i + 1, f + d);
	    } else {
		carry(i, V.end(), V);
	    }
	}
    }
}

void multiply_to (vector<double>::const_iterator i1,
		  vector<double>::iterator o1,
		  vector<double>::iterator o2,
		  vector<double>& V)
{
    while (o1 < o2) {
	*o1++ += *i1;
    }
    adjust_numbers(o1, o2, V);
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
	
    double y1 = log10(y);
    for_each(logs.begin(), logs.end(),
	     [y1](double &x) {
    	    x += y1;
    	});
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
    vector<double> V(logs.size() * y.logs.size());
    
#pragma omp parallel for
    for (unsigned i = 0; i < logs.size(); i++) {
	for (unsigned j = 0; j < y.logs.size(); j++) {
	    V[i * logs.size() + j] = logs[i] + y.logs[j];
	}
    }
    return *this;
}

SmallNumber operator* (const SmallNumber& x, const SmallNumber& y)
{
    if (x.logs.empty() || y.logs.empty())
	return SmallNumber();

    vector<double> V(x.logs.size() * y.logs.size());
#pragma omp parallel for
    for (unsigned i = 0; i < x.logs.size(); i++) {
	for (unsigned j = 0; j < y.logs.size(); j++) {
	    V[i * y.logs.size() + j] = x.logs[i] + y.logs[j];
	}
    }
    return SmallNumber(V);
}

const SmallNumber& SmallNumber::operator += (const SmallNumber& y)
{
    add_to(y.logs.begin(), y.logs.end(),
	   logs.begin(), logs.end(), logs);
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
    SmallNumber Y(y)
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
