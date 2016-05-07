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
	logs.push_back(log(x));
}

SmallNumber::SmallNumber(const vector<double>& logs)
    : logs(logs) {
}

SmallNumber::SmallNumber(const SmallNumber& x)
    : logs(x.logs) {
}

const SmallNumber& SmallNumber::operator= (double y)
{
    assert(y >= 0);
    if (y == 0)
	logs.clear();
    else 
	logs.push_back(log(y));
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
	
    double y1 = log(y);
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
	    V[i * x.logs.size() + j] = x.logs[i] + y.logs[j];
	}
    }
    return SmallNumber(V);
}

const SmallNumber& SmallNumber::operator += (double y)
{
    assert(y >= 0);
    
    if (y == 0) return *this;
    logs.insert(logs.end(), log(y));
    return *this;
}
    
SmallNumber operator + (const SmallNumber& x, double y)
{
    SmallNumber z(x);
    z += y;
    return z;
}

const SmallNumber& SmallNumber::operator += (const SmallNumber& y)
{
    logs.insert(logs.end(), y.logs.begin(), y.logs.end());
    // for_each(y.logs.begin(), y.logs.end(),
    // 	     [&](double element){
    // 		 logs.push_back(element);
    // 	     });
    return *this;
}
    
SmallNumber operator + (const SmallNumber& x, const SmallNumber& y)
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
		entry[i][j] += X(i, k) * Y(k, i);
	    }
	}
    }
    return *this;
}
