#include <cmath>
#include <cassert>
#include <algorithm>
#include "ExtremeNumber2.hpp"


using namespace std;

ExtremeNumber::ExtremeNumber(void)
    : sign(0)
{
}

ExtremeNumber::ExtremeNumber(double x)
{
    operator = (x);
}

ExtremeNumber::ExtremeNumber(const ExtremeNumber& x)
    :sign(x.sign), mylog(x.mylog)
{
}

const ExtremeNumber& ExtremeNumber::operator = (double x)
{
    if (x == 0) {
	sign = 0;
    } else if (x > 0) {
	sign = 1;
	mylog = log10(x);
    } else {
	sign = -1;
	mylog = log10(-x);
    }
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator = (const ExtremeNumber &y)
{
    sign = y.sign;
    mylog = y.mylog;
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator *= (double y)
{
    return operator *= (ExtremeNumber(y));
}

const ExtremeNumber& ExtremeNumber::operator *= (const ExtremeNumber& y)
{
    if (sign == 0) return *this;
    if (y.sign == 0) {
	sign = 0;
	return *this;
    }
    if (sign == y.sign) {
	sign = 1;
    } else {
	sign = -1;
    }
    mylog += y.mylog;
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator /= (double y)
{
    assert(y != 0);
    if (sign == 0) return *this;
    ExtremeNumber x(1/y);
    return operator *= (x);
}

const ExtremeNumber& ExtremeNumber::operator /= (const ExtremeNumber& y)
{
    assert(y.sign != 0);
    if (sign == 0) return *this;
    if (sign != y.sign)
	sign = -1;
    else
	sign = 1;
    mylog -= y.mylog;
    return *this;
}

const ExtremeNumber& ExtremeNumber::operator += (double y)
{
    return operator += (ExtremeNumber(y));
}

const ExtremeNumber& ExtremeNumber::operator += (const ExtremeNumber& y)
{
    if (y.sign == 0) return *this;
    if (sign == 0) {
	operator = (y);
	return *this;
    }
    if (sign == y.sign) {
	double f = max(floor(mylog), floor(y.mylog));
	mylog = log10(pow(10, mylog - f) + pow(10, y.mylog - f)) + f;
	return *this;
    }
    if (mylog == y.mylog) {
	sign = 0;
	return *this;
    }
    return operator -= (-y);
}

const ExtremeNumber& ExtremeNumber::operator -= (double x)
{
    return operator -= (ExtremeNumber(x));
}

const ExtremeNumber& ExtremeNumber::operator -= (const ExtremeNumber &x)
{
    if (x.sign == 0) return *this;
    if (sign == 0) {
	operator = (x);
	sign = -sign;
	return *this;
    }
    if (sign == x.sign) {
	if (mylog == x.mylog) {
	    sign = 0;
	    return *this;
	}
	if (mylog > x.mylog) {
	    double f = floor(x.mylog);
	    mylog = log10(pow(10, mylog - f) - pow(10, x.mylog - f)) + f;
	    return *this;
	} else {
	    sign = -sign;
	    double f = floor(mylog);
	    mylog = log10(-pow(10, mylog - f) + pow(10, x.mylog - f)) + f;
	    return *this;
	}
    }
    if (sign == 1)
	return operator += (-x);
    else {
	sign = 1;
	operator += (x);
	sign = -1;
	return *this;
    }
}

const ExtremeNumber& ExtremeNumber::operator ^= (double u)
{
    assert((sign == 1) || (sign == 0 && u >= 0));
    if (u == 0) {
	sign = 1;
	mylog = 0;
	return *this;
    }
    if (sign == 0) return *this;
    mylog *= u;
    return *this;
}

bool ExtremeNumber::operator < (const ExtremeNumber& x) const
{
    if (sign < x.sign) return true;
    if (sign > x.sign) return false;
    return mylog < x.mylog;
}

double ExtremeNumber::operator () (void) const
{
    if (sign == 0) return 0;
    double x = pow(10, mylog - floor(mylog));
    return sign == 1 ? x : -x;
}

double ExtremeNumber::operator () (long n) const
{
    if (sign == 0) return 0;
    double x = pow(10, mylog - n);
    return sign == 1 ? x : -x;
}

ExtremeNumber::operator double() const
{
    if (mylog == 0) return 0;
    return pow(10, mylog);
}

long ExtremeNumber::power(void) const
{
    return (long)floor(mylog);
}

ExtremeNumber operator * (const ExtremeNumber& x, double y)
{
    ExtremeNumber z(x);
    z *= y;
    return z;
}

ExtremeNumber operator * (const ExtremeNumber& x, const ExtremeNumber& y)
{
    ExtremeNumber z(x);
    z *= y;
    return z;
}

ExtremeNumber operator / (const ExtremeNumber& x, double y)
{
    ExtremeNumber z(x);
    z /= y;
    return z;
}

ExtremeNumber operator / (const ExtremeNumber& x, const ExtremeNumber& y)
{
    ExtremeNumber z(x);
    z /= y;
    return z;
}

ExtremeNumber operator + (const ExtremeNumber &x, const ExtremeNumber& y)
{
    ExtremeNumber z(x);
    z += y;
    return z;
}
    
ExtremeNumber operator + (const ExtremeNumber &x, double y)
{
    ExtremeNumber z(x);
    z += y;
    return z;
}

ExtremeNumber operator - (const ExtremeNumber &x)
{
    ExtremeNumber y(x);
    if (y.sign != 0) y.sign = -y.sign;
    return y;
}

ExtremeNumber operator - (const ExtremeNumber &x, const ExtremeNumber &y)
{
    ExtremeNumber z(x);
    z -= y;
    return z;
}

ExtremeNumber operator - (const ExtremeNumber &x, double y)
{
    ExtremeNumber z(x);
    z -= y;
    return z;
}

ExtremeNumber operator ^ (const ExtremeNumber &x, double u)
{
    ExtremeNumber z(x);
    return z ^= u;
}

ExtremeNumber abs(const ExtremeNumber &x)
{
    ExtremeNumber y(x);
    y.sign = 1;
    return y;
}

double log10(const ExtremeNumber &x)
{
    assert(x.sign == 1);
    return x.mylog;
}

double log(const ExtremeNumber &x)
{
    return log10(x) * log(10);
}

// ExtremeNumber log10(const ExtremeNumber &x)
// {
//     assert(x.sign == 1);
//     ExtremeNumber y = x;
//     y.mylog = log10(x.mylog);
//     return y;
// }

ExtremeNumber exp(const ExtremeNumber &x)
{
    ExtremeNumber R;
    long p = x.power();
    long inc = p > 0 ? 1 : -1;
    double q = p > 0 ? 10 : 0.1;
    double y = 1;
    for (long i = 0; i != p; i += inc) {
	y *= q;
    }
    R.mylog = x() * y / log(10);
    return R;
}
