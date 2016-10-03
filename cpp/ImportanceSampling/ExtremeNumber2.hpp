#ifndef ExtremeNumber_hpp
#define ExtremeNumber_hpp

class ExtremeNumber
{
public:
    int sign;
    double mylog;

    ExtremeNumber(void);
    ExtremeNumber(double x);
    ExtremeNumber(const ExtremeNumber& x);

    const ExtremeNumber& operator= (double y);
    const ExtremeNumber& operator= (const ExtremeNumber &y);

    const ExtremeNumber& operator *= (double y);
    const ExtremeNumber& operator *= (const ExtremeNumber& y);

    const ExtremeNumber& operator /= (double y);
    const ExtremeNumber& operator /= (const ExtremeNumber& y);
    
    const ExtremeNumber& operator += (double y);
    const ExtremeNumber& operator += (const ExtremeNumber& y);

    const ExtremeNumber& operator -= (double x);
    const ExtremeNumber& operator -= (const ExtremeNumber &x);

    const ExtremeNumber& operator ^= (double u);
    
    bool operator < (const ExtremeNumber& x) const;
    
    double operator() (void) const;
    double operator() (long n) const;
    long power(void) const;
};

ExtremeNumber operator * (const ExtremeNumber& x, double y);
ExtremeNumber operator * (const ExtremeNumber& x, const ExtremeNumber& y);
ExtremeNumber operator / (const ExtremeNumber& x, double y);
ExtremeNumber operator / (const ExtremeNumber& x, const ExtremeNumber& y);
ExtremeNumber operator + (const ExtremeNumber& x, double y);
ExtremeNumber operator + (const ExtremeNumber& x, const ExtremeNumber& y);
ExtremeNumber operator - (const ExtremeNumber &x);
ExtremeNumber operator - (const ExtremeNumber &x, const ExtremeNumber &y);
ExtremeNumber operator - (const ExtremeNumber &x, double y);
ExtremeNumber operator ^ (const ExtremeNumber &x, double u);

ExtremeNumber abs(const ExtremeNumber &x);
double log10(const ExtremeNumber &x);
double log(const ExtremeNumber &x);
double exp(const ExtremeNumber &x);


#endif
