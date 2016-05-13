#include <armadillo>

#if defined DEBUG
#define container_t vector
#include <vector>
#else
#define container_t list
#include <list>
#endif


using namespace std;
using namespace arma;

class ExtremeNumber
{
protected:
    const ExtremeNumber& erase_small(void);

public:
    container_t<double> logs;

    ExtremeNumber(void);
    ExtremeNumber(double x);
    ExtremeNumber(const container_t<double>& logs);
    ExtremeNumber(const ExtremeNumber& x);
    ~ExtremeNumber(void) {logs.clear();}

    const ExtremeNumber& operator= (double y);
    const ExtremeNumber& operator= (const ExtremeNumber &y);

    const ExtremeNumber& operator *= (double y);
    friend ExtremeNumber operator * (const ExtremeNumber& x, double y);

    const ExtremeNumber& operator /= (double y);
    friend ExtremeNumber operator / (const ExtremeNumber& x, double y);

    const ExtremeNumber& operator *= (const ExtremeNumber& y);
    friend ExtremeNumber operator * (const ExtremeNumber& x, const ExtremeNumber& y);
    
    friend ExtremeNumber operator + (const ExtremeNumber& x, double y);
    const ExtremeNumber& operator += (double y);
    
    friend ExtremeNumber operator + (const ExtremeNumber& x, const ExtremeNumber& y);
    const ExtremeNumber& operator += (const ExtremeNumber& y);

    ExtremeNumber operator ^ (unsigned long n) const;

    double nth_root (unsigned long n) const;

    bool operator < (const ExtremeNumber& x) const;

    long leading_power(void) const;
    double comptify(long power) const;
    double comptify(void) const;

    friend double log10(const ExtremeNumber& x);
    friend double log(const ExtremeNumber& x);
};



class XMatrix
{
public:
   
    vector< vector<ExtremeNumber> > entry;

public:

    XMatrix(void);
    XMatrix(const XMatrix& M);
    XMatrix(unsigned int m, unsigned int n);
    XMatrix(double** data, unsigned int m, unsigned int n);
    ExtremeNumber& operator() (unsigned i, unsigned j);
    const ExtremeNumber& operator() (unsigned i, unsigned j) const;
    const XMatrix& operator *= (const XMatrix& Y);
    friend XMatrix operator * (const XMatrix& X, const XMatrix& Y);
    const XMatrix& operator += (const XMatrix& Y);
    friend XMatrix operator + (const XMatrix& X, const XMatrix& Y);
    const XMatrix& operator = (const XMatrix& Y);
    const XMatrix& operator ^ (unsigned n);
    XMatrix operator ~ (void);
    mat comptify(long *power);
};
