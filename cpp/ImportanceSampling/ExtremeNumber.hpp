#include <list>
#include <armadillo>

using namespace std;
using namespace arma;

class ExtremeNumber
{
// private:
//     friend void add_to (list<double>::const_iterator i1,
// 			ExtremeNumber& V);
    
//     friend void multiply_to (list<double>::const_iterator i,
// 			     const ExtremeNumber& V);

public:
    list<double> logs;

    ExtremeNumber(void);
    ExtremeNumber(double x);
    ExtremeNumber(const list<double>& logs);
    ExtremeNumber(const ExtremeNumber& x);
    ~ExtremeNumber(void) {logs.clear();}

    const ExtremeNumber& operator= (double y);
    const ExtremeNumber& operator= (const ExtremeNumber &y);

    const ExtremeNumber& erase_small(void);

    const ExtremeNumber& operator *= (double y);
    friend ExtremeNumber operator * (const ExtremeNumber& x, double y);

    const ExtremeNumber& operator *= (const ExtremeNumber& y);
    friend ExtremeNumber operator * (const ExtremeNumber& x, const ExtremeNumber& y);
    
    friend ExtremeNumber operator + (const ExtremeNumber& x, double y);
    const ExtremeNumber& operator += (double y);
    
    friend ExtremeNumber operator + (const ExtremeNumber& x, const ExtremeNumber& y);
    const ExtremeNumber& operator += (const ExtremeNumber& y);

    double comptify(int power);
    int leading_power(void);
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
    mat comptify(int *power);
};
