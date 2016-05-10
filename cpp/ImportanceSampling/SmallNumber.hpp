#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class SmallNumber
{
// private:
//     friend void add_to (vector<double>::const_iterator i1,
// 			SmallNumber& V);
    
//     friend void multiply_to (vector<double>::const_iterator i,
// 			     const SmallNumber& V);

public:
    vector<double> logs;

    SmallNumber(void);
    SmallNumber(double x);
    SmallNumber(const vector<double>& logs);
    SmallNumber(const SmallNumber& x);

    const SmallNumber& operator= (double y);
    const SmallNumber& operator= (const SmallNumber &y);

    const SmallNumber& erase_small(void);

    const SmallNumber& operator *= (double y);
    friend SmallNumber operator * (const SmallNumber& x, double y);

    const SmallNumber& operator *= (const SmallNumber& y);
    friend SmallNumber operator * (const SmallNumber& x, const SmallNumber& y);
    
    friend SmallNumber operator + (const SmallNumber& x, double y);
    const SmallNumber& operator += (double y);
    
    friend SmallNumber operator + (const SmallNumber& x, const SmallNumber& y);
    const SmallNumber& operator += (const SmallNumber& y);
};

class XMatrix
{
public:
   
    vector< vector<SmallNumber> > entry;

public:

    XMatrix(void);
    XMatrix(const XMatrix& M);
    XMatrix(unsigned int m, unsigned int n);
    XMatrix(double** data, unsigned int m, unsigned int n);
    SmallNumber& operator() (unsigned i, unsigned j);
    const SmallNumber& operator() (unsigned i, unsigned j) const;
    const XMatrix& operator *= (const XMatrix& Y);
    friend XMatrix operator * (const XMatrix& X, const XMatrix& Y);
    const XMatrix& operator += (const XMatrix& Y);
    friend XMatrix operator + (const XMatrix& X, const XMatrix& Y);
    const XMatrix& operator = (const XMatrix& Y);
    const XMatrix& operator ^ (unsigned n);
    XMatrix operator ~ (void);
};
