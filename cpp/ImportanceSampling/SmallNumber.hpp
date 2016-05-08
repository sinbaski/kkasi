#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

class SmallNumber
{
private:
    friend void add_to (vector<double>::const_iterator i1,
			vector<double>::iterator o1,
			vector<double>::iterator o2,
			vector<double>& V);

    friend void add_to (vector<double>::const_iterator i1,
			vector<double>::const_iterator i1,
			vector<double>::iterator o1,
			vector<double>::iterator o2,
			vector<double>& V);

    friend void multiply_to (vector<double>::const_iterator i1,
			     vector<double>::iterator o1,
			     vector<double>::iterator o2,
			     vector<double>& V);

    friend void multiply_to (vector<double>::const_iterator i1,
			     vector<double>::const_iterator i1,
			     vector<double>::iterator o1,
			     vector<double>::iterator o2,
			     vector<double>& V);

public:
    vector<double> logs;

    SmallNumber(void);
    SmallNumber(double x);
    SmallNumber(const vector<double>& logs);
    SmallNumber(const SmallNumber& x);

    const SmallNumber& operator= (double y);
    const SmallNumber& operator= (const SmallNumber &y);

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
    XMatrix(int m, int n);
    SmallNumber& operator() (unsigned i, unsigned j);
    const SmallNumber& operator() (unsigned i, unsigned j) const;
    const XMatrix& operator *= (const XMatrix& Y);
    friend XMatrix operator * (const XMatrix& X, const XMatrix& Y);
};
