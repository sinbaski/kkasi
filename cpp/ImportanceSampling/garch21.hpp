#ifndef XXIE_GARCH21_H
#define XXIE_GARCH21_H

#include "GarchPQ.hpp"

typedef array<double, 2> funval;
// typedef vector<funval> emp_fun;

double interpolate_fun
(double arg, const vector<funval> &fun);

class Garch21 : public GarchPQ
{
public:
    Garch21(const vector<double> &alpha, const vector<double> &beta,
	    double tail_index_sup = 5.0);
    const double tail_index;
    double quantile(double u, double angle) const;
    double draw_z2(void) const;
    double draw_z2(double) const;
    void simulate_sample_path(double threshold, unsigned int n,
			      vector<vec> &path) const;

    double Lyapunov(size_t n, size_t iterations);

    double right_eigenfunction
    (double index, vector<funval> &eigenfunction) const;

    double M_fun(double theta);
    double b_fun(double theta);
    
protected:
    vector<double> pool;
    // vector<funval> eigenfunction;

};



#endif
