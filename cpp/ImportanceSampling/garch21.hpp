#ifndef XXIE_GARCH21_H
#define XXIE_GARCH21_H

#include "GarchPQ.hpp"

class Garch21 : public GarchPQ
{
public:
    Garch21(const vector<double> &alpha, const vector<double> &beta,
	    double tail_index_sup = 5.0);
    const double tail_index;
    double r(double) const;
    double quantile(double u, double angle) const;
    double draw_z2(void) const;
    double draw_z2(double) const;
    void simulate_sample_path(double threshold, unsigned int n,
			      vector<vec> &path) const;

protected:
    vector<double> pool;
    typedef array<double, 2> funval;
    vector<funval> eigenfunction;
    void right_eigenfunction(unsigned n);
};

#endif
