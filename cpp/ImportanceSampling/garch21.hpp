#ifndef XXIE_GARCH21_H
#define XXIE_GARCH21_H

#include "GarchPQ.hpp"
typedef array<double, 2> funval;
// typedef vector<funval> emp_fun;

double interpolate_fun
(double arg, const vector<funval> &fun);

class Garch21 : public GarchPQ
{
protected:
    const size_t nbr_eigenfunctions;
    const size_t nbr_eigenfunction_points;

    vector<double> pool;
    vector<mat> A_matrices;

    struct eigenfunction {
	double kappa;
	double lambda_kappa;
	vector<funval> *r_kappa;
    };

    vector<eigenfunction> eigenfunctions;
    void compute_eigenfunctions(void);
    double compute_M(void);

public:
    Garch21(const vector<double> &alpha, const vector<double> &beta,
	    double tail_index_sup = 5.0);
    ~Garch21(void);

    // const double tail_index;
    double tail_index;
    const double M;
    vector<funval> r_xi;
    
    double quantile(double u, double angle) const;
    double draw_z2(void) const;
    double draw_z2(double) const;
    void simulate_path(vector<vec> &path) const;
    pair<double, size_t> sample_estimator(const vec &V0, double u);
    vector<double> estimate_prob(double u, size_t);

    double Lyapunov(size_t n, size_t iterations);

    double right_eigenfunction
    (double index, vector<funval> &eigenfunction) const;

};



#endif
