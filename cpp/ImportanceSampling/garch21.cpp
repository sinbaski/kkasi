#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include "garch21.hpp"

using namespace std;
using namespace arma;

extern size_t num_iter;

double sum_norm(vec& arg)
{
    double s = 0;
    for (size_t i = 0; i < arg.size(); i++) {
	s += arg(i);
    }
    return s;
}

garch21::F_Set::F_Set(garch21 *markov)
{
    this->markov = markov;
}

bool garch21::F_Set::includes(const array<double, 2> &arg) const
{
    double x = arg[0];
    double eta = exp(arg[1]);

    bool ret = x >= x_interval[0] && x <= x_interval[1];
    ret = ret && eta >= eta_interval[0] && eta <= eta_interval[1];

    return ret;
}

garch21::garch21(array<double, 4> &params)
    :F(this), nu(this)
{
    a0 = params[0];
    a1 = params[1];
    a2 = params[2];
    b1 = params[3];
    double t = 1 + (1 + a1) / (1 - a1) / 2;
    double b = t * pow(b1, 2);
    // rho = 1/2 * (1 + b1^2/b)
    double rho = 0.5 * (1 + 1/t);
    F.b = b;
    F.rho = rho;
    F.x_interval[0] = (a1 + rho)/(1 + a1);
    F.x_interval[1] = 1;
    F.eta_interval[0] = pow(b1, 2)/rho;
    F.eta_interval[1] = b;

    C[0] = 0.99;
    C[1] = 1;

    double p = C[0];
    double q = C[1];

    // compute the constants of nu distribution
    {

	double t3 = pow(b1, 0.2e1);
	double t12 = (p * a2 * b1 - p * t3 - a2 * b1 + b) / p / (a1 * b1 + a2);
	double chi1 = gsl_ran_chisq_pdf(t12, 1);

	double t10 = pow(a2, 0.2e1);
	double t17 = (a1 * b1 + a2) * (1 - rho) * b / (a1 + 0.1e1) / (b * rho * a1 + p * a2 * b1 - p * t10 + t10);
	double chi2 = gsl_ran_chisq_pdf(t17, 1);
	nu.c1 = chi1 * chi2;
    }
    assert(nu.c1 > 0);

    {	
	double t7 = q * rho;
	double t8 = pow(a2, 0.2e1);
	double t13 = pow(b1, 0.2e1);
	double t19 = 0.1e1 / (a2 * b1 * t7 + t13 * a1 - t8 * t7 + t8 * rho) / (0.1e1 + a1) / q * rho * (0.1e1 - rho) / 0.2e1;
	nu.c2 = t19;
    }
    assert(nu.c2 > 0);

    {
	nu.c3 = b * b - pow(b1, 4)/rho/rho;
    }
    assert(nu.c3 > 0);

    /* Computes delta, the nu measure of F */
    {
	double t3 = pow(a2, 0.2e1);
	double t7 = b * a1;
	double t9 = log(q * a2 * b1 - t3 * q + t3 + t7);
	double t10 = t9 * t3;
	double t11 = q * rho;
	double t15 = b1 * q * rho;
	double t21 = pow(b1, 0.2e1);
	double t22 = t21 * a1;
	double t24 = 0.1e1 / rho;
	double t26 = log(-t24 * (-a2 * b1 * t11 + t3 * t11 - t3 * rho - t22));
	double t27 = t26 * t3;
	double t40 = pow(a1, 0.2e1);
	double t45 = -t24 / t40 / (0.1e1 + a1) / q * (t15 * t26 * a2 - t15 * t9 * a2 + t11 * t10 - rho * t10 - t11 * t27 + rho * t27 + rho * t7 - t22) * (-0.1e1 + rho) * nu.c1;
	nu.delta = t45;
    }
    assert(nu.delta > 0);
}

double garch21::kernel_density(const array<double, 2> &arg, double x0) const
{
    assert(x0 > 0 && x0 <= 1);

    double x = arg[0];
    double eta = exp(arg[1]);

    double t1 = eta * a1;
    double t3 = a2 * b1;
    double t5 = b1 * b1;
    double t14 = (x * eta + t1 * x + t3 * x0 - t5 * x0 - t1 - t3) / x0 / (a1 * b1 + a2);
    double chi1 = gsl_ran_chisq_pdf(t14, 1);


    double t2 = eta * (1 - x);
    double t7 = a1 * a1;
    double t9 = a2 * a2;
    double t16 = t2 * (a1 * b1 + a2) / (x * eta * a1 + a2 * b1 * x0 - t2 * t7 - t9 * x0 + t9);
    double chi2 = gsl_ran_chisq_pdf(t16, 1);

    t1 = eta * eta;
    double t4 = a1 * a1;
    t5 = eta * t4;
    t9 = a2 * a2;
    double t15 = t1 / x0 / (x * eta * a1 + a2 * b1 * x0 + t5 * x - t9 * x0 - t5 + t9);

    return chi1 * chi2 * t15;
}

inline bool garch21::C_includes(array<double, 2> arg)
{
    return arg[0] >= C[0] && arg[0] <= C[1];
}


array<double, 2> garch21::simulate_path(void)
{
    random_device dev;
    array<double, 2> X = nu.draw(dev);
    // Place the MC initially on the unit sphere.
    X[1] = 0;
    bool reg = false;
    uniform_real_distribution<double> unif;
    size_t n = 0;
    double s = 0;
    do {
	if (C_includes(X)) {
	    double u = unif(dev);
	    reg = u < nu.delta;
	    if (! reg) {
		X = forward(dev, X[0], false);
		n++;
		s += X[1];
	    } else {
		X = nu.draw(dev);
		n++;
		s += X[1];		
	    }
	} else {
	    X = forward(dev, X[0]);
	    n++;
	    s += X[1];
	}
    } while(! reg);
    array<double, 2> ret;
    ret[0] = n;
    ret[1] = s;
    return ret;
}

double tail_index_fun(double alpha, void* param)
{
    vector<double> *stats = (vector<double> *)param;
    double n = (double)stats->size();
    unsigned long counter = 0;
    sort(stats->begin(), stats->end());
    double expected = accumulate(stats->begin(), stats->end(), 0.0,
				 [alpha, &counter](double average, double s)
				 {
				     counter++;
				     if (average == 0.0)
					 return s;
				     double greater, x = s * alpha, y;
				     if (average < x) {
					 greater = x;
					 y = exp(average - greater);
				     } else {
					 greater = average;
					 y = exp(x - greater);
				     }
				     double ret = log(1 + y) + greater;
				     return ret;
				 });
    double ret = expected  - log(n);
    return ret;
}

double garch21::compute_tail_index(size_t beg_line, size_t end_line)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub = 100, xi = -1;

    vector<double> stats;

    char name[64];
    sprintf(name, "%.4f_%.4f_%.4f_stats.txt", a1, a2, b1);

    ifstream infile;
    double x;
    size_t n = 0;
    infile.open(name, ios::in);
    while (infile >> x) {
	if (n >= beg_line && n < end_line)
	    stats.push_back(x);
	n++;
    }
    infile.close();

// #pragma omp parallel for
    for (size_t i = 0; i < num_iter; i++) {
	stats.push_back(simulate_path()[1]);
    }

    ofstream outfile;
    outfile.open(name, ios::out | ios::app);
    for_each(stats.rbegin(), stats.rbegin() + num_iter,
	     [&](double x)
	     {
		 outfile << x << endl;
	     });
    outfile.close();

    double a;
    double bounds[2];
    if (tail_index_fun(a = 2, &stats) > 0) {
	while (tail_index_fun(a = a - 1, &stats) > 0);
	bounds[0] = a;
	bounds[1] = a + 1;
    } else {
	while (tail_index_fun(a = a + 1, &stats) < 0);
	bounds[0] = a - 1;
	bounds[1] = a;
    }

    F.function = tail_index_fun;
    F.params = &stats;

    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);

    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 1.0e-6, 0);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) {
	xi = -1;
    }
    return xi;
}

garch21::nu_dist::nu_dist(garch21 *markov)
{
    this->markov = markov;
}

double garch21::nu_dist::eta_draw(URNG& dev)
{
    bool accepted = false;
    uniform_real_distribution<double> unif;
    double proposed;
    do {
	proposed = eta_proposal_draw(dev);
	double d1 = eta_marginal_density(proposed);
	double d2 = eta_proposal_density(proposed);
	accepted = unif(dev) <= d1 / d2 * delta / c1 / c2 / c3;
    } while (! accepted);
    return proposed;
}

double garch21::nu_dist::eta_marginal_density(double eta)
{
    double a1 = markov->a1;
    double a2 = markov->a2;
    double b1 = markov->b1;
    double q = markov->C[1];
    double rho = markov->F.rho;

    double t4 = pow(a2, 0.2e1);
    double t18 = (double) (1 - rho) / (0.1e1 + a1) / q / (q * a2 * b1 + eta * a1 - t4 * q + t4) / delta * c1 * eta;
    return t18;
}

double garch21::nu_dist::eta_proposal_density(double eta)
{
    return 2 * eta / c3;
}

double garch21::nu_dist::eta_proposal_draw(URNG& dev)
{
    uniform_real_distribution<double> unif(0.0, 1.0);
    double y = unif(dev);

    double b1 = markov->b1;
    double rho = markov->F.rho;

    return sqrt(c3 * y + pow(b1, 4) / pow(rho, 2));
}

/*
  double a1 = markov->a1;
  double a2 = markov->a2;
  double b1 = markov->b1;
 */
double garch21::nu_dist::density(const array<double, 2> &arg) const
{
    if (! markov->F.includes(arg)) return 0;

    double eta = exp(arg[1]);

    double a1 = markov->a1;
    double a2 = markov->a2;
    double b1 = markov->b1;
    double q = markov->C[1];

    double t3 = eta * eta;
    double t4 = pow(a2, 0.2e1);
    double t14 = 0.1e1 / q / (q * a2 * b1 + eta * a1 - t4 * q + t4) * t3 / delta * c1;
    return t14;
}

array<double, 2> garch21::nu_dist::draw(URNG& dev)
{
    array<double, 2> ret;
    uniform_real_distribution<double> unif(0.0, 1.0);
    uniform_real_distribution<double> Ux(markov->F.x_interval[0],
					 markov->F.x_interval[1]);
    ret[0] = Ux(dev);
    ret[1] = log(eta_draw(dev));
    assert(markov->F.includes(ret));
    return ret;
}

array<double, 2> garch21::forward(URNG& dev, double x0, bool orig)
{
    // move forward according to the normal transition kernel if orig = 1
    assert(x0 > 0);
    
    array<mat, 2> A;
    chi_squared_distribution<double> chi2;
    vec V(2);
    array<double, 2> ret;

    V(0) = x0;
    V(1) = 1 - x0;

    for (size_t i = 0; i < A.size(); i++) {
	A[i].set_size(2, 2);
	double z2 = chi2(dev);
	A[i](0, 0) = a1 * z2 + b1;
	A[i](0, 1) = a2;
	A[i](1, 0) = z2;
	A[i](1, 1) = 0;
    }
    V = A[1] * A[0] * V;

    double s = sum_norm(V);
    ret[0] = V(0) / s;
    ret[1] = log(s);

    bool accepted = orig;
    uniform_real_distribution<double> unif;
    while (! accepted) {
	double d1 = kernel_density(ret, x0);
	double d2 = nu.delta * nu.density(ret);

	accepted = unif(dev) < 1 - d2/d1;
	if (! accepted)
	    ret = forward(dev, x0);
    }
    return ret;
}

