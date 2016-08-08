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
    double eta = arg[1];

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
    
    {	    
	double t1 = q * rho;
	double t2 = pow(a2, 0.2e1);
	double t7 = pow(b1, 0.2e1);
	nu.c2 = -0.1e1 / (-t1 * a2 * b1 + t1 * t2 - rho * t2 - a1 * t7) / q * rho;
    }
    
    {
	double t5 = b * b;
	double t7 = pow(b1, 0.2e1);
	double t8 = t7 * t7;
	double t10 = rho * rho;
	nu.c3 = (1 - rho) / (a1 + 1) * (t5 * b - t8 * t7 / t10 / rho) / 0.3e1;
    }

    /* Computes delta, the nu measure of F */
    {
	double t5 = a1 * a1;
	double t6 = 0.1e1 / t5;
	double t7 = a2 * a2;
	double t8 = t6 * t7;
	double t10 = t6 * a2;
	double t13 = 0.1e1 / q;
	double t15 = t13 / a1;
	double t16 = b * b;
	double t19 = t13 * t6;
	double t22 = t7 * t7;
	double t24 = 0.1e1 / t5 / a1;
	double t25 = t22 * t24;
	double t26 = q * t7;
	double t27 = q * a2;
	double t31 = log(b * a1 + t27 * b1 - t26 + t7);
	double t32 = t31 * q;
	double t34 = t7 * a2;
	double t35 = t34 * t24;
	double t39 = t24 * t31;
	double t40 = b1 * b1;
	double t50 = 0.1e1 / rho;
	double t56 = t40 * t40;
	double t57 = rho * rho;
	double t72 = log(-(-t27 * rho * b1 + t26 * rho - a1 * t40 - t7 * rho) * t50);
	double t73 = t22 * t72;
	double t76 = t34 * t72;
	double t91 = t8 * b - t10 * b1 * b + t15 * t16 / 0.2e1 - t19 * t7 * b + t25 * t32 - 0.2e1 * t35 * t32 * b1 + t26 * t39 * t40 - 0.2e1 * t25 * t31 + 0.2e1 * t35 * t31 * b1 + t13 * t22 * t39 - t8 * t40 * t50 + t10 * t40 * b1 * t50 - t15 * t56 / t57 / 0.2e1 + t19 * t7 * t40 * t50 - t73 * q * t24 + 0.2e1 * t76 * q * b1 * t24 - t26 * t72 * t40 * t24 + 0.2e1 * t73 * t24 - 0.2e1 * t76 * b1 * t24 - t73 * t13 * t24;
	double t92 = (0.1e1 - rho) / (a1 + 0.1e1) * t91;
	nu.delta = nu.c1 * t92;
    }
}

double garch21::kernel_density(const array<double, 2> &arg, double x0) const
{
    assert(x0 > 0 && x0 < 1);

    double x = arg[0];
    double eta = arg[1];

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
    X[1] = 1;
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
		s += log(X[1]);
	    }
	} else {
	    X = forward(dev, X[0]);
	    n++;
	    s += log(X[1]);
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
    double expected = 0;
    double greatest = *max_element(stats->begin(), stats->end());
    double n = (double)stats->size();
    expected = accumulate(stats->begin(), stats->end(), 0.0,
			  [=](double average, double xi)
			  {
			      return average + exp(alpha * (xi - greatest));
			  });
    double ret = log(expected) + alpha * greatest - log(n);
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
    sprintf(name, "%3.2f_%3.2f_%3.2f_stats.txt", a1, a2, b1);

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

#pragma omp parallel for
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
    for (a = 2; a > 0 && tail_index_fun(a, &stats) > 0; a -= 1);
    if (a > 0) {
    	bounds[0] = a;
    } else {
    	printf("%s %.2f.\n", "lower bound less than ", a);
    	return -1;
    }
    for (a = 2; a < ub && tail_index_fun(a, &stats) < 0; a += 1);
    if (a < ub) {
    	bounds[1] = a;
    } else {
    	printf("%s %.2f.\n", "Upper bound larger than ", (double)ub);
    	return -1;
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

double garch21::nu_dist::proposal_density(array<double, 2> arg)
{
    return arg[1] * arg[1] / c3;
}

array<double, 2> garch21::nu_dist::proposal_draw(URNG& dev)
{
    array<double, 2> ret;
    uniform_real_distribution<double> unif(0.0, 1.0);
    double y = unif(dev);
    
    double b1 = markov->b1;
    double a1 = markov->a1;
    double rho = markov->F.rho;
    
    double t1 = pow(b1, 0.2e1);
    double t2 = t1 * t1;
    double t4 = rho * rho;
    double t10 = pow((3 * c3 * y * t4 * rho) + t2 * t1, 0.1e1 / 0.3e1);
    ret[1] = 0.1e1 / rho * t10;

    uniform_real_distribution<double> U((a1 + rho)/(1 + a1), 1.0);
    ret[0] = U(dev);

    return ret;
}

/*
  double a1 = markov->a1;
  double a2 = markov->a2;
  double b1 = markov->b1;
 */
double garch21::nu_dist::density(const array<double, 2> &arg) const
{
    if (! markov->F.includes(arg)) return 0;

    double eta = arg[1];

    double a1 = markov->a1;
    double a2 = markov->a2;
    double b1 = markov->b1;
    double q = markov->C[1];

    double t1 = eta * eta;
    double t2 = pow(a2, 0.2e1);
    double t11 = t1 / (q * a2 * b1 + eta * a1 - q * t2 + t2) / q;

    return c1 * t11 / delta;
}

array<double, 2> garch21::nu_dist::draw(URNG& dev)
{
    array<double, 2> proposed;
    uniform_real_distribution<double> unif(0.0, 1.0);
    double d1, d2, u;
    do {
	proposed = proposal_draw(dev);
	d2 = proposal_density(proposed);
	d1 = density(proposed);
	u = unif(dev);
    } while (u > d1 / d2 / c1 / c2 / c3 * delta);
    assert(markov->F.includes(proposed));
    return proposed;
}

array<double, 2> garch21::forward(URNG& dev, double x0, bool orig)
{
    // move forward according to the normal transition kernel if orig = 1
    array<mat, 2> A;
    chi_squared_distribution<double> chi2;
    vec V(2);
    array<double, 2> ret;

    V(0) = x0;
    V(1) = 1 - x0;

//    V.print();
    for (size_t i = 0; i < A.size(); i++) {
	A[i].set_size(2, 2);
	double z2 = chi2(dev);
	A[i](0, 0) = a1 * z2 + b1;
	A[i](0, 1) = a2;
	A[i](1, 0) = z2;
	A[i](1, 1) = 0;
    }
    // A[0].print();
    // A[1].print();
    
    V = A[1] * A[0] * V;
//    V.print();
    
    ret[1] = sum_norm(V);
    ret[0] = V(0) / ret[1];

    bool accepted = orig;
    while (! accepted) {
	double d1 = kernel_density(ret, x0);
	double d2 = nu.delta * nu.density(ret);
	uniform_real_distribution<double> unif;
	accepted = unif(dev) < 1 - d2/d1;
	if (! accepted)
	    ret = forward(dev, x0);
    }
    return ret;
    
}

