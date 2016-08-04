#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
#include <armadillo>
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

    C[0] = 0.5;
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
	    
	double t2 = pow(a2, 0.2e1);
	double t6 = pow(b1, 0.2e1);
	nu.c2 = 0.1e1 / q / (q * a2 * b1 - q * t2 + a1 * t6 + t2);

	double t5 = b * b;
	double t7 = pow(b1, 0.2e1);
	double t8 = t7 * t7;
	t10 = rho * rho;
	nu.c3 = (1 - rho) / (a1 + 1) * (t5 * b - t8 * t7 / t10 / rho) / 0.3e1;
    }

    /* Computes delta, the nu measure of F */
    {
	double t1 = b1 * b1;
	double t2 = a2 * a2;
	double t3 = t1 * t2;
	double t4 = q * a1;
	double t7 = t1 * b1 * a2;
	double t9 = q * b1;
	double t10 = t9 * a2;
	double t11 = q * t2;
	double t12 = t1 * a1;
	double t14 = log(t10 - t11 + t12 + t2);
	double t15 = q * t14;
	double t16 = t2 * a2;
	double t17 = t16 * b1;
	double t20 = q * q;
	double t21 = t20 * t14;
	double t28 = log(b * rho * a1 + t10 - t11 + t2);
	double t29 = q * t28;
	double t32 = a1 * rho;
	double t34 = b * t2;
	double t36 = rho * rho;
	double t41 = log(b * a1 + t10 - t11 + t2);
	double t42 = t20 * t41;
	double t43 = t2 * t2;
	double t44 = t43 * t36;
	double t54 = log(-(q * rho * t2 - t9 * rho * a2 - rho * t2 - t12) / rho);
	double t55 = t20 * t54;
	double t57 = t41 * t43;
	double t58 = q * t36;
	double t61 = t54 * t43;
	double t64 = b * b;
	double t66 = a1 * a1;
	double t67 = t66 * t36;
	double t69 = t20 * t28;
	double t79 = -t34 * a1 * t36 + t28 * t64 * t67 - t41 * t64 * t67 + 0.2e1 * t15 * t17 - 0.2e1 * t21 * t17 - 0.2e1 * t29 * t17 + 0.2e1 * t69 * t17 + t21 * t3 + t21 * t43 + 0.2e1 * t29 * t43 + t3 * t32 + t3 * t4 - t69 * t3 + t34 * t32 + t57 * t36 - t7 * t4 + t42 * t44 - t55 * t44 - 0.2e1 * t57 * t58 + 0.2e1 * t61 * t58;
	double t82 = t1 * t1;
	double t91 = b * a2 * b1;
	double t92 = t4 * t36;
	double t94 = t4 * rho;
	double t98 = t17 * t36;
	double t101 = t3 * t36;
	double t107 = t9 * t36;
	double t117 = 0.2e1 * t41 * t16 * t107 - 0.2e1 * t54 * t16 * t107 - t14 * t82 * t66 + t54 * t82 * t66 - t3 * a1 + t42 * t101 - t55 * t101 + t14 * t43 - 0.2e1 * t15 * t43 - t28 * t43 - t3 * t94 + t34 * t92 - t34 * t94 - t61 * t36 - 0.2e1 * t42 * t98 - t69 * t43 + 0.2e1 * t55 * t98 + t7 * t94 - t91 * t92 + t91 * t94;
	nu.delta = -(t79 + t117) / q / t66 / a1 / (a1 + 0.1e1) / t36 / 0.2e1;
    }
}

double garch21::kernel_density(array<double, 2> arg, double x0)
{
    assert(x0 > 0 && x0 < 1);
    assert(state_space_includes(arg));

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

bool garch21::state_space_includes(array<double, 2> arg)
{
    double t2 = pow(b1, 0.2e1);
    double t8 = (arg[1] * a1 + t2) / (a1 + 0.1e1) / arg[1];

    return arg[1] > t2 && arg[0] > t8 && arg[0] < 1;
}

inline bool garch21::C_includes(array<double, 2> arg)
{
    return arg[0] >= C[0] && arg[0] <= C[1];
}


array<double, 2> garch21::simulate_path(void)
{
    array<double, 2> X = nu.draw();
    // Place the MC initially on the unit sphere.
    X[1] = 1;
    bool reg = false;
    uniform_real_distribution<double> unif;
    size_t n = 0;
    double es = X[1];
    do {
	if (C_includes(X)) {
	    double u = unif(dev);
	    reg = u < nu.delta;
	    if (! reg) {
		X = forward(X[0], false);
		n++;
		es *= X[1];
	    }
	} else {
	    X = forward(X[0]);
	    n++;
	    es *= X[1];
	}
    } while(! reg);
    array<double, 2> ret;
    ret[0] = n;
    ret[1] = es;
    return ret;
}

double tail_index_fun(double alpha, void* param)
{
    static vector<double> stats(num_iter);
    static int counter = 0;
    garch21 *markov = (garch21 *)param;
    if (counter == 0) {
//#pragma omp parallel for
	for (size_t i = 0; i < stats.size(); i++) {
	    stats[i] = markov->simulate_path()[1];
	}
	counter++;
    }
    return accumulate(stats.begin(), stats.end(), 0.0,
    			      [alpha](double average, double eta)
    			      {
    				  return average + pow(eta, alpha);
    			      });
}

double garch21::compute_tail_index(void)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub = 30, xi = -1;
    
    double a;
    double bounds[2];
    for (a = 2; a > 0 && tail_index_fun(a, this) > 0; a -= 1);
    if (a > 0) {
    	bounds[0] = a;
    } else {
    	printf("%s %.2f.\n", "lower bound less than ", a);
    	return -1;
    }
    for (a = 2; a < ub && tail_index_fun(a, this) < 0; a += 1);
    if (a < ub) {
    	bounds[1] = a;
    } else {
    	printf("%s %.2f.\n", "Upper bound larger than ", (double)ub);
    	return -1;
    }

    F.function = tail_index_fun;
    F.params = this;

    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);

    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 1.0e-6, 0);
	// if (status == GSL_SUCCESS)
	//     cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) {
	// cout << "The Brent algorithm did not converge after " << max_iter
	//      << " iterations." << endl;
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
    return pow(arg[1], 2) / c3;
}

array<double, 2> garch21::nu_dist::proposal_draw(void)
{
    array<double, 2> ret;
    uniform_real_distribution<double> unif(0.0, 1.0);
    double y = unif(markov->dev);
    
    double b1 = markov->b1;
    double a1 = markov->a1;
    double rho = markov->F.rho;
    
    double t1 = pow(b1, 0.2e1);
    double t2 = t1 * t1;
    double t4 = rho * rho;
    double t10 = pow((3 * c3 * y * t4 * rho) + t2 * t1, 0.1e1 / 0.3e1);
    ret[1] = 0.1e1 / rho * t10;

    uniform_real_distribution<double> U((a1 + rho)/(1 + a1), 1.0);
    ret[0] = U(markov->dev);

    return ret;
}

/*
  double a1 = markov->a1;
  double a2 = markov->a2;
  double b1 = markov->b1;
 */
double garch21::nu_dist::density(array<double, 2> arg)
{
    double x = arg[0];
    double eta = arg[1];

    if (x < markov->F.x_interval[0] || x > markov->F.x_interval[1])
	return 0.0;
    if (eta < markov->F.eta_interval[0] || eta > markov->F.eta_interval[1])
	return 0.0;

    double a1 = markov->a1;
    double a2 = markov->a2;
    double b1 = markov->b1;
    double q = markov->C[1];

    double t1 = eta * eta;
    double t4 = eta * x;
    double t5 = pow(a1, 0.2e1);
    double t9 = pow(a2, 0.2e1);
    double t15 = t1 / q / (q * a2 * b1 - eta * t5 - q * t9 + t4 * t5 + t4 * a1 + t9);

    return c1 * t15 / delta;
}

array<double, 2> garch21::nu_dist::draw(void)
{
    array<double, 2> proposed;
    uniform_real_distribution<double> unif(0.0, 1.0);
    double d1, d2, u;
    do {
	proposed = proposal_draw();
	d2 = proposal_density(proposed);
	d1 = density(proposed);
	u = unif(markov->dev);
    } while (u > d1 / d2 / c1 / c2 / c3);
    return proposed;
}

array<double, 2> garch21::forward(double x0, bool orig)
{
    // move forward according to the normal transition kernel if orig = 1
    array<mat, 2> A;
    chi_squared_distribution<double> chi2;
    vec V(2);
    array<double, 2> ret;

    V(0) = x0;
    V(1) = 1 - x0;

    V.print();
    for (size_t i = 0; i < A.size(); i++) {
	double z2 = chi2(dev);
	A[i](0, 0) = a1 * z2 + b1;
	A[i](0, 1) = a2;
	A[i](1, 0) = z2;
	A[i](1, 1) = 0;
    }
    A[0].print();
    A[1].print();
    
    V = A[1] * A[0] * V;
    V.print();
    
    ret[1] = sum_norm(V);
    ret[0] = V(0) / ret[1];

    bool accepted = orig;
    while (! accepted) {
	double d1 = kernel_density(ret, x0);
	double d2 = nu.delta * nu.density(ret);
	uniform_real_distribution<double> unif;
	accepted = unif(dev) < 1 - d2/d1;
	if (! accepted)
	    ret = forward(x0);
    }
    return ret;
    
}

