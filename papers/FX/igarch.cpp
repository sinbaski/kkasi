#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;

#define SIZE 4
extern double estimateLambda(const vector<double>& a,
			     const vector<double>& b,
			     double theta,
			     unsigned long N,
			     unsigned long K,
			     double &sd);

extern double find_root(const vector<double>& a,
			const vector<double>& b,
			unsigned long N,
			unsigned long K,
			double bounds[2]);

struct f_params
{
    double kappa;
    double alpha[2];
    double rho;
};

double f(double *x, size_t dim, void *params)
{
    struct f_params *par = (struct f_params *)params;
    double kappa = par->kappa;
    double rho = par->rho;
    double alpha[2];
    copy(par->alpha, par->alpha + 2, alpha);

    double r = x[0];
    double theta = x[1];

    double t1 = cos(theta);
    double t2 = t1 * t1;
    double t4 = r * r;
    double t7 = pow(t4 * t2 * alpha[0] - alpha[0] + 0.1e1, kappa);
    double t9 = rho * rho;
    double t11 = sqrt(-t9 + 1);
    double t13 = sin(theta);
    double t20 = t13 * t13;
    double t27 = pow(t4 * (0.2e1 * alpha[1] * rho * t13 * t1 * t11 + alpha[1] * t9 * t2
		    - alpha[1] *  t9 * t20 + alpha[1] * t20) + 0.1e1 - alpha[1], kappa);
    double t29 = exp(-t4 / 0.2e1);
    double t31 = t29 * t27 * t7 * r;
    return t31;
 }

double f_integral(double kappa, double rho, double* alpha, double eps, double* err)
{

    gsl_rng_env_setup();
    gsl_rng *gen = gsl_rng_alloc(gsl_rng_default);
    struct f_params params;
    
    params.kappa = kappa;
    params.rho = rho;
    copy(alpha, alpha + 2, params.alpha);
    
    gsl_monte_function F = { &f, 2, &params};
    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(2);

    double M = 10;
    double integral = 0, inc, dM = 5;
    do {
	double res;
	double lim1[] = {0, 0};
	double lim2[] = {M, 2 * M_PI};
    
	// Warm up...
	gsl_monte_vegas_integrate (&F, lim1, lim2, 2, 10000, gen, state,
				   &res, err);
	do {
	    gsl_monte_vegas_integrate (&F, lim1, lim2, 2, 100000, gen, state,
				       &res, err);
	} while (fabs (gsl_monte_vegas_chisq (state) - 1.0) > 0.5);
	inc = res - integral;
	integral = res;
	M += dM;
    } while (inc / integral > eps);
    gsl_monte_vegas_free(state);
    return integral / 2 / M_PI;
}

struct func_par
{
    double rho;
    double *alpha;
};

double func(double theta, void *p)
{
    struct func_par* par = (struct func_par*) p;
    double err;
    return f_integral(theta, par->rho, par->alpha, 1.0e-4, &err) - 1;
}

struct corr_func_par
{
    double kappa;
    double *alpha;
};

double corr_func(double rho, void *p)
{
    struct corr_func_par* par = (struct corr_func_par*) p;
    double err;
    return f_integral(par->kappa, rho, par->alpha, 1.0e-4, &err) - 1;
}

void map_tail_index(array<double, 4>& field)
{
    unsigned long N = 1000;
    unsigned long K = 4000;

    for (unsigned i = 0; i < field.size(); i++) {
	vector<double> alpha(1, (double)(i + 1) * 0.05);
	vector<double> beta(1, 1 - alpha[0]);
	double sd;
	double theta = 0.5;
	double Lambda = estimateLambda(alpha, beta, theta, N, K, sd);
	double inc = Lambda < 0 ? 1 : -1;
	double bounds[2];
	while(theta < 10 && Lambda * inc < 0) {
	    inc *= 2;
	    Lambda = estimateLambda(alpha, beta, theta += inc, N, K, sd);
	}
	if (theta >= 10 && Lambda * inc < 0) {
	    field[i] = GSL_POSINF;
	    return;
	}
	bounds[0] = min(theta, theta - inc);
	bounds[1] = max(theta, theta - inc);
	field[i] = find_root(alpha, beta, N, K, bounds);
    }
}

double product_tail_index(double alpha[2], double rho)
{
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    gsl_root_fsolver *solver;
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);

    double xi = 0.5;
    double err;
    // first find an interval in which the solution lies
    struct func_par par = {
	rho, alpha
    };
    double integral = f_integral(xi, rho, alpha, 1.0e-4, &err);
    double inc = integral > 1 ? -0.05 : 0.05;
    double bounds[2];
    while (inc * (integral - 1) < 0) {
	inc *= 2;
	integral = f_integral(xi += inc, rho, alpha, 1.0e-4, &err);
    }
    bounds[0] = min(xi, xi - inc);
    bounds[1] = max(xi, xi - inc);

    F.function = func;
    F.params = &par;
    gsl_root_fsolver_set(solver, &F,
			 min(bounds[0], bounds[1]),
			 max(bounds[0], bounds[1]));
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    assert(status == GSL_SUCCESS);
    gsl_root_fsolver_free(solver);
    return xi;
}

void map_product_tail_index(array<array<double, SIZE>, SIZE>& field, double rho)
{
    for (unsigned i = 0; i < field.size(); i++) {
#pragma omp parallel for schedule(dynamic)
	for (unsigned j = i; j < field[0].size(); j++) {
	    double alpha[2];
	    alpha[0] = (double)(i + 1) * (double)1/SIZE;
	    alpha[1] = (double)(j + 1) * (double)1/SIZE;
	    field[i][j] = field[j][i] = product_tail_index(alpha, 0.5);
	}
    }
}

int main(int argc, char *argv[])
{
    array<array<vector<double>, SIZE>, SIZE> field;
    array<double, SIZE> uni_indices;
 
    map_tail_index(uni_indices);
    for (unsigned i = 0; i < field.size(); i++) {
        for (unsigned j = 0; j < field[0].size(); j++) {
	    double xi = min(min(uni_indices[i], uni_indices[j])/2, 10.0);
	    
//	    double rho;
	    double alpha[2];
	    double err;
	    alpha[0] = (double)(i + 1) * (double)1/SIZE;
	    alpha[1] = (double)(j + 1) * (double)1/SIZE;
	    array<double, 20> expectaions;
	    for (unsigned k = 0; k < expectaions.size(); k++) {
		expectaions[k] = f_integral(xi, k * 0.1 - 1, alpha, 1.0e-4, &err);
		if (k > 0 && (expectaions[k] - 1) * (expectaions[k-1] - 1) < 0) {
		    field[i][j].push_back(k);
		}
	    }
	    
        }
    }    

    for (unsigned i = 0; i < field.size(); i++) {
        for (unsigned j = 0; j < field[0].size(); j++) {
            //printf("%8.4f", field[i][j]);
    	    printf("%6ld", field[i][j].size());
        }
    	printf("\n");
    }    
    return 0;
}
