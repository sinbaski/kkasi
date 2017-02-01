#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <array>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

using namespace std;

#define SIZE 10

double g(double z, void *params)
{
    vector<double>* par = (vector<double>*)params;
    double rho = par->at(0);
    double xi = par->back();
    return pow(abs(sin(M_PI * z) + rho), 2*xi);
}

double g_integral(vector<double> *params)
{
    gsl_integration_workspace *w 
	= gsl_integration_workspace_alloc (1000);
  
    double result, error;
    gsl_function G;
    G.function = &g;
    G.params = params;    
    gsl_integration_qag(&G, -1, 1, 0, 1.0e-4, 1000, GSL_INTEG_GAUSS61,
			w, &result, &error);
    return result;
}

double g_integral_func(double kappa, void *params)
{
    vector<double> *par = (vector<double>*)params;
    vector<double> V(*par);
    V.push_back(kappa);
    // product of alpha_1 and alpha_2
    double a = par->at(1);
    double r = g_integral(&V) -
	2 / gsl_sf_gamma(2 * kappa + 1) / pow(a, kappa);
    return r;
}

double find_tail_index(double rho, double alpha_prod, double bounds[2])
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi;
    vector<double> params(2);
    params[0] = rho;
    params[1] = alpha_prod;
    F.function = g_integral_func;    
    F.params = &params;
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    assert(status == GSL_SUCCESS);
    return xi;
}

double f(double xi, void* param)
{
    double alpha = *(double *)param;
    double t1 = gsl_sf_gamma(xi + 0.5) * pow(2, xi) / sqrt(M_PI);
    t1 =  pow(t1, -1.0/xi);
    return t1 - alpha;
}

double arch_tail_index(double alpha)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi;
    F.function = f;    
    F.params = &alpha;

    xi = 2;
    double t = f(xi, &alpha);
    if (t > 0) {
	do {
	    xi *= 2;
	    t = f(xi, &alpha);
	} while (t > 0);
	ub = xi;
	lb = xi / 2;
    } else {
	do {
	    xi /= 2;
	    t = f(xi, &alpha);
	} while (t < 0);
	ub = xi * 2;
	lb = xi;
    }
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, lb, ub);
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    assert(status == GSL_SUCCESS);
    return xi;
}

int main(int argc, char *argv[])
{
    for (double rho = 0.1; rho <= 1; rho += 0.1) {
    	for (double a = 0.09; a < 0.99; a += 0.1) {
    	    double kappa = 2;
    	    vector<double> params({rho, a});
    	    double bounds[2];
    	    if (g_integral_func(kappa, &params) < 0) {
    		while (g_integral_func(kappa, &params) < 0) {
    		    kappa *= 2;
    		}
    		bounds[0] = kappa / 2;
    		bounds[1] = kappa;
    	    } else {
    		while (g_integral_func(kappa, &params) > 0)
    		    kappa /= 2;
    		bounds[0] = kappa;
    		bounds[1] = kappa * 2;
    	    }
    	    kappa = find_tail_index(rho, a, bounds);
    	    printf("%6.3f", kappa * 2);
    	}
    	printf("\n");
    }
   return 0;

    // array<double, 10> uni_indices;
    // int i = 0;
    // for (double a = 0.09; a < 0.99; a += 0.1) {
    // 	printf("%6.3f", a);
    // 	uni_indices[i++] = arch_tail_index(sqrt(a));
    // }
    // cout << endl;
    // for_each(uni_indices.begin(), uni_indices.end(),
    // 	     [](double x){
    // 		 printf("%6.3f", x);
    // 	     });
    // printf("\n");
    // return 0;
}
