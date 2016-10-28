#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <gsl/gsl_math.h>
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
    return pow(sin(M_PI * z) + rho, 2*xi);
}

double g_integral(vector<double> *params)
{
    gsl_integration_workspace *w 
	= gsl_integration_workspace_alloc (1000);
  
    double result, error;
    gsl_function G;
    G.function = &g;
    G.params = params;    
    gsl_integration_qag(&G, -1, 1, 0, 1.0e-7, 1000, w, &result, &error);
    return result;
}

double g_integral_func(double kappa, void *params)
{
    vector<double> *par = (vector<double>*)params;
    vector<double> V(*par);
    V.push_back(kappa);
    // product of alpha_1 and alpha_2
    double a = par->at(1);
    g_integral(V) - 2 / gsl_sf_gamma(2 * kappa + 1) / pow(a, kappa);
}

double find_tail_index(double rho, doube alpha_prod, double bounds[2])
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

int main(int argc, char *argv[])
{
    for (double rho = 0.1; rho <= 1; rho++) {
	for (double a = 0.09; a < 0.99; a++) {
	    double xi = 0.1;
	    while ()
	}
    }


    // for (unsigned i = 0; i < field.size(); i++) {
    //     for (unsigned j = 0; j < field[0].size(); j++) {
    //         //printf("%8.4f", field[i][j]);
    // 	    printf("%6ld", field[i][j].size());
    //     }
    // 	printf("\n");
    // }    
    return 0;
}
