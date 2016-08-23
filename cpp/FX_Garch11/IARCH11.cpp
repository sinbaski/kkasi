#include <math.h>
#include <array>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>

#define ARRAY_SIZE(X) sizeof(X)/sizeof((X)[0])

double f (double x, void *params)
{
    double (*par)[2] = (double (*)[2])params;
    double t1 = sin(M_PI * x) + fabs((*par)[0]);
    return pow(t1 * t1, (*par)[1]);
}

double f_integral(double (*params)[2])
{
    gsl_integration_workspace * w 
	= gsl_integration_workspace_alloc (1000);
  
    double result, error;
    gsl_function F;
    F.function = &f;
    F.params = params;
    gsl_integration_qag(&F, -1, 1, 0, 1.0e-7, 1000,
			GSL_INTEG_GAUSS61,
			w, &result, &error);
    gsl_integration_workspace_free (w);
    return result;
}

double target_fun(double xi, void *par)
{
    double rho = *(double *)par;
    double X[2];
    X[0] = rho;
    X[1] = xi;
    double u = f_integral(&X);
    double v = 2/gsl_sf_gamma(2 * xi + 1);
    return  u - v;
}

double tail_index_fun(double rho)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub, xi;

    F.function = target_fun;
    F.params = &rho;
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, 0.5, 1.5);
    
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 1.0e-6, 0);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) return -1;
    return xi;
}

int main(int argc, char **argv)
{
    // gsl_integration_workspace * w 
    // 	= gsl_integration_workspace_alloc (1000);
  
    // double result, error;
    // gsl_function F;
    // double params[] = {0.1, 0.8};
    // F.function = &f;
    // F.params = &params;

    // gsl_integration_qag(&F, -1, 1, 0, 1.0e-7, 1000, GSL_INTEG_GAUSS61,
    // 			w, &result, &error);
    // printf ("result          = % .18f\n", result);
    // printf ("estimated error = % .18f\n", error);
    // printf ("intervals       = %zu\n", w->size);

    // gsl_integration_workspace_free (w);

    const size_t N = 1000;
    double X[N + 1];
    double Y[N + 1];
    double dx = 1/(double)N;
    X[0] = 0;
    Y[0] = 1;
    for (size_t i = 1; i < ARRAY_SIZE(X); i++) {
	X[i] = i * dx;
	Y[i] = tail_index_fun(X[i]);
	printf("%e\t%e\n", X[i], Y[i]);
    }
    return 0;
}
