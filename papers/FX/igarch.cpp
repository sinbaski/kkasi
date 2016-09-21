#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <random>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

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
    memcpy(alpha, par->alpha, sizeof(alpha));

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

double f_integral(double ub, double kappa, double rho, double* alpha, double* err)
{
    double res;
    gsl_rng_env_setup();
    gsl_rng *gen = gsl_rng_alloc(gsl_rng_default);

    struct f_params params;
    
    params.kappa = kappa;
    params.rho = rho;
    memcpy(params.alpha, alpha, sizeof(params.alpha));
    
    gsl_monte_function F = { &f, 2, &params};
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
    double lim1[] = {0, 0};
    double lim2[] = {ub, 2 * M_PI};
    
    gsl_monte_vegas_integrate (&F, lim1, lim2, 2, 10000, gen, s,
                               &res, err);
    do {
        gsl_monte_vegas_integrate (&F, lim1, lim2, 2, 100000, gen, s,
                                   &res, err);
        // printf ("result = % .6f sigma = % .6f "
        //         "chisq/dof = %.1f\n", res, *err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
    gsl_monte_vegas_free(s);
    return res;
    
}

int main(int argc, char *argv[])
{
    double err;
    double alpha[] = {1, 1};
    double s = f_integral(10, 0.8, 1, alpha, &err);
    printf("%e    %e\n", s, err);
    return 0;
}
