#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <random>

using namespace std;

#define SAMPLE_SIZE 400000

double pool[SAMPLE_SIZE];

double fun(double xi, void *param)
{
    double *coef =(double *)param;
    double r;
    r = accumulate(pool, pool + SAMPLE_SIZE, (double)0,
		   [xi, coef](double s, double x) {
		       return s + pow(coef[1] * x + coef[2], xi)/SAMPLE_SIZE;
		   }
	);
    return r - 1;
}

double tail_index(double coef[3], double interval[2])
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub, xi = -1;
    F.function = fun;
    F.params = coef;

    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, interval[0], interval[1]);

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

int main(int argc, char *argv[])
{
    double coef[3];
    random_device gen;
    chi_squared_distribution<double> dist(1.0);
    
    coef[0] = 1.0e-7;
    coef[1] = stod(argv[1]);
    coef[2] = stod(argv[2]);

//#pragma omp parallel for
    for (int i = 0; i < SAMPLE_SIZE; i++) {
	pool[i] = dist(gen);
    }
 
   double bounds[] = {1, 6};
   double a;
   for (a = 2; a > 0 && fun(a, coef) > 0; a -= 1);
   if (a > 0) {
       bounds[0] = a;
   } else {
       printf("%s %.2f.\n", "lower bound less than ", a);
       return 0;
   }
   for (a = 2; a < 10 && fun(a, coef) < 0; a += 1);
   if (a < 10) {
       bounds[1] = a;
   } else {
       printf("%s %.2f.\n", "Upper bound larger than ", (double)10);
       return 0;
   }
   double xi = tail_index(coef, bounds);
   cout << xi << endl;
}
