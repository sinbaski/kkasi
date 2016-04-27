#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include "garch1D.hpp"

using namespace std;
using namespace arma;

#define A_DIM 4

void print_matrix(const mat& M)
{
    for (unsigned i = 0; i < M.n_rows; i++) {
	for (unsigned j = 0; j < M.n_cols; j++) {
	    printf("% 12.2e", M(i,j));
	}
	printf("\n");
    }
    printf("\n");
}

/**
   alpha: a vector of size 2+ that contain the coefficients of past squared returns.
   beta: a vector of size 1+ that contain the coefficients of past variances.
   This notation is consistent with Mikosch and Davis'es paper. on GARCH(p,q)
 */
mat& gen_rand_matrices(const vector<double> &alpha,
		       const vector<double> &beta,
		       mat& A)
{
    static chi_squared_distribution<double> dist;
    static random_device gen;
    
    double z = dist(gen);
    A.row(0) = rowvec({alpha[0] * z + beta[0], beta[1], alpha[1], alpha[2]});
    A.row(1) = rowvec({1, 0, 0, 0});
    A.row(2) = rowvec({z, 0, 0, 0});
    A.row(3) = rowvec({0, 0, 1, 0});
    return A;
}

double norm_pow_prod_mat(const vector<mat>& matrices, double alpha)
{
    mat x(A_DIM, A_DIM);
    x.eye();
    for (size_t i=0; i < matrices.size(); i++) {
	x *= matrices[i];
    }
    // double y = norm(x);
    // vec S = svd(x);
    // double* z = S.memptr();
    return pow(norm(x), alpha);
}

double mean_norm_pow_prod_mat(const vector< vector<mat> >& matrices,
			      double alpha)
{
    size_t p = matrices.size();
    double x = 0;
    size_t i;
    for (x = 0, i = 0; i < p; i++) {
	x += norm_pow_prod_mat(matrices[i], alpha)/p;
    }
    return x;
}

// double test(void)
// {
//     static normal_distribution<double> dist;
//     // static default_random_engine gen;
//     static random_device gen;
//     return dist(gen);
// }

double Lambda_func(double xi, unsigned int n, unsigned int p,
		   const vector<double> &alpha, const vector<double> &beta)
{
    vector< vector<mat> > matrices;
    mat M(A_DIM, A_DIM);
    matrices.resize(p);
    for (unsigned i = 0; i < p; i++) {
    	matrices[i].resize(n);
    	for (unsigned j = 0; j < n; j++) {
    	    matrices[i][j] = gen_rand_matrices(alpha, beta, M);
    	}
    }
    double x = mean_norm_pow_prod_mat(matrices, xi);
    return log(x)/double(n);
}

double LHS_func(double xi, void *par)
{
    Garch1D<double> *garch11 = (Garch1D<double> *)par;
    return garch11->moment_func(xi) - 1;
}

struct minimizer_func_par
{
    double beta;
    Garch1D<double> *garch11;
};

double M_minimizer_func(double alpha, void *par)
{
    static double normalizer = 0;
    static bool flag = true;

    struct minimizer_func_par *p = (struct minimizer_func_par *)par;
    if (flag && abs(p->beta) > 1.0e-3)
	normalizer = p->garch11->moment_func(p->beta, 0, 0);
    else if (flag)
	normalizer = 1;
    flag = false;
    double lambda = p->garch11->moment_func(alpha, p->beta, normalizer);
    // if (alpha < 1)
    // 	return garch11->a0 / pow(1 - lambda, 1/alpha);
    // else
    // 	return garch11->a0 / (1 - pow(lambda, 1/alpha));
    return lambda;
}

int main(int argc, char* argv[])
{
    /**
     * Simulate GARCH(1,1) processes.
     */
    double coef[] = {
	1.0e-7, 0.11, 0.88
    };
    Garch1D<double> garch11(coef[0], coef[1], coef[2]);
    // double moment = stod(argv[1]);
    // cout << "E A^(" <<  moment << ") = "
    // 	 << garch11.moment_func(moment)
    // 	 << endl;

    /**
     * Find the tail index
     */
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double xi, lb, ub;
    F.function = &LHS_func;
    F.params = &garch11;
    
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, 1.8, 1.9);

    printf ("using %s method\n", 
	    gsl_root_fsolver_name (solver));

    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0001, 0);
	if (status == GSL_SUCCESS)
	    cout << "Converged: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) {
	cout << "The Brent algorithm did not converge after " << max_iter
	     << " iterations." << endl;
	return status;
    }
    garch11.set_shift_par(xi);

    printf("lambda(-%.4f) = %.4f\n", xi, garch11.moment_func(-xi, 0, 0));

    /**
     * Find the set C = [-M, M]
     * Choose the smallest M_beta(alpha) given beta.
     * beta is either 0 or -xi. Choose the larger M in these two cases.
     */
    gsl_min_fminimizer *minimizer = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    struct minimizer_func_par par;
    par.beta = -garch11.get_shift_par();
    par.garch11 = &garch11;

    F = {
	M_minimizer_func,
	&par
    };
    double arginf;

    gsl_min_fminimizer_set(
	minimizer, &F, 3, 2, 3.5
	);
    iter = 0;
    do {
	iter++;
	status = gsl_min_fminimizer_iterate(minimizer);
	arginf = gsl_min_fminimizer_x_minimum(minimizer);
	lb = gsl_min_fminimizer_x_lower (minimizer);
	ub = gsl_min_fminimizer_x_upper (minimizer);

	status = gsl_min_test_interval (lb, ub, 0.0001, 0.0);

      if (status == GSL_SUCCESS)
	  printf ("Minimizer Converged: alpha = %.4f, "
		  "lambda = %.4f\n", arginf,
		  gsl_min_fminimizer_f_minimum(minimizer));
    } while (status == GSL_CONTINUE && iter < max_iter);
    if (iter == max_iter && status == GSL_CONTINUE) {
	printf ("Minimizer didn't converge.\n");
    }
    

    gsl_min_fminimizer_free(minimizer);
    
    /**
     * Evaluate the Lambda(.) function by simulation
     */
    // unsigned n = stoi(argv[1]);
    // unsigned p = stoi(argv[2]);
    // double xi = stod(argv[3]);
    // vector<double> alpha({0.6, 0.1, 0.1});
    // vector<double> beta({0.1, 0.05});

    // for (int i=0; i < stoi(argv[4]); i++) {
    // 	double x = Lambda_func(xi, n, p, alpha, beta);
    // 	cout<< x << endl;
    // }


    /**
     * Test random number generation.
     */
    // for (int i = 0; i < stoi(argv[1]); i++) {
    // 	cout << test() << endl;
    // }
    
    return 0;
}
