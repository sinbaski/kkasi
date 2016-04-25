#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
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

int main(int argc, char* argv[])
{
    /**
     * Simulate GARCH(1,1) processes.
     */
    double coef[] = {
	0.11, 0.88
    };
    Garch1D<double> garch11(coef[0], coef[1]);
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
