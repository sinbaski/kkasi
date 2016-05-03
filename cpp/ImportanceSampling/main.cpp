#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
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


struct lambda_minimizer_data
{
    double beta;
    double alpha_bounds[2];
    double alpha;
    double inf;
};

struct minimizer_func_par
{
    double beta;
    const Garch1D<double> *garch11;
};

double lambda_minimizer_func(double alpha, void *par)
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

int find_inf_lambda (
    struct lambda_minimizer_data *mdata,
    int data_size,
    const Garch1D<double>* garch11
    )
{
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    gsl_min_fminimizer *minimizer = gsl_min_fminimizer_alloc(
	gsl_min_fminimizer_brent
	);

    for (int i = 0; i < data_size; i++) {
	struct minimizer_func_par par;
	gsl_function F;
	par.beta = mdata[i].beta;
	par.garch11 = garch11;
	F = {
	    lambda_minimizer_func,
	    &par
	};
	gsl_min_fminimizer_set(
	    minimizer, &F, mdata[i].alpha,
	    mdata[i].alpha_bounds[0], 
	    mdata[i].alpha_bounds[1]
	    );
	iter = 0;
	do {
	    iter++;
	    status = gsl_min_fminimizer_iterate(minimizer);
	    mdata[i].alpha = gsl_min_fminimizer_x_minimum(minimizer);
	    lb = gsl_min_fminimizer_x_lower (minimizer);
	    ub = gsl_min_fminimizer_x_upper (minimizer);

	    status = gsl_min_test_interval (lb, ub, 0.0001, 0.0);

	    if (status == GSL_SUCCESS) {
		mdata[i].inf = gsl_min_fminimizer_f_minimum(minimizer);
		printf ("Minimizer Converged for beta = %.4f:"
			" alpha = %.4f, "
			"lambda = %.4f\n", mdata[i].beta,
			mdata[i].alpha, mdata[i].inf
			);
	    }
	} while (status == GSL_CONTINUE && iter < max_iter);
	if (status != GSL_SUCCESS) {
	    printf ("Minimizer didn't converge for beta 0%.4f.\n",
		    mdata[i].beta);
	    gsl_min_fminimizer_free(minimizer);
	    break;
	}
    }
    gsl_min_fminimizer_free(minimizer);
    return status;
}

int main(int argc, char* argv[])
{
    /**
     * Simulate GARCH(1,1) processes.
     */
    double coef[] = {
	1.0e-7, 0.11, 0.88
    };
    Garch1D<double> garch11(coef[0], coef[1], coef[2], stol(argv[1]));
    // double moment = stod(argv[1]);
    // cout << "E A^(" <<  moment << ") = "
    // 	 << garch11.moment_func(moment)
    // 	 << endl;

    printf("lambda(-%.4f) = %.4f\n", garch11.xi, garch11.moment_func(-garch11.xi, 0, 0));

    /**
     * Find M. The set C = [-M, M]
     */
    double M = -numeric_limits<double>::infinity();
    struct lambda_minimizer_data mdata[] = {
	// beta = -xi
	{ 0, {2, 3.5}, 2.75, 0},
	// beta = 0
	{0, {0, 0}, 0, 0}
    };
    mdata[0].beta = -garch11.xi;
    mdata[1].beta = 0;
    mdata[1].alpha_bounds[1] = garch11.xi;
    mdata[1].alpha = garch11.xi/2;
    find_inf_lambda(mdata, 2, &garch11);
    for (int i = 0; i < 2; i++) {
	double x;
	if (mdata[i].alpha <= 1) {
	    x = garch11.a0 / pow(1 - mdata[i].inf, 1/mdata[i].alpha);
	} else {
	    x = garch11.a0 / (1 - pow(mdata[i].inf, 1/mdata[i].alpha));
	}
	M = max(M, x);
    }
    garch11.M = M;

    /**
     * Simulate the process
     */
    double u = stod(argv[3]);
    vector<double> sim_stat(stoi(argv[2]));
#pragma omp parallel for
    for (vector<double>::iterator i = sim_stat.begin(); i < sim_stat.end(); i++) {
	uniform_real_distribution<double> unif;
	random_device alice;
	double V = garch11.init_quantile(unif(alice));
	long Nu = 0;
	vector<double> As;
	// As.clear();
	int measure_index = Garch1D<double>::SHIFTED_M;
	int status = 0;
	while (status <= 2) {
	    switch(status) {
		double A;
	    case 0: //before exceeding u
		do {
		    A = garch11.quantile_func(unif(alice), measure_index);
		    As.push_back(A);
		    V = A * V + garch11.a0;
		} while (V < u && V > M);
		if (V <= M) { // Restart the process
		    As.clear();
		    Nu = 0;
		    V = garch11.init_quantile(unif(alice));
		    continue;
		} else {
		    Nu++;
		    status = 1;
		    // Sample from the original measure from now on
		    measure_index = Garch1D<double>::ORIG_M;
		}
		break;
	    case 1: // u exceeded, still outside C
		A = garch11.quantile_func(unif(alice), measure_index);
		V = A * V + garch11.a0;
		if (V <= M) {
		    status = 2;
		} else if (V > u) {
		    Nu++;
		}
		break;
	    case 2: // u exceeded, has returned to C
		*i = accumulate(As.begin(), As.end(), (double)Nu,
			       [&garch11](double s, double a) {
				   return s * pow(a, -garch11.xi);
			       }
		    );
		status++;
	    }
	}
    }

    double estimator[] = {0, 0};
    estimator[0] = accumulate(sim_stat.begin(), sim_stat.end(), estimator[0]);
    estimator[0] *= garch11.stationary_prob(garch11.M)/(double)sim_stat.size();
    
    estimator[1] = accumulate(sim_stat.begin(), sim_stat.end(), estimator[1],
			      [=](double s, double x) {
				  return s + pow(x - estimator[0], 2)/(double)sim_stat.size();
			      });
    estimator[1] = sqrt(estimator[1]);
    printf("M = %.4e, measure(C) = %.4f\n", M, garch11.stationary_prob(M));
    printf("P(V > %.2f) = %.4e( %.4e ), std/mean = %.2f\n", u,
	   estimator[0], estimator[1], estimator[1]/estimator[0]);

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
