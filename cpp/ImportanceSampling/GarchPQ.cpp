#include <stdio.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace arma;

random_device gen;


mat& gen_rand_matrix(const vector<double> &alpha,
		     const vector<double> &beta,
		     double z2, mat &M)
{
    unsigned p = beta.size(), q = alpha.size() - 1;
    unsigned n = p + q - 1;

    M.zeros(n, n);
    M(0, 0) = alpha[1] * z2 + beta[0];
    for (unsigned i = 1; i < p; i++) {
	M(0, i) = beta[i];
    }
    for (unsigned i = p; i < n; i++) {
	M(0, i) = alpha[i - p + 2];
    }
    for (unsigned i = 1; i < p; i++) {
	M(i, i - 1) = 1;
    }
    if (n > 1) M(p, 0) = z2;
    for (unsigned i = p + 1; i < n; i++) {
	M(i, i - 1) = 1;
    }
    return M;
}

double estimateLambda(const vector<double>& a,
		      const vector<double>& b,
		      double theta,
		      unsigned long N,
		      unsigned long K,
		      double &sd)
{

    uniform_real_distribution<double> unif;
    chi_squared_distribution<double> chi2;
    // vector<double> beta(N, 0);
    vector<vec> E(K);
    vector< array<double, 2> > alpha(K);
#pragma omp parallel for schedule(dynamic) shared(gen, unif)
    for (unsigned i = 0; i < K; i++) {
	E[i].set_size(a.size() + b.size() - 2);
	for_each(E[i].begin(), E[i].end(),
		 [&](double &x)
		 {
		     x = unif(gen);
		 });
    }
    double Lambda = 0;
    sd = 0;
    for (unsigned j = 0; j < N; j++) {
	double beta = 0;
	vector<mat> A(K);
	vector<double> Q(K);
	// transform(alpha.begin(), alpha.end(), Q.begin(),
	// 	  [&s](const array<double, 2> &a) {
	// 	      return s += a[0];
	// 	  });
	if (j == 0) {
#pragma omp parallel for
	    for (unsigned k = 0; k < K; k++) {
		alpha[k][0] = 1;
		Q[k] = k + 1;
	    }
	    beta = K;
	} else {
	    for (unsigned k = 0; k < K; k++) {
		if (k == 0) {
		    Q[k] = alpha[k][0];
		} else {
		    Q[k] = Q[k-1] + alpha[k][0];
		}
		beta += alpha[k][0];
		alpha[k][0] = alpha[k][1];
	    }
	}

#pragma omp parallel for schedule(dynamic) shared(gen, chi2, beta)
	for (unsigned k = 0; k < K; k++) {
	    double z2;

#pragma omp critical
	    z2 = chi2(gen);
	    gen_rand_matrix(a, b, z2, A[k]);
	    alpha[k][1] = pow(norm(A[k] * E[k], "inf"), theta);
	}
	vector<vec> E_prev(K);
	copy(E.begin(), E.end(), E_prev.begin());

#pragma omp parallel for shared(gen, unif)
	for (unsigned k = 0; k < K; k++) {
	    double U;

#pragma omp critical
	    U = unif(gen) * beta;
	    unsigned l = upper_bound(Q.begin(), Q.end(), U) - Q.begin();
	    vec V = A[k] * E_prev[l];
	    E[k] = V / norm(V, "inf");
	}

	double lbt = log(beta/K);
	Lambda += lbt / N;
	sd += pow(lbt, 2) / N;
	if (j == N) {
	    accumulate(alpha.begin(), alpha.end(), beta = 0,
		       [](double s, const array<double, 2>& a) {
			   return s + a[1];
		       });
	    lbt = log(beta/K);
	    Lambda += lbt / N;
	    sd += pow(lbt, 2)/N;
	}
    }
    sd = sqrt(sd - pow(Lambda, 2));
    return Lambda;
}

struct func_par
{
    const vector<double>& a;
    const vector<double>& b;
    unsigned long N;
    unsigned long K;
};

double func(double theta, void *p)
{
    struct func_par* par = (struct func_par*) p;
    double sd;
    return estimateLambda(par->a, par->b, theta, par->N, par->K, sd);
}

double find_root(const vector<double>& a,
		 const vector<double>& b,
		 unsigned long N,
		 unsigned long K,
		 double bounds[2])
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi;
    struct func_par par = {
	a, b, N, K
    };
    F.function = func;
    F.params = &par;
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
	if (status == GSL_SUCCESS)
	    cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    if (status != GSL_SUCCESS) {
	cout << "The Brent algorithm did not converge after " << max_iter
	     << " iterations." << endl;
	xi = -1;
    }
    return xi;
}

int main(int argc, char*argv[])
{
    // // madeup
    // vector<double> alpha({1.0e-7, 0.11, 0});
    // vector<double> beta({0.88});
    // // DAX
    // vector<double> alpha({1.0e-7, 0.02749864, 0.04228535});
    // vector<double> beta({0.8968533});
    // FTSE100
    vector<double> alpha({1.0e-7, 0.10623464, 0.02904907});
    vector<double> beta({0.7829784});

    double Lambda;
    
    cout << "alpha[0]= "  << alpha[0] << ", alpha[1]=" <<
    	alpha[1] << ", alpha[2]=" << alpha[2] <<
    	", beta[1]=" << beta[0] << endl;
    cout << "N = " << argv[1] << endl;
    cout << "K = " << argv[2] << endl;

    unsigned long N = stoul(argv[1]), K = stoul(argv[2]);
    
    double bounds[2];
    bool flag = false;
    for (double nu = stod(argv[3]); nu <= stod(argv[4]); nu += 0.1) {
	double sd;
	Lambda = estimateLambda(alpha, beta, nu, N, K, sd);
	printf("%.4f    % .4f    %.4f    %.4f\n", nu, Lambda, sd, sd/abs(Lambda));
	if (!flag && Lambda > 0) {
	    bounds[1] = nu;
	    bounds[0] = nu - 0.1;
	    flag = true;
	}
    }
    // double xi = find_root(alpha, beta, N, K, bounds);
    // cout << "Lambda(" << xi << ") = 0" << endl;
    return 0;
}


