#include <stdio.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>

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
		      unsigned long K)
{

    uniform_real_distribution<double> unif;
    chi_squared_distribution<double> chi2;
    vector<double> beta(K, 0);
    vector<vec> E(K);
    for_each(
	E.begin(), E.end(),
	[&](vec &e){
	    e.set_size(2);
	    for_each(e.begin(), e.end(),
		     [&](double &x) {
			 x = unif(gen);
		     });
	});
    for (unsigned j = 1; j < N; j++) {
	vector<double> alpha(K);
	vector<mat> A(K);
	double s = 0;
//#pragma omp parallel for schedule(dynamic) shared(gen, chi2, s)
	for (unsigned k = 1; k < K; k++) {
	    double z2;

//#pragma omp critical
	    z2 = chi2(gen);
	    gen_rand_matrix(a, b, z2, A[k]);
	    alpha[k] = pow(norm(A[k] * E[k], "inf"), theta);
//#pragma omp atomic
	    s += alpha[k];
	}
        beta[j] = s;

	vector<double> Q(K);
	partial_sum(alpha.begin(), alpha.end(), Q.begin());

//#pragma omp parallel for shared(gen, unif)
	for (unsigned k = 1; k < K; k++) {
	    double U;

//#pragma omp critical
	    U = unif(gen);
	    U = U * s;
	    unsigned l = upper_bound(Q.begin(), Q.end(), U) - Q.begin();
	    vec V = A[k] * E[l];
	    E[k] = V / norm(V, "inf");
	}	
    }

    double Lambda = accumulate(beta.begin(), beta.end(), 0.0,
	       [K, N](double s, double x)
	       {
	           return s + log(x/K) / N;
	       }
    );
    return Lambda;
}


int main(int argc, char*argv[])
{
    vector<double> alpha({1.0e-7, 0.6, 0.001});
    // vector<double> alpha({1.0e-7, stod(argv[3])});
    vector<double> beta({0.005});
    double Lambda;
    
    cout << "alpha[0]= "  << alpha[0] << ", alpha[1]=" <<
    	alpha[1] << ", alpha[2]=" << alpha[2] <<
    	", beta[1]=" << beta[0] << endl;
    cout << "N = " << argv[1] << endl;
    cout << "K = " << argv[2] << endl;

    // double nu = stod(argv[1]);
    // Lambda = estimateLambda(
    // 	alpha, beta, nu, stoul(argv[2]),
    // 	stoul(argv[3]), bounds);
    // printf("%e    %e\n", nu, Lambda);

    double nu = 1.05;
    Lambda = estimateLambda(
	alpha, beta, nu, stoul(argv[1]), stoul(argv[2]));
    printf("%e    %e\n", nu, Lambda);

    // for (double nu = stod(argv[1]); nu < 2; nu += 0.1) {
    // 	Lambda = estimateLambda(
    // 	    alpha, beta, nu, stoul(argv[2]), stoul(argv[3]));
    // 	printf("%e    %e\n", nu, Lambda);
    // }
    

    // printf("Lambda(%s) = %.4e\n", argv[1], Lambda);
    // printf("lambda(xi)^n = %.4fE%+ld (%.4fE+%ld, %.4fE+%ld)\n\n",
    // 	   bounds[1].first, bounds[1].second,
    // 	   bounds[0].first, bounds[0].second,
    // 	   bounds[2].first, bounds[2].second
    // 	);
    return 0;
}

