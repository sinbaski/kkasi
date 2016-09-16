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
		      unsigned long K,
		      double &sd)
{

    uniform_real_distribution<double> unif;
    chi_squared_distribution<double> chi2;
    // vector<double> beta(N, 0);
    vector<vec> E(K);
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
	vector<double> alpha(K);
	vector<mat> A(K);
	double beta = 0;

#pragma omp parallel for schedule(dynamic) shared(gen, chi2, beta)
	for (unsigned k = 0; k < K; k++) {
	    double z2;

#pragma omp critical
	    z2 = chi2(gen);
	    gen_rand_matrix(a, b, z2, A[k]);
	    alpha[k] = pow(norm(A[k] * E[k], "inf"), theta);
#pragma omp atomic
	    beta += alpha[k];
	}

	vector<double> Q(K);
	partial_sum(alpha.begin(), alpha.end(), Q.begin());
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
    }
    sd = sqrt(sd - pow(Lambda, 2));
    return Lambda;
}


int main(int argc, char*argv[])
{
    // vector<double> alpha({1.0e-7, 0.6, 0.001});
    // vector<double> beta({0.005});
    vector<double> alpha({1.0e-7, 0.11, 1.0e-16});
    vector<double> beta({0.88});

    double Lambda;
    
    cout << "alpha[0]= "  << alpha[0] << ", alpha[1]=" <<
    	alpha[1] << ", alpha[2]=" << alpha[2] <<
    	", beta[1]=" << beta[0] << endl;
    cout << "N = " << argv[1] << endl;
    cout << "K = " << argv[2] << endl;

    for (double nu = stod(argv[3]); nu < stod(argv[4]); nu += 0.1) {
	double sd;
	Lambda = estimateLambda(
	    alpha, beta, nu, stoul(argv[1]), stoul(argv[2]), sd);
	printf("%.4f    % .4f    %.4f    %.4f\n", nu, Lambda, sd, sd/abs(Lambda));
    }
    
    return 0;
}

