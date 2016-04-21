#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <random>

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

int main(int argc, char* argv[])
{
    unsigned n = stoi(argv[1]);
    unsigned p = stoi(argv[2]);
    double xi = stod(argv[3]);
    vector<double> alpha({0.6, 0.1, 0.1});
    vector<double> beta({0.1, 0.05});

    for (int i=0; i < stoi(argv[4]); i++) {
    	double x = Lambda_func(xi, n, p, alpha, beta);
    	cout<< x << endl;
    }

    // for (int i = 0; i < stoi(argv[1]); i++) {
    // 	cout << test() << endl;
    // }
    
    return 0;
}
