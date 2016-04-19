#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <vector>
#include <random>


using namespace std;
using namespace arma;

void print_matrix(const mat& M)
{
    for (unsigned i = 0; i < M.n_rows; i++) {
	for (unsigned j = 0; j < M.n_cols; j++) {
	    printf("% 12.2e", M(i,j));
	}
	printf("\n");
    }
}

mat& gen_rand_matrices(vector<double> &alpha, vector<double> &beta, mat& A)
{
    static normal_distribution<double> ndist;
    static default_random_engine gen(time(NULL));

    z = pow(ndist(gen), 2);
    A.row(0) = rowvec({alpha[0] * z + beta[0], beta[1]});
    A.row(1) = rowvec({z, 0});
    return M;
}

double norm_pow_prod_mat(const vector<mat>& matrices, double alpha)
{
    mat x(4, 4);
    x.eye();
    for (size_t i=0; i < matrices.size(); i++) {
	x *= matrices[i];
    }
    double y = norm(x, 2);
    return pow(y, alpha);
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

int main(int argc, char* argv[])
{
    unsigned n = stoi(argv[1]);
    unsigned p = stoi(argv[2]);
    double alpha = stod(argv[3]);
    // double epsilon = 1.0e-4;
    // double x;
    vector< vector<mat> > matrices;
    
    matrices.resize(p);
    for (unsigned i = 0; i < p; i++) {
	matrices[i].resize(n);
	for (unsigned j = 0; j < n; j++) {
	    mat M(2, 2);
	    vector<double> alpha({0.7});
	    vector<double> beta({0.2, 0.05});
	    matrices[i][j] = gen_rand_matrices(alpha, beta, M);
	}
    }
    double x = mean_norm_pow_prod_mat(matrices, alpha);
    cout<< "Estimated lambda(" << alpha << ") = "
	<< pow(x, 1/double(n))
	<< endl;
    return 0;
}
