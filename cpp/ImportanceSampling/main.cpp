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

mat& gen_rand_matrices(double a, double b, mat& M)
{
    static normal_distribution<double> ndist;
    static default_random_engine gen(time(NULL));

    M.row(0) = rowvec({b, 0, a, 0});
    M.row(1) = rowvec({0, b, 0, a});
    M.row(2) = M.row(0) * ndist(gen);
    M.row(3) = M.row(1) * ndist(gen);
    return M;
}

double norm_pow_prod_mat(const vector<mat>& matrices, double alpha)
{
    mat x(4, 4);
    x.eye();
    for (size_t i=0; i < matrices.size(); i++) {
	x *= matrices[i];
    }
    return pow(norm(x, 2), alpha);
}

double mean_norm_pow_prod_mat(const vector< vector<mat> >& matrices, double alpha)
{
    size_t p = matrices.size();
    double x;
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
	    mat M(4, 4);
	    double a = 0.5;
	    double b = 0.4;
	    matrices[i][j] = gen_rand_matrices(a, b, M);
	}
    }
    cout<< "Estimated Lambda(" << alpha << ") = " << log(mean_norm_pow_prod_mat(matrices, alpha))/n << endl;
    return 0;
}
