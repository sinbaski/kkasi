#include <algorithm>
#include <armadillo>
#include <cmath>
#include <functional>
#include <iterator>
#include <assert.h>
#include <iostream>
#include <random>
#include <vector>

using namespace std;
using namespace arma;

random_device randev;

class Lindep
{
protected:
public:
    mat h; // coefficients' matrix
    double alpha = 1.5; // the tail index
    Lindep(mat &h);
    virtual mat H(unsigned s, unsigned n_rows, unsigned n_cols);
    virtual mat M(unsigned s, unsigned n_rows, unsigned n_cols);
    virtual mat K(unsigned s1, unsigned s2, unsigned n_rows, unsigned n_cols);
    virtual mat C(mat &X, unsigned s);
    virtual mat& sim_X(mat &X);
};

Lindep::Lindep(mat &h)
    :h(h)
{
}

mat Lindep::H(unsigned s, unsigned n_rows, unsigned n_cols)
{
    mat T = zeros(n_cols, n_rows);
    if (s >= h.n_cols) return T.t();
    mat g = h.t();
    for (unsigned i = 0; i < h.n_rows; i++) {
	double *p = g.colptr(i);
	double *q = T.colptr(i);
    	copy(p + s, h.n_cols - s < T.n_cols ? p + h.n_cols : p + s + T.n_cols, q);
    }
    return T.t();
}

mat Lindep::M(unsigned s, unsigned n_rows, unsigned n_cols)
{
    mat H0 = H(0, n_rows, n_cols);
    mat Hs = H(s, n_rows, n_cols);
    return H0 * Hs.t();
}

mat Lindep::K(unsigned s1, unsigned s2, unsigned n_rows, unsigned n_cols)
{
    mat T = zeros(n_rows, n_cols);
    for (unsigned s = s1; s <= s2; s++) {
	mat m = M(s, n_rows, n_cols);
	T += m * m.t();
    }
    return T;
}

mat Lindep::C(mat &X, unsigned s)
{
    
}

mat& Lindep::sim_X(mat &X)
{
    unsigned p = X.n_rows;
    unsigned n = X.n_cols;
    student_t_distribution<double> stud(alpha);
    mat Z(p + h.n_rows, n + h.n_cols);

    for (unsigned i = 0; i < Z.n_rows; i++) {
#pragma omp parallel for
	for (unsigned j = 0; j < Z.n_cols; j++) {
	    Z(i, j) = stud(randev);
	}
    }
    mat g = h.t();
    for (unsigned i = 0; i < h.n_rows; i++) {
	double *ptr = g.colptr(i);
	reverse(ptr, ptr + h.n_cols);
    }
    g = g.t();
    for (unsigned i = 0; i < p; i++) {
#pragma omp parallel for
	for (unsigned j = 0; j < n; j++) {
	    X(i, j) = sum(sum(Z.submat(i, j, i + h.n_rows - 1, j + h.n_cols - 1) % g));
	}
    }
    return X;
}

int main(int argc, char *argv[])
{
    vec c = {1, 0.5, 0.25};
    vec d = {1, 0.3, 0.09};
    // unsigned LEN = c.n_elem;
    
    // vec c_bar(LEN);
    // for (unsigned s = 0; s < LEN; s++) {
    // 	vec temp = zeros(LEN);
    // 	vec temp2(LEN);
    // 	copy(c.begin() + s, c.end(), temp.begin());
    // 	transform(
    // 	    c.begin(), c.end(), temp.begin(), temp2.begin(),
    // 	    std::multiplies<double>()
    // 	    );
    // 	c_bar[s] = accumulate(temp2.begin(), temp2.end(), double(0));
    // }
    mat h = kron(d, c.t());
    // mat D = d * d.t();

    Lindep lindep(h);
    vec val;
    mat modal;

    mat K00 = lindep.K(0, 0, 3, 3);
    cout << "matrix K00" << endl;
    K00.print();

    eig_sym(val, modal, K00);
    cout << "eigenvalues of  K00" << endl;
    val.print();

    for (unsigned i = 0; i < modal.n_cols; i++) {
	double sign = modal(0, i) >= 0 ? 1 : -1;
	transform(modal.begin_col(i), modal.end_col(i), modal.begin_col(i),
		  [=](double x) {
		      return x * sign;
		  });
    }
    cout << "eigenvectors of  K00" << endl;
    modal.print();

    mat X(200, 2000);
    lindep.sim_X(X);
    cout << "Matrix X:" << endl;
    X.submat(0, 0, 10, 10).print();

    return 0;
}
