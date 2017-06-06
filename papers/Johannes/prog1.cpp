#include <algorithm>
#include <armadillo>
#include <cmath>
#include <functional>
#include <iterator>
#include <assert.h>
#include <iostream>
#include <random>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;
using namespace arma;

random_device randev;

class Lindep
{
protected:
public:
    mat h; // coefficients' matrix
    double alpha = 1.5; // the tail index
    double beta = 0.5;
    mat *X;
    
    Lindep(mat &h);
    mat &set_X(mat *X);
    virtual mat H(unsigned s, unsigned n_rows, unsigned n_cols);
    virtual mat M(unsigned s, unsigned n_rows, unsigned n_cols);
    virtual mat K(unsigned s1, unsigned s2, unsigned n_rows, unsigned n_cols);
    virtual mat C(unsigned s);
    virtual mat& sim_X(mat &X);
};

Lindep::Lindep(mat &h)
    :h(h), X(NULL)
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

mat Lindep::C(unsigned s)
{
    return X->submat(0, 0, X->n_rows - 1, X->n_cols - s - 1) * X->submat(0, s, X->n_rows - 1, X->n_cols - 1).t();
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
    // cout << "Matrix Z:" << endl;
    // Z.print();

    // cout << "Matrix h:" << endl;
    // h.print();

    mat g = h.t();
    for (unsigned i = 0; i < h.n_rows; i++) {
	double *ptr = g.colptr(i);
	reverse(ptr, ptr + h.n_cols);
    }
    g = g.t();

    // cout << "Matrix g:" << endl;
    // g.print();
    for (unsigned i = 0; i < p; i++) {
#pragma omp parallel for
	for (unsigned j = 0; j < n; j++) {
	    X(i, j) = sum(sum(Z.submat(i, j, i + h.n_rows - 1, j + h.n_cols - 1) % g));
	}
    }
    return X;
}

mat &Lindep::set_X(mat *X)
{
    this->X = X;
    return *X;
}

void case_iid(void)
{
    mat h({1});
    
#pragma omp parallel for
    for (unsigned n = 1000; n < 40000; n+=1000) {
	Lindep lindep(h);
	unsigned p = ceil(pow(n, lindep.beta));
	mat X(p, n);
	lindep.sim_X(X);
	lindep.set_X(&X);

	// cout << "matrix X" <<endl;
	// X.submat(0, 0, 9, 12).print();

	mat C0 = lindep.C(0);
	mat T = C0;

	unsigned s2 = 5;
	for (unsigned s = 1; s <= s2; s++) {
	    T += lindep.C(s);
	}
	double a_np = gsl_cdf_tdist_Pinv(
	    1 - 1/(double)(X.n_cols * X.n_rows),
	    lindep.alpha);
	double y = norm(T - C0)/a_np/a_np;
	// cout << "matrix T" << endl;
	// T.print();

	// cout << "matrix C0" << endl;
	// C0.print();

#pragma omp critical
	printf("%5u\t%e\n", n, y);
    }
}

void case_email(void)
{
    vec c = {1, 0.5, 0.25};
    vec d = {1, 0.3, 0.09};
    mat h = kron(d, c.t());

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
}

int main(int argc, char *argv[])
{
    case_iid();
    return 0;
}
