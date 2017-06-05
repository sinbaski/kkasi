#include <algorithm>
#include <armadillo>
#include <cmath>
#include <functional>
#include <iterator>
#include <assert.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace arma;

class Lindep
{
protected:
    vector<vector<double>> h;
public:
    Lindep(mat &h);
    virtual mat H(unsigned s, unsigned n_rows, unsigned n_cols);
    virtual mat M(unsigned s, unsigned n_rows, unsigned n_cols);
    virtual mat K(unsigned s1, unsigned s2, unsigned n_rows, unsigned n_cols);
};

Lindep::Lindep(mat &h)
{
    // this->h.resize(h.n_rows);
    // for (unsigned i = 0; i < h.n_rows; i++) {
    // 	this->h[i].resize(h.n_cols);
    // 	copy(h.begin_row(i), h.end_row(i), this->h[i].begin());
    // }
}

mat Lindep::H(unsigned s, unsigned n_rows, unsigned n_cols)
{
    mat T = zeros(n_rows, n_cols);
    if (s >= h[0].size()) return T;
    // for (unsigned i = 0; i < h.size(); i++) {
    // 	copy(
    // 	    next(h[i].begin(), s),
    // 	    h[i].size() - s < T.n_cols ? h[i].end() : next(h[i].begin(), s + T.n_cols),
    // 	    T.begin_row(i)
    // 	    );
    // }
    return T;
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

// mat K(unsigned s1, unsigned s2, mat &D, vec &c_bar)
// {
//     vec squared(s2 - s1 + 1);
//     transform(
// 	c_bar.begin() + s1, c_bar.begin() + s2 + 1,
// 	squared.begin(),
// 	[&](double x) {
// 	    return x*x;
// 	});
//     double a = accumulate(squared.begin(), squared.end(), double(0));
//     return a * D;
// }

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
    // vec val;
    // mat modal;

    // mat K00 = lindep.K(0, 0, 3, 3);
    // cout << "matrix K00" << endl;
    // K00.print();

    // eig_sym(val, modal, K00);
    // cout << "eigenvalues of  K00" << endl;
    // val.print();

    // cout << "eigenvectors of  K00" << endl;
    // modal.print();


    return 0;
}
