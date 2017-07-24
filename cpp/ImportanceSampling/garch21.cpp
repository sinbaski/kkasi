#include "garch21.hpp"
#include <cmath>
#include <ctgmath>
#include <gsl/gsl_math.h>
#include <algorithm>

Garch21::Garch21(const vector<double> &alpha, const vector<double> &beta,
		 double tail_index_sup)
    :GarchPQ(alpha, beta),
     tail_index(find_tail_index(vector<double>({0, tail_index_sup}).data()))
{
    random_device randev;
    normal_distribution ndist(0, 1);
    right_eigenfunction(4000u);
    pool.resize(100000u);
    for (auto i = pool.begin(); i != pool.end(); ++i) {
	*i = ndist(randev);
    }
}

double Garch21::quantile_shifted(double u, double angle) const
{
    mat A;
    vec x(2);
    x(0) = cos(angle);
    x(1) = sin(angle);
    
}

double Garch21::draw_z(bool shift) const
{
    random_device randev;
    normal_distribution ndist(0, 1);
    if (! shift) {
	return ndist(randev);
    }
    
}

double Garch21::r(double angle) const
{
    funval val;
    val[0] = angle;
    val[1] = 0;
    auto j = upper_bound(eigenfunction.cbegin(), eigenfunction.cend(),
			 val, [](const funval &x, const funval &y) {
			     return x[0] < y[0];
			 });
    auto i = eigenfunction.cbegin();
    advance(i, distance(i, j) - 1);
    double slope = (j->at(1) - i->at(1))/(j->at(0) - i->at(0));
    return i->at(1) + slope * (angle - i->at(0));
}

void Garch21::right_eigenfunction(unsigned n)
{
    double ang_max = atan(1/alpha[1]);
    double flag;
    vec angles(n);
    mat P(n, n);
    cx_vec eigenval(n);
    cx_mat eigenmat(n, n);
    double c1 = pow(pow(alpha[2], 2) + pow(beta[0], 2), 0.5);
    double ang1 = atan(beta[0] / alpha[2]);
    
    do {
	vec eigenfun(n);
	vec tangents(n);
	vec cosines(n);

	double d = ang_max/(n + 1);
	P.set_size(n, n);
	angles.set_size(n);
	for (unsigned i = 0; i < n; i++) {
	    angles[i] = d * (i + 1);
	    tangents[i] = tan(angles[i]);
	    cosines[i] = cos(angles[i]);
	}

	for (unsigned i = 0; i < n; i++) {
	    for (unsigned j = 0; j < n; j++) {
		double t1 = sin(angles[i] + ang1);
		double t2 = sqrt(tangents[i] * alpha[2] + beta[0]);
		double t3 = (alpha[2] * tangents[i] + beta[0]);
		t3 *= tangents[j] / (1 - tangents[j] * alpha[1]) / 2;
		t3 = exp(-t3);

		double t5 = sqrt(2 * M_PI * tangents[j]);
		double t6 = pow(cosines[j], tail_index + 2);
		double t7 = pow(1 - tangents[j] * alpha[1],
				tail_index + 1.5);
		P(i, j) = pow(c1 * t1, tail_index) * t2 * t3 / t5 / t6 / t7;
		P(i, j) *= d;
	    }
	}
	flag = det(P - eye(n, n));
    } while (flag > 5.0e-3 && (n *= 2));

    eig_gen(eigenval, eigenmat, P);

    vec values = real(eigenval);
    mat vectors = real(eigenmat);
    auto p =
    	min_element(values.begin(), values.end(),
    		    [](double x, double y) {
    			return fabs(x - 1) < fabs(y - 1);
    		    });
    vec V = vectors.col(p - values.begin());
    eigenfunction.resize(n);
    for (unsigned int i = 0; i < n; ++i) {
	eigenfunction[i][0] = angles(i);
	eigenfunction[i][1] = V(i);
    }
}



int main(int argc, char *argv[])
{
    return 0;
}
