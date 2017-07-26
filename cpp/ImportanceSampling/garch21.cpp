#include "garch21.hpp"
#include <cmath>
#include <ctgmath>
#include <gsl/gsl_math.h>
#include <algorithm>
#include <cassert>

Garch21::Garch21(const vector<double> &alpha, const vector<double> &beta,
		 double tail_index_sup)
    :GarchPQ(alpha, beta),
     tail_index(find_tail_index(vector<double>({0, tail_index_sup}).data()))
{
    random_device randev;
    chi_squared_distribution<double> chi2;
    right_eigenfunction(4000u);
    pool.resize(10000u);
    for (auto i = pool.begin(); i != pool.end(); ++i) {
	*i = chi2(randev);
    }
    sort(pool.begin(), pool.end());
}

double Garch21::quantile(double u, double angle) const
{
    assert(u > 0 && u < 1);
    vector<double> edf = pool;
    vec x = {cos(angle), sin(angle)};
    static vector<mat> matrices;
    
    if (matrices.empty()) {
#pragma omp parallel for
	for (unsigned int i = 0; i < pool.size(); ++i) {
	    gen_rand_matrix(pool[i], matrices[i]);
	}
    }
    
#pragma omp parallel for
    for (unsigned int i = 0; i < pool.size(); ++i) {
	vec y = matrices[i] * x;
	double len = norm(y);
	double theta = acos(y[0]/len);
	double inc =
	    pow(len, tail_index) * r(theta) / r(angle) / pool.size();
	edf[i] = inc;
    }
    for (auto i = edf.begin() + 1; i != edf.end(); ++i) {
	*i = *i + *(i - 1);
    }
    auto p = upper_bound(edf.begin(), edf.end(), u);
    if (p == edf.end()) {
	size_t n = pool.size();
	double edf_slope = (edf[n-1] - edf[n-2])/(pool[n-1] - pool[n-2]);
	return pool[n-1] + (u - edf[n-1])/edf_slope;
    } else if (p == edf.begin()) {
	double edf_slope = (edf[1] - edf[0])/(pool[1] - pool[0]);
	return pool[0] + (u - edf[0])/edf_slope;
    } else {
	size_t n = distance(edf.begin(), p);
	double edf_slope = (edf[n] - edf[n-1])/(pool[n] - pool[n-1]);
	return pool[n-1] + (u - edf[n-1])/edf_slope;
    }
}

double Garch21::draw_z2(void) const
{
    random_device randev;
    chi_squared_distribution<double> chi2;
    return chi2(randev);
}

/**
 * Draw from the shifted conditional distribution
 */
double Garch21::draw_z2(double angle) const
{
    random_device randev;
    uniform_real_distribution<double> unif(0, 1);
    return quantile(unif(randev), angle);
}

void Garch21::simulate_sample_path
(double threshold, unsigned int n, vector<vec> path) const
{
    
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
    if (j == eigenfunction.cend()) {
	size_t n = eigenfunction.size();
	double slope =
	    (eigenfunction[n-1][1] - eigenfunction[n-2][1])/
	    (eigenfunction[n-1][0] - eigenfunction[n-2][0]);
	return eigenfunction[n-1][1] +
	    (angle - eigenfunction[n-1][0])*slope;
    } else if (j == eigenfunction.cbegin()) {
	double slope =
	    (eigenfunction[1][1] - eigenfunction[0][1])/
	    (eigenfunction[1][0] - eigenfunction[0][0]);
	return eigenfunction[0][1] + (angle - eigenfunction[0][0])*slope;
    } else {
	auto i = eigenfunction.cbegin();
	advance(i, distance(i, j) - 1);
	double slope = (j->at(1) - i->at(1))/(j->at(0) - i->at(0));
	return i->at(1) + slope * (angle - i->at(0));
    }
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
