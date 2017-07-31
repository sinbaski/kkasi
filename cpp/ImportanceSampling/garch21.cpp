#include "garch21.hpp"
#include <cmath>
#include <ctgmath>
#include <gsl/gsl_math.h>
#include <algorithm>
#include <cassert>
#include "XMatrix.hpp"

double interpolate_fun
(double arg, const vector<funval> &fun)
{
    funval val;
    val[0] = arg;
    val[1] = 0;
    auto j = upper_bound(fun.cbegin(), fun.cend(),
			 val, [](const funval &x, const funval &y) {
			     return x[0] < y[0];
			 });
    if (j == fun.cend()) {
	size_t n = fun.size();
	double slope =
	    (fun[n-1][1] - fun[n-2][1])/
	    (fun[n-1][0] - fun[n-2][0]);
	return fun[n-1][1] +
	    (arg - fun[n-1][0])*slope;
    } else if (j == fun.cbegin()) {
	double slope =
	    (fun[1][1] - fun[0][1])/
	    (fun[1][0] - fun[0][0]);
	return fun[0][1] + (arg - fun[0][0])*slope;
    } else {
	size_t n = distance(fun.begin(), j);
	double slope = (fun[n][1] - fun[n-1][1])/(fun[n][0] - fun[n-1][0]);
	return fun[n-1][1] + (arg - fun[n-1][0]) * slope;
    }
}

Garch21::Garch21(const vector<double> &alpha, const vector<double> &beta,
		 double tail_index_sup)
    :GarchPQ(alpha, beta),
     tail_index(find_tail_index(vector<double>({0, tail_index_sup}).data()))
{
    random_device randev;
    chi_squared_distribution<double> chi2;
    pool.resize(10000u);
    for (auto i = pool.begin(); i != pool.end(); ++i) {
	*i = chi2(randev);
    }
    sort(pool.begin(), pool.end());
}

double Garch21::quantile(double u, double angle) const
{
    assert(u > 0 && u < 1);
    vector<funval> equantile(pool.size());
    vec x = {cos(angle), sin(angle)};
    static vector<mat> matrices;
    
    if (matrices.empty()) {
#pragma omp parallel for
	for (unsigned int i = 0; i < pool.size(); ++i) {
	    gen_rand_matrix(pool[i], matrices[i]);
	}
    }

    vector<funval> r(pool.size());
    right_eigenfunction(tail_index, r);
#pragma omp parallel for
    for (unsigned int i = 0; i < pool.size(); ++i) {
	vec y = matrices[i] * x;
	double len = norm(y);
	double theta = acos(y[0]/len);
	double inc =
	    pow(len, tail_index) *
	    interpolate_fun(theta, r) /
	    interpolate_fun(angle, r) /
	    pool.size();
	equantile[i][1] = pool[i];
	equantile[i][0] = inc;
    }
    for (auto i = equantile.begin() + 1; i != equantile.end(); ++i) {
	i->at(0) = i->at(0) + (i - 1)->at(0);
    }
    return interpolate_fun(u, equantile);
    // auto p = upper_bound(edf.begin(), edf.end(), u);
    // if (p == edf.end()) {
    // 	size_t n = pool.size();
    // 	double edf_slope = (edf[n-1] - edf[n-2])/(pool[n-1] - pool[n-2]);
    // 	return pool[n-1] + (u - edf[n-1])/edf_slope;
    // } else if (p == edf.begin()) {
    // 	double edf_slope = (edf[1] - edf[0])/(pool[1] - pool[0]);
    // 	return pool[0] + (u - edf[0])/edf_slope;
    // } else {
    // 	size_t n = distance(edf.begin(), p);
    // 	double edf_slope = (edf[n] - edf[n-1])/(pool[n] - pool[n-1]);
    // 	return pool[n-1] + (u - edf[n-1])/edf_slope;
    // }
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
(double threshold, unsigned int n, vector<vec> &path) const
{
}

void Garch21::right_eigenfunction
(double index, vector<funval> &eigenfunction) const
{
    double ang_max = atan(1/alpha[1]);
    double flag;
    size_t n = eigenfunction.size();
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
		double t6 = pow(cosines[j], index + 2);
		double t7 = pow(1 - tangents[j] * alpha[1],
				index + 1.5);
		P(i, j) = pow(c1 * t1, index) * t2 * t3 / t5 / t6 / t7;
		P(i, j) *= d;
	    }
	}
	flag = det(P - eye(n, n));
    } while (flag > 5.0e-3 && (n *= 2));

    eig_gen(eigenval, eigenmat, P);

    cout << "eigenvalues: " << endl;
    for (auto i = eigenval.begin(); i != eigenval.end(); ++i) {
	cout << *i << endl;
    }

    vec values = real(eigenval);
    mat vectors = real(eigenmat);
    auto p =
    	min_element(values.begin(), values.end(),
    		    [](double x, double y) {
    			return fabs(x - 1) < fabs(y - 1);
    		    });
    vec V = vectors.col(p - values.begin());
    for (unsigned int i = 0; i < n; ++i) {
	eigenfunction[i][0] = angles(i);
	eigenfunction[i][1] = V(i);
    }
}

double Garch21::Lyapunov(size_t n, size_t iterations)
{
    vector<XMatrix> matrices(n * iterations);
    random_device randev;
    chi_squared_distribution<double> chi2;
    for (size_t i = 0; i < matrices.size(); i++) {
	double z2 = chi2(randev);
	mat M(2,2);
	gen_rand_matrix(z2, M);
	matrices[i] = XMatrix(M);
    }
    vector<double> samples(iterations);
    cout << "Samples of top Lyapunov exponent:" << endl;
    for (size_t i = 0; i < iterations; i++) {
	XMatrix M = matrices[i * n];
	for (size_t j = 1; j < n; j++) {
	    M *= matrices[i * n + j];
	}
	long power;
	mat Pi_n = M.comptify(&power);
	samples[i] = log(norm(Pi_n)) + log(power);
	cout<< samples[i] << endl;
    }
    double gamma = accumulate(samples.begin(), samples.end(), 0.0,
			      [iterations](double s, double a) {
				  return s + a/iterations;
			      });
    return gamma;
}

// double Garch21::G_fun(double phi, double theta, double m)
// {
//     size_t n = 10000u;
//     vector<funval> r_phi(n);
//     vector<funval> r_theta(n);
//     right_eigenfunction(phi, r_phi);
//     right_eigenfunction(theta, r_theta);
//     double r_phi_sup =
// 	max_element(r_phi.begin(), r_phi.end(),
// 		    [](const funval &a, const funval &b) {
// 			return a[1] < b[1];
// 		    })->at(1);
//     double r_theta_sup =
//     	max_element(r_theta.begin(), r_theta.end(),
//     		    [](const funval &a, const funval &b) {
//     			return a[1] < b[1];
//     		    })->at(1);
//     return r_phi_sup + r_theta_sup;
// }

int main(int argc, char *argv[])
{
    // DJIA
    vector<double> alpha({3.374294e-06, 0.061577928, 0.12795424});
    vector<double> beta({0.6610499});

    Garch21 model(alpha, beta);
    double theta = 0.5;
    vector<funval> r_theta(100);
    model.right_eigenfunction(theta, r_theta);
    return 0;
}
