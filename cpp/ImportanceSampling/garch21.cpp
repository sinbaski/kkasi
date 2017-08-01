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
     tail_index(find_tail_index(vector<double>({1, tail_index_sup}).data()))
{
    random_device randev;
    chi_squared_distribution<double> chi2;
    // Estimate the tail index
    pool.resize(10000u);
    for (auto i = pool.begin(); i != pool.end(); ++i) {
	*i = chi2(randev);
    }
    sort(pool.begin(), pool.end());
    // compute right eigenfunctions
    eigenfunctions.resize(400);
    double dx = (tail_index - 0.1)/eigenfunctions.size();
#pragma omp parallel for
    for (unsigned i = 0; i < eigenfunctions.size(); i++) {
	double kappa, lambda_kapp;
	eigenfunctions[i].r_kappa = new vector<funval>(500);
	kappa = i * dx + 0.1;
	lambda_kappa = right_eigenfunction(
	    kappa, eigenfunctions[i]->r_kappa
	    );
	eigenfunctions[i].kappa = kappa;
	eigenfunctions[i].lambda_kappa = lambda_kappa;
    }
}

Garch21::~Garch21(void)
{
    // for (auto p = eigenfunctions.begin();
    // 	 p != )
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

double Garch21::right_eigenfunction
(double index, vector<funval> &eigenfunction) const
{
    double ang_max = atan(1/alpha[1]);
    size_t n = eigenfunction.size();
    vec angles(n);
    mat P(n, n);
    cx_vec eigenval(n);
    cx_mat eigenmat(n, n);
    double c1 = pow(pow(alpha[2], 2) + pow(beta[0], 2), 0.5);
    double ang1 = atan(beta[0] / alpha[2]);
    
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

    eig_gen(eigenval, eigenmat, P);

    vec values = real(eigenval);
    mat vectors = real(eigenmat);
    for (unsigned int i = 0; i < n; ++i) {
	eigenfunction[i][0] = angles(i);
	eigenfunction[i][1] = vectors(i, 0);
    }
    return values(0);
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

double Garch21::b_fun(double theta)
{
    auto p = find_if(
	eigenfunctions.begin(), eigenfunctions.end(),
	[theta](const eigenfunction &fun) {
	    return fun.kappa == theta;
	});

    if (!eigenfunctions.size() || p == eigenfunctions.end()) {
	double omega = alpha[0];
	vector<funval> *r_theta = new vector<funval>(400);
	double lambda_theta =
	    right_eigenfunction(theta, r_theta);
	eigenfunctions.push_back(
	    {theta, lambda_theta, r_theta}
	    );
	assert(lambda_theta < 1);
	return (1 + lambda_theta)/2;
    } else {
	return (1 + p->lambda_kappa)/2;
    }
}

double Garch21::M_fun(double theta)
{
    return 0;
}

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
