#include <cmath>
#include <random>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "garch21.hpp"

using namespace std;

garch21::garch21(array<double, 4> &params)
    :params(params), Hdist(params)
{
    copy(params.begin(), params.end(), this->params.begin());
}

/* Prepare a chain trajectory that has reached stationarity.
   Assume the chain has reached stationarity after 1000 steps.
   
 */
garch21::h_dist::h_dist(array<double, 4> params)
    :params(params)
{
    uniform_real_distribution<double> unif;
    int i = 1, n = 1000;

    /* Start at a random position in the state space */
    array<double, 2> T = MH_proposal_draw();
    double omg = density(T) / MH_proposal_density(T);
    traj.push(T);
    do {
	array<double, 2> newstate = MH_proposal_draw();
	double ratio = density(newstate) / MH_proposal_density(newstate);
	if (unif(dev) <= ratio / omg) {
	    traj.push(newstate);
	    omg = ratio;
	    i++;
	}
    } while (i < n);
}

array<double, 2> garch21::h_dist::MH_proposal_draw(void)
{
    array<double, 2> final;
    gamma_distribution<double> Gamma(1/2, 4);
    final[1] = Gamma(dev);

    uniform_real_distribution<double> Unif((params[1] + pow(params[3], 2)/final[1])/(1 + params[1]), 1);
    final[0] = Unif(dev);

    return final;
}

/**
   The proposal density function q(x, eta)
   q(x, eta) = gamma_pdf(eta, 1/2, 4) * unif_pdf(x, (a1 + b1^2 / eta) / (1 + a1), 1)
 */
double garch21::h_dist::MH_proposal_density(array<double, 2> arg)
{
    double zeta = (arg[1]/params[3] - params[3])/(1 + params[1]);
    double den = gsl_ran_gamma_pdf(zeta, 1/2, 4);
    den /= 1 - (params[1] + pow(params[3], 2)/arg[1])/(1 + params[1]);
    return den;
}

array<double, 2> garch21::h_dist::draw(void)
{
    static double omg = density(traj.back()) / MH_proposal_density(traj.back());
    static uniform_real_distribution<double> unif;
    array<double, 2> T = MH_proposal_draw();
    double ratio;
    while (unif(dev) > (ratio = density(T) / MH_proposal_density(T))/omg) {
	T = MH_proposal_draw();
    }
    omg = ratio;
    traj.push(T);
    traj.pop();
    return T;
}

double garch21::h_dist::density(array<double, 2> arg)
{
    double chi_arg = arg[1] * (1 - arg[0]) * (params[1] * params[3] + params[2]);
    double a = params[1] * (arg[0] * params[1] + arg[0] - params[1]) * arg[1];
    chi_arg /= a + params[2] * params[3];
    double J = pow(arg[1], 2) / (a + pow(params[2], 2));
    return gsl_ran_chisq_pdf(chi_arg, 1) * J;
}

double garch21::compute_tail_index(void)
{
    array<double, 2> sample = Hdist.draw();
    return sample[1];
}
