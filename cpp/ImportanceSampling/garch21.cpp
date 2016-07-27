#include <cmath>
#include <random>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "garch21.hpp"

using namespace std;

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


garch21::F_Set::F_Set(garch21 *markov)
    :markov(markov)
{
}

void garch21::F_Set::initialize(double b, double rho)
{
    this->b = b;
    this->rho = rho;
}

array<double, 2> garch21::F_Set::Fx(void)
{
    array<double, 2> ntv;
    ntv[0] = (markov->params[1] + rho)/(1 + markov->params[1]);
    ntv[1] = 1;
    return ntv;
}

array<double, 2> garch21::F_Set::Feta(void)
{
    array<double, 2> ntv;
    ntv[0] = pow(markov->params[3], 2)/rho;
    ntv[1] = b;
    return ntv;
}

garch21::C_transition_kernel::C_transition_kernel(array<double, 4> params, garch21* markov)
    :params(params), Hdist(params), markov(markov)
{
}

array<double, 2> garch21::C_transition_kernel::draw(double x0)
{
    assert(x0 > 0 && x0 < 1);    
    return array<double, 2>();
}

bool garch21::C_transition_kernel::check_K(array<double, 2> arg, double x0)
{
    bool cond = arg[1] > pow(params[3], 2);
    cond = cond && arg[0] > (params[1] + pow(params[3], 2)/arg[1]) / (1 + params[1]);
    cond = cond && arg[0] < 1;
    return cond;
}

double garch21::C_transition_kernel::density(array<double, 2> arg, double x0)
{
    assert(x0 > 0 && x0 < 1);
    assert(check_K(arg, x0));
    double den = markov->kernel_density(arg, x0);
    return den;
}

garch21::garch21(array<double, 4> &params)
    :params(params), ctk(params, this), F(this)
{
    double t = 1 + (1 + params[1]) / (1 - params[1]) / 2;
    double b = t * pow(params[3], 2);
    double rho = 1/2 * (1 + 1/t);
    F.initialize(b, rho);
    C[0] = 0.5;
    C[1] = 1;

    delta = -(b * params[2] * params[3] * params[1] * C[1] * rho + b * pow(params[2], 0.2e1) * params[1] * rho + b * pow(params[2], 0.2e1) * params[1] * C[1] * rho * rho + C[1] * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) - log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[1], 0.2e1) * pow(params[3], 0.4e1) - 0.2e1 * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) - C[1] * C[1] * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) * rho * rho * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) - 0.2e1 * pow(params[2], 0.3e1) * params[3] * C[1] * rho * rho * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) - 0.2e1 * C[1] * C[1] * pow(params[2], 0.3e1) * params[3] * rho * rho * log(rho) + C[1] * C[1] * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) * rho * rho * log(rho) + 0.2e1 * pow(params[2], 0.3e1) * params[3] * C[1] * rho * rho * log(rho) - log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) + C[1] * rho * params[1] * params[2] * pow(params[3], 0.3e1) - C[1] * rho * params[1] * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) - C[1] * C[1] * log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) + log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) * rho * rho + 0.2e1 * C[1] * log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) - C[1] * C[1] * pow(params[2], 0.4e1) * rho * rho * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) + 0.2e1 * pow(params[2], 0.4e1) * C[1] * rho * rho * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) + C[1] * C[1] * pow(params[2], 0.4e1) * rho * rho * log(rho) - 0.2e1 * pow(params[2], 0.4e1) * C[1] * rho * rho * log(rho) + log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) - pow(params[2], 0.4e1) * rho * rho * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) - pow(params[1], 0.2e1) * pow(params[3], 0.4e1) * log(rho) + pow(params[1], 0.2e1) * pow(params[3], 0.4e1) * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) + pow(params[2], 0.4e1) * rho * rho * log(rho) + rho * params[1] * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) - C[1] * params[1] * params[2] * pow(params[3], 0.3e1) + C[1] * params[1] * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) - 0.2e1 * C[1] * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.3e1) * params[3] + C[1] * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) + 0.2e1 * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[3], 0.2e1) * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.3e1) * params[3] + 0.2e1 * C[1] * C[1] * pow(params[2], 0.3e1) * params[3] * rho * rho * log(-C[1] * rho * pow(params[2], 0.2e1) + C[1] * rho * params[2] * params[3] + rho * pow(params[2], 0.2e1) + pow(params[3], 0.2e1) * params[1]) - params[1] * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) + C[1] * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) * rho * rho - 0.2e1 * C[1] * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.3e1) * params[3] * rho * rho + C[1] * C[1] * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) * rho * rho + 0.2e1 * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.3e1) * params[3] * C[1] * rho * rho - b * pow(params[2], 0.2e1) * params[1] * rho * rho - b * params[2] * params[3] * params[1] * C[1] * rho * rho - b * pow(params[2], 0.2e1) * params[1] * C[1] * rho - 0.2e1 * log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * pow(params[2], 0.4e1) * C[1] * rho * rho + log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * b * b * pow(params[1], 0.2e1) * rho * rho + 0.2e1 * C[1] * C[1] * log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * pow(params[2], 0.3e1) * params[3] - C[1] * C[1] * log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * pow(params[2], 0.2e1) * pow(params[3], 0.2e1) - log(-C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + b * params[1] + pow(params[2], 0.2e1)) * b * b * pow(params[1], 0.2e1) * rho * rho - 0.2e1 * C[1] * log(b * rho * params[1] - C[1] * pow(params[2], 0.2e1) + C[1] * params[2] * params[3] + pow(params[2], 0.2e1)) * pow(params[2], 0.3e1) * params[3]) * pow(params[1], -0.3e1) / (0.1e1 + params[1]) / C[1] * pow(rho, -0.2e1) / 0.2e1;

    double t1 = (params[2] * params[3] - pow(params[3], 2))/(params[1] * params[3] + params[2]);
    t1 += (-params[2] * params[3] + F.b)/ C[0] / (params[1] * params[3] + params[2]);
    t1 = gsl_ran_chisq_pdf(t1, 1);
}

double garch21::kernel_density(array<double, 2> arg, double x0)
{
    assert(x0 > 0 && x0 < 1);
    assert(check_K(arg, x0));

    double t1 = arg[1] * (1 - arg[0]);
    double t2 = arg[1] - t1;
    double t3 = params[1] * params[3] + params[2];
    double u = t2 - params[1] * t1;
    double v = t1 * t3;
    v /= t2 * params[1] - (t1 + x0 - 1) * pow(params[1], 2) + x0 * params[2] * params[3];

    double a = (params[2] * params[3] - pow(params[3], 2)) / t3;
    a += (u - params[2] * params[3]) / t3 / x0;
    a = gsl_ran_chisq_pdf(a, 1);

    double b = gsl_ran_chisq_pdf(v, 1);

    double J = pow(arg[1], 2) / x0;
    J /= arg[1] * params[1] * ((params[1] + 1) * arg[0] - params[1]) +
	params[2] * (params[2] - x0 * (params[2] - params[2]));

    return a * b * J;
}

double garch21::compute_tail_index(void)
{
    return 0;
}
