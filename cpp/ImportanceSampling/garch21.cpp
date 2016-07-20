#include <cmath>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "garch21.hpp"


array<double, 2> garch21::h_dist::MH_proposal_draw(random_device &dev, const array<double, 2> &init)
{
    array<double, 2> final;
    double a = pow(params[3], 2);
    gamma_distribution<double> Gamma(1/2, 4);
    final[1] = Gamma(dev);

    uniform_real_distribution<double> Unif((params[1] + pow(params[3], 2)/final[1])/(1 + params[1]), 1);
    final[0] = Unif(dev);

     return final;
}

double garch21::h_dist::MH_proposal_density(array<double, 2> arg)
{
    double zeta = (arg[1]/params[3] - params[3])/(1 + params[1]);
    double den = gsl_ran_gamma_pdf(zeta, 1/2, 4);
    den /= 1 - (params[1] + pow(params[3], 2)/arg[1])/(1 + params[1]);
    return den;
}

array<double, 2> garch21::h_dist::MH_target_draw(random_device &dev, const array<double, 2> &init)
{
    
}

double garch21::h_dist::MH_target_density(array<double, 2> arg)
{
}
