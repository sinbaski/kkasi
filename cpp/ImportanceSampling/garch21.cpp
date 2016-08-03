#include <cmath>
#include <random>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "garch21.hpp"

using namespace std;

garch21::F_Set::F_Set(garch21 *markov)
{
    this->markov = markov;
}

// garch21::C_transition_kernel::C_transition_kernel(array<double, 4> params, garch21* markov)
//     :params(params), Hdist(params), markov(markov)
// {
// }

// array<double, 2> garch21::C_transition_kernel::draw(double x0)
// {
//     assert(x0 > 0 && x0 < 1);    
//     return array<double, 2>();
// }

// bool garch21::C_transition_kernel::check_K(array<double, 2> arg, double x0)
// {
//     bool cond = arg[1] > pow(b1, 2);
//     cond = cond && arg[0] > (a1 + pow(b1, 2)/arg[1]) / (1 + a1);
//     cond = cond && arg[0] < 1;
//     return cond;
// }

// double garch21::C_transition_kernel::density(array<double, 2> arg, double x0)
// {
//     assert(x0 > 0 && x0 < 1);
//     assert(check_K(arg, x0));
//     double den = markov->kernel_density(arg, x0);
//     return den;
// }

garch21::garch21(array<double, 4> &params)
    :F(this), nu(this)
{
    a0 = params[0];
    a1 = params[1];
    a2 = params[2];
    b1 = params[3];
    double t = 1 + (1 + a1) / (1 - a1) / 2;
    double b = t * pow(b1, 2);
    // rho = 1/2 * (1 + b1^2/b)
    double rho = 1/2 * (1 + 1/t);
    F.b = b;
    F.rho = rho;
    F.x_interval[0] = (a1 + rho)/(1 + a1);
    F.x_interval[1] = 1;
    F.eta_interval[0] = pow(b1, 2)/rho;
    F.eta_interval[1] = b;

//    double p = C[0] = 0.5;
    double q = C[1] = 1;

    /* Computes delta, the nu measure of F */
    double t1 = b1 * b1;
    double t2 = a2 * a2;
    double t3 = t1 * t2;
    double t4 = q * a1;
    double t7 = t1 * b1 * a2;
    double t9 = q * b1;
    double t10 = t9 * a2;
    double t11 = q * t2;
    double t12 = t1 * a1;
    double t14 = log(t10 - t11 + t12 + t2);
    double t15 = q * t14;
    double t16 = t2 * a2;
    double t17 = t16 * b1;
    double t20 = q * q;
    double t21 = t20 * t14;
    double t28 = log(b * rho * a1 + t10 - t11 + t2);
    double t29 = q * t28;
    double t32 = a1 * rho;
    double t34 = b * t2;
    double t36 = rho * rho;
    double t41 = log(b * a1 + t10 - t11 + t2);
    double t42 = t20 * t41;
    double t43 = t2 * t2;
    double t44 = t43 * t36;
    double t54 = log(-(q * rho * t2 - t9 * rho * a2 - rho * t2 - t12) / rho);
    double t55 = t20 * t54;
    double t57 = t41 * t43;
    double t58 = q * t36;
    double t61 = t54 * t43;
    double t64 = b * b;
    double t66 = a1 * a1;
    double t67 = t66 * t36;
    double t69 = t20 * t28;
    double t79 = -t34 * a1 * t36 + t28 * t64 * t67 - t41 * t64 * t67 + 0.2e1 * t15 * t17 - 0.2e1 * t21 * t17 - 0.2e1 * t29 * t17 + 0.2e1 * t69 * t17 + t21 * t3 + t21 * t43 + 0.2e1 * t29 * t43 + t3 * t32 + t3 * t4 - t69 * t3 + t34 * t32 + t57 * t36 - t7 * t4 + t42 * t44 - t55 * t44 - 0.2e1 * t57 * t58 + 0.2e1 * t61 * t58;
    double t82 = t1 * t1;
    double t91 = b * a2 * b1;
    double t92 = t4 * t36;
    double t94 = t4 * rho;
    double t98 = t17 * t36;
    double t101 = t3 * t36;
    double t107 = t9 * t36;
    double t117 = 0.2e1 * t41 * t16 * t107 - 0.2e1 * t54 * t16 * t107 - t14 * t82 * t66 + t54 * t82 * t66 - t3 * a1 + t42 * t101 - t55 * t101 + t14 * t43 - 0.2e1 * t15 * t43 - t28 * t43 - t3 * t94 + t34 * t92 - t34 * t94 - t61 * t36 - 0.2e1 * t42 * t98 - t69 * t43 + 0.2e1 * t55 * t98 + t7 * t94 - t91 * t92 + t91 * t94;
    nu.delta = -(t79 + t117) / q / t66 / a1 / (a1 + 0.1e1) / t36 / 0.2e1;
}

double garch21::kernel_density(array<double, 2> arg, double x0)
{
    // assert(x0 > 0 && x0 < 1);
    // assert(check_K(arg, x0));

    // double t1 = arg[1] * (1 - arg[0]);
    // double t2 = arg[1] - t1;
    // double t3 = a1 * b1 + a2;
    // double u = t2 - a1 * t1;
    // double v = t1 * t3;
    // v /= t2 * a1 - (t1 + x0 - 1) * pow(a1, 2) + x0 * a2 * b1;

    // double a = (a2 * b1 - pow(b1, 2)) / t3;
    // a += (u - a2 * b1) / t3 / x0;
    // a = gsl_ran_chisq_pdf(a, 1);

    // double b = gsl_ran_chisq_pdf(v, 1);

    // double J = pow(arg[1], 2) / x0;
    // J /= arg[1] * a1 * ((a1 + 1) * arg[0] - a1) +
    // 	a2 * (a2 - x0 * (a2 - a2));

    // return a * b * J;
    return 0;
}

double garch21::compute_tail_index(void)
{
    return 0;
}

garch21::nu_dist::nu_dist(garch21 *markov)
{
    this->markov = markov;
    double a1 = markov->a1;
    double a2 = markov->a2;
    double b1 = markov->b1;

    double p = markov->F.x_interval[0];
    double q = markov->F.x_interval[1];
    double b = markov->F.b;
    double rho = markov->F.rho;

    double t3 = pow(b1, 0.2e1);
    double t12 = (p * a2 * b1 - p * t3 - a2 * b1 + b) / p / (a1 * b1 + a2);
    double chi1 = gsl_ran_chisq_pdf(t12, 1);

    double t10 = pow(a2, 0.2e1);
    double t17 = (a1 * b1 + a2) * (1 - rho) * b / (a1 + 0.1e1) / (b * rho * a1 + p * a2 * b1 - p * t10 + t10);
    double chi2 = gsl_ran_chisq_pdf(t17, 1);

    double t2 = pow(a2, 0.2e1);
    double t6 = pow(b1, 0.2e1);
    t10 = 0.1e1 / q / (q * a2 * b1 - q * t2 + a1 * t6 + t2);
    c1 = chi1 * chi2 * t10;

    double t5 = b * b;
    double t7 = pow(b1, 0.2e1);
    double t8 = t7 * t7;
    t10 = rho * rho;
    c2 = (1 - rho) / (a1 + 1) * (t5 * b - t8 * t7 / t10 / rho) / 0.3e1;
}

double garch21::nu_dist::proposal_density(array<double, 2> arg)
{
    return pow(arg[1], 3) / c2 / 3;
}

array<double, 2> garch21::nu_dist::proposal_draw(void)
{
    array<double, 2> ret;
    uniform_real_distribution<double> unif(0.0, 1.0);
    double y = unif(markov->dev);
    
    double b1 = markov->b1;
    double a1 = markov->a1;
    double rho = markov->F.rho;
    
    double t1 = pow(b1, 0.2e1);
    double t2 = t1 * t1;
    double t4 = rho * rho;
    double t10 = pow((3 * c2 * y * t4 * markov->F.rho) + t2 * t1, 0.1e1 / 0.3e1);
    ret[1] = 0.1e1 / rho * t10;

    uniform_real_distribution<double> U((a1 + rho)/(1 + a1), 1.0);
    ret[0] = U(markov->dev);

    return ret;
}

/*
  double a1 = markov->a1;
  double a2 = markov->a2;
  double b1 = markov->b1;
 */
double garch21::nu_dist::density(array<double, 2> arg)
{
    double a1 = markov->a1;
    double a2 = markov->a2;
    double b1 = markov->b1;
    double q = markov->C[1];

    double x = arg[0];
    double eta = arg[1];
    double t1 = eta * eta;
    double t4 = eta * x;
    double t5 = pow(a1, 0.2e1);
    double t9 = pow(a2, 0.2e1);
    double t15 = t1 / q / (q * a2 * b1 - eta * t5 - q * t9 + t4 * t5 + t4 * a1 + t9);

    return c1 * t15 / delta;
}

array<double, 2> garch21::nu_dist::draw(void)
{
    array<double, 2> proposed;
    uniform_real_distribution<double> unif(0.0, 1.0);
    double d1, d2, u;
    do {
	proposed = proposal_draw();
	d2 = proposal_density(proposed);
	d1 = density(proposed);
	u = unif(markov->dev);
    } while (u > d1 / d2 / c1 / c2);
    return proposed;
}
