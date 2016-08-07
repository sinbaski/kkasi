#ifndef XXIE_GARCH21_H
#define XXIE_GARCH21_H
#include <array>
#include <queue>
#include <random>

using namespace std;

class garch21
{
protected:
    double tail_index;
    random_device dev;
public:
    class F_Set
    {
    public:
	garch21 *markov;
	double b;
	double rho;
	array<double, 2> x_interval;
	array<double, 2> eta_interval;
	F_Set(garch21 *markov);
	bool includes(const array<double, 2> &arg) const;
    };

    class nu_dist
    {
    public:
	/* The nu measure of the F set */
	double delta;
	double c1, c2, c3;
	garch21 *markov;
	array<double, 2> draw(void) ;
	double density(const array<double, 2> &arg) const;
	array<double, 2> proposal_draw(void);
	double proposal_density(array<double, 2> arg);
	nu_dist(garch21 *markov);
    };

    double a0, a1, a2, b1;
    array<double, 2> C;
    F_Set F;
    nu_dist nu;
    garch21(array<double, 4> &params);
    inline bool C_includes(array<double, 2> arg) ;
    double kernel_density(const array<double, 2> &arg, double x0) const;
    array<double, 2> forward(double x0, bool orig = true) ;
    array<double, 2> simulate_path(void) ;
    double compute_tail_index(void) ;
};

#endif
