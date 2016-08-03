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
    /* Draw a random vector from the proposal distribution */
    // class h_dist
    // {
    // private:
    // 	/* a0, a1, a2, b1 */
    // 	array<double, 4> params;
    // 	/* Trajectory of the chain after reaching stationarity */
    // 	queue< array<double, 2> > traj;
    // 	random_device dev;
    // public:
    // 	array<double, 2> MH_proposal_draw(void);
    // 	double MH_proposal_density(array<double, 2> arg);
    // 	/* Un-normalized density function */
    // 	double density(array<double, 2> arg);
    // 	array<double, 2> draw(void);
    // 	h_dist(array<double, 4> params);
    // };

    /* The transition kernel inside C when the chain does NOT re-generate */
    class C_transition_kernel
    {
    private:
	/* a0, a1, a2, b1 */
	double a0, a1, a2, b1;
	h_dist Hdist;
	garch21 *markov;
    public:
	array<double, 2> draw(double x0);
	double density(array<double, 2> arg, double x0);
	bool check_cond(array<double, 2> arg, double x0);
	C_transition_kernel(array<double, 4> params, garch21* markov);
    };

    class F_Set
    {
    public:
	garch21 *markov;
	double b;
	double rho;
	array<double, 2> x_interval;
	array<double, 2> eta_interval;
	F_Set(garch21 *markov);
    };

    class nu_dist
    {
    public:
	/* The nu measure of the F set */
	double delta;
	double c1, c2;
	garch21 *markov;
	double draw(void);
	double density(array<double, 2> arg);
	nu_dist(garch21 *markov);
    };
    
    /* a0, a1, a2, b1 */
    double a0, a1, a2, b1;
    array<double, 2> C;
    F_Set F;
    nu_dist nu;
    C_transition_kernel ctk;
    garch21(array<double, 4> &params);
    double kernel_density(array<double, 2> arg, double x0);
    double compute_tail_index(void);
};

#endif
