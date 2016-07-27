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
    class h_dist
    {
    private:
	/* a0, a1, a2, b1 */
	array<double, 4> params;
	/* Trajectory of the chain after reaching stationarity */
	queue< array<double, 2> > traj;
	random_device dev;
    public:
	array<double, 2> MH_proposal_draw(void);
	double MH_proposal_density(array<double, 2> arg);
	/* Un-normalized density function */
	double density(array<double, 2> arg);
	array<double, 2> draw(void);
	h_dist(array<double, 4> params);
    };

    /* The transition kernel inside C when the chain does NOT re-generate */
    class C_transition_kernel
    {
    private:
	/* a0, a1, a2, b1 */
	array<double, 4> params;
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
	double b;
	double rho;
	garch21 *markov;
	array<double, 2> Fx(void);
	array<double, 2> Feta(void);
	F_Set(garch21 *markov);
	void initialize(double b, double rho);
    };

    class nu_dist
    {
    public:
	
    };
    
    /* a0, a1, a2, b1 */
    array<double, 4> params;
    double delta;
    array<double, 2> C;
    F_Set F;
    C_transition_kernel ctk;
    garch21(array<double, 4> &params);
    double kernel_density(array<double, 2> arg, double x0);
    double compute_tail_index(void);
};

#endif
