#ifndef XXIE_GARCH21_H
#define XXIE_GARCH21_H
#include <array>
#include <queue>
#include <random>

using namespace std;

class garch21
{
protected:
    /* a0, a1, a2, b1 */
    array<double, 4> params;
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
    h_dist Hdist;
    garch21(array<double, 4> &params);
    double compute_tail_index(void);
};

#endif
