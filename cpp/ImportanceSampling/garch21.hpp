#ifndef XXIE_GARCH21_H
#define XXIE_GARCH21_H
#include <array>

using namespace std;

class garch21
{
protected:
    /* a0, a1, a2, b1 */
    double params[4];
public:
    /* Draw a random vector from the proposal distribution */
    class h_dist
    {
    public:
	array<double, 2> MH_proposal_draw(random_device &dev, const array<double, 2> &init);
	double MH_proposal_density(array<double, 2> arg);
	array<double, 2> MH_target_draw(random_device &dev, const array<double, 2> &init);
	double MH_target_density(array<double, 2> arg);
    };
    
};

#endif
