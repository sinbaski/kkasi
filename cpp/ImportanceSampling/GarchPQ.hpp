#ifndef garchPQ_HPP
#define garchPQ_HPP

#include <random>
#include <vector>
#include <cmath>
#include <ctgmath>
#include <armadillo>

using namespace std;
using namespace arma;

class GarchPQ
{
public:
    const vector<double> alpha;
    const vector<double> beta;

    // double M;
    // double xi;

    GarchPQ(const vector<double> &alpha, const vector<double> &beta);

    mat& gen_rand_matrix(double z2, mat &M) const;

    /*
     * Computes the value of Lambda(theta)
     */
    double estimateLambda(
	double theta, double &sd,
	unsigned long N = 80,
	unsigned long K = 10000
	) const;
	    
    double find_tail_index(
	double bounds[2],
	unsigned long N = 80,
	unsigned long K = 10000
	) const;

protected:
    // mutable random_device gen;
    mutable mt19937_64 gen;
};

#endif
