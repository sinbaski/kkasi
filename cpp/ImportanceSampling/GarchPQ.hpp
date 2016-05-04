#ifndef garch1D_HPP
#define garch1D_HPP

#include <random>
#include <vector>
#include <cmath>
#include <ctgmath>
#include <armadillo>

using namespace std;
using namespace arma;

class Garch1D
{
public:
    const vector<double> alpha;
    const vector<double> beta;

    const static int ORIG_M = 1;
    const static int SHIFTED_M = 2;

    double M;
    double xi;

    GarchPQ(vector<double> &alpha, vector<double> &beta, unsigned long sample_size);

    /**
     * The original or the shifted distribution of A
     */
    double moment_func(double moment, double measure_shift = 0, double normalizer = 0) const;
    double density_func(double x, int measure_index) const;
    double dist_func(double x, int measure_index) const;
    double quantile_func(double u, int measure_index) const;

    /**
     * Stationary distribution of V
     */
    double stationary_quantile(double u) const;
    double stationary_prob(double x) const;

    /*
      Initial distribution
     */
    double init_quantile(double u) const;
    double init_prob(double x) const;


    int find_tail_index(void);

protected:
    random_device gen;
    gamma_distribution<double> dist;
    Mat<double> pool;
    vector<pool_data> stationary;
};

#endif
