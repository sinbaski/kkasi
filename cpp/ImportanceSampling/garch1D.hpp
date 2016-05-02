#ifndef garch1D_HPP
#define garch1D_HPP

#include <random>
#include <vector>
#include <cmath>
#include <ctgmath>
#include <armadillo>

using namespace std;
using namespace arma;

#define SAMPLE_SIZE 1000000

template <typename T>
class Garch1D
{
public:
    /**
     * E |X|^t = 2^((t + 1)/2) * Gamma((t+1)/2)
     */
    struct pool_data {
	T quantile;
	T prob;
	bool operator< (const pool_data& rhs) const;
    };

    const T a0;
    const T a;
    const T b;

    const static int ORIG_M = 1;
    const static int SHIFTED_M = 2;

    T M;
    T xi;
    int measure_index;

    Garch1D(T a, T b, T a0);

    /**
     * The original or the shifted distribution of A
     */
    T moment_func(T moment, T measure_shift = 0, T M = 0) const;
    T density_func(T x) const;
    T dist_func(T x) const;
    T quantile_func(T u) const;

    /**
     * Stationary distribution of V
     */
    T stationary_quantile(T u) const;
    T stationary_prob(T x) const;

    /*
      Initial distribution
     */
    T init_quantile(T u) const;
    T init_prob(T x) const;

protected:
    Mat<T> pool;
    vector<pool_data> stationary;
    gamma_distribution<T> dist;
    random_device gen;

    int find_tail_index(void);
};

template <typename T>
struct LHS_func_par {
    const Garch1D<T> *garch11;
};
    
template <typename T>
inline T LHS_func(T alpha, void *par) {
const Garch1D<T> *garch11 = ((LHS_func_par<T> *)par)->garch11;
return garch11->moment_func(alpha) - 1;
}
  
#endif
