#ifndef garch1D_HPP
#define garch1D_HPP

#include <random>
#include <vector>
#include <cmath>
#include <ctgmath>

using namespace std;

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
	inline bool operator< (const pool_data& rhs) const {
	    if (isnormal(quantile) && isnormal(rhs.quantile))
		return quantile < rhs.quantile;
	    else
		return prob < rhs.prob;
	}
    };
    const T a0;
    const T a;
    const T b;

    T M;

    Garch1D(T a, T b, T a0);

    void set_shift_par(T shift_par);
    inline T get_shift_par(void) {
	return shift_par;
    };

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
    T shift_par = 0;
    T normalizer = 1;
    vector<pool_data> pool;
    vector<pool_data> stationary;
    gamma_distribution<T> dist;
    random_device gen;
};

  
#endif
