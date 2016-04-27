#ifndef garch1D_HPP
#define garch1D_HPP

#include <random>
#include <vector>

using namespace std;

template <typename T>
class Garch1D
{
protected:
    T shift_par = 0;
    T a0;
    T a;
    T b;
    T normalizer = 1;
    T M;
    vector<T> quantiles;
    vector<T> probs;
    gamma_distribution<T> dist;
    random_device gen;
    
public:
    /**
     * E |X|^t = 2^((t + 1)/2) * Gamma((t+1)/2)
     */

    T moment_func(T moment, T measure_shift = 0, T M = 0);
    Garch1D(T a, T b, T a0);
    void set_shift_par(T shift_par);
    T get_shift_par(void) {
	return shift_par;
    };
    void set_M(T M) {
	this-> M = M;
    }
    T M_func(T alpha);
    T density_func(T x);
    T dist_func(T x);
    T quantile_func(T u);
};

    
#endif
