#include <time.h>
#include <cmath>
#include <algorithm>

#include "garch1D.hpp"

using namespace std;

template <typename T>
T Garch1D<T>::moment_func(T moment) {
    T n = (T)quantiles.size();
    T mean = 0;
    // for (typename vector<T>::iterator i = quantiles.begin(); i < quantiles.end(); i++) {
    // 	mean += pow(*i, moment)/n;
    // }
    if (normalizer <= 0) {
	mean = accumulate(
	    quantiles.begin(),
	    quantiles.end(),
	    mean,
	    [=](T s, T x) {
		return s + pow(x, moment)/n;
	    }
	    );
    } else {
	mean = accumulate(
	    quantiles.begin(),
	    quantiles.end(),
	    mean,
	    [=](T s, T x) {
		return s + pow(x, moment + shift_par)/n;
	    }
	    ) / normalizer;
    }
    return mean;
}

template <typename T>
void Garch1D<T>::set_shift_par(T shift_par)
{
    this->shift_par = shift_par;
    normalizer = moment_func(shift_par);

    probs.resize(quantiles.size());
    typename vector<T>::iterator i, j;
    T s;
    T n = (T)quantiles.size();
    for (i = quantiles.begin(), j = probs.begin(), s = 0;
	 i < quantiles.end(); i++, j++) {
    	*j = s;
    	s += pow(*i, shift_par)/n/normalizer;
    }
}

template <typename T>
Garch1D<T>::Garch1D(T a, T b)
    :dist(0.5, 2*a)
{
    this->a = a;
    this->b = b;
	
    // shift_par = find_tail_index(a, b);

    // T C = norm_moment(shift_par);
    size_t n = 1000000;

    quantiles.resize(n);
    for (size_t i = 0; i < n; i++) {
	quantiles[i] = dist(gen) + b;
    }
    sort(quantiles.begin(), quantiles.end());
	
}

template<typename T>
T Garch1D<T>::density_func(T x) {
    T y;
    y = pow(x, shift_par) / sqrt(x - b) * exp(-(x-b)/(2*a));
    y = y / sqrt(2*M_PI*a);
    y = y/normalizer;
    return y;
}
	
    
template<typename T>
T Garch1D<T>::dist_func(T x) {
    static T C = moment_func(shift_par);
	
    typename vector<T>::iterator ub = upper_bound(quantiles.begin(), quantiles.end(), x);
    if (ub == quantiles.begin())
	return 0;
    T x0 = *(ub - 1);
    T s = *(probs.begin() + (ub - quantiles.begin() - 1)) +
	pow(x0, shift_par) * density_func(x0) * (x - x0);
    return s / C;
}
    
template<typename T>
T Garch1D<T>::quantile_func(T u) {
    typename vector<T>::iterator ub = upper_bound(probs.begin(), probs.end(), u);
    T u0 = *(ub - 1);
    T q0 = *(quantiles.begin() + (ub - probs.begin() - 1));
    return (u - u0)/density_func(q0) + q0;
}

template class Garch1D<double>;


