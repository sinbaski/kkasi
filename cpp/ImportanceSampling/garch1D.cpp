#include <time.h>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_min.h>
#include "garch1D.hpp"

using namespace std;

template <typename T>
T Garch1D<T>::moment_func(T moment, T measure_shift, T M)
{
    T n = (T)quantiles.size();
    T e = 1.0e-3;
    T m1, m2;
    bool flag = false;
    if (abs(moment) <= e) {
	return 1;
    } else if (abs(measure_shift) <= e) {
	m1 = 0;
	m2 = 1;
    } else if (M > e) {
	m1 = 0;
	m2 = M;
    } else {
	m1 = m2 = 0;
	flag = true;
    }

    typename vector<T>::iterator i;
    for (i = quantiles.begin(); i < quantiles.end(); i++) {
	m1 += pow(*i, moment + measure_shift)/n;
	if (flag) m2 += pow(*i, measure_shift)/n;
    }
    return m1/m2;
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
Garch1D<T>::Garch1D(T a0, T a, T b)
    :dist(0.5, 2*a)
{
    this->a0 = a0;
    this->a = a;
    this->b = b;
	
    size_t n = 1000000;

    quantiles.resize(n);
    for (size_t i = 0; i < n; i++) {
	quantiles[i] = dist(gen) + b;
    }
    sort(quantiles.begin(), quantiles.end());
	
}

template<typename T>
T Garch1D<T>::M_func(T alpha)
{
    T lambda = moment_func(alpha, shift_par, normalizer);
    if (alpha < 1) {
	return a0 / pow(1 - lambda, 1/alpha);
    } else {
	return a0 / (1 - pow(lambda, 1/alpha));
    }
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
    typename vector<T>::iterator ub = upper_bound(quantiles.begin(), quantiles.end(), x);
    if (ub == quantiles.begin())
	return 0;
    T x0 = *(ub - 1);
    T s = *(probs.begin() + (ub - quantiles.begin() - 1)) +
	pow(x0, shift_par) * density_func(x0) * (x - x0);
    return s / normalizer;
}
    
template<typename T>
T Garch1D<T>::quantile_func(T u) {
    typename vector<T>::iterator ub = upper_bound(probs.begin(), probs.end(), u);
    T u0 = *(ub - 1);
    T q0 = *(quantiles.begin() + (ub - probs.begin() - 1));
    return (u - u0)/density_func(q0) + q0;
}

template class Garch1D<double>;


