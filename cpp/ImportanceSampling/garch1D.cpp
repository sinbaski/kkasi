#include <time.h>
#include <cmath>
#include <algorithm>
#include "garch1D.hpp"

using namespace std;

template <typename T>
T Garch1D<T>::moment_func(T moment, T measure_shift, T M) const
{
    T n = (T)pool.size();
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

    typename vector<pool_data>::const_iterator i;
    for (i = pool.begin(); i < pool.end(); i++) {
	m1 += pow(i->quantile, moment + measure_shift)/n;
	if (flag) m2 += pow(i->quantile, measure_shift)/n;
    }
    return m1/m2;
}

template <typename T>
void Garch1D<T>::set_shift_par(T shift_par)
{
    bool flag = abs(shift_par) > 1.0e-3;
    this->shift_par = shift_par;
    normalizer = flag ? moment_func(shift_par) : 1;

    typename vector<pool_data>::iterator i;
    T s = 0;
    T n = (T)pool.size();
    i = pool.begin();
    do {
    	s += flag ? pow(i->quantile, shift_par)/n/normalizer : 1/n;
	(i++)->prob = s;
    } while (i < pool.end());
}

template <typename T>
Garch1D<T>::Garch1D(T a0, T a, T b)
    :dist(0.5, 2*a), a0(a0), a(a), b(b)
{
    /*
      Set up the distribution of A
     */
    typename vector<pool_data>::iterator i;
    size_t n = 1000000;

    pool.resize(n);
    for (i = pool.begin(); i < pool.end(); i++) {
	i->quantile = dist(gen) + b;
	i->prob = nan("");
    }
    sort(pool.begin(), pool.end());

    /*
      Set up the stationary distribution 
     */
    stationary.resize(10000);
    size_t N = 2000;
    for (i = stationary.begin(); i < stationary.end(); i++) {
	double Vn = a0;
	for (int step = 0; step < N; step++) {
	    Vn = Vn * (dist(gen) + b) + a0;
	}
	i->quantile = Vn;
	i->prob = nan("");
    }
    sort(stationary.begin(), stationary.end());
    T k = (T)stationary.size();
    for (n = 0, i = stationary.begin(); n < stationary.size(); n++, i++) {
	i->prob = ((T)n+1)/k;
    }
}

template <typename T>
T Garch1D<T>::stationary_quantile(T u) const
{
    pool_data e;
    e.prob = u;
    e.quantile = nan("");

    typename vector<pool_data>::const_iterator ub =
	upper_bound(pool.begin(), pool.end(), e);
    return (ub - 1)->quantile;
}

template <typename T>
T Garch1D<T>::stationary_prob(T x) const
{
    pool_data e;
    e.prob = nan("");
    e.quantile = x;

    typename vector<pool_data>::const_iterator ub =
	upper_bound(pool.begin(), pool.end(), e);
    return (ub - 1)->prob;
}

template <typename T>
T Garch1D<T>::init_quantile(T u) const
{
    T u1 = stationary_prob(M);
    return stationary_quantile(u * u1);
}

template <typename T>
T Garch1D<T>::init_prob(T x) const
{
    return stationary_prob(x)/stationary_prob(M);
}

template<typename T>
T Garch1D<T>::density_func(T x) const 
{
    T y;
    y = pow(x, shift_par) / sqrt(x - b) * exp(-(x-b)/(2*a));
    y = y / sqrt(2*M_PI*a);
    y = y/normalizer;
    return y;
}
	
template<typename T>
T Garch1D<T>::dist_func(T x) const
{
    pool_data e;
    e.quantile = x;
    e.prob = nan("");

    typename vector<pool_data>::const_iterator ub =
	upper_bound(pool.begin(), pool.end(), e);
    if (ub == pool.begin())
	return 0;
    T x0 = (ub - 1)->quantile;
    T s = (ub - 1)->prob + pow(x0, shift_par) * density_func(x0) * (x - x0);
    return s / normalizer;
}
    
template<typename T>
T Garch1D<T>::quantile_func(T u) const
{
    pool_data e;
    e.prob = u;
    e.quantile = nan("");

    typename vector<pool_data>::const_iterator ub =
	upper_bound(pool.begin(), pool.end(), e);
    T u0 = (ub - 1)->prob;
    T q0 = (ub - 1)->quantile;
    return (u - u0)/density_func(q0) + q0;
}

template class Garch1D<double>;


