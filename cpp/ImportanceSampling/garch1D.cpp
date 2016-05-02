#include <time.h>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
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

    typename Mat<T>::const_iterator i;
    for_each(pool.begin_col(0), pool.end_col(0),
	     [&](T x) {
		 m1 += pow(x, moment + measure_shift)/n;
		 if (flag) m2 += pow(x, measure_shift)/n;
	     });
    return m1/m2;
}

// template <typename T>
// void Garch1D<T>::set_shift_par(T shift_par)
// {
//     bool flag = abs(shift_par) > 1.0e-3;
//     this->shift_par = shift_par;
//     normalizer = flag ? moment_func(shift_par) : 1;

//     typename vector<pool_data>::iterator i;
//     T s = 0;
//     T n = (T)pool.size();
//     i = pool.begin();
//     do {
//     	s += flag ? pow(i->quantile, shift_par)/n/normalizer : 1/n;
// 	(i++)->prob = s;
//     } while (i < pool.end());
// }

template<typename T>
int Garch1D<T>::find_tail_index(void)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    struct LHS_func_par<T> param = {this};
    F.function = &LHS_func;
    F.params = &param;
    
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, 1.8, 1.9);

    printf ("using %s method\n", 
	    gsl_root_fsolver_name (solver));

    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0001, 0);
	if (status == GSL_SUCCESS)
	    cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);

    if (status != GSL_SUCCESS) {
	cout << "The Brent algorithm did not converge after " << max_iter
	     << " iterations." << endl;
    }
    return status;
}

template <typename T>
Garch1D<T>::Garch1D(T a0, T a, T b)
    :a0(a0), a(a), b(b), xi(nan("")),
     measure_index(Garch1D<T>::ORIG_M),
     dist(0.5, 2*a), pool(SAMPLE_SIZE, 3)
{
    /*
      Set up the distribution of A
     */
    for_each(pool.begin_col(0), pool.end_col(0),
	     [this, b](T &x) {
		 x = this->dist(this->gen) + b;
	     });
    sort(pool.begin_col(0), pool.end_col(0));
    
    /*
      Set up the stationary distribution 
     */
    stationary.resize(10000);
    size_t N = 2000;
    for (typename vector<pool_data>::iterator i = stationary.begin();
	 i < stationary.end(); i++) {
	double Vn = a0;
	for (unsigned int step = 0; step < N; step++) {
	    Vn = Vn * (dist(gen) + b) + a0;
	}
	i->quantile = Vn;
	i->prob = nan("");
    }
    sort(stationary.begin(), stationary.end());
    T k = (T)stationary.size();
    int n = 1;
    for_each(stationary.begin(), stationary.end(),
	     [&n, k](pool_data &d) {
		 d.prob = (T)n/k;
		 n++;
	     }
	);
    /*
      The value of xi is set when the tail index is found.
    */
    find_tail_index();

    /*
      Set up the probs in both the original and the shifted measure
    */
    T s0 = 0, s1 = 0;
    T x0 = 1/(double)n;
    int i = 0;
    do {
    	s0 += x0;
	s1 += pow(pool[i, 0], xi) * x0;
	pool[i, 1] = s0;
	pool[i, 2] = s1;
    } while (++i < pool.n_cols);
    
}

template <typename T>
bool Garch1D<T>::pool_data::operator< (const Garch1D<T>::pool_data & rhs) const
{
    if (isnormal(quantile) && isnormal(rhs.quantile))
	return quantile < rhs.quantile;
    else
	return prob < rhs.prob;
    
}

template <typename T>
T Garch1D<T>::stationary_quantile(T u) const
{
    pool_data e;
    e.prob = u;
    e.quantile = nan("");

    typename vector<pool_data>::const_iterator ub =
	upper_bound(stationary.begin(), stationary.end(), e);
    if (ub == stationary.begin())
	return stationary.begin()->quantile;
    else
	return (ub - 1)->quantile;
}

template <typename T>
T Garch1D<T>::stationary_prob(T x) const
{
    pool_data e;
    e.prob = nan("");
    e.quantile = x;

    typename vector<pool_data>::const_iterator ub =
	upper_bound(stationary.begin(), stationary.end(), e);
    if (ub == stationary.begin())
	return 0;
    else
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

/*
  Density function of the shifted distribution
 */
template<typename T>
T Garch1D<T>::density_func(T x) const 
{
    T y;
    y = 1 / sqrt(x - b) * exp(-(x-b)/(2*a));
    if (measure_index == SHIFTED_M)
	y *= pow(x, xi);
    y = y / sqrt(2*M_PI*a);
    return y;
}
	
template<typename T>
T Garch1D<T>::dist_func(T x) const
{
    T q0, u0;
    int l;
    typename Mat<T>::const_col_iterator i =
	upper_bound(pool.begin_col(0), pool.end_col(0), x);
    if (i == pool.begin_col(0)) {
	q0 = pool[0, 0];
	u0 = pool[0, measure_index];
    } else {
	l = 1;
	q0 = *(i-l);
	u0 = pool[i-pool.begin_col(0)-l, measure_index];
    }
    return u0 + density_func(q0) * (x - q0);
}
    
template<typename T>
T Garch1D<T>::quantile_func(T u) const
{
    T q0, u0;
    typename Mat<T>::const_col_iterator i =
	upper_bound(pool.begin_col(measure_index),
		    pool.end_col(measure_index), u);
    if (i == pool.begin_col(measure_index)) {
	q0 = pool[0, 0];
	u0 = pool[0, measure_index];
    } else {
	u0 = *(i - 1);
	q0 = pool[i-pool.begin_col(measure_index)-1, 0];
    }

    return q0 + (u - u0)/density_func(q0);
}

template class Garch1D<double>;


