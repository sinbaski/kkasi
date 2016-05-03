#include <time.h>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "garch1D.hpp"

using namespace std;

template <typename T>
T Garch1D<T>::moment_func(T moment, T measure_shift, T normalizer) const
{
    T n = (T)pool.n_rows;
    T e = 1.0e-3;
    T m1, m2;
    bool flag = false;
    if (abs(moment) <= e) {
	return 1;
    } else if (abs(measure_shift) <= e) {
	m1 = 0;
	m2 = 1;
    } else if (normalizer > e) {
	m1 = 0;
	m2 = normalizer;
    } else {
	m1 = m2 = 0;
	flag = true;
    }

    for_each(pool.begin_col(0), pool.end_col(0),
	     [&](T x) {
		 m1 += pow(x, moment + measure_shift)/n;
		 if (flag) m2 += pow(x, measure_shift)/n;
	     });
    return m1/m2;
}

template <typename T>
struct LHS_func_par {
    const Garch1D<T> *garch11;
};

template<typename T>
T LHS_func(T xi, void *par)
{
    const Garch1D<T> *garch11 = ((LHS_func_par<T> *)par)->garch11;
    return garch11->moment_func(xi) - 1;
}

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
    F.function = LHS_func<T>;
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
	status = gsl_root_test_interval(lb, ub, 1.0e-6, 0);
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
Garch1D<T>::Garch1D(T a0, T a, T b, unsigned long sample_size)
    :a0(a0), a(a), b(b), xi(nan("")), dist(0.5, 2*a)
{
    /*
      Set up the distribution of A
     */
    pool.set_size(sample_size, 3);
    for_each(pool.begin_col(0), pool.end_col(0),
	     [this, b](T &x) {
		 x = this->dist(this->gen) + b;
	     });
    sort(pool.begin_col(0), pool.end_col(0));

    find_tail_index();
    
    /*
      Set up the stationary distribution 
     */
    size_t N = 2000;
    stationary.resize(10000);
#pragma omp parallel for
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
    printf("Left end point of stationary dist = %.4e\n", stationary[0].quantile);

    T k = (T)stationary.size();

#pragma omp parallel for
    for (unsigned int i = 0; i < stationary.size(); i++) {
	stationary[i].prob = (T)(i+1)/k;
    }

    /*
      Set up the probs in both the original and the shifted measure
    */
    unsigned int i = 0;
    T s0, s1, x0;
    s0 = s1 = 0;
    x0 = 1/(double)pool.n_rows;
    for (i = 0, s0 = 0, s1 = 0, x0 = 1/(double)pool.n_rows; i < pool.n_rows; i++) {
    	s0 += x0;
	s1 += pow(pool(i, 0), xi) * x0;
	pool(i, ORIG_M) = s0;
	pool(i, SHIFTED_M) = s1;
    }
    
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
    static T u1 = nan("");
    if (!isnormal(u1))
	u1 = stationary_prob(M);
    return stationary_quantile(u * u1);
}

template <typename T>
T Garch1D<T>::init_prob(T x) const
{
    static T uM = nan("");
    if (!isnormal(uM))
	uM = stationary_prob(M);
    return stationary_prob(x)/uM;
}

/*
  Density function of the shifted distribution
 */
template<typename T>
T Garch1D<T>::density_func(T x, int measure_index) const 
{
    T y;
    y = 1 / sqrt(x - b) * exp(-(x-b)/(2*a));
    if (measure_index == SHIFTED_M)
	y *= pow(x, xi);
    y = y / sqrt(2*M_PI*a);
    return y;
}
	
template<typename T>
T Garch1D<T>::dist_func(T x, int measure_index) const
{
    T q0, u0;
    typename Mat<T>::const_col_iterator i =
	upper_bound(pool.begin_col(0), pool.end_col(0), x);
    if (i == pool.begin_col(0)) {
	q0 = pool(0, 0);
	u0 = pool(0, measure_index);
    } else {
	q0 = *(i-1);
	u0 = pool(i-pool.begin_col(0)-1, measure_index);
    }
    return u0 + density_func(q0, measure_index) * (x - q0);
}
    
template<typename T>
T Garch1D<T>::quantile_func(T u, int measure_index) const
{
    T q0, u0;
    typename Mat<T>::const_col_iterator i =
	upper_bound(pool.begin_col(measure_index),
		    pool.end_col(measure_index), u);
    if (i == pool.begin_col(measure_index)) {
	q0 = pool(0, 0);
	u0 = pool(0, measure_index);
    } else {
	u0 = *(i - 1);
	q0 = pool(i-pool.begin_col(measure_index)-1, 0);
    }

    return q0 + (u - u0)/density_func(q0, measure_index);
}

template class Garch1D<double>;


