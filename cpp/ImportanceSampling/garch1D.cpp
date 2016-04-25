#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <cmath>
#include <algorithm>

#include "garch1D.hpp"

using namespace std;

template<typename T>
class Garch1D
{
protected:
    T xi;
    T a;
    T b;
    set<T> Qset;
    set<T> Pset;
    normal_distribution<T> ndist;
    random_device gen;
    
public:
    /**
     * E |X|^t = 2^((t + 1)/2) * Gamma((t+1)/2)
     */
    T norm_moment(T t) {
	return 2^((t + 1)/2) * tgamma((t+1)/2);
    }

    T mean_A_to_xi (void) {
	T s = accumulate(Qset.begin(), Qset.end(), 0,
		   [=](T s, T x) {
		       size_t n = this->Qset.size();
		       return s += (this->a * x^2 + this->b)/n;
		   });
	return s;
    }
	    
    
    Garch1D(T a, T b) {
	this->a = a;
	this->b = b;
	
	// xi = find_tail_index(a, b);

	// T C = norm_moment(xi);
	size_t n = 1000000;

	Qset.resize(n);
	for (size_t i = 0; i < n; i++) {
	    Qset.insert(ndist(gen));
	}
	
    // 	Pset.resize(n);
    // 	typename set<T>::iterator i;
    // 	T s;
    // 	for (i = Qset.begin(), s = 0; i < Qset.end(); i++) {
    // 	    Pset.insert(s);
    // 	    s += pow(*i, xi)/n/C;
    // 	}
    // }

    // T density_func(T x) {
    // 	T y;
    // 	y = exp(-x/(2*b)) * pow(x, xi - 1/2);
    // 	y = y / sqrt(2*M_PI*b);
    // 	y = y/norm_moment(xi);
    // 	return y
    }
	
    
    T dist_func(T x) {
    	static T C = norm_moment(xi);
	
    	typename set<T>::iterator ub = Qset.upper_bound(x);
    	if (ub == Qset.begin())
    	    return 0;
    	T x0 = *(ub - 1);
    	T s = *(Pset.begin() + ub - Qset.begin() - 1) + pow(x0, xi) * gsl_ran_ugaussian_pdf(x0) * (x - x0);
    	return s / C;
    }
    
    T quantile_func(T u) {
	typename set<T>::iterator ub = Pset.upper_bound(u);
	T u0 = *(ub - 1);
	T q0 = *(Qset.begin() + ub - Pset.begin() - 1);

	return (u - u0)/gsl_ran_ugaussian_pdf(q0) + q0;
    }
};

    
