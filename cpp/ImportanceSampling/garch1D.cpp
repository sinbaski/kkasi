#include <time.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <set>
#include <algorithm>


using namespace std;

template<typename T>
class Garch1D
{
protected:
    T alpha;
    T a;
    T b;
    set<T> Qset;
    set<T> Pset;
    normal_distribution<T> ndist;
    default_random_engine gen;
    
public:
    /**
     * E |X|^t = 2^((t + 1)/2) * Gamma((t+1)/2)
     */
    T norm_moment(T t) {
	return 2^((t + 1)/2) * tgamma((t+1)/2);
    }
    
    Garch1D(T alpha, T a, T b)
	:gen(time(NULL)) {

    	T C = norm_moment(alpha);
	size_t n = 1000000;

	this->alpha = alpha;
	this->a = a;
	this->b = b;

	Qset.resize(n);
	Pset.resize(n);

	for (size_t i = 0; i < n; i++) {
	    Qset.insert(ndist(gen));
	}
	
	typename set<T>::iterator i;
	T s;
	for (i = Qset.begin(), s = 0; i < Qset.end(); i++) {
	    Pset.insert(s);
	    s += pow(*i, alpha)/n/C;
	}
    }

    // T density_func(T x) {
    // 	T y;
    // 	y = exp(-x/(2*b)) * pow(x, alpha - 1/2);
    // 	y = y / sqrt(2*M_PI*b);
    // 	y = y/norm_moment(alpha);
    // 	return y;
    // }
	
    
    T dist_func(T x) {
    	static T C = norm_moment(alpha);

    	typename set<T>::iterator ub = Qset.upper_bound(x);
    	if (ub == Qset.begin())
    	    return 0;
    	T x0 = *(ub - 1);
    	T s = *(Pset.begin() + ub - Qset.begin() - 1) + pow(x0, alpha) * gsl_ran_ugaussian_pdf(x0) * (x - x0);
    	return s / C;
    }
    
    T quantile_func(T u) {
	typename set<T>::iterator ub = Pset.upper_bound(u);
	T u0 = *(ub - 1);
	T q0 = *(Qset.begin() + ub - Pset.begin() - 1);

	return (u - u0)/gsl_ran_ugaussian_pdf(q0) + q0;
    }
};

    
