#include <time.h>
#include <cmath>
#include <algorithm>
#include <set>

using namespace std;

template<typename T>
class Garch1D
{
protected:
    T alpha;
    T a;
    T b;
    multiset<T> pool;
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
	sizt_t n = 1000000;

	this->alpha = alpha;
	this->a = a;
	this->b = b;

	pool.resize(n);
	for (size_t i = 0; i < n; i++) {
	    pool.insert(ndist(gen));
	}
    }

    T density_func(T x) {
	T y;
	y = exp(-x/(2*b)) * pow(x, alpha - 1/2);
	y = y / sqrt(2*M_PI*b);
	y = y/norm_moment(alpha);
	return y;
    }
	
    
    T dist_func(T x) {
	static T C = norm_moment(alpha);
	multiset<T>::const_iterator ub = pool.upper_bound(x);
	if (ub == pool.begin())
	    return 0;
	T x0 = *(ub - 1);
	sizt_t n = ub - pool.begin();
	T s = accumulate(pool.begin(), ub, 0,
			 [&n](T s, T y) {
			     return x + pow(y, alpha)/n;
			 }
	    );
	pow(x0, alpha)/C * 
	return s;
    }
    
    T quantile_func(T q) {
	return 0;
    }
};

    
