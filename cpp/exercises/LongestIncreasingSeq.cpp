#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <iterator>
#include <cstdlib>
#include <utility>
#include <type_traits>

using namespace std;
//random_device gen;
unsigned gen(void)
{
    return (unsigned)rand();
}

// bool identical_seq(vector<unsigned> *u, vector<unsigned> *v)
// {
//     if (u->size() != v->size()) return false;
//     vector<unsigned> T(u->size());
//     set_symmetric_difference(
// 	u->begin(), u->end(),
// 	v->begin(), v->end(),
// 	T.begin()); 
//     return T.empty();
// }

bool comp_seq(vector<unsigned> *u, vector<unsigned> *v)
{
    if (u->back() < v->back())
	return true;
    if (u->back() > v->back())
	return false;
    if (u->size() > v->size())
	return true;
    if (u->size() < v->size())
	return false;
    return lexicographical_compare(
	u->begin(), u->end(),
	v->begin(), v->end()
	);
}

class SeqComp
{
public:
    bool operator() (vector<unsigned> *u, vector<unsigned> *v) const {
	return comp_seq(u, v);
    }
    static bool comp_size(vector<unsigned> *u, vector<unsigned> *v) {
	return u->size() < v->size();	
    }
};

typedef set<vector<unsigned> *, SeqComp> MySet;

void dump_results(const MySet &U)
{
    for_each(U.begin(), U.end(),
	     [](vector<unsigned> *v) {
		 for_each(v->begin(), v->end(),
			  [](unsigned x) {
			      printf("%u, ", x);
			  });
		 printf("\n");
	     });
    printf("\n");
}

void extend_it_maybe(const vector<unsigned> *u, unsigned val, MySet &U)
{
    vector<unsigned> *v = new vector<unsigned>(*u);
    v->push_back(val);
    auto result = U.insert(v);
    if (!result.second) delete v;
}

template<class Ri>
void find_seq(Ri first, Ri last, unsigned depth, MySet &U)
{
    // static_assert(is_base_of<forward_iterator_tag, Ri>::value, "A random access iterator is needed here.");
    if (distance(first, last) == 1) {
	vector<unsigned> *v = new vector<unsigned>(first, last);
	U.insert(v);
	return;
    }
    MySet W;
    Ri p = prev(last);
    find_seq(first, prev(last), depth + 1, W);
    U = W;
    
    vector<unsigned> temp(p, next(p));
    auto i =
	lower_bound(W.begin(), W.end(), &temp,
		    [](vector<unsigned> *u, vector<unsigned> *v) {
			return u->back() < v->back();
		    });
    auto j = max_element(W.begin(), i, SeqComp::comp_size);
    while (distance(j, i) > 0) {
	extend_it_maybe(*j, *p, U);
	j = find_if(next(j), i,
		    [&](vector<unsigned> *v) {
			return v->size() == (*j)->size();
		    });
    };
    vector<unsigned> *v = new vector<unsigned>(p, last);
    if (!U.insert(v).second) delete v;

    auto bounds =
    	equal_range(U.begin(), U.end(), &temp,
    		    [](vector<unsigned> *u, vector<unsigned> *v) {
    			return u->back() < v->back();
    		    });
    auto k =
	upper_bound(bounds.first, bounds.second,
		    *bounds.first,
    		    [](vector<unsigned> *u, vector<unsigned> *v) {
    			return u->size() > v->size();
    		    });
    while (distance(k, U.end()) > 0) {
	auto i =
	    find_if(next(k), U.end(),
		    [&](vector<unsigned> *v) {
			bool b = (v->back() >= (*bounds.first)->back());
			b = (b && (v->size() < (*bounds.first)->size()));
			return b;
		    });
	delete *k;
	U.erase(k);
	k = i;
    }

    if (depth == 0) {
	auto i =
	    max_element(U.begin(), U.end(),
			[](vector<unsigned>* u, vector<unsigned>* v) {
			    return u->size() < v->size();
			});
	if (distance(U.begin(), i) > 0) {
	    for_each(U.begin(), i,
		     [](vector<unsigned> *v) {
			 delete v;
		     });
	    U.erase(U.begin(), i);
	}
	for (i = U.begin(); distance(i, U.end()) > 0;) {
	    auto j =
		find_if(next(i), U.end(),
			[&i](vector<unsigned> *u) {
			    return (*i)->size() == u->size();
			});
	    for_each(next(i), j,
		     [](vector<unsigned> *v) {
			 delete v;
		     });
	    if (distance(i, j) > 1)
		i = U.erase(next(i), j);
	    else i = j;
	}
    }
    printf("At depth %u: #%u\n", depth, *p);
    dump_results(U);
}

int main(int argc, char *argv[])
{
    vector<unsigned> V(stoi(argv[1]));
    MySet U;

    srand(stoi(argv[2]));
    generate(V.begin(), V.end(),
	     [](void) {
		 return gen() % 10;
	     });
    for_each(V.begin(), V.end(),
	     [](unsigned x) {
		 printf("%u, ", x);
	     });
    printf("\n\n");
    find_seq(V.begin(), V.end(), 0, U);
    for_each(U.begin(), U.end(),
    	     [](vector<unsigned> *seq) {
		 for_each(seq->begin(), seq->end(),
			  [](unsigned x) {
			      printf("%u, ", x);
			  });
		 printf("\n");
		 delete seq;
    	     });
    return 0;
}
