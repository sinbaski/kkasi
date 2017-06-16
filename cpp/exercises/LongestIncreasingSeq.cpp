#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <list>
#include <iterator>
#include <cstdlib>
#include <utility>

using namespace std;
//random_device gen;
unsigned gen(void)
{
    return (unsigned)rand();
}

void dump_results(const list< vector<unsigned> *> &U)
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

bool identical_seq(vector<unsigned> *u, vector<unsigned> *v)
{
    if (u->size() != v->size()) return false;
    vector<unsigned> T(u->size());
    set_symmetric_difference(u->begin(), u->end(), v->begin(), v->end(), T.begin());
    return T.empty();
}

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
	v->begin(), v->end(),
	);
}

template<class Ri>
void find_seq(Ri first, Ri last, unsigned depth, list< vector<unsigned> *> &U)
{
    if (distance(first, last) == 1) {
	vector<unsigned> *v = new vector<unsigned>(first, last);
	U.push_back(v);
	return;
    }
    list< vector<unsigned> *> W;
    unsigned extended = 0;
    find_seq(first, prev(last), depth + 1, W);
    Ri p = prev(last);
    for (auto i = W.begin(); i != W.end();) {
	if ((*i)->back() < *p) {
	    vector<unsigned> *v = new vector<unsigned>(*(*i));
	    v->push_back(*p);
	    // check for a range of equal last  elem and equal length
	    
	    U.push_back(v);
	    extended++;
	    if (depth == 0) {
		delete *i;
		i = W.erase(i);
	    } else {
		U.push_back(*i);
		i++;
	    }
	} else {
	    U.push_back(*i);
	    i++;
	}
    }
    if (!extended) U.push_back(new vector<unsigned>(p, p+1));
    U.sort(comp_seq);

    printf("At depth %u: #%u\n", depth, *p);
    dump_results(U);

    for (auto i = U.begin(); i != U.end();) {
	auto j =
	    upper_bound(i, U.end(), *i,
			[&](vector<unsigned> *u, vector<unsigned> *v) {
			    return u->back() < v->back();
			});
	if (distance(i, j) == 1) {
	    i = next(i);
	    continue;
	}
	auto k =
	    upper_bound(i, j, *i,
			[&](vector<unsigned> *u, vector<unsigned> *v) {
			    return u->size() > v->size();
			});
	if (distance(k, j)) {
	    for_each(k, j,
		     [](vector<unsigned> *v) {
			 delete v;
		     });
	    U.erase(k, j);	    
	}
	for (auto l = next(i); l != k;) {
	    if (identical_seq(*i, *l)) {
		delete *l;
		l = U.erase(l);
	    }
	    l = next(l);
	}
	i = j;
    }
    for (auto i = U.begin(); i != U.end();) {
	auto j =
	    upper_bound(i, U.end(), *i,
			[](vector<unsigned> *u, vector<unsigned> *v) {
			    return u->back() < v->back();
			});
	if (j == U.end()) break;
	auto k =
	    upper_bound(j, U.end(), *j,
			[](vector<unsigned> *u, vector<unsigned> *v) {
			    return u->back() < v->back();
			});
	if ((*j)->size() < (*i)->size()) {
	    for_each(j, k,
		     [](vector<unsigned> *v) {
			 delete v;
		     });
	    i = U.erase(j, k);
	}
	i = j;
    }
    if (depth == 0) {
	   U.sort(
		[](vector<unsigned> *u, vector<unsigned> *v) {
		    return u->size() > v->size();
		});
	   auto i = upper_bound(next(U.begin()), U.end(),
				*U.begin(),
				[](vector<unsigned> *u, vector<unsigned> *v) {
				    return u->size() > v->size();
				});
	   for_each(i, U.end(),
		    [](vector<unsigned> *v) {
			delete v;
		    });
	   U.erase(i, U.end());
    }
    dump_results(U);
}

int main(int argc, char *argv[])
{
    vector<unsigned> V(stoi(argv[1]));
    list< vector<unsigned> *> U;

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
