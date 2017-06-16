#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <functional>

using namespace std;
random_device gen;

// void merge_them(vector<unsigned> &V, const vector<unsigned> &pos)
// {
//     for (auto i = pos.begin() + 1; i != pos.end(); i++) {
// 	auto j = next(i) == pos.end() ? V.end() : V.begin() + *next(i);
// 	inplace_merge(V.begin(), V.begin() + *i, j);
//     }
// }
bool comp_vec(const vector<unsigned> *v, const vector<unsigned> *w)
{
    return v->size() < w->size();    
}

void merge_them(vector<vector<unsigned>*> &L, vector<unsigned> &V)
{
    unsigned n = accumulate(
	L.begin(), L.end(), (unsigned)0,
	[](unsigned m, const vector<unsigned> *v) {
	    return m + v->size();
	});
    V.resize(n);
    make_heap(L.begin(), L.end(), comp_vec);
    pop_heap(L.begin(), L.end(), comp_vec);
    auto j = prev(L.end());
    auto k = V.begin();
    copy((*j)->begin(), (*j)->end(), k);
    k += (*j)->size();
    for (unsigned c = L.size() - 1; c > 0; c--) {
	pop_heap(L.begin(), j);
	copy((*prev(j))->begin(), (*prev(j))->end(), k);
	auto m = k + (*prev(j))->size();
	inplace_merge(V.begin(), k, m);
	k = m;
	j = prev(j);
    }
}

int main(void)
{
    unsigned n = 0;
    unsigned k = 8;
    // vector<unsigned> pos(k);
    // for_each(V.begin(), V.end(),
    // 	     [](int &i) {
    // 		 i = gen() % 10;
    // 	     });
    // for_each(pos.begin() + 1, pos.end(),
    // 	     [=](unsigned &u) {
    // 		 u = (gen() % (n-1)) + 1;
    // 	     });
    // pos[0] = 0;

    // sort(pos.begin(), pos.end());
    // for (auto i = pos.begin(); i != pos.end(); i++) {
    // 	auto j = next(i) == pos.end() ? V.end() : V.begin() + *next(i);
    // 	sort(V.begin() + *i, j);
    // }
    // for_each(V.begin(), V.end(),
    // 	     [](int &x) {
    // 		 cout << x << ", ";
    // 	     });
    // cout << endl;
    // merge_them(V, pos);
    vector<vector<unsigned> *> L(k);
    for (auto i = L.begin(); i != L.end(); i = next(i)) {
	unsigned m = gen() % 10;
	vector<unsigned> *p = new vector<unsigned>(m);
	generate(p->begin(), p->end(),
		 [&](void) {
		     return gen() % 10;
		 });
	n += m;
	*i = p;
    }
    vector<unsigned> V(n);
    merge_them(L, V);
    for_each(V.begin(), V.end(),
	     [](int &x) {
		 cout << x << ", ";
	     });
    cout << endl;
    for_each(L.begin(), L.end(),
	     [](vector<unsigned> *p) {
		 delete p;
	     });
    return 0;
}
