#include <iostream>
#include <algorithm>
#include <random>

using namespace std;
random_device gen;

int main(void)
{
    unsigned n = 20;
    unsigned k = 8;
    vector<int> V(n);
    for_each(V.begin(), V.end(),
	     [&](int &x) {
		 x = gen() % 10;
		 cout << x << ", ";
	     });
    cout << endl;

    make_heap(V.begin(), V.end());
    auto iter = V.end();
    for (unsigned i = 0; i < k; i++) {
	pop_heap(V.begin(), iter);
	cout << *(--iter) << ", ";
    }
    cout << endl;
    return 0;
}
