#include <algorithm>
#include <iostream>
#include <random>

using namespace std;

random_device gen;

int main(void)
{
    unsigned n = 10;
    vector<int> V(n);
    vector<int> W;
    for_each(V.begin(), V.end(),
	     [](int &i) {
		 i = gen() % 5;
		 cout << i << ", ";
	     });
    cout << endl;
    sort(V.begin(), V.end());
    auto i = V.begin();
    while (i != V.end()) {
	auto j = upper_bound(i, V.end(), *i);
	if (j - i > floor(n/4)) {
	    W.push_back(*i);
	}
	i = j;
    }
    for_each(W.begin(), W.end(),
	     [](int x) {
		 cout << x << ", ";
	     });
    cout << endl;
}
