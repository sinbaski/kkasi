#include <iostream>
#include <random>
#include <algorithm>

using namespace std;
random_device gen;

int find_mode(vector<unsigned int> &V)
{
    sort(V.begin(), V.end());
    unsigned int mode, n;
    auto i = V.begin();
    mode = *i;
    n = 0;
    while (i != V.end()) {
	unsigned y = *i;
	auto j = find_if(i, V.end(),
			 [=](unsigned x) {
			     return x != y;
			 });
	if (j - i > n) {
	    n = j - i;
	    mode = y;
	}
	i = j;
    }
    return mode;
}

int main(void)
{
    vector<unsigned int> V(20), W(20);
    for_each(V.begin(), V.end(),
	     [](unsigned int &a) {
		 a = gen() % 5;
	     });
    W = V;
    sort(W.begin(), W.end());
    for_each(W.begin(), W.end(),
	     [](unsigned i) {
		 cout << i << ", ";
	     });
    cout << endl;
    unsigned mode = find_mode(V);
    cout << "The mode is " << mode << endl;
    return 0;
}
