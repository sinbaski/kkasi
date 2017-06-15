#include <algorithm>
#include <set>
#include <iostream>

using namespace std;
bool check_it(const set<unsigned> &S1, const set<unsigned> &S2, unsigned x)
{
    for (auto i = S1.begin(); i != S1.end(); i = next(i)) {
	if (binary_search(S2.begin(), S2.end(), x - *i))
	    return true;
    }
    return false;
}

int main(void)
{
    unsigned n = 10;
    unsigned x = 60;
    set<unsigned> S1, S2;
    for (unsigned i = 0; i < n; i++) {
	S1.insert(i);
	S2.insert(i + 10);
    }
    cout << check_it(S1, S2, x) << " for " << x << endl;
    return 0;
}

