#include <algorithm>
#include <numeric>
#include <vector>
#include <list>
#include <set>
#include <bits/stdc++.h>

using namespace std;

typedef vector<long>::iterator Iter;

long getWays(long n,  const Iter first, const Iter last){
    // Complete this function
    if (distance(first, last) == 1) {
	return n % *first == 0 ? 1 : 0;
    }
    long l = *prev(last);
    long k = n / l;
    long result = 0;
    for (long i = 0; i <= k; i++) {
	long w = getWays(n - i * l, first, prev(last));
	result += w;
    }
    return result;
}

int main() {
    int n;
    int m;
    cin >> n >> m;
    vector<long> c(m);
    for(int c_i = 0; c_i < m; c_i++){
       cin >> c[c_i];
    }
    // Print the number of ways of making change for 'n' units using coins having the values given by 'c'
    long ways = getWays(n, c.begin(), c.end());
    cout << ways << endl;
    return 0;
}
