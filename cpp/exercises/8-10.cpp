#include <algorithm>
#include <iostream>
#include <random>
#include <cstdlib>
#include <set>
#include <list>
#include <vector>
#include <assert.h>

using namespace std;

class SetPartition
{
public:
    multiset<int> S1;
    multiset<int> S2;
    typedef vector<int>::iterator T;
    
    SetPartition(T first1, T last1, T first2, T last2) {
	S1.insert(first1, last1);
	S2.insert(first2, last2);
    }
    SetPartition(T first1, T last1) {
	S1.insert(first1, last1);
    }

    int get_diff(SetPartition &sp) {
	return (accumulate(S1.begin(), S1.end(), 0) -
		accumulate(S2.begin(), S2.end(), 0));
    };

    multiset<int> &get_set(bool small_one) {
	int a = accumulate(S1.begin(), S1.end(), 0);
	int b = accumulate(S2.begin(), S2.end(), 0);
	if (small_one && a <= b)
	    return S1;
	else if (small_one && a > b)
	    return S2;
	else if (!small_one && a <= b)
	    return S1;
	else
	    return S2;
    }

    void print(void) {
	cout << "Set 1:" << endl;
	for_each(S1.begin(), S1.end(),
		 [](int x) {
		     cout << x << ", ";
		 });
	cout << endl;
	for_each(S2.begin(), S2.end(),
		 [](int x) {
		     cout << x << ", ";
		 });
	cout << endl;
    }
};

bool operator< (SetPartition &P, SetPartition &Q)
{
    return
	lexicographical_compare(
	    P.get_set(true).begin(),
	    P.get_set(true).end(),
	    Q.get_set(true).begin(),
	    Q.get_set(true).end());
}

typedef set<SetPartition *> setpar;

template<class T>
void find_partitions(T first, T last, setpar &partitions, int diff)
{
    assert(diff >= 0);
    if (first == last) return;
    if (distance(first, last) == 1) {
	if (*first == diff) {
	    SetPartition *sp =
		new SetPartition(first, last);
	    partitions.insert(sp);
	}
	return;
    }
    if (distance(first, last) == 2) {
	if (abs(*first - *next(first)) == diff) {
	    SetPartition *sp =
		new SetPartition(first, next(first),
				 next(first), last);
	    partitions.insert(sp);
	}
	return;
    }
    // when the set contains more than 2 elem.
    setpar P1, P2;
    int d = diff;
    int z = *prev(last);
    if (z <= -diff) {
	find_partitions(first, prev(last), P1, -z-d);
	for_each(P1.begin(), P1.end(),
		 [z](SetPartition *P) {
		     P->get_set(false).insert(z);
		 });
	find_partitions(first, prev(last), P2, d-z);
	for_each(P2.begin(), P2.end(),
		 [z](SetPartition *P) {
		     P->get_set(false).insert(z);
		 });
    } else if (z > -d && z < d) {
	find_partitions(first, prev(last), P1, d-z);
	for_each(P1.begin(), P1.end(),
		 [z](SetPartition *P) {
		     P->get_set(false).insert(z);
		 });
	find_partitions(first, prev(last), P2, d+z);
	for_each(P2.begin(), P2.end(),
		 [z](SetPartition *P) {
		     P->get_set(true).insert(z);
		 });
    } else {// z >= d
	find_partitions(first, prev(last), P1, d+z);
	for_each(P1.begin(), P1.end(),
		 [z](SetPartition *P) {
		     P->get_set(true).insert(z);
		 });
	find_partitions(first, prev(last), P2, z - d);
	for_each(P2.begin(), P2.end(),
		 [z](SetPartition *P) {
		     P->get_set(true).insert(z);
		 });
    }
    partitions.insert(P1.begin(), P1.end());
    partitions.insert(P2.begin(), P2.end());
}

int main(int argc, char *argv[])
{
    // vector<int> numbers(stoi(argv[1]));
    vector<int> numbers({1, 5, 8, 11, 9});
    setpar partitions;
    srand(stoi(argv[2]));

    // generate(numbers.begin(), numbers.end(),
    // 	     [](void) {
    // 	     	 return rand() % 10;
    // 	     });
    for_each(numbers.begin(), numbers.end(),
	     [](int x) {
		 cout << x << ", ";
	     });
    cout << endl;
    find_partitions(numbers.begin(), numbers.end(), partitions, 0);
    auto i = partitions.begin();
    unsigned c = 0;
    while (distance(i, partitions.end()) > 0) {
	cout << "Partition " << c++ << endl;
	(*i)->print();
	delete *i;
	i = next(i);
    }
}
