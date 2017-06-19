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
	    return S2;
	else
	    return S1;
    }

    void print(void) {
	printf("(");
	for_each(S1.begin(), S1.end(),
		 [](int x) {
		     printf("%d, ", x);
		 });
	printf(") + (");
	for_each(S2.begin(), S2.end(),
		 [](int x) {
		     printf("%d, ", x);
		 });
	printf(")\n");
    }
};

class comp_ppt
{
public:
    bool operator() (SetPartition *P, SetPartition *Q) {
	return
	    lexicographical_compare(
		P->get_set(true).begin(),
		P->get_set(true).end(),
		Q->get_set(true).begin(),
		Q->get_set(true).end());
    }
};

typedef set<SetPartition *, comp_ppt> setpar;

void merge_setpar(setpar &S, setpar &P)
{
    for ( auto i = P.begin(); i != P.end(); advance(i, 1))
	if (!S.insert(*i).second) delete *i;
}

template<class T>
void find_partitions(const T first, const T last,
		     setpar &partitions, int diff)
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
	if (abs(*first + *next(first)) == diff) {
	    SetPartition *sp =
		new SetPartition(first, last);
	    partitions.insert(sp);
	}
	return;
    }
    // when the set contains more than 2 elem.
    int d = diff;
    int z = *prev(last);
    auto iter = prev(last);
    // printf("z = %d, d = %d:\n", z, d);
    // printf("From previous round:\n");
    if (d >= z) {
	setpar P;
	find_partitions(first, iter, P, d - z);
	printf("z=%d, d=%d:\n", z, d - z);
	for_each(P.begin(), P.end(),
		 [](SetPartition *sp) {
		     sp->print();
		 });
	for(auto i = P.begin(); i != P.end(); i++) {
	    (*i)->get_set(false).insert(z);
	}
	merge_setpar(partitions, P);
    } else {
	setpar P;
	find_partitions(first, iter, P, z - d);
	printf("z=%d, d=%d:\n", z, z - d);
	for_each(P.begin(), P.end(),
		 [](SetPartition *sp) {
		     sp->print();
		 });
	for(auto i = P.begin(); i != P.end(); i++) {
	    (*i)->get_set(true).insert(z);
	}
	merge_setpar(partitions, P);
    }
    if (d + z >= 0) {
	setpar P;
	find_partitions(first, iter, P, d + z);
	printf("z=%d, d=%d:\n", z, d + z);
	for_each(P.begin(), P.end(),
		 [](SetPartition *sp) {
		     sp->print();
		 });
	for(auto i = P.begin(); i != P.end(); i++) {
	    (*i)->get_set(true).insert(z);
	}
	merge_setpar(partitions, P);
    } else {
	setpar P;
	find_partitions(first, iter, P, -d - z);
	printf("z=%d, d=%d:\n", z, -d - z);
	for_each(P.begin(), P.end(),
		 [](SetPartition *sp) {
		     sp->print();
		 });
	for(auto i = P.begin(); i != P.end(); i++) {
	    (*i)->get_set(false).insert(z);
	}
	merge_setpar(partitions, P);
    }
    // printf("current candidates:\n");
    // for_each(partitions.begin(), partitions.end(),
    // 	     [](SetPartition *sp) {
    // 		 sp->print();
    // 	     });
}

int main(int argc, char *argv[])
{
    // vector<int> numbers(stoi(argv[1]));
    // vector<int> numbers({1, 5, 8, 11, 9});
    vector<int> numbers({-4, 8, 4, -7, 7, 7, 12, 3});
    // vector<int> numbers({0, 7, 4, 9, 3, 3, 2, 5, 1});
    // vector<int> numbers({7, 4, 9, 3, 3, 2, 5, 1}); // sum to 17
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
    if (partitions.empty()) {
	cout << "No partition is found." << endl;
	return 0;
    }
    auto i = partitions.begin();
    unsigned c = 0;
    while (distance(i, partitions.end()) > 0) {
	cout << "Partition " << c++ << endl;
	(*i)->print();
	delete *i;
	advance(i, 1);
//	cout << endl;
    }
}
