#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

random_device gen;

struct bill;

struct bill
{
    unsigned name;
    unsigned id;

    class equal_name
    {
    public:
    	bool operator() (const bill &a, const bill &b) const {
    	    return a.name == b.name;
    	}
    };

    class hash
    {
    public:
	size_t operator() (const bill &b) const {
	    // unsigned n = sizeof(size_t) >> 1;
	    // return (b.name << n) | b.id;
	    return b.name;
	}
    };
    // static bool equal_name (const bill &a, const bill &b) {
    // 	return a.name == b.name;
    // }
    static bool smaller_id (const bill &a, const bill &b) {
	return a.id < b.id;
    }
    // static hash(const bill &b) {
    // 	unsigned n = sizeof(size_t) >> 1;
    // 	return (b.name << n) | b.id;
    // }
};

// class comp_bills
// {
// public:
//     bool operator() (const bill &a, const bill &b) const {
// 	return a.name == b.name;
//     }
// };

// class hash_bill
// {
// public:
//     size_t operator() (const bill &b) const {
// 	unsigned n = sizeof(size_t) >> 1;
// 	return (b.name << n) | b.id;
//     }
// };

typedef unordered_multiset<bill, bill::hash, bill::equal_name> myset;
// typedef unordered_multiset<bill, bill::hash> myset;

void find_didntpay(const myset &bills,
		   const myset &checks,
		   myset &didntpay) 
{
    vector<bill> U(bills.size());
    auto i = unique_copy(bills.begin(), bills.end(), U.begin(), bill::equal_name());
    U.resize(U.end() - i);
    for (auto i = U.begin(); i != U.end(); i = next(i)) {
	unsigned b = bills.count(*i);
	unsigned c = checks.count(*i);
	if (b > c) {
	    auto B = bills.equal_range(*i);
	    auto C = checks.equal_range(*i);
	    if (C.first == C.second) {
	    	didntpay.insert(C.first, C.second);
	    } else {
		vector<bill> VB(B.first, B.second);
		vector<bill> VC(C.first, C.second);
	    	sort(VB.begin(), VB.end(), bill::smaller_id);
	    	sort(VC.begin(), VC.end(), bill::smaller_id);
	    	vector<bill> V(b - c);
	    	set_difference(VC.begin(), VC.end(), VB.begin(), VB.end(), V.begin(), bill::smaller_id);
	    	didntpay.insert(V.begin(), V.end());
	    }
	}
    }
}

int main(void)
{
    unsigned n = 10, m=3;
    vector<bill> V(n);
    // vector<unsigned> V(n);
    for (unsigned i = 0; i < n; i++) {
    	V[i].name = gen() % m;
    	V[i].id = gen() % n;
	// V[i] = i;
    }
    myset checks(V.begin(), V.end());
    myset bills = checks;
    myset didntpay;
    bills.insert({0, 1});
	// bill({0, 2}), bill({1, 2}));
    find_didntpay(bills, checks, didntpay);
    for_each(didntpay.begin(), didntpay.end(),
    	     [](const bill &b) {
    		 cout << b.name << "did't pay his bill No." << b.id << endl;
    	     });
}
