#include <list>
#include <algorithm>
#include <iostream>
#include <random>

using namespace std;

random_device gen;

// enum color {
//     red, blue, yellow
// };

class item {
public:
    double number;
    unsigned color;
};

int main(void)
{
    list<item> V(5);
    for_each(V.begin(), V.end(),
	     [](item &i) {
		 i.number = gen();
		 i.color = gen() % 3;
	     });
    V.sort([](item &a, item &b) {
	    return a.number < b.number;
	});
    for_each(V.begin(), V.end(),
    	     [](item &i) {
    		 cout << "(" << i.number << ", " << i.color << ")" << endl;
    	     });
    auto iter = V.begin();
    auto pos = V.begin();
    unsigned color = 0;
    do {
	auto i =
	    find_if(iter, V.end(),
		    [=](item &k) {
			return k.color == color;
		    });
	if (i == pos) {
	    iter = pos = next(pos);
	} else if (i != V.end()) {
	    item c = *i;
	    iter = V.erase(i);
	    pos = next(V.insert(pos, c));
	} else {
	    if (++color == 2) break;
	    iter = pos;
	}
    } while (1);
    cout << endl;
    for_each(V.begin(), V.end(),
	     [](item &i) {
		 cout << "(" << i.number << ", " << i.color << ")" << endl;
	     });

}
