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
    for_each(V.begin(), V.end(),
    	     [](item &i) {
    		 cout << "(" << i.number << ", " << i.color << ")" << endl;
    	     });
    V.sort([](item &a, item &b) {
	    return a.number < b.number;
	});
    // sort(V.begin(), V.end(),
    // 	 [](item &a, item &b) {
    // 	     return a.number < b.number;
    // 	 });
    auto iter = V.begin();
    auto pos = V.begin();
    unsigned color = 0;
    do {
	auto i =
	    find_if(iter, V.end(),
		    [=](item &i) {
			return i.color == color;
		    });
	if (i != V.end()) {
	    item c = *i;
	    iter = V.erase(i);
	    pos = V.insert(pos, c);
	    pos = next(pos);
	} else {
	    if (++color == 2) break;
	    iter = pos;
	}
    } while (1);
    // auto k = partition(V.begin(), V.end(),
    // 		       [](item &i) {
    // 			   return i.color == red;
    // 		       });
    // partition(k, V.end(),
    // 	      [](item &i) {
    // 		  return i.color == blue;
    // 	      });
    for_each(V.begin(), V.end(),
	     [](item &i) {
		 cout << "(" << i.number << ", " << i.color << ")" << endl;
	     });

}
