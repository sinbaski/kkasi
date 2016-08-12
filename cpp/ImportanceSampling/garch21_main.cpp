#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <stdio.h>
#include "garch21.hpp"

size_t num_iter;

int main(int argc, char **argv)
{
    array<double, 4> params;
    size_t lines[2];

    transform(argv+1, argv + params.size() + 1, params.begin(),
	      [](char *str)
	      {
	      	  return stod(str);
	      });
    num_iter = stoi(argv[params.size() + 1]);
    transform(argv + params.size() + 2, argv + argc, lines,
	      [](char *str)
	      {
		  return stoul(str);
	      });
    garch21 model(params);

    cout << "Tail index = " << model.compute_tail_index(lines[0], lines[1]) << endl;
    // random_device dev;
    // for (size_t i = 0; i < 4000; i++) {
    // 	array<double, 2> r = model.forward(dev, 0.995, false);
    // 	printf("%e\t%e\n", r[0], r[1]);
    // }
}
