#include <algorithm>
#include <random>
#include <string>
#include "garch21.hpp"

size_t num_iter;

int main(int argc, char **argv)
{
    array<double, 4> params;
    random_device dev;

    transform(argv+1, argv + params.size() + 1, params.begin(),
	      [](char *str)
	      {
	      	  return stod(str);
	      });
    garch21 model(params);
    num_iter = stoi(argv[argc - 1]);
    model.compute_tail_index();
}
