#include <algorithm>
#include <random>
#include <string>
#include "garch21.hpp"

int main(int argc, char **argv)
{
    array<double, 4> params;
    random_device dev;

    transform(argv+1, argv + argc, params.begin(),
	      [](char *str)
	      {
	      	  return stod(str);
	      });
    garch21 model(params);
    
    model.compute_tail_index();
}
