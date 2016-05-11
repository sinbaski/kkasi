CC := g++
srcdir := .
objdir := obj

## srcs := $(wildcard $(srcdir)/*.cpp)
srcs := ./test.cpp ./ExtremeNumber.cpp
objs := $(notdir $(patsubst %.cpp, %.o, $(srcs)))
objs := $(addprefix $(objdir)/, $(objs))

ifeq ($(debug),1)
	CFLAGS := -std=c++11 -c -Wall -Wno-unused-variable -Werror -g3 -fopenmp
else
	CFLAGS := -std=c++11 -c -Wall -Werror -O3 -fopenmp
endif

LDFLAGS := -lm -larmadillo -lgsl -lgslcblas -fopenmp

$(objs): $(objdir)/%.o: $(srcdir)/%.cpp
	@if [ ! -e $(objdir) ]; then mkdir $(objdir); fi
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: main
main: $(objs)
##	$(CC) $(filter %.o, $^) $(LDFLAGS) -o $@
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: all
all: main

.PHONY: clean
clean:
	$(RM) $(objdir)/*.o
	$(RM) main
