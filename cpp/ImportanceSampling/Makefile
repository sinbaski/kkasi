CC := g++
srcdir := .
objdir := obj

## srcs := $(wildcard $(srcdir)/*.cpp)
## srcs := ./test.cpp ./ExtremeNumber2.cpp ./XMatrix.cpp
## srcs := ./GarchPQ.cpp
## srcs := garch21.cpp garch21_main.cpp
# srcs := garch1D.cpp main.cpp
# objs := $(notdir $(patsubst %.cpp, %.o, $(srcs)))
# objs := $(addprefix $(objdir)/, $(objs))

ifeq ($(debug),1)
	CFLAGS := -std=c++11 -c -Wall -g3 -D DEBUG -D GARCH_PQ_ALONE
	LDFLAGS := -lm -larmadillo -lgsl -lgslcblas
else
	CFLAGS := -std=c++11 -c -Wall -Wno-unused-variable -O3 -fopenmp -D GARCH_PQ_ALONE
	LDFLAGS := -lm -larmadillo -lgsl -lgslcblas -fopenmp	
endif

# $(objs): $(objdir)/%.o: $(srcdir)/%.cpp
# 	@if [ ! -e $(objdir) ]; then mkdir $(objdir); fi
# 	$(CC) $(CFLAGS) $^ -o $@


$(objdir)/%.o :%.cpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: garch1d
garch1d: $(objdir)/main.o $(objdir)/garch1D.o
	$(CC) $(LDFLAGS) $^ -o $@

.PHONY: garchpq
garchpq: $(objdir)/GarchPQ.o
	$(CC) $(LDFLAGS) $^ -o $@

# .PHONY: main
# main: $(objs)
# ##	$(CC) $(filter %.o, $^) $(LDFLAGS) -o $@
# 	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: all
all: garch1d garchpq

.PHONY: clean
clean:
	$(RM) $(objdir)/*.o
	$(RM) garch1d garchpq


