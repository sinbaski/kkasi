CC := g++
srcdir := .
objdir := obj

## srcs := $(wildcard $(srcdir)/*.cpp)
## srcs := ./test.cpp ./ExtremeNumber2.cpp ./XMatrix.cpp
srcs := ./GarchPQ.cpp
## srcs := garch21.cpp garch21_main.cpp
# srcs := garch1D.cpp main.cpp
## srcs := right_eigenfunction.cpp
objs := $(notdir $(patsubst %.cpp, %.o, $(srcs)))
objs := $(addprefix $(objdir)/, $(objs))

ifeq ($(debug),1)
	CPPFLAGS := -std=c++11 -c -Wall -g3 -D DEBUG
	CPPFLAGS += -fopenmp
	LDFLAGS := -lm -larmadillo -lgsl -lgslcblas
	LDFLAGS += -fopenmp
else
	CPPFLAGS := -std=c++11 -c -Wall -Wno-unused-variable -O3
	CPPFLAGS += -fopenmp
	LDFLAGS := -lm -larmadillo -lgsl -lgslcblas -fopenmp	
endif

ifeq ($(main), 1)
	CPPFLAGS += -D MAIN_ALGO
endif


$(objdir)/%.o: $(srcdir)/%.cpp
	@if [ ! -e $(objdir) ]; then mkdir $(objdir); fi
	$(CC) $(CPPFLAGS) $^ -o $@

.PHONY: garch1d
garch1d: $(objdir)/main.o $(objdir)/garch1D.o
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: garchpq
garchpq: obj/GarchPQ.o
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: garch21
garch21: $(objdir)/GarchPQ.o $(objdir)/garch21.o $(objdir)/XMatrix.o $(objdir)/ExtremeNumber2.o
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: MyTailEst
MyTailEst: $(objs)
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: right_eigenfunction
right_eigenfunction: $(objs)
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: all
all: right_eigenfunction

.PHONY: clean
clean:
	$(RM) $(objdir)/*.o



