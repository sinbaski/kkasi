CC := g++
srcdir := .
objdir := obj

ifeq ($(debug),1)
	CFLAGS := -std=c++11 -c -Wall -Wno-unused-variable
	CFLAGS += -g3 -fopenmp -D DEBUG
	CFLAGS += -fpermissive 
	LDFLAGS := -lm -lgsl -lgslcblas -fopenmp -larmadillo
	LDFLAGS += -L/usr/local/lib -lnlopt
else
	CFLAGS := -std=c++11 -c -Wall -O3 -fopenmp -larmadillo
	CFLAGS += -fpermissive 
	LDFLAGS := -lm -lgsl -lgslcblas -fopenmp
	LDFLAGS += -L/usr/local/lib -lnlopt
endif

$(objdir)/%.o :%.cpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: prepare
prepare:
	@if [ ! -e $(objdir) ]; then mkdir $(objdir); fi

.PHONY: igarch
igarch: $(objdir)/igarch.o
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: arch
arch: $(objdir)/arch.o
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: implied_tail_index
implied_tail_index: $(objdir)/implied_tail_index.o
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: all
all: prepare igarch

.PHONY: clean
clean:
	$(RM) $(objdir)/*.o
	$(RM) igarch


