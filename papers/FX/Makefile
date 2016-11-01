CC := g++
srcdir := .
objdir := obj

ifeq ($(debug),1)
	CFLAGS := -std=c++11 -c -Wall -Wno-unused-variable -g3 -fopenmp -D DEBUG
	LDFLAGS := -lm -lgsl -lgslcblas -fopenmp
else
	CFLAGS := -std=c++11 -c -Wall -O3 -fopenmp
	LDFLAGS := -lm -lgsl -lgslcblas -fopenmp	
endif

$(objdir)/%.o :%.cpp
	$(CC) $(CFLAGS) $< -o $@

.PHONY: prepare
prepare:
	@if [ ! -e $(objdir) ]; then mkdir $(objdir); fi

.PHONY: igarch
igarch: $(objdir)/igarch.o
	$(CC) $(LDFLAGS) $^ -o $@

.PHONY: arch
arch: $(objdir)/arch.o
	$(CC) $(LDFLAGS) $^ -o $@

.PHONY: all
all: prepare igarch

.PHONY: clean
clean:
	$(RM) $(objdir)/*.o
	$(RM) igarch


