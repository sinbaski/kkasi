CC := g++
CFLAGS := -c -g3 -Wall -std=c++11
LDFLAGS := -lgsl -lgslcblas -lm

src := find_garch11_tail.cpp
obj := $(patsubst %.cpp, %.o, $(src))

$(obj): %.o: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

garch_tail: $(obj)
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: compile
compile: $(obj)

.PHONY: clean
clean:
	rm -f $(obj)
	rm garch_tail

.PHONY: all
all: garch_tail
