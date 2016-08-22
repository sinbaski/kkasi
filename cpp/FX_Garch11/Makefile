CC := g++
CFLAGS := -c -g3 -Wall -std=c++11 -fopenmp
LDFLAGS := -lgsl -lgslcblas -lm -fopenmp

src := IARCH11.cpp
obj := $(patsubst %.cpp, %.o, $(src))

$(obj): %.o: %.cpp
	$(CC) $(CFLAGS) $^ -o $@

main: $(obj)
	$(CC) $^ $(LDFLAGS) -o $@

.PHONY: compile
compile: $(obj)

.PHONY: clean
clean:
	rm -f $(obj)
	rm main

.PHONY: all
all: main
