CC := g++
srcdir := .
objdir := obj

srcs := $(wildcard $(srcdir)/*.cpp)
objs := $(notdir $(patsubst %.cpp, %.o, $(srcs)))
objs := $(addprefix $(objdir)/, $(objs))

## debug := 0

ifeq ($(debug),1)
	CFLAGS := -std=c++11 -c -Wall -Werror -g3 -Og
else
	CFLAGS := -std=c++11 -c -Wall -Werror -O2
endif

LDFLAGS := -lm -larmadillo

$(objs): $(objdir)/%.o: $(srcdir)/%.cpp
	@if [ ! -e $(objdir) ]; then mkdir $(objdir); fi
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: main
main: $(objs)
	$(CC) $(filter %.o, $^) $(LDFLAGS) -o $@

.PHONY: all
all: main

.PHONY: clean
clean:
	$(RM) $(objdir)/*.o
	$(RM) main
