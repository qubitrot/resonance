CXX=clang++

SRCDIR = ./src
OBJDIR = ./obj

MAGMADIR ?= /opt/magma
CUDADIR  ?= /opt/cuda

SRC     = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/json/*.cpp) $(wildcard $(SRCDIR)/profiler/*.cpp)
OBJS    = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))
CFLAGS  = -O3 -std=c++11 -Wall -Wextra -DNDEBUG
LDFLAGS = -pthread -llapack 
INCLUDE = -I/usr/local/include -I./src -I$(MAGMADIR)/include -I$(CUDADIR)/include -I$(SRCDIR)/profiler

res:$(OBJS)
	$(CXX) -o $@ ${OBJS} $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(INCLUDE) $(CFLAGS) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/json/%.cpp
	$(CXX) -c $(INCLUDE) $(CFLAGS) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/profiler/%.cpp
	$(CXX) -c $(INCLUDE) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f obj/*.o res
