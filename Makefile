CC=gcc
CXX=g++

SRCDIR = ./src
OBJDIR = ./obj

SRC     = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/json/*.cpp)
OBJS    = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.cpp=.o)))
CFLAGS  = -O3 -std=c++11 -Wall -Wextra -pedantic
LDFLAGS = -pthread -llapack
INCLUDE = -I/usr/local/include -I./src

res:$(OBJS)
	$(CXX) -o $@ ${OBJS} $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(INCLUDE) $(CFLAGS) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/json/%.cpp
	$(CXX) -c $(INCLUDE) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f obj/*.o res
