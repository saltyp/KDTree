##
## Makefile
##

CC = g++ -std=c++17 # compiler to use

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CXXFLAGS = -g -Wall -Wextra #compile flags
# LDFLAGS = -g # flags for linker (eg -lm to link to math library ); for debugging, -g has to be used both in compilation & linkage

# ----------------------------------
# Files that are part of the project
# ----------------------------------
PROG = a.out #executable
MAINAPP = test-harness.cpp
# CLASSES = KDTree.h Point.h BoundedPQueue.h #dependencies so that make reruns if any of these change
MAINOBJ = $(MAINAPP:.cpp=.o)

# ----------------------------------
#linkage
$(PROG) : $(MAINOBJ)
	$(CC) $(CXXFLAGS) -o $@ $^

#compilation
$(MAINOBJ) : $(MAINAPP)
	$(CC) $(CXXFLAGS) -c $^

#dependency
# $(MAINAPP) : $(CLASSES)

# ----------------------------------
clean :
	rm -i *.out *.o *.gch
