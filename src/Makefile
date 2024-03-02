##
## Makefile
#  Should compile as : 
# g++ -std=c++17  -g -Wall -Wextra  -c evil-hangman.cpp
# g++ -g -Wall -Wextra    -c -o wordFamily.o wordFamily.cpp
# g++ -g -Wall -Wextra    -c -o hangmanGame.o hangmanGame.cpp
# g++ -std=c++17  -g -Wall -Wextra  -o playgame evil-hangman.o wordFamily.o hangmanGame.o
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
PROG = codetest #executable
MAINAPP = test-harness.cpp
MAINOBJ = $(MAINAPP:.cpp=.o)

# ----------------------------------
#linkage
$(PROG) : $(MAINOBJ)
	$(CC) $(CXXFLAGS) -o $@ $^

#compilation
$(MAINOBJ) : $(MAINAPP)
	$(CC) $(CXXFLAGS) -c $^

# ----------------------------------
clean :
	rm -i *.out *.o *.h.gch