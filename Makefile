#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#

# Executables to build
EXEC += viewer
EXEC += subgraph
EXEC += shortest_path
EXEC += mass_spring
EXEC += test_nodes
EXEC += test_edges
EXEC += poisson
EXEC += mtl_test
EXEC += space_search_test

# Get the shell name to determine the OS
UNAME := $(shell uname)

# Define the C++ compiler to use

ifeq ($(UNAME), Linux)
  CXX := g++ -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP 
endif
ifeq ($(UNAME), Darwin)
  CXX := g++ -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP 
endif

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) -MP

# Define any directories containing header files
#   To include directories use -Ipath/to/files
INCLUDES += -I.
INCLUDES += -I./MTL
INCLUDES += -I./CME212
INCLUDES += -I./thrust
 
# Define CXX compile flags
CXXFLAGS += -std=c++11 -O3 -funroll-loops -W -Wall -Wextra #-Wfatal-errors
#CXXFLAGS += -DTHRUST_DEVICE_SYSTEM = THRUST_DEVICE_SYSTEM_OMP

# Define any directories containing libraries
#   To include directories use -Lpath/to/files
LDFLAGS +=

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
ifeq ($(UNAME), Linux)
  LDLIBS += -lSDL -lGL -lGLU
endif
ifeq ($(UNAME), Darwin)
  LDLIBS += -L/usr/local/lib -lSDLmain -lSDL -Wl,-framework,Cocoa,-framework,OpenGL
endif

##################
# The following part of the makefile defines generic rules; it can be used to
# build any executable just by changing the definitions above.
#
#   $^: the name of the prereqs of the rule
#   $<: the name of the first prereq of the rule
#   $@: the name of the target of the rule
##################

# 'make' - default rule
all: $(EXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): % : %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -c -o $@ $<

# Extra dependencies for executables
#   Nothing here

# 'make clean' - deletes all .o files, exec, and dependency files
clean:
	-$(RM) *.o $(EXEC)
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)
