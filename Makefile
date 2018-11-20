# Declaration of variables
CC = g++
CC_FLAGS = -std=c++0x -I ${GEMSHOME}/gmml/includes/ -I includes/
RFLAGS = -Wl,-rpath,${GEMSHOME}/gmml/lib/
LFLAGS = -I ${GEMSHOME}/gmml/includes/ -L ${GEMSHOME}/gmml/lib/ -lgmml
 
# File names
BUILD = build/
EXEC = bin/gp_builder
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o $(EXEC)
 
# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $(BUILD)/$@
 
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
