#By Yohanna White 1/4/18
#Edited by Davis Templeton 2/5/2018

#g++ -std=c++0x -I$GEMSHOME/gmml/includes/* -L$GEMSHOME/gmml//bin/ -Wl,-rpath,$GEMSHOME/gmml/bin/ *.cpp -lgmml -o gp_builder

CC = g++
CFLAGS = -std=gnu++11 -I${GEMSHOME}/gmml/includes
RFLAGS = -Wl,-rpath,${GEMSHOME}/gmml/bin/
LFLAGS = -I${GEMSHOME}/gmml/includes -L${GEMSHOME}/gmml/bin -lgmml
COMPILE = $(CC) $(CFLAGS) -c
LINK = $(LFLAGS)
RUNTIME = $(RFLAGS)
#DEBUG = -g 2
RM = rm -f

SRC = ./src
INC = ./includes
BUILD = ./build
BIN = ./bin
all: $(BIN)/main $(BUILD)/main.o $(BUILD)/aromaticResidue.o $(BUILD)/chPiPair.o 

run: $(BIN)/main
	$(BIN)/main

# Notice the change in the order of the $(CC) and the $(LINK) variables in the command.
$(BIN)/main: $(BUILD)/main.o $(BUILD)/aromaticResidue.o $(BUILD)/chPiPair.o
	$(CC) $(BUILD)/main.o $(BUILD)/aromaticResidue.o $(BUILD)/chPiPair.o $(LINK) $(RUNTIME) -o $(BIN)/main

$(BUILD)/main.o: $(SRC)/main.cpp
	$(COMPILE) $(SRC)/main.cpp -o $(BUILD)/main.o

$(BUILD)/aromaticResidue.o: $(SRC)/aromaticResidue.cpp $(INC)/aromaticResidue.hpp
	$(COMPILE) $(SRC)/aromaticResidue.cpp -o $(BUILD)/aromaticResidue.o

$(BUILD)/chPiPair.o: $(SRC)/chPiPair.cpp $(INC)/chPiPair.hpp
	$(COMPILE) $(SRC)/chPiPair.cpp -o $(BUILD)/chPiPair.o

clean:
	rm $(BUILD)/*.o
	rm $(BIN)/main

# Doesn't work
# g++ -I/home/korin/WORKSPACE/CCRC/public/gems/gmml/includes \
# 	-L/home/korin/WORKSPACE/CCRC/public/gems/gmml/bin -lgmml \
# 	./build/main.o ./build/aromatic.o ./build/aliphaticH.o ./build/chPiPair.o ./build/aromAliphContainer.o \
# 	-o ./bin/main
#
# WORKS!
# g++ ./build/main.o ./build/aromatic.o ./build/chPiPair.o ./build/aliphaticH.o \
# 	-I/home/korin/WORKSPACE/CCRC/public/gems/gmml/includes -L/home/korin/WORKSPACE/CCRC/public/gems/gmml/bin -lgmml \
# 	-o main
