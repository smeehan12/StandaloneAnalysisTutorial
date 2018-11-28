SOURCES=NTupleMaker.cc 
EXECUTABLE=NTuple
OBJECTS=$(SOURCES:.cpp=.o)
FASTJET=/Users/meehan/work/fastjet-install/

CC=g++
CFLAGS=-c -I$(ROOTSYS)/include 
LDFLAGS=-L$(ROOTSYS)/lib $(shell root-config --glibs) $(shell root-config --cflags) `$(FASTJET)/bin/fastjet-config --cxxflags --libs --plugins`

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *~
