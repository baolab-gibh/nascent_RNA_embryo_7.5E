CXX          = g++
CXXFLAGS     = -g -Wall -std=c++11
INCLUDE_DIR  = -I/usr/include -I./src/htslib -I./src
LIB_DIR      = -L/usr/lib -L/usr/lib/x86_64-linux-gnu

LIBS         = -lz -ldeflate -llzma -lbz2 -lm -lcurl -lcrypto -lpthread
HTSLIB       = -L./src/htslib
EXECUTABLE   = slamtk

# OBJECTS      = $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) build/libhts.a
OBJECTS      := build/BamReader.o build/SamRecord.o build/CmdLineOpts.o build/SamAuxTag.o build/DNASequence.o build/FastqQuality.o build/CigarString.o build/libhts.a


all: $(EXECUTABLE)


$(EXECUTABLE): main.cpp $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(INCLUDE_DIR) $(HTSLIB) $(LIBS)

build/%.o: src/%.cpp
	$(CXX) -o $@ $< -c ${CXXFLAGS} ${INCLUDE_DIR}

build/libhts.a:
	mkdir -p build
	cd src/htslib && ./configure && make lib-static
	cp src/htslib/libhts.a build/libhts.a

test:
	./$(EXECUTABLE) -n 3 data/test.bam

clean:
	rm -rf $(EXECUTABLE) build/*.o

clean-all:
	rm -rf $(EXECUTABLE) build
	cd src/htslib && make clean


.PHONY: clean all
