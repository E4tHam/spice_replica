
LIBS = -L$(MATLABROOT)/extern/bin/glnxa64 -L$(MATLABROOT)/bin/glnxa64
INCLUDES = -Isrc -I$(MATLABROOT)/extern/include -Ilib/eigen -Ilib/json/single_include
FLAGS = -std=c++11 $(LIBS) $(INCLUDES)

HEADERS = src/analysis.h src/circuit.h src/dc.h src/helper.h src/matlab.h src/tran.h

main: main.o analysis.o circuit.o dc.o helper.o matlab.o tran.o
	g++ $(FLAGS) $^ -leng -lmx -o $@

main.o: main.cpp $(HEADERS)
	g++ $(FLAGS) -c $<

%.o: src/%.cpp $(HEADERS)
	g++ $(FLAGS) -c $<

clean:
	rm -rf *.o main out.json
