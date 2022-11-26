
LIBS = -L$(MATLABROOT)/extern/bin/glnxa64 -L$(MATLABROOT)/bin/glnxa64
INCLUDES = -Isrc -I$(MATLABROOT)/extern/include -Ilib/eigen -Ilib/json/single_include
FLAGS = -std=c++11 $(LIBS) $(INCLUDES)

main: main.o circuit.o circuit_interface.o matlab.o
	g++ $(FLAGS) $^ -leng -lmx -o $@

main.o: main.cpp src/circuit.h
	g++ $(FLAGS) -c main.cpp

circuit.o: src/circuit.cpp src/circuit.h src/circuit_interface.h
	g++ $(FLAGS) -c src/circuit.cpp

circuit_interface.o: src/circuit_interface.cpp src/circuit.h src/circuit_interface.h
	g++ $(FLAGS) -c src/circuit_interface.cpp

matlab.o: src/matlab.cpp src/matlab.h
	g++ $(FLAGS) -c src/matlab.cpp

clean:
	rm -rf *.o main out.json
