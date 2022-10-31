
main: main.o circuit.o circuit_interface.o
	g++ $^ -o $@

main.o: main.cpp circuit.h
	g++ -I lib/eigen -c main.cpp

circuit.o: circuit.cpp circuit.h circuit_interface.h
	g++ -I lib/eigen -c circuit.cpp

circuit_interface.o: circuit_interface.cpp circuit.h circuit_interface.h
	g++ -I lib/eigen -I lib/json/single_include -c circuit_interface.cpp

clean:
	rm -rf *.o main out.json
