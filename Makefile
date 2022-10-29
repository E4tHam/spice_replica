
main: main.o circuit.o circuit_loader.o
	g++ $^ -o $@

main.o: main.cpp circuit.h
	g++ -I lib/eigen -c main.cpp

circuit.o: circuit.cpp circuit.h circuit_loader.h
	g++ -I lib/eigen -c circuit.cpp

circuit_loader.o: circuit_loader.cpp circuit.h circuit_loader.h
	g++ -I lib/eigen -I lib/json/single_include/nlohmann -c circuit_loader.cpp

clean:
	rm -rf *.o main
