SRC := ../src/metropolis.hpp ../src/rng.hpp ../src/histogram.hpp

all: assert metropolis

assert: %: %.cpp $(SRC)
	g++ -g -Wall -Wpedantic assert.cpp -o assert && ./assert

metropolis: %: %.cpp $(SRC)
	g++ -O3 -ffast-math -Wall -Wpedantic metropolis.cpp -o metropolis && ./metropolis

clean:
	rm -I assert metropolis
