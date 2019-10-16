SRC := ../src/metropolis.hpp ../src/rng.hpp ../src/histogram.hpp ../src/synchrotron.hpp
BIN := assert metropolis synchrotron
FLAGS := -O3 -ffast-math -std=c++17 -Wall -Wpedantic 

all: $(BIN)

assert: %: %.cpp $(SRC)
	g++ -g -Wall -Wpedantic $< -o $@ && ./$@

metropolis: %: %.cpp $(SRC)
	g++ $(FLAGS) $< -o $@ && ./$@

synchrotron: %: %.cpp $(SRC)
	g++ $(FLAGS) $< -o $@ && ./$@

clean:
	rm -I $(BIN)
