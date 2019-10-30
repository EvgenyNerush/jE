SRC    := ../src/metropolis.hpp ../src/rng.hpp ../src/histogram.hpp ../src/radiation.hpp
BIN    := assert metropolis synchrotron_radiation
FLAGS  := -std=c++17 -Wall -Wpedantic -I ../range-v3/include
FLAGS3 := -O3 -ffast-math $(FLAGS)

all: $(BIN)

assert: %: %.cpp $(SRC)
	g++ -g $(FLAGS) $< -o $@ && ./$@

metropolis: %: %.cpp $(SRC)
	g++ $(FLAGS3) $< -o $@ && ./$@

synchrotron_radiation: %: %.cpp $(SRC)
	g++ $(FLAGS3) $< -o $@ && ./$@

clean:
	rm -I $(BIN)
