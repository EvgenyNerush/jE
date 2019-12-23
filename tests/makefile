SRC    := ../src/metropolis.hpp ../src/rng.hpp ../src/histogram.hpp ../src/radiation.hpp ../src/proposal_density.hpp
BIN    := assert metropolis synchrotron_radiation proposal_density bks_low_chi bks_high_chi
FLAGS  := -std=c++17 -fconcepts -Wall -Wpedantic -fdiagnostics-color=always -I ../range-v3/include
FLAGS3 := -O3 -ffast-math $(FLAGS)

.PHONY: all clean

all: $(BIN)

assert: assert.cpp $(SRC)
	g++ -g $(FLAGS) $< -o $@ && ./$@

%: %.cpp $(SRC)
	g++ $(FLAGS3) $< -o $@ && ./$@
clean:
	rm -I $(BIN)
