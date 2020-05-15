SRC    := ../src/metropolis.hpp ../src/rng.hpp ../src/histogram.hpp ../src/radiation.hpp ../src/proposal_density.hpp
FLAGS  := -std=c++17 -fconcepts -Wall -Wpedantic -fdiagnostics-color=always -I ../range-v3/include
FLAGS3 := -O3 -ffast-math $(FLAGS)

.PHONY: all clean

all: bks_dipole_fig

dipole: %: %.cpp $(SRC)
	g++ $(FLAGS) $< -o $@ && ./$@

clean:
	rm -I bks_electron_data/*.txt bks_electron bks_electron.pdf
