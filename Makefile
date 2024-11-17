build:
	@g++ -o ljp_omp main.cpp -std=c++11 -fopenmp -O3

clean:
	@echo "Cleaning..."
	@rm -f ljp*

run:
	@chmod +x run.sh
	@./run.sh

diff-in:
	@diff -w m.in mols.in

all: build run

.PHONY: build clean run diff-in all