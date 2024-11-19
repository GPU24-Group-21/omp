build:
	@g++ -o ljp main.cpp -std=c++17 -fopenmp -O3
	@g++ -o validator validator.cpp -std=c++17

clean:
	@echo "Cleaning..."
	@rm -f ljp* validator
	@rm -rf output

run:
	@chmod +x run.sh
	@./run.sh

run-output:
	@chmod +x run.sh
	@./run.sh -v

run-serial:
	@chmod +x run.sh
	@./run.sh -c

run-omp:
	@chmod +x run.sh
	@./run.sh -o

validate:
	@chmod +x validate.sh
	@./validate.sh

all: build run

.PHONY: build clean run all validate run-serial run-omp run-output