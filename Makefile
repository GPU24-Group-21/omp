build:
	@g++ -o ljp main.cpp -std=c++11 -fopenmp -O3
	@g++ -o validator validator.cpp -std=c++11

clean:
	@echo "Cleaning..."
	@rm -f ljp* validator
	@rm -rf output

run:
	@chmod +x run.sh
	@./run.sh

run-serial:
	@chmod +x run.sh
	@./run.sh -c

run-omp:
	@chmod +x run.sh
	@./run.sh -o

validate:
	@chmod +x validate.sh
	@./validate.sh

all: build run validate

.PHONY: build clean run all validate run-serial run-omp