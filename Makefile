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

run-serial-output:
	@chmod +x run.sh
	@./run.sh -c -v

run-omp-output:
	@chmod +x run.sh
	@./run.sh -o -v

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
all-output: build run-output

help:
	@echo "Usage: make [target]"
	@echo "Targets:"
	@echo "  build: Compile the source code"
	@echo "  clean: Remove the compiled files"
	@echo "  run: Run the program"
	@echo "  run-output: Run the program and output the result"
	@echo "  run-serial: Run the program in serial mode"
	@echo "  run-omp: Run the program in OpenMP mode"
	@echo "  run-serial-output: Run the program in serial mode and output the result"
	@echo "  run-omp-output: Run the program in OpenMP mode and output the result"
	@echo "  validate: Validate the output"
	@echo "  all: Compile and run the program"
	@echo "  all-output: Compile and run the program and output the result"

.PHONY: build clean run all validate run-serial run-omp run-output run-serial-output run-omp-output all-output help