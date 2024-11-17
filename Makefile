build:
	@g++ -o ljp main.cpp -std=c++11 -fopenmp -O3
	@g++ -o validator validator.cpp -std=c++11

clean:
	@echo "Cleaning..."
	@rm -f ljp*

run:
	@chmod +x run.sh
	@./run.sh

diff-in:
	@diff -w m.in mols.in

validate:
	@chmod +x validate.sh
	@./validate.sh

all: build run validate

.PHONY: build clean run diff-in all validate