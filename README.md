# ModMolecular Dynamics Simulation (OMP Parallelized Version)

This version is a simple implement the Molecular Dynamics Simulation in OMP compared to the sequential version.
We want compared this wtih cuda to see the performance difference.

## How to run

```bash
# build the program
make

# automatic build and run (will run the faster mode, no output step file)
make all

# or you can only run serial or omp version
make run-serial
make run-omp

# if you want to output the step file, you can run the following command to run both serial and omp version
make run-output # run both serial and omp version with output step file
make run-omp-output # run only omp version with output step file
make run-serial-output # run only serial version with output step file

# or you can run the prog directly, adjest the size
./ljp config.in [size] [0: cpu, 1: omp] [0: not output step file, 1: output (optional)]

# to validate the output
make validate

# to clean the output
make clean
```

## Configuration

The configuration file is `config.in` which contains the following parameters:

```
deltaT  0.005    // time step
density	0.8      // density
stepAvg 10       // average steps to output the summary report, reduce this will speed up the simulation
stepLimit   100  // total steps to loop
temperature 1.0  // temperature
```

The `run.sh` have a series of test cases to run the program with different size and different version.
We currenlty set the max size to 400 which is 160000 atoms, in cpu version, it will take a long time to run(hour unit).
You can change the size in the `run.sh` to test the performance of the program.

or you can run the program directly with the following command:

```bash
./ljp config.in 100 0 // run the cpu version with 100 atoms and no output step file

./ljp config.in 100 1 // run the omp version with 100 atoms and output step file

./ljp config.in 100 0 1 // run the cpu version with 100 atoms and output step file

./ljp config.in 100 1 1 // run the omp version with 100 atoms and output step file
```

## Output

All the ouput files are stored in the `output` directory. Each will group by cpu or omp version and the size of the system.
Each output file will contain the position of each atom in each step and the energy, temperature, pressure, and timing of each step.

Output folder Structures:
```
output
|── cpu
 |── <size>
   |── <step>.out  // output file for each step including position of each atom
   ...
   └── final // summary of the steps based on the stepAvg and timing
└── omp
 |── <size>
   |── <step>.out // output file for each step including position of each atom
   ...
   └── final // summary of the steps based on the stepAvg and timing
```