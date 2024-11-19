# ModMolecular Dynamics Simulation (OMP Parallelized Version)

This version is a simple implement the Molecular Dynamics Simulation in OMP compared to the sequential version.
We want compared this wtih cuda to see the performance difference.

## How to run

```bash
# automatic run
$ make all

# or you can only run serial or omp version
$ make run-serial
$ make run-omp

# or you can run the prog directly, adjest the size
$ ./ljp config.in [size] [0: cpu, 1: omp] [0: not output step file, 1: output]

# to validate the output
$ make validate
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

## Output

All the ouput files are stored in the `output` directory.

Structures:
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