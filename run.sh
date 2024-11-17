prog="./ljp_omp"
infile="Rap_2_LJP.in"

# find any .in file and set it as the input file
if [ -f *.in ]; then
    infile=$(ls *.in)
fi
if [ ! -d "output" ]; then
    mkdir output
fi

clear

# Run the program with the given arguments
echo "Running $prog with input file $infile\n"
$prog $infile 0 0