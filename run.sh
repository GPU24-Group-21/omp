#!/bin/bash
prog="./ljp"
infile="config.in"

# find any .in file and set it as the input file
if [ -f *.in ]; then
    infile=$(ls *.in)
fi

mkdir -p output/cpu
mkdir -p output/omp
rm -rf output/cpu/*
rm -rf output/omp/*

clear

# read the arguments -c for cpu, -o for omp if not, run both

series="10 20 40 80 100 200 400"

verbose=0

if [ "$1" == "-v" ]; then
    verbose=1
fi

if [ "$1" == "-c" ]; then
echo "----------------Run Series Verion----------------"
# run cpu version
for size in $series; do
    # create sub folder
    mkdir -p output/cpu/$size
    rm -rf output/cpu/$size/*
    # run cpu version
    echo -n "Running CPU version($size x $size)"
    $prog $infile $size 0 $verbose > "output/cpu/$size/final"
    echo " - $(grep '^\[Seq Time\]' output/cpu/$size/final)"
done
elif [ "$1" == "-o" ]; then
# run omp version
echo "----------------Run OMP Verion----------------"
for size in $series; do
    # create sub folder
    mkdir -p output/omp/$size
    rm -rf output/omp/$size/*
    # run cpu version
    echo -n "Running OMP version($size x $size)"
    $prog $infile $size 1 $verbose > "output/omp/$size/final"
    echo " - $(grep '^\[OMP Time\]' output/omp/$size/final)"
done
else
echo "----------------Run Series Verion----------------"
# run cpu version
for size in $series; do
    # create sub folder
    mkdir -p output/cpu/$size
    rm -rf output/cpu/$size/*
    # run cpu version
    echo -n "Running CPU version($size x $size)"
    $prog $infile $size 0 $verbose > "output/cpu/$size/final"
    echo " - $(grep '^\[Seq Time\]' output/cpu/$size/final)"
done
# run omp version
echo "----------------Run OMP Verion----------------"
for size in $series; do
    # create sub folder
    mkdir -p output/omp/$size
    rm -rf output/omp/$size/*
    # run cpu version
    echo -n "Running OMP version($size x $size)"
    $prog $infile $size 1 $verbose > "output/omp/$size/final"
    echo " - $(grep '^\[OMP Time\]' output/omp/$size/final)"
done
fi