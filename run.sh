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

series="10 20 40 80 100 200 300 400 500"

if [ "$1" == "-c" ]; then
echo "----------------Run Series Verion----------------"
# run cpu version
for size in $series; do
    # create sub folder
    mkdir -p output/cpu/$size
    rm -rf output/cpu/$size/*
    # run cpu version
    echo "Running CPU version($size x $size)"
    $prog $infile $size 0 > "output/cpu/$size/final"
done
elif [ "$1" == "-o" ]; then
# run omp version
echo "----------------Run OMP Verion----------------"
for size in $series; do
    # create sub folder
    mkdir -p output/omp/$size
    rm -rf output/omp/$size/*
    # run cpu version
    echo "Running OMP version($size x $size)"
    $prog $infile $size 1 > "output/omp/$size/final"
done
else
echo "----------------Run Series Verion----------------"
# run cpu version
for size in $series; do
    # create sub folder
    mkdir -p output/cpu/$size
    rm -rf output/cpu/$size/*
    # run cpu version
    echo "Running CPU version($size x $size)"
    $prog $infile $size 0 > "output/cpu/$size/final"
done
# run omp version
echo "----------------Run OMP Verion----------------"
for size in $series; do
    # create sub folder
    mkdir -p output/omp/$size
    rm -rf output/omp/$size/*
    # run cpu version
    echo "Running OMP version($size x $size)"
    $prog $infile $size 1 > "output/omp/$size/final"
done
fi