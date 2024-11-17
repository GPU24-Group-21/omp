#!/bin/bash
prog="./ljp"
infile="Rap_2_LJP.in"

# find any .in file and set it as the input file
if [ -f *.in ]; then
    infile=$(ls *.in)
fi

mkdir -p output/cpu
mkdir -p output/omp
rm -rf output/cpu/*
rm -rf output/omp/*

clear

series="10 20 40 80 160 320 640"

# run cpu version
for size in $series; do
    # create sub folder
    mkdir -p output/cpu/$size
    rm -rf output/cpu/$size/*
    # run cpu version
    echo "Running CPU version($size x $size)"
    $prog $infile $size 0 > "output/cpu/$size/final"
done

# run gpu version
echo "------------------------------------"
for size in $series; do
    # create sub folder
    mkdir -p output/omp/$size
    rm -rf output/omp/$size/*
    # run cpu version
    echo "Running OMP version($size x $size)"
    $prog $infile $size 1 > "output/omp/$size/final"
done