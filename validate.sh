#!/bin/bash

output_folder="output"

if [ ! -d $output_folder ]; then
    echo "Output folder does not exist"
    exit 1
fi

if [ ! -d "$output_folder/cpu" ] || [ ! -d "$output_folder/omp" ]; then
    echo "Output folder does not contain the cpu and omp subfolders"
    exit 1
fi

# loop through the output/cpu and output/omp folder together
for cpu_folder in $output_folder/cpu/*; do
    folder_name=$(basename $cpu_folder)
    omp_folder="${cpu_folder/cpu/omp}"
    if [ ! -d "$omp_folder" ]; then
        echo "Missing corresponding omp folder for $cpu_folder"
        exit 1
    fi
    # check if the final output files exist in both cpu and omp folders
    if [ ! -f "$cpu_folder/final" ] || [ ! -f "$omp_folder/final" ]; then
        echo "Missing final output file in $cpu_folder or $omp_folder"
        exit 1
    fi
    # loop through the final output files in cpu and omp folders with .out extension
    for cpu_file in $cpu_folder/*.out; do
        omp_file="${cpu_file/cpu/omp}"
        if [ ! -f "$omp_file" ]; then
            echo "Missing corresponding omp file for $cpu_file"
            exit 1
        fi
        
        # run the validation prog, if not exit with error code
        ./validator $cpu_file $omp_file
        if [ $? -ne 0 ]; then
            echo "Validation failed for $cpu_file and $omp_file"
            exit 1
        fi
    done
done

echo "Validation successful"
