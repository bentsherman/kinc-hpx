#!/bin/bash

module purge
module add use.own
module add hpx/1.4.1

BSIZE_VALUES="256 512 1024 2048 4096"
NP_VALUES="1 2 4 8 16"
EMX_FILE="Yeast.emx.txt"
CMX_FILE="Yeast.cmx.txt"

printf "method\tbsize\tnp\truntime\n"

# run kinc-omp
for BSIZE in ${BSIZE_VALUES}; do
    for NP in ${NP_VALUES}; do
        ./kinc-omp \
            ${BSIZE} \
            ${NP} \
            ${EMX_FILE} \
            ${CMX_FILE}
    done
done

# run kinc-hpx
for BSIZE in ${BSIZE_VALUES}; do
    for NP in ${NP_VALUES}; do
        ./kinc-hpx \
            --bsize ${BSIZE} \
            --hpx:threads ${NP} \
            --infile ${EMX_FILE} \
            --outfile ${CMX_FILE}
    done
done
