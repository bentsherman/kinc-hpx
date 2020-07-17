#!/bin/bash

module purge
module add anaconda3/5.1.0
module add cuda-toolkit/10.1.168
module add julia/1.1.1

source activate myenv

# define parameters
N_TRIALS=10
EMX_FILE_CPU="Yeast-100.emx.txt"
EMX_FILE_GPU="Yeast-1000.emx.txt"
CMX_FILE_CPU="Yeast-100.cmx.txt"
CMX_FILE_GPU="Yeast-1000.cmx.txt"

printf "implementation\thardware\truntime\n"

# test KINC/C++
for N in $(seq ${N_TRIALS}); do
    printf "c++\tcpu\t"
    env time -f "%e" ~/workspace/kinc-c++/kinc.sh serial 1 ${EMX_FILE_CPU}
done

for N in $(seq ${N_TRIALS}); do
    printf "c++\tgpu\t"
    env time -f "%e" ~/workspace/kinc-c++/kinc.sh cuda   1 ${EMX_FILE_GPU}
done

# test KINC/Numba
for N in $(seq ${N_TRIALS}); do
    printf "numba\tcpu\t"
    env time -f "%e" python ~/workspace/kinc-numba/kinc-numba.py ${EMX_FILE_CPU} ${CMX_FILE_CPU} 0
done

for N in $(seq ${N_TRIALS}); do
    printf "numba\tgpu\t"
    env time -f "%e" python ~/workspace/kinc-numba/kinc-numba.py ${EMX_FILE_GPU} ${CMX_FILE_GPU} 1
done

# test KINC/Julia
for N in $(seq ${N_TRIALS}); do
    printf "julia\tcpu\t"
    env time -f "%e" julia ~/workspace/kinc-julia/kinc-julia.jl ${EMX_FILE_CPU} ${CMX_FILE_CPU} 0
done

# for N in $(seq ${N_TRIALS}); do
#     printf "julia\tgpu\t"
#     env time -f "%e" julia ~/workspace/kinc-julia/kinc-julia.jl ${EMX_FILE_GPU} ${CMX_FILE_GPU} 1
# done
