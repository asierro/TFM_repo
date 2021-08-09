#!/bin/bash

bash gen_hostfile.sh

module load SIESTA/4.1.5-intel-2019b

echo 'Loaded'

unset I_MPI_HYDRA_BOOTSTRAP
unset I_MPI_PMI_LIBRARY

export NPROCS=`wc -l < hostfile`

mpirun -np $NPROCS --hostfile hostfile siesta RUN_0_0.fdf

# RUN AS: nohup bash {filename}.sh &
