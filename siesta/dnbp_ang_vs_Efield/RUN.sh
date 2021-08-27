#!/bin/bash

np001=24  # Number of cores in atlas-001
np002=24   # Number of cores in atlas-002

>hostfile
if [ $np001 -gt 0 ]
then
    for ((i=1;i<=np001;i++))
    do
        echo 'atlas-001'>> hostfile
    done
fi
if [ $np002 -gt 0 ]
then
    for ((i=1;i<=np002;i++))
    do
        echo 'atlas-002' >> hostfile
    done
fi

module load SIESTA/4.1.5-intel-2019b
echo 'Loaded SIESTA'

unset I_MPI_HYDRA_BOOTSTRAP
unset I_MPI_PMI_LIBRARY

export NPROCS=`wc -l < hostfile`
mpirun -np $NPROCS --hostfile hostfile siesta RUN.fdf >& RUN.out

# RUN AS: nohup bash RUN.sh &
