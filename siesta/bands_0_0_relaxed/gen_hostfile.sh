#!/bin/bash

np001=24  # Number of cores in atlas-001
np002=0   # Number of cores in atlas-002

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
