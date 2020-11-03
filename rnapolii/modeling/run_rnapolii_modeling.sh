#!/bin/bash

# Usage:
# ./run_rnapolii_modeling.py out_dir n n_steps

py=/home/saltzberg/swr/anaconda3/bin/python

out_dir=$1
n=$2
n_steps=$3

for x in $( seq 1 $n) 
do
    $py ./modeling.py $out_dir $x $n_steps & > ./output_files/$out_dir$x.out
done
