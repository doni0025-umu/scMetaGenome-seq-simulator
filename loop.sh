#!/bin/bash

for depth in 0.001 0.1 0.5 1 5 10 50 100
do
    bash sim_to_bgzip_runner.sh -n "cpynm10_b-caccae_b-uniformis" -d $depth -r
done