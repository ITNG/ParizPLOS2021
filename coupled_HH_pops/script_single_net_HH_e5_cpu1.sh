#!/bin/bash
for dI in `seq 0 3`
    do
    for deeout in `seq 0 59`
        do
        gfortran -m64 -O3 hh_mod.f95 main_hh.f95 -o a1.out
        echo $deeout $dI | ./a1.out 
        rename "s/ //g" *.txt
        nohup matlab -nodisplay <rho_creator.m> a1.dat 
        rm rho*.txt
    done
done
