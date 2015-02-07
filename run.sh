#!/bin/bash
source /usr/lib/root/bin/thisroot.sh

cd Debug
make clean
make

time ./micromegas-analysis > /localscratch/praktikum/output/cutStats.txt
