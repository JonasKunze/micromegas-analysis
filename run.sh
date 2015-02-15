#!/bin/bash
source /usr/lib/root/bin/thisroot.sh

cd Debug
make clean
make

mkdir /localscratch/praktikum/output/
time ./micromegas-analysis > /localscratch/praktikum/output/cutStats.txt
