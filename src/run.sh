#!/bin/bash
source /usr/lib/root/bin/thisroot.sh
make
time ./MMPlots > /localscratch/praktikum/output/cutStats.txt
