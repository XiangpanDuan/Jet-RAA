#!/bin/bash
##############################
# Author: Xiangpan Duan
# Date: 21-10-2024
##############################


rm -f JetCtl.dat JetData.dat
cp ../../../VISHNew_Code/results_5/JetCtl.dat ./
cp ../../../VISHNew_Code/results_5/JetData.dat ./

# Change output file name in vish_gen.f90
make clean
make
nohup ./generate.exe &
