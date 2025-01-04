#!/bin/bash
##############################
# Author: Xiangpan Duan
# Date: 21-10-2024
##############################


rm -f JetCtl.dat JetData.dat
cp ../../../../VISHNew/results/JetCtl.dat  ./
cp ../../../../VISHNew/results/JetData.dat ./

# Change output file name in Calculate.for
make clean
make
./Calculate.exe
