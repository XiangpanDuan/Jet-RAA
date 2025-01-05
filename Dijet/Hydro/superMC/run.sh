#!/bin/bash
##############################
# Initial Condition
# Author: Xiangpan Duan
# Date: 3-10-2024
##############################

rm -f superMC.e
cd src
make clean
rm -f superMC.e
make
cd ..

rm -rf data
mkdir data
cp src/superMC.e ./
nohup ./superMC.e &

# cp ./data/sdAvg_order_2_block.dat ../VISHYDRO/VISHNew/Initial/
