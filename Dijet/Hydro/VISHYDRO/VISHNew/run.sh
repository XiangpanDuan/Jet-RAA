#!/bin/bash
##############################
# Author: Xiangpan Duan
# Date: 21-10-2024
##############################


rm -f Initial/*.dat
cp ../../superMC/data/sdAvg_order_2_block.dat Initial/
mv Initial/sdAvg_order_2_block.dat Initial/InitialSd.dat


rm -f VISHNew.e
cd src
# rm -rf VISHNew.e obj
make distclean
make
cd ..

rm -rf results
mkdir results

cp src/VISHNew.e ./
./VISHNew.e


# mkdir build
# # sed -i -e '41i    set (CMAKE_Fortran_FLAGS "-w")\' CMakeLists.txt
# cd build
# cmake ..
# make
# cp src/VISHNew.e ../
# cd ..
# rm -rf results
# mkdir results
# ./VISHNew.e
# rm -rf build
