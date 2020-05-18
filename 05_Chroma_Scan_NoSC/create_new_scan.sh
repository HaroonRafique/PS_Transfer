#!/bin/bash

cp -r 01_01 03_01
cp -r 01_02 03_02
cp -r 01_03 03_03
cp -r 01_04 03_04
cp -r 01_05 03_05
cp -r 01_06 03_06
cp -r 01_07 03_07
cp -r 01_08 03_08
cp -r 01_09 03_09
cp -r 01_10 03_10
cp -r 01_11 03_11

cd 03_01
./clean_all.sh
cd ../03_02
./clean_all.sh
cd ../03_03
./clean_all.sh
cd ../03_04
./clean_all.sh
cd ../03_05
./clean_all.sh
cd ../03_06
./clean_all.sh
cd ../03_07
./clean_all.sh
cd ../03_08
./clean_all.sh
cd ../03_09
./clean_all.sh
cd ../03_10
./clean_all.sh
cd ../03_11
./clean_all.sh
cd ../
