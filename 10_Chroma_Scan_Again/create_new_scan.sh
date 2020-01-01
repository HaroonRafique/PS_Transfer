#!/bin/bash
folder_in=01
folder_out=03

cp -r ${folder_in}_01 ${folder_out}_01
cp -r ${folder_in}_02 ${folder_out}_02
cp -r ${folder_in}_03 ${folder_out}_03
cp -r ${folder_in}_04 ${folder_out}_04
cp -r ${folder_in}_05 ${folder_out}_05
cp -r ${folder_in}_06 ${folder_out}_06
cp -r ${folder_in}_07 ${folder_out}_07
cp -r ${folder_in}_08 ${folder_out}_08
cp -r ${folder_in}_09 ${folder_out}_09
cp -r ${folder_in}_10 ${folder_out}_10
cp -r ${folder_in}_11 ${folder_out}_11

cd ${folder_out}_01
./clean_all.sh
cd ../${folder_out}_02
./clean_all.sh
cd ../${folder_out}_03
./clean_all.sh
cd ../${folder_out}_04
./clean_all.sh
cd ../${folder_out}_05
./clean_all.sh
cd ../${folder_out}_06
./clean_all.sh
cd ../${folder_out}_07
./clean_all.sh
cd ../${folder_out}_08
./clean_all.sh
cd ../${folder_out}_09
./clean_all.sh
cd ../${folder_out}_10
./clean_all.sh
cd ../${folder_out}_11
./clean_all.sh
cd ../
