#!/bin/bash

#read -p 'Enter value of Q : ' Q
#read -p 'Enter value of tau : ' tau

echo ------------------------------------
echo ------------------------------------

#echo Finding similar pairs within edit distance $tau and q-gram $Q ...
#../GSim/gsim ../data/data_2k.txt $Q $tau -v >> ../Output/GSim_$Q-$tau.txt
#rm aids_result.txt
#../GSim/gsim ../data/data_test.txt 4 6 -v >> ../Output/GSim_2-6.txt

../GSim/gsim ../data/data_2k.txt 4 1 -v >> ../Output/GSim_4-1.txt
../GSim/gsim ../data/data_2k.txt 4 2 -v >> ../Output/GSim_4-2.txt
../GSim/gsim ../data/data_2k.txt 4 3 -v >> ../Output/GSim_4-3.txt
../GSim/gsim ../data/data_2k.txt 4 4 -v >> ../Output/GSim_4-4.txt
../GSim/gsim ../data/data_2k.txt 4 5 -v >> ../Output/GSim_4-5.txt
../GSim/gsim ../data/data_2k.txt 4 6 -v >> ../Output/GSim_4-6.txt


../GSim/gsim ../data/data_5k.txt 4 1 -v >> ../Output/GSim_4-1.txt
../GSim/gsim ../data/data_5k.txt 4 2 -v >> ../Output/GSim_4-2.txt
../GSim/gsim ../data/data_5k.txt 4 3 -v >> ../Output/GSim_4-3.txt
../GSim/gsim ../data/data_5k.txt 4 4 -v >> ../Output/GSim_4-4.txt
../GSim/gsim ../data/data_5k.txt 4 5 -v >> ../Output/GSim_4-5.txt
../GSim/gsim ../data/data_5k.txt 4 6 -v >> ../Output/GSim_4-6.txt



echo ------------------------------------
echo ------------------------------------

echo Complete


