#!/bin/bash

read -p 'Enter value of Q : ' Q
read -p 'Enter value of tau : ' tau

echo ------------------------------------
echo ------------------------------------

echo Finding similar pairs within edit distance $tau and q-gram $Q ...
../GSim/gsim ../data/data_test2.txt $Q $tau -v >> ../Output/GSim_$Q-$tau.txt


echo ------------------------------------
echo ------------------------------------

echo Complete