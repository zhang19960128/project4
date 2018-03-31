#!/bin/bash
declare -a temp=(5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80)
j=5;
for i in ${temp[@]}
do
    echo $i
    sed -i -e "s/double t=$j/""double t=$i/g" main.cpp
    g++ -std=c++11 main.cpp atom.cpp -o out
    ./out >>ener.txt
    j=$i;
done
