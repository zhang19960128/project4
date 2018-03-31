#!/bin/bash
declare -a temp=(5 5.1 10 10.1 15 15.1 20 20.1 25 25.1 30 30.1 35 35.1 40 40.1 45 45.1 50 50.1 55 55.1 60 60.1 65 65.1 70 70.1 75 75.1 80 80.1)
j=5;
for i in ${temp[@]}
do
    echo $i
    sed -i -e "s/double t=$j/""double t=$i/g" main.cpp
    g++ -std=c++11 main.cpp atom.cpp -o out
    ./out >>energradient.txt
    j=$i;
done
