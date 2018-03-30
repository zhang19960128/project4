#include<iostream>
#include "parameter.h"
#include "atom.h"
#include <cmath>
#include <vector>
#include <list>
#include <fstream>
int main(){
	//***********************initialize the system********units SI**************************//
    int size=20;
    int N=size*size;
    double delta_t=2*1e-14;//units s, typically this is about 0.01 ps.
    std::vector<atom> atomall(N);
    int temp;
    double r_verlet=1.2*r_cut;
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    std::vector<double> totalstress;
    std::vector<double> initial(2,0);
    for(size_t i=0;i<size;i++)
     for(size_t j=0;j<size;j++){
        temp=i*size+j;
        atomall[temp].setx(i*r_min);
        atomall[temp].sety(j*r_min);
        atomall[temp].setr(0.1*r_min);
        atomall[temp].setm(39.948*1e-3/NA);
        atomall[temp].setv(initial);
        atomall[temp].setf(initial);
       }
    //**********************end initialize ****************************************//
    cool(delta_t,r_verlet,atomall);
    updatetensor(atomall);
    //ntsimu(delta_t,r_verlet,80,atomall,10000);
    for(size_t i=0;i<size*size;i++){
        atomall[i].printinfo();
    }
}
