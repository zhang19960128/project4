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
    double delta_t=0.2*1e-14;//units s, typically this is about 0.02 ps.
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
      double t=5;
      ntsimu(delta_t,r_verlet,t,atomall,30000);
      int steps=50000;
      std::cout<<steps<<" "<<autocorrelate(delta_t,r_verlet,steps,5,atomall)<<std::endl;
}
