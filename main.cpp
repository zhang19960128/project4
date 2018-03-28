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
    double delta_t=1e10;//units s
    double maxdis=0.0;
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
    updatelist(atomall,r_cut);
    int count=0;
    double t_before=0.0;
    double t_end=0.0;
    double tc=10;
    do{
        t_before=t_end;
        /*force has already been updated on the updateallposition*/
        maxdis=updateallposition(atomall,delta_t);
        //****verletrun contains velocity update and force update****//
        verletrun(delta_t,atomall);
        updatelist(atomall,r_cut);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        count++;
        t_end=temperature(atomall);
        if(t_end>tc){
         freeze(atomall);
         if(tc>2)  tc=tc-1;
        }
        std::cout<<count<<" temp= "<<t_end<<std::endl;
    }while(std::fabs(t_end-t_before)>0.01 || tc>1);
}
