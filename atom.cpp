#include "atom.h"
#include "parameter.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
void atom::setx(double x1){
	x=x1;
}
void atom::sety(double x2){
	y=x2;
}
void atom::setr(double ra){
	radius=ra;
}
void atom::setstress_tensor(std::vector<double> a){
        stresstensor=a;
}
void atom::setm(double m){
   mass=m;
}
std::vector<double> atom::getstress(){
	return stresstensor;
}
void atom::setv(std::vector<double> v){
        speed=v;
}
void atom::setf(std::vector<double> f){
   force=f;
}
double atom::getx(){
	return x;
}
double atom::gety(){
	return y;
}
double atom::getr(){
	return radius;
}
void atom::printneighbor(){
		std::cout<<"the neighbors are:"<<std::endl;
		for(std::list<int>::iterator a=neighbor.begin();a!=neighbor.end();a++){
			std::cout<<*a<<" "<<std::endl;
		}
	}
void atom::printstress(){
	std::cout<<stresstensor[0]<<" "<<stresstensor[1]<<" "<<stresstensor[2]<<" "<<stresstensor[3]<<std::endl;
}
void atom::printinfo(){
	std::cout<<x<<" "<<y<<" "<<speed[0]<<" "<<speed[1]<<" "<<force[0]<<" "<<force[1]<<std::endl;
}
void atom::printforce(){
    std::cout<<"F(x):"<<force[0]<<std::endl;
    std::cout<<"F(y):"<<force[1]<<std::endl;
}
double distance(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	return sqrt(r);
}
double dxx(atom& one,atom& two){
	return one.x-two.x;
}
double dyy(atom& one,atom& two){
	return one.y-two.y;
}
double potential(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	r=sqrt(r);
	if(r<r0){
		return 4*eps*(pow(sigma/r,12)-pow(sigma/r,6));
	}
	else if(r>r_cut){
		return 0;
	}
	else{
		return A*pow(r-r_cut,3)+B*pow(r-r_cut,2);
	}
}
double allpotential(std::vector<atom>& allatom){
    int size=allatom.size();
    double sum=0.0;
    double temp=0.0;
    for(size_t i=0;i<size;i++){
       for(std::list<int>::iterator a=allatom[i].neighbor.begin();a!=allatom[i].neighbor.end();a++){
        temp=potential(allatom[i],allatom[*a]);
        sum=sum+temp;
       }
    }
    return sum/2;
}
atom& atom::operator =(atom& one){
	this->mass=one.mass;
	this->x=one.x;
	this->y=one.y;
	this->speed=one.speed;
  this->force=one.force;
	return *this;
}
double allener(std::vector<atom>& allatom){
    double poten=allpotential(allatom);
    double kinet=0.0;
    for(size_t i=0;i<allatom.size();i++){
        kinet=kinet+0.5*(allatom[i].mass)*(allatom[i].speed[0]*allatom[i].speed[0]+allatom[i].speed[1]*allatom[i].speed[1]);
    }
    return kinet+poten;
}
//the force exerted on one by two
std::vector<double> str_tensor(atom& one,atom& two){
    std::vector<double> a(4,0);
    double r=distance(one,two);
    double deri=0;
    if(r<r0){
        deri=24*eps/r*(-2*pow(sigma/r,12)+pow(sigma/r,6));
    }
    else if(r<r_cut){
        deri=3*A*(r-r_cut)*(r-r_cut)+2*B*(r-r_cut);
    }
    else{
        deri=0.0;
    }
    a[0]=deri*(one.x-two.x)*(one.x-two.x)/r;
    a[1]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[2]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[3]=deri*(one.y-two.y)*(one.y-two.y)/r;
    return a;
}
std::vector<double> cal_force(atom& one,atom& two){
   std::vector<double> a(2,0);
   double r=distance(one,two);
   double deri=0;
   if(r<r0){
      deri=24*eps/r*(-2*pow(sigma/r,12)+pow(sigma/r,6));
   }
   else if(r<r_cut){
      deri=3*A*(r-r_cut)*(r-r_cut)+2*B*(r-r_cut);
   }
   else{
      deri=0.0;
   }
   a[0]=-1*deri*(one.x-two.x)/r;
   a[1]=-1*deri*(one.y-two.y)/r;
   return a;
}
std::vector<double>& operator +=(std::vector<double>& one,std::vector<double>& two){
   for(size_t i=0;i<one.size();i++){
      one[i]=one[i]+two[i];
   }
   return one;
}
void atom::updateforce(std::vector<atom>& atomall){
   double r=0;
   std::vector<double> allforce(3,0);
   std::vector<double> temp(3,0);
   for(std::list<int>::iterator a=neighbor.begin();a!=neighbor.end();a++){
      temp=cal_force(*this,atomall[*a]);
      allforce+=temp;
   }
   this->force=allforce;
}
double atom::updateposition(double delta_t){
   double deltax=speed[0]*delta_t+delta_t*delta_t/2/mass*force[0];
   double deltay=speed[1]*delta_t+delta_t*delta_t/2/mass*force[1];
   x=x+deltax;
   y=y+deltay;
   return sqrt(deltax*deltax+deltay*deltay);
}
double hydropressure(std::vector<atom> atomall){
    int size=atomall.size();
    double xx=0.0;
    double yy=0.0;
    for(size_t i=0;i<size;i++){
        xx=xx+atomall[i].stresstensor[0];
        yy=yy+atomall[i].stresstensor[3];
    }
    return xx+yy;
}
void updatelist(std::vector<atom>& input,double r_set){
    int size=input.size();
        for(size_t i=0;i<size;i++){
	    input[i].neighbor.clear();
	   }
	for(size_t i=0;i<size;i++)
		for(size_t j=0;j<size;j++){
                    if(i==j) continue;
               	  if(distance(input[i],input[j])<r_set && distance(input[i],input[j])>0.0){
		     input[i].neighbor.push_back(j);
        }
   }
}
double updateallposition(std::vector<atom>& atomall,double delta_t){
    double maxdis=0.0;
    double tempdis=0.0;
    for(size_t i=0;i<atomall.size();i++){
      tempdis=atomall[i].updateposition(delta_t);
      if(tempdis>maxdis){
         maxdis=tempdis;
      }
   }
    return maxdis;
}
double verletrun(double delta_t,std::vector<atom>& allatom){
   int size=allatom.size();
   double maxdis=0.0;
   double tempdis=0.0;
   std::vector<double> f_before(2,0.0);
   for(size_t i=0;i<size;i++){
      f_before=allatom[i].force;
      tempdis=allatom[i].updateposition(delta_t);
      if(tempdis>maxdis){
        maxdis=tempdis;
      }
      allatom[i].updateforce(allatom);
      allatom[i].speed[0]=allatom[i].speed[0]+delta_t/2/allatom[i].mass*(allatom[i].force[0]+f_before[0]);
      allatom[i].speed[1]=allatom[i].speed[1]+delta_t/2/allatom[i].mass*(allatom[i].force[1]+f_before[1]);
   }
   return maxdis;
}
double verletrun_standard(double delta_t,std::vector<atom>& allatom){
   int size=allatom.size();
   double maxdis=0.0;
   double tempdis=0.0;
   std::vector<double> f_before(2,0.0);
   maxdis=updateallposition(allatom,delta_t);
   for(size_t i=0;i<size;i++){
      f_before=allatom[i].force;
//      tempdis=allatom[i].updateposition(delta_t);
//      if(tempdis>maxdis){
//         maxdis=tempdis;
//      }
      allatom[i].updateforce(allatom);
      allatom[i].speed[0]=allatom[i].speed[0]+delta_t/2/allatom[i].mass*(allatom[i].force[0]+f_before[0]);
      allatom[i].speed[1]=allatom[i].speed[1]+delta_t/2/allatom[i].mass*(allatom[i].force[1]+f_before[1]);
   }
   return maxdis;
}
void updatetensor(std::vector<atom>& atomall){
    std::vector<double> temp(4,0);
    std::vector<double> all(4,0);
    int size=atomall.size();
    for(size_t i=0;i<size;i++){
        for(size_t k=0;k<4;k++){
            all[k]=0.0;
        }
        for(std::list<int>::iterator a=atomall[i].neighbor.begin();a!=atomall[i].neighbor.end();a++){
            temp=str_tensor(atomall[i],atomall[*a]);
            for(size_t k=0;k<4;k++){
                all[k]=all[k]+temp[k];
            }
        }
    atomall[i].setstress_tensor(all);
    }
}
std::ostream& operator<<(std::ostream& os,atom& output){
		os<<"for atom position ("<<output.x<<","<<output.y<<")"<<std::endl;
		os<<"stress tensor(xx,xy,yx,yy): "<<output.stresstensor[0]<<" "<<output.stresstensor[1]<<" "<<output.stresstensor[2]<<" "<<output.stresstensor[3]<<std::endl;
		return os;
}
std::fstream& operator<<(std::fstream& os,atom& output){
		os<<output.x<<" "<<output.y<<" "<<output.stresstensor[0]<<" "<<output.stresstensor[1]<<" "<<output.stresstensor[2]<<" "<<output.stresstensor[3];
		return os;
}
int count(std::vector<atom> all,atom& input,double r){
    int c=0;
    int size=all.size();
    for(int i=0;i<size;i++){
        if(distance(all[i],input)<r){
            c++;
        }
    }
    return c-1;
}
void print_radial_dis(double r_start,double r_stop,std::vector<atom>& atomall,std::string name){
	 // double r_start=0.0000001;
  //  double r_stop=15;
    int N=10000;
    int size=atomall.size();
    std::vector<double> ra_dis(N,0.0);
    std::vector<double> ra_dis_all(N,0.0);
    std::vector<double> rinter(N,0.0);
    double r_delta=(r_stop-r_start)/N;
    double r_inter=0.0;
    int count_old=0;
    int count_new=0;
    int count_delta=0;
    std::fstream radis;
    radis.open(name,std::fstream::out);
    for(size_t j=0;j<size;j++){
    for(size_t i=0;i<N;i++){
        r_inter=i*r_delta+r_start;
        rinter[i]=r_inter;//record the radius.
        count_new=count(atomall,atomall[j],r_inter);
        count_delta=count_new-count_old;
        ra_dis[i]=count_delta/2/Pi/r_inter/r_delta;
        count_old=count_new;
       // radis<<r_inter<<" "<<ra_dis[i]<<std::endl;
    }
    std::cout<<j<<std::endl;
    ra_dis_all+=ra_dis;
    }
    for(size_t j=0;j<size;j++){
        ra_dis_all[j]=ra_dis_all[j]/size;
    }
    for(size_t i=0;i<N;i++){
        radis<<rinter[i]<<" "<<ra_dis_all[i]<<std::endl;
    }
    radis.close();
}
double temperature(std::vector<atom>& allatom){
   double sum=0;
   double vx;
   double vy;
   int size=allatom.size();
   for(size_t i=0;i<size;i++){
      vx=allatom[i].speed[0];
      vy=allatom[i].speed[1];
      sum=sum+0.5*allatom[i].mass*(vx*vx+vy*vy);
   }
   return sum/(Kb)/size;
}
void freeze(std::vector<atom>& allatom){
   int size=allatom.size();
   for(size_t i=0;i<size;i++){
      allatom[i].speed[0]=0.0;
      allatom[i].speed[1]=0.0;
   }
}
void settemp(double t,std::vector<atom>& allatom){
  int size=allatom.size();
  double temp=temperature(allatom);
  for(size_t i=0;i<size;i++){
   for(size_t j=0;j<2;j++){
      allatom[i].speed[j]=allatom[i].speed[j]*sqrt(t/temp);
   }
  }
}
void cool(double delta_t,double r_verlet,std::vector<atom>& atomall){
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    int count=0;
    double temp_before=0.0;
    double temp_now=0.0;
    double maxdis=0.0;
    int neg_count=5;
    double p_before=0.0;
    double p_end=0.0;
    updatelist(atomall,r_cut);
    do{
        p_before=p_end;
        temp_before=temp_now;
        /*force has already been updated on the updateallposition*/
        //****verletrun contains velocity update and force update****//
        maxdis=verletrun(delta_t,atomall);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        count++;
        temp_now=temperature(atomall);
        if((temp_now-temp_before)< 0){
            neg_count++;
            if(temp_now < 2) break;
            if(neg_count>2){
               freeze(atomall);
               neg_count=0;
            }
            temp_now=temperature(atomall);
        }
        p_end=allpotential(atomall);
    }while(std::fabs((p_end-p_before)/p_end)>1e-9);
}
void equilibrium(double delta_t,double r_verlet,std::vector<atom>& atomall){
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    int count=0;
    double temp_before=0.0;
    double temp_now=0.0;
    double maxdis=0.0;
    int neg_count=5;
    updatelist(atomall,r_cut);
    do{
        temp_before=temp_now;
        //****verletrun contains velocity update and force update****//
        maxdis=verletrun(delta_t,atomall);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        count++;
        temp_now=temperature(atomall);
        std::cout<<temp_now<<std::endl;
    }while(1);
}
void ntsimu(double delta_t,double r_verlet,double t,std::vector<atom>& atomall,int steps){
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    int count=0;
    double temp_before=0.0;
    double temp_now=0.0;
    double maxdis=0.0;
		double screen=10.0;
    int neg_count=5;
    updatelist(atomall,r_cut);
		settemp(t,atomall);
    do{
        count++;
        temp_before=temp_now;
        //****verletrun contains velocity update and force update****//
        maxdis=verletrun_standard(delta_t,atomall);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        count++;
        temp_now=temperature(atomall);
				std::cout<<temp_now<<std::endl;
				if(std::fabs(temp_now-t)>screen){
        settemp(t,atomall);
				}
    }while(count<steps);
}
void nesimu(double delta_t,double r_verlet,int steps,std::vector<atom>& atomall){
   double r_shell_initial=r_verlet-r_cut;
   double r_shell=r_shell_initial;
    int count=0;
    double temp_now=0.0;
    double maxdis=0.0;
    double totalenergy=0.0;
    updatelist(atomall,r_cut);
    do{
        /*force has already been updated on the updateallposition*/
        //****verletrun contains velocity update and force update****//
        maxdis=verletrun_standard(delta_t,atomall);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        count++;
        std::cout<<"all energy: "<<allener(atomall)<<"potential energy: "<<allpotential(atomall)<<std::endl;
    }while(count<steps);
}
double autocorrelate(double delta_t,double r_verlet,int steps,double t,std::vector<atom>& atomall){
    int runtime=steps*2;
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    int size=atomall.size();
    double temp_before=0.0;
    double temp_now=0.0;
    double maxdis=0.0;
		double screen=1.0;
    updatelist(atomall,r_cut);
    std::vector<double> a(runtime,0.0);
    std::vector<std::vector<double>> allx(size,a);
    std::vector<std::vector<double>> ally(size,a);
    for(size_t i=0;i<runtime;i++){
        temp_before=temp_now;
        //****verletrun contains velocity update and force update****//
        maxdis=verletrun_standard(delta_t,atomall);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
        temp_now=temperature(atomall);
				if(std::fabs(temp_now-t)>screen){
             settemp(t,atomall);
				}
        for(size_t j=0;j<size;j++){
            allx[j][i]=atomall[j].speed[0];
            ally[j][i]=atomall[j].speed[1];
            }
         }
    std::vector<double> correx(size,0.0);
    std::vector<double> correy(size,0.0);
    for(size_t i=0;i<size;i++){
       for(size_t j=0;j<steps;j++){
         correx[i]=correx[i]+allx[i][j]*allx[i][j+steps];
         correy[i]=correy[i]+ally[i][j]*ally[i][j+steps];
       }
    correx[i]=correx[i]/steps;
    correy[i]=correy[i]/steps;
    }
    double sumx=0.0;
    double sumy=0.0;
    for(size_t i=0;i<size;i++){
      sumx=sumx+correx[i];
      sumy=sumy+correy[i];
    }
    return (sumx/size+sumy/size)/2;
}
void diffusive(double delta_t,double r_verlet,double t,std::vector<atom>& atomall,int steps){
    double r_shell_initial=r_verlet-r_cut;
    double r_shell=r_shell_initial;
    int count=0;
    double maxdis=0.0;
		double screen=10.0;
		double temp_now=0.0;
		double dall=0.0;
		std::fstream diff;
		diff.open("diffusive.txt",std::fstream::out);
		std::vector<atom> initial(atomall.size());
		for(size_t i=0;i<atomall.size();i++){
			initial[i]=atomall[i];
		}
    updatelist(atomall,r_cut);
    do{
        count++;
        //****verletrun contains velocity update and force update****//
        maxdis=verletrun_standard(delta_t,atomall);
        if(maxdis*2>r_shell){
            updatelist(atomall,r_verlet);
            r_shell=r_shell_initial;
        }
        else{
            r_shell=r_shell-2*maxdis;
        }
				dall=0.0;
				for(size_t i=0;i<atomall.size();i++){
					dall=dall+dxx(initial[i],atomall[i])*dxx(initial[i],atomall[i])+dyy(initial[i],atomall[i])*dyy(initial[i],atomall[i]);
				}
				diff<<count<<" "<<dall/atomall.size()<<std::endl;
    }while(count<steps);
}
