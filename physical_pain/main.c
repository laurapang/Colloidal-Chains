//Date: 07.07.15 
//Description: 18 particle simulation, 4 different species,
//Name: Ramin Khajeh
#include "main.h"

int main(){
    
   traj = fopen("myTraj_0721.xyz","w");
   InitialSet();
   t=0;
   while(t<time){
      t = t + 1;
      newcal();
      note();      
      for(a=0;a<N-1;a++){
          for(b=a+1;b<N;b++){
              CalDist();
              CalEdepth();
              CalForce();
          }
      }
      SumForces();
      Renew();
      }
   return 0;
}

