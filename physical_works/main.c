//Date: 07.21.15 
//Description: 18 particle simulation, 4 different species, 6p simulation modified to 18
//this one seems to work properly
//Name: Ramin Khajeh
#include "main.h"

int main(){
   traj = fopen("myTraj_0721_18p_modified6p.xyz","w");
   stateFile = fopen("state_0721.txt","w");
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
      State();
      Renew();
      }
   return 0;
}

