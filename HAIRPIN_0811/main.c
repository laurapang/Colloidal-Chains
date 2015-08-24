//Date: 08.22.2015
//Description: Brownian Simulation of RNA self-assembly
//Harvard School of Engineering and Applied Sciences, Summer 2015 REU Program

//Names: 
//       Ramin Khajeh | raminkh@berkeley.edu
//       Laura Pan    | laura.pangx@gmail.com 
//Advisor: 
//       Professor Michael P. Brenner
//       Hidenori Tanaka

#include "main.h"

int main(int argc, char **argv){                
                                                
   //traj = fopen("myTraj_0723_2.xyz","w");     
   //stateFile = fopen("state_0723_2.txt","w");
   //bond = fopen("bond)9723_2.txt","w");

   char* filename = argv[1];
   openFiles(filename);
   char* TempArray = argv[2];
   sscanf(TempArray, "%lf",&T);
 
   T = 80.0 + 25.0*T;
   zi = 9.42E-12;
   Dif = (Kb*T)/zi;
   tau = (diam*diam)/Dif;
   E = 7.0*Kb*(298.0);
   double time = (sim_time*tau*.000002)/h;
   //printf("%lf\t%.10lf\n",T, time);
   
   start = clock();
   InitialSet();
   NewVerletList();
   t=0;
   while(t<time){
      t = t + 1;
      newcal();
      note();      
      if(t%Vtime==0){
        NewVerletList();
      }
      IniCA();
      oneBond();
      State();
      for(a=0;a<N;a++){
         int c;
         for(c=0;c<nlist[a];c++){            
            b=list[a][c];
            CalDist();
            CalEdepth();
            CalForce();
         }
      }
      // for(a=0;a<N-1;a++){
      //     for(b=a+1;b<N;b++){
      //         CalDist();
      //         CalEdepth();
      //         CalForce();
      //     }
      // }
      SumForcesV();
      //State();
      Renew();
   }
   stop = clock();
   //double timeEllapse = (double)(stop-start)/(CLOCKS_PER_SEC);
   //printf("%lf\n",timeEllapse);
   //fprintf(stateFile,"%lf\n",timeEllapse);
   return 0;
}

