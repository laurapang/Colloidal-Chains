//Date: 07.21.15 
//Description: 18 particle simulation, 4 different species, 6p simulation modified to 18
//this one seems to work properly
//Name: Ramin Khajeh
#include "main.h"

int main(int argc, char **argv){

   // traj = fopen("myTraj_0723_2.xyz","w");
   // stateFile = fopen("state_0723_2.txt","w");
   // bond = fopen("bond)9723_2.txt","w");
  char* filename = argv[1];
    openFiles(filename);
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
    double timeEllapse = (double)(stop-start)/(CLOCKS_PER_SEC);
    printf("%lf\n",timeEllapse);
    fprintf(stateFile,"%lf\n",timeEllapse);
   return 0;
}

