//Date: 08.22.2015
//Description: Brownian Simulation of RNA self-assembly
//Harvard School of Engineering and Applied Sciences, Summer 2015 REU Program

//Names: 
//       Ramin Khajeh | raminkh@berkeley.edu
//       Laura Pang   | laura.pangx@gmail.com 
//Advisor: 
//       Professor Michael P. Brenner
//       Hidenori Tanaka

#include "main.h"

int main(int argc, char **argv){  						//compilatiion of program in terminal:
														//gcc -o EXEFILE.o main.c func.c -lm
														//./EXEFILE.o OUTPUTFILENAME DIFCOEFFICIENT  				                                                           
   //traj = fopen("myTraj_0723_2.xyz","w");     
   //stateFile = fopen("state_0723_2.txt","w");
   //bond = fopen("bond)9723_2.txt","w");

   char* filename = argv[1];							//READS IN OUTPUT FILENAME AND OPENS FILES
   openFiles(filename);
   char* TempArray = argv[2];							//READS IN DIFFUSION COEFFICIENT
   sscanf(TempArray, "%lf",&T);
 
   T = 80.0 + 25.0*T;									//DIMENSIONALIZED CALCULATIONS
   zi = 9.42E-12;
   Dif = (Kb*T)/zi;
   tau = (diam*diam)/Dif;
   E = 7.0*Kb*(298.0);
   double time = (sim_time*tau*.000002)/h;
   //printf("%lf\t%.10lf\n",T, time);
   
   start = clock();										//INITIAL CLOCK TO KEEP TRACK OF HOW LONG THIS PROGRAM RUNS
   InitialSet();										//INITIALIZES BASE VARIABLES
   NewVerletList();										//INITIALIZES VERLET LIST
   t=0;
   while(t<time){
      t = t + 1;
      newcal();											//BROWNIAN CALCULATION OF NEXT POSITION
      note();      										//NOTES DOWN THE TRAJECTORY
      if(t%Vtime==0){									//HOW OFTEN THE VERLET LIST IS RESET
        NewVerletList();
      }
      IniCA();											//INITIALIZES A CONNECTIVITY ARRAY
      oneBond();										//**METHOD NEEDS TO BE REMADE - ONLY BOND TO ONE OTHER ELEMENT AT A TIME
      State();											//NOTES DOWN THE CURRENT CONNECTIVITY OF ELEMENTS
      for(a=0;a<N;a++){
         int c;
         for(c=0;c<nlist[a];c++){            			//ITERATES THROUGH THE VERLET LIST
            b=list[a][c];
            CalDist();									//CALCULATES THE DISTANCE BETWEEN PARTICLE A AND B
            CalEdepth();								//CALCULATES THE ENERGY DEPTH
            CalForce();									//CALCULATES THE FORCE
         }	
      }
      // for(a=0;a<N-1;a++){
      //     for(b=a+1;b<N;b++){
      //         CalDist();
      //         CalEdepth();
      //         CalForce();
      //     }
      // }
      SumForcesV();										//SUMS THE INDIVIDUAL FORCES (CALCULATED IN CALFORCE)
      //State();
      Renew();											//UPDATES THE VARIABLES
   }
   stop = clock();										//ENDS THE CLOCK THAT KEEPS TRACK OF HOW LONG THE PROGRAM RUNS
   //double timeEllapse = (double)(stop-start)/(CLOCKS_PER_SEC);
   //printf("%lf\n",timeEllapse);
   //fprintf(stateFile,"%lf\n",timeEllapse);
   return 0;
}

