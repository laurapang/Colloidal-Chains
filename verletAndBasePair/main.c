/*--- Brownian Dynamics simulation of 2D Replicating colloidal cluster with the Morse potential and the forward Euler time-step ---*/
#include "main.h"
int main(int argc, char **argv){

  char* difArray = argv[1];
  sscanf(difArray, "%lf",&Dif);
  Dif = .5+(Dif*0.02);
  printf("%lf\n",Dif);
char* filename = argv[2];
   

    openFiles(filename);


        InitialSet();
    /* time integration loop */
        t=0;
        newcal();               // given Xold , FXpold, FbX (so as Y) calculate Xnew based on forward Euler time-step algorithm
        note();
        NewVerletList();
    while(t<time){    
        t = t + 1;
        newcal();               // given Xold , FXpold, FbX (so as Y) calculate Xnew based on forward Euler time-step algorithm
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
            
            
        CalDist();              // Given configuration of particles calculate DX, R[a][b] between each pair of particles.
        CalEdepth();            // calculate depth of potential energy depending on its situation
        CalForce();             // calculate force acting between pairs of particle
            
        }
    }


        
        SumForces();   
        Renew();
        
    }
    fclose(traj);
    return 0;
            }




