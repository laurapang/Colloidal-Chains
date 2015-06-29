/*--- Brownian Dynamics simulation of 2D Replicating colloidal cluster with the Morse potential and the forward Euler time-step ---*/
#include "main.h"
int main(){

traj = fopen("0629_08trajChain.xyz","w");
Amatrix = fopen("0629_08Amatrix.txt","w");
Product = fopen("0629_08Product.txt","w");
	int i;
	Dif = .42;
	for ( i = 0; i<6; i++){
		Dif+=(0.02);
        InitialSet();
    /* time integration loop */
        t=0;
    while(t<time){    t = t + 1;
        newcal();               // given Xold , FXpold, FbX (so as Y) calculate Xnew based on forward Euler time-step algorithm
        note();                 // write down Xnew we've got on a file

		//fprintf(traj,"%d\n",1);
    /*--- calculations of pairwise quantities (distance, interactions etc...) ---*/
            for(a=0;a<N-1;a++){
                for(b=a+1;b<N;b++){
                    CalDist();              // Given configuration of particles calculate DX, R[a][b] between each pair of particles.
                    CalEdepth();            // calculate depth of potential energy depending on its situation
                    CalForce();             // calculate force acting between pairs of particle
                    
                }
            }
       
        SumForces();
        AdMatrix();
        Renew();
        
      }
     printf("%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
}
    //printf("%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
    fclose(traj);
    return 0;
            }




