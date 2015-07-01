/*--- Brownian Dynamics simulation of 2D Replicating colloidal cluster with the Morse potential and the forward Euler time-step ---*/
#include "main.h"
int main(int argc, char **argv){
	char* difArray = argv[1];
	sscanf(difArray, "%lf",&Dif);
	Dif = 0.6+(Dif*0.02);
	printf("%lf\n",Dif);

	char* filename = argv[2];


	//runMain(fileName);

	char trajName[100]; char amatrixName[100]; char productName[100]; 
	strcpy(trajName,filename); 
	strcat(trajName, "trajChain.xyz");
	strcpy(amatrixName,filename); 
	strcat(amatrixName, "Amatrix.txt");
	strcpy(productName,filename); 
	strcat(productName, "Product.txt");

	printf("%s\t%s\t%s\n",trajName,amatrixName,productName);
	traj = fopen(trajName,"w");
	Amatrix = fopen(amatrixName,"w");
	Product = fopen(productName,"w");
	BP = fopen("myBP.txt","a");
	int i;
	//Dif = dif;
	// for ( i = 0; i<6; i++){
	// 	Dif+=(0.02);
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
        
     //  }
     // printf("%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
}
    printf("%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
    fprintf(BP,"%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
    //fclose(traj);
    return 0;
}




