/*---Brownian Dynamics simulation of 2D colloidal cluster with the Morse potential and the forward Euler time-step---*/
#include "main.h"
void InitialSet(){
// first, we create and open files to note our result
/* set initial conditions */
IniConf();                  // IniCond(Initial Condition) => here we set initial condition of particles
Roff = 1+(4/Rho) ;          // Roff(Cut-off distance of pariwise force)   Roff  = 1+(4/30) = 1.133
IniSV();                    // Species vector
IniIM();   					// set 9 by 9 interaction matrix
count_a = 0; 
count_b=0;        
count_c=0;
count_d=0;
count_e=0;         
}



//%%%%%%%%%%%%%%%%% Initial Species Vector %%%%%%%%%%%%%%%%%//
// IniSpVec => Initial Species Vector; we assign initial "Species" for each particle.
void IniSV(){
    int i;
    for(i=0;i<N;i++){
        SV[i]=0;
    }
             }

//%%%%%%%%%%%%%%%%% Interaction Matrix %%%%%%%%%%%%%%%%%%%%//
/* Interaction Matrix between species */
void IniIM(){
    IM[0][0]= 1; 
                // IM[0][1]= 1; IM[0][2]=1;
                //  IM[1][1]= 0; IM[1][2]=1;
                //               IM[2][2]=0;
    

    
// We complete unfilled components of the "Interaction Matrix"
        int i,j;
    for(i=0;i<8;i++){
        for(j=i+1;j<9;j++){
                   IM[j][i]=IM[i][j];
                           }
                     }
              }


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// IniConf => Initial Configuration; function to set initial condition //
void IniConf(){
/* first we set initial configurations, so that particles make chain */
    int i,j;
    for(i=0;i<N;i++){
        rold[i][0]=i;
        for(j=1;j<Dim;j++){
            rold[i][j]=0;
        }
    }
    

    // rold[0][0]= 0 ; rold[0][1]= -.707; rold[0][2] = 0;
    // rold[1][0]= .707 ; rold[1][1]= 0; rold[1][2] = 0;
    // rold[2][0]= 0 ; rold[2][1]= 0; rold[2][2] = -.707;
    // rold[3][0]= 0 ; rold[3][1]= .707; rold[3][2] = 0;
    // rold[4][0]= -.707 ; rold[4][1]= 0; rold[4][2] = 0;
    // rold[5][0]= 0 ; rold[5][1]= 0; rold[5][2] = .707;


	   // int i,j;
    // for(i=0;i<N;i++){
    //     rold[i][0]=i%3;
    //     for(j=1;j<Dim;j++){
    //         rold[i][j]=0;
    //     }
    // }

   	// rold[3][1]=1;
   	// rold[4][1]=1;
   	// rold[5][1]=1;

    }




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// given rold and fpold, return rnew
void newcal(){
    
    int i,j;
    for(i=0;i<N;i++){        // calculate positions of next time step
        for(j=0;j<Dim;j++){
        rnew[i][j]  =  rold[i][j]  + fpold[i][j]*h + sqrt(2.0*Dif*h)*RandNormal();  // over damped Langevin equation
        rnew[i][j]  =  rnew[i][j] - round(rnew[i][j]/L)*L; // making boundaries
// apply periodic boundary condition: if rnew[i][j]<L/2 -> rnew[i][j]=rnew[i][j]
//                                  : if rnew[i][j]>L/2 -> rnew[i][j]=rnew[i][j]-L
                           }
                     }
              }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void printtest(){

fprintf(traj,"%d",3231231);
}

void note(){
    int i,j;

    if(t%dat == 0){
        fprintf(traj, "%d\n%s\n",N,"empty");
        for(i=0;i<N;i++){
            fprintf(traj,"%d\t",1);
            for(j=0;j<Dim;j++){
                    fprintf(traj,"%f\t",rnew[i][j]);
                               }
                    fprintf(traj,"\n");
                         }
                   }
            }






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// function to calculate distance R[a][b] between particles 'a' and 'b' and force which acts between them
void CalDist(){
    int i;
            Rsq[a][b]  = 0;
    for(i=0;i<Dim;i++){
            D[a][b][i] = rnew[a][i]  - rnew[b][i];
            D[a][b][i] = D[a][b][i]  - round(D[a][b][i]/L)*L;
            Rsq[a][b]  = Rsq[a][b]   + D[a][b][i]*D[a][b][i];
                       }
    
            R[a][b]  = sqrt(Rsq[a][b]);
            R[b][a]  = R[a][b];
               }





            
// Calculate potential depth between a and b//
// given SpVec, CIN, CA, let's calculate potential

/* if "a" and "b" are connected or not => if "CA[a]=b (CA[b])=a" or not */
/* if particles are inside the same cluster or not => if "CIA[a]=CIA[b]" or not" */




void CalEdepth(){
    /*--- if two particles are next to each other ---*/
    if(abs(a-b)==1){
        Edep = 5.0;         // for particles to be chain
                    }
    
    // --- if two "non neighbor" particles with SV[a] and SV[b] are interacting with each other they interact with ---
     else {
        Edep = IM[SV[a]][SV[b]];
         }
    // Edep = IM[SV[a]][SV[b]];
                 }


// CalForce Calculate Forces acting between particle "a" and "b" //
void CalForce(){
    int i;
    if(Edep > 0.5){
                if(Roff < R[a][b]){
                                F[a][b] = 0;
                                   }
                
                else if (1< R[a][b]){
                    F[a][b] = (-2.0*Rho*Edep*E)*(exp(-Rho*(R[a][b]-1.0))*(1-exp(-Rho*(R[a][b]-1.0)))-((exp(4.0)-1.0)/exp(8.0)));
                                     }
                
                else if (R[a][b]<=1){
                    F[a][b] = (-2.0*Rho*Edep*E)*(pow(M, 2.0)*Rho*(R[a][b]-1.0)-((exp(4.0)-1.0)/exp(8.0)));
                                     }
        
                   }
            
            
              else{
                     if(1. < R[a][b]){
                                F[a][b] = 0;
                                      }
                
                else if (R[a][b]<=1 ){
                    F[a][b] = (-2.0*Rho*E)*(pow(M, 2.0)*Rho*(R[a][b]-1.0)-((exp(4.0)-1.0)/exp(8.0)));
                                      }
                  
                   }
    
    for(i=0;i<Dim;i++){
            f[a][b][i] =  F[a][b]*(D[a][b][i]/R[a][b]);
            f[b][a][i] =  -f[a][b][i];
                       }
    
                 }




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// function to calculate forces when we are not dealing with the verlet algorithm

void SumForces(){
    int i,j,k,l,m;
    for(i=0;i<N;i++){
        for(j=0;j<Dim;j++){
            fpnew[i][j]=0.0;
            
        }
    }
    
    
    for(k=0;k<N;k++){
        for(l=0;l<N;l++){
            for(m=0;m<Dim;m++){
                fpnew[k][m]=fpnew[k][m]+f[k][l][m];
            }
        }
    }
}






//%%%%%%%%%%%%%%%%% Generate Uniforml Random Numbers %%%%%%%%%%%%%%%%%%%//
// generate uniformly distributed random numbers between 0 to 1 (can we improve this random number generation??) //
double Uniform(){
    static int x=10;
    int i=1103515245, j=12345, k=2147483647;
    x = (i*x + j)&k;
    
    return ((double)x+1.0) / ((double)k+2.0);
}

//%%%%%%%%%%%%%%%%% Generate Normally Distributed Numbers %%%%%%%%%%%%%%%%%//
// return numbers with random normal distribution //
double RandNormal(){
    double x=sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
    return x;
}




void AdMatrix(){
    if(t%dat == 0){
        BoNT = 0;
        BoNP = 1;
        int a,b,c,d;
        for(a=0;a<N;a++){
            BoN[a] = 0;
            for(b=0;b<N;b++){
                if (R[a][b] < BoL && a!=b){ // if distance between particle "a" and "b" is closer than distance "BoL", we substitute "1" in the AM[a][b]
                    AM[a][b] = 1;
                }
                else{
                    AM[a][b] = 0;
                }
                BoN[a] = BoN[a] + AM[a][b]; // calculate total number of bonds, particle "a" has.
            }
            BoNT = BoNT + BoN[a];  // calculate total number of bonds *2
            BoNP = BoNP * BoN[a];  // calculate all product of bonds
        }
        
        for(c=0;c<N;c++){          // this is roop to note adjacency matrix, if the cluster has more than seven bonds and no particle apart //
            for(d=0;d<N;d++){
                if(/*BoNT > 13 &&*/ BoNP !=0 ){fprintf(Amatrix,"%d\t",AM[c][d]);

            }
                else{fprintf(Amatrix, "0\t");}
            }
            fprintf(Amatrix,"\n");
        }
        fprintf(Product,"%d\n",BoNP);
    	countStruct(BoNP);
        // if(BoNP != 0){
        //     TApart = 0;
        // }
        
    }
}

void countStruct(int BoNP){
         if (BoNP==20736){
            count_a++;
         }
         if (BoNP == 20250){
            count_b++;
        }
        if (BoNP == 21600){
            count_c++;
        
        }
        if (BoNP == 25600){
            count_d++;
        }
        
        if (BoNP == 24000){
            count_e++;
        }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Renew(){
    int i,j;
    for(i=0;i<N;i++){
        for(j=0;j<Dim;j++){
        rold[i][j]   = rnew[i][j]  ;
        fpold[i][j]  = fpnew[i][j];
                          }
                    }
             }


// void runMain(char filename[]){
// 	char trajName[100]; char amatrixName[100]; char productName[100];
// 	strcpy(trajName,filename); 
// 	strcat(trajName, "trajChain.xyz");
// 	strcpy(amatrixName,filename); 
// 	strcat(amatrixName, "Amatrix.txt");
// 	strcpy(productName,filename); 
// 	strcat(productName, "Product.txt");

// 	printf("%s\t%s\t%s\n",trajName,amatrixName,productName);
// traj = fopen(trajName,"w");
// Amatrix = fopen(amatrixName,"w");
// Product = fopen(productName,"w");
// 		int i;
// 	//Dif = dif;
// 	// for ( i = 0; i<6; i++){
// 	// 	Dif+=(0.02);
//         InitialSet();
//     /* time integration loop */
//         t=0;
//     while(t<time){    t = t + 1;
//         newcal();               // given Xold , FXpold, FbX (so as Y) calculate Xnew based on forward Euler time-step algorithm
//         note();                 // write down Xnew we've got on a file

// 		//fprintf(traj,"%d\n",1);
//     /*--- calculations of pairwise quantities (distance, interactions etc...) ---*/
//             for(a=0;a<N-1;a++){
//                 for(b=a+1;b<N;b++){
//                     CalDist();              // Given configuration of particles calculate DX, R[a][b] between each pair of particles.
//                     CalEdepth();            // calculate depth of potential energy depending on its situation
//                     CalForce();             // calculate force acting between pairs of particle
                    
//                 }
//             }
       
//         SumForces();
//         AdMatrix();
//         Renew();
        
//      //  }
//      // printf("%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
// }
//     printf("%f\t%d\t%d\t%d\t%d\t%d\n",Dif,count_a,count_b,count_c, count_d,count_e);
// }






































