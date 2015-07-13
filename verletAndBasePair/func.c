/*---Brownian Dynamics simulation of 2D colloidal cluster with the Morse potential and the forward Euler time-step---*/
#include "main.h"

void InitialSet(){
// first, we create and open files to note our result
// traj     =   fopen("0713_07trajChain.xyz","w");
// bond = fopen("0713_07bond.txt","w");
// state = fopen("0713_07state.txt","w");
/* set initial conditions */
IniConf();                  // IniCond(Initial Condition) => here we set initial condition of particles
Roff = 1+(4/Rho) ;          // Roff(Cut-off distance of pariwise force)   Roff  = 1+(4/30) = 1.133
IniSV();                    // Species vector
IniIM();                    // set 9 by 9 interaction matrix
}

void openFiles(char* filename){
  //runMain(fileName);

 char trajName[100]; char amatrixName[100]; char bondName[100]; char stateName[100];
  strcpy(trajName,filename); 
  strcat(trajName, "trajChain.xyz");
  strcpy(stateName,filename); 
  strcat(stateName, "state.txt");
  strcpy(bondName,filename); 
  strcat(bondName, "bond.txt");

  printf("%s\t%s\t%s\n",trajName, stateName, bondName);
  traj = fopen(trajName,"w");
  state = fopen(stateName,"w");
  bond = fopen(bondName,"w");
}


//%%%%%%%%%%%%%%%%% Initial Species Vector %%%%%%%%%%%%%%%%%//
// IniSpVec => Initial Species Vector; we assign initial "Species" for each particle.
void IniSV(){
    // int i;
    // for(i=0;i<N;i++){
    //     SV[i]=0;
    // }
  int i;
    char seq[] = "GCGTTGCTTCTCCAACGC";
    //char seq[] = "AAAAAAAAATTTTTTTTT";
    for(i=0;i<sizeof(seq);i++){
        char base =seq[i];
        switch(base){
            case 'A': SV[i]=0; break;
            case 'T': SV[i]=1; break;
            case 'C': SV[i]=2; break;
            case 'G': SV[i]=3; break;
        }
        //SV[i]=i%4;
    }   
    for (i = 0; i<N;i++){
        printf("%d  ",SV[i]);
    }printf("\n");
    //scp filename.txt lpang@login.rc.fas.harvard.edu:
    printf("%s\n",seq);
             }

//%%%%%%%%%%%%%%%%% Interaction Matrix %%%%%%%%%%%%%%%%%%%%//
/* Interaction Matrix between given pairs of particles */
void IniIM(){
      IM[0][0]=0; IM[0][1]=1; IM[0][2]=0; IM[0][3]=0;
               IM[1][1]=0; IM[1][2]=0; IM[1][3]=0;
                           IM[2][2]=0; IM[2][3]=1;
                                       IM[3][3]=0;
                 int i,j;
    for(i=0;i<8;i++){
        for(j=i+1;j<9;j++){
                   IM[j][i]=IM[i][j];
                           }
                     }
             }
void State(){
    if(t%dat == 0){
        int count_a=0;int count_b=0;int count_c=0;int count_d=0;int count_e=0;int count_f=0;
        if (CA[0]==17){
            count_a++;
        }
        if (CA[1]==16){
            count_b++;
        }
        if (CA[2]==15){
            count_c++;
        }
        if (CA[3]==14){
            count_d++;
        }
        if (CA[4]==13){
            count_e++;
        }
        if (CA[5]==12){
            count_f++;
        }
        int tot = count_a+count_b+count_c+count_d+count_e+count_f;
        fprintf(state,"CAstate: %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",t, tot, count_a, count_b,count_c,count_d,count_e,count_f);
        
    }
}
void oneBond(){
  for (a=0; a<N;a++){
            double minBond = 1.05;
            //fprintf(bond,"a= %d\t: nlist[a]=%d\t",a,nlist[a]);
            int c; 
            for(c=0;c<nlist[a];c++){            
            b=list[a][c];
            //fprintf(bond, "%d\t%lf\t",b,R[a][b]);
                if (R[a][b]<minBond&&a!=b&&!(b==(a-1)||b==(a+1))){
                    minBond=R[a][b]; CA[a]=b;
                }
            }

        }
        
        int j; int k;
        for (j=0; j<N; j++){
            k = CA[j];
            if (k!=-1){
                if (CA[k]!=-1){
                if (CA[k]!=j){
                    //fprintf(bond,"%d\t%d\t%d\t%d\t%d\n",t,j,CA[j],k,CA[k]);
                    double amin = R[j][CA[j]];
                    double bmin = R[k][CA[k]];
                    if (amin<bmin)
                        CA[k]=-1;
                    else
                        CA[j]=-1;
                }}
            }
        }
        if (t%dat==0){
        fprintf(bond,"2Bonded: ");
        int cc;
        for (cc = 0; cc<N;cc++){
            fprintf(bond,"%d\t%lf\t",CA[cc],R[cc][CA[cc]]);
        }fprintf(bond,"\n");
        

        }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// IniConf => Initial Configuration; function to set initial condition //
void IniConf(){
/* first we set initial configurations, so that particles make chain */
    int i,j;
    // for(i=0;i<N;i++){
    //     rold[i][0]=i;
    //     for(j=1;j<Dim;j++){
    //         rold[i][j]=0;
    //     }
    // }


    // for(i=0;i<N;i++){
    //     for(j=0;j<Dim;j++){
    //         rold[i][j]=0;
    //     }
    // }
    // rold[1][1]=-1; rold [2][1]=-2;
    // rold[3][0]=1; rold[3][1]=-2; rold[4][0]=1; rold[4][1]=-1; rold[5][0]=1;
    // rold[6][0]=2; rold[7][0]=2; rold[7][1]=-1; rold[8][0]=2; rold[8][1]=-2;

    //int i; int j;
        for(i=0;i<N;i++){
        //rold[i][0]=i;
        for(j=0;j<Dim;j++){
            rold[i][j]=0;
        }
    }
    
    for (i = 0; i<N/2; i++){
      rold[i][0]=i;
        rold[i][1]=0;
      
    }
    for (i=N/2; i<N; i++){
      rold[i][0]=(N-1-i);
        rold[i][1]=1;
      
    }
    }






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// given rold and fpold, return rnew
void newcal(){
    
    int i,j;
    for(i=0;i<N;i++){        // calculate positions of next time step
        for(j=0;j<Dim;j++){
        rnew[i][j]  =  rold[i][j]  + fpold[i][j]*h + sqrt(2.0*Dif*h)*RandNormal();  // over damped Langevin equation
        rnew[i][j]  =  rnew[i][j] - round(rnew[i][j]/L)*L;
// apply periodic boundary condition: if rnew[i][j]<L/2 -> rnew[i][j]=rnew[i][j]
//                                  : if rnew[i][j]>L/2 -> rnew[i][j]=rnew[i][j]-L
                           }
                     }
              }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
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



// NewVerletList => Given rnew[PIN][Dim] of every particles, create verlet list list[][]
void NewVerletList(){
    int i,j,k,l;
    for(i=0;i<N;i++){
        nlist[i]=0;
    }
    for(j=0;j<N-1;j++){
        for(k=j+1;k<N;k++){
            Rsq[j][k] = 0;
            
            for(l=0;l<Dim;l++){
                D[j][k][l] = rnew[j][l]  - rnew[k][l];
                D[j][k][l] = D[j][k][l]  - round(D[j][k][l]/L)*L;
                
                
                Rsq[j][k]  = Rsq[j][k] + D[j][k][l]*D[j][k][l];
            
            }
            R[j][k]  = sqrt(Rsq[j][k]);
            
            R[k][j]  = R[j][k];
            if(R[j][k]<Rv){
                list[j][nlist[j]]=k;
                list[k][nlist[k]]=j;
                nlist[j]=nlist[j]+1;
                nlist[k]=nlist[k]+1;
            }
            
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
        Edep = 5.0;
                    }
    
    /*--- if two "non neighbor" particles with SV[a] and SV[b] are interacting with each other they interact with ---*/
    else {
      if (CA[a]==b||CA[a]==-1)
        Edep = IM[SV[a]][SV[b]];
      else
         Edep=0;
           }
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
// function to calculate potential force Fp
void SumForces(){
    int i,j,k,l,m,n;
    for(i=0;i<N;i++){
            for(j=0;j<Dim;j++){
                fpnew[i][j]=0.0;
                              }
                     }

    for(k=0;k<N;k++){
        for(l=0;l<nlist[k];l++){
            n=list[k][l];
                for(m=0;m<Dim;m++){
                fpnew[k][m]=fpnew[k][m]+f[k][n][m];
                                  }
                         }
                     }
                 }




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// function to calculate forces when we are not dealing with the verlet algorithm

void nonVSumForces(){
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


























//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%For Future Use%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// Connection array remembers which attached monomers and particles inside clusters are connected. if there is no connection between attatched monomer and cluster it is -1//

void IniCA(){
    int i;
    for(i=0;i<N;i++){
        CA[i]= -1; /* -1 stands for no connection */ // if a particle has bond between "monomer" and "parents inside cluster" CA[a]= "PIN of the pair"
    }
}

//%%%%%%%% Initial Cluster Identification Number %%%%%%%%%%%//
// Cluster identification number remembers inside which cluster does a particles belongs to.
void IniCIN(){
    int i;
    for(i=0;i<N;i++){
        CIN[i]= 0; /* CIN[a]= "cluster number" when inside parent or catalysts, CIN[a]= -1 when it is monomer　*/
    }
}



// Species Vector is modified only when "type 0" attaches to an "unoccupied" particle inside cluster//

void ModifySV(){
    if(R[a][b]<=BoL){
        if(SV[a]==0/*a is type 0 <=> unoccupied */ && CIN[b]!= -1 /* b is inside cluster */ && CA[b] == -1 /* b is unoccupied */ ){
            if(SV[b]>4){
                SV[a]=SV[b]-4;
            }
            else{
                SV[a]=SV[b]+4;
            }
            
            /* make the species for the one of complementally particle*/
        }
        else if(SV[b]==0 /* b is type 0 */ && CIN[a]!= -1 /* a is inside cluster */ && CA[a] == -1 /* a is unoccupied */ ){
            if(SV[a]>4){
                SV[b]=SV[a]-4;
            }
            else{
                SV[b]=SV[a]+4;
            }
            
            
            /* make the species for the one of complementally particle*/
        }
    }
}

//Modify Connection Array when "unoccupied" "monomer and cluster" comes togeter or "connected" particles fallen apart
void ModifyCA(){
    if(R[a][b]<=BoL){
        if(CA[a]==-1 && CA[b] ==-1/*both unoccupied*/){
            if(CIN[a]==-1 && CIN[b] != -1){
                CA[a] = b, CA[b]= a;
            }
            else if(CIN[b] ==-1 && CIN[a] != -1){
                CA[a] = b, CA[b]= a;
            }
        }
    }
    
    else if(R[a][b]>BoL && CA[a] == b /*particles are connected*/){
        CA[a] = -1, CA[b] = -1;
    }
}

































