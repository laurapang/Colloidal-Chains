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

//This function sets initial conditions of the system.
void InitialSet(){
   IniConf();                 
   Roff = (1+(4/Rho))*diam; //defining the cut-off length (above which there is no interaction) based on Holmes-Cerfon paper, section 5
   IniSV();                 
   IniIM();
}
   
//initializes the file name of the output
//DURING COMPILATION, when user runs "./EXEFILE.o OUTPUT-TEST 1" in terminal, output files will be:
//OUTPUT-TESTtrajChain.xyz and OUTPUT-TESTstate.txt
void openFiles(char* filename){
   //runMain(fileName);
 
   char trajName[100]; 
   char amatrixName[100];  
   //char bondName[100]; 
   char stateName[100];
   strcpy(trajName,filename); 
   strcat(trajName, "trajChain.xyz");
   strcpy(stateName,filename); 
   strcat(stateName, "state.txt");
   //strcpy(bondName,filename); 
   //strcat(bondName, "bond.txt");
 
   //printf("%s\t%s\t%s\n",trajName, stateName, bondName);
   traj = fopen(trajName,"w");
   stateFile = fopen(stateName,"w");
   //bond = fopen(bondName,"w");
}
//-----------------------------------------------------------
//This function constructs a Species Vector
void IniSV(){
   int i;
   // for (i=0; i<N;i++)
   //    SV[i]=0;
   // for (i=0;i<9;i++){
   //    SV[i]=i;
   //    SV[i+9]=8-i;
   // }
   char seq[] = "GCGTTGCTTCTCCAACGC";
   for(i=0;i<sizeof(seq);i++){
       char base =seq[i];
       switch(base){
           case 'A': SV[i]=0; break;      //converting DNA bases into a numbered species
           case 'T': SV[i]=1; break;
           case 'C': SV[i]=2; break;
           case 'G': SV[i]=3; break;
       }
   }   
   //for (i = 0; i<N;i++){
   //    printf("%d  ",SV[i]);
   //}printf("\n");
   //printf("%s\n",seq);
}
void IniIM(){
   //SINGLE SPECIES
   //IM[0][0]=1;

   //4 SPECIES LIKE DNA
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

void IniConf(){
   int i,j;
   for(i=0;i<N;i++){
      rold[i][0]=(double)i*diam;
      for(j=1;j<Dim;j++){
         rold[i][j]=(double)0.0;
      }
   }

   // for (i=9; i<18; i++){
   //    rold[i][0]=(double)(17-i)*diam;
   //    rold[i][1]=1*diam;
   // }
}
    
void newcal(){
   int i,j;
   for(i=0;i<N;i++){
      for(j=0;j<Dim;j++){
         rnew[i][j] = rold[i][j] + (fpold[i][j]*Kb*T*h)/(zi*diam) + sqrt(2.0*Dif*h)*RandNormal();
         rnew[i][j] = rnew[i][j] - round(rnew[i][j]/L)*L;
      }
   }
}

void note(){
   int i,j;
   if(t%(dat*10) == 0){// && t>200000000){
      fflush(stateFile);
      fprintf(traj, "%d\n%s\n",N,"empty");
      for(i=0;i<N;i++){
         fprintf(traj,"%d\t",1);
         for(j=0;j<Dim;j++){
            fprintf(traj,"%.10lf\t",rnew[i][j]/diam);
         }
         fprintf(traj,"\n");
      }
   }
}

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
// if "a" and "b" are connected or not => if "CA[a]=b (CA[b])=a" or not 
// if particles are inside the same cluster or not => if "CIA[a]=CIA[b]" or not" 

//------------------------------------------------------
//This function sets the depth of the Morse potential. 
void CalEdepth(){
   if(abs(a-b)==1){                     //If two particles are neighbors, set the potential depth to 5.0
      Edep = 5.0;
   }
   else {
      if (CA[a]==b||CA[a]==-1)          //
         Edep = 2.5*IM[SV[a]][SV[b]];
      else
         Edep = 0;
   }
}
//calculates the forces between 2 given particles, a and b
void CalForce(){
   int i;
   if(Edep > 0.5){  //if the energy depth is large enough (low Edep is no attraction)
      if(Roff/diam < R[a][b]/diam){ //if the Radius between a and b is too large (greater than the cut off distance Roff), there is no force
         F[a][b] = 0;
      }
      else if (1 < R[a][b]/diam){ //if the Radius is close enough for a force to be applied, calculate with modified morse potential equation
          F[a][b] = (-2.0*Rho*Edep*(E/(Kb*T)))*(exp(-Rho*(R[a][b]/diam-1.0))*(1-exp(-Rho*(R[a][b]/diam-1.0)))-((exp(4.0)-1.0)/exp(8.0)));
      }
      else if (R[a][b]/diam <= 1){ //if the Radius is too close (R<1 means that the particles are overlapping), particles need to repel to prevent overlap
          F[a][b] = (-2.0*Rho*Edep*(E/(Kb*T)))*(pow(M, 2.0)*Rho*(R[a][b]/diam-1.0)-((exp(4.0)-1.0)/exp(8.0)));
      }
   }

   else{ //if there is no attraction between the particles
      if(1 < R[a][b]/diam){
         F[a][b] = 0;
      }
      else if (R[a][b]/diam <= 1){ //if the particles are overlapping, they also need to repel
         F[a][b] = (-2.0*Rho*(E/(Kb*T)))*(pow(M, 2.0)*Rho*(R[a][b]/diam-1.0)-((exp(4.0)-1.0)/exp(8.0)));
      }
   }

   for(i=0;i<Dim;i++){
      f[a][b][i] =  F[a][b]*(D[a][b][i]/R[a][b]);
      f[b][a][i] =  -f[a][b][i];
   }
}
//OLD METHOD (NOT USED) that sums forces NOT based on the verlet list
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
//NEW method that sums forces based on verlet list
void SumForcesV(){
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


//Generates a uniformly distributed random number
double Uniform(){
   static int x=10;
   int i=1103515245, j=12345, k=2147483647;
   x = (i*x + j)&k;
   return ((double)x+1.0) / ((double)k+2.0);
}
//generates a number within normal distribution curve
double RandNormal(){
   double x=sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
   return x;
}
//hard-coded state function for the OxDNA hairpin example
void State(){
    if(t%dat == 0){
       fflush(traj);
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
        fprintf(stateFile,"%d\t",tot);
        int k;
        for (k=0;k<N;k++){
           fprintf(stateFile,"%d\t",CA[k]);
        }
        fprintf(stateFile,"\n");
        //fprintf(stateFile,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",t, tot, count_a, count_b,count_c,count_d,count_e,count_f);
        
    }
}
//updates the variables for next calculation: the old position is set to the "new" position and old force is set to "new" force
void Renew(){
   int i,j;
   for(i=0;i<N;i++){
      for(j=0;j<Dim;j++){
         rold[i][j] = rnew[i][j]  ;
         fpold[i][j] = fpnew[i][j];
      }
   }
}

//calculates verlet list
//nlist[i] is the number of particles within a certain radius to "i"
//in example, if i is near particle a and b
//then nlist[i] = 2, and list[i][0]=a and list[i][1]=b
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
//THIS FUNCTION IS FLAWED*****************************
//ensures that the particle should only be bonded to one other particle at a time
//refer to bintree.c and bintree.h for future work perhaps - bintree.c can generate sorted array
//can sort the radius between a and b, and then designate connectivity starting based on the smallest distances
void oneBond(){
  for (a=0; a<N;a++){
            double minBond = 1.05*diam;
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
        //i commented out this part of the code to save time
        //if (t%dat==0){
        //fprintf(bond,"2Bonded: ");
        //int cc;
        //for (cc = 0; cc<N;cc++){
        //    fprintf(bond,"%d\t%lf\t",CA[cc],R[cc][CA[cc]]/diam);
        //}fprintf(bond,"\n");
        //

        //}
}
void IniCA(){
    int i;
    for(i=0;i<N;i++){
        CA[i]= -1; /* -1 stands for no connection */ // if a particle has bond between "monomer" and "parents inside cluster" CA[a]= "PIN of the pair"
    }
}
