//07.20.15
//Ramin Khajeh
#include "main.h"

void InitialSet(){
   IniConf();
   Roff = (1+(4/Rho))*diam;
   IniSV();
   IniIM();
}
   
void IniSV(){
   int i;
   for (i=0;i<9;i++){
      SV[i]=i;
      SV[i+9]=8-i;
   }
}
void IniIM(){
   int i,j;
   for (i=0;i<9;i++){
      for (j=0;j<9;j++){
         if (i==j){
            IM[i][j] = 1;
         }
         else{
            IM[i][j] = 0;
         }
      }
   }
   IM[6][6]= 0;
   IM[7][7]= 0;
   IM[8][8]= 0;
}

void IniConf(){
   int i,j;
   for(i=0;i<N;i++){
      rold[i][0]=(double)i*diam;
      for(j=1;j<Dim;j++){
         rold[i][j]=(double)0.0;
      }
   }
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
   if(t%dat == 0){
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

void CalEdepth(){
   if(abs(a-b)==1){
      Edep = 5.0;
   }
   else {
      Edep = IM[SV[a]][SV[b]];
   }
}

void CalForce(){
   int i;
   if(Edep > 0.5){
      if(Roff/diam < R[a][b]/diam){
         F[a][b] = 0;
      }
      else if (1 < R[a][b]/diam){
          F[a][b] = (-2.0*Rho*Edep*(E/(Kb*T)))*(exp(-Rho*(R[a][b]/diam-1.0))*(1-exp(-Rho*(R[a][b]/diam-1.0)))-((exp(4.0)-1.0)/exp(8.0)));
      }
      else if (R[a][b]/diam <= 1){
          F[a][b] = (-2.0*Rho*Edep*(E/(Kb*T)))*(pow(M, 2.0)*Rho*(R[a][b]/diam-1.0)-((exp(4.0)-1.0)/exp(8.0)));
      }
   }

   else{
      if(1 < R[a][b]/diam){
         F[a][b] = 0;
      }
      else if (R[a][b]/diam <= 1){
         F[a][b] = (-2.0*Rho*(E/(Kb*T)))*(pow(M, 2.0)*Rho*(R[a][b]/diam-1.0)-((exp(4.0)-1.0)/exp(8.0)));
      }
   }

   for(i=0;i<Dim;i++){
      f[a][b][i] =  F[a][b]*(D[a][b][i]/R[a][b]);
      f[b][a][i] =  -f[a][b][i];
   }
}

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

double Uniform(){
   static int x=10;
   int i=1103515245, j=12345, k=2147483647;
   x = (i*x + j)&k;
   return ((double)x+1.0) / ((double)k+2.0);
}

double RandNormal(){
   double x=sqrt(-2.0*log(Uniform()))*sin(2.0*M_PI*Uniform());
   return x;
}

void State(){
   if (t%dat==0){
      b1=0;b2=0;b3=0;b4=0;b5=0;b6=0;
      if (R[0][17] <BoL){
         b1++;
      }
      if (R[1][16] <BoL){
         b2++;
      }
      if (R[2][15] <BoL){
         b3++;
      }
      if (R[3][14] <BoL){
         b4++;
      }
      if (R[4][13] <BoL){
         b5++;
      }
      if (R[5][12] <BoL){
         b6++;
      }
      bTotal = b1+b2+b3+b4+b5+b6;
      fprintf(stateFile,"%d\t%d\t%d\t%d\t%d\t%d\t%d\n",b1,b2,b3,b4,b5,b6,bTotal);
     }
}

void Renew(){
   int i,j;
   for(i=0;i<N;i++){
      for(j=0;j<Dim;j++){
         rold[i][j] = rnew[i][j]  ;
         fpold[i][j] = fpnew[i][j];
      }
   }
}
