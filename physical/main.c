//Date: 07.16.15 
//Description: dimensionalization
//Name: Ramin Khajeh
#include "main.h"

int main(){
   traj = fopen("myTraj01_low.xyz","w");
   InitialSet();
   t=0;
   note();//double aa = 1E-9; double bb = 2E-9; double add = aa+bb; double times = aa*bb; double divide = bb/aa;

   while(t<time){
      t = t + 1;//dt;
      newcal();
      note();      
      for(a=0;a<N-1;a++){
      	//printf("%d\t",a);
          for(b=a+1;b<N;b++){
              CalDist();
              //printf("%.2lf\t",R[a][b]);
              CalEdepth();
              CalForce();
          }
          //printf("\n");
      }
      SumForces();
      Renew();
   } 

   //printf("%lf\t%lf\t%.30lf\n",T,E,Dif);
   return 0;
}

