#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define Dim   3
#define tau   2.0E-9
#define diam  1.0E-9
#define T     200.0
#define Kb    1.38E-23
#define h     0.000002*tau
#define time  100000000
#define zi    1.0*(Kb*T*tau)/(diam*diam) 
#define N     6    
#define L     64.0
#define E     7.0*Kb*T
#define dat   1000
#define Vtime 100
#define M     2 
#define Dif   1.0*((diam*diam)/tau)
#define Rho   30.0
#define BoL   1.05*diam
#define IE0   0.0
#define Rv    2.0

void    InitialSet();
void    IniSV();
void    IniIM();
void    IniConf();
void    newcal();
void    note();
void    CalDist();
void    CalEdepth();
void    CalForce();
void    SumForces();
void    Renew();

int     SV[N];
int     NumCl;
int     t; 
int     a,b;
double  Uniform();
double  RandNormal();
double  Edep;
double  VL[N][100];
double  rnew[N][Dim];
double  rold[N][Dim];
double  fpnew[N][Dim];
double  fpold[N][Dim];
double  F[N][N];
double  f[N][N][Dim];
double  fb[N][Dim];
double  D[N][N][Dim];
double  Rsq[N][N];
double  R[N][N];
double  IM[1][1];
double  Roff;
int     list[N][N];
int     nlist[N];
int BoNT;
int BoN[N];
int BoNP;
int AM[N][N];
FILE *RandNum;
FILE *traj;
