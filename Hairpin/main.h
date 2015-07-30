#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define Dim   3
//#define tau   2.0E-9
#define diam  1.0E-9
//#define T     334.0
#define Kb    1.38E-23
//#define h     0.000002*tau
#define time  750000000// 1 microsecond
//#define zi    1.022*(Kb*T*tau)/(diam*diam) 
#define N     18    
#define L     64.0
//#define E     7.0*Kb*(298.0)
#define dat   20000
#define Vtime 100
#define M     2 
//#define Dif   1.0*((diam*diam)/tau)
#define Rho   30.0
#define BoL   1.05*diam
#define IE0   0.0
#define Rv    2.0*diam

double tau;
double T;
double h;
double zi;
double Dif;
double E;

void    InitialSet();
void    openFiles(char* filename);
void    SumForcesV();
void    IniSV();
void    IniIM();
void    IniConf();
void    newcal();
void    NewVerletList();  
void    note();
void    CalDist();
void    CalEdepth();
void    CalForce();
void    SumForces();
void    Renew();
void    State();
void    IniCA(); 
void	oneBond();
void	State();
int     CA[N];   
int     SV[N];
int     NumCl;
int     t; 
int     a,b;
int     b1,b2,b3,b4,b5,b6;
int     bTotal;
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
double  IM[9][9];
double  Roff;
int     list[N][N];
int     nlist[N];
int BoNT;
int BoN[N];
int BoNP;
int AM[N][N];

int count_a, count_b, count_c, count_d, count_e, count_f;
FILE *RandNum;
FILE *traj;
FILE *stateFile;
FILE *bond;

clock_t start;
clock_t stop; 
