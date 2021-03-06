#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define Dim   3
#define T     300.0
#define Kb    1.38E-23 
#define zi    4E-12//9.42E-12
#define diameter 1.0E-9
#define h     4.0E-15//2E-6 
#define time  1000000//2.0E-8
//#define dt    2E-9
#define N     6
#define L     64.0
#define E     7.0*(Kb*T)
#define dat   1000//2E-6
#define Vtime 100
#define M     2 
#define Dif   (Kb*T)/zi 
#define Rho   30.0E9
#define BoL   1.05E-9
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
