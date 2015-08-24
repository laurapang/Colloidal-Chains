//Date: 08.22.2015
//Description: Brownian Simulation of RNA self-assembly
//Harvard School of Engineering and Applied Sciences, Summer 2015 REU Program

//Names: 
//       Ramin Khajeh | raminkh@berkeley.edu
//       Laura Pan    | laura.pangx@gmail.com 
//Advisor: 
//       Professor Michael P. Brenner
//       Hidenori Tanaka

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define Dim   3 					//Dimension of simulation
#define diam  1.0E-9				//Length scale of the system (diameter of each sphere)
#define Kb    1.38E-23				//Boltzman's constant
#define h     4.0E-15				//Integration time-step in Langevin equation
#define sim_time  10000000000		//Simulation time (dimensionless)
#define N     18    				//Number of particles in the simulation
#define L     64.0 					//Size of the simulation box
#define dat   100000 				//Sampling frequency (sample every "dat" times)
#define Vtime 100 					//
#define M     2 					//The factor 'm' in parabolic potential (potential that is respnsible for avoiding overlaps)
#define Rho   30.0 					//
#define BoL   1.05*diam   			//Bond length (distane below which two particles are considered "bonded")
#define IE0   0.0  					//
#define Rv    2.0*diam  			//

double tau;							//Time scale of the system
double T;							//Temperature in Kelvin
double zi;							//Drag coefficient in Langevin equation
double Dif;							//Diffusion Coefficient 
double E;							//
void    InitialSet();				//Sets initial conditions of the system
void    openFiles(char* filename);	//?Opens necessary files
void    SumForcesV();				//?Sums forces 
void    IniSV();					//Fills in the species array SV[] 
void    IniIM();					//Assigns interaction amongst species
void    IniConf();					//Sets initial configuration of particles
void    newcal();					//
void    NewVerletList();  			//
void    note();						//Writes to trajectory file
void    CalDist();					//Calculates distance between parir of particles
void    CalEdepth();				//Calculates the potential depth amongst a pair of particles
void    CalForce();					//
void    SumForces();				//
void    Renew();					//?Update everything
void    State();					//Writes the state of 
void    IniCA(); 					//
void	oneBond();					//
int     CA[N];   					//
int     SV[N];						//
int     NumCl;						//
int     t; 							//
int     a,b;						//
int     b1,b2,b3,b4,b5,b6;			//
int     bTotal;						//
double  Uniform();					//
double  RandNormal();				//
double  Edep;						//
double  VL[N][100];					//
double  rnew[N][Dim];				//
double  rold[N][Dim];				//	
double  fpnew[N][Dim];				//
double  fpold[N][Dim];				//
double  F[N][N];					//
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
