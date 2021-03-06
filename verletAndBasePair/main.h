/*--- Brownian Dynamics Simulation of self-replicating colloidal cluster with the Morse potential using the forward Euler time-step---*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
/*--- Define constants ---*/
#define Dim   3           // dimension of the system
#define h     0.000002    // time step of this simulation (we should also think about what is the possible largest time step, which doesn't change physics)
#define time  100000    // total time-steps of this simulation
#define N     18          // total number of particles involved in this simulation is N*(N-1)+12
#define L     1000000000.0// size of the simulation box (we use square(for 2D) or cube(for 3D) "L" each side.)
#define E     7.0         // "unit" depth of the Morse potential (interaction potential)
#define dat   1000        // time interval: taking trajectory data
#define Vtime 100         // time interval: renewing the verlet list
#define M     2           // coefficient for harmonic potential
//#define Dif   1.2         // Diffusion coefficient
#define Rho   30.0        // probability measure (described in Miranda PNAS paper)
#define BoL   1.05        // we say bond is created when distance between two particles is less than BoL
#define IE0   0.0         //
#define Rv    2.0         // radius for verlet algorithm

/*----Declare functions----*/
void    InitialSet();
void    IniSV();          // function to initially set species vector
void    IniIM();          // function to set Interaction Matrix
void    IniConf();        // function to set initial configuration
void    newcal();         // calculate positions and velocities of next time step for every paritlces
void    note();           // every 'dat' time, write down informations of positions of every particles
void    NewVerletList();  // create new verlet list
void    CalDist();        // calculate distance between every possible pairs of particles
void    CalEdepth();      // calculate depth of potential energy given SV,CIN,CA
void    CalForce();       // calculate force acting between every possible pairs
void    SumForces();      // function to sum-up all forces acting between particles
void    nonVSumForces();  // SumForces for non verlet case
void    CalRandomForce(); // function to calculate gaussian distributed random force
void    Renew();          // put new variables into old variables preparing for integration at next time step
/*-----Files------*/
FILE    *traj;
FILE 	*bond;
FILE    *state;

double Dif;
int     SV[N];            // species vector (vector to remember which particle is assigned to be which species)
int     NumCl;            // Number of Clusters => total number of clusters (parents or catalysts)
int     t;                // global variable to make time roop
int     a,b;
double  Uniform();        // generate uniform random numbers distributing between 0 to 1
double  RandNormal();     // generate random normal distributed numbers with mean=0, sigma(variance)=1
double  Edep;
double  VL[N][100];       // matrix to store verlet list
double  rnew[N][Dim];
double  rold[N][Dim];     // positions of particles(notation: -old=> input for integration,    -new=> output of integration)
double  fpnew[N][Dim];
double  fpold[N][Dim];    // sum of potential forces acting on single particle
double  F[N][N];
double  f[N][N][Dim];     // each force acting between a pair of particle
double  fb[N][Dim];       // random white noise force
double  D[N][N][Dim];     // Di[a][b][0]=rnew[a][0]-rnew[b][0]
double  Rsq[N][N];        // R square
double  R[N][N] ;         // sqrt(DX^2+DY^2)
double  IM[9][9];         // Interaction Matrix
double  Roff ;            // Roff(Cut-off distance of pariwise force)   Roff  = 1+(4/30) = 1.133
int     list[N][N];       // verlet list -> list[i][] represents PIN of particles inside verlet radius
int     nlist[N];         // nlist remembers how many particles are within the verlet radius


















/*--- below are function you may use in future ---*/

void    IniCIN();         // function to set initial cluster identification number
void    IniCA();          // function to set initial connection array
void    ModifySV();       // modify species vector when type 0 particle attaches to cluster
void    ModifyCA();       // modify Connection vector
int     CIN[N];           // cluster identification number (if a particle is not inside a cluster, -1)
int     CA[N];            // connection array
int     NCA[6];           // remember 6 monomers attached to a parent and catalyst
