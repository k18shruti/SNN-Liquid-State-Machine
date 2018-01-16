//header file

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define C 0.2E-6
#define C_op 1E-9
#define tau_op 1.5E-3
#define gL 1E-6
#define VT -55E-3
#define EL -62E-3
#define Rp 3E-3			//refractory period
#define tau 6E-3
#define X 10			//3D arrangement of the NMC block
#define Y 5
#define Z 20
#define N X*Y*Z                 //no. of neurons
#define simT  1.0		//Total simulation duration
#define M 10000			//no. of time steps
#define max_epochs 1000		//max no. of training iterations

#define max_spks 1667		//maximum no. of spikes within a given duration
#define N_pattern 1		//no. of trials/spike patterns

extern int spk_gen(int, int, double, int, double);
extern int nmc_str();

