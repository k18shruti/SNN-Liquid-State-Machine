//generate Poisson distributed spike times for input and desired
//lambda_ip: avg spike rate in the input spike stream
//lambda_des: avg spike rate in the desired spike stream
//T: length/duration of the spike streams
//dt: time step used for simulation

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"snn.h"

#define L 1000

//global variables to store the input and desired spike streams
extern int inp_spk[N_pattern][M], des_spk[N_pattern][M];

int spk_gen(int lambda_ip, int lambda_des, double T, int seed, double dt)
{
    float des_time[N_pattern][max_spks], inp_time[N_pattern][max_spks];
    float success_ip, success_des, unif_ip[L], unif_des[L];
    float poisson_ip[L], poisson_des[L], t_ms_ip[L], t_ms_des[L];
    int r[N_pattern], p[N_pattern];		//actual no. of spikes in inp and desired
    int i, j, k, l, rand_ip, rand_des;
    float sum_tms_ip=0, sum_tms_des=0;

    success_ip = (float)lambda_ip*T;
    success_des = (float)lambda_des*T;
    //actual no. of spike times to be created:
    rand_ip = success_ip*10;
    rand_des = success_des*10;
    printf("T=%f\n",T);
    printf("lambda_des=%d, lambda_ip=%d, max_spks=%d\n",lambda_des, lambda_ip, max_spks);
	
    srand(seed);

    for(k=0; k<N_pattern; k++) {
	r[k]=0; p[k]=0;
	sum_tms_des=0; sum_tms_ip=0;
	while(r[k]<1 && p[k]<1) {	  //while there are no desired no. of spikes in the streams
	    for(j=0; j<max_spks; j++) {
		des_time[k][j]=0;
		inp_time[k][j]=0;
	    }

            for(i=0; i<L; i++) {
		//for iput spikes:
	        unif_ip[i] = (float)rand()/(float)RAND_MAX;
  	        poisson_ip[i] = -log(1-unif_ip[i])/lambda_ip;
	        t_ms_ip[i] = round(poisson_ip[i]*1000)/1000;	//precision upto 0.1ms
		//for desired spikes:
		unif_des[i] = (float)rand()/(float)RAND_MAX;
		poisson_des[i] = -log(1-unif_des[i])/lambda_des;
		t_ms_des[i] = round(poisson_des[i]*1000)/1000;
            }
	
	    inp_time[k][0] = t_ms_ip[0] + 1E-3;

	    //1st spike to have non-zero time
	    i=0;
	    j=1;
	    while(inp_time[k][j-1]<T) {
		i++;
		if(t_ms_ip[i]<Rp) {
		    continue;}
		else {
		    //sum up the times
		    sum_tms_ip = 0;
		    for(l=0; l<i; l++) {
			sum_tms_ip = sum_tms_ip + t_ms_ip[l];
		    }


		    inp_time[k][j] = sum_tms_ip;
		    r[k]++;
		    j++;
		}
	    }	

            //create the desired spike times:
	    i=0; j=1;

	    while(sum_tms_des < (inp_time[k][0]+Rp)) {
		i++;
		sum_tms_des = 0;
		for(l=0; l<i; l++) {
		    sum_tms_des = sum_tms_des + t_ms_des[l];
		}
	    }
	    des_time[k][0] = sum_tms_des;

	    while(des_time[k][j-1]<T) {
		i++;
		if(t_ms_des[i]<Rp) {
		    continue;
		}
		else {
		    sum_tms_des=0;
		    for(l=0; l<i; l++) {
			sum_tms_des = sum_tms_des + t_ms_des[l];
		    }
		    des_time[k][j] = sum_tms_des;
		    p[k]++;
		    j++;
		}
	    }
	}	//end of while loop
	printf("p=%d, r=%d, k=%d\n",p[k],r[k],k);
    }	//end of k loop

    

    for(i=0; i<N_pattern; i++) {
	k=0; l=0;
	for(j=0; j<M; j++) {
	    for(k=0; k<r[i]; k++) {
	    if(j==(int)(inp_time[i][k]/dt)) {
		inp_spk[i][j]=1;
	    }
	    }
	
	    for(l=0; l<p[i]; l++) {
	    if(j==(int)(des_time[i][l]/dt)) {
		des_spk[i][j]=1;
	    }
	    }
	}
    }	
}	//end of function
