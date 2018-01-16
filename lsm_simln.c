//Simulation of the Recurrent spiking neural network
//Integrate and fire neuron C model:
//LIF model with refractory period of the neuron
//With spike train as stimulus to neurons
//N neurons, N input synapses
//N neurons, with connectivity info as per ReSuMe's NMC
//fanins stored as a linked list per neuron
//1 input signal triggering the nmc block
//Inlcude the excitatory and inhibitory weights
//Include axonal delays

//The input and desired spike streams are created inthe spk_gen 
//the probalistic connections of the network are also created in 
//C for different NMC sizes in function nmc_str

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
//#include<conio.h>	//Used while running on Dev CPP

#include"snn.h"

//structure to store the fanin info of each NMC neuron
struct node {
    int val;
    double w_nmc, I_syn, D_nmc;
    struct node *next;
}*start[N];
   
//function for time keeping
long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
    return elapsed;
}

int spk_gen(int, int, double, int, double);
int nmc_str();

//need to declare them as global as they were too large and gave segmentation fault
double Vm[N][M];
int inp_spk[N_pattern][M], Spk[N][M];
int Spk_op[M];
int Spk_op_saved[max_epochs][M];
int des_spk[N_pattern][M], e[M];
//ds to store the fanins of each of the 'N' nmc neurons
//array of linked lists
struct node *fanin[N], *current[N];
int nmc_conn[N][N], op_conn[N], nmc_conn_exc[N][N];
int inp_conn[N];
double w_NMC[N][N], w_IP[N];
double d_nmc[N][N], d_ip[N];
double D_ip[N];			//exp LUT for the delays
double del_w[N], ci[N][M], d_hat[N][M];
double norm_dh[M];

double Vm_op[M], w_op[N];
double tdiff[max_spks], op_time[max_spks], des_time[N_pattern][max_spks];		//stores spike times

double max_tdiff;
int iterations[N_pattern];

int main()
{
    int i=0, j,k=0;     //no. of neurons
    int training_on=1, p=0, q=0; //p: no. of o/p spikes q: no. of des spikes
    double T=simT;         //total sim time (in sec)
    double dt=0.0001;     //time step = 0.1ms;
    double k1, k2;              //RK method variables
    FILE *FS;              //files to store NMC spike times
    FILE *F_ispk;               //file to store input spike trains
    
    double Iapp[N], I_op=0, weights[N];
    double ref_time[N], decay, I_in_prev[N], I_nmc_syn[N];
    clock_t t1, t2, t3, t4, t5;
    long elapsed;
    int req_iter;
    double ref_time_op=0, decay_op;
    double r = 6.4006E-10, decay1;	//learning rate
    int lambda_ip=20, lambda_des=20;	//spike rates
    int seed =27, t=0;
       
    printf("Running sim for %d time steps\n", M);

    FS=fopen("nmc_spikes.dat","w");
    if (FS==NULL) {
        printf("Can't open nmc_spik.dat for writing\n");
    }
    
    F_ispk = fopen("inp_spks.dat","w"); // read mode
    if(F_ispk == NULL)
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    
    printf("done opening the files for write\n");

    srand(seed);

    //populate the nmc connectivity matrix and conn types
    //input and output conenctions:
    nmc_str();

    printf("Created the SNN connectivity\n");

    //Randomize the weights & delays:
    for(j=0; j<N; j++) {
	//w_ip in range 4e-9 to 8e-9
	w_IP[j] = 4E-9 + 4E-9*((double)rand()/RAND_MAX);
	d_ip[j] = ((double)rand()/RAND_MAX)*(2E-3);
	//input synapse delay in range:
	w_op[j]=0;
	for(i=0; i<N; i++) {
	    //w_NMC in range: 2e-9 to 2e-12, -1e-12 to -1e-9
	    //d_nmc in range: [0,2E-3]
	    if(nmc_conn_exc[j][i]==1) {
		w_NMC[j][i]=2E-12+((double)rand()/RAND_MAX)*1.998E-9;
	    }
	    else if(nmc_conn_exc[j][i]==-1) {
		w_NMC[j][i]=-1E-9+((double)rand()/RAND_MAX)*(-9.99E-10);
	    }

	    d_nmc[j][i]=2E-3*((double)rand()/RAND_MAX);
	}
    }

    decay = exp(-dt/tau);  //exp decay of syn current, use -lm for gcc compile!
    decay1 = exp(-dt/1E-3);	//exp decay of mem potential, tauL'=1ms
    
    decay_op = exp(-dt/tau_op);


    //clear ref_time
    for(j=0; j<N; j++) {
        ref_time[j]=0;
		weights[j] = 40*w_IP[j];
		D_ip[j] = exp(d_ip[j]/tau);
        I_in_prev[j]=0;
        I_nmc_syn[j]=0;
        Vm[j][0]=EL;
    }

    //get the fan in info per NMC neuron
    //store in a linked list
    for(i=0; i<N; i++) {  //for each col 
        for(j=0; j<N; j++) {    //for each row in that col
            if(nmc_conn[j][i]==1) {
                fanin[i]=(struct node *)malloc(sizeof(struct node));                                  
                fanin[i]->val=j;
                fanin[i]->I_syn=0;
                fanin[i]->w_nmc=50*w_NMC[j][i];	//weights	
				fanin[i]->D_nmc=exp(d_nmc[j][i]/tau);	//exp of delays
                fanin[i]->next=NULL;
                if(start[i]==NULL)         //start uninitialised
                {
                    start[i]=fanin[i];     //initialise start to 1st node
                    current[i]=fanin[i];   //current to 1st node
                }
                else
                {
                    current[i]->next=fanin[i];
                    current[i]=fanin[i];
                }
            }
        }
    }

    //read the input spikes and clear the nmc spikes
    for(i=0; i<M; i++) {
        //fscanf(F_ispk,"%d,",&inp_spk[0][i]);
	//fscanf(FD,"%d,",&des_spk[0][i]);
	Spk_op[i]=0;
        for(j=0; j<N; j++)
            Spk[j][i]=0;     
    }

    //////////////////////////////////////////////////////
    /////////////Generate the spike patterns//////////////
    // Practical cases: spikes would be read from a file, or a sensor
    // then need to change the file permissions to read for inp and des spike streams
    //p: non zeros in desired, r: non zeros in input spike train
    spk_gen(lambda_ip, lambda_des, T, seed, dt);
    //printf("inp_spk[190]=%d, des_spk[500]=%d\n",inp_spk[0][190],des_spk[0][500]);
    //this function creates the reqd no. of spike streams for
    //in and des (global variables) for the given duration 'T'	
    //////////////////////////////////////////////////////

    printf("Simulating with N=%d, for T=%f sec\n",N,T);

    for(t=0; t<N_pattern; t++) {
		printf("expt. no. %d\n",t);
        q=0;
    
        //clear the nmc and output spikes:
        for(i=0; i<M; i++) {
            Spk_op[i]=0;
            for(j=0; j<N; j++)
                Spk[j][i]=0;     
        }

        //clear the synaptic currents and neuron potentials:
        for(j=0; j<N; j++) {
            ref_time[j]=0;
            I_in_prev[j]=0;
            I_nmc_syn[j]=0;
            Vm[j][0]=EL;
            fanin[j]=start[j];
            while(fanin[j]!=NULL)
            {
                fanin[j]->I_syn = 0;
                fanin[j]=fanin[j]->next;
            }
        }

        //CPU time required for computation
        t1 = clock();
        printf("Total no. of time steps in simulation= %f\n",(T/dt));

        i=0;
        ////START THE SIMULATION
        //would result in next mem([j+1][0] location being corrupted with (i+1)th value, at the end of i
        //NMC simulation:
        while(i<(M-1)) {
            for(j=0; j<N; j++) {
                if (ref_time[j]<i) {
                    ref_time[j]=0;     //jth neuron is not in refracory period               
                }
            
            //Compute the synaptic Current for each neuron in the network
            //composed of 2 components: input spikes and recerrent synaptic current
            //1. input spikes
                if(inp_spk[t][i]==1) { 
                    I_in_prev[j] = (I_in_prev[j] + weights[j]*D_ip[j])*inp_conn[j];
                }
                
                I_in_prev[j]=I_in_prev[j]*decay;

                //2. current due to synaptic spikes, from within the nmc
                fanin[j]=start[j];                     
                while(fanin[j]!=NULL) {
                    k=fanin[j]->val;      
                    if(Spk[k][i]==1) {
                        fanin[j]->I_syn = fanin[j]->I_syn + fanin[j]->w_nmc*fanin[j]->D_nmc;
                        I_nmc_syn[j] = I_nmc_syn[j] + fanin[j]->I_syn;
                    }
                    
                    fanin[j]->I_syn = fanin[j]->I_syn*decay;
                    fanin[j]=fanin[j]->next;
                }
                I_nmc_syn[j] = I_nmc_syn[j]*decay;
                Iapp[j]=I_in_prev[j] + I_nmc_syn[j];
                       
                
                //RK 2nd order method of solution for membrane potential
                k1 = (-gL*(Vm[j][i]-EL)+Iapp[j])/C;
                k2 = (-gL*((Vm[j][i]+dt*k1)-EL)+Iapp[j])/C;
                Vm[j][i+1] = Vm[j][i]+(dt*(k1+k2)/2)*(ref_time[j]==0);

                //check for spike
                if(Vm[j][i+1]>=VT) {
                    Vm[j][i+1]=EL;
                    ref_time[j]=i+Rp/dt;
                    Spk[j][i+1]=1;
                } else {
                    Spk[j][i+1]=0;   
                }            
                I_nmc_syn[j]=0;
            }   //close for(j) loop
            i++;
        }   //close while(i) loop
        printf("End of NMC simulation\n");

        t2 = clock();
        elapsed = timediff(t1,t2);
        printf("Elapsed time for NMC simulation: %ld ms\n", elapsed);

        //save the nmc spikes
        for(j=0; j<N; j++) {
            for(i=0; i<M; i++) {
                fprintf(FS,"%d,",Spk[j][i]);
            }
            fprintf(FS,"\n");
        } 

        t3 = clock();
        elapsed = timediff(t2,t3);
        printf("Elapsed time for NMC saving: %ld ms\n", elapsed);

    }	//end of t loop
    
    printf("Saving the input spikes\n");
    for(i=0; i<M; i++) {
	fprintf(F_ispk,"%d,",inp_spk[t-1][i]);
    }
    
    //getch(); //used when running the code on devC++
    return(0);
}
