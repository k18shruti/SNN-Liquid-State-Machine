//Function to generate the connectivity within the Neural Micro-circuit (NMC) block
//N: total number of neurons, within the NMC
//X, Y, Z: (in header snn.h) dimensions of the network

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include"snn.h"

extern int nmc_conn[N][N], nmc_conn_exc[N][N];
extern int op_conn[N], inp_conn[N];

int nmc_str()
{
    int i, j, k;
    int exc_nmc[N];
    double lambda=4.8, c, pr_conn; 
    int x1[N], yy[N], z1[N];
    static int xd[N][N], yd[N][N], zd[N][N];
    static int dist[N][N];

    //input connections; 30% prob of connections:
    for(i=0; i<N; i++) {
        inp_conn[i]=round((double)rand()/(double)RAND_MAX-0.2);
    }
    printf("created inp connections\n");

    for(k=0; k<Z; k++) {
	for(j=0; j<Y; j++) {
	    for(i=0; i<X; i++) {
		x1[i+(j*X)+(k*X*Y)]=i;
                yy[i+(j*X)+(k*X*Y)]=j;
                z1[i+(j*X)+(k*X*Y)]=k;
	    }
	}
    }
    
    for(k=0; k<N; k++) {
	for(j=0; j<N; j++) {
	    xd[k][j]=x1[k]-x1[j];
	    yd[k][j]=yy[k]-yy[j];
	    zd[k][j]=z1[k]-z1[j];
	}
    }
    
    for(k=0; k<N; k++) {
	for(j=0; j<N; j++) {
	    dist[k][j]=xd[k][j]*xd[k][j]+yd[k][j]*yd[k][j]+zd[k][j]*zd[k][j];
	}
    }

    for(k=0; k<N; k++) {
        exc_nmc[k] = 2*round((double)rand()/(double)RAND_MAX+0.3)-1;
        for(j=0; j<N; j++) {
            nmc_conn_exc[k][j]=exc_nmc[k];
        }
    }
    
    for(i=0; i<N; i++) {
	for(j=0; j<N; j++) {
	    if(exc_nmc[i]==1 && exc_nmc[j]==1)
		c=0.3;
	    else if(exc_nmc[i]==1 && exc_nmc[j]==-1)
		c=0.2;
	    else if(exc_nmc[i]==-1 && exc_nmc[j]==1)
		c=0.4;
	    else if(exc_nmc[i]==-1 && exc_nmc[j]==-1)
		c=0.1;
	    pr_conn = c*exp(-dist[i][j]/lambda);
	    if(i==j) {
		nmc_conn[i][j]=0;
	    } else {
  	        nmc_conn[i][j] = round((double)rand()/(double)RAND_MAX-(0.5-pr_conn));
	    }
	}
    }
    printf("created nmc connections\n");

    for(k=0; k<N; k++) {
        op_conn[k]=round((double)rand()/(double)RAND_MAX-(0.5-(0.75-0.05*exc_nmc[k])));
    }
    printf("created op connections\n");

return(0);
}

