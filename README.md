# SNN-Liquid-State-Machine
A simple implementation of Liquid State machine using Leaky integrate and fire spiking neurons


The folder consists of C-files for simulation of a spiking neural network, they can be compiled and run using gcc.

The neurons are simple leaky integrate-and-fire models.
The synaptic current is modeled by a single decaying expoential kernel.

For all simulations, the time step (dt) is kept at 0.1ms.

The main file is lsm_simln.c and 2 files for the functions used are spk_gen.c & nmc_str.c.

1. spk_gen.c: Generates the input and desired (used only when running a supervised learning algorithm) spike streams at the specified average rate.

2. nmc_str.c: Creates the connectivity matrix for the LSM based neural micro circuit.

lsm_simln.c: This simulates a network of N recurrently connected neurons. The header file (snn.h) has the network parameters which can be modified.
   N: no. of neurons in the network.
   simT: duration of simulation
   M: no. of time steps (simT/dt)
   X,Y,Z: no. of neurons in each of the 3 directions within the NMC block.

To compile: gcc lsm_simln.c spk_gen.c nmc_str.c -lm -o lsm
To run: ./lsm
