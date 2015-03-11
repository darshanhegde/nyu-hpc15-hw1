//
//  jacobi-mpi.c
//  Does Jacobi iterations to solve laplace equation in 1-d using MPI
//
//  Created by Darshan Hegde on 3/10/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


void print_solution(double* u_k, int n, int N, int mpirank, int mpisize){
    int i, index;
    for (i=1; i<n-1; i++) {
        index = (N/mpisize)*mpirank + i;
        fprintf(stderr, "value at %d is %f \n", index, u_k[i]);
    }
}

// Performs 1 iteration of jacobi iteration assuming
// all the required values are in place.
// modify kernel to suit the fact
void jacobi_kernel(double* u_k, int n, int N){
    double a_inv = 0.5/pow((double)N+1, 2);
    double h_inv = pow((double)N+1, 2);
    double u_prev=0, u_next=0;
    int i;
    u_prev = u_k[0];
    for (i=1; i<n-1; i++) {
        u_next = a_inv * (1+ h_inv*(u_prev+u_k[i+1]));
        u_prev = u_k[i];
        u_k[i] = u_next;
    }
}

// Performs Jacobi iterations using jacobi_kernel. It's main job is to
// perform boundary communications.
void jacobi_iteration(double* u_k, int n, int N, int mpirank, int mpisize, int max_iter){
    int iter=0, tag1=25, tag2=35;
    int prev, next;
    double first, last;
    
    MPI_Request reqs[4];
    MPI_Status stats[4];

    for (iter=0; iter<max_iter; iter++) {
        // Requires a barrier synchronization
        MPI_Barrier(MPI_COMM_WORLD);
        
        // exhange boundary values
        first = u_k[1];
        last = u_k[n-2];
        prev = mpirank-1;
        next = mpirank+1;
        if (mpirank == 0) {
            prev = mpisize-1;
            first = 0;
        }
        if (mpirank+1 == mpisize){
            next = 0;
            last = 0;
        }
        
        // Requires a barrier synchronization
        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Irecv(u_k, 1, MPI_DOUBLE, prev, tag1, MPI_COMM_WORLD, &reqs[0]);
        MPI_Irecv(u_k+n-1, 1, MPI_DOUBLE, next, tag2, MPI_COMM_WORLD, &reqs[1]);
        
        MPI_Isend(&first, 1, MPI_DOUBLE, prev, tag2, MPI_COMM_WORLD, &reqs[2]);
        MPI_Isend(&last, 1, MPI_DOUBLE, next, tag1, MPI_COMM_WORLD, &reqs[3]);
        
        MPI_Waitall(4, reqs, stats);
        
        // Requires a barrier synchronization
        MPI_Barrier(MPI_COMM_WORLD);
        
        // perform 1 step of jacobi
        jacobi_kernel(u_k, n, N);
        
        // Just to keep boundry values consistent
        if (mpirank == 0) {
            u_k[0] = 0;
        }
        if (mpirank+1 == mpisize) {
            u_k[n-1] = 0;
        }
    }
}

int main(int argc, char** argv){
    int mpisize, mpirank;
    long N, max_iter;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    
    if (argc != 3) {
        fprintf(stderr, "Please provide following arguments: <N: number of grid points> <max_iter: maximum number of iterations>\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    N = atol(argv[1]);
    max_iter = atol(argv[2]);
    
    if ((N%mpisize != 0) || (N>mpisize)) {
        fprintf(stderr, "N <number of grid points> must be a multiple of p <number of processors> \n");
    }
    
    int n = (N/mpisize)+2;
    double* u_k = (double*) malloc(n*sizeof(double));
    
    int i;
    for (i=0; i<n; i++) {
        u_k[i] = 0;
    }
    
    jacobi_iteration(u_k, n, N, mpirank, mpisize, max_iter);
//    print_solution(u_k, n, N, mpirank, mpisize);
    
    MPI_Finalize();
    
    return 0;
}