#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "util.h"

main(int argc, char *argv[])  {
    int mpisize, mpirank, N;
    int* message= (int*) malloc(500000*sizeof(int));
    int tag = 25;
    
    MPI_Status status;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    
    if (argc != 2) {
        fprintf(stderr, "Function needs number of passes as input arguments!\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    N = atol(argv[1]);
    int prev = (mpirank==0) ? mpisize-1 : mpirank-1;
    int next = (mpirank+1==mpisize) ? 0: mpirank+1;
    
    int iter;
    
    fprintf(stderr, "Thread %d reporting for duty. \n", mpirank);
    
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    for(iter=0; iter<N; iter++){
        if (mpirank==0) {
            if (iter==0) {
                MPI_Send(message, 500000, MPI_INT, next, tag, MPI_COMM_WORLD);
            } else {
                MPI_Recv(message,  500000, MPI_INT, prev, tag, MPI_COMM_WORLD, &status);
                MPI_Send(message, 500000, MPI_INT, next, tag, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(message,  500000, MPI_INT, prev, tag, MPI_COMM_WORLD, &status);
            MPI_Send(message, 500000, MPI_INT, next, tag, MPI_COMM_WORLD);
        }
    }
    if (mpirank==0) {
        MPI_Recv(message,  500000, MPI_INT, prev, tag, MPI_COMM_WORLD, &status);
    }
    
    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    
    fprintf(stderr, "Thread %d done. \n", mpirank);
    
    printf("Time elapsed is %f seconds.\n", elapsed);
    
    printf("Data bandwidth is: %f GB/s\n", 2*N*mpisize*5E5*sizeof(int)/elapsed/1E9);
    
    MPI_Finalize();
    
    
  	return 0;
}