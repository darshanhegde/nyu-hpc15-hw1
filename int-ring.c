#include <mpi.h>
#include <stdio.h>
#include <assert.h>

main(int argc, char *argv[])  {
    int mpisize, mpirank, N;
    long message=0;
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
    int next = (mpirank==mpisize-1) ? 0: mpirank+1;
    
    int iter;
    
    for(iter=0; iter<N; iter++){
        if (mpirank==0) {
            if (iter==0) {
                MPI_Send(&message, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&message,  1, MPI_INT, prev, tag, MPI_COMM_WORLD, &status);
                message += mpirank;
                MPI_Send(&message, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(&message,  1, MPI_INT, prev, tag, MPI_COMM_WORLD, &status);
            message += mpirank;
            MPI_Send(&message, 1, MPI_INT, next, tag, MPI_COMM_WORLD);
        }
    }
    
    if (mpirank==0) {
        MPI_Recv(&message,  1, MPI_INT, prev, tag, MPI_COMM_WORLD, &status);
    }
    
    // check for correctness of implementation
    long ans = ((mpisize*(mpisize-1))/2)*(N-1) + (mpirank*(mpirank+1))/2;
    assert(message==ans);
    
    MPI_Finalize();
  	return 0;
}