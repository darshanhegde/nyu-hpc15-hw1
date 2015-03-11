EXECS=int-ring jacobi-mpi data-ring
MPICC?=mpicc

all: ${EXECS}

int-ring: int-ring.c
	${MPICC} -o int-ring int-ring.c

jacobi-mpi: jacobi-mpi.c
	${MPICC} -o jacobi-mpi jacobi-mpi.c

data-ring: data-ring.c
	${MPICC} -o data-ring data-ring.c

clean:
	rm ${EXECS}
