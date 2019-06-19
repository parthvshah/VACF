#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	MPI_Init(NULL, NULL);
	int wSize;
	MPI_Comm_size(MPI_COMM_WORLD, &wSize);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char procName[MPI_MAX_PROCESSOR_NAME];
	int nameLen;
	MPI_Get_processor_name(procName, &nameLen);

	printf("HW from %s, rank %d of %d.\n", procName, rank, wSize);

	return 0;
}
