#include <mpi.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>  
#include <iostream>

using namespace std;

int main (int argc, char **argv)
{
	int rank, size;
	char Val; 
	int i;
	double t0=0.0,tf=0.0;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (rank == 0)
		Val = 'F';
	
	t0 = MPI_Wtime();
	
	for (i=0;i<size;i++)
	{
		if (rank == size-1)
		{
			MPI_Send(&Val,1,MPI_CHAR,0,rank,MPI_COMM_WORLD);
			MPI_Recv(&Val,1,MPI_CHAR,rank-1,rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else if (rank == 0)
		{
			MPI_Send(&Val,1,MPI_CHAR,rank+1,rank,MPI_COMM_WORLD);
			MPI_Recv(&Val,1,MPI_CHAR,size-1,size-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Send(&Val,1,MPI_CHAR,rank+1,rank,MPI_COMM_WORLD);
			MPI_Recv(&Val,1,MPI_CHAR,rank-1,rank-1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
	}
	tf = MPI_Wtime();
	
	printf("Latencia: %1.6f ns\n", (tf-t0));
	
	MPI_Finalize();
}