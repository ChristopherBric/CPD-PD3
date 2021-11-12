#include <mpi.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;

int main (int argc, char **argv)
{
    int rank, size;
    int N, i, j, K;
    double t0=0.0,t1=0.0,t2=0.0,tf=0.0;
    int  **A, *v, *x, *F;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2)
    {
        printf("Se necesita como minimo 2 procesadores\n");
        MPI_Abort(MPI_COMM_WORLD,-1);
        return 0;
    }

	if (rank==0)
	{
		cin >> N;
		if (rank==0) printf("Matriz %d x %d\n",N,N);
	}
	MPI_Bcast(&N,1,MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
    A = new int *[(N*N)/size];
    v = new int [N];

    if (rank==0)
    {
        A[0] = new int [N*N];
        for(int i=1;i<N;i++)
        {
            A[i]=A[i-1]+ N;
        }
        
        x = new int [N];

        srand(time(0));
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++)
            {
                A[i][j]= (rand()%100)+1;
            }
            v[i] = rand()%100;
        }
    }

    F = new int [(N*N)/size];

    MPI_Barrier(MPI_COMM_WORLD);

    // Comunicacion
    t0 = MPI_Wtime();

    MPI_Bcast(&v,N,MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(A[0],(N*N)/size,MPI_INT,F,(N*N)/size,MPI_INT,0,MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

    t1 = MPI_Wtime();

    int SF[(N*N)/size]={};

    for(int i=0;i<(N*N)/size;i++)
    {
        SF[i]=0;
        for(int j=0;j<N;j++)
        {
            SF[i] += F[(N*i)+j]*v[j];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    t2 = MPI_Wtime();

    MPI_Gather(&SF,(N*N)/size,MPI_INT,A,(N*N)/size,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    tf = MPI_Wtime();

    if(rank==0)
    {
        printf("Tiempo Com. Inicial: %1.6f ms\n",1000*(t1-t0));
        printf("Tiempo Comp.: %1.6f ms\n",1000*(t2-t1));
        printf("Tiempo Com. Final: %1.6f ms\n",1000*(tf-t2));
    }

    MPI_Finalize();
}