#include <mpi.h>
#include <mkl.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>

#define world MPI_COMM_WORLD

#define matrix_size 8100

double MatrixFormula(int i, int j)
{
	return sin((double) (3.0 * i + 7.0 * j + 1.0));
	//return -1.0/((1.0 + i + j)*(1.0 + i + j));
};

int main(int argc, char ** argv)
{
	int mpi_world_size_sq;
	int mpi_size;
	int myid;
	int myid_I;
	int myid_J;

	int mtx_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(world, &mpi_world_size_sq);
	MPI_Comm_rank(world, &myid);

	mpi_size = sqrt(mpi_world_size_sq);
	mtx_size = matrix_size / mpi_size;

    double* mtx = new double[mtx_size * mtx_size];

	if (mpi_size * mpi_size != mpi_world_size_sq || mtx_size * mpi_size != matrix_size)
	{
		if (myid == 0)
		{
			printf("The number of processes is not a full square or matrix size is not divisible by process count, exiting...");
		}; 

		MPI_Finalize();
		return -2;
	};

	myid_I = myid % mpi_size;
	myid_J = myid / mpi_size;

	char cN = 'N';
	double done = 1.0;
	double dzero = 0.0;
	double dmone = -1.0;
	int ione = 1;
	int izero = 0;

	double first_max;
    double err = 1.0;

	int maxi, maxj;

	double *recvbuf1 = new double[mpi_world_size_sq];
	int *recvbuf2 = new int[2*mpi_world_size_sq];
	int *sendbuf = new int[2];

    double* sendrecvRow = new double[mtx_size];
    double* sendrecvCol = new double[mtx_size];
    

	int istart;
	int jstart;

	istart = myid_I * mtx_size;
	jstart = myid_J * mtx_size;

	for(int i = 0; i < mtx_size; i++)
	{
		for(int j = 0; j < mtx_size; j++)
		{
			mtx[i + j * mtx_size] = MatrixFormula(istart + i, jstart + j);
		};
	};

    int k;
    k = 0;

    while(err > 1.0e-5)
    {

		//searching for the maximum element of matrix
		double max = 0;
        for(int i = 0; i < mtx_size; i++)
        {
			for(int j = 0; j < mtx_size; j++)
	        {
				if(fabs(mtx[i + j * mtx_size]) >= fabs(max))
				{
				max = mtx[i + j * mtx_size];
				sendbuf[0] = i + istart;
				sendbuf[1] = j + jstart;
				}
			}
		}

        MPI_Allgather(sendbuf, 2, MPI_INT, recvbuf2, 2, MPI_INT, world);

        MPI_Allgather(&max, 1, MPI_DOUBLE, recvbuf1, 1, MPI_DOUBLE, world);

		max = 0;

        for(int i = 0; i < mpi_world_size_sq; i++)
        {
			if(fabs(recvbuf1[i]) >= fabs(max))
			{
				max = recvbuf1[i];
				maxi = recvbuf2[2*i];
				maxj = recvbuf2[2*i+1];
            }
        }

		// saving first maximum
        k++;
		if(k == 1)
			 first_max = max;


		//message part 
        int maxProc_i = maxi/mtx_size;
		int maxProc_j = maxj/mtx_size;

		if (myid_I == maxProc_i) 
		{
			for (int i = 0; i < mtx_size; ++i) 
			{
				sendrecvRow[i] = mtx[maxi % mtx_size + i * mtx_size];
			}
			MPI_Bcast(sendrecvRow, mtx_size, MPI_DOUBLE, myid, world);
		}
		MPI_Barrier(world);


		if (myid_J == maxProc_j)
		{
			for (int i = 0; i < mtx_size; ++i)
			{
				sendrecvCol[i] = mtx[i + (maxj % mtx_size) * mtx_size];
			}
			MPI_Bcast(sendrecvCol, mtx_size, MPI_DOUBLE, myid, world);
		}
		MPI_Barrier(world);

		// changing the matrix
        for(int i = 0; i < mtx_size; i++)
		{
			for(int j = 0; j < mtx_size; j++)
			{
				mtx[i + j * mtx_size] -= sendrecvRow[j] * sendrecvCol[i] / max;
			}
        }
     
        err = fabs(max) / fabs(first_max);

		if (myid == 0)
		{
			printf("\n");
			printf("MAX: %lf; indexes: %d, %d; err: %lf.\n", fabs(max), maxi, maxj, err);
			printf("\n");
		}
   }

        delete[] sendbuf;
        delete[] recvbuf1;
        delete[] recvbuf2;
        delete[] sendrecvRow;
        delete[] sendrecvCol;
        delete[] mtx;
	
	MPI_Finalize();
};
