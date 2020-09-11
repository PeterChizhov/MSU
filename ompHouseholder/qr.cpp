#include <omp.h>
#include <iostream>
#include <mkl.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include "reflector.h"

#define block_size 16

void QR_Decompose(int n, const double * A, double * Q, double * R)
{

/*
	printf("%d NUMTHREADS\n", omp_get_num_threads());
	
	#pragma omp master
	printf("%d THREADNUM\n", omp_get_thread_num());
*/
	

	int ione = 1;
	int izero = 0;
	double dzero = 0.0;

	/* Initialization */

	#pragma omp for
	for(int i = 0; i < n; i++)
	{
		dcopy(&n, &dzero, &izero, Q + i * n, &ione);
		Q[i + n * i] = 1.0;
	};

	#pragma omp for
	for(int i = 0; i < n; i++)
	{
		dcopy(&n, A + i * n, &ione, R + i * n, &ione);
	};

	/* Decomposition */
/*
	//non-block
	
	static double * reflector;
	#pragma omp single
	reflector = new double[n];
	#pragma omp barrier


	//cycle over reflectors
	for(int i = 0; i < n; i++)
	{
	#pragma omp single
		//computing reflector vector
		GetReflector(n - i, R + i * n + i, reflector);
#pragma omp barrier		
		#pragma omp for
		for(int j = i; j < n; j++)
		{
			//multiplication: R
			ApplyReflector(n - i, reflector, R + j * n + i);
		};
		#pragma omp barrier

		#pragma omp for
		for(int j = 0; j < n; j++)
		{
			//multiplication: Q
			ApplyReflector(n - i, reflector, Q + j * n + i);
		};
		#pragma omp barrier

	};


#pragma omp single
	delete[] reflector;
*/	

		
	//block version
	if (n % block_size != 0)
	{
		printf("Error: Invalid block partitioning.\n");
	};
	int blocks = n / block_size;


static	double * reflector;

#pragma omp single
	//reflector = new double[n * block_size];
    reflector = new double[n * n];
#pragma omp barrier

	for(int ib = 0; ib < blocks; ib++)
	{
		// ldr = leading dimension row (otstup do sled element in the row of each block);
		int ldr = n - ib * block_size;

		//computing reflector block
		for(int i = 0; i < block_size; i++)
		{
			// the index in Matrix R corresponding to i index inside the block; 
			int index = ib * block_size + i;

		#pragma omp single	
			//reflectors in the block cant be made parallel

{			//GetReflector(n - index, R + index * n + index, reflector + i * ldr + i);
		    GetReflector(n - index, R + index * n + index, reflector + n * block_size * ib + i * ldr + i);
			//dcopy(&i, &dzero, &izero, reflector + i * ldr, &ione);
			dcopy(&i, &dzero, &izero, reflector + n * block_size * ib + i * ldr, &ione);
			//  applying every i reflector to the whole block, so next reflectors would be fine;

			for(int j = i; j < block_size; j++)
			{
				int jndex = ib * block_size + j;

				ApplyReflector(n - index, reflector + n * block_size * ib + i * ldr + i, R + jndex * n + index);
			}
}			
		};
		#pragma omp barrier

		static int NumberOfThreads = omp_get_num_threads();
		int id = omp_get_thread_num();
		for (int i = 0; ib + 1 + id + i*NumberOfThreads < blocks; ++i) {
			//first shift = (ib + 1)* block_size * n + ib * block_size
			ApplyReflectorArray(ldr,
								block_size,
								reflector + n * block_size * ib,
								ldr,
								R + (ib + 1 + id)* block_size * n + ib * block_size + i * NumberOfThreads * n * block_size,
								n);
		}
		#pragma omp barrier
		/*#pragma omp for
		for(int jb = ib + 1; jb < blocks; jb++)
		{
			ApplyReflectorArray(ldr, block_size, reflector + n * block_size * ib, ldr, R + jb * block_size * n + ib * block_size, n);
		};
		#pragma omp barrier*/

		//#pragma omp for
		//for(int jb = 0; jb < blocks; jb++)
		//{
		//	ApplyReflectorArray(ldr, block_size, reflector, ldr, Q + jb * block_size * n + ib * block_size, n);
		//};
		//#pragma omp barrier
	};

	for (int ib = 0; ib < blocks; ib++) 
	{
		int ldr = n - ib * block_size;
        #pragma omp for
		for (int jb = 0; jb < blocks; jb++)
		{
			ApplyReflectorArray(ldr, block_size, reflector + n * block_size * ib, ldr, Q + jb * block_size * n + ib * block_size, n);
		};
        #pragma omp barrier
	}

#pragma omp single
delete[] reflector;
#pragma omp barrier		

	//transposing...
#pragma omp single
{
	double * tmp = new double[n * n];
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			tmp[i + j * n] = Q[j + i * n];
		};
	};

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			Q[i + j * n] = tmp[i + j * n];
		};
	};
	delete[] tmp;

}
};

void QR_Decompose_MKL(int n, const double * A, double * Q, double * R)
{
	
	int sz2 = n * n;
	int ione = 1;
	int info;

	double * tau = new double[n];
	int lwork = 64 * n;
	double * work = new double[lwork];

	dcopy(&sz2, A, &ione, R, &ione);

	info = 257;
	dgeqrf(&n, &n, R, &n, tau, work, &lwork, &info);
	if (info != 0) printf("dgeqrf failed!\n");

	dcopy(&sz2, R, &ione, Q, &ione);
	info = 257;
	dorgqr(&n, &n, &n, Q, &n, tau, work, &lwork, &info);
	if (info != 0) printf("dorgqr failed!\n");

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < i; j++)
		{
			R[i + j * n] = 0.0;
		};
	};

	delete[] tau;
	delete[] work;
};

void PrintMatrix(int n, double * A)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			printf("%2.5e ", A[i + j * n]);
		};
		printf("\n");
	};
};

int main(int argc, char ** argv)
{
	int n = 4096;
	double * A = new double[n * n];
	double * Q = new double[n * n];
	double * R = new double[n * n];

	double time;

	srand(123);

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			//A[i + j * n] = 1.0 / (1.0 + i + j);
			A[i + j * n] = -1.0 + 2.0 * rand() / (double) RAND_MAX;
		};
	};

	if (n <= 16)
	{
		printf("Matrix A: \n");
		PrintMatrix(n, A);	
	};

	//time = omp_get_wtime();
	clock_t vremya = clock();
	
    #pragma omp parallel 
	QR_Decompose(n, A, Q, R);
	//QR_Decompose_MKL(n, A, Q, R);

	vremya = clock() - vremya;
	//time += omp_get_wtime();

	//printf("Decomposition time: %f seconds.\n", time);
	printf("Decomposition time: %f seconds.\n", ((float)vremya)/CLOCKS_PER_SEC);
	if (n <= 16)
	{
		printf("Matrix Q: \n");
		PrintMatrix(n, Q);
	};

	if (n <= 16)
	{
		printf("Matrix R: \n");
		PrintMatrix(n, R);
	};

	//checks
	double dzero = 0;
	double done = 1;
	double dmone = -1;
	int ione = 1;
	char cN = 'N';
	char cT = 'T';
	int sz2 = n * n;

	double nrm = dnrm2(&sz2, A, &ione);
	dgemm(&cN, &cN, &n, &n, &n, &dmone, Q, &n, R, &n, &done, A, &n);
	double errmult = dnrm2(&sz2, A, &ione) / nrm;

	dgemm(&cN, &cT, &n, &n, &n, &done, Q, &n, Q, &n, &dzero, A, &n);
	for(int i = 0; i < n; i++)
	{
		A[i * n + i] -= 1.0;
	};
	double errQ = dnrm2(&sz2, A, &ione);

	printf("A - QR relative residual: %e\n", errmult);
	printf("QQ^T - I residual: %e\n", errQ);

	delete[] A;
	delete[] Q;
	delete[] R;

	int fuckCoach = 1;
	std::cin >> fuckCoach;
};
