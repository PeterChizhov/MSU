#include <omp.h>

#include <mkl.h>
#include <cmath>
#include "reflector.h"

void GetReflector(int sz, const double * vector, double * reflector)
{
	double nrm2;
	int ione = 1;

	nrm2 = ddot(&sz, vector, &ione, vector, &ione);
	dcopy(&sz, vector, &ione, reflector, &ione);

	if (sz != 1)
	{
		if (vector[0] >= 0)
		{
			nrm2 -= reflector[0] * reflector[0];
			reflector[0] -= sqrt(nrm2 + reflector[0] * reflector[0]);
			nrm2 += reflector[0] * reflector[0];
		}
		else
		{
			nrm2 -= reflector[0] * reflector[0];
			reflector[0] += sqrt(nrm2 + reflector[0] * reflector[0]);
			nrm2 += reflector[0] * reflector[0];
		};
	};
	nrm2 = 1.0 / sqrt(nrm2);

	dscal(&sz, &nrm2, reflector, &ione);
};

void ApplyReflector(int sz, const double * reflector, double * vec)
{
	double scalar;
	int ione = 1;

	scalar = - 2.0 * ddot(&sz, vec, &ione, reflector, &ione);
	daxpy(&sz, &scalar, reflector, &ione, vec, &ione);
};

void ApplyReflectorArray(int sz, int block_sz, const double * reflector, int ldr, double * block, int ldblock)
{
	double * scalars;
	double dm2 = -2.0;
	double done = 1.0;
	double dzero = 0.0;
	int ione = 1;
	char cN = 'N';
	char cT = 'T';

	scalars = new double[block_sz];

	for(int i = 0; i < block_sz; i++)
	{
		dgemv(&cT, &sz, &block_sz, &done, block, &ldblock, reflector + i * ldr, &ione, &dzero, scalars, &ione);
		dger(&sz, &block_sz, &dm2, reflector + i * ldr, &ione, scalars, &ione, block, &ldblock);
	};

	delete[] scalars;
};
