#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "random.h"

#define boolean         int
#define RIGHT           -1
#define WRONG           0

double u[98], c, cd, cm;
int i97, j97;
boolean test = WRONG;

void rmarin(int ij, int kl)
{
	/*
	 C This is the initialization routine for the random number generator RANMAR()
	 C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
	 C                                                      0 <= KL <= 30081
	 C The random number sequences created by these two seeds are of sufficient
	 C length to complete an entire calculation with. For example, if sveral
	 C different groups are working on different parts of the same calculation,
	 C each group could be assigned its own IJ seed. This would leave each group
	 C with 30000 choices for the second seed. That is to say, this random
	 C number generator can create 900 million different subsequences -- with
	 C each subsequence having a length of approximately 10^30.
	 C
	 C Use IJ = 1802 & KL = 9373 to test the random number generator. The
	 C subroutine RANMAR should be used to generate 20000 random numbers.
	 C Then display the next six random numbers generated multiplied by 4096*4096
	 C If the random number generator is working properly, the random numbers
	 C should be:
	 C           6533892.0  14220222.0  7275067.0
	 C           6172232.0  8354498.0   10633180.0
	 */
	int i, j, k, l, ii, jj, m;
	double s, t;

	if (ij<0 || ij>31328 || kl<0 || kl>30081)
	{
		printf("The first random number seed must have a value between 0 and 31328.");
		printf("The second seed must have a value between 0 and 30081.\n");
		exit(1);
	}

	i = (ij/177)%177 + 2;
	j = ij%177 + 2;
	k = (kl/169)%178 + 1;
	l = kl%169;

	for (ii=1; ii<=97; ii++)
	{
		s = 0.0;
		t = 0.5;
		for (jj=1; jj<=24; jj++)
		{
			m = (((i*j)%179)*k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l + 1) % 169;
			if ((l*m)%64 >= 32)
				s += t;
			t *= 0.5;
		}
		u[ii] = s;
	}

	c = 362436.0 / 16777216.0;
	cd = 7654321.0 / 16777216.0;
	cm = 16777213.0 / 16777216.0;

	i97 = 97;
	j97 = 33;

	test = RIGHT;
}

void ranmar(double rvec[], int len)
/* len random numbers are placed in rvec[0..len-1] */
/*
 C This is the random number generator proposed by George Marsaglia in
 C Florida State University Report: FSU-SCRI-87-50
 C It was slightly modified by F. James to produce an array of pseudorandom
 C numbers.
 */
{
	int ivec;
	double uni;

	if (test==WRONG)
	{
		printf("ERROR: Call the init routine rmarin() before calling ranmar().\n");
		exit(1);
	}
	for (ivec=0; ivec<len; ivec++)
	{
		uni = u[i97] - u[j97];
		if (uni < 0.0)
			uni += 1.0;
		u[i97] = uni;
		i97--;
		if (i97==0)
			i97 = 97;
		j97--;
		if (j97==0)
			j97 = 97;
		c -= cd;
		if (c<0.0)
			c += cm;
		uni -= c;
		if (uni<0.0)
			uni += 1.0;
		rvec[ivec] = uni;
	}
}

/* Gaussian random number generator */
double gaussran()
{
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;
	double rand_num[2];

	if (iset == 0)
	{
		do
		{
			ranmar(rand_num, 2);
			v1 = 2.0 * rand_num[0] - 1.0;
			v2 = 2.0 * rand_num[1] - 1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1 * fac;
		iset = 1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}
//---------------------------------------------------------


