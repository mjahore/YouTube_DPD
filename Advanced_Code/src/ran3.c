/*
 * ran3.c - This is the ran3() function as presented in Numerical Recipes in C. It is sufficient
 *          for DPD simulation in terms of an extremely long period for the RNG. idum is the seed
 *          and should be set initially to a large, negative integer. ran3() returns a random number
 *          between 0 and 1 from a uniform distribution. 
 */

// The following 4 lines will be global variables.
int inext,inextp;
long ma[56];    
int iff=0;     
int iseed;    

#define __RAN3__
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(int *idum)
{
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=labs(MSEED-labs(*idum));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

