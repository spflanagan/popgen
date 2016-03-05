#pragma once
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// The following code has to do with the random number
// generator. This is the Mersenne Twister prng.
// The reference is:
// M. Matsumoto and T. Nishimura,
// Mersenne Twister: A 623-Dimensionally Equidistributed
// Uniform Pseudo-Random Number Generator,
// ACM Transactions on Modeling and Computer Simulation,
// Vol. 8, No. 1, January 1998, pp. 3-30
// Copyright corresponding to Mersenne Twister Code:
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. The names of its contributors may not be used to endorse or promote
// products derived from this software without specific prior written
// permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)
#define My_PI 3.14159265358979323846

static unsigned long mt[N]; /* the array for the state vector */
static int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */
double Two2the36 = 4294967296.0;
static double rnZ = 0;

/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long int seed)
{
	/* setting initial seeds to mt[N] using the generator Line 25 of Table 1 in	[KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102] */
	mt[0] = seed & 0xffffffff;
	for (mti = 1; mti<N; mti++)
		mt[mti] = (69069 * mt[mti - 1]) & 0xffffffff;
}

inline double genrand()
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0, MATRIX_A };
	/* mag01[x] = x * MATRIX_A for x=0,1 */
	if (mti >= N) { /* generate N words at one time */
		int kk;
		if (mti == N + 1) /* if sgenrand() has not been called, */
			sgenrand(4357); /* a default initial seed is used */
		for (kk = 0; kk<N - M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
		}

		for (; kk<N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

	return ((double)y / Two2the36); /* reals */
	//return y; /* for integer generation */
	// This should give a range of [0,1); for [0,1]
	// use Two2the36-1.
}

inline int randnum(int Max) // returns a value between 0 and Max-1.
{
	double rnms;
	double dbrnum;
	int irnum;
	rnms = genrand();
	dbrnum = floor(rnms * Max);
	irnum = dbrnum;
	return irnum;
}

// The following function is a standard routine for generating a random number from
// a normal distribution.
inline double randnorm(double rnmean, double rnstd)
{
	double S, v1, v2, X2;
	if (rnZ != 0)
	{
		X2 = rnZ;
		rnZ = 0;
	}
	else
	{
		do {
			v1 = 2.0 * genrand() - 1.0;
			v2 = 2.0 * genrand() - 1.0;
			S = v1 * v1 + v2 * v2;
		} while (S >= 1.0);
		S = sqrt((-2.0 * log(S)) / S);
		X2 = v1 * S;
		rnZ = v2 * S;
	}
	return rnmean + X2 * rnstd;
}

// The function below will generate random numbers from a bivariate
// normal distribution. The random numbers will be in the global
// variables bivN1 and bivN2. The std deviations of the two
// distributions
// are given by sigmaX and sigmaY, and the correlation is rho.
// These bivariate normal routines are from the GNU free software database
// (C) 2000 James Theiler and Brian Gough.

inline void randbivnorm(double sigmaX, double sigmaY, double rho, double &bivN1, double &bivN2)
{
	double bivU, bivV, bivR2, bivScale;

	do
	{
		// choose x and y from a uniform square (-1, -1) to (1, 1)
		bivU = 2 * genrand() - 1;
		bivV = 2 * genrand() - 1;
		// see if it is in the unit circle
		bivR2 = bivU * bivU + bivV * bivV;
	} while (bivR2 > 1.0 || bivR2 == 0);
	bivScale = sqrt(-2.0 * log(bivR2) / bivR2);
	bivN1 = sigmaX * bivU * bivScale;
	bivN2 = sigmaY * (rho * bivU + sqrt(1 - rho*rho) * bivV) * bivScale;
} // end of randbivnorm

inline double bivnormpdf(const double bnX, const double bnY, const double sigmaX, const double sigmaY, const double rho)
{
	double bnU = bnX / sigmaX;
	double bnV = bnY / sigmaY;
	double bnC = 1 - rho*rho;
	double bnP = (1 / (2 * My_PI * sigmaX * sigmaY * sqrt(bnC))) * exp(-(bnU * bnU - 2 * rho * bnU * bnV + bnV * bnV) / (2 * bnC));
	return bnP;
} // end of bivnormpdf

//took these from Adam's simulation
inline int poissonrand(double Pmean) // Returns a Poisson-distributed random number
{
	double PL = exp(-1 * Pmean);
	int Pkk = 0;
	double Ppp = 1;
	do
	{
		Pkk++;
		Ppp = Ppp*genrand();
	} while (Ppp > PL);
	return Pkk - 1;
}  // end of poissonrand

inline int positiveroundnorm(double PRmean, double PRstddev)  // Normal distribution rounded to nearest integer, with negatives converted to zeros
{
	double prnnumber, prnfloor;
	int prnint;
	prnnumber = randnorm(PRmean, PRstddev);
	if (prnnumber <= 0)
		return 0;
	prnfloor = floor(prnnumber);
	if (prnnumber - prnfloor < 0.5)
		prnint = static_cast<int>(prnfloor);
	else
		prnint = static_cast<int>(prnfloor + 1);
	return prnint;
}


double lowerincompletegamma(double shape, double x)
{
	double sum = 0;
	double term = 1.0 / shape;
	double n = 1;

	while (term != 0 && n < 100) {
		sum = sum + term;
		term = term*(x / (shape + n));
		n++;
	}

	return pow(x, shape)*exp(-1 * x)*sum;
}

double gammln(double xx) // From numerical recipes in C
{
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146, -86.50532032941677, 24.01409824083091,
		-1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
	int j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005*ser / x);
}

double chisqr(double teststat, double df)
{
	double numerator;
	double denominator;

	numerator = lowerincompletegamma(df / 2, teststat / 2);
	denominator = exp(gammln(df / 2));

	return numerator / denominator;
}


//this function was written by Sarah P. Flanagan
//17 December 2015
//It reads in a file containing known allele frequency distributions
//and returns an integer for which allele index value to use
using namespace std;
int d_empirical_afs(int num_alleles, string afs_file_name)
{
	//read in the known distribution
	ifstream afs_file;
	string line;
	vector<double> intervals;
	vector<double> densities;
	vector<int> count_intervals;
	double t, tt, last_interval, last_prob, sum;
	int j, jj, jjj, int_index, count, return_allele;
	bool found;

	afs_file.open(afs_file_name);
	while (getline(afs_file, line))
	{
		stringstream temp;
		temp.str(line);
		temp >> t >> tt;
		intervals.push_back(t);
		densities.push_back(tt);
	}
	//need to divide up the known distribution for the correct number of alleles
	vector<double> allele_probs;
	int num_allele_intervals = intervals.size() / num_alleles;
	if (num_allele_intervals*num_alleles < intervals.size())
		num_allele_intervals++;
	for (j = 0; j < num_alleles; j++)
		allele_probs.push_back(0);
	int_index = count = 0;
	count_intervals.push_back(0);
	for (j = 0; j < intervals.size(); j++)
	{
		if (count >= num_allele_intervals)
		{
			count_intervals.push_back(0);
			int_index++;
			count = 0;
		}
		allele_probs[int_index] = allele_probs[int_index] + densities[j];
		count_intervals[int_index]++;
		count++;
	}
	sum = 0;
	for (j = 0; j < num_alleles; j++)
	{
		sum = sum + allele_probs[j];
		allele_probs[j] = allele_probs[j] / count_intervals[j];
	}
	for (j = 0; j < num_alleles; j++)
		allele_probs[j] = allele_probs[j] / sum;

	//now need to draw a number from that distribution and know which allele it should be.
	found = false;
	while (!found)
	{
		count = randnum(num_alleles);
		t = genrand();
		if (count*t < allele_probs[count])
		{
			return_allele = count;
			found = true;
		}
	}
	return return_allele;
}
