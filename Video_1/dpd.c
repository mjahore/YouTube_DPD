/***********************************************************************************
 *
 * dpd.c - Example DPD simulation of a polymer in solution
 *         from YouTube video series available at:
 *         https://www.youtube.com/playlist?list=PLx9F_EFL1lVaAqnF-k091LjhjGvC0pLOq
 *
 * Author: Michael J. A. Hore (hore@case.edu)
 *
 * Repository: https://www.github.com/mjahore/YouTube_DPD
 *
 * License: This code is licensed under GPL v3. Please see
 *          LICENSE for more details.
 ***********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795f

// Random Number Generator
int inext,inextp;      // Ran3
long ma[56];           // Ran3
int iff=0;             // Ran3
int iseed;             // Random seed
#include "ran3.c"      // ran3() from Numerical Recipes in C

// Global Variables
double *x, *y, *z;        // Coordinates w/ PBC
double *rx, *ry, *rz;     // Coordinates w/o PBC
double *vx, *vy, *vz;     // Velocities (current time step)
double *vxp, *vyp, *vzp;  // Velocities (prev. time step)
double *fx, *fy, *fz;     // Forces (current time step)
double *fxp, *fyp, *fzp;  // Forces (prev. time step)
int    *p_type;           // Particle type or species
double Lx, Ly, Lz;        // System dimensions
double rho;               // Fluid density
double poly_phi;          // Polymer volume fraction
int    N;                 // Polymer deg. polym.
double b;                 // Bond length
double b_avg;             // Average bond length from sim.
double T;                 // Specified value of kBT
double U;                 // Potential energy
double K;                 // Kinetic energy
double a_ii;              // Conservative int. strength (like particles)
double a_ij;              // Conservative int. strength (unlike particles)
double sigma, gamma1;     // Random and dissipative strengths.
double kH;                // Bond spring const.
double dt;                // Integration time step.
int    time_equil;        // # of equilbration time steps.
int    time_prod;         // # of production time steps.
int    t;                 // Current time step
int    n_dpd;             // Total number of DPD particles in sim.
int    n_poly;            // Number of polymers in the simulation.

// Function prototypes.
void initialize_system(void);
void read_parameters(char *filename);
void write_xyz(char *filename, int append);

int main(int argc, char *argv[]) {

	// Cosmetic stuff.
	printf("\n");

	// Set seed value for ran3()
	iseed = -time(0);

	// Read parameters from disk.
	read_parameters("param.cfg");

	// Create a random configuration of particles.
	initialize_system();

	// Write out the polymer coordinates to initial_config.xyz.
	write_xyz("initial_config.xyz", 0);

	// Calculate forces:

	// Move particles/update velocities with velocity-Verlet algorithm.

	// Calculate some stuff (Rg, T, p, etc.)

	// End!

	// More cosmetic stuff.
	printf("\n");

	free(x);
	free(y);
	free(z);
	free(rx);
	free(ry);
	free(rz);
	free(vx);
	free(vy);
	free(vz);
	free(vxp);
	free(vyp);
	free(vzp);
	free(fx);
	free(fy);
	free(fz);
	free(fxp);
	free(fyp);
	free(fzp);
	free(p_type);
	return 0;
}

void initialize_system(void) {
	int i,j,k;
	double phi, theta;
	double vxt, vyt, vzt;

	// Figure out how many DPD particles are in the system:
	n_dpd = (int) (rho * Lx * Ly * Lz);

	// Calculate the number of polymers in the system.
	n_poly = (int) (poly_phi * n_dpd / N);

	// Write some output to screen
	printf("There are %d DPD particles and %d polymers with chain length %d.\n", n_dpd, n_poly, N);

	// Allocate memory for our x, y,z, etc. arrays.
	x       = (double*) malloc(n_dpd * sizeof(double));
	y       = (double*) malloc(n_dpd * sizeof(double));
	z       = (double*) malloc(n_dpd * sizeof(double));
	rx      = (double*) malloc(n_dpd * sizeof(double));
	ry      = (double*) malloc(n_dpd * sizeof(double));
	rz      = (double*) malloc(n_dpd * sizeof(double));
	vx      = (double*) malloc(n_dpd * sizeof(double));
	vy      = (double*) malloc(n_dpd * sizeof(double));
	vz      = (double*) malloc(n_dpd * sizeof(double));
	vxp     = (double*) malloc(n_dpd * sizeof(double));
	vyp     = (double*) malloc(n_dpd * sizeof(double));
	vzp     = (double*) malloc(n_dpd * sizeof(double));
	fx      = (double*) malloc(n_dpd * sizeof(double));
	fy      = (double*) malloc(n_dpd * sizeof(double));
	fz      = (double*) malloc(n_dpd * sizeof(double));
	fxp     = (double*) malloc(n_dpd * sizeof(double));
	fyp     = (double*) malloc(n_dpd * sizeof(double));
	fzp     = (double*) malloc(n_dpd * sizeof(double));
	p_type  = (int*)    malloc(n_dpd * sizeof(int));

	// Put the polymers into the system.
	vxt = 0.0;
	vyt = 0.0;
	vzt = 0.0;
	k = 0;

	for(i=0; i<n_poly; i++) {
		// Assign a random position in the system.
		x[k] = ran3(&iseed) * Lx;
		y[k] = ran3(&iseed) * Ly;
		z[k] = ran3(&iseed) * Lz;

		// Assign a velocity to this particle.
		theta = ran3(&iseed) * PI;
		phi   = ran3(&iseed) * 2.0 * PI;
		vx[k] = sqrt(T) * cos(phi) * sin(theta);
		vy[k] = sqrt(T) * sin(phi) * sin(theta);
		vz[k] = sqrt(T) * cos(theta);
		vxt  += vx[k];
		vyt  += vy[k];
		vzt  += vz[k];

		// Assign the particle type.
		p_type[k] = 1;

		k++;
	
		for (j=1; j<N; j++) {
			// Assign a random position in the system.
			theta = ran3(&iseed) * PI;
			phi   = ran3(&iseed) * 2.0 * PI;
			x[k]  = x[k-1] + sqrt(b) * cos(phi) * sin(theta);
			y[k]  = y[k-1] + sqrt(b) * sin(phi) * sin(theta);
			z[k]  = z[k-1] + sqrt(b) * cos(theta);

			// Assign a velocity to this particle.
			theta = ran3(&iseed) * PI;
			phi   = ran3(&iseed) * 2.0 * PI;
			vx[k] = sqrt(T) * cos(phi) * sin(theta);
			vy[k] = sqrt(T) * sin(phi) * sin(theta);
			vz[k] = sqrt(T) * cos(theta);
			vxt  += vx[k];
			vyt  += vy[k];
			vzt  += vz[k];

			// Assign the particle type.
			p_type[k] = 1;

			k++;
		}
	}

	// Generate or place the solvent particles into the system.
	for(i=k; i<n_dpd; i++) {
		// Random position!
		x[k] = ran3(&iseed) * Lx;
		y[k] = ran3(&iseed) * Ly;
		z[k] = ran3(&iseed) * Lz;

		// Assign a velocity to this particle.
		theta = ran3(&iseed) * PI;
		phi   = ran3(&iseed) * 2.0 * PI;
		vx[k] = sqrt(T) * cos(phi) * sin(theta);
		vy[k] = sqrt(T) * sin(phi) * sin(theta);
		vz[k] = sqrt(T) * cos(theta);
		vxt  += vx[k];
		vyt  += vy[k];
		vzt  += vz[k];

		// Set particle type.
		p_type[k] = -1;

		k++;
	}

	// Set the non-PBC arrays.
	memcpy(rx, x, n_dpd*sizeof(double));
	memcpy(ry, y, n_dpd*sizeof(double));
	memcpy(rz, z, n_dpd*sizeof(double));

	// Subtract off the average velocity to eliminate any drift.
	vxt /= (double) n_dpd;
	vyt /= (double) n_dpd;
	vzt /= (double) n_dpd;
	for (i=0; i<n_dpd; i++) {
		vx[i] -= vxt;
		vy[i] -= vyt;
		vz[i] -= vzt;
	}

	// BOOM WE ARE DONE!
	return;
}

void read_parameters(char *filename) {
	FILE *param_file;
	char junk[256];

	param_file = fopen(filename, "r");
	if (param_file == NULL) {
		printf("Parameter file %s not found.\n", filename);
		exit(-1);
	}

	// Read in the comment line.
	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &Lx);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &Ly);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &Lz);
	
	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &rho);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &poly_phi);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%d\n", &N);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &b);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &T);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &a_ii);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &a_ij);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &kH);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &sigma);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%lf\n", &dt);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%d\n", &time_equil);

	fgets(junk, 256, param_file);
	fscanf(param_file, "%d\n", &time_prod);

	// Define gamma1 from sigma.
	gamma1 = pow(sigma,2.0) / (2.0 * T);
	sigma  = sigma * sqrt(3.0/dt);
	return;
}

void write_xyz (char* filename, int append) {
	FILE *output;
	int i;

	printf("\nWriting coordinates to %s\n\n", filename);

	if (append) {
		output = fopen(filename, "aw");
	} else {
		output = fopen(filename, "w");
	}

	fprintf(output, "%d\n", n_poly*N);
	fprintf(output, "# Snapshot from time t = %d\n", t);
	for(i=0; i<n_poly*N; i++) {
		fprintf(output, "C %3.6lf %3.6lf %3.6lf\n", x[i], y[i], z[i]);
	}
	fclose(output);

	return;
}
