/***********************************************************
 *
 * dpd.c - Example DPD simulation of a polymer in solution
 *         from YouTube video series available at:
 *         https://www.youtube.com/
 *
 * Author: Michael J. A. Hore (hore@case.edu)
 *
 * License: This code is licensed under GPL v3. Please see
 *          LICENSE for more details.
 ***********************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// Global Variables 
double *x, *y, *z;       // Coordinates w/ PBC
double *rx, *ry, *rz;    // Coordinates w/o PBC
double *vx, *vy, *vz;    // Velocities (current time step)
double *vxp, *vyp, *vzp; // Velocities (previous time step)
double *fx, *fy, *fz;    // Forces (current time step)
double *fxp, *fyp, *fzp; // Forces (previous time step)
int    *p_type;          // Particle type
double Lx, Ly, Lz;       // System dimensions
double rho;              // Fluid density
double poly_phi;         // Polymer volume fraction
int    N;                // Polymer deg. of polymerization
double b;                // Polymer bond length
double T;                // Temperature (actually, kB*T)
double a_ii;             // Same-type interaction
double a_ij;             // Different-type interaction
double sigma, gamma;     // Random and dissipative force strengths
double kH;               // Harmonic potential strength.
double dt;               // Integration time step.
int    time_equil;       // Number of equilibration time steps.
int    time_prod;        // Number of production time steps.
int    t;                // Current time step.
int    n_dpd;            // Total number of particles.
int    n_poly;           // Total number of polymers.

// Random Number Generator
int inext,inextp;      // Ran3
long ma[56];           // Ran3
int iff=0;             // Ran3
int iseed;             // Random seed
#include "ran3.c"      // ran3() from Numerical Recipes in C

#define PI 3.1415926535897932384626433832795f

// Function prototypes:
void   calc_forces(void);
void   calc_Rg(void);
double calc_p(void);
double calc_T(void);
void   initialize_system(void);
void   make_list(void);
void   read_parameters(char *filename);
void   velocity_verlet(void);
void   write_xyz(char *filename);

int main(int argc, char *argv[]) {
	// Cosmetic stuff.
	printf("\n");

	// Set ran3() random seed based on current time.
	iseed=-time(0);

	// Read simulation parameters.
	read_parameters("param.cfg");

	// Just some info to make sure things are okay...
	printf("Simulating a %lf x %lf x %lf system with density p = %lf.\n", Lx, Ly, Lz, rho);

	// Set up the initial configuration of the system
	initialize_system();

	// Cosmetic stuff.
	printf("\n");

	// Free the memory that was allocated in initialize_system().
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
	return 0;
}

void initialize_system(void) {
	int i,j,k;
	double phi, theta;     // Spherical coordinates.
	double vxt, vyt, vzt;  // Average velocity in x,y,z directions.

	// Calculate the number of DPD particles in the simulation
	// box.
	n_dpd = (int) (rho * Lx * Ly * Lz);

	// Calculate the number of polymers in the simulation box.
	n_poly = (int) (poly_phi * n_dpd / N);

	printf("There are %d DPD particles and %d polymers with chain length %d.\n", n_dpd, n_poly, N);

	// Now, allocate memory in the arrays that store particle information
	// (x, y, z, fx, fy, ...). This memory is freed at the end of the
	// main() portion of the code.
	x   = (double*) malloc(n_dpd * sizeof(double));
	y   = (double*) malloc(n_dpd * sizeof(double));
	z   = (double*) malloc(n_dpd * sizeof(double));
	rx  = (double*) malloc(n_dpd * sizeof(double));
	ry  = (double*) malloc(n_dpd * sizeof(double));
	rz  = (double*) malloc(n_dpd * sizeof(double));
	vx  = (double*) malloc(n_dpd * sizeof(double));
	vy  = (double*) malloc(n_dpd * sizeof(double));
	vz  = (double*) malloc(n_dpd * sizeof(double));
	vxp = (double*) malloc(n_dpd * sizeof(double));
	vxp = (double*) malloc(n_dpd * sizeof(double));
	vyp = (double*) malloc(n_dpd * sizeof(double));
	fx  = (double*) malloc(n_dpd * sizeof(double));
	fy  = (double*) malloc(n_dpd * sizeof(double));
	fz  = (double*) malloc(n_dpd * sizeof(double));
	fxp = (double*) malloc(n_dpd * sizeof(double));
	fyp = (double*) malloc(n_dpd * sizeof(double));
	fzp = (double*) malloc(n_dpd * sizeof(double));

	// Place the n_poly polymers at random positions in the
	// system. Give them a random velocity based on the value
	// of kBT that is specified.
	vxt = 0.0;
	vyt = 0.0;
	vzt = 0.0;
	k   = 0;
	for (i=0; i<n_poly; i++) {
		// Random position:
		x[k]  = ran3(&iseed) * Lx;
		y[k]  = ran3(&iseed) * Ly;
		z[k]  = ran3(&iseed) * Lz;

		// vx^2 + vy^2a + vz^2 = kBT, so... pick a random direction and
		// give the particle a velocity with a magnitude of T. Theta and
		// phi are two random angles in spherical coordinates (i.e., the 		// polar and azimuthal angles).
		theta = ran3(&iseed) * PI;
		phi   = ran3(&iseed) * 2.0 * PI;
		vx[k] = sqrt(T) * cos(phi) * sin(theta);
		vy[k] = sqrt(T) * sin(phi) * sin(theta);
		vz[k] = sqrt(T) * cos(theta);
		vxt  += vx[k];
		vyt  += vy[k];
		vzt  += vz[k];

		// Pick a particle type. Let's call all monomers
		// type 1.
		p_type[k] = 1;

		k++;

		// Do the remaining monomers.
		for (j=1; j<N; j++) {
			// Find a random direction:
			theta = ran3(&iseed) * PI;
			phi   = ran3(&iseed) * 2.0 * PI;
			
			// The position of this monomer should be the position
			// of the previous monomer + the bond length, oriented
			// in a random direction.
			x[k]  = x[k-1] + sqrt(b) * cos(phi) * sin(theta);
			y[k]  = y[k-1] + sqrt(b) * sin(phi) * sin(theta);
			z[k]  = z[k-1] + sqrt(b) * cos(theta);

			// Same principles for the velocity here:
			theta = ran3(&iseed) * PI;
			phi   = ran3(&iseed) * 2.0 * PI;
			vx[k] = sqrt(T) * cos(phi) * sin(theta);
			vy[k] = sqrt(T) * sin(phi) * sin(theta);
			vz[k] = sqrt(T) * cos(theta);
			vxt  += vx[k];
			vyt  += vy[k];
			vzt  += vz[k];

			// Particle type.
			p_type[k] = 1;
			k++;
		}
	}

	// Place the solvent particles randomly into the system, and assign all
	// of the associated quantities (velocity,etc.) to each particle.
	for (i=k; i<n_dpd; i++) {
		// Random position.
		x[k]  = ran3(&iseed) * Lx;
		y[k]  = ran3(&iseed) * Ly;
		z[k]  = ran3(&iseed) * Lz;

		// Random velocity at T.
		theta = ran3(&iseed) * PI;
		phi   = ran3(&iseed) * 2.0 * PI;
		vx[k] = sqrt(T) * cos(phi) * sin(theta);
		vy[k] = sqrt(T) * sin(phi) * sin(theta);
		vz[k] = sqrt(T) * cos(theta);
		vxt  += vx[k];
		vyt  += vy[k];
		vzt  += vz[k];
		
		// Particle type. Let's call solvent -1.
		p_type[k] = -1;

		k++;
	}

	// Subtract the average velocity from  *all* particles to eliminate 
        // any drift that may be present.
	vxt /= (double) n_dpd;
	vyt /= (double) n_dpd;
	vzt /= (double) n_dpd;
	for (i=0; i<n_dpd; i++) {
		vx[i] -= vxt;
		vy[i] -= vyt;
		vz[i] -= vzt;
	}
	
	// Allocate the cell list for processing interactions.
	return;
}


void read_parameters(char *filename) {
	FILE *param_file;
	char junk[256]; 

	param_file = fopen(filename, "r");
	if (param_file == NULL) {
		printf("Parameter file %s not found.\n", filename);
		exit(1);
	}
	
	// Read in the "junk" comment from the file.
	fgets(junk, 256, param_file);

	// Read in the value that follows.
	fscanf(param_file, "%lf\n", &Lx);

	// Repeat the same process for the remaining lines of the file.
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

	// Define gamma from sigma, since they are related through the
	// fluctuation-dissipation theorem:
	gamma = pow(sigma, 2.0) / (2.0 * T);
	sigma = sigma * sqrt(3.0 / T);
	return;
}