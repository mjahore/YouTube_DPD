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

#include "cell_list.h" // Cell list variables/parameters.

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
void   calc_forces(void);
double calc_Rg(void);
void   initialize_system(void);
void   make_list(void);
void   read_parameters(char *filename);
void   write_xyz(char *filename, int append);

int main(int argc, char *argv[]) {
	double Rg_rms;

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
	make_list();
	calc_forces();

	// Calculate initial rms Rg.
	Rg_rms = calc_Rg();

	printf("Rg: %2.2lf r_{c}\n", Rg_rms);

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
	free(cell_head);
	free(cell_list);
	return 0;
}

void calc_forces(void) {
	int i,j;
	int idx_i, idx_j;
	double dx, dy, dz, dr;
	double f_net, f_comp; 

	// Initialize forces
	memset(fx, 0.0, n_dpd*sizeof(double));
	memset(fy, 0.0, n_dpd*sizeof(double));
	memset(fz, 0.0, n_dpd*sizeof(double));

	// Reset energy.
	U     = 0.0;
	K     = 0.0;
	b_avg = 0.0;

	// Apply the harmonic force between all adjacent monomers.
	for (i=0; i<n_poly; i++) {
		for(j=0; j<N-1; j++) {
			idx_i = i*N + j;
			idx_j = idx_i + 1;

			// How apart are monomers idx_i and idx_j?
			dx = rx[idx_i] - rx[idx_j];
			dy = ry[idx_i] - ry[idx_j];
			dz = rz[idx_i] - rz[idx_j];

			dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

			b_avg += dr;

			// Calculate the magnitude of the harmonic force.
			f_net = -kH * (1.0 - b/dr);

			// Add to our potential energy.
			U += 0.5 * kH * pow((dr-b), 2);

			// Add each component of the force to the associated
			// force arrays.
			f_comp     = f_net * dx;
			fx[idx_i] += f_comp;
			fx[idx_j] -= f_comp;

			f_comp     = f_net * dy;
			fy[idx_i] += f_comp;
			fy[idx_j] -= f_comp;

			f_comp     = f_net * dz;
			fz[idx_i] += f_comp;
			fz[idx_j] -= f_comp;
		}
	}


	// Process all cells & their neighbors and calculate the conservative
	// random, and dissipative forces (NEXT TIME!).

	return;
}

double calc_Rg(void) {
	int i,j,k;
	int idx_i, idx_j;
	double dx, dy, dz, dr2;
	double Rg;

	Rg = 0.0;

	for(i=0; i<n_poly; i++) {
		for(j=0; j<N; j++) {
			idx_i = i*N + j;
			for(k=j+1; k<N; k++) {
				idx_j = i*N + k;

				dx  = rx[idx_i] - rx[idx_j];
				dy  = ry[idx_i] - ry[idx_j];
				dz  = rz[idx_i] - rz[idx_j];
				dr2 = pow(dx,2) + pow(dy,2) + pow(dz,2);

				Rg += dr2;
			}
		}
	}

	// We need to divide Rg by N^2
	Rg /= pow(N,2);

	// Divide by the number of polymers.
	Rg /= n_poly;

	Rg = sqrt(Rg); // Calculate root mean-squared.

	return Rg;
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
		vx[k] = sqrt(3.0*T) * cos(phi) * sin(theta);
		vy[k] = sqrt(3.0*T) * sin(phi) * sin(theta);
		vz[k] = sqrt(3.0*T) * cos(theta);
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
			x[k]  = x[k-1] + b * cos(phi) * sin(theta);
			y[k]  = y[k-1] + b * sin(phi) * sin(theta);
			z[k]  = z[k-1] + b * cos(theta);

			// Assign a velocity to this particle.
			theta = ran3(&iseed) * PI;
			phi   = ran3(&iseed) * 2.0 * PI;
			vx[k] = sqrt(3.0*T) * cos(phi) * sin(theta);
			vy[k] = sqrt(3.0*T) * sin(phi) * sin(theta);
			vz[k] = sqrt(3.0*T) * cos(theta);
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
		vx[k] = sqrt(3.0*T) * cos(phi) * sin(theta);
		vy[k] = sqrt(3.0*T) * sin(phi) * sin(theta);
		vz[k] = sqrt(3.0*T) * cos(theta);
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

	// Allocate the cell list for processing interactions.
	nx     = (int) (Lx/CELL_SIZE);
	ny     = (int) (Ly/CELL_SIZE);
	nz     = (int) (Lz/CELL_SIZE);
	n_cell = nx*ny*nz;

	cell_head = (int*) malloc(n_cell*sizeof(int));
	cell_list = (int*) malloc(n_dpd*sizeof(int));

	// BOOM WE ARE DONE!
	return;
}

void make_list(void) {
	int i, cell_id, tmp_id;

	// Initialize our linked list.
	memset(cell_head, -1, n_cell*sizeof(int));
	memset(cell_list, -1, n_dpd*sizeof(int));

	// Loop through all DPD particles and place them after
	// applying the PBC.
	for(i=0; i<n_dpd; i++) {
		// Apply PBC:
		if (x[i] > Lx) {
			x[i] -= Lx;
		} else if (x[i] < 0) {
			x[i] += Lx;
		}

		if (y[i] > Ly) {
			y[i] -= Ly;
		} else if (y[i] < 0) {
			y[i] += Ly;
		}
		
		if (z[i] > Lz) {
			z[i] -= Lz;
		} else if (z[i] < 0) {
			z[i] += Lz;
		}

		// Calculate the cell to which particle i belongs:
		cell_id = (int) (x[i]/CELL_SIZE) + (int) (y[i]/CELL_SIZE)*nx + (int) (z[i]/CELL_SIZE)*ny*nx;
	
		tmp_id = cell_head[cell_id];
		cell_head[cell_id] = i;
		cell_list[i] = tmp_id;
	}
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
