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
double b_avg;            // Average bond length (from simulation)
double T;                // Temperature (actually, kB*T)
double U;                // Total system potential energy.
double K;                // Total system kinetic energy.
double a_ii;             // Same-type interaction
double a_ij;             // Different-type interaction
double sigma, gamma1;    // Random and dissipative force strengths (gamma is taken)
double kH;               // Harmonic potential strength
double dt;               // Integration time step
int    time_equil;       // Number of equilibration time steps
int    time_prod;        // Number of production time steps
int    t;                // Current time step
int    n_dpd;            // Total number of particles
int    n_poly;           // Total number of polymers

// Random Number Generator
int inext,inextp;      // Ran3
long ma[56];           // Ran3
int iff=0;             // Ran3
int iseed;             // Random seed
#include "ran3.c"      // ran3() from Numerical Recipes in C

#include "cell_list.h" // Cell list for finding interacting pairs

#define PI 3.1415926535897932384626433832795f
#define CALC_FREQ 10
#define XYZ_FREQ  100

// Function prototypes:
void   calc_forces(void);
double calc_Rg(void);
double calc_p(void);
double calc_T(void);
void   initialize_system(void);
void   make_list(void);
void   read_parameters(char *filename);
void   velocity_verlet(double lambda);
void   write_xyz(char *filename, int append);

int main(int argc, char *argv[]) {
	FILE *output;
	double kBT_temp;

	// Cosmetic stuff.
	printf("\n");

	// Set ran3() random seed based on current time.
	iseed=-time(0);

	// Read simulation parameters.
	read_parameters("param.cfg");

	// Just some info to make sure things are okay...
	printf("Simulating a %lf x %lf x %lf system with density rho = %lf.\n", Lx, Ly, Lz, rho);

	// Set up the initial configuration of the system
	initialize_system();

	// Write out our initial configuration:
	write_xyz("initial.xyz", 0);

	// Calculate the initial forces before entering the main loop
	// of the simulation.
	make_list();
	calc_forces();

	// Equilibrate the system.
	output = fopen("dpd_data.dat", "w");
	fprintf(output, "# timestep\t kBT\t p\t U\t K\t b_avg\t Rg\n");
	fclose(output);
	for(t=0; t<time_equil; t++) {
		velocity_verlet(0.5);

		// Output to disk at regular intervals.
		if (t % CALC_FREQ == 0) {
			kBT_temp = calc_T();
			printf("(%d) kBT = %1.5lf\t p = %lf\t U = %1.5lf\t K = %1.5lf\t Rg: %1.5lf\n", t, kBT_temp, calc_p(), U, K, calc_Rg()); 
			output = fopen("dpd_data.dat", "aw");
			fprintf(output, "%d   %lf   %lf   %lf   %lf   %lf   %lf\n", t, kBT_temp, calc_p(), U, K, b_avg, calc_Rg());
			fclose(output);
		}
	}

	// Run production steps.
	for(t=time_equil; t<(time_equil+time_prod); t++) {
		velocity_verlet(0.5);

		// Output to disk at regular intervals.
		if (t % CALC_FREQ == 0) {
			kBT_temp = calc_T();
			printf("(%d) kBT = %1.5lf\t p = %lf\t U = %1.5lf\t K = %1.5lf\t Rg: %1.5lf\n", t, kBT_temp, calc_p(), U, K, calc_Rg()); 
			output = fopen("dpd_data.dat", "aw");
			fprintf(output, "%d   %lf   %lf   %lf   %lf   %lf   %lf\n", t, kBT_temp, calc_p(), U, K, b_avg, calc_Rg());
			fclose(output);
		}

		if (t % XYZ_FREQ == 0) {
			// Write coordinates to .xyz file for viewing.
			write_xyz("trajectory.xyz", 1);
		}
	}

	// Simulation complete.
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
	free(p_type);
	free(cell_head);
	free(cell_list);
	return 0;
}

// Here we calculate the harmonic force between all adjacent monomers, and then calculate
// the three DPD forces for all particles that are within a distance r_{c}.
//
void calc_forces (void) {
	int i, j, k, l, idx_i, idx_j;
	int ni, nj, nk;
	int cell0, cell1;
	int p_ij;
	double dx, dy, dz, dr;
	double x_trans, y_trans, z_trans;
	double f_net, f_comp, fC, fR, fD;
	double a, theta_ij, r_dot_v, wr, wd;

	// Initialize the forces.
	memset(fx, 0.0, n_dpd*sizeof(double));
	memset(fy, 0.0, n_dpd*sizeof(double));
	memset(fz, 0.0, n_dpd*sizeof(double));

	// Reset energy.
	U     = 0.0;
	K     = 0.0;
	b_avg = 0.0;

	// Apply the harmonic force between all adjacent monomers.
	for (i=0; i<n_poly; i++) {
		// j < N-1 is due to the fact that the last monomer has no adjacent monomer in the
		// forward direction!
		for (j=0; j<N-1; j++) {
			idx_i = i*N+j;
			idx_j = idx_i + 1;

			// How far apart are the adjacent monomers?
			dx = rx[idx_i] - rx[idx_j];
			dy = ry[idx_i] - ry[idx_j];
			dz = rz[idx_i] - rz[idx_j];
			dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));

			// Track average bond length.
			b_avg += dr;

			// Harmonic potential has no cutoff, so regardless of
			// the value of dr, we compute the force.
			// U = 1/2 * kH * (dr - b)^2.
			// The b/dr term comes from pre-computing part of the projection
			// of f_net into the 3 directions. 
			f_net = -kH * (1.0 - b/dr); 
	
			// Add to our potential energy:
			U += 0.5 * kH * pow((b - dr), 2);

			// Do each component of the force. Newton's 3rd law!
			f_comp     = f_net*dx;
			fx[idx_i] += f_comp;
			fx[idx_j] -= f_comp;

			f_comp     = f_net*dy;
			fy[idx_i] += f_comp;
			fy[idx_j] -= f_comp;

			f_comp     = f_net*dz;
			fz[idx_i] += f_comp;
			fz[idx_j] -= f_comp;
		}
	}

	b_avg /= (n_poly * (N - 1));
	
	// Now we process *all* particles in the system to compute the three
	// DPD forces. Since all forces are pair-wise we really only need to
	// examine 1/2 of the particles in the system. 
	for(k=0; k<nz; k++) {
		for(j=0; j<ny; j++) {
			for(i=0; i<nx; i++) {
				// Calculate this cell's index.
				cell0 = i + j*nx + k*nx*ny;

			        // Loop through the 14 cells we need to examine (13 neighbors + cell0)
			  	for(l=0; l<14; l++) {
					// Who is in this cell?
					idx_i = cell_head[cell0];
					
					// Compute the index of the cell we're looking at.
					ni = i + nix[l];
					nj = j + niy[l];
					nk = k + niz[l];	

					// If we're on a boundary, we need to change our
					// indices to account for the PBC.
					x_trans = 0.0;
					y_trans = 0.0;
					z_trans = 0.0;
					if (ni > nx-1) {
						ni      = 0;
				       		x_trans = -Lx;
				 	} else if (ni < 0) {
						ni      = nx-1;
						x_trans = Lx;
					}	

					if (nj > ny-1) {
						nj      = 0;
				       		y_trans = -Ly;
				 	} else if (nj < 0) {
						nj      = ny-1;
						y_trans = Ly;
					}	

					if (nk > nz-1) {
						nk      = 0;
				       		z_trans = -Lz;
				 	} else if (nk < 0) {
						nk      = nz-1;
						z_trans = Lz;
					}		
	
					// Get the top of our neighboring cell.
					cell1 = ni + nj*nx + nk*nx*ny;

					// Loop through the contents of cell0 and process with
					// all of cell1 until we reach the bottom!
					idx_j = cell_head[cell1];
					if (cell0 == cell1) idx_j = cell_list[idx_i];

			        	while (idx_i > -1) {
						while (idx_j > -1) {
							dx = x[idx_i] - x[idx_j] + x_trans;
							dy = y[idx_i] - y[idx_j] + y_trans;
							dz = z[idx_i] - z[idx_j] + z_trans;
	
							dr = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
	
							// If the particles are closer than dr = r_{c}, calculate
							// the forces.
							if (dr < 1.0) {
								p_ij = p_type[idx_i]*p_type[idx_j];
	
								// Choose interaction strength for conservative force.
								switch(p_ij) {
									case 1:
										a = a_ii;
										break;
									case -1:
										a = a_ij;
										break;
								}
					
								// Weighting terms for (r)andom and (d)isspative forces.
								wr = 1.0/dr - 1.0;
								wd = pow(wr, 2);
							
								// Random number for random force, zero mean, unit variance.
								theta_ij = 2.0*ran3(&iseed) - 1.0;
	
								// Inner product between r_ij and  v_ij.
								r_dot_v = dx*(vx[idx_i] - vx[idx_j]) + dy*(vy[idx_i] - vy[idx_j]) + dz*(vz[idx_i] - vz[idx_j]);
								// Add to our potential energy:
								U += 0.5*a*pow((1.0-dr), 2);
	
								// Compute our forces.
								fC = a * wr;
								fR = sigma * wr * theta_ij;
								fD = -gamma1 * wd * r_dot_v;
	
								f_net = fC + fR + fD;
	
								// Compute the components.
								f_comp = f_net * dx;
								fx[idx_i] += f_comp;
								fx[idx_j] -= f_comp;
	
								f_comp = f_net * dy;
								fy[idx_i] += f_comp;
								fy[idx_j] -= f_comp;
	
								f_comp = f_net * dz;
								fz[idx_i] += f_comp;
								fz[idx_j] -= f_comp;
							} // end if (dr < 1)
							// Get next particle in the cell.
							idx_j = cell_list[idx_j];
						} // end while (idx_j > -1)
						idx_j = cell_head[cell1];
						idx_i = cell_list[idx_i]; 
						if (cell0 == cell1) idx_j = cell_list[idx_i];
				    	} // end while (idx_i > -1)
				} // end for l = [0, 14);		
			} // end for i
		} // end for j
	} // end for k

	return;
}

// Here, we calculate p^2 = p (dot) p and assume the mass of
// all particles to be m = 1 (not strictly necessary).
double calc_p(void) {
	int i;
	double p, px, py, pz;

	px = 0;	
	py = 0;
	pz = 0;
	for (i=0; i<n_dpd; i++) {
		px += vx[i];
		py += vy[i];
		pz += vz[i];
	}
	p = pow(px,2) + pow(py,2) + pow(pz,2);

	return p;
}

// Here, we calculate the root mean-squared radius of gyration of the
// polymers.
double calc_Rg(void) {
	int i,j,k;
	int idx_i, idx_j;
	double dx, dy, dz, dr2;
	double Rg;

	Rg = 0.0;

	// There are a few ways to calculate Rg. I'm using the simplified
	// pair-wise expression derived in Polymer Physics (Rubinstein & Colby)
	// Equation 2.48.
	//
	// You could do it according to the formal definition (Eq. 2.44) or
	// you could compute the gyration tensor and find the eigenvalues.
	for (i=0; i<n_poly; i++) {
		for (j=0; j<N; j++) {
			idx_i = i*N + j;
			for (k=j+1; k<N; k++) {
				idx_j = i*N + k;
				dx    = rx[idx_i] - rx[idx_j];
				dy    = ry[idx_i] - ry[idx_j];
				dz    = rz[idx_i] - rz[idx_j];
				dr2   = pow(dx,2) + pow(dy,2) + pow(dz,2);

				Rg += dr2;
			}
		}
	}

	// We need to divide Rg by N^2, and by the number of polymers.
	Rg /= pow(N,2);
	Rg /= n_poly;

	// And take square root since above we calculated Rg^2.
	Rg = sqrt(Rg);

	return Rg;
}


// Here, we calculate the temperature T from the kinetic energy K and
// the degrees of freedom in the system (3 * number of particles - 3).
//
double calc_T(void) {
	int i, deg_freedm;
	double kBT = 0;

	deg_freedm = 3 * n_dpd - 3;

	// This is the explicit calculation.
	/*
	for (i=0; i<n_dpd; i++) {
		kBT += pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2);
	}
	kBT /= deg_freedm;
	*/
	// This is how you'd calculate it directly from the kinetic energy.
	 kBT = 2.0 * K / deg_freedm;
	return kBT;
}

// Here, we create the random initial condition for the simulation. Polymers and solvent
// particles are placed randomly into the system along with random velocities. The average
// velocity is subtracted from all particles to eliminate any net drift that may be present.
//
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
	x      = (double*) malloc(n_dpd * sizeof(double));
	y      = (double*) malloc(n_dpd * sizeof(double));
	z      = (double*) malloc(n_dpd * sizeof(double));
	rx     = (double*) malloc(n_dpd * sizeof(double));
	ry     = (double*) malloc(n_dpd * sizeof(double));
	rz     = (double*) malloc(n_dpd * sizeof(double));
	vx     = (double*) malloc(n_dpd * sizeof(double));
	vy     = (double*) malloc(n_dpd * sizeof(double));
	vz     = (double*) malloc(n_dpd * sizeof(double));
	vxp    = (double*) malloc(n_dpd * sizeof(double));
	vyp    = (double*) malloc(n_dpd * sizeof(double));
	vzp    = (double*) malloc(n_dpd * sizeof(double));
	fx     = (double*) malloc(n_dpd * sizeof(double));
	fy     = (double*) malloc(n_dpd * sizeof(double));
	fz     = (double*) malloc(n_dpd * sizeof(double));
	fxp    = (double*) malloc(n_dpd * sizeof(double));
	fyp    = (double*) malloc(n_dpd * sizeof(double));
	fzp    = (double*) malloc(n_dpd * sizeof(double));
	p_type = (int*)    malloc(n_dpd * sizeof(int));

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
		// phi are two random angles in spherical coordinates (i.e., the 		
		// azimuthal and polar angles).
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

	// Copy the coordinates to our non-PBC arrays.
	memcpy(rx, x, n_dpd*sizeof(double));
	memcpy(ry, y, n_dpd*sizeof(double));
	memcpy(rz, z, n_dpd*sizeof(double));

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
	nx     = (int) (Lx/CELL_SIZE);
	ny     = (int) (Ly/CELL_SIZE);
	nz     = (int) (Lz/CELL_SIZE);
	n_cell = nx*ny*nz;

	cell_head = (int*) malloc(n_cell*sizeof(int));
	cell_list = (int*) malloc(n_dpd*sizeof(int));

	return;
}

// Here we apply the periodic boundary condition (PBC) and place all particles
// in the system into our linked list so that we can quickly determine all interacting
// particle pairs within a distance of r_{c}.
//
void make_list(void) {
	int i, cell_id, tmp_id;

	// Clear out the cell list. -1 means empty cell.
	memset(cell_head, -1, n_cell*sizeof(int));
	memset(cell_list, -1, n_dpd*sizeof(int));

	// Loop through all DPD particles and place them.
	for (i=0; i<n_dpd; i++) {

		// Apply periodic boundary condition.
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

		// See where the particle is, and add to the linked list. Everybody goes somewhere!
		cell_id = (int) (x[i]/CELL_SIZE) + (int) (y[i]/CELL_SIZE)*nx + (int)(z[i]/CELL_SIZE)*nx*ny;
		tmp_id             = cell_head[cell_id];
		cell_head[cell_id] = i;
		cell_list[i]       = tmp_id;
	}

	return;
}

// Here we read in the parameters for our simulation from the config file we specify in filename.
//
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
	gamma1 = pow(sigma, 2.0) / (2.0 * T);
	sigma  = sigma * sqrt(3.0 / dt);
	return;
}

// Below is the implementation of the velocity-Verlet
// algorithm as detailed in Groot & Warren.
//
void velocity_verlet(double lambda) {
	int i;
	double dt2 = pow(dt,2);

	// We need to store the current values of the forces
	// and velocities to be used in the 2nd half of the 
	// velocity-Verlet algorithm.
	memcpy(fxp, fx, n_dpd*sizeof(double));
	memcpy(fyp, fy, n_dpd*sizeof(double));
	memcpy(fzp, fz, n_dpd*sizeof(double));
	memcpy(vxp, vx, n_dpd*sizeof(double));
	memcpy(vyp, vy, n_dpd*sizeof(double));
	memcpy(vzp, vz, n_dpd*sizeof(double));

	// Integrate equations of motion to update position & velocity.
	// Lambda is a parameter from Groot & Warren that accounts
	// for some effects of the stochastic processes. It's usually
	// set to 0.5.
	for (i=0; i<n_dpd; i++) {
		// Update positions (w/ PBC)
		x[i] += dt*vx[i] + 0.5*dt2*fx[i];
		y[i] += dt*vy[i] + 0.5*dt2*fy[i];
		z[i] += dt*vz[i] + 0.5*dt2*fz[i];

		// Update positions (w/o PBC)
		rx[i] += dt*vx[i] + 0.5*dt2*fx[i];
		ry[i] += dt*vy[i] + 0.5*dt2*fy[i];
		rz[i] += dt*vz[i] + 0.5*dt2*fz[i];

		// Predict velocities (will be corrected in 2nd half
		// of velocity-Verlet algorithm)
		vx[i] += dt*lambda*fx[i];
		vy[i] += dt*lambda*fy[i];
		vz[i] += dt*lambda*fz[i];
	}

	// Now we calculate the forces again:
	make_list();
	calc_forces();

	// Correct velocities and calculate kinetic energy
	for (i=0; i<n_dpd; i++) {
		vx[i] = vxp[i] + 0.5*dt*(fx[i]+fxp[i]);
		vy[i] = vyp[i] + 0.5*dt*(fy[i]+fyp[i]);
		vz[i] = vzp[i] + 0.5*dt*(fz[i]+fzp[i]);
	
		K += 0.5*(pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
	}

	return;
}

// Finally, here's where we write the polymer coordinates (only!) to disk in 
// .xyz format so that VMD, Ovito, etc. can read them to render our snapshot.
//
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
	for (i=0; i<n_poly*N; i++) {
		fprintf(output, "C %3.6lf %3.6lf %3.6lf\n", x[i], y[i], z[i]);
	}
	fclose(output);
	return;
} 
