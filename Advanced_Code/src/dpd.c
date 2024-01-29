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
#define MAX_COMPONENTS 2
#define BUF_SIZE 65536
#define CALC_FREQ 100
#define OUT_FREQ  1000

// Global Variables 
double *x, *y, *z;              // Coordinates w/ PBC
double *rx, *ry, *rz;           // Coordinates w/o PBC
double *vx, *vy, *vz;           // Velocities (current time step)
double *vxp, *vyp, *vzp;        // Velocities (previous time step)
double *fx, *fy, *fz;           // Forces (current time step)
double *fxp, *fyp, *fzp;        // Forces (previous time step)
int    *p_type;                 // Particle type
int    *p_comp;                 // Particle component
double Lx, Ly, Lz;              // System dimensions
double rho;                     // Fluid density
double phi[MAX_COMPONENTS];     // Polymer volume fraction
int    N[MAX_COMPONENTS];       // Polymer deg. of polymerization
double b0[MAX_COMPONENTS];      // Polymer bond length (harmonic or r_eq)
double b1[MAX_COMPONENTS];      // Polymer bond length (FENE - r_max)
double b_avg[MAX_COMPONENTS];   // Average bond length (from simulation)
double T;                       // Temperature (actually, kB*T)
double U;                       // Total system potential energy.
double K;                       // Total system kinetic energy.
double *a;                      // Conservative interaction strengths.
double sigma, gamma1;           // Random and dissipative force strengths (gamma is taken)
double kH[MAX_COMPONENTS];      // Harmonic potential strength
double dt;                      // Integration time step
int    time_equil;              // Number of equilibration time steps
int    time_prod;               // Number of production time steps
int    t;                       // Current time step
int    n_dpd;                   // Total number of particles
int    n_poly[MAX_COMPONENTS];  // Total number of polymers

// Architecture configuration
char   arch_names[MAX_COMPONENTS][BUF_SIZE];
char   arch_files[MAX_COMPONENTS][BUF_SIZE];

// Random Number Generator
#include "ran3.c"      // ran3() from Numerical Recipes in C

#include "nlist.h"     // Cell list for finding interacting pairs
                       // and for intrachain bonding.

// Prototypes for main functions.
void   calc_forces(int prod);
void   calc_intrachain(int comp, int idx);
double calc_Rg(int comp);
double calc_p(void);
double calc_T(void);
void   initialize_system(void);
int    line_count(FILE *fp);
int    make_arch(int comp, int idx, double phi_comp, char *filename);
void   make_list(void);
void   read_interactions(char *filename);
void   read_parameters(char *filename);
void   velocity_verlet(int prod, double lambda);
void   write_xyz(int comp, char *filename, int append);

// Prototypes for forces
double f_harmonic(double, double, double);   

int main(int argc, char *argv[]) {
	FILE *data_output;
	double T_calc, p_calc, Rg_rms;

	// Cosmetic stuff.
	printf("\n");
	printf("dpd.c by Michael J. A. Hore, Case Western Reserve University\n");
	printf("Copyright (c) 2024 Michael J. A. Hore <hore@case.edu>\n\n");
	printf("Version: %s  (Compiled: %s)\n\n", VERSION, BUILD_INFO);
	
	// We pass the parameter file as an input.
	if (argc < 2) {
		printf("\nUsage: ./gDPD [PARAMETER FILE]\n\n");
		return -1;
	}

	// Set ran3() random seed based on current time.
	iseed=-time(0);

	// Read simulation parameters.
	read_parameters(argv[1]);

	// Set up the initial configuration of the system
	initialize_system();

	// Write out our initial configuration:
	write_xyz(0, "initial0.xyz", 0);
	write_xyz(1, "initial1.xyz", 0);

	// Calculate the initial forces before entering the main loop
	// of the simulation.
	calc_forces(0);

	// Equilibrate the system.
	for(t=0; t<time_equil; t++) {
		velocity_verlet(0,0.5);

		// Output to disk at regular intervals.
		if (t % CALC_FREQ == 0) {
			Rg_rms = calc_Rg(0);
			T_calc = calc_T();
			p_calc = calc_p();

			printf("(%d)\t kBT: %2.5lf   |p|: %2.5lf    <Rg^2>^(1/2): %2.5lf\n", t, T_calc, p_calc, Rg_rms);

			// Output temperature to disk.
			data_output = fopen("temperature.dat", "aw");
			fprintf(data_output, "%d\t%lf\n", t, T_calc);
			fclose(data_output);

			// Output momentum to disk.
			data_output = fopen("momentum.dat", "aw");
			fprintf(data_output, "%d\t%lf\n", t, p_calc);
			fclose(data_output);

			// Output Rg to disk.
			data_output = fopen("radius_gyration.dat", "aw");
			fprintf(data_output, "%d\t%lf\t%lf\n", t, Rg_rms, b_avg[0]);
			fclose(data_output);

			// Output kinetic/potential energy.
			data_output = fopen("energy.dat", "aw");
			fprintf(data_output, "%d\t%lf\t%lf\n", t, K, U);
			fclose(data_output);
		}

		if (t % OUT_FREQ == 0) {
			write_xyz(0,"trajectory0.xyz", 1);
			write_xyz(1,"trajectory1.xyz", 1);
		}
	}

	// Run production steps.
	for(t=time_equil; t<(time_equil+time_prod); t++) {
		velocity_verlet(1,0.5);

		// Output to disk at regular intervals.
		if (t % CALC_FREQ == 0) {
			Rg_rms = calc_Rg(0);
			T_calc = calc_T();
			p_calc = calc_p();

			printf("(%d)\t kBT: %2.5lf   |p|: %2.5lf    <Rg^2>^(1/2): %2.5lf\n", t, T_calc, p_calc, Rg_rms);

			// Output temperature to disk.
			data_output = fopen("temperature.dat", "aw");
			fprintf(data_output, "%d\t%lf\n", t, T_calc);
			fclose(data_output);

			// Output momentum to disk.
			data_output = fopen("momentum.dat", "aw");
			fprintf(data_output, "%d\t%lf\n", t, p_calc);
			fclose(data_output);

			// Output Rg to disk.
			data_output = fopen("radius_gyration.dat", "aw");
			fprintf(data_output, "%d\t%lf\t%lf\n", t, Rg_rms, b_avg[0]);
			fclose(data_output);

			// Output kinetic/potential energy.
			data_output = fopen("energy.dat", "aw");
			fprintf(data_output, "%d\t%lf\t%lf\n", t, K, U);
			fclose(data_output);
		}

		if (t % OUT_FREQ == 0) {
			write_xyz(0,"trajectory0.xyz", 1);
			write_xyz(1,"trajectory1.xyz", 1);
		}
	}

	// Simulation complete.
	// Cosmetic stuff.
	printf("\n");

	// Free the memory that was allocated in initialize_system().
	free(a);
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
	free(p_comp);
	free(cell_head);
	free(cell_list);
	free(bond_list);
	return 0;
}

//
// calc_intrachain - Calculates pair-wise bonding interactions using bond_list[].
//                   This routine should work regardless of the polymer architecture
//                   because all bonds are defined by the indices contained within
//                   bond_list[] without making any assumptions about adjacent monomers.
//
void calc_intrachain (int comp, int i) {
	int j;
	double dx, dy, dz, dr;
	double f_net, f_comp;

	// This intrachain bonding routine should work regardless of the polymer
	// architecture because all bonds are defined by the indices contained within
	// bond_list[] and not any assumptions about adjacent monomers.
	j  = bond_list[i];
	if (j > -1) {
		dx           = rx[i] - rx[j];
		dy           = ry[i] - ry[j];
		dz           = rz[i] - rz[j];
		dr           = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
		b_avg[comp] += dr;

		// We use the harmonic potential for bonding. Replace
		// with a suitable function to use FENE, etc.
		f_net  = f_harmonic(kH[comp], b0[comp], dr);
	
		// Add the components of the force to the force arrays.
		f_comp = f_net * dx;
		fx[i] += f_comp;
		fx[j] -= f_comp;
			
		f_comp = f_net * dy;
		fy[i] += f_comp;
		fy[j] -= f_comp;

		f_comp = f_net * dz;
		fz[i] += f_comp;
		fz[j] -= f_comp;
	}
	return;
}

// Here we calculate the harmonic force between all adjacent monomers, and then calculate
// the three DPD forces for all particles that are within a distance r_{c}.
//
void calc_forces (int prod) {
	int i, j, k, l, idx_i, idx_j;
	int ni, nj, nk;
	int cell0, cell1;
	int p_ij;
	double dx, dy, dz, dr;
	double x_trans, y_trans, z_trans;
	double f_net, f_comp, fC, fR, fD;
	double theta_ij, r_dot_v, wr, wd;


	// Place our particles into a cell list, and process any
	// intrachain interactions. Average bond length b_avg is
	// also calculated within make_list().
	make_list();

	for (i=0; i<MAX_COMPONENTS; i++) {
		b_avg[i] /= (n_poly[i] * (N[i] - 1));
	}
	
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
								// Determine interaction type (below):
								p_ij = p_type[idx_i]*p_type[idx_j]*prod;
	
								// Weighting terms for (r)andom and (d)isspative forces.
								wr = 1.0/dr - 1.0;
								wd = pow(wr, 2);
							
								// Random number for random force, zero mean, unit variance.
								theta_ij = 2.0*ran3(&iseed) - 1.0;
	
								// Inner product between r_ij and  v_ij.
								r_dot_v = dx*(vx[idx_i] - vx[idx_j]) + dy*(vy[idx_i] - vy[idx_j]) + dz*(vz[idx_i] - vz[idx_j]);
								// Add to our potential energy:
								U += 0.5*a[p_ij]*pow((1.0-dr), 2);
	
								// Compute our forces.
								fC = a[p_ij] * wr;
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
double calc_Rg(int comp) {
	int i,j,k;
	int idx_i, idx_j;
	int offset;
	double dx, dy, dz, dr2;
	double Rg;

	Rg = 0.0;

	// Determine offset for arrays.
	offset = 0;
	for (i=0; i<comp; i++) {
		offset += n_poly[i]*N[i];
	}

	// There are a few ways to calculate Rg. I'm using the simplified
	// pair-wise expression derived in Polymer Physics (Rubinstein & Colby)
	// Equation 2.48.
	//
	// You could do it according to the formal definition (Eq. 2.44) or
	// you could compute the gyration tensor and find the eigenvalues.
	for (i=offset; i<n_poly[comp]; i++) {
		for (j=0; j<N[comp]; j++) {
			idx_i = i*N[comp] + j + offset;
			for (k=j+1; k<N[comp]; k++) {
				idx_j = i*N[comp] + k + offset;
				dx    = rx[idx_i] - rx[idx_j];
				dy    = ry[idx_i] - ry[idx_j];
				dz    = rz[idx_i] - rz[idx_j];
				dr2   = pow(dx,2) + pow(dy,2) + pow(dz,2);

				Rg += dr2;
			}
		}
	}

	// We need to divide Rg by N^2, and by the number of polymers.
	Rg /= pow(N[comp],2);
	Rg /= n_poly[comp];

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

	// Just calculate directly from kinetic energy and the equipartition 
	// theorem.
	 kBT = 2.0 * K / deg_freedm;

	return kBT;
}

//
// f_harmonic - Force from harmonic potential (normalized by dr).
//
double f_harmonic(double k, double r0, double dr) {
	double f_mag;
	
	// Magnitude of the force, normalized by dr.
	f_mag = -k * (1.0 - r0/dr);

	// U is a global potential energy variable. Add to it.
	U += 0.5 * k * pow((dr - r0), 2.0);

	return f_mag;
}


// Here, we create the random initial condition for the simulation. Polymers and solvent
// particles are placed randomly into the system along with random velocities. The average
// velocity is subtracted from all particles to eliminate any net drift that may be present.
//
void initialize_system(void) {
	int i,j,k;
	double sc_phi, sc_theta;     // Spherical coordinates.
	double vxt, vyt, vzt;  // Average velocity in x,y,z directions.

	// Calculate the number of DPD particles in the simulation
	// box.
	n_dpd = (int) (rho * Lx * Ly * Lz);

	// Now, allocate memory in the arrays that store particle information
	// (x, y, z, fx, fy, ...). This memory is freed at the end of the
	// main() portion of the code.
	x         = (double*) malloc(n_dpd * sizeof(double));
	y         = (double*) malloc(n_dpd * sizeof(double));
	z         = (double*) malloc(n_dpd * sizeof(double));
	rx        = (double*) malloc(n_dpd * sizeof(double));
	ry        = (double*) malloc(n_dpd * sizeof(double));
	rz        = (double*) malloc(n_dpd * sizeof(double));
	vx        = (double*) malloc(n_dpd * sizeof(double));
	vy        = (double*) malloc(n_dpd * sizeof(double));
	vz        = (double*) malloc(n_dpd * sizeof(double));
	vxp       = (double*) malloc(n_dpd * sizeof(double));
	vyp       = (double*) malloc(n_dpd * sizeof(double));
	vzp       = (double*) malloc(n_dpd * sizeof(double));
	fx        = (double*) malloc(n_dpd * sizeof(double));
	fy        = (double*) malloc(n_dpd * sizeof(double));
	fz        = (double*) malloc(n_dpd * sizeof(double));
	fxp       = (double*) malloc(n_dpd * sizeof(double));
	fyp       = (double*) malloc(n_dpd * sizeof(double));
	fzp       = (double*) malloc(n_dpd * sizeof(double));
	p_type    = (int*)    malloc(n_dpd * sizeof(int));
	p_comp    = (int*)    malloc(n_dpd * sizeof(int));
	bond_list = (int*)    malloc(n_dpd * sizeof(int));

	// Put in component 1 using the defined architecture.
	j = 0;
	k = 0;
	for (i=0; i<MAX_COMPONENTS; i++) {
		j  = make_arch(i, j, phi[i], arch_files[i]);
		k += j;
	}

	// Rescale n_dpd based on the actual number of DPD particles
	// placed into the system.
	n_dpd = k;
	printf("Simulating %d DPD particles in system of %3.2lf x %3.2lf x %3.2lf (rho = %1.2lf)\n", n_dpd, Lx, Ly, Lz, rho);

	// Copy the coordinates to our non-PBC arrays.
	memcpy(rx, x, n_dpd*sizeof(double));
	memcpy(ry, y, n_dpd*sizeof(double));
	memcpy(rz, z, n_dpd*sizeof(double));

	// Place the n_poly polymers at random positions in the
	// system. Give them a random velocity based on the value
	// of kBT that is specified.
	vxt = 0.0;
	vyt = 0.0;
	vzt = 0.0;
	for (i=0; i<n_dpd; i++) {
		// vx^2 + vy^2a + vz^2 = kBT, so... pick a random direction and
		// give the particle a velocity with a magnitude of T. Theta and
		// phi are two random angles in spherical coordinates (i.e., the 		
		// azimuthal and polar angles).
		sc_theta = ran3(&iseed) * PI;
		sc_phi   = ran3(&iseed) * 2.0 * PI;
		vx[i]    = sqrt(3.0*T) * cos(sc_phi) * sin(sc_theta);
		vy[i]    = sqrt(3.0*T) * sin(sc_phi) * sin(sc_theta);
		vz[i]    = sqrt(3.0*T) * cos(sc_theta);
		vxt     += vx[i];
		vyt     += vy[i];
		vzt     += vz[i];
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
	nx     = (int) (Lx/CELL_SIZE);
	ny     = (int) (Ly/CELL_SIZE);
	nz     = (int) (Lz/CELL_SIZE);
	n_cell = nx*ny*nz;

	cell_head = (int*) malloc(n_cell*sizeof(int));
	cell_list = (int*) malloc(n_dpd*sizeof(int));

	return;
}

//
// Gets number of lines in file.
//
int line_count(FILE *fp) {
	int i;
	int lines = 0;
	char buf[BUF_SIZE];

	for(;;) {
		size_t res = fread(buf, 1, BUF_SIZE, fp);

        	for(i = 0; i < res; i++) {
            		if (buf[i] == '\n') lines++;
		}

        	if (feof(fp)) break;
    	}
	
	// Read in architecture description.
	fseek(fp, 0, SEEK_SET);

	return lines;
}


//
// Create polymers using architecture defined in filename.
//
int make_arch(int comp, int idx, double phi_comp, char *filename) {
	FILE *arch_in;
        char buf[BUF_SIZE];
	char type;
	int *template_coord, *template_ptype;
	int  bond_id, block_length, block_type;
       	int i,j,k;
	int templ_idx;
        int arch_lines;
	double dx, dy, dz;
	double sc_theta, sc_phi;

	// Open our architecture file and make sure it's there.
	arch_in = fopen(filename, "r");
	if (!arch_in) {
		printf("ERROR: Cannot find architecture file %s.\n\n", filename);
		exit(-1);
	}
	arch_lines = line_count(arch_in) - 5;

	// Discard the comment line.
	fgets(buf, BUF_SIZE, arch_in);

	// Get the name of this architecture for output to disk.
	fgets(arch_names[comp], BUF_SIZE, arch_in);
	arch_names[comp][strcspn(arch_names[comp], "\n")] = 0;

	// Read in the degree of polymerization and
	// bonding potential parameters.
	fscanf(arch_in, "%d\n", &N[comp]);
	fscanf(arch_in, "%lf\n", &kH[comp]);
	fscanf(arch_in, "%lf %lf\n", &b0[comp], &b1[comp]);

	// Read in the molecule's template.
	template_coord = (int*) malloc(N[comp]*sizeof(int));
	template_ptype = (int*) malloc(N[comp]*sizeof(int));
	printf("Creating component %d (%s - %s): N = %d, kH = %3.2lf, b0 = %1.2lf, b1 = %1.2lf\n", comp,arch_names[comp],filename,N[comp],kH[comp],b0[comp],b1[comp]);
	templ_idx = 0;
	for (i=0; i<arch_lines; i++) {
		fscanf(arch_in, "%c %d %d %d\n", &type, &bond_id, &block_length, &block_type);
		switch (type) {
			// Place junction randomly (r) in the system.
			case 'r':
				template_coord[templ_idx] = -1;
				template_ptype[templ_idx] = block_type;
				templ_idx++;
				break;
			// Create a block of monomers.
			case 'b':
				template_coord[templ_idx] = bond_id;
				template_ptype[templ_idx] = block_type;
				templ_idx++;
				for (j=1; j<block_length; j++) {
					template_coord[templ_idx] = templ_idx - 1;
					template_ptype[templ_idx] = block_type;
					templ_idx++;
				}
				break;
			default:
				printf("Invalid architecture in %s (type = %c).\n\n", filename,type);
				fclose(arch_in);
				exit(-1);
				break;

		}
	}
	fclose(arch_in);

	/*
	for (i=0; i<N[comp]; i++) {
		printf("%d %d\n", template_coord[i], template_ptype[i]);
	}
	exit;
	*/

	// Place the polymers throughout the system using the template.
	n_poly[comp] = phi_comp * n_dpd / N[comp];
	printf("\t\t+ Adding %d molecules ('%s') to the system.\n", n_poly[comp], arch_names[comp]);
	
	k = 0;
	for (i=0; i<comp; i++) {
		k += n_poly[i]*N[i];
	}

	for (i=0; i<n_poly[comp]; i++) {
		for (j=0; j<N[comp]; j++) {
			templ_idx = j - template_coord[j];
			if (template_coord[j] < 0) {
				// -1 => Random location in system.
				x[k]  = ran3(&iseed) * Lx;
				y[k]  = ran3(&iseed) * Ly;
				z[k]  = ran3(&iseed) * Lz;
			
				// -1 => No intrachain bond processing.
				bond_list[k] = -1;
			} else {
				// Generate a random displacement.
				sc_theta = ran3(&iseed) * PI;
				sc_phi   = 2.0 * ran3(&iseed) * PI;
				dx = b0[comp] * cos(sc_phi) * sin(sc_theta);
				dy = b0[comp] * sin(sc_phi) * sin(sc_theta);
				dz = b0[comp] * cos(sc_theta);

				// Add the coordinates for this particle.
				x[k] = x[k-templ_idx] + dx;
				y[k] = y[k-templ_idx] + dy;
				z[k] = z[k-templ_idx] + dz;
		
				// Add the bond partner for intrachain interactions.
				bond_list[k] = k - templ_idx;
			}

			// Set the particle type
			p_type[k] = template_ptype[j];

			// Set the particle component.
			p_comp[k] = comp;			

			k++;
		}
	}

	printf("\n");
	free(template_coord);
	free(template_ptype);
	return (n_poly[comp]*N[comp]);
}

// Here we apply the periodic boundary condition (PBC) and place all particles
// in the system into our linked list so that we can quickly determine all interacting
// particle pairs within a distance of r_{c}.
//
// For efficiency, we also process intrachain interactions within this function since we
// are already looping over all DPD particles.
//
void make_list(void) {
	int i, cell_id, tmp_id;

	// Clear out the cell list. -1 means empty cell.
	memset(cell_head, -1, n_cell*sizeof(int));
	memset(cell_list, -1, n_dpd*sizeof(int));

	// Initialize the forces.
	memset(fx, 0.0, n_dpd*sizeof(double));
	memset(fy, 0.0, n_dpd*sizeof(double));
	memset(fz, 0.0, n_dpd*sizeof(double));

	// Reset energy.
	U     = 0.0;
	K     = 0.0;
	
	for (i=0; i<MAX_COMPONENTS; i++) {
		b_avg[i] = 0.0;	
	}

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

		// Consider any intrachain interactions that may be present.
		calc_intrachain(p_comp[i],i);
	}

	return;
}

//
// read_interactions - 
//
void read_interactions(char *filename) {
	FILE *inter_in;
	int inter_size;
	int *inter_idx;
	double *inter_val;
	char junk[BUF_SIZE];
	char buf[BUF_SIZE];
	int i,j,max_idx;

	// Open our architecture file and make sure it's there.
	inter_in = fopen(filename, "r");
	if (!inter_in) {
		printf("ERROR: Cannot find interaction strengths file %s.\n\n", filename);
		exit(-1);
	}
	
	// Determine the number of lines in the file.
	inter_size = line_count(inter_in) - 2;


	// Discard the comment line.
	fgets(buf, BUF_SIZE, inter_in);

	// Read the name of these interactions.
	fgets(buf, BUF_SIZE, inter_in);
	buf[strcspn(buf, "\n")] = 0;

	// Output information for user:
	printf("Reading interaction configuration from %s (%s).\n\n", filename, buf);

	// Read interactions.
	max_idx   = 0;
	inter_idx = (int*)    malloc(inter_size*sizeof(int));
	inter_val = (double*) malloc(inter_size*sizeof(double));
	for (i=0; i<inter_size; i++) {
		fscanf(inter_in, "%d %lf\n", &inter_idx[i], &inter_val[i]);

		if (inter_idx[i] > max_idx) {
			max_idx = inter_idx[i];
		}
	}
	max_idx++;

	a = (double*) malloc((max_idx) * sizeof(double));
	memset(a, 0.0, max_idx*sizeof(double));
	printf("Interaction table:\n");
	printf("---------------------------------------\n");
	for (i=0; i<inter_size; i++) {
		a[inter_idx[i]] = inter_val[i];
		printf("IDX: %d\t\t a: %3.2lf\n", inter_idx[i], a[inter_idx[i]]);
	}
	printf("\n");
	printf("Thermal noise (sigma): %2.2lf\n", sigma);
	printf("Dissipation (gamma)  : %2.2lf\n", pow(sigma, 2)/(2*T));
	printf("kBT                  : %2.2lf\n", T);
	printf("---------------------------------------\n\n");


	printf("\n");
	free(inter_idx);
	free(inter_val);
	return;
}

// Here we read in the parameters for our simulation from the config file we specify in filename.
//
void read_parameters(char *filename) {
	FILE *param_file;
	int i;
	char junk[BUF_SIZE]; 

	param_file = fopen(filename, "r");
	if (param_file == NULL) {
		printf("Parameter file %s not found.\n", filename);
		exit(1);
	}
	
	// Read in the "junk" comment from the file.
	fgets(junk, BUF_SIZE, param_file);

	// Read in the value that follows.
	fscanf(param_file, "%lf\n", &Lx);

	// Repeat the same process for the remaining lines of the file.
	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%lf\n", &Ly);
	
	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%lf\n", &Lz);

	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%lf\n", &rho);

	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%lf\n", &T);

	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%lf\n", &sigma);

	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%lf\n", &dt);

	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%d\n", &time_equil);

	fgets(junk, BUF_SIZE, param_file);
	fscanf(param_file, "%d\n", &time_prod);

	printf("Reading configuration of %d components.\n", MAX_COMPONENTS);
	for (i=0; i<MAX_COMPONENTS; i++) {
		// Read the volume fraction of this component.
		fgets(junk, BUF_SIZE, param_file);
		fscanf(param_file, "%lf\n", &phi[i]);

		// Read the architecture definition of this component, and strip newline
		// characters.
		fgets(junk, BUF_SIZE, param_file);
		fgets(arch_files[i], BUF_SIZE, param_file);
		arch_files[i][strcspn(arch_files[i], "\n")] = 0;

		printf("\t\t - Component %d: %s (phi = %lf)\n", i, arch_files[i], phi[i]);
	}
	printf("\n");

	// Read in the interaction configuration filename.
	fgets(junk, BUF_SIZE, param_file);
	fgets(junk, BUF_SIZE, param_file);
	junk[strcspn(junk, "\n")] = 0;

	// Read and assign the interactions strengths.	
	read_interactions(junk);

	// Define gamma from sigma, since they are related through the
	// fluctuation-dissipation theorem:
	gamma1 = pow(sigma, 2.0) / (2.0 * T);
	sigma  = sigma * sqrt(3.0 / dt);
	return;
}

// Below is the implementation of the velocity-Verlet
// algorithm as detailed in Groot & Warren.
//
void velocity_verlet(int prod, double lambda) {
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
	calc_forces(prod);

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
void write_xyz (int comp, char* filename, int append) {
	FILE *output;
	int i, offset;
	char p_char;

	printf("\nWriting coordinates to %s\n", filename);
	if (append) {
		output = fopen(filename, "aw");
	} else {
		output = fopen(filename, "w");
	}
	
	// Determine offset in arrays.
	offset = 0;
	for (i=0; i<comp; i++) {
		offset += n_poly[i] * N[i];
	}

	fprintf(output, "%d\n", n_poly[comp]*N[comp]);
	fprintf(output, "# Snapshot from time t = %d\n", t);
	for (i=offset; i<offset+n_poly[comp]*N[comp]; i++) {
		p_char = 'x';
		switch(p_type[i]) {
			case 1:
				p_char = 'C';
				break;
			case 3:
				p_char = 'S';
				break;
			default:
				p_char = 'H';
				break;
		}	
		fprintf(output, "%c %3.6lf %3.6lf %3.6lf\n", p_char, x[i], y[i], z[i]);
	}
	fclose(output);
	return;
} 
