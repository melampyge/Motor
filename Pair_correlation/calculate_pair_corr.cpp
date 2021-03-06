
/* calculate intermediate scattering function */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/calculate_pair_corr.cpp ${upath}/read_write.cpp -lhdf5 -fopenmp -o paircorr
// ./paircorr out.h5 ${path}/Pair_corr.txt

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "omp.h"
#include "../Utility/read_write.hpp"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void compute_pair_corr (double **x, double **y, double *gr, int max_distance, int nsteps, int ncells, double lx, double ly) {
  /* calculate pair correlation function */

  int *cnt_distances = new int[max_distance];
  for (int j = 0; j < max_distance; j++) {
    cnt_distances[j] = 0;
  }
  
  for (int step = 0; step < nsteps; step++) {
    
    cout << "step / nsteps: " << step << " / " << nsteps << endl;
    
    for (int i = 0; i < ncells-1; i++) {
      
      for (int j = i+1; j < ncells; j++) {
	
	// calculate the distance between the particles
	
	double dx = x[step][j] - x[step][i];
	dx = min_img_dist(dx, lx);
	double dy = y[step][j] - y[step][i];
	dy = min_img_dist(dy, ly);
	double dr = sqrt(dx*dx + dy*dy);
	
	// increase the histogram depending on the distance-based bin
	
	int bin = inearbyint(dr);
	if (bin < max_distance) {
	  cnt_distances[bin] += 2;
	}
	
      } // particle j
      
    } //particle i
    
  } // timesteps
  
  // normalization
  
  double area = lx*ly;
  double number_density = ncells/area;
  double jfactor = 1./(2*pi*ncells*number_density);
  double delr = 1.;
  for (int j = 1; j < max_distance; j++) {
    double factor = jfactor/(delr*j);
    cnt_distances[j] /= nsteps;
    gr[j] = factor*cnt_distances[j];
  }

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

  // get the file name by parsing
  
  char *filename = argv[1];
  cout << "Calculating pair correlation function of the following file: \n" << filename << endl;

  // read in general simulation data
  
  int nsteps, nbeads, nsamp, ncells, nbpf;
  nsteps = nbeads = nsamp = ncells = nbpf = 0;
  double lx, ly, dt, density, kappa, km, pa, pp, bl, sigma;
  lx = ly = dt = density = kappa = km = pa = pp = bl = sigma = 0.;
  
  read_sim_data(filename, nsteps, nbeads, nsamp, ncells, nbpf, lx, ly, dt, density, kappa, km, pa, pp, bl, sigma);
  
  // print simulation information
  
  cout << "nsteps = " << nsteps << endl;
  cout << "ncells = " << ncells << endl;
  
  // read in array data
  
  int nbpc[ncells];
  for (int i = 0; i < ncells; i++) nbpc[i] = nbpf;
  //read_integer_array(filename, "/cell/nbpc", nbpc);
 
  // read in the cell position data all at once
  
  /* the position data is stored in the following format:
  (nsteps, 2, ncells)
  the data will be loaded as follows:
  (nsteps, ncells) in x and y separately 
  */
  
  double **x = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) x[i] = new double[ncells];
  
  double **y = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) y[i] = new double[ncells];
  
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < ncells; j++) {
      x[i][j] = 0.; y[i][j] = 0.;
    }
  }
  
  read_all_pos_data(filename, x, y, nsteps, ncells, "/cells/comu");
  
  // allocate the arrays
  
  int max_distance = static_cast<int>(lx/2.);
  double *gr = new double[max_distance];
  for (int j = 0; j < max_distance; j++) {
    gr[j] = 0;
  }
  
  // calculate the pair correlation function
  
  compute_pair_corr(x, y, gr, max_distance, nsteps, ncells, lx, ly);  
  
  // write the computed data
  
  char *outfilepath = argv[2];
  cout << "Writing pair correlation function to the following file: \n" << outfilepath << endl;  
  write_1d_analysis_data(gr, max_distance, outfilepath);
  
  // deallocate the arrays
  
  for (int i = 0; i < nsteps; i++) {
    delete [] x[i];
    delete [] y[i];
  }
  delete [] x;
  delete [] y;
  
  return 0;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
