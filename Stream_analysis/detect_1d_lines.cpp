
/* detect 1 dimensional lines for motors */

// COMPILATION AND RUN COMMANDS:
// g++ -Wl,-rpath=$HOME/hdf5/lib -L$HOME/hdf5/lib -I$HOME/hdf5/include ${spath}/detect_1d_lines.cpp ${upath}/read_write.cpp -lhdf5 -o detect_lines
// ./detect_lines

//////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include "omp.h"
#include "../Utility/read_write.hpp"
#include "../Utility/kdtree.hpp"
#include "../Utility/cartesian.hpp"

#define pi M_PI

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void populate_extended_pos_array(vector<Point<2> > &points, const int npoints, const int ncells, const int step, double **x, double **y, double lx, double ly) {
  /* populate an extended array with all the image points */
  
  for (int i = 0; i < ncells; i++) {
    points[i].x[0] = x[step][i];
    points[i].x[1] = y[step][i];
    
    points[i+ncells].x[0] = x[step][i] - lx;
    points[i+ncells].x[1] = y[step][i] + ly;
    
    points[i+2*ncells].x[0] = x[step][i];
    points[i+2*ncells].x[1] = y[step][i] + ly;
    
    points[i+3*ncells].x[0] = x[step][i] + lx;
    points[i+3*ncells].x[1] = y[step][i] + ly;
    
    points[i+4*ncells].x[0] = x[step][i] + lx;
    points[i+4*ncells].x[1] = y[step][i];
    
    points[i+5*ncells].x[0] = x[step][i] + lx;
    points[i+5*ncells].x[1] = y[step][i] - ly;
    
    points[i+6*ncells].x[0] = x[step][i];
    points[i+6*ncells].x[1] = y[step][i] - ly;
    
    points[i+7*ncells].x[0] = x[step][i] - lx;
    points[i+7*ncells].x[1] = y[step][i] - ly;
    
    points[i+8*ncells].x[0] = x[step][i] - lx;
    points[i+8*ncells].x[1] = y[step][i];
  }
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {
  
  // get the file name by parsing
  
  char *filename = argv[1];
  cout << "Detecting 1D lines of the following file: \n" << filename << endl;
  
  // read in general simulation data
  
  int nsteps, nbeads, nsamp, nfils, nbpf;
  nsteps = nbeads = nsamp = nfils = nbpf = 0;
  double lx, ly, dt, density, kappa, km, pa, bl, sigma;
  lx = ly = dt = density = kappa = km = pa = bl = sigma = 0.;
  
  read_sim_data(filename, nsteps, nbeads, nsamp, nfils, nbpf,
                lx, ly, dt, density, kappa, km, pa, bl, sigma);
  
  // print simulation information
  
  cout << "nsteps = " << nsteps << endl;
  cout << "nfils = " << nfils << endl;
  
  // read in the filament position data all at once
  
  /* the position data is stored in the following format:
   (nsteps, 2, nfils)
   the data will be loaded as follows:
   (nsteps, nfils) in x and y separately
   */
  
  double **x = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) x[i] = new double[nfils];
  
  double **y = new double*[nsteps];
  for (int i = 0; i < nsteps; i++) y[i] = new double[nfils];
  
  for (int i = 0; i < nsteps; i++) {
    for (int j = 0; j < nfils; j++) {
      x[i][j] = 0.; y[i][j] = 0.;
    }
  }
  
  read_all_pos_data(filename, x, y, nsteps, nfils, "/positions/xi");
  
  // get the image positions in the central box for all time steps at once
  
  get_img_pos(x, y, nsteps, ncells, lx, ly);

  // write the computed data
  
//  char *outfilepath = argv[2];
//  cout << "Writing detected 1D lines to the following file: \n" << outfilepath << endl;
//  write_single_analysis_data(order_param, outfilepath);
  
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
