
#pragma once

#include <string>
#include <fstream>
#include "/usr/users/iff_th2/duman/hdf5/include/hdf5.h"
//#include "/homec/jiff26/jiff2610/hdf/include/hdf5.h"

/* wrapper to read integer data from hdf5 file --note that buffer needs to be an array of size 1 for single entries-- */
int read_integer_data (hid_t file, std::string path_in_file, int *buffer);

/* wrapper to read double data from hdf5 file --note that buffer needs to be an array of size 1 for single entries-- */
double read_double_data (hid_t file, std::string path_in_file, double *buffer);
  
/* wrapper to read integer array data from hdf5 file --note that buffer needs to be the array size-- */
void read_integer_array (std::string filename, std::string path_in_file, int *buffer);

/* wrapper to read double array data from hdf5 file --note that buffer needs to be the array size-- */
void read_double_array (std::string filename, std::string path_in_file, double *buffer);
  
/* read general simulation data in hdf5 format */
void read_sim_data (std::string filename, int &nsteps, int &nbeads, int &nsamp, int &nfils, int &nbpf, double &lx, double &ly, double &dt, double &density, double &kappa, double &km, double &pa, double &bl, double &sigma);

/* read the position data at a single timestep */  
void read_single_pos_data (int step, hid_t dataset, hid_t dataspace, double **x, double **y, int natoms);
  
/* read position data in hdf5 format all at once */
void read_all_pos_data (std::string filename, double **x, double **y, int nsteps, int natoms, std::string datapath);

/* write analysis data with single parameter to the outfile */
void write_single_analysis_data (double data, std::string outpath);
  
/* write the 1d analysis data to the outfile */
void write_1d_analysis_data (double *x, int ndata, std::string outpath);
  
/* write the 2d analysis data to the outfile */
void write_2d_analysis_data (double *x, double *y, double ndata, std::string outpath);
    
