
""" Read hdf5 data of cells and write some data
        IMPORTANT NOTE:
        All trajectory information is assumed UNWRAPPED!"""

### example command line arguments: 
###

##############################################################################

import argparse
import numpy as np
import os
import h5py
import misc_tools

##############################################################################
    
def read_h5_file(folder):
    """ read data from hdf5 file"""
    
    ### file path
    
    fpath = folder + 'out.h5'
    assert os.path.exists(fpath), "The out.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### positions of beads
    
    xu = np.array(fl['/beads/xu'], dtype=np.float32)
    cid = np.array(fl['/beads/cid'], dtype=np.float32)
    
    ### cell information
    
    comu = np.array(fl['/cells/comu'], dtype=np.float32)
    pol = np.array(fl['/cells/pol'], dtype=np.float32)
    nbpc = np.array(fl['/cells/nbpc'], dtype=np.float32)
    
    ### simulation information
    
    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    ncells = fl['/info/ncells'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]
    
    ### simulation parameters
    
    eps = fl['/param/eps'][...]
    rho = fl['/param/rho'][...]
    fp = fl['/param/fp'][...]
    areak = fl['/param/areak'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]
    
    fl.close()
    
    ### generate classes to submerge data
    
    sim = misc_tools.Simulation(lx, ly, dt, nsteps, ncells, nbeads, nsamp, nbpc, \
                     eps, rho, fp, areak, bl, sigma)
    cells = misc_tools.Cells(comu, pol, nbpc, sim)
    beads = misc_tools.Beads(xu, cid)
    
    return sim, cells, beads

##############################################################################
    
def read_h5_file_fils(folder):
    """ read data from hdf5 file for filament information"""
    
    ### file path
    
    fpath = folder + 'out_fil.h5'
    assert os.path.exists(fpath), "The out_fil.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')
    
    ### cell information
    
    xu = np.array(fl['/positions/xi'], dtype=np.float32)
    #xi = np.array(fl['/positions/x'], dtype=np.float32)
    #xi = [xt.T  for xt in xi[:]]
    pol = np.array(fl['/positions/ori'], dtype=np.float32)
    pol = np.array([xt.T for xt in pol[:]])
    
    ### simulation information
    
    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nfils = fl['/info/nfils'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]
    nbpf = fl['/info/nbpf'][...]
    
    ### simulation parameters
    
    density = fl['/param/density'][...]
    kappa = fl['/param/kappa'][...]
    km = fl['/param/km'][...]
    pa = fl['/param/pa'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]
    
    fl.close()
    
    ### generate classes to submerge data
    
    sim = misc_tools.Simulation(lx, ly, dt, nsteps, nfils, nbeads, nsamp, nbpf, \
                     density, kappa, km, pa, bl, sigma)
    fils = misc_tools.Cells(xu, pol, nbpf, sim)
    
    return sim, fils
    
##############################################################################
    
def read_sim_info(folder):
    """ read simulation info from hdf5 file"""
    
    ### file path
    
    fpath = folder + 'out_fil.h5'
    assert os.path.exists(fpath), "The out_fil.h5 file does NOT exist for " + fpath
    fl = h5py.File(fpath, 'r')    
    
    ### simulation information
    
    lx = fl['/info/box/x'][...]
    ly = fl['/info/box/y'][...]
    dt = fl['/info/dt'][...]
    nsteps = fl['/info/nsteps'][...]
    nfils = fl['/info/nfils'][...]
    nbeads = fl['/info/nbeads'][...]
    nsamp = fl['/info/nsamp'][...]
    nbpf = fl['/info/nbpf'][...]
    
    ### simulation parameters
    
    density = fl['/param/density'][...]
    kappa = fl['/param/kappa'][...]
    km = fl['/param/km'][...]
    pa = fl['/param/pa'][...]
    bl = fl['/param/bl'][...]
    sigma = fl['/param/sigma'][...]
    
    fl.close()
    
    ### generate classes to submerge data
    
    sim = misc_tools.Simulation(lx, ly, dt, nsteps, nfils, nbeads, nsamp, nbpf, \
                     density, kappa, km, pa, bl, sigma)
    
    return sim
    
##############################################################################
        
def write_2d_analysis_data(x, y, savebase, savefolder, sim):
    """ write 2d analysis data to the corresponding file
    EXAMPLE:
    savebase: /usr/users/iff_th2/duman/Motor/DATA/
    savefolder: MSD"""
    
    ### create the path
    
    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)
    fpath = base + savefolder + "_density_" + str(sim.density) + "_kappa_" + \
        str(sim.kappa) + "_km_" + str(sim.km) + "_panti_" + str(sim.pa) + ".txt"  
  
    ### write the data
    
    fl = open(fpath, 'w')
    
    N = len(x)
    for j in range(N):
        fl.write(str(x[j]) + '\t\t' + str(y[j]) + '\n')
    
    fl.close()

    return
 
##############################################################################

def read_single_analysis_data(f):
    """ read single analysis data"""
    
    data = np.loadtxt(f, dtype=np.float64)

    return data  

##############################################################################    
    
def read_2d_analysis_data(f):
    """ read 2d analysis data"""
    
    data = np.transpose(np.loadtxt(f, dtype=np.float64))
    x = data[0]
    y = data[1]

    return x, y   

##############################################################################

def gen_folders(rho, kappa, km, pa, analysis, dbase, analysisdbase):
    """ data and analysis folders generator"""
    
    path1 = 'density_' + + str(rho) + "_kappa_" + \
        str(kappa) + "_km_" + str(km) + "_panti_" + str(pa)
    path2 = analysis + '_density_' + + str(rho) + "_kappa_" + \
        str(kappa) + "_km_" + str(km) + "_panti_" + str(pa) + '.txt'  
    datafolder = dbase + path1 + '/'
    analysisfile = analysisdbase + path2  

    return datafolder, analysisfile
    
##############################################################################
    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")        
    args = parser.parse_args()    
    read_h5_file(args.folder)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
