
""" Combine multiple dump files from different restart instances into a single hdf5 file"""

### example command line arguments: 
###    -fl=/homea/ias2/duman/Cells_in_LAMMPS/ -t -d=0.8 -e=1.0 -f=0.5 -a=10.0 
###         -dt=0.001 -ns=50000 -b=0.5 -s=1.0 -nc=5000

##############################################################################

import argparse
import numpy as np
import os
import h5py

##############################################################################

def read_contextual_info():
    """ read the contextual information provided by the user"""
    
    ### get the data folder and the last timestep info
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", help="Folder containing data")
    parser.add_argument("-nfl", "--newfolder", help="Folder to put the data in")    
    parser.add_argument("-t", "--tstep", nargs="?", const="100000000", \
                            type=int, help="The last time step that is being searched for")
    parser.add_argument("-d", "--density", type=float, help="Packing fraction of the system")  
    parser.add_argument("-k", "--kappa", type=float, help="Bending rigidity")
    parser.add_argument("-m", "--km", type=float, help="Motor strength")
    parser.add_argument("-pa", "--panti", type=float, 
                        help="Probability of motor attachment in anti-orientation")
    parser.add_argument("-pp", "--ppar", type=float, 
                        help="Probability of motor attachment in parallel-orientation")    
    parser.add_argument("-dt", "--timestep", type=float, help="Timestep of the simulation")
    parser.add_argument("-ns", "--nsamp", type=int, help="Sampling rate of data")
    parser.add_argument("-b", "--bl", type=float, help="Bond length of the simulation")    
    parser.add_argument("-s", "--sigma", type=float, help="Lennard Jones length")     
    parser.add_argument("-nc", "--ncells", type=int, help="Number of cells")    
    parser.add_argument("-nb", "--nbpc", type=int, help="Number of beads per bond")              
    args = parser.parse_args()
    
    ### generate folder path
    
    folder = args.folder + 'density_' + str(args.density) + '/pa_' + str(args.panti) + \
        '/pp_' + str(args.ppar) 
    newfolder = args.newfolder + 'density_' + str(args.density) + '/pa_' + str(args.panti) + \
        '/pp_' + str(args.ppar)
    fpath = folder + '/out.dump'
    assert os.path.exists(fpath), "\nOUT.DUMP DOES NOT EXIST FOR: " + folder 

    print fpath

    ### determine the total number of beads and the box size
    
    fl = open(fpath, 'r')
    fl.readline()
    fl.readline()
    fl.readline()
    line = fl.readline()
    line = line.split()
    nbeads = int(line[0])
    
    fl.readline()
    line = fl.readline()
    line = line.split()
    lx = float(line[1])
    line = fl.readline()
    line = line.split()
    ly = float(line[1])

    fl.close()
    
    ### total number of steps
    
    nsteps = args.tstep/args.nsamp
    nsteps += 1
            
    return folder, newfolder, nbeads, nsteps, lx, ly, args

##############################################################################
    
def get_number_of_snaps(f, nbeads):
    """ determine the number of snapshots in the file"""
    
    os.system('wc -l ' + f + ' > tmp.txt')
    ifile = open('tmp.txt')
    line = ifile.readline()
    line = line.split()
    nlines = int(line[0])
    nsnaps = nlines/(nbeads+9)
    ifile.close()
    os.system('rm tmp.txt')
    
    return nsnaps

##############################################################################
    
def read_pos(fl, x, nbeads, nsnaps, lx, ly, checked, tstep_cnt, T, mid):
    """ read the position data from the file and return the last timestep at the end"""

    ### read the positions unique per each tstep in a single dump file

    already_checked = False
    for snap in range(nsnaps):
        
        ### read the headers to check the uniqueness of current tstep
        
        fl.readline()
        line = fl.readline()
        line = line.split()
        tstep = int(line[0])
        
        ### finish if the last tstep is exceeded already
        
        if tstep > T:
            return tstep_cnt, tstep
        
        ### make sure the current tstep is unique
        
        if tstep not in checked:
            checked.append(tstep)
            tstep_cnt += 1
            already_checked = False
        else:
            print "ALREADY CHECKED " + str(tstep)             
            already_checked = True
    
        ### read the remaining part of the header part
        
        for j in range(7):
            fl.readline()
         
        ### read the positions per bead if the tstep is unique
        
        for j in range(nbeads):
            line = fl.readline()
            if already_checked:
                continue
            line = line.split()
            if len(line) < 2:
                print line
                continue
            bid = int(line[0]) - 1
            mid[bid] = int(line[1]) - 1
            ix = float(line[6])
            iy = float(line[7])
            x[tstep_cnt, 0, bid] = float(line[3]) + ix*lx
            x[tstep_cnt, 1, bid] = float(line[4]) + iy*ly

        print tstep_cnt, tstep

    return tstep_cnt, tstep

##############################################################################
    
def read_pos_from_dump_files(folder, nbeads, ncells, nsteps, T, lx, ly):
    """ read the position data of each dump file until the last tstep is reached"""

    ### generate file path and the total number of snapshots in the file
    
    x = np.zeros((nsteps, 2, nbeads), dtype=np.float32)
    mid = np.zeros((nbeads), dtype=np.int32)
    tstep_cnt = -1
    checked = []
    tstep = 0
    fpath = folder + '/out.dump'
    assert os.path.exists(fpath), "out dump file does NOT exist for: " + fpath
    fl = open(fpath, 'r')
    tstep_cnt, tstep = read_pos(fl, x, nbeads, nsteps, lx, ly, checked, tstep_cnt, T, mid)
    fl.close()
        
    return x, mid
    
##############################################################################
    
def write_h5_file(folder, x, d, mid, com, nbeads, nsteps, nbpf, lx, ly, args):
    """ write data to hdf5 file"""
    
    ### file path
    
    os.system("mkdir -p " + folder)
    fpath = folder + '/out.h5'
    fl = h5py.File(fpath, 'w')
    
    ### positions of beads
    
    bead = fl.create_group('beads')
    bead.create_dataset('xu', (nsteps, 2, nbeads), data=x, dtype=np.float32, compression='gzip') 
    bead.create_dataset('cid', data=mid) 
    #bead.create_dataset('pol', (nsteps, args.nbeads-1), data=d, dtype=np.float32, compression='gzip')
    
    
    ### cell information
    
    cell = fl.create_group('cells')
    cell.create_dataset('comu', (nsteps, 2, args.ncells), data=com, dtype=np.float32, compression='gzip') 
    
    ### simulation information
    
    info = fl.create_group('info')
    box = info.create_group('box')
    box.create_dataset('x', data=lx)
    box.create_dataset('y', data=ly)
    info.create_dataset('dt', data=args.timestep)
    info.create_dataset('nsteps', data=nsteps)
    info.create_dataset('ncells', data=args.ncells)
    info.create_dataset('nbeads', data=nbeads)
    info.create_dataset('nsamp', data=args.nsamp)
    cell.create_dataset('nbpf', data=nbpf)
    
    ### simulation parameters
    
    param = fl.create_group('param')
    param.create_dataset('density', data=args.density)
    param.create_dataset('kappa', data=args.kappa)
    param.create_dataset('km', data=args.km)
    param.create_dataset('pa', data=args.pa)
    param.create_dataset('pp', data=args.pa)    
    param.create_dataset('bl', data=args.bl)
    param.create_dataset('sigma', data=args.sigma)
    
    fl.close()
    
    return

##############################################################################
    
def calculate_com_of_cells(xu, nsteps, nbpc, args):   
    """ calculate the center of mass of cells"""
    
    print "Calculating center of mass of filaments"
    
    com = np.zeros((nsteps, 2, args.ncells), dtype=np.float32)
    
    k = 0
    for j in range(args.ncells):
        com[:, :, j] = np.mean(xu[:, :, k:k+nbpc], axis=2)
        k += nbpc
    
    return com

##############################################################################
    
def calculate_com_of_cells_one_liner(xu, mid, nsteps, nbeads, nbpc):   
    """ calculate the center of mass of cells"""
    
    splitted = np.split(xu, np.cumsum(nbpc)[:-1], axis=2)
    r = np.array([np.mean(sfil, axis=2) for sfil in splitted])
    com = np.swapaxes(np.swapaxes(r, 0, 1), 1, 2)
    
    return com

##############################################################################
 
#def calculate_orientations(xu, nbeads, ncells, nsteps, nbpf, lx, ly):
#    """ calculate bond orientations"""
#    
#    print "Calculating bond orientations"
#    
#    d = np.zeros((nsteps, nbeads-1), dtype=np.float32)
#    k = 0
#    for tstep in range(nsteps):
#        for j in range(ncells):
#            for l in range(nbpf):
#                dx = x[tstep, k+1, ]
#                d[tstep, k] =
#                k += 1
#        
#    
#    return d
    
##############################################################################
    
def main():

    folder, newfolder, nbeads, nsteps, lx, ly, args = read_contextual_info()
    xu, mid = read_pos_from_dump_files(folder, nbeads, args.ncells, nsteps, args.tstep, lx, ly)
    #d = calculate_orientations(xu, nbeads, args.ncells, nsteps, args.nbpc, lx, ly)
    com = calculate_com_of_cells(xu, nsteps, args.nbpc, args)
    write_h5_file(newfolder, xu, mid, com, nbeads, nsteps, args.nbpc, lx, ly, args)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################