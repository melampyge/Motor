
""" Plot the overlap function per simulation"""

### example command line arguments: 
###    -e=1.0 -f=0.5 -k=10.0
###    -fl=/local/duman/SIMULATIONS/Cells_in_LAMMPS/density_0.8/ 
###         -sb=/usr/users/iff_th2/duman/Cells_in_LAMMPS/PLOTS/
###             -sf=Overlap_fnc

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import read_write
import misc_tools 
   
##############################################################################
     
def plot_data(x, y, sim, savebase, savefolder, save_eps):
    """ plot the data as a function of the chosen parameter"""

    ### set general plot properties

    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)
    #downlim = -1
    #uplim = sim.lx/4.
    num_ticks = 5
    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
                          
    ### plot 

    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    label = '$\\phi=$' + str(sim.density) + '$,km=$' + str(sim.km)
    line0 = ax0.plot(x/sim.tau_D, y, \
                     linewidth=2.0, label=label)
    
    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels
        
    ax0.set_xlabel("$\\Delta t/\\tau_{D}$", fontsize=40)
    ax0.set_ylabel("$Q(\\Delta t)$", fontsize=40)

    ### limits

    #ax0.set_xlim((-1, 15))
    #ax0.set_ylim((downlim, uplim))
    
    ### ticks
    
    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

#    ax0.legend(bbox_to_anchor=(1.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
#        prop={'size': 20}, mode="expand", frameon=False)
    
    ### save

    savepath1 = base + savefolder + "_density_" + str(sim.density) + "_kappa_" + str(sim.kappa) + \
        "_km_" + str(sim.km) + '_panti_' + str(sim.pa) + ".png"
    if save_eps:
        savepath2 = base + savefolder + "_density_" + str(sim.density) + "_kappa_" + str(sim.kappa) + \
        "_km_" + str(sim.km) + '_panti_' + str(sim.pa) + ".eps"            
    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    if save_eps:
        plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                       
        
    return
        
##############################################################################

def main():
    
    ### get the data folder
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--density", type=float, help="Density")
    parser.add_argument("-ka", "--kappa", type=float, help="Bending rigidity")
    parser.add_argument("-km", "--km", type=float, help="Motor spring strength")
    parser.add_argument("-p", "--pa", type=float, help="Probability when filaments are antialigned")    
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/Motor/Data/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/Motor/Data/")    
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/Motor/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/Motor/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="Overlap_fnc", \
                        help="Specific folder for saving, as in Overlap_fnc") 
    parser.add_argument("-s","--save_eps", action="store_true", help="Decide whether to save in eps or not")            
    args = parser.parse_args()
    
    ### read general information from the folder
    
    path1 = 'density_' + str(args.density) + '/kappa_' + str(args.kappa) + \
        '/km_' + str(args.km) + '/panti_' + str(args.pa) + '/' 
    datafolder = args.folder + path1
    sim = read_write.read_sim_info(datafolder)
    
    ### read the data
    
    path2 = '_density_' + str(args.density) + '_kappa_' + str(args.kappa) + \
        '_km_' + str(args.km) + '_panti_' + str(args.pa) 
    analysisdatabase = '/usr/users/iff_th2/duman/Motor/DATA/'
    analysisdatabase += args.savefolder + '/'
    analysisdata = analysisdatabase + args.savefolder + path2 + '.txt'
    x, y = read_write.read_2d_analysis_data(analysisdata)
        
    ### plot the data in the given time window
    
    plot_data(x, y, sim, args.savebase, args.savefolder, args.save_eps)
    
    return
    
##############################################################################

if __name__ == '__main__':
    main()    
    
##############################################################################
