
""" A collection of helper functions and classes to do post-processing analysis
        IMPORTANT NOTE:
        All trajectory information is assumed UNWRAPPED!"""

##############################################################################

import numpy as np
import math

############################################################################

def nearbyint(x):
    """ Round to the nearby integer"""
    
    if x >= 0:
        return math.floor(x+0.5)
    else:
        return math.floor(x-0.5)

############################################################################
    
def min_img_dist(x1, x2, lx):
    """ compute the minimum image distance between two positions"""
    
    dx = x2 - x1 
    return dx-nearbyint(dx/lx)*lx

############################################################################

def get_img(x, lx):
    """ get the image position in the central box 
    -- can be numpy array or single pos"""
    
    return x-np.floor(x/lx)*lx

##############################################################################

class Beads:
    """ data structure for storing particle/bead information"""
    
    def __init__(self, x, cid, d):
        
        ### assign bead positions
        
        self.xu = x
        
        ### assign cell indices to beads
        
        self.cid = cid
        
        ### assign orientations to bonds
        
        #self.pol = d
        
        return
        
    ##############
    
    def get_img_pos(self, lx):
        """ get the image positions of beads in the central box"""
        
        self.xi = get_img(self.xu, lx)
        
        return

##############################################################################

class Cells:
    """ data structure for storing cell information"""
    
    def __init__(self, x, d, nbpf, sim):
        
        ### assign bead positions
        
        self.xu = x
        
        ### assign number of beads info
        
        self.nbpc = nbpf
        
        return

    ##############
    
    def get_img_pos(self, lx):
        """ get the image positions of cells in the central box"""
        
        self.xi = get_img(self.xu, lx)
        
        return
        
##############################################################################

class Simulation:
    """ data structure for storing general simulation information"""

    def __init__(self, lx, ly, dt, nsteps, nfils, nbeads, nsamp, nbpf, \
                 density, kappa, km, pa, pp, bl, sigma):
        
        self.lx = lx
        self.ly = ly
        self.dt = dt
        self.nsteps = nsteps
        self.nbeads = nbeads
        self.nsamp = nsamp
        self.nfils = nfils
        self.nbpc = nbpf
        self.density = float("{:.1f}".format(float(density)))
        self.kappa = float("{:.1f}".format(float(kappa)))
        self.km = float("{:.1f}".format(float(km)))
        self.pa = float("{:.1f}".format(float(pa)))
        self.pp = float("{:.1f}".format(float(pp)))
        self.bl = bl
        self.sigma = sigma
        
        ### normalize certain variables
          
        self.dt *= self.nsamp
        
        ### define more simulation parameters
        
        self.kT = 1
        self.gamma_n = 1
        self.length = self.bl*self.nbpf
        self.tau_D = self.length**2 * self.nbpf / 4. / self.kT
#        if self.fp == 0.:
#            self.tau_A = 0.0
#        else:
#            self.tau_A = 2 * self.r_avg * self.gamma_n / self.fp
        
        return
        
##############################################################################

class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
        return
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ### increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ### get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##############################################################################
        
