#!/usr/bin/env python

# -----------------------
# FDAP fitting script 1.0
# -----------------------
#
# REQUIREMENTS: Python 2.7, Numpy, Scipy, mpi4py, OpenMPI
#
# DESCRIPTION:
# 	The routine includes three main blocks. The first block is the
# 	numerical Laplace inversion of arbitrary functions (modified from J. Valsa
# 	and Brancik, 1998). The algorithm is relatively simple and produces accurate
# 	transforms of functions of different complexity. The second block of the
# 	routine preforms the least square analysis. For this to be done, a simple
# 	brute-force mapping procedure is used to approximately minimize the chi^2-
# 	value in the parametric space (k*_on, k_off). The minimization procedure
# 	consists in ordered discretization of the parametric space and subsequent
# 	finding a vector of parameters that minimizes the chi^2-value. Furthermore,
# 	mpi4py is used to reduce the total calculation time. The third block of the
# 	routine computes 95% confidence intervals for the estimated parameters (the
# 	parametric vector computed in the second block). All necessary formulae were
# 	taken from D.M. Bates and D.G. Watts, "Nonlinear Regression Analysis and Its
# 	Applications", 1988. A very similar algorithm is implemented in Matlab. For
# 	more details, please see the Supporting Material for Igaev et al., 2014
#
# HOW IT WORKS:
# 	1) FDAP curves (each having been normalized to the unity) must be saved as
# 	   separate *.dat files (e.g., curve_1.dat, curve_2.dat, etc.) and put together
# 	   with the script into a folder
#	2) Choose parameters you want to have in the section "Defining constants and ranges"
# 	   (diffusion constant Df, half-length of the activation area R, vIni_kon, vEnd_kon,
# 	   vIni_koff, and vEnd_koff). Note that setting, e.g., vIni_kon to -2 and vEnd_kon
# 	   to 2 will create a range of k*_on values like [0.01, 100] discretized with a
# 	   vStep = 0.01 increment of a log scale
# 	3) Run "./residual_analysis_kon_koff.py" to get a basic usage example
# 	4) For massive processing use a simple auxiliary script "FDAP_automatic_analysis.bash"
# 	5) The python script will produce *.txt files contaning fit parameters for each curve.
# 	   These files can be further processed either manually in Scidavis, Origin, etc. or
# 	   by applying additional bash/tcl/python scripts
#
# TODOs:
# 	- Confidence interval calculations for other kinetic models
# 	- Brute-force mapping requires a sufficient number of CPUs (optimal values range
# 	  between 8 and 30). A more elegant method (e.g., Gauss-Newton or any of quasi-Newton
# 	  methods) has to be implemented.
# 	- Implement multiple curve processing directly in the python script 
# 	
#
# AUTHORS: Maxim Igaev (maxim.igaev at biologie.uni-osnabrueck.de):
# 		- FDAP formulae
# 		- Brute-force mapping
# 		- Confidence intervals
#	   Frederik Suendermann (frederik.suendermann at biologie.uni-osnabrueck.de):
#		- Code parallelization
# 		- Inverse Laplace transform
# 		- Data output
#


import numpy as np
import scipy as sc
import math
import scipy.special as scs
import scipy.io as sio
import sys
import copy as cp
from numpy import linalg
from numpy.linalg import inv
from scipy.stats import t
from mpi4py import MPI
from optparse import OptionParser


#---------------------------------------------------------------------
# Multiplication of matrices of different sizes
# Needed for confidence interval caclulation
#
def multmat(a,b):
    m=len(a)
    n=(a[0])
    k=len(b)
    res=[]
    if len(b) != 1:
        p=len(b[0])-1
    else:
        p=0

    if len(b)==0:
        print('bad size!')
    #elif n!=k:
    #    print('bad size!')
    else:
        n=k

        for q in range(m):
            res.append([0])
        for q in range(m):
            for w in range(p):
                res[q].append(0)

        for i in range(m):
            for j in range(p+1):
                for r in range(n):
                    res[i][j] = a[i][r]*b[r][j] + res[i][j]
    return res
#
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Inverse Laplace transform class
# Parameters a, ns, nd must not be changed
#
class TransformClass:
    def __init__(self,tini,tend,nnt,a=6,ns=20,nd=19):
        self.a = a
        self.ns = ns
        self.nd = nd
        self.tini = tini
        self.tend = tend
        self.nnt = nnt
        self.alfa = np.array([(np.float(a)+(n-1.0)*sc.pi*1j) for n in xrange(1,(ns+1+nd)+1)])
        self.beta = np.array([(-math.exp(a)*math.pow(-1.0,n)) for n in xrange(1,(ns+1+nd)+1)])
        self.stepSize = ((self.tend-self.tini)/self.nnt)
        self.radt = np.array([self.tini+(i*self.stepSize) for i in xrange(self.nnt+1)])

    def invlap(self,FF,kon,koff):
#        alfa = cp.deepcopy(self.alfa)
        beta = cp.deepcopy(self.beta)
        ft = np.zeros((self.nnt+1)) # return value
#        if tini==0:
#            print "tini appears to be 0!!!"
#            radt=radt[1:]
        N = np.array([i+1 for i in xrange(self.nd)])
        bdif = sc.cumsum(scs.gamma(self.nd+1)/scs.gamma(self.nd+2-N)/scs.gamma(N))
        bdif = bdif[::-1]
        bdif = bdif/math.pow(2,self.nd)
        beta[self.ns+1:self.ns+1+self.nd] *= bdif
        beta[0] /= 2
        for kt in xrange(0,self.nnt+1):
            tt = self.radt[kt]
            s = self.alfa/tt
            bt = beta/tt
            btF = 0
            btF = np.array(bt * FF(s,kon,koff))
            ft[kt] = btF.real.sum()
        return ft
#
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Setting up MPI environment
#
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print "%i is up" % rank
comm.Barrier()
#
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Seting up command line usage
#
parser = OptionParser()
parser.add_option("-d","--diff_const",action="store",dest="diff_const",help="value of diffusion constant",default="NONE")
parser.add_option("-r","--half_area",action="store",dest="half_area",help="value of half-area",default="NONE")
parser.add_option("-i","--input_curve",action="store",dest="input_curve",help="name of input datafile",default="NONE")
parser.add_option("-s","--input_std",action="store",dest="input_std",help="name of input datafile",default="NONE")
parser.add_option("-m","--input_model",action="store",dest="input_model",help="name of the model",default="NONE")
parser.add_option("-o","--output",action="store",dest="out",help="name of output file",default="NONE")
opt,args = parser.parse_args()

# Checking for correct parameters
if opt.input_curve=="NONE":
    print "ATTENTION: No valid input filename given"
    print ""
    print ""
    print "FDAP Fitting Script. 2013-2014. Maxim Igaev, Frederik Suendermann"
    print "Further reading: Igaev et al. 2014. Biophysical Journal, 107(11): 2567-78"
    print ""
    print "Usage example:"
    print "mpiexec -n 8 python residual_analysis_kon_koff.py -i curve.dat -s std_error.dat -m full_model -o output"
    print "   -n:   number of CPUs"
    print "   -d:   value of diffusion constant"
    print "   -r:   value of half-area"
    print "   -i:   normalized input curve (one column *.dat file)"
    print "   -sd:   standard error file (one column *.dat file); NOT MANDATORY"
    print "   -m:   kinetic model (full_model, hybrid, eff, pure, reaction, reaction_pure)"
    print "   -o:   output file's basename (*.txt file)"
    sys.exit()

if opt.input_std=="NONE":
    print ""
    print "WARNING: No standard deviation file given. Proceeding without it..."
    print "You will need to compute confidence intervals manually"
    std = 1.0
else:
    # Loading standard deviation values from a file
    std = np.loadtxt(opt.input_std,dtype=float)

if opt.input_model=="NONE":
    print "ATTENTION: No model indicated"
    print ""
    print ""
    print "FDAP Fitting Script. 2013-2014. Maxim Igaev, Frederik Suendermann"
    print "Further reading: Igaev et al. 2014. Biophysical Journal, 107(11): 2567-78"
    print ""
    print "Usage example:"
    print "mpiexec -n 8 python residual_analysis_kon_koff.py -i curve.dat -s std_error.dat -m full_model -o output"
    print "   -n:   number of CPUs"
    print "   -d:   value of diffusion constant"
    print "   -r:   value of half-area"
    print "   -i:   normalized input curve (one column *.dat file)"
    print "   -s:   standard error file (one column *.dat file); NOT MANDATORY"
    print "   -m:   kinetic model (full_model, hybrid, eff, pure, reaction, reaction_pure)"
    print "   -o:   output file's basename (*.txt file)"
    sys.exit()

if opt.out=="NONE":
    print "ATTENTION: No valid output filename given"
    print ""
    print ""
    print "FDAP Fitting Script. 2013-2014. Maxim Igaev, Frederik Suendermann"
    print "Further reading: Igaev et al. 2014. Biophysical Journal, 107(11): 2567-78"
    print ""
    print "Usage example:"
    print "mpiexec -n 8 python residual_analysis_kon_koff.py -i curve.dat -s std_error.dat -m full_model -o output"
    print "   -n:   number of CPUs"
    print "   -d:   value of diffusion constant"
    print "   -r:   value of half-area"
    print "   -i:   normalized input curve (one column *.dat file)"
    print "   -s:   standard error file (one column *.dat file); NOT MANDATORY"
    print "   -m:   kinetic model (full_model, hybrid, eff, pure, reaction, reaction_pure)"
    print "   -o:   output file's basename (*.txt file)"
    sys.exit()
#
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Defining constans and ranges
#

# Loading FDAP curve
FDAP_experiment = np.loadtxt(opt.input_curve,dtype=float)

# Diffusion constant and activation area's half-length
Df = float(opt.diff_const)
R = float(opt.half_area)

# Defining the (k*_on, k_off)-space
vStep = 0.1 # discretization step

vIni_kon = -5   # initial value
vEnd_kon = +5 # end value
vN_kon = np.int(((vEnd_kon-vIni_kon)/vStep)+1) # number of points

vIni_koff = -5 # initial value
vEnd_koff = +5 # end value
vN_koff = np.int(((vEnd_koff-vIni_koff)/vStep)+1) # number of points

# k*_on and k_off linear arrays
vDegree_kon = np.array([(vIni_kon+(i*vStep)) for i in xrange(vN_kon)])
vDegree_koff = np.array([(vIni_koff+(i*vStep)) for i in xrange(vN_koff)])

# Final remapping of k*_on and k_off onto a log scale
kon = np.power(10,vDegree_kon)
koff = np.power(10,vDegree_koff)

# Time variables for looping (must coincide with the experimental procedure)
tini = 0.001 # initial time in seconds (0 must be avioded to exclude division by zero)
tend = 112   # final time in seconds
timeSteps = 112 # number of time steps within the above interval
#
#---------------------------------------------------------------------


#---------------------------------------------------------------------
# Theoretical FDAP curves
#
def FullModel(s,kon,koff):
    return ((1.0/(1.0+kon/koff))*(1.0+kon/(s+koff))*(1.0/s-1.0/2.0/s/sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)*(1.0-sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df))))+(kon/koff/(1.0+kon/koff))/(s+koff))

def EffectiveDiffusion(s,kon,koff):
    return (1.0/s-1.0/2.0/s/sc.sqrt(R*R*s*(1.0+kon/koff)/Df)*(1.0-sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/koff)/Df))))

def ReactionDominant_Pure(s,kon,koff):
    return (koff/(kon+koff))*((1.0/s-1.0/2.0/s/sc.sqrt(R*R*s/Df)*(1.0-sc.exp(-2.0*sc.sqrt(R*R*s/Df))))) + (kon/(kon+koff))*(1.0/(s+koff))

def ReactionDominant(s,kon,koff):
    return ((kon/(kon+koff))/(s+koff))

def PureDiffusion(s,kon,koff):
    return (1.0/s-1.0/2.0/s/sc.sqrt(R*R*s/Df)*(1.0-sc.exp(-2.0*sc.sqrt(R*R*s/Df))))

def HybridModel(s,kon,koff):
    return ((koff/(s+koff))*(1.0/s-1.0/2.0/s/sc.sqrt(R*R*kon*s/Df/(s+koff))*(1.0-sc.exp(-2.0*sc.sqrt(R*R*kon*s/Df/(s+koff)))))+1.0/(s+koff))
#
#---------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# STARTING THE MAIN ROUTINE
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Ensuring the same package size
dataSize = np.int(round(len(vDegree_kon)*len(vDegree_koff)/size))+1
partMatrix = {}

# Initializing inverse Laplase transform class
tc = TransformClass(tini,tend,timeSteps)

# Starting distributed computing
job = 0l
for j in xrange(len(vDegree_koff)):
    for i in xrange(len(vDegree_kon)):
        job += 1
        if job%size == rank:
            e = np.zeros(timeSteps)
            
            if opt.input_model=="full_model":            
                FDAP_full_model = tc.invlap(FullModel,kon[i],koff[j])
                FDAP_full_model[0] = 1
	        e = (FDAP_full_model - FDAP_experiment)/std
            elif opt.input_model=="eff":
                FDAP_eff = tc.invlap(EffectiveDiffusion,kon[i],koff[j])
                FDAP_eff[0] = 1
	        e = (FDAP_eff - FDAP_experiment)/std
            elif opt.input_model=="pure":
                FDAP_pure = tc.invlap(PureDiffusion,kon[i],koff[j])
                FDAP_pure[0] = 1
	        e = (FDAP_pure - FDAP_experiment)/std
            elif opt.input_model=="hybrid":
                FDAP_hybrid = tc.invlap(HybridModel,kon[i],koff[j])
                FDAP_hybrid[0] = 1
	        e = (FDAP_hybrid - FDAP_experiment)/std
            elif opt.input_model=="reaction":
                FDAP_reaction = tc.invlap(ReactionDominant,kon[i],koff[j])
                FDAP_reaction[0] = 1
	        e = (FDAP_reaction - FDAP_experiment)/std
            elif opt.input_model=="reaction_pure":
                FDAP_reaction_pure = tc.invlap(ReactionDominant_Pure,kon[i],koff[j])
                FDAP_reaction_pure[0] = 1
	        e = (FDAP_reaction_pure - FDAP_experiment)/std

            partMatrix[(j,i)] = (e*e).sum()

if not(len(partMatrix)==dataSize):
    print "I'm the %i th and I'm done!" % rank
    for addLine in xrange(len(partMatrix),dataSize):
        partMatrix[('X','X',addLine)] = np.random.rand()

#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#---------------------------------------------------------------------
# Gathering distributed results
#

data = None
data = comm.gather(partMatrix,root=0)
if rank==0:
    print ""
    print "I am the master"
    print "and I found %i result sets" %len(data)
    rez = np.zeros((len(vDegree_kon),len(vDegree_koff)))

    for i in xrange(len(data)):
        print "\tset %i has a length of %i " % (i,len(data[i]))
#        print "and the following keys:"
#        print data[i].keys()
    print "Done with the calculations... Saving the data now..."
    print ""

    for pack in data:
        for k in pack.keys():
            if not(k[0]=='X'):
                (cKon,cKoff) = k
                rez[cKoff,cKon] = pack[k]
    
    # Idices transposition for kon,koff
    (ikon,ikoff) = np.unravel_index(rez.argmin(),rez.shape)
    coordinates5 = np.argwhere(rez<=5.0)
    coordinates01 = np.argwhere(rez<=0.1)
    print "Exponent values of the minimal k*_on and k_off are %f and %f respectively" %(ikon*vStep-5,ikoff*vStep-5)
    #sio.savemat(opt.out+"_residuals.mat",{'rez':rez})
    #sio.savemat(opt.out+"_coord.mat",{'coordinates':coordinates})
    #np.savetxt(opt.out+"_residuals.txt", rez*vStep)
    np.savetxt(opt.out+"_coordinates5.txt", coordinates5*vStep)
    np.savetxt(opt.out+"_coordinates01.txt", coordinates01*vStep)


    #---------------------------------------------------------------------
    # Recording the best FDAP fit + writing the curve
    #
    #FDAP = tc.invlap(FullModel,kon[ikon],koff[ikoff])
    #FDAP[0] = 1
    #fh = open(opt.out+"_best_fit.txt","w")
    #for i in xrange(len(FDAP)):
    #    fh.write("%f\n" %FDAP[i])
    #fh.close()
    #
    #---------------------------------------------------------------------


    #---------------------------------------------------------------------
    # Calculating residuals for the best FDAP fit + writing the residuals
    #
    #e = (FDAP - FDAP_experiment)/std
    #fh = open(opt.out+"_residuals.txt","w")
    #for i in xrange(len(e)):
    #    fh.write("%f\n" %e[i])
    #fh.close()
    #
    #---------------------------------------------------------------------


    #---------------------------------------------------------------------
    # Jacobian matrix and error estimation (manual calculation)
    # Works only for fitting with FullModel so far
    #

    # Derivatives
    def JJ_full_model_kon(s,kon,koff):
        return -koff*(-kon-koff+koff*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + kon*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 2.0*koff*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 2.0*kon*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) - 2.0*s + 2.0*s*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)))/4.0/s/(s+koff)/(kon+koff)/(kon+koff)/sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)

    def JJ_full_model_koff(s,kon,koff):
        return kon*(-koff**2.0-2.0*s**2.0-4.0*s*koff-kon*koff + kon*koff*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 2.0*koff*koff*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 2.0*s*s*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 4.0*s*koff*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + koff*koff*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 2.0*s*kon*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) + 2.0*kon*koff*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)*sc.exp(-2.0*sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)) - 2*s*kon)/4.0/s/(s+koff)/(s+koff)/(kon+koff)/(kon+koff)/sc.sqrt(R*R*s*(1.0+kon/(s+koff))/Df)

    # Confidence intervals
    J_kon = tc.invlap(JJ_full_model_kon,kon[ikon],koff[ikoff])
    J_koff = tc.invlap(JJ_full_model_koff,kon[ikon],koff[ikoff])
    J = zip(J_kon,J_koff)

    df = len(FDAP_experiment) - 2
    s2 = (linalg.norm(e))**2/df
    v = s2*linalg.inv(multmat(np.transpose(J),J))
    error = np.absolute(sc.stats.t.isf(0.95,df))*sc.sqrt(np.diag(v))
    #
    #---------------------------------------------------------------------


    #---------------------------------------------------------------------
    # Writing the result file
    #
    fh = open(opt.out+"_parameters.txt","w")
    fh.write("bound = %f +- %f\n" %(100.0 - 100.0/(1.0 + kon[ikon]/koff[ikoff]),100.0*(error[0]/koff[ikoff] - kon[ikon]*error[1]/koff[ikoff]/koff[ikoff])/(1.0 + kon[ikon]/koff[ikoff])/(1.0 + kon[ikon]/koff[ikoff])))
    fh.write("k_on = %f +- %f\n" %(kon[ikon],error[0]))
    fh.write("exponent of k_on = %f\n" %(ikon*vStep-vEnd_kon))
    fh.write("k_off = %f +- %f\n" %(koff[ikoff],error[1]))
    fh.write("exponent of k_off = %f\n" %(ikoff*vStep-vEnd_koff))
    fh.write("chi2-value = %f\n" %rez[ikon,ikoff])
    fh.close()
    #
    #---------------------------------------------------------------------

#
#---------------------------------------------------------------------

    print "...done\nEXIT"
#    print rez
