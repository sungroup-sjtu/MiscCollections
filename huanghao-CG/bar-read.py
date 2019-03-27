#!/home/gongzheng/local/bin/apy
# coding=utf-8

import os, sys, pymbar
import numpy as np
from pymbar import timeseries
from pymbar.bar import BAR
from optparse import OptionParser

#================================================================================
# CONSTANTS AND PARAMETERS
#===============================================================================

kB = 1.381*6.02214/1000.0 # Boltzmann's constant (kJ/mol/K).
relative_tolerance = 1e-10 # Convergence criterion of the energy estimates for BAR and MBAR.
initialize_with = 'BAR' # The initial MBAR free energy guess is set to BAR.
# methods = ['BAR', 'MBAR'] # Free energy estimation methods.
methods = ['BAR'] # Free energy estimation methods.

#=================================================================================
# INPUT OPTIONS
#==================================================================================

parser = OptionParser()
parser.add_option('-t', '--temperature', dest = 'temperature', help = "Temperature in K. Default: 298.15 K.", default = 298.15, type=float)
parser.add_option('-u', '--units', dest = 'units', help = 'Units to report energies: \'kJ\', \'kcal\', and \'kBT\'. Default: \'kcal\'', default = 'kcal')
parser.add_option('-v', '--verbose', dest = 'verbose', help = 'Verbose option for BAR and MBAR. Default: False.', default = False, action = 'store_true')

(options, args) = parser.parse_args()
units = options.units
verbose = options.verbose

#====================================================================================
# MAIN
#===================================================================================

# Define units the energies will be reported in.
beta = 1./(kB*options.temperature)
if units == 'kJ':
    beta_report = beta
    units = '(kJ/mol)'
elif units == 'kcal':
    beta_report = 4.184*beta
    units = '(kcal/mol)'
elif units == 'kBT':
    beta_report = 1
    units = '(k_BT)'
else:
    parser.error('\nThe only options: \'kJ\', \'kcal\', \'kBT\'')

# Read lambda list
l_list = []
for dir in os.listdir(os.getcwd()):
    if not (dir.startswith('bar-') and os.path.isdir(dir)):
        continue
    lam = float(dir.split('-')[1])
    l_list.append(lam)
l_list.sort()
K = len(l_list)
maxn = len(np.genfromtxt('bar-%.2f/bar-%.2f-dl.txt' %(l_list[0],l_list[0])))
u_klt = np.zeros([K,K,maxn], np.float64)      # u_klt[k,m,t] is the reduced potential energy of snapshot t of state k evaluated at state m

for k in range(K):
    if k != K-1:
        data = np.genfromtxt('bar-%.2f/bar-%.2f-dl.txt' %(l_list[k],l_list[k]))
        u_klt[k,k+1,:] = data[:,1]*beta*4.184
    if k != 0:
        data = np.genfromtxt('bar-%.2f/bar-%.2f-ndl.txt' %(l_list[k],l_list[k]))
        u_klt[k,k-1,:] = data[:,1]*beta*4.184

#======================================================================================
# The functions we'll be using presently.
#=====================================================================================
def getNkandUkln():
    # u_kln = u_klt
    # N_k = [maxn]*K
    # return (N_k, u_kln)

    """Identifies uncorrelated samples and updates the arrays of the reduced potential energy and dhdlt retaining data entries of these samples only."""
    u_kln = np.zeros([K,K,maxn], np.float64) # u_kln[k,m,n] is the reduced potential energy of uncorrelated sample index n from state k evaluated at state m
    N_k = np.zeros(K, int) # N_k[k] is the number of uncorrelated samples from state k
    g = np.zeros(K,float) # autocorrelation times for the data
    print "Number of correlated and uncorrelated samples:\n\n%8s %10s %12s %12s" % ('Lambda', 'N', 'N_k', 'N/N_k')
    for k in range(K):
        if k == 0:
            g[k] = timeseries.statisticalInefficiency(u_klt[k,k+1,:])
            indices = np.array(timeseries.subsampleCorrelatedData(u_klt[k,k+1,:])) # indices of uncorrelated samples
        else:
            g[k] = timeseries.statisticalInefficiency(u_klt[k,k-1,:])
            indices = np.array(timeseries.subsampleCorrelatedData(u_klt[k,k-1,:]))
        N = len(indices) # number of uncorrelated samples
        N_k[k] = N # Store the number of uncorrelated samples from state k.
        for l in range(K):
            u_kln[k,l,0:N] = u_klt[k,l,indices]
        print "%6.2f %12s %12s %12.2f" % (l_list[k], maxn, N_k[k], g[k])
    print ''
    return (N_k, u_kln)

def estimatewithMBAR(reltol=relative_tolerance, regular_estimate=False):
    """Computes the MBAR free energy given the reduced potential and the number of relevant entries in it."""
    MBAR = pymbar.MBAR(u_kln, N_k, verbose = verbose, method = 'adaptive', relative_tolerance = reltol, initialize = initialize_with)
    # Get matrix of dimensionless free energy differences and uncertainty estimate.
    (Deltaf_ij, dDeltaf_ij) = MBAR.getFreeEnergyDifferences(uncertainty_method='svd-ew')

    if (verbose): 
        print "The overlap matrix is..."
        O = MBAR.computeOverlap('matrix')
        for k in range(K):
            line = ''
            for l in range(K):
                line += ' %5.2f ' % O[k, l]
            print line
    if regular_estimate:
        return (Deltaf_ij, dDeltaf_ij)
    return (Deltaf_ij[0,K-1]/beta_report, dDeltaf_ij[0,K-1]/beta_report)

#========================================================================================
# Calculate free energies with different methods.
#=======================================================================================

# Obtain uncorrelated samples
N_k, u_kln = getNkandUkln()
# Estimate free energy difference with MBAR -- all states at once.
if 'MBAR' in methods:
    print "Estimating the free energy change with MBAR..."
    Deltaf_ij, dDeltaf_ij = estimatewithMBAR(regular_estimate=True)

print ("Estimating the free energy change with %s..." % ', '.join(methods)).replace(', MBAR', '')
df_allk = list(); ddf_allk = list()

for k in range(K-1):
    df = dict(); ddf = dict()
    w_F = u_kln[k,k+1,0:N_k[k]] - u_kln[k,k,0:N_k[k]] 
    w_R = u_kln[k+1,k,0:N_k[k+1]] - u_kln[k+1,k+1,0:N_k[k+1]] 

    for name in methods:
        if name == 'BAR':
            (df['BAR'], ddf['BAR']) = pymbar.bar.BAR(w_F, w_R, relative_tolerance=relative_tolerance, verbose = verbose)        
        if name == 'MBAR':
            (df['MBAR'], ddf['MBAR']) =  Deltaf_ij[k,k+1], dDeltaf_ij[k,k+1]

    df_allk = np.append(df_allk,df)
    ddf_allk = np.append(ddf_allk,ddf)

#=========================================================================================
# All done with calculations, now summarize and print stats.
#=========================================================================================   

dF  = dict.fromkeys(methods, 0)
ddF = dict.fromkeys(methods, 0)

for name in methods:
     if name == 'MBAR':
         dF['MBAR']  =  Deltaf_ij[0, K-1]
         ddF['MBAR'] = dDeltaf_ij[0, K-1]

     else:
         for k in range(0, K-1):
             dF[name] += df_allk[k][name]
             ddF[name] += (ddf_allk[k][name])**2
         ddF[name] = np.sqrt(ddF[name])

# Display results.
def printLine(str1, str2, d1=None, d2=None):
    """Fills out the results table linewise."""
    print str1,
    text = str1
    for name in methods:
        if d1 == 'plain':
            print str2,
            text += ' ' + str2
        if d1 == 'name':
            print str2 % (name, units),
            text += ' ' + str2 % (name, units)
        if d1 and d2:
            print str2 % (d1[name]/beta_report, d2[name]/beta_report),
            text += ' ' + str2 % (d1[name]/beta_report, d2[name]/beta_report)
    print ''
    outtext.append(text + '\n')
    return

outtext = []
printLine(12*'-', 21*'-', 'plain')
printLine('%-12s' % '   States', '%9s %-11s', 'name')
printLine(12*'-', 21*'-', 'plain')
for k in range(K-1):
    printLine('%5.2f -> %3.2f' % (l_list[k], l_list[k+1]), '%10.3f +- %6.3f', df_allk[k], ddf_allk[k])
printLine(12*'-', 21*'-', 'plain')
printLine('%9s    ' %'TOTAL:', '%10.3f +- %6.3f', dF, ddF)

