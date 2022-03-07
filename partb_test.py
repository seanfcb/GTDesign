import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from tabulate import tabulate
from CompressibleFlowFunctions.Isentropic import *
from scipy.optimize import *


#####################################################################################
#                             Functions in use                                      #
#####################################################################################


def eqn6(M1,C1,T01,gamma,Rs):
    return C1/np.sqrt(T01) - np.sqrt(gamma*Rs)*M1*(1+(gamma-1)/2*M1**2)**(-0.5)

#####################################################################################
#                              Input parameters (from A)                            #
#                             [P] = kPa, [T] = K                                    #
#####################################################################################

N     = 6300 #RPM
gamma = 1.4
mdot  = 15/2.205
Rs    = 287
P03   = 79.032
P04   = 237.096
T03   = 348.395
T04   = 500.426

##For simplicity, we denote station 3 and 4 from part A as stations 1 and 2 for part B
T01 = T03
T02 = T04
P01 = P03
P02 = P04

#####################################################################################
#                              Input parameters (from B statement)                  #
#####################################################################################

HTR = np.linspace(0.4,0.7,100) #Hub to tip ratio, typical values for centrifugal compressors
alpha1 = [20*np.pi/180]#np.linspace(5,35,100)*np.pi/180 #given alpha1 angle range
r1h = 0.194852941 #coordinate values from WebplotDigitizer
r1T = 0.268382353 #coordinate values from WebplotDigitizer
CT_ratio = r1h/r1T #Ratio of impeller inlet hub radius to turbine hub radius
r1T = 50.64758741/100 #Value calculated in part C, m
#r1h  = CT_ratio*r1T #Impeller hub radius
r1h  = 15/100
r1sh = r1h/HTR
Mrel1sh = 0*r1sh
U1sh = 0*r1sh
U1m  = 0*r1sh
C1x  = 0*r1sh
##Begin solutionf
fig = plt.figure()

for i in alpha1: #iterating for every alpha value
    rcount = 0
    for j in r1sh: #iterating through range of inlet shroud radii
        U1m[rcount]   = (j+r1h)/2*2*np.pi/60*N
        rho01 = P01/(Rs/1000*T01) #assume rho1 = rho01 for first iteration
        rho1  = rho01
        A_inlet = np.pi*(j**2-r1h**2)
        delrho  = 1

        while delrho > 1e-5:
            C1x[rcount] = mdot/(rho1*A_inlet)
            spec = N*2*np.pi/60*np.sqrt(C1x[rcount]*A_inlet)*(1/(P02-P01))**(3/4)
            print(spec)
            C1  = C1x[rcount]/np.cos(i)
            M1  = bisect(eqn6,0,0.99999,args=(C1,T01,gamma,Rs))
            c1  = C1/M1

            rhoi = rho1
            rho1 = rho01/(1+(gamma-1)/2*M1**2)**(1/(gamma-1))

            delrho = np.abs(rhoi - rho1)

        U1sh[rcount] = 2*np.pi*j*N/60
        W1sh = np.sqrt(C1x[rcount]**2+(U1sh[rcount]-C1x[rcount]*np.tan(i))**2)
        Mrel1sh[rcount] = W1sh/c1
        rcount = rcount + 1
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.plot(C1x/U1sh,Mrel1sh, label = 'alpha = {}'.format(i*180/np.pi))
    plt.xlabel(r"$C_{1x}/U_{1sh}$")
    plt.ylabel(r"$M_{rel,1sh}$")
    fig.savefig('Mrel1sh_vs_r1sh.png', dpi=fig.dpi)
