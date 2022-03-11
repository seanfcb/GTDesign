import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from tabulate import tabulate
from CompressibleFlowFunctions.Isentropic import *
from scipy.optimize import *


#####################################################################################
#                             Functions in use                                      #
#####################################################################################
tab1 = []
def table1(*args):
    tab1.append(args)

def eqn6(M1,C1,T01,gamma,Rs):
    return C1/np.sqrt(T01) - np.sqrt(gamma*Rs)*M1*(1+(gamma-1)/2*M1**2)**(-0.5)

#####################################################################################
#                              Input parameters (from A)                            #
#                             [P] = kPa, [T] = K                                    #
#####################################################################################

N     = 13700 #RPM
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

HTR = np.linspace(0.4,0.8,100) #Hub to tip ratio, typical values for centrifugal compressors
alpha1 = [10.96*np.pi/180] #np.linspace(5,35, int((35-5)/5+1))*np.pi/180 #given alpha1 angle range[20*np.pi/180]
r1h = 0.24920391555607968 #coordinate values from WebplotDigitizer
r1T = 0.3481778511616934 #coordinate values from WebplotDigitizer
CT_ratio = r1h/r1T #Ratio of impeller inlet hub radius to turbine hub radius
#CT_ratio = 0.87778
r1T = 23.218/100 #Value calculated in part C, m
r1h  = CT_ratio*r1T #Impeller hub radius
#r1h  = 15/100
r1sh = r1h/HTR
Mrel1sh = 0*r1sh
U1sh = 0*r1sh
U1m  = 0*r1sh
C1x  = 0*r1sh
spec = 0*r1sh
#####################################################################################
#                              Begin solution                                       #
#####################################################################################
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
            #spec[rcount] = N*2*np.pi/60*np.sqrt(C1x[rcount]*A_inlet)*(1/(P02-P01))**(3/4)
            spec[rcount] = N*2*np.pi/60*np.sqrt(C1x[rcount]*A_inlet)*(rho1/(P02-P01))**(3/4)
            #print(spec)
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
    rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
    rc('text', usetex=True)
    plt.plot(C1x/U1sh,Mrel1sh, label = r'$\alpha_1 = {}$'.format(round(i*180/np.pi,3)))
    plt.xlabel(r"$C_{1x}/U_{1sh}$")
    plt.ylabel(r"$M_{rel,1sh}$")
    plt.legend()

#########TESTING INDENT AS OF HERE
fig.savefig('Mrel1sh_vs_C1x.png', dpi=fig.dpi)
#plt.show()
Mach_index = np.argmin(Mrel1sh)
alpha1  = alpha1[0]
HTR     = HTR[Mach_index]
C1x     = C1x[Mach_index]
Mrel1sh = Mrel1sh[Mach_index]
r1sh    = r1sh[Mach_index]
U1sh    = U1sh[Mach_index]
spec    = spec[Mach_index]

#####################################################################################
#                                 Hub properties                                    #
#####################################################################################
U1h    = np.pi/180*N*r1h
W1h    = np.sqrt(C1x**2 + (U1h-C1x*np.tan(alpha1))**2)
beta1h = np.arccos(C1x/W1h)*180/np.pi

#####################################################################################
#                                 Mean properties                                   #
#####################################################################################
r1m     = (r1sh+r1h)/2
U1m     = np.pi/180*N*r1m
W1m     = np.sqrt(C1x**2 + (U1m-C1x*np.tan(alpha1))**2)
beta1m  = np.arccos(C1x/W1m)*180/np.pi


#####################################################################################
#                                 Shroud properties                                 #
#####################################################################################
U1sh    = U1h*r1sh/r1h
W1sh    = np.sqrt(C1x**2 + (U1sh-C1x*np.tan(alpha1))**2)
beta1sh = np.arccos(C1x/W1sh)*180/np.pi



table1("Location", "alpha1","beta1", "i1", "N (RPM)","Specific Speed", "Mrel1sh", "r1 (m)", "C1x","C1","U1","W1")
table1("Hub",round(alpha1*180/np.pi,4),beta1h,alpha1*180/np.pi-beta1h,N,spec,Mrel1sh,r1h,C1x,C1,U1h,W1h)
table1("Mean",round(alpha1*180/np.pi,4),beta1m,alpha1*180/np.pi-beta1m,N,spec,Mrel1sh,r1m,C1x,C1,U1m,W1m)
table1("Shroud",round(alpha1*180/np.pi,4),beta1sh,alpha1*180/np.pi-beta1sh,N,spec,Mrel1sh,r1sh,C1x,C1,U1sh,W1sh)
print(tabulate(tab1))
