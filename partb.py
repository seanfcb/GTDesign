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
tab2 = []
def table1(*args):
    tab1.append(args)
def table2(*args):
    tab2.append(args)

def eqn6(M1,C1,T01,gamma,Rs):
    return C1/np.sqrt(T01) - np.sqrt(gamma*Rs)*M1*(1+(gamma-1)/2*M1**2)**(-0.5)

#####################################################################################
#                              Input parameters (from A)                            #
#                             [P] = kPa, [T] = K                                    #
#####################################################################################

Cpa   = 1005
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

HSR = np.linspace(0.4,0.8,100) #Hub to tip ratio, typical values for centrifugal compressors
alpha1 = [(25)*np.pi/180] #np.linspace(5,35, int((35-5)/5+1))*np.pi/180 #given alpha1 angle range[20*np.pi/180] ###10.96
r1h = 0.24920391555607968 #coordinate values from WebplotDigitizer
r1T = 0.3481778511616934 #coordinate values from WebplotDigitizer
CT_ratio = r1h/r1T #Ratio of impeller inlet hub radius to turbine hub radius
#CT_ratio = 0.87778
r1T = 23.218/100 #Value calculated in part C, m
r1h  = CT_ratio*r1T #Impeller hub radius
#r1h  = 20/100
r1sh = r1h/HSR
Mrel1sh = 0*r1sh
U1sh = 0*r1sh
U1m  = 0*r1sh
C1x  = 0*r1sh
C1   = 0*r1sh
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
            C1[rcount]  = C1x[rcount]/np.cos(i)
            M1  = bisect(eqn6,0,0.99999,args=(C1[rcount],T01,gamma,Rs))
            c1  = C1[rcount]/M1

            rhoi = rho1
            rho1 = rho01/(1+(gamma-1)/2*M1**2)**(1/(gamma-1))

            delrho = np.abs(rhoi - rho1)

        U1sh[rcount] = 2*np.pi*j*N/60
        W1sh = np.sqrt(C1x[rcount]**2+(U1sh[rcount]-C1x[rcount]*np.tan(i))**2)
        Mrel1sh[rcount] = W1sh/c1
        rcount = rcount + 1
    rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
    rc('text', usetex=True)
    plt.plot(C1x,Mrel1sh, label = r'$\alpha_1 = {}$'.format(round(i*180/np.pi,3)))
    plt.xlabel(r"$C_{1x}/U_{1sh}$")
    plt.ylabel(r"$M_{rel,1sh}$")
    plt.legend()

#fig.savefig('Mrel1sh_vs_C1x.png', dpi=fig.dpi)
#plt.show()
Mach_index = np.argmin(Mrel1sh)
alpha1  = alpha1[0]
HSR     = HSR[Mach_index]
C1x     = C1x[Mach_index]
C1      = C1[Mach_index]
Mrel1sh = Mrel1sh[Mach_index]
r1sh    = r1sh[Mach_index]
#U1sh    = U1sh[Mach_index]
spec    = spec[Mach_index]

#####################################################################################
#                                 Hub inlet properties                              #
#####################################################################################
U1h         = 2*np.pi/60*N*r1h
W1h         = np.sqrt(C1x**2 + (U1h-C1x*np.tan(alpha1))**2)
beta1h      = np.arccos(C1x/W1h)*180/np.pi #relative flow angle
ih_min      = 5
beta1h_star = beta1h - ih_min #blade metal angle

#####################################################################################
#                                 Mean inlet properties                             #
#####################################################################################
r1m         = (r1sh+r1h)/2
U1m         = 2*np.pi/60*N*r1m
W1m         = np.sqrt(C1x**2 + (U1m-C1x*np.tan(alpha1))**2)
beta1m      = np.arccos(C1x/W1m)*180/np.pi #relative flow angle
im_min      = 2
beta1m_star = beta1m - im_min #blade metal angle

#####################################################################################
#                                 Shroud inlet properties                           #
#####################################################################################
U1sh         = U1h*r1sh/r1h
W1sh         = np.sqrt(C1x**2 + (U1sh-C1x*np.tan(alpha1))**2)
beta1sh      = np.arccos(C1x/W1sh)*180/np.pi #relative flow angle
ish_min      = 1
beta1sh_star = beta1h - ish_min #blade metal angle

#####################################################################################
#                                 Display inlet properties                          #
#####################################################################################

# table1("Location", "alpha1","beta1", "i1", "N (RPM)","Specific Speed", "Mrel1sh", "r1 (m)", "C1x","C1","U1","W1")
# table1("Hub",round(alpha1*180/np.pi,4),beta1h,alpha1*180/np.pi-beta1h,N,spec,Mrel1sh,r1h,C1x,C1,U1h,W1h)
# table1("Mean",round(alpha1*180/np.pi,4),beta1m,alpha1*180/np.pi-beta1m,N,spec,Mrel1sh,r1m,C1x,C1,U1m,W1m)
# table1("Shroud",round(alpha1*180/np.pi,4),beta1sh,alpha1*180/np.pi-beta1sh,N,spec,Mrel1sh,r1sh,C1x,C1,U1sh,W1sh)
table1("Location", "alpha1","beta1", "Blade metal angle", "Incidence", "N (RPM)","Specific Speed", "Mrel1sh", "r1 (m)", "C1x","C1","U1","W1")
table1("Hub",round(alpha1*180/np.pi,4),beta1h,beta1h_star,ih_min,N,spec,Mrel1sh,r1h,C1x,C1,U1h,W1h)
table1("Mean",round(alpha1*180/np.pi,4),beta1m,beta1m_star,im_min,N,spec,Mrel1sh,r1m,C1x,C1,U1m,W1m)
table1("Shroud",round(alpha1*180/np.pi,4),beta1sh,beta1sh_star,ish_min,N,spec,Mrel1sh,r1sh,C1x,C1,U1sh,W1sh)

print(tabulate(tab1))

#####################################################################################
#                                 Impeller exit conditions                          #
#####################################################################################
W_HPC= 1.0398*1e6/mdot
U2_i = 1250*np.sqrt(T01*1.8/518.7)/3.281
U2   = [U2_i]
r2   = 60*U2[0]/(N*2*np.pi)
L    = r2 - r1h #impeller axial length
##Redefining values
# r2 = r1h/np.linspace(0.4,0.6,int(3000+1))
# U2 = N*2*np.pi/60*r2
sig = 1-np.pi*0.63/(16+16)
beta2_star = np.linspace(0,40,int(40+1))*np.pi/180#[10*np.pi/180]#
table2("Crat","Wrat", "beta2*","Ct2","W2","Ct2i","Wt2i","W2i","U2")
for j in U2:
    for i in beta2_star:
        Ct2  = (Cpa*(T02-T01))/j #+U1m*W1m
        Ct2i = Ct2/sig
        Wt2i = (j - Ct2i)
        W2i  = Wt2i/np.sin(i)
        W2   = W2i
        Wrat = W2/W1sh#np.abs(W2/W1sh)
        if Wrat > 0.5 and Wrat < 0.6:
            Cr2 = np.sqrt(W2**2-Wt2i**2)
            Crat = Cr2/C1x
            if Crat > 0.8 and Crat < 1 and Wrat > 0.5 and Wrat < 0.6:
                table2(round(np.abs(Crat),5),round(Wrat,5), round(i*180/np.pi,3),round(Ct2,3),round(W2,3),round(Ct2i,3),round(Wt2i,3),round(W2i,3),round(j,3))
print(tabulate(tab2))

############Taking velocities and angles from table2
U2   = 418.907
Crat = 0.88916
Wrat = 0.50841
beta2_star = 12
Ct2  = 364.737
W2   = 144.886
Ct2i = 388.784
Wt2i = 30.123
W2i  = 144.886

###########Calculating remaining velocities and angles in outlet triangles
##Real flow triangle
Wt2    = U2 - Ct2
Cr2    = np.sqrt(W2**2-Wt2**2)
C2     = np.sqrt(Cr2**2+Ct2**2)
a2     = np.sqrt(gamma*Rs*T02)
M2     = C2/a2
alpha2 = np.arctan(Ct2/Cr2)*180/np.pi
beta2  = np.arctan(Wt2/Cr2)*180/np.pi
##Ideal flow triangle
Cr2i    = Crat*C1x
C2i     = np.sqrt(Cr2i**2+Ct2i**2)
alpha2i = np.arctan(Ct2i/Cr2i)*180/np.pi

phi   = Cr2/U2 #Flow coefficient
psi   = W_HPC/(sig*U2**2)#Aerodynamic loading
delta = beta2_star - beta2#Deviation

#####################################################################################
#                                 Impeller blockage                                 #
#####################################################################################
T2   = T02*(1+(gamma-1)/2*M2**2)**(-1)
P2   = P02*(T2/T02)**(gamma/(gamma-1))#Static pressure at 2
rho2 = P2/(Rs/1000*T2)
B2_star_aero  = (0.12+0.17)/2
B2_star_blade = (0.03+0.06)/2
B2_star       = B2_star_aero+B2_star_blade
CD2           = 1-B2_star
b2            = mdot/(rho2*Cr2*2*np.pi*r2*CD2)
