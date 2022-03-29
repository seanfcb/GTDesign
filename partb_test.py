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
tab3 = []
def table1(*args):
    tab1.append(args)
def table2(*args):
    tab2.append(args)
def table3(*args):
    tab3.append(args)

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
#fig = plt.figure()

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
    # rc('font', **{'family': 'serif', 'serif': ['Times New Roman']})
    # rc('text', usetex=True)
    # plt.plot(C1x,Mrel1sh, label = r'$\alpha_1 = {}$'.format(round(i*180/np.pi,3)))
    # plt.xlabel(r"$C_{1x}/U_{1sh}$")
    # plt.ylabel(r"$M_{rel,1sh}$")
    # plt.legend()

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
beta2_star = -np.linspace(0,40,int(40+1))*np.pi/180#[10*np.pi/180]#
table2("Crat","Wrat", "beta2*","Ct2","W2","Ct2i","Wt2i","W2i","U2")
for j in U2:
    for i in beta2_star:
        Ct2  = (Cpa*(T02-T01))/j #+U1m*W1m
        Ct2i = Ct2/sig
        Wt2i = (j - Ct2i)
        W2i  = Wt2i/np.sin(i)
        W2   = W2i
        Wrat = np.abs(W2/W1sh)#W2/W1sh#
        #print(Wrat)
        if Wrat > 0.5 and Wrat < 0.6:
            Cr2 = np.sqrt(W2**2-Wt2i**2)
            Crat = Cr2/C1x
            #print(Wrat,Crat,i*180/np.pi)

            if Crat > 0.8 and Crat < 1 and Wrat > 0.5 and Wrat < 0.6:
                table2(round(np.abs(Crat),5),round(Wrat,5), round(i*180/np.pi,3),round(Ct2,3),round(W2,3),round(Ct2i,3),round(Wt2i,3),round(W2i,3),round(j,3))
#print(tabulate(tab2))

############Taking velocities and angles from table2
# U2   = 418.907
# Crat = 0.88916
# Wrat = 0.50841
# beta2_star = 12
# Ct2  = 364.737
# W2   = 144.886
# Ct2i = 388.784
# Wt2i = 30.123
# W2i  = 144.886

U2   = 418.907
Crat = 0.88916
Wrat = 0.50841
beta2_star = -12
Ct2  = 364.737
W2   = -144.886
Ct2i = 388.784
Wt2i = 30.123
W2i  = -144.886


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
T2            = T02*(1+(gamma-1)/2*M2**2)**(-1)
P02           = P01*(1+0.88*(T02/T01-1))**(gamma/(gamma-1))#Total pressure at impeller exit
P2            = P02*(T2/T02)**(gamma/(gamma-1))#Static pressure at 2
rho2          = P2/(Rs/1000*T2)
B2_star_aero  = (0.17)
B2_star_blade = (0.06)
B2_star       = B2_star_aero+B2_star_blade
CD2           = 1-B2_star
b2            = mdot/(rho2*Cr2*2*np.pi*r2*CD2)
A2            = 2*np.pi*r2*b2

#####################################################################################
#                                Diffuser inlet                                     #
#####################################################################################

r3r2 = 1.05 #r3/r2, [1.05;1.10]. Smallest value selected to minimize gaspath length
r3   = r3r2*r2
b3   = b2
A3   = 2*np.pi*r3*b3
P03  = 0.99*P02
T03  = T02
rho03 = P03/(Rs/1000*T03)
#B3_star = 0.125 #Assumption


Ct3  = -Ct2*r2/r3

def station3(Cr3,Ct3,P02,P03,Cpa,Rs,mdot,rho2,r2,r3,Cr2):
    C3 = np.sqrt(Cr3**2+Ct3**2)
    T3 = T03-C3**2/(2*Cpa)
    P3 = P03*(T3/T03)**(gamma/(gamma-1))
    rho3 = P3/(Rs/1000*T3)
    #print(C3)

    # print(Cr2*r2*rho2/(r3*rho3))
    return Cr3 - Cr2*r2*rho2/(r3*rho3)

Cr3 = bisect(station3,0.0001,600, args=(Ct3,P02,P03,Cpa,Rs,mdot,rho2,r2,r3,Cr2)) #Iterating on the previous function to determine optimal Cr3
C3 = np.sqrt(Cr3**2+Ct3**2)
T3 = T03-C3**2/(2*Cpa)
P3 = P03*(T3/T03)**(gamma/(gamma-1))
rho3 = P3/(Rs/1000*T3)

def Mach3(M,T03,gamma,Rs):
    return C3/np.sqrt(T03) - np.sqrt(gamma*Rs)*M*(1+(gamma-1)/2*M**2)**(-0.5)
M3 = bisect(Mach3,0.00001,0.9999999,args=(T03,gamma,Rs))
alpha3 = -np.arctan(Cr3/Ct3)*180/np.pi
i3 = -1 #degrees
alpha3star = alpha3+i3

table3("r3","b3", "A3","P03","P3","T03","T3","rho03","rho3","Ct3","Cr3","C3","M3","alpha3","i3","alpha3*")
table3(round(r3,3),round(b3,3),round(A3,3),round(P03,3),round(P3,3),round(T03,3),round(T3,3),round(rho03,3),round(rho3,3),round(Ct3,3),round(Cr3,3),round(C3,3),round(M3,3),round(alpha3,3),round(i3,3),round(alpha3star,3))
#print(tabulate(tab3))

#####################################################################################
#                                Diffuser throat                                    #
#####################################################################################
AR_star =0.9
M_star  = 1
b_star  = b3
w_star  = b_star*AR_star*1.0125

M4      = 0.1 #Diffuser exit condition

##From the isentropic flow relations
Tstar_ratio   = (1+((gamma-1)/2)*1**2)**(-1)
Pstar_ratio   = (1+((gamma-1)/2)*1**2)**(-gamma/(gamma-1))
Rhostar_ratio = (1+((gamma-1)/2)*1**2)**(-1/(gamma-1))
Astar_ratio   = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*1**2)**((gamma+1)/(2*(gamma-1)))/1

T4_ratio   = (1+((gamma-1)/2)*M4**2)**(-1)
P4_ratio   = (1+((gamma-1)/2)*M4**2)**(-gamma/(gamma-1))
Rho4_ratio = (1+((gamma-1)/2)*M4**2)**(-1/(gamma-1))
A4_ratio   = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M4**2)**((gamma+1)/(2*(gamma-1)))/M4
A3_ratio   = ((gamma+1)/2)**(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M3**2)**((gamma+1)/(2*(gamma-1)))/M3

T4Tstar     = T4_ratio/Tstar_ratio
P4Pstar     = P4_ratio/Pstar_ratio
rho4rhostar = Rho4_ratio/Rhostar_ratio
A4Astar     = A4_ratio/Astar_ratio

# #Exit static conditions
# P4 = P04*P4_ratio
# T4 = T04*T4_ratio
# rho4 = P4/(Rs/1000*T4)
#From Slater
P4 = 3*P01
P04 = P4/P4_ratio
T04 = T03
T4 = T04*T4_ratio
rho4 = P4/(Rs/1000*T4)

#Throat conditions
T0star  = T03
Tstar   = T0star*Tstar_ratio
P0star  = P03#Pstar/Pstar_ratio
Pstar   = P0star*Pstar_ratio#P4/P4Pstar


phip   = (Pstar-P3)/(P03-P3) #LE to throat static pressure recovery coefficient
Bstar_throat = 0.04
#Astar = mdot*np.sqrt(T0star*Rs/1000/gamma)/(1*(1+(gamma-1)/2*1**2)**((gamma+1)/(2-2*gamma))*P0star*(1-Bstar_throat))

Astar =A3/A3_ratio/(1-Bstar_throat)#Including blockage, assuming the total conditions at 3 are the same as the total conditions at choke

#Astar= mdot*1.0125/P0star*np.sqrt(Rs/1000*T0star/gamma)*np.sqrt((gamma+1)/2)**((gamma+1)/(gamma-1))/(1-Bstar_throat)

#####################################################################################
#                                Diffuser exit                                      #
#####################################################################################
r4r3 = 2 #Radius ratio
r4 = r3*r4r3
A4  = A4_ratio*Astar
Nv = 21
phis= 360/Nv

#Throat to exit static pressure recovery
Cp = (P4-Pstar)/(P0star-Pstar)

omega = 100


X = np.linspace(0,w_star,10)#np.linspace(0,w_star,10)
Y1 = -np.tan((90-(alpha3star-phis+omega/2))*np.pi/180)*X +r3*(np.tan((90-(alpha3star-phis+omega/2))*np.pi/180)*np.sin(phis*np.pi/180)+np.cos(phis*np.pi/180))
V1cl = -np.tan((90-alpha3star)*np.pi/180)*X+r3
V2cl = -np.tan((90-(alpha3star-phis))*np.pi/180)*X+r3*(np.tan((90-(alpha3star-phis))*np.pi/180)*np.sin(phis*np.pi/180)+np.cos(phis*np.pi/180))
V2ss = -np.tan((90-(alpha3star-phis+omega/2))*np.pi/180)*X+r3*(np.tan((90-(alpha3star-phis+omega/2))*np.pi/180)*np.sin(phis*np.pi/180)+np.cos(phis*np.pi/180))

fig2 = plt.figure()
#plt.plot(X,Y1)
# plt.plot(X,V1cl,"b--")
plt.plot(X,V2cl,"r--")
plt.plot(X,V2ss,"g-")

circle1 = plt.Circle((0,r3),w_star,color='r', fill=False)
plt.gca().add_patch(circle1)
plt.show()

####Because Mexit = 0.1 and the area ratio is outside the chart bounds
twotheta = 12
Lwstar   = 25

eta_ts = T01/(T4-T01)*((P4/P01)**((gamma-1)/gamma)-1)

print("Impeller inlet")
print(tabulate(tab1))
tab2 = []
table2('Triangle',"Crat","Wrat", "beta2*",'beta2',"Ct2","W2","Wt2","W2i","U2")
table2('Actual',Crat,Wrat,beta2_star,beta2,round(Ct2,3),round(W2,3),round(Wt2,3),round(W2,3),round(U2,3))
table2('Ideal',Crat,Wrat, beta2_star,beta2,round(Ct2i,3),round(W2i,3),round(Wt2i,3),round(W2i,3),round(U2,3))
print("Impeller exit")
print(tabulate(tab2))
print("Station 3")
print(tabulate(tab3))
print("Choke")
tab2 = []
table2('Parameter','AR','P0*','P*','T0*','T*','A*','PHIp')
table2('Value',round(AR_star,3),round(P0star,3),round(Pstar,3),round(T0star,3),round(Tstar,3),round(Astar,3),round(phip,3))
print(tabulate(tab2))
print("Exit conditions")
tab2 = []
table2('Parameter','P04','P4','T04','T4','A4','Cp','2Theta','L/w*','eta_ts')
table2('Value',round(P04,3),round(P4,3),round(T04,3),round(T4,3),round(A4,3),round(Cp,3),round(twotheta,3),round(Lwstar,3),round(eta_ts,3))
print(tabulate(tab2))
