import numpy as np
from tabulate import tabulate
from CompressibleFlowFunctions.Isentropic import *
from scipy.optimize import *


exec(open('partA.py').read())

def Asearch(A,M,mdot,Po,Rs,To,gamma):
    return delta_mass_stag(M,mdot,Po,Rs,To,gamma,A)

#######################################################
#                     Vane properties                 #
#######################################################
AR_vane    = 0.7
Chord_v    = 3.999    #centimeters
ttoc_vane  = 0.172
Zwie_vane  = np.linspace(0.7,0.8,11) #Zwiefel number
mdot_vane  = m_coola #mass flow after the vanes
M_inlet    = 0.125 #inlet Mach number
T_inlet    = T_from_Tratio(T05a,gamma2,M_inlet)
C_insound  = np.sqrt(gamma2*287*T_inlet)
C1         = C_insound*M_inlet #Inlet velocity vector
alpha1     = -10*np.pi/180

#######################################################
#                    Blade properties                 #
#######################################################
AN2        = 4.5e10/60/60*2*np.pi/(1100*12)**2 ##Non-dimensional AN2 value
AR_blade   = 1.45
Chord_b    = 1.619 #centimeters
ttoc_blade = 0.2026
Zwie_blade = np.linspace(0.85,0.95,11) #Zwiefel number
mdot_blade = m_coolb #mass flow rate at stage exit
M_exit     = np.linspace(0.3,0.45,16)
T_exit     = T_from_Tratio(T06b,gamma2,M_exit) #Static temperature exiting the stage
alpha3 = np.linspace(-5,30,36)

#######################################################
#                    Blade solution                   #
#######################################################

C_exsound  = np.sqrt(gamma2*287*T_exit) #Sound speed exiting the stage
C3         = np.multiply(M_exit,C_exsound)

A_bexit = []

HT_blade   = 1/np.sqrt(AN2/np.pi+1)
Span_blade = AR_blade*Chord_b/100 #Dividing by 100 to convert cm to m
clearance  = np.linspace(0.9,1.5,11)/100

Umax = 1100/3.281
Rtip = Span_blade/(1-HT_blade)
Rhub = Rtip*HT_blade
RPM  = Umax/Rhub*60/(2*np.pi)

for x in clearance: #Calculate a range of appropriate exit areas
    A_bexit = np.pi*(Rtip**2*(1+clearance)-Rhub**2)



U_hub = Umax
U_tip = RPM*Rtip*2*np.pi/60
U_mid = RPM*(Rhub+Span_blade/2)*2*np.pi/60

#######################################################
#                     Vane solution                   #
#######################################################
Ca_in = C1*np.cos(alpha1)
# rho05a = P05/(0.287*T05a)
# rho5a  = rho05a*(1+(gamma2-1)/2*M_inlet**2)**(-1/(gamma2-1))

# A_vane = m_htot/(rho5a*Ca_in)

Span_vane = AR_vane*Chord_v/100
Rvtip  = Rhub
Rvhub  = Rvtip - Span_vane
# Rvhub  = Rhub
# Rvtip  = Rvhub+Span_vane
A_vane = np.pi*(Rvtip**2-Rvhub**2)
M_mid  = bisect(delta_mass_stag, 0.000001,1, args=(m_coola*1000,P05*1000,287,T05,gamma2,A_vane))
T_mid  = T_from_Tratio(T05,gamma2,M_mid)
Cs_mid = np.sqrt(gamma2*287*T_mid)
Ca_mid = M_mid*Cs_mid

M_out = []

for x in A_bexit:
    try:
        M_out = np.append(M_out, bisect(delta_mass_stag, 0.000001,1, args=(m_coolb*1000,P06a*1000,287,T06b,gamma2,x)))
    except:
        pass
Cs_out = np.sqrt(gamma2*287*T06b)
Ca_out = M_out*Cs_out

C3_mid = np.sqrt(Ca_out[0]**2+(Ca_out[0]*np.tan(alpha3*np.pi/180)+U_mid)**2)
Mout_mid = C3_mid/Cs_out
print(Mout_mid)
# A_vane    = bisect(Asearch,0.000001, 999, args=(M_inlet,m_htot*1000,P05*1000,0.287,T05a,gamma2))
# Rvtip     = Span_vane - A_vane/(2*np.pi*Span_vane)

#
# Ca_in  = m_htot/(rho5a*A_vane)
# Ca_mid =
# Ca_out =
