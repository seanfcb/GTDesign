### TO RUN THIS CODE, CHANGE EXTENSION TO .PY

import numpy as np
from tabulate import tabulate

tab1 = []
tab2 = []
def table1(*args):
    tab1.append(args)
def table2(*args):
    tab2.append(args)
table1("Stage", "T0 (K)", "P0 (kPa)", "Work or Q (kJ)")


y1 = 19.4
y2 = 14.17

y = y1 + (12496.8-12000)*(y2-y1)/(14000-12000)

########################################################
#                    Input parameters                  #
########################################################
T0inf       = 273.15-56.6 ##Temperature in Kelvin
P0inf       = y ## kPa
Cinf        = 295.1 ## m/s
Minf        = 0.8 #Free stream Mach number
gamma1      = 1.4 #gamma for air
gamma2      = 1.33333 #gamma combustion products
scoop       = 0.9 #Scoop factor (Vin/Vinf)
BPR         = 5 #Bypass ratio
m_core      = 15/2.205 #15 lbm/sec to kg/s
Cpa         = 1.005
Cpg         = 1.148
eta_nozzle  = 0.93
eta_exhaust = 1.4/100

########################################################
#              Fan and compressor                      #
########################################################

eta_fan   = 0.91
eta_boost = 0.855
eta_hpc   = 0.845
PR_fan    = 1.5
PR_fanb   = 3
PR_hpc    = 3

########################################################
#                    Combustion parameters             #
########################################################
FAR      = 0.02 #Fuel air ratio
LHV      = 17200*2.326 #Heat of combustion BTU/lb to kJ/kg
eta_comb = 0.99 #Combustion efficiency
Ploss    = 1.8/100 #Combustion pressure loss

########################################################
#                    Turbine parameters                #
########################################################
eta_hpt  = 0.85
hpt_vane = 3/100
hpt_disk = 1.65/100
ITD      = 0.6/100

eta_lpt  = 0.89
########################################################
#                   Inlet conditions                   #
########################################################


Vinf = Minf*Cinf
Vin  = scoop*Vinf

#Assuming M scales with the Scoop factor
Min = scoop*Minf

P01 = P0inf*(1+eta_fan*(Vin**2)/(2*Cpa*T0inf*1000))**(gamma1/(gamma1-1))
T01 = T0inf + (Vin**2)/(2*Cpa*1000)
table1("Inlet", round(T01,2), round(P01,2), "Not applicable")




########################################################
#                   Cold stream 1-->2                  #
########################################################

m_cold = m_core*BPR
P02 = PR_fan * P01
T02 = T01 + T01/eta_fan*(PR_fan**((gamma1-1)/gamma1)-1)
W12 = m_cold*(T02-T01)#/eta_fan
table1("Fan (1-->2)", round(T02,2), round(P02,2), round(W12,2))

########################################################
#                   Hot stream  1-->3                  #
########################################################

m_hot = m_core
P03   = PR_fanb * P01
T03   = T02 + T02/(eta_boost)*((P03/P02)**((gamma1-1)/gamma1)-1)
W13   = m_hot*Cpa*(T03-T02)#/eta_boost
W23   = m_hot*(T03-T02)
table1("Low pressure compressor (2-->3)", round(T03,2), round(P03,2), round(W13,2))


########################################################
#                   Hot stream  3-->4                  #
########################################################

P04   = PR_hpc * P03
T04   = T03 + T03/(eta_hpc)*((PR_hpc)**((gamma1-1)/gamma1)-1)
W34   = m_hot*Cpa*(T04-T03)/eta_hpc
table1("High pressure compressor (3-->4)", round(T04,2), round(P04,2), round(W34,2))

########################################################
#           Hot stream  Combustion stage  4-->5a        #
########################################################

m_a    = m_cold + m_hot
m_fuel = FAR*m_hot
m_htot = m_fuel + m_hot*(1-9/100) #Removing 9% of the air flow for cooling


#T05a = (m_fuel*LHV + m_hot*Cpa*T04)/(m_htot*Cpg)
T05a = (FAR*LHV*eta_comb+Cpa*T04)/((1+FAR*eta_comb)*Cpg)
P05  = P04*(1-Ploss)
Q45  = eta_comb*LHV*m_fuel
table1("Combustion (Q, not W)  (4-->5a)", round(T05a,2), round(P05,2), round(Q45,2))

########################################################
#           Hot stream  Cooling stage  5a-->5               #
########################################################
m_ha    = 0.03*m_hot
m_coola = m_htot+m_ha
T05     = (m_htot*Cpg*T05a+m_ha*Cpa*T04)/((m_coola)*Cpg)

table1("Adding Vane cooling  (5a-->5)", round(T05,2), round(P05,2), round(Q45,2))



########################################################
#           Hot stream  HPT stage  5-->6a               #
########################################################

# T06a = T05 - (Cpa/(eta_hpt*Cpg))*(T04-T03)
WHPC = Cpa*(T04-T03)
T06a = T05 - WHPC*m_hot/(0.99*m_coola*Cpg)
P06a = P05*(1-(1/eta_hpt)*(1-T06a/T05))**(gamma2/(gamma2-1))
W56 = m_coola*Cpg*(T05-T06a)*eta_hpt
table1("High pressure turbine (5-->6a)", round(T06a,2), round(P06a,2), round(W56,2))

########################################################
#           Hot stream  HPT Disk  Cooling  6a-->6b     #
########################################################
m_hb    = 0.0265*m_hot
m_coolb = m_coola + m_hb
T06b    = (m_coola*Cpg*T06a + m_hb*Cpa*T04)/(m_coolb*Cpg)
table1("High pressure turbine disk cooling (6a-->6b)", round(T06b,2), round(P06a,2), round(W56,2))

########################################################
#           Hot stream  ITD LOSS  6b-->6               #
########################################################
T06    = T06b
P06    = P06a*(1-ITD)
table1("High pressure turbine ITD Loss (6b-->6)", round(T06,2), round(P06,2), round(W56,2))

########################################################
#           Hot stream  LPT stage  6-->7a               #
########################################################

T07a = T06 - ((W12)+(W23))/(0.99*m_coolb*Cpg)#(Cpa/(eta_lpt*Cpg))*(T03-T01)
P07  = P06*(1-(1/eta_lpt)*(1-T07a/T06))**(gamma2/(gamma2-1))
W67  = m_coolb*Cpg*(T06-T07a)*eta_lpt
table1("Low pressure turbine (6-->7a)", round(T07a,2), round(P07,2), round(W67,2))

########################################################
#           Hot stream  LPT Vane Cooling  7a-->7       #
########################################################
m_hc    = m_hot*0.021
m_coolc = m_coolb+m_hc
T07     = (m_coolb*Cpg*T07a+m_hc*Cpa*T04)/(m_coolc*Cpg)
table1("Low pressure turbine Disk Cooling (7a-->7)", round(T07,2), round(P07,2), round(W67,2))

########################################################
#           Hot stream  Exhaust  7-->8                #
########################################################
P08 = P07*(1-eta_exhaust)
T08 = T07

########################################################
#           Hot stream Nozzle      8-->9               #
########################################################
T09 = T08
P9  = P0inf
T9  = T09 * (1-eta_nozzle*(1-(P9/P08)**((gamma2-1)/gamma2)))
V9  = np.sqrt((T09-T9)*2*Cpg*1000)
Fh  = m_coolc*V9

########################################################
#                 Thrust calculations                  #
########################################################
##Cold stream
dT_cold = eta_nozzle*T02*(1-(P01/P02)**((gamma1-1)/gamma1))
C8      = np.sqrt(2*Cpa*dT_cold*1000)
Fc      = m_cold*C8
table2("Cold thrust (N)", round(Fc,2))
##Hot stream
dT_hot  = eta_nozzle*T07*(1-(P01/P07)**((gamma2-1)/gamma2))
C9      = np.sqrt(2*Cpg*dT_hot*1000)
#Fh      = m_htot*C9
table2("Hot thrust (N)", round(Fh,2))

##Total engine thrust
F_engine = Fc + Fh
table2("Total engine thrust (N)", round(F_engine,2))

########################################################
#                Specific fuel consumption             #
########################################################

SFC = FAR*m_hot*3600/F_engine
table2("Specific fuel consumption (kg/h N)", round(SFC,4))

print(tabulate(tab1))
print(tabulate(tab2))
