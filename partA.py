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
T0inf      = 273.15-56.6 ##Temperature in Kelvin
P0inf      = y ## kPa
Cinf       = 295.1 ## m/s
Minf       = 0.8 #Free stream Mach number
gamma1     = 1.4 #gamma for air
gamma2     = 1.33333 #gamma combustion products
scoop      = 0.9 #Scoop factor (Vin/Vinf)
BPR        = 5 #Bypass ratio
m_core     = 15/2.205 #15 lbm/sec to kg/s
Cpa        = 1.005
Cpg        = 1.148
eta_nozzle = 0.93

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

m_cold = m_core*BPR/(BPR+1)
P02 = PR_fan * P01
T02 = T01 + T01/eta_fan*(PR_fan**((gamma1-1)/gamma1)-1)
W12 = m_cold*Cpa*(T02-T01)/eta_fan
table1("Fan (1-->2)", round(T02,2), round(P02,2), round(W12,2))

########################################################
#                   Hot stream  1-->3                  #
########################################################

m_hot = m_core/(BPR+1)
P03   = PR_fanb * P01
T03   = T02 + T02/(eta_boost)*((PR_fanb-PR_fan)**((gamma1-1)/gamma1)-1)
W13   = m_hot*Cpa*(T03-T01)/eta_boost
table1("Low pressure compressor (fan + boost) (1-->3)", round(T03,2), round(P03,2), round(W13,2))


########################################################
#                   Hot stream  3-->4                  #
########################################################

P04   = PR_hpc * P03
T04   = T03 + T03/(eta_hpc)*((PR_hpc)**((gamma1-1)/gamma1)-1)
W34   = m_hot*Cpa*(T04-T03)/eta_hpc
table1("High pressure compressor (3-->4)", round(T04,2), round(P04,2), round(W34,2))

########################################################
#           Hot stream  Combustion stage  4-->5        #
########################################################

m_fuel = FAR*m_hot
m_htot = m_fuel + m_hot

T05 = (m_fuel*LHV + m_hot*Cpa*T04)/(m_htot*Cpg)
P05 = P04*(1-Ploss)
Q45 = eta_comb*LHV*m_fuel
table1("Combustion (Q, not W)  (4-->5)", round(T05,2), round(P05,2), round(Q45,2))

########################################################
#           Hot stream  HPT stage  5-->6               #
########################################################

T06 = T05 - (Cpa/(eta_hpt*Cpg))*(T04-T03)
P06 = P05*(1-(1/eta_hpt)*(1-T06/T05))**(gamma2/(gamma2-1))
W56 = m_htot*Cpg*(T05-T06)*eta_hpt
table1("High pressure turbine (5-->6)", round(T06,2), round(P06,2), round(W56,2))

########################################################
#           Hot stream  LPT stage  6-->7               #
########################################################

T07 = T06 - (Cpa/(eta_lpt*Cpg))*(T03-T01)
P07 = P06*(1-(1/eta_lpt)*(1-T07/T06))**(gamma2/(gamma2-1))
W67 = m_htot*Cpg*(T06-T07)*eta_lpt
table1("Low pressure turbine (6-->7)", round(T07,2), round(P07,2), round(W67,2))

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
Fh      = m_htot*C9
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
