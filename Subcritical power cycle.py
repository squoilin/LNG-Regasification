# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:26:24 2021

@author: hayua
"""
import matplotlib.pyplot as plt

from pyswarm import pso
"Imports"
import CoolProp.CoolProp as CP
import numpy as np

from misc.utils import NoStdStreams
from components.heat_exchangers import hx
       
def model(x):
    
    P_ev=x[0]
    P_cd=x[1]
    p = np.zeros(10, dtype=float)
    T = np.zeros(10, dtype=float)
    x = np.zeros(10, dtype=float)
    h = np.zeros(10, dtype=float)
    s = np.zeros(10, dtype=float)
    
    
    # Cycle parameters:
    fluid ='R1270'
    p_low=P_cd
    p_high=P_ev
    DELTAT_sh = 5              #degreee of superheating
    DELTAT_sc = 5              #degree of subcooling
    DELTAT_pinch=10              # pinch point temperature in the condenser /K
    epsilon_s = 0.6            # efficiency of pump and turbine 
    
    # Heat sink parameters:
    T_hsin_in =273.15-160        # Natural gas inlet temperature/K -160 celsius
    
    
    # Heat source parameters:
    T_hsor_in=273.15 + 200       # waste heat inlet temperature /K  200 celsius
    cp_w = 100000                # Heat capactity flow rate of both waste heat and LNG 100kW/C

    
    #saturate liquid state at the condensation pressure:
    p[0] = p_low
    s[0] = CP.PropsSI('S','Q', 0, 'P', p_low, fluid)
    T[0] = CP.PropsSI('T','Q', 0, 'P', p_low, fluid)
    h[0] = CP.PropsSI('H','Q', 0, 'P', p_low, fluid)
    
    # pump inlet:
    p_su_cp = p_low
    T_su_cp = CP.PropsSI('T','Q', 0, 'P', p_low, fluid)-DELTAT_sc
    h_su_cp = CP.PropsSI('H','T', T_su_cp, 'P', p_su_cp, fluid)
    s_su_cp = CP.PropsSI('S','T', T_su_cp, 'P', p_su_cp, fluid)
    s[1] = s_su_cp
    T[1] = T_su_cp
    h[1] = h_su_cp
    p[1] = p_su_cp
    
    #pump outlet:
    p_ex_cp = p_high
    h_ex_cp_s = CP.PropsSI('H','S', s_su_cp, 'P', p_ex_cp, fluid)
    h_ex_cp = h_su_cp - (h_su_cp - h_ex_cp_s)/epsilon_s
    T_ex_cp = CP.PropsSI('T','H', h_ex_cp, 'P', p_ex_cp, fluid)
    s_ex_cp = CP.PropsSI('S','H', h_ex_cp, 'P', p_ex_cp, fluid)
    
    p[2] = p_high
    s[2] = s_ex_cp
    T[2] = T_ex_cp
    h[2] = h_ex_cp
    
    #From statepoint 2 to state ponit 3, we assume seawater is used to heat to state point 3 at 0 Celsius:
    p[3] = p_high
    T[3] = 273.15
    s[3] = CP.PropsSI('S','T', T[3] , 'P', p_high, fluid)
    h[3] = CP.PropsSI('H','T', T[3] , 'P', p_high, fluid)
    
    #Saturate liquid 
    p[4] = p_high
    T[4] = CP.PropsSI('T','Q', 0, 'P', p_high, fluid)
    s[4] = CP.PropsSI('S','Q', 0, 'P', p_high, fluid)
    h[4] = CP.PropsSI('H','Q', 0, 'P', p_high, fluid)
    
    #staturate vapor 
    p[5] = p_high
    T[5] = CP.PropsSI('T','Q', 1, 'P', p_high, fluid)
    s[5] = CP.PropsSI('S','Q', 1, 'P', p_high, fluid)
    h[5] = CP.PropsSI('H','Q', 1, 'P', p_high, fluid)
    
    #superheated vapor
    p[6] = p_high
    T[6] = T[5]+DELTAT_sh
    s[6] = CP.PropsSI('S','T', T[6], 'P', p_high, fluid)
    h[6] = CP.PropsSI('H','T', T[6], 'P', p_high, fluid)
    
    #isentropic turbine outlet enthalpy 
    h_ex_turb_s = CP.PropsSI('H','S', s[6], 'P', p_low, fluid)
    
    #real turbine outlet 
    p[7] = p_low
    h[7] = h[6] - (h[6] - h_ex_turb_s)*epsilon_s
    T[7] = CP.PropsSI('T','H', h[7], 'P', p_low, fluid)
    s[7] = CP.PropsSI('S','H', h[7], 'P', p_low, fluid)
    
    #staturate vapor in the condensation process 
    p[8] = p_low
    s[8] = CP.PropsSI('S','Q', 1, 'P', p_low, fluid)
    T[8] = CP.PropsSI('T','Q', 1, 'P', p_low, fluid)
    h[8] = CP.PropsSI('H','Q', 1, 'P', p_low, fluid)
    
    
    # Mass flowrate of working fluid determined by the condenser side
    M_dot1 = cp_w*(T[8]-DELTAT_pinch-T_hsin_in)/ (h[8] - h[1]) # mass flowrate calculated by LNG end
    
    # initial mass flowrate value 
    
    M_dot=M_dot1
    
    #calculate the pinch temperature for both condenser and evaporator
    # Temperature profile in the condenser and evaporator
    print("For the condenser:")
    s_cd,T_cd,T_hf,pinch_cd = hx(fluid,100,M_dot,h[7],h[1],p_low,p_low,T_hsin_in,cp_w)
    print("For the evaporator:")
    s_ev,T_ev,T_cf,pinch_ev = hx(fluid,100,M_dot,h[3],h[6],p_high,p_high,T_hsor_in,cp_w)
 
    #if the pinch temperature difference of the evaporator is less than 10, then we have to reduce the mass flowrate of the working fluid, then recalculate the temperature 
    # for bothe evaporator and condenser. However, the pinch temperature of the condenser is greater than 10 definitely if mass flowrate is decreased  
   
    while pinch_ev<10:
        
        M_dot-=0.5   
        print("For the condenser:")
        s_cd,T_cd,T_hf,pinch_cd = hx(fluid,100,M_dot,h[7],h[1],p_low,p_low,T_hsin_in,cp_w)
        print("For the evaporator:")
        s_ev,T_ev,T_cf,pinch_ev = hx(fluid,100,M_dot,h[3],h[6],p_high,p_high,T_hsor_in,cp_w)
    
    #main results
    
    W_tur=(h[6]-h[7])*M_dot
    W_pump=(h[2]-h[1])*M_dot
    W_net=W_tur-W_pump
    
    heatld_eva=M_dot*(h[6]-h[3])
    T_hsor_out=T_hsor_in-heatld_eva/cp_w
    heatld_cond=M_dot*(h[7]-h[1])
    T_hsin_out=T_hsin_in+heatld_cond/cp_w
    heatld_cond_sup=M_dot*(h[7]-h[8])
    
    D=CP.PropsSI('D','T|gas',T[7],'P',p[7],fluid)
    V=M_dot/D
    v=V/W_net
    
    print("The temperature of each state")
    print(T)
    print("The pressure of each state (Pa) ")
    print(p)
    print("The enthalpy of each state (J/kg)")
    print(h)
    print("The entropy of each state (J/kg/K)")
    print(s)
    print ("The mass flow rate of working fluid is %.2f kg/s" %(M_dot))
    print("The condensation temperature is %.2f C" %(T[8]-273.15))
    print("The heat load of evaproator is %.2f kW and the outlet temperature of heat source is %.2f C" %(heatld_eva/1000,T_hsor_out-273.15))
    print("The heat load of condenser is %.2f kW and the outlet temperature of heat sink is %.2f C" %(heatld_cond/1000,T_hsin_out-273.15))
    print("The turbine power output is %.2f kW" %(W_tur/1000))
    print("The pump power input is %.2f kW" %(W_pump/1000))
    print("The net power output is %.2f kW" %(W_net/1000))
    print("The volume flowrate of the turbine outlet is %.2f m3/s" %(V))

    print("For the condenser:")
    s_cd,T_cd,T_hf,pinch_cd = hx(fluid,200,M_dot,h[7],h[1],p_low,p_low,T_hsin_in,cp_w)
    plt.plot(T_cd)
    plt.plot(T_hf)
    plt.show()
    # # Temperature profile in the evaporator
    print("For the evaporator:")
    s_ev,T_ev,T_cf,pinch_ev = hx(fluid,200,M_dot,h[3],h[6],p_high,p_high,T_hsor_in,cp_w)
    plt.plot(T_ev)
    plt.plot(T_cf)
    plt.show()
    return -W_net,V,v

# lb = [1000000, 100000]  # lower bound for propane 
# ub = [3825000, 1500000] # upper bound for propane

lb = [1000000, 100000]  # lower bound for C3H6
ub = [4095000, 1500000] # upper bound for C3H6

# xopt, fopt = pso(model, lb, ub,swarmsize=50, maxiter=100)
#xopt= [3681703.96, 100000]  # optimal results obtained with propane as the working fluid
xopt= [3731628.34, 100000] #optimal results obtained with C3H6 as the working fluid

net_power, total_volume, specific_volume= model(xopt)
# a= CP.PropsSI('S','T', 50, 'P', 3500000, 'propane')

#%%

# from CoolProp.Plots import PropertyPlot,
# from CoolProp.Plots.SimpleCycles import StateContainer
# import pickle
# import os

# cache_plot = False
# if cache_plot:
#     filename = fluid + '.p'
#     if os.path.isfile(filename):   # Load previously saved plot
#         pp = pickle.load(open( filename, "rb" ) )
#     else:
#         states = StateContainer()
#         states_hf = StateContainer()
#         states_cf = StateContainer()
#         pp = PropertyPlot('HEOS::'+fluid, 'TS')
#         with NoStdStreams():
#             pp.calc_isolines()
#         with open(filename, 'wb') as f: 
#             pickle.dump(pp, f) 
# else:
#     states = StateContainer()
#     states_hf = StateContainer()
#     states_cf = StateContainer()
#     pp = PropertyPlot('HEOS::'+fluid, 'TS')
#     with NoStdStreams():
#         pp.calc_isolines()
        
    
# for i in range(3):
#     states[i,'T'] = T[i]
#     states[i,"S"] = s[i]
    
# for i,Tx in enumerate(T_cd):
#     states.append({'T':Tx,'S':s_cd[i]})

# for i in range(4,len(T)):
#     states.append({'T':T[i],'S':s[i]})  
# states.append({'T':T[1],'S':s[1]})    # for some reasons, the second point needs to be repeated to close the cycle
    
# for i,Tx in enumerate(T_hf):
#     states_hf.append({'T':Tx,'S':s_cd[i]})
    
# for i,Tx in enumerate(T_cf):
#     states_cf.append({'T':Tx,'S':s_ev[i]})

# with NoStdStreams():
#     pp.draw_process(states,line_opts={'color':'green'})
#     pp.draw_process(states_hf,line_opts={'color':'red', 'linestyle':'dashed'})
#     pp.draw_process(states_cf,line_opts={'color':'blue', 'linestyle':'dashed'})

# pp.show()




