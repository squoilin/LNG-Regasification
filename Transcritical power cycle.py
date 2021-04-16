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
    p = np.zeros(8, dtype=float)
    v = np.zeros(8, dtype=float)
    T = np.zeros(8, dtype=float)
    x = np.zeros(8, dtype=float)
    h = np.zeros(8, dtype=float)
    s = np.zeros(8, dtype=float)
    
    # Cycle parameters:
    fluid ='propane'
    # P_ev = 6000000  #pa
    # P_cd = 200000  #pa
    p_low=P_cd
    p_high=P_ev
    
    DELTAT_sh = 5
    DELTAT_sc = 5
    epsilon_s = 0.6
    
    # Heat sink parameters:
    T_hsin_in =273.15-160           # Natural gas inlet temperature /K
    DELTAT_cd=10                    # pinch point temperature of cooling water /K
    # Q_dot_cd=10000                  # heat capacity flowrates of cooling water kW/K
    
    # Heat source parameters:
    T_hsor_in=273.15 + 200        # cooling water inlet temperature /K
    DELTAT_pinch=10                    # pinch point temperature of cooling water /K
    
    # get fluid properties:
    #T_crit = CP.PropsSI("Tcrit",fluid)
    #p_crit = CP.PropsSI("Pcrit",fluid)
    T_low = CP.PropsSI('T','Q', 0.5, 'T', P_ev, fluid)
    T_high = CP.PropsSI('T','Q', 0.5, 'T', P_cd, fluid)
    #h_crit = CP.PropsSI('H','T', T_crit, 'P', p_crit, fluid)
    #s_crit = CP.PropsSI('S','T', T_crit, 'P', p_crit, fluid)
    
    #Evaporator outlet:
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
    
    p[4] = p_high
    T[4]=190+273.15
    s[4] = CP.PropsSI('S','T', T[4] , 'P', p_high, fluid)
    h[4] = CP.PropsSI('H','T', T[4] , 'P', p_high, fluid)
    
    
    p[5] = p_low
    h_ex_turb_s = CP.PropsSI('H','S', s[4], 'P', p_low, fluid)
    h[5] = h[4] - (h[4] - h_ex_turb_s)*epsilon_s
    T[5] = CP.PropsSI('T','H', h[5], 'P', p_low, fluid)
    s[5] = CP.PropsSI('S','H', h[5], 'P', p_low, fluid)
    
    
    p[6] = p_low
    s[6] = CP.PropsSI('S','Q', 1, 'P', p_low, fluid)
    T[6] = CP.PropsSI('T','Q', 1, 'P', p_low, fluid)
    h[6] = CP.PropsSI('H','Q', 1, 'P', p_low, fluid)
    
    p[7] = p[0]
    s[7] = s[0]
    T[7] = T[0]
    h[7] = h[0] 
    
    # heat sink:
    cp_w = 100000  # Heat capactity flow rate of both hot water and LNG 100kW/C
    M_dot1 = cp_w*(T[6]-DELTAT_pinch-T_hsin_in)/ (h[6] - h[1]) # mass flowrate calculated by LNG end
    

    M_dot=M_dot1
    W_tur=(h[4]-h[5])*M_dot
    W_pump=(h[2]-h[1])*M_dot
    W_net=W_tur-W_pump
    
    heatld_eva=M_dot*(h[4]-h[3])
    T_hsor_out=T_hsor_in-heatld_eva/cp_w
    heatld_cond=M_dot*(h[5]-h[1])
    T_hsin_out=T_hsin_in+heatld_cond/cp_w
    heatld_cond_sup=M_dot*(h[5]-h[6])
    
    # print("The temperature of each state")
    # print(T)
    # print("The pressure of each state (Pa) ")
    # print(p)
    # print("The enthalpy of each state (J/kg)")
    # print(h)
    # print("The entropy of each state (J/kg/K)")
    # print(s)
    # print ("The mass flow rate of working fluid is %.2f kg/s" %(M_dot))
    # print("The condensation temperature is %.2f C" %(T[6]-273.15))
    # print("The heat load of evaproator is %.2f kW and the outlet temperature of heat source is %.2f C" %(heatld_eva/1000,T_hsor_out-273.15))
    # print("The heat load of condenser is %.2f kW and the outlet temperature of heat sink is %.2f C" %(heatld_cond/1000,T_hsin_out-273.15))
    # print("The turbine power output is %.2f kW" %(W_tur/1000))
    # print("The pump power input is %.2f kW" %(W_pump/1000))
    # print("The net power output is %.2f kW" %(W_net/1000))
    # Heat source:
    # print(heatld_cond_sup)
    
    # Temperature profile in the condenser:
    # print("For the condenser:")
    s_cd,T_cd,T_hf,pinch_cd = hx(fluid,100,M_dot,h[5],h[1],p_low,p_low,T_hsin_in,cp_w)
    # plt.plot(T_cd)
    # plt.plot(T_hf)
    # plt.show()
    # # # Temperature profile in the evaporator
    # print("For the evaporator:")
    s_ev,T_ev,T_cf,pinch_ev = hx(fluid,100,M_dot,h[3],h[4],p_high,p_high,T_hsor_in,cp_w)
    # plt.plot(T_ev)
    # plt.plot(T_cf)
    # plt.show()
    
    while pinch_ev<10: # if the pinch temperature is less than 5, we just decrease the target temperature of working fluid
        
        T[4]-=2   #step size is chosen as 2 celsius
        s[4] = CP.PropsSI('S','T', T[4] , 'P', p_high, fluid)
        h[4] = CP.PropsSI('H','T', T[4] , 'P', p_high, fluid)
        h_ex_turb_s = CP.PropsSI('H','S', s[4], 'P', p_low, fluid)
        h[5] = h[4] - (h[4] - h_ex_turb_s)*epsilon_s
        T[5] = CP.PropsSI('T','H', h[5], 'P', p_low, fluid)
        s[5] = CP.PropsSI('S','H', h[5], 'P', p_low, fluid)
        
        M_dot=M_dot1
        W_tur=(h[4]-h[5])*M_dot
        W_pump=(h[2]-h[1])*M_dot
        W_net=W_tur-W_pump
        
        heatld_eva=M_dot*(h[4]-h[3])
        T_hsor_out=T_hsor_in-heatld_eva/cp_w
        heatld_cond=M_dot*(h[5]-h[1])
        T_hsin_out=T_hsin_in+heatld_cond/cp_w
        heatld_cond_sup=M_dot*(h[5]-h[6])
        s_cd,T_cd,T_hf,pinch_cd = hx(fluid,100,M_dot,h[5],h[1],p_low,p_low,T_hsin_in,cp_w)
        s_ev,T_ev,T_cf,pinch_ev = hx(fluid,100,M_dot,h[3],h[4],p_high,p_high,T_hsor_in,cp_w)
        
    print("The temperature of each state")
    print(T)
    print("The pressure of each state (Pa) ")
    print(p)
    print("The enthalpy of each state (J/kg)")
    print(h)
    print("The entropy of each state (J/kg/K)")
    print(s)
    print ("The mass flow rate of working fluid is %.2f kg/s" %(M_dot))
    print("The condensation temperature is %.2f C" %(T[6]-273.15))
    print("The heat load of evaproator is %.2f kW and the outlet temperature of heat source is %.2f C" %(heatld_eva/1000,T_hsor_out-273.15))
    print("The heat load of condenser is %.2f kW and the outlet temperature of heat sink is %.2f C" %(heatld_cond/1000,T_hsin_out-273.15))
    print("The turbine power output is %.2f kW" %(W_tur/1000))
    print("The pump power input is %.2f kW" %(W_pump/1000))
    print("The net power output is %.2f kW" %(W_net/1000))

    print("For the condenser:")
    s_cd,T_cd,T_hf,pinch_cd = hx(fluid,100,M_dot,h[5],h[1],p_low,p_low,T_hsin_in,cp_w)
    plt.plot(T_cd)
    plt.plot(T_hf)
    plt.show()
    # # Temperature profile in the evaporator
    print("For the evaporator:")
    s_ev,T_ev,T_cf,pinch_ev = hx(fluid,100,M_dot,h[3],h[4],p_high,p_high,T_hsor_in,cp_w)
    plt.plot(T_ev)
    plt.plot(T_cf)
    plt.show()
    
    D=CP.PropsSI('D','T|gas',T[5],'P',p[5],fluid)
    V=M_dot/D
    v=V/W_net
    
    from CoolProp.Plots import PropertyPlot
    from CoolProp.Plots.SimpleCycles import StateContainer
    import pickle
    import os
    
    cache_plot = False
    if cache_plot:
        filename = fluid + '.p'
        if os.path.isfile(filename):   # Load previously saved plot
            pp = pickle.load(open( filename, "rb" ) )
        else:
            states = StateContainer()
            states_hf = StateContainer()
            states_cf = StateContainer()
            pp = PropertyPlot('HEOS::'+fluid, 'HEOS::R134a','TS', tp_limits='DEF')
            pp.set_axis_limits([-0.5, 3, -200, 530])
            with NoStdStreams():
                pp.calc_isolines()
            with open(filename, 'wb') as f: 
                pickle.dump(pp, f) 
    else:
        states = StateContainer()
        states_hf = StateContainer()
        states_cf = StateContainer()
        pp = PropertyPlot('HEOS::'+fluid, 'TS')
        pp.set_axis_limits([0, 3.5, 100, 530])
        with NoStdStreams():
            pp.calc_isolines()
            
        
    for i in range(4):
        states[i,'T'] = T[i]
        states[i,"S"] = s[i]
        
    for i,Tx in enumerate(T_ev):
        states.append({'T':Tx,'S':s_ev[i]})
    
    for i in range(4,len(T)):
        states.append({'T':T[i],'S':s[i]})  
    states.append({'T':T[1],'S':s[1]})    # for some reasons, the second point needs to be repeated to close the cycle
        
    for i,Tx in enumerate(T_hf):
        states_hf.append({'T':Tx,'S':s_cd[i]})
        
    for i,Tx in enumerate(T_cf):
        states_cf.append({'T':Tx,'S':s_ev[i]})
    
    with NoStdStreams():
        pp.draw_process(states,line_opts={'color':'green'})
        pp.draw_process(states_hf,line_opts={'color':'blue', 'linestyle':'dashed'})
        pp.draw_process(states_cf,line_opts={'color':'red', 'linestyle':'dashed'})
    
    pp.show()

    
    return -W_net,V,v

lb = [8129000, 650000] #lower bound with CO2 as the working fluid
ub = [15000000,3000000]  #upper bound with CO2 as the working fluid

# lb = [4675000, 100000] #lower bound with propane as the working fluid
# ub = [10000000,1500000]  #upper bound with propane as the working fluid

# xopt, fopt = pso(model, lb, ub,swarmsize=50, maxiter=100) 
# xopt=[11569607.3265, 1161490.312] # optimum value with evaporation pressure upper bound 150 bar with CO2 as the working fluid

xopt=[5863765.83, 100000]  #optimum value with evaporation pressure upper bound 150 bar with propane as the working fluid

netpower,Volume_total, Volume_specific=model(xopt)

# xopt=[10000000, 998484] currently obtained best results. 
#%%

# from CoolProp.Plots import PropertyPlot
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




