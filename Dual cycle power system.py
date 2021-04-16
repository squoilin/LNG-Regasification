# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:26:24 2021

@author: hayua
"""
import matplotlib.pyplot as plt
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

from pyswarm import pso
"Imports"
import CoolProp.CoolProp as CP
import numpy as np

from components.heat_exchangers import hx
       
def model(x):
    
    P_ev=x[0]
    P_cd=x[1]
    P_ev_top=x[2]
    P_cd_top=x[3]
    
    p = np.zeros(10, dtype=float)
    v = np.zeros(10, dtype=float)
    T = np.zeros(10, dtype=float)
    x = np.zeros(10, dtype=float)
    h = np.zeros(10, dtype=float)
    s = np.zeros(10, dtype=float)
    
        
    p_top = np.zeros(10, dtype=float)
    v_top = np.zeros(10, dtype=float)
    T_top = np.zeros(10, dtype=float)
    x_top = np.zeros(10, dtype=float)
    h_top = np.zeros(10, dtype=float)
    s_top = np.zeros(10, dtype=float)
    

    #transfer variables to the PropsSI 
    p_low=P_cd
    p_high=P_ev
    p_low_top=P_cd_top
    p_high_top=P_ev_top


    # Cycle parameters:
    fluid ='ethane'
    fluid_top ='CO2'
    DELTAT_sh = 5 #degree of superheating 
    DELTAT_sc = 5 #degree of subcooling 
    epsilon_s = 0.6 #isentropic efficiency of pump and turbine 
    
    # Heat sink parameters:
    T_hsin_in =273.15-160         # Natural gas inlet temperature /K
    # Heat source parameters:
    T_hsor_in=273.15 + 200        # waste heat  inlet temperature /K
    DELTAT_pinch=10               # pinch point temperature of heat exchangers /K
    cp_w = 100000                 # Heat capactity flow rate of both hot water and LNG 100kW/C
    
    
    # get fluid properties:
    #saturate liquid state in condensation process:
    p[0] = p_low
    s[0] = CP.PropsSI('S','Q', 0, 'P', p_low, fluid)
    T[0] = CP.PropsSI('T','Q', 0, 'P', p_low, fluid)
    h[0] = CP.PropsSI('H','Q', 0, 'P', p_low, fluid)
    
    # pump inlet:'
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
    #staturate liquid state at evaporation pressure 
    p[4] = p_high
    T[4] = CP.PropsSI('T','Q', 0, 'P', p_high, fluid)
    s[4] = CP.PropsSI('S','Q', 0, 'P', p_high, fluid)
    h[4] = CP.PropsSI('H','Q', 0, 'P', p_high, fluid)
    
    if T[4]>273.15-DELTAT_sh :
        print("The evaporation pressure of bottom cycle is too high, please reduce the evaporation pressure")
        return
    
    #saturate gas state at the evaporation pressure   
    p[5] = p_high
    T[5] = CP.PropsSI('T','Q', 1, 'P', p_high, fluid)
    s[5] = CP.PropsSI('S','Q', 1, 'P', p_high, fluid)
    h[5] = CP.PropsSI('H','Q', 1, 'P', p_high, fluid)
    
    
      
    p[6] = p_high
    T[6] = T[3]
    s[6] = CP.PropsSI('S','T', T[6], 'P', p_high, fluid)
    h[6] = CP.PropsSI('H','T', T[6], 'P', p_high, fluid)
    
    h_ex_turb_s = CP.PropsSI('H','S', s[6], 'P', p_low, fluid)
    
    p[7]=p_low
    h[7] = h[6] - (h[6] - h_ex_turb_s)*epsilon_s
    T[7] = CP.PropsSI('T','H', h[7], 'P', p_low, fluid)
    s[7] = CP.PropsSI('S','H', h[7], 'P', p_low, fluid)
    
    p[8] = p_low
    s[8] = CP.PropsSI('S','Q', 1, 'P', p_low, fluid)
    T[8] = CP.PropsSI('T','Q', 1, 'P', p_low, fluid)
    h[8] = CP.PropsSI('H','Q', 1, 'P', p_low, fluid)
    
    # heat sink:
    M_dot1 = cp_w*(T[8]-DELTAT_pinch-T_hsin_in)/ (h[8] - h[1]) # mass flowrate calculated by LNG end
    

    M_dot=M_dot1
    W_tur=(h[6]-h[7])*M_dot
    W_pump=(h[2]-h[1])*M_dot
    W_net=W_tur-W_pump
    
    heatld_eva=M_dot*(h[6]-h[2])
    heatld_cond=M_dot*(h[7]-h[1])
    T_hsin_out=T_hsin_in+heatld_cond/cp_w
    heatld_cond_sup=M_dot*(h[7]-h[8])
    
    D=CP.PropsSI('D','T|gas',T[7],'P',p[7],fluid)
    V=M_dot/D
    v=V/W_net
    
    s_cd,T_cd,T_hf,pinch_cd = hx(fluid,100,M_dot,h[7],h[1],p_low,p_low,T_hsin_in,cp_w)
   
    #top cycle
    T_hsin_in_top= T_hsin_out  #first assiagn the outlet LNG temperature of bottom cycle to be the inlet temperature of top cycle 
    T_hsor_in_top= T_hsor_in   # the heat source inlet temperature is just the heat source inlet temperature
    p_top[0] = p_low_top
    s_top[0] = CP.PropsSI('S','Q', 0, 'P', p_low_top, fluid_top)
    T_top[0] = CP.PropsSI('T','Q', 0, 'P', p_low_top, fluid_top)
    h_top[0] = CP.PropsSI('H','Q', 0, 'P', p_low_top, fluid_top)
    
    # pump inlet:
    p_su_cp_top = p_low_top
    T_su_cp_top = CP.PropsSI('T','Q', 0, 'P', p_low_top, fluid_top)-DELTAT_sc
    h_su_cp_top = CP.PropsSI('H','T', T_su_cp_top, 'P', p_su_cp_top, fluid_top)
    s_su_cp_top = CP.PropsSI('S','T', T_su_cp_top, 'P', p_su_cp_top, fluid_top)
    s_top[1] = s_su_cp_top
    T_top[1] = T_su_cp_top
    h_top[1] = h_su_cp_top
    p_top[1] = p_su_cp_top
    
    #pump outlet:
    p_ex_cp_top = p_high_top
    h_ex_cp_s_top = CP.PropsSI('H','S', s_su_cp_top, 'P', p_ex_cp_top, fluid_top)
    h_ex_cp_top = h_su_cp_top - (h_su_cp_top - h_ex_cp_s_top)/epsilon_s
    T_ex_cp_top = CP.PropsSI('T','H', h_ex_cp_top, 'P', p_ex_cp_top, fluid_top)
    s_ex_cp_top = CP.PropsSI('S','H', h_ex_cp_top, 'P', p_ex_cp_top, fluid_top)
    
    p_top[2] = p_high_top
    s_top[2] = s_ex_cp_top
    T_top[2] = T_ex_cp_top
    h_top[2] = h_ex_cp_top
    
    #From statepoint 2 to state ponit 3, we assume seawater is used to heat to state point 3 at 0 Celsius:
    p_top[3] = p_high_top
    T_top[3] = 273.15
    s_top[3] = CP.PropsSI('S','T', T_top[3] , 'P', p_high_top, fluid_top)
    h_top[3] = CP.PropsSI('H','T', T_top[3] , 'P', p_high_top, fluid_top)
    
    
    #assume that waste heat can heat the working fluid to the maximum temperature 190 K.
    p_top[4] = p_high_top
    T_top[4]=190+273.15
    s_top[4] = CP.PropsSI('S','T', T_top[4] , 'P', p_high_top, fluid_top)
    h_top[4] = CP.PropsSI('H','T', T_top[4] , 'P', p_high_top, fluid_top)
    
    p_top[5] = p_low_top
    h_ex_turb_s_top = CP.PropsSI('H','S', s_top[4], 'P', p_low_top, fluid_top)
    h_top[5] = h_top[4] - (h_top[4] - h_ex_turb_s_top)*epsilon_s
    T_top[5] = CP.PropsSI('T','H', h_top[5], 'P', p_low_top, fluid_top)
    s_top[5] = CP.PropsSI('S','H', h_top[5], 'P', p_low_top, fluid_top)
    
    
    p_top[6] = p_low_top
    s_top[6] = CP.PropsSI('S','Q', 1, 'P', p_low_top, fluid_top)
    T_top[6] = CP.PropsSI('T','Q', 1, 'P', p_low_top, fluid_top)
    h_top[6] = CP.PropsSI('H','Q', 1, 'P', p_low_top, fluid_top)
    
    M_dot_top = cp_w*(T_top[6]-DELTAT_pinch-T_hsin_in_top)/ (h_top[6]-h_top[1]) # mass flowrate calculated by LNG end
    

    s_ev_top,T_ev_top,T_cf_top,pinch_ev_top = hx(fluid_top,100,M_dot_top,h_top[3],h_top[4],p_high_top,p_high_top,T_hsor_in_top,cp_w)
    
    s_cd_top,T_cd_top,T_hf_top,pinch_cd_top = hx(fluid_top,100,M_dot_top,h_top[5],h_top[1],p_low_top,p_low_top,T_hsin_in_top,cp_w)


    while pinch_ev_top<10: # if the pinch temperature is less than 5, we just decrease the target temperature of working fluid
        
        T_top[4]-= 2
        p_top[4] = p_high_top
        s_top[4] = CP.PropsSI('S','T', T_top[4] , 'P', p_high_top, fluid_top)
        h_top[4] = CP.PropsSI('H','T', T_top[4] , 'P', p_high_top, fluid_top)
        
        
        p_top[5] = p_low_top
        h_ex_turb_s_top = CP.PropsSI('H','S', s_top[4], 'P', p_low_top, fluid_top)
        h_top[5] = h_top[4] - (h_top[4] - h_ex_turb_s_top)*epsilon_s
        T_top[5] = CP.PropsSI('T','H', h_top[5], 'P', p_low_top, fluid_top)
        s_top[5] = CP.PropsSI('S','H', h_top[5], 'P', p_low_top, fluid_top)
        
        
        M_dot_top = cp_w*(T_top[6]-DELTAT_pinch-T_hsin_in_top)/ (h_top[6] - h_top[1])
        
        s_ev_top,T_ev_top,T_cf_top,pinch_ev_top = hx(fluid_top,100,M_dot_top,h_top[3],h_top[4],p_high_top,p_high_top,T_hsor_in_top,cp_w)
    
        s_cd_top,T_cd_top,T_hf_top,pinch_cd_top = hx(fluid_top,100,M_dot_top,h_top[5],h_top[1],p_low_top,p_low_top,T_hsin_in_top,cp_w)
        
        
    W_tur_top=(h_top[4]-h_top[5])*M_dot_top
    W_pump_top=(h_top[2]-h_top[1])*M_dot_top
    W_net_top=W_tur_top-W_pump_top
    
    heatld_eva_top=M_dot_top*(h[4]-h[3])
    T_hsor_out_top=T_hsor_in_top-heatld_eva_top/cp_w
    heatld_cond_top=M_dot_top*(h_top[5]-h_top[1])
    T_hsin_out_top=T_hsin_in_top+heatld_cond_top/cp_w
    heatld_cond_sup_top=M_dot_top*(h_top[6]-h_top[5])
    
    
    D_top=CP.PropsSI('D','T|gas',T_top[5],'P',p_top[5],fluid_top)
    V_top=M_dot_top/D_top
    v_top=V_top/W_net_top
    
    #overall performance 
    W_net_overall=W_net+W_net_top
    

    print("The temperature of each state")
    print(T)
    print("The pressure of each state (Pa) ")
    print(p)
    print("The enthalpy of each state (J/kg)")
    print(h)
    print("The entropy of each state (J/kg/K)")
    print("For the bottom cycle")
    print ("The mass flow rate of working fluid is %.2f kg/s" %(M_dot))
    print("The condensation temperature of is %.2f C" %(T[8]-273.15))
    print("The heat load of evaproator is %.2f kW " %(heatld_eva/1000))
    print("The heat load of condenser is %.2f kW and the outlet temperature of heat sink is %.2f C" %(heatld_cond/1000,T_hsin_out-273.15))
    print("The turbine power output is %.2f kW" %(W_tur/1000))
    print("The pump power input is %.2f kW" %(W_pump/1000))
    print("The net power output is %.2f kW" %(W_net/1000))
    
    
    print("For the top cycle")
    print ("The mass flow rate of working fluid is %.2f kg/s" %(M_dot_top))
    print("The condensation temperature is %.2f C" %(T_top[6]-273.15))
    print("The heat load of evaproator is %.2f kW and the outlet temperature of heat source is %.2f C" %(heatld_eva_top/1000,T_hsor_out_top-273.15))
    print("The heat load of condenser is %.2f kW and the outlet temperature of heat sink is %.2f C" %(heatld_cond_top/1000,T_hsin_out_top-273.15))
    print("The turbine power output is %.2f kW" %(W_tur_top/1000))
    print("The pump power input is %.2f kW" %(W_pump_top/1000))
    print("The net power output is %.2f kW" %(W_net_top/1000))
    
    
    print("For the bottom cycle condenser:")
    s_cd,T_cd,T_hf,pinch_cd = hx(fluid,200,M_dot,h[7],h[1],p_low,p_low,T_hsin_in,cp_w)
    plt.plot(T_cd)
    plt.plot(T_hf)
    plt.show()
    
    print("For the top cycle condenser:")
    s_cd_top,T_cd_top,T_hf_top,pinch_cd_top = hx(fluid_top,100,M_dot_top,h_top[5],h_top[1],p_low_top,p_low_top,T_hsin_in_top,cp_w)
    plt.plot(T_cd_top)
    plt.plot(T_hf_top)
    plt.show()
    # # Temperature profile in the evaporator
    print("For the top cycle evaporator:")
    s_ev_top,T_ev_top,T_cf_top,pinch_ev_top = hx(fluid_top,100,M_dot_top,h_top[3],h_top[4],p_high_top,p_high_top,T_hsor_in_top,cp_w)
    plt.plot(T_ev_top)
    plt.plot(T_cf_top)
    plt.show()
    return -W_net_overall, V, v, V_top, v_top 

lb = [1000000, 100000, 8000000, 2000000]
ub = [1600000, 500000, 15000000,5000000]

# xopt, fopt = pso(model, lb, ub,swarmsize=50, maxiter=100) # activate this line when we perform optimization 

# xoptfound= [2000000, 100000, 10000000,2530130]  #upper bound is 100 bar with ethane as the working fluid
xoptfound= [2000000, 100000, 14801194,4048317]    #upper bound is 150 bar with ethane as the working fluid

# xoptfound= [1600000, 100000, 13823939,2767688]    #upper bound is 150 bar with ethylene as the working fluid
# xoptfound= [1600000, 100000, 15000000,5000000]  #upper bound is 150 bar with R116 as the working fluid


netpower,total_volume, spefic_volume, total_volume_top, spefic_volume_top = model(xoptfound)

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




