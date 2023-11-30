
"""
Created on Mon Aug  7 13:20:51 2023

@author: bikashadhikari
"""

import gurobipy as gp
from gurobipy import GRB
import numpy as np


def MPC(x0, N, ref,Q, A, B, C, D, t):
    # INPUTS: 
        # x0 - initial state
        # N - prediction horizon
        # ref - reference input signal
        # A, B, C, D - state space matrix of low pass filter
        # t - total simulation time
        
    # OUTPUTS:
        # u_opt - optimal control 
        
    
    len_ref = len(t)
    u_mpc = []
    j = 0
    # MPC Uniform
    # xx = []
    last_print = 0
    while j <= len_ref:
        progress = round((j/len_ref)*100);
        if(last_print != progress) and (progress % 10 == 0):
            print("Progress: %d%%" % progress)
            last_print = progress
        
        # Initialize gurobi optimizer
        m = gp.Model("MPC")
        u = m.addMVar(N, vtype = GRB.INTEGER, name = "u", lb = np.min(Q), ub = np.max(Q))  # control variable
        x = m.addMVar((2*(N+1),1), vtype = GRB.CONTINUOUS, lb= -GRB.INFINITY, ub = GRB.INFINITY, name = "x") # state variable        
            
        # Objective function
        obj = 0                  # initialize objective function
        for i in range(0,N):
            k = 2*i
            con = u[i]-ref[j]
            e_t = C @ x[k:k+2] + D*con
            obj = obj + e_t*e_t
            
        # Constraints
        m.addConstr(x[0:2,:] == x0)  # initial condition constraint
        for i in range(0,N):
            k = 2*i
            st =x[k:k+2] 
            con = u[i]-ref[j]  # current control 
            f_value = A @ st + B * con #function value
            st_next = x[k+2: k+4]
            m.addConstr(st_next == f_value)
        
        m.update
        m.setObjective(obj, GRB.MINIMIZE) # 
        # m.params.MIPFocus= 1 #  find feasible solutions quickly
        # m.params.MIPFocus= 2 # if there is not difficulty in finding the feasible solution but we want to focus
                                # on providing the optimiality of the solutions
#        m.params.MIPFocus= 3 # if you want your run to imporve the best objective bound
#        m.tune() # perform an automated search for parameter setting that imporve performance
        
        m.Params.LogToConsole = 0 # 1: show the log in console ; 0: does not show the log in console

        m.optimize() # optimize
        allvars = m.getVars()
        values = m.getAttr("X", allvars)
        values = np.array(values) 
            
        control_u = values[0:N]  # optimal control variables
        states_x = values[N:len(values)] 

        # optimal state variavles # Store the first value of the optimal quantizer
        u_opt = [control_u[0]]  # moving horizon control
        u_mpc = np.append(u_mpc, u_opt)  # store the optimal value
        
        # apply the new optimal control to the system
        x0_new = A @ x0 + B*(u_opt - ref[j])
        #this updates states is the new initial condition for the next prediction horizon
        x0 = x0_new # set the new initial condition as the prediction using optimal control. 
        j += 1 ;
    return u_mpc

        
