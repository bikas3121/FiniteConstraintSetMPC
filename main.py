#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:25:01 2023

@author: bikashadhikari
"""
import numpy as np
from scipy  import signal
from datetime import datetime
import os
# from pathlib import Path
import csv
# import matplotlib.pyplot as plt


# custom modules
import DirectQantization as dq
import gurobiMPC as gMPC  # solve the MPC using gurobi without INL feedback
import gurobiMPCINL as gMPCINL  # solve the MPC using gurobi with INL feedback
from cfMPC import ClosedFormMPC # solve the closed from MPC problem
from signalProcessing import signalProcessing
import read_measr as rM

# %% Function to add INL to the  quantizer level
def add_INL(u_opt, Q, INL):
    u_opt = np.round(u_opt)
    INL_dict = dict(zip(Q, INL))
    u_direct_INL = []
    for i in u_opt:
        u_op_INL = 0
        INL_i = INL_dict[i]
        u_op_INL = i + INL_i
        u_direct_INL = np.append(u_direct_INL, u_op_INL)
    return u_direct_INL


# %% MPC Prediction Horiozn 
N = 2   # prediction horizon 

# %% Reference signal generation
Fs = 1e6  # sampling frequency
Ts = 1/Fs; # sampling rate
t_end = 1; # time vector duration
t = np.arange(0,t_end,Ts)  # time vector or time samples 

# %% Quantizer Parameters
l = np.longdouble(1.0) # normalized step
nob = 4  # number of bits
lvl = 2**nob  # number of quantizer levels
Q =  np.arange(0, int(lvl),1)  #quantizer levels


# %% Filter parameters 
Fc = 10e3 # cutoff frequency
Wn = Fc / (Fs / 2)
# Butterworth Filter properties 
b, a = signal.butter(2, Wn, 'low')  # tf numerator and denominator
butter = signal.dlti(*signal.butter(2, Wn,'low'))
w, h = signal.dimpulse(butter, n = 10)  #w - angular frequency, h= frequency response
h = h[0]  
A, B, C, D = signal.tf2ss(b,a) # Transfer function to StateSpace


# %% Reference Signal (signal to be recovered)
fs = 1000 # reference/ input signal frequency  , 1kHz
As = lvl  # amplitude of the signal scaled to the quantizer levels
ws = 2*np.pi*fs
t_new = np.arange(0, t_end+0.1,Ts)
ref = 1/2*((As-1)*np.sin(ws*t_new)) + (As-1)/2  # ref signal
Tc = len(t)   # truncation length
# %% Measured INL
level, INL, DNL  = rM.read_INL_measurements('Measured_INL/measurements_2023-09-10.csv', nob) 
# INL = np.round(2*INL)


# %%  Nearest neighbor(Direct)  Quantization 
u_direct = dq.DirectQuantization(Q, ref[0:Tc]) # Direct quantization with uniform quatizer levels
# INL - add INL to the quantizer levels after quantization of the reference signal i.e., replicate the behavior of the DAC. 
u_direct_INL = add_INL(u_direct, Q, INL)

# %% MPC setup - Initial conditions
x0 = [1,1]   # initial condition

# %% Moving Horizon Closed-Form Solution Parameters
mhoq_cf = ClosedFormMPC(N, A, B, C, h)
# Transformation matrices and vector quantizer levels 
# See Goodwin papaer for the references. 
Tau = mhoq_cf.Tau   # mhoq transformation matrix 1
Psi = mhoq_cf.Psi   # mhoq transformation matrix 2
Un = mhoq_cf.vectorQuantizationLevels(Q)   # vector quantization levels
Un_tilde = np.matmul(Psi, Un)   # MHOQ, transformed vector quantization levels

u_mhoq = mhoq_cf.cfQuantizer(x0, Tc, Un_tilde, ref)
u_mhoq = u_mhoq[0:Tc]   
u_mhoq_INL = add_INL(u_mhoq, Q, INL)

# %% MPC- Gurobi  without INL
x0 = np.reshape(x0,(2,1))   # reshaping into column vector
u_mpc = gMPC.MPC(x0, N, ref, Q, A, B, C, D, t) #solve the MPC problem with Gurobi. 
u_mpc = u_mpc[0:Tc]
u_mpc_INL = add_INL(u_mpc, Q, INL)  # Quantization using Non-unifrom DAC

# %% MPC Gurobi  with INL Feedback
nu_mpc_trun, nu_mpc_INL = gMPCINL.MPC(x0, N, ref, A, B, C, D, t, Q, INL)
nu_mpc_trun = nu_mpc_trun[0:Tc]
nu_mpc_INL = nu_mpc_INL[0:Tc]
nu_mpc_trun = add_INL(nu_mpc_trun, Q, INL)
ref = ref[0:Tc]

# %% Signal Processing
sgp = signalProcessing(N, t, ref, b, a)
filt_reference = sgp.referenceFilter      # Filtered reference signal#
filt_u_direct, error_direct, var_direct = sgp.signalFilter(u_direct) # Filtered Direct Quantizer by Unifrom DAC 
filt_u_direct_INL, error_direct_INL, var_direct_INL = sgp.signalFilter(u_direct_INL) # #  Direct Quantizer by Non Unifrom DAC #
filt_u_mhoq, error_mhoq, var_mhoq = sgp.signalFilter(u_mhoq)  # Ideal Uniform MPC using unifrom DAC
filt_u_mhoq_INL, error_mhoq_INL, var_mhoq_INL = sgp.signalFilter(u_mhoq_INL)  # Ideal Uniform MPC using unifrom DAC
filt_u_mpc, error_mpc, var_mpc = sgp.signalFilter(u_mpc)  # Ideal Uniform MPC using unifrom DAC
filt_u_mpc_INL, error_mpc_INL, var_mpc_INL = sgp.signalFilter(u_mpc_INL)  # Ideal Unifrom MPC using Non-unifrom DAC (added INL) 
filt_nu_mpc_INL, error_nu_mpc_INL, var_nu_mpc_INL = sgp.signalFilter(nu_mpc_INL)  # Non ideal Uniform MPC using Non-unifrom DAC
filt_nu_mpc_trun, error_nu_mpc_trun, var_nu_mpc_trun = sgp.signalFilter(nu_mpc_trun)  # Non ideal Uniform MPC truncated (i.e., before added INL)
    

# %% Write data into file
# Filename Convention: Understanding filename
# UF_Sim = unfiltered simulation data
# UF_Sim_bit_signalfrequency_cutoffrequency_samplingfrequency.csv
# frequency in kHz


# create file name with time stamp and other information
current_datetime = datetime.now().strftime("%y%h%d_%H-%M")
# Unfiltered signal
uf_folder_name = "UF-"
uf_folder_to_save_file = uf_folder_name + str(int(nob)) + str("_")+ str(int(fs))  +str("_")+ str(int(Fc)) + str("_")+ str(int(Fs))+ str("_")+current_datetime
uf_str_file_name = str(uf_folder_to_save_file)
uf_file_name =uf_str_file_name+".csv" 


# Filtered signal
f_folder_name = "F-"
f_folder_to_save_file = f_folder_name + str(int(nob)) + str("_")+ str(int(fs))  +str("_")+ str(int(Fc)) + str("_")+ str(int(Fs))+ str("_")+current_datetime
f_str_file_name = str(f_folder_to_save_file)
f_file_name =f_str_file_name+".csv" 

# File saving path
# path = r'SimulationData'
if not os.path.exists('SimulationData'):
    os.makedirs('SimulationData')
path = r'SimulationData'

# Unfiltered data
headerlist = ['Time', 'Reference', 'Direct(Uniform)', 'Direct(INL)', 'Uniform MPC', 'Uniform MPC+(INL)', 'NonUniform MPC(INL Feedback)', 'MPC Truncated']
datalist = zip(t, ref, u_direct, u_direct_INL, u_mhoq, u_mhoq_INL, u_mpc, u_mpc_INL, nu_mpc_INL, nu_mpc_trun)
with open(os.path.join(path,uf_file_name),'w') as f1:
    writer = csv.writer(f1, delimiter ='\t')
    writer.writerow(headerlist)
    writer.writerows(datalist)

# Filtered data
filt_datalist = zip(t, filt_reference, filt_u_direct, filt_u_direct_INL, filt_u_mhoq, filt_u_mhoq_INL, filt_u_mpc, filt_u_mpc_INL, filt_nu_mpc_INL, filt_nu_mpc_trun)
# write the filtered result to plot 
with open(os.path.join(path,f_file_name),'w') as f2:
    writer = csv.writer(f2, delimiter ='\t')
    writer.writerow(headerlist)
    writer.writerows(filt_datalist)



# %%
# # Plot of noise power
# leg_bar = ['D', 'D-INL','cf-MPC', 'cf-MPC-INL', 'U-MPC', 'U-MPC-INL', 'NU-MPC','NU-MPC-Trun']
# bars = [var_direct, var_direct_INL,var_mhoq, var_mhoq_INL, var_mpc, var_mpc_INL, var_nu_mpc_INL, var_nu_mpc_trun ] 
# # bars = [var_INL_direct, var_mpc1, var_mpc_INL] 
# x_pos = [0,1,2,3,4,5,6,7]  # x position of bars
# # bar plot noise power
# plt.figure()
# plt.bar(x_pos, bars, width = 0.9 )
# plt.xticks([i  for i in range(len(bars))],leg_bar[0:len(leg_bar)])
# for i in range(8):
#     plt.text(x = x_pos[i], y= bars[i] + 0.00015 , s= bars[i], size = 10)
# plt.subplots_adjust(bottom= 0.1, top = 1.5)
# plt.ylabel('Error Variance')
# plt.show()
