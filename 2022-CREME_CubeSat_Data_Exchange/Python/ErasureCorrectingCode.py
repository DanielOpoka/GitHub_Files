#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 12:44:18 2022

@author: d.opoka
"""

"""
##############################################################################
                    TELECOMMUNICATIONS PROTOCOL
##############################################################################
"""

import numpy as np  
import matplotlib.pyplot as plt  
from matplotlib.pyplot import figure
from datetime import datetime 
from crc import CrcCalculator, Configuration
import random
import zfec
import yaml
from yaml.loader import SafeLoader
from matplotlib.backends.backend_pdf import PdfPages

# %%

"""
Transfer Frames are data structures that contain many Space Packets,
When the model is in the Bad state, the Transfer Frame will lose a Space Packet. 

Binary Erasure Channel:
This is an abstract model that uses an erasure probability (epsilon = p) 
calculated from the astrodynamics of SimulationOfTumbling.

2-State Markov Chain using the Gilbert Elliot Model:
It uses 2 states:
    1 - Good state (G): packet is received
    2 - Bad state (B): packet is lost
It is made up of 2 probabilities, p (from G to B) and r (from B to G).
"""

'-------------------------------- FUNCTIONS ----------------------------------'
# Generator for Space Packets with a defined length with Random Content in bytearrays from hexadecimal
def SP_generator(len_SP_cont):
    # Primary Header
    header_PSP = bytearray.fromhex("080ac40b002d1003190002")
    
    space_packet = header_PSP
    # Secondary Header - Adding CUC time
    now = datetime.utcnow()
    TAI_ref = datetime(1958,1,1)
    duration = now - TAI_ref
    elapsed_seconds = int(duration.total_seconds())-1800 # artificially decreasing the packet time for tests reasons
    CUC_time = elapsed_seconds.to_bytes(4, 'big')
    header_SSP = CUC_time
    space_packet += header_SSP
    
    # Content of XX bytes
    cont_SP = "%00x" % random.randrange(16**(len_SP_cont * 2))
    while (len(cont_SP) != (len_SP_cont * 2)):
        cont_SP = "0" + cont_SP
    
    space_packet += bytearray.fromhex(cont_SP)
    
    # Add CCSDS CRC
    checksum = crc_calculator.calculate_checksum(space_packet)
    checksum_bytes = checksum.to_bytes(2, 'big')
    space_packet += checksum_bytes
    
    return space_packet

# Transfer Frame with Random Content
def TM_frame_generator(len_TF, nb_SP, len_SP_cont):    
    # Transfer Frame = Header + n Space Packets + (Idle Packet)
    # TM_frame_list = []      
    # Transfer Frame Header
    header_TMTF = bytearray.fromhex("3fc1e30b1800")
    len_header_TMTF = len(header_TMTF)

    # TM_frame_list.append(header_TMTF)
    TM_frame = header_TMTF
    
    # Space Packet(s) with Random Content:
    # Space Packet = Primary Header + Secondary Header + Content + CRC
    len_SP = 0
    for ind in range(nb_SP):
        space_packet = SP_generator(len_SP_cont)
        
        # TM_frame_list.append(space_packet)
        TM_frame += space_packet
    
        len_SP += len(space_packet)
    
    # Idle Packet Header
    header_IP = bytearray.fromhex("07ffffff0416")
    # Idle Packet Content
    idle_packet = header_IP
    idle_packet_end = bytearray.fromhex("0100c000")
    
    len_IP = len_TF - (len_header_TMTF + len_SP + len(header_IP) + len(idle_packet_end))
    
    cont_IP = ""
    for ind in range(len_IP * 2): 
        cont_IP += "a"
    
    idle_packet += bytearray.fromhex(cont_IP)
    idle_packet += idle_packet_end
    
    # TM_frame_list.append(idle_packet)
    TM_frame += idle_packet
  
    return TM_frame


# Binary Erasure Channel with dictionaries
def BEC(eps, data_sent, bool_txt):
    
    if type(data_sent) == list:
        data_rec = []                                 # Received Data Unit (list)
    elif type(data_sent) == dict:
        data_rec = {}                                 # Received Data Unit (dictionary)
        
    state = (random.random() > eps)                   # State [True / False]
    data_sent_state = []                              # Sent Data Unit [bits] 
    
    for elem in data_sent:
        data_sent_state.append(state)
        
        if state == True:
                if type(data_sent) == list:
                    data_rec.append(elem)
                elif type(data_sent) == dict:
                    data_rec[elem] = data_sent[elem]
                    
        state = (random.random() > eps)

    # Number of Received Data Units [#]
    nb_data_rec = np.count_nonzero(data_sent_state)  
    
    # Theoretical packet loss [%]
    packets_loss_the = eps
    # Actual packet loss [%]
    packets_loss_act = (1 - (nb_data_rec / len(data_sent)))
    # Percentage error between actual and theoretical packet loss [%]
    per_error = ((abs(packets_loss_the - packets_loss_act) / packets_loss_the) * 100)    
    
    if bool_txt == True:
        print("BINARY ERASURE CHANNEL")
        print("packet_sent = " + str(len(data_sent)))
        print('packet_rec = ' + str(nb_data_rec) + '\n')  
    
        print("packet_loss_the = " + str(packets_loss_the))
        print("packet_loss_act = " + str(packets_loss_act))    
        print('per_error = ' + str(per_error) + '% \n')  

    return data_rec, packets_loss_the, packets_loss_act, per_error

# zfec encoding + decoding sent through the BEC
def zfec_BEC(nb_k, nb_n, eps, data_k, bool_txt):
    # nb_k ->  number of useful data units [#]
    # nb_n -> number of total data units (useful + redundant) [#]
    # eps -> erasure probability []
    # data_k -> useful data [bits]
    # bool_txt -> to show text
    
    enc = zfec.Encoder(nb_k, nb_n)
    ind_n = list(range(nb_k, nb_n))
    data_n = enc.encode(data_k, ind_n)
    data_n = data_k + data_n
    
    data_n_dict = {}
    for ind, elem in enumerate(data_n):
        data_n_dict[ind] = elem
    data_n_rec, plt_zfc, pla_zfc, per_zfc = BEC(eps, data_n_dict, bool_txt)
    
    if ((len(data_n_rec) < nb_k) == False):   
        while (len(data_n_rec) != nb_k):
            d = data_n_rec
            (k := next(iter(d)), d.pop(k))
        
        dec = zfec.Decoder(nb_k, nb_n)
        data_n_dec = dec.decode(list(data_n_rec.values()), list(data_n_rec.keys()))
        
        if (data_n_dec == data_k):
            pla_zfc = 0
        
        data_rec = data_n_dec
    else:
        packets_rec = 0
        data_rec = list(data_n_rec.values())
        for ind0 in range(len(data_n_rec)):
            for ind1 in range(len(SP_k)):
               if (data_rec[ind0] == data_k[ind1]): 
                   packets_rec += 1
    
        # Actual packet loss [%]
        pla_zfc = (1 - (packets_rec / len(SP_k)))
        # Percentage error between actual and theoretical packet loss [%]
        per_zfc = ((abs(plt_zfc - pla_zfc) / plt_zfc) * 100)
        
        data_rec = data_n_rec
    
    if (bool_txt == True):
        print("BINARY ERASURE CHANNEL")
        if ((len(data_n_rec) < nb_SP_k) == False):
            print('*** Message was correctly decoded! :) ***')
            print("packet_loss_act = " + str(pla_zfc)) 
        else:
            print('*** Decoding cannot be done :( *** \n')
            
            print('However... some correct Space Packets were received:')
            print("packet_k_sent = " + str(len(SP_k)))
            print('packet_n_rec = ' + str(packets_rec) + '\n')  
        
            print("packet_loss_the = " + str(plt_zfc))
            print("packet_loss_act = " + str(pla_zfc))    
            print('per_error = ' + str(per_zfc) + '% \n')  
            
    return data_rec, plt_zfc, pla_zfc, per_zfc

# 2-State Markov Chain
def MC2S(p, r, TF_list, bool_txt):
    'TF must be a list'
    TF_rec_list = []                               # Received Transfer Frames  
    state = True                              # State [True / False]
    packets = []                              # Transmitted Packets [] 
    
    for elem in TF_list:
        packets.append(state)
        
        if state == True:
            state = (random.random() > p)
            TF_rec_list.append(elem)
        elif state == False:
            state = (random.random() > (1 - r))
      
    packets_rec = np.count_nonzero(packets)  
    
    # Theoretical packet loss [%]
    packets_loss_the = 1 - (r / (p + r))  
    # Actual packet loss [%]
    packets_loss_act = (1 - (packets_rec / len(TF_list)))
    # Percentage error between actual and theoretical packet loss [%]
    per_error = ((abs(packets_loss_the - packets_loss_act) / packets_loss_the) * 100)    
    
    if bool_txt == True:
        print("2-STATE MARKOV CHAIN")
        print("packet_sent = " + str(len(TF_list)))
        print('packet_rec = ' + str(packets_rec) + '\n')  
    
        print("packet_loss_the = " + str(packets_loss_the))
        print("packet_loss_act = " + str(packets_loss_act))    
        print('per_error = ' + str(per_error) + '% \n')  

    return TF_rec_list, packets_loss_the, packets_loss_act, per_error

# For determining the Antenna Intermittency from Received Power
def det_ant_int(ant_clr_dwl_PRX, lag_time):
    ant_int = np.zeros(sim_tele)            # Antenna Intermittence [-1 / 0 / 1]
    
    # Time Left for Transmission [s]
    if (ant_clr_dwl_PRX[0] == 0):
            trm_time = 0
    else:
        trm_time = lag_time        
             
    for ind_tele, elem_tele in enumerate(ant_clr_dwl_PRX):
        if (elem_tele > 0):
            ant_int[ind_tele] += 1
            trm_time = lag_time 
    
        if (elem_tele == 0):
            trm_time -= sim_tstp
    
        if (elem_tele == 0) and (trm_time >= 0):
            ant_int[ind_tele] -= 1
    
    return ant_int

def calc_R_eps(ant_int):
    # Total of time for communication [s]
    time_tot = np.count_nonzero(ant_int) * sim_tstp
    # Time of connection for communication [s]
    time_con = np.sum(ant_int) * sim_tstp         
    
    R_est = (time_con / time_tot)
    eps_est = 1 - R_est
    
    return R_est, eps_est

# Calculating the mean Fading Time as per the Antenna Intermittence's Lag Time
def calc_fading(ant_int, lag_time):
    count_fading = 0
    bool_fading = False
    
    ind_fading = []
    time_fading = []
    
    for ind_fad, elem_fad in enumerate(ant_int):
        
        if (elem_fad == -1):
            count_fading += 1 * sim_tstp
        else:
            bool_fading = True
            
        if (bool_fading == True):
            if (count_fading >= (lag_time - sim_tstp)):
                count_fading = 0
            ind_fading.append(ind_fad)    
            time_fading.append(count_fading)
            count_fading = 0
            bool_fading = False
            
        
    ind_fading = np.array(ind_fading)
    time_fading = np.array(time_fading)
    
    time_fading_mean = np.sum(time_fading) / np.count_nonzero(time_fading)    
          
    return ind_fading, time_fading, time_fading_mean

def closest_frac(ref):
    dem = 20
    num = dem - 1
    frac = num / dem
    
    while (frac >= ref):
        dem -= 1
        num -= 1
        
        frac = num / dem
        
    dem += 1
    num += 1
        
    return num, dem

# Determine index/position of minimum value
def min_ind(est, real):
    ind = np.abs(est - real)
    ind = np.where(ind == np.amin(ind))
    
    return ind

def subplot_labels(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


# %%

# For controlling the experiment's "randomness"
random.seed(0)

'--------------------------------- YAML INPUT FILE ---------------------------'
with open('data/input.yaml') as f:
    data = yaml.load(f, Loader=SafeLoader)


# %%
'---------------------------- DOWNLINK PARAMETERS ----------------------------'
sim_tini = data['SIM']['int_time']                 # Beginning of time [s]
sim_tend = data['SIM']['end_time']                 # Ending of time [s]
sim_tstp = data['SIM']['time_stp']                 # Time Step [s]
# Number of elements for time [#]
sim_tele = round((sim_tend - sim_tini) / sim_tstp)              
# Time [s]
sim_time = np.linspace(sim_tini, sim_tend, sim_tele)

orb_per = data['ORBIT']['per'] * 60                # Orbit period/duration [s]
ind_1op = min_ind(sim_time, orb_per)               # Index of 1 rrbit period [#]
nb_orb = round(sim_tend / orb_per)                 # Number of orbits per day [#] 

dwl_spd = 85 * 1000                                # Bit-rate [bps]
b_daily = 67.6 * (10 ** 6)                         # Daily Download length [bits]
dwl_time = b_daily / dwl_spd                       # Time to download daily data [s]    

'Parameters for Creation of Transfer Frame with randomized Space Packet'
len_TF = 1115                                      # Length of Transfer Frame [bytes / octets] 
len_header_TMTF = 6                                # Length of TM Transfer Frame Header [bytes / octets]
len_IP_min = 12                                    # Minimum size of an Idle Packet [bytes / octets]
len_SP = 64                                        # Length of SP [bytes / octets]                            
len_SP_cont = len_SP - 11 - 4 - 2                  # Length of SP content [bytes / octets]    

width = 16
poly=0x1021
init_value=0xffff
final_xor_value=0x00
reverse_input=False
reverse_output=False

configuration = Configuration(width, poly, init_value, final_xor_value, reverse_input, reverse_output)  
use_table = True
crc_calculator = CrcCalculator(configuration, use_table)

time_TF = (len_TF * 8) / dwl_spd                   # Time to download 1 Transfer Frame [s]

time_con_min = 16 * 60                             # Minimum connection time [s]
time_con_exp = 73                                  # Experimental time of connection [s]
time_dis_exp = 17                                  # Experimental time of disconnection [s]
time_con_max = 25 * 60                             # Maximum connection time [s]

R_real_min = dwl_time / time_con_max               # Minimum experimental coding rate []
R_real_max = dwl_time / time_con_min               # Maximum experimental coding rate []

eps_real_min = 1 - R_real_max                      # Minimum experimental erasure probability []               
eps_real_max = 1 - R_real_min                      # Maximum experimental erasure probability []

R_real_ran = np.array([(2 / 3), (3 / 4), (5 / 6), (7 / 8)])
R_real_num = np.array([2, 3, 5, 7])
R_real_dem = np.array([3, 4, 6, 8])
eps_real_ran = 1 - R_real_ran

# Importing communication intermettincy data for modelling: 
# Cleared Antenna Downlink Received Power 
# data_str = 'data/ant_clr_gd=15-t_end=86400.0s-t_stp=0.1s.dat'
data_str = 'data/ant_clr_avg-t_end=86400.0s-t_stp=0.1s.dat'
ant_clr = np.loadtxt(data_str)
ant_clr_vis = ant_clr[0]
ant_clr_dwl_PRX = ant_clr[1]

R_vis = np.zeros(nb_orb)
time_con = np.zeros(nb_orb)
time_vis = np.zeros(nb_orb)

for ind_orb in range(nb_orb):
    ind_ini = ind_orb * ind_1op[0][0]
    ind_end = (ind_orb + 1) * ind_1op[0][0]
    
    time_vis[ind_orb] = np.count_nonzero(ant_clr_vis[ind_ini:ind_end]) * sim_tstp
    time_con[ind_orb] = np.count_nonzero(ant_clr_dwl_PRX[ind_ini:ind_end]) * sim_tstp

    if (time_con[ind_orb] > 0):
        R_vis[ind_orb] = time_con[ind_orb] / time_vis[ind_orb]

time_con_tot = np.sum(time_con)

# Mean Simulated Coding Rate
R_vis_mean = np.sum(R_vis) / np.count_nonzero(R_vis)
# Mean Simulated Erasure Probability
eps_vis_mean = 1 - R_vis_mean


# %%
'--------------------------- DATA SIZE STUDY ------------------------'
# Studying all relevant information about the size of the data as per the size
# of the Space Packet.
len_SP_ran = np.array([32, 64, 128])               # Range of SP lengths [bytes]
time_SP_ran = (len_SP_ran * 8) / dwl_spd           # Range of Time to send SP [s]

# X -> Range of Lengths of SP, Y -> Range of Coding Rates
# Number of Space Packets in a Transfer as per the size of the Space Packet [#]
nb_SP_TF = np.zeros([len_SP_ran.shape[0], R_real_ran.shape[0] + 1]) 
# Number of Transfer Frames as per the size of the Space Packet [#]
nb_TF_SP = np.zeros([len_SP_ran.shape[0], R_real_ran.shape[0] + 1]) 
# Daily Data sizeas per the size of the Space Packet [bits]
b_SP_daily = np.zeros([len_SP_ran.shape[0], R_real_ran.shape[0] + 1])

for ind_len_SP, elem_len_SP in enumerate(len_SP_ran):
    # Number of Space Packets per Transfer Frame [#]
    nb_SP = int(np.floor((len_TF - len_header_TMTF - len_IP_min) / elem_len_SP))   
    # Number of Total Space Packets for b_daily [#]
    nb_SP_tot = np.ceil(b_daily / (elem_len_SP * 8))
    # Daily number of bits in Space Packets [bits]
    nb_SP_TF[ind_len_SP][-1] = nb_SP
    nb_TF_SP[ind_len_SP][-1] = round(nb_SP_tot / nb_SP)
    b_SP_daily[ind_len_SP][-1] = nb_TF_SP[ind_len_SP][-1] * (len_TF * 8)

    for ind_R, R_real in enumerate(R_real_ran):
        # Number of total Space Packets (coded) [#]                     
        nb_SP_n = np.floor(nb_SP / R_real_dem[ind_R]) * R_real_dem[ind_R]
        nb_SP_n = int(nb_SP_n)
        # Number of useful Space Packets (coded) [#]  
        nb_SP_k = int(nb_SP_n * R_real)
        
        nb_SP_TF[ind_len_SP][ind_R] = nb_SP_k
        nb_TF_SP[ind_len_SP][ind_R]  = np.ceil(nb_SP_tot / nb_SP_k)
        b_SP_daily[ind_len_SP][ind_R] = nb_TF_SP[ind_len_SP][ind_R] * (len_TF * 8)                                                          
        

# %%
'-------------------------- BEC PARAMETERS CALCULATIONS ----------------------'
# Coding Rate (k / n); k = useful, n = useful + redundant 
# Theretical real coding rate []
R_real_exp = time_con_exp / (time_con_exp + time_dis_exp)
# Theretical real  []
eps_real_exp = time_dis_exp / (time_con_exp + time_dis_exp)

ref = R_real_exp
num, dem = closest_frac(ref)
R_real = num / dem

lag_time = 75                            # Final lag time of satellite comm. [s]

bool_txt = False

ant_int_lag = det_ant_int(ant_clr_dwl_PRX, lag_time)
ind_fading, time_fading, time_fading_mean = calc_fading(ant_int_lag, lag_time)
ant_int_fad = det_ant_int(ant_clr_dwl_PRX, time_fading_mean)

R_est_fad, eps_est_fad = calc_R_eps(ant_int_fad)


# %%
'----------------------------- zfec OF SPACE PACKETS -------------------------'
# If the fading times are in the order of a few milliseconds,
# encoding the Space Packets works best.

# Fading Time for Space Packets [s]
time_fading_SP = time_fading_mean / time_SP_ran

# Number of Space Packets per Transfer Frame (uncoded) [#]
nb_SP = int(np.floor((len_TF - len_header_TMTF - len_IP_min) / len_SP))     

'For Binary Erasure Channel (BEC_list)'
# For Simulation
nb_TF = 1000
nb_elem = 100
sim_eps_SP = np.linspace(0, 0.5, nb_elem)

# Number of coding rate [#]
nb_R = R_real_ran.shape[0]

pla_zfc_BEC_SP_stat = np.zeros([nb_R, nb_elem])        # Statistics of Actual Packet Loss through BEC (zfec-encoded) 
pla_zfc_BEC_SP = np.zeros(nb_TF)                       # Actual Packet Loss through BEC (zfec-encoded)

for ind_R, R_real in enumerate(R_real_ran):
    # Number of total Space Packets (coded) [#]                     
    nb_SP_n = np.floor(nb_SP / R_real_dem[ind_R]) * R_real_dem[ind_R]
    nb_SP_n = int(nb_SP_n)
    # Number of useful Space Packets (coded) [#]  
    nb_SP_k = int(nb_SP_n * R_real)

    for ind_eps, eps in enumerate(sim_eps_SP):
        for ind_TF in range(nb_TF):
            bool_txt = False
            if (bool_txt == True):
                print("----------------- nb_TF = " + str(ind_TF) + " -----------------")    
            
            # Generated Min. SP for zfec encoding
            SP_k = []
            for ind in range(nb_SP_k):
                space_packet = SP_generator(len_SP_cont)
                SP_k.append(space_packet)
             
            # plt_zfc = Theoretical % of Packet Loss of zfec-encoded Stream
            # pla_zfc = Actual % of Packet Loss of zfec-encoded Stream
            # per_zfc = % Error between Theoretical (true value) and Actual (estimated value)
            TF_rec, plt_zfc, pla_zfc, per_zfc = zfec_BEC(nb_SP_k, nb_SP_n, eps, SP_k, bool_txt)
            
            # TFs_rec_BEC_list_list.append(TF_rec_list) 
            pla_zfc_BEC_SP[ind_TF] = pla_zfc

        pla_zfc_BEC_SP_stat[ind_R][ind_eps] = np.mean(pla_zfc_BEC_SP)
'''
# Plotting Parameters
plt.close('all')
nb_fig = 0
tit_fontsize = 18
txt_fontsize = 14
figsize_X = 8
figsize_Y = 6
dpi_ind = 100

figure()
'''
# %%
'--------------------------- zfec OF TRANSFER FRAMES -------------------------'
# If the fading times are in the order of a few seconds,
# encoding the Transfer Frames works best.

bool_txt = False

# Number Faded (Erased) Transfer Frames [#]
nb_fad_TF = int(np.ceil(time_fading_mean / time_TF)) 
# Number of Useful Transfer Frames to mitigate the Faded Transfer Frames [#]
nb_TF_k = np.zeros(R_real_ran.shape[0])
# Number of Total Transfer Frames to mitigate the Faded Transfer Frames [#]
nb_TF_n = np.zeros(R_real_ran.shape[0])

for ind_R, R_real in enumerate(R_real_ran):
    nb_TF_k[ind_R] = int(nb_fad_TF * R_real_num[ind_R]) 
    nb_TF_n[ind_R] = nb_TF_k[ind_R] + nb_fad_TF

'For Binary Erasure Channel (BEC_list)'
# For Simulation
nb_sim = 100
nb_elem = 100
sim_eps_TF = np.linspace(0, 0.5, nb_elem)

pla_zfc_BEC_TF_stat = np.zeros([nb_R, nb_elem])        # Statistics of Actual Packet Loss through BEC (zfec-encoded) 
pla_zfc_BEC_TF = np.zeros(nb_sim)                       # Actual Packet Loss through BEC (zfec-encoded)

nb_TF_k_zfc = int(nb_TF_k[0])
nb_TF_n_zfc = int(nb_TF_n[0])

for ind_eps, eps in enumerate(sim_eps_TF):
    for ind_sim in range(nb_sim):
        bool_txt = False
        if (bool_txt == True):
            print("----------------- nb_sim = " + str(ind_sim) + " -----------------")    
        
        # Generated Min. SP for zfec encoding
        TM_frame_k = []
        for ind_TF_zfc in range(nb_TF_k_zfc):
            TM_frame_k.append(TM_frame_generator(len_TF, nb_SP, len_SP_cont))
         
        # plt_zfc = Theoretical % of Packet Loss of zfec-encoded Stream
        # pla_zfc = Actual % of Packet Loss of zfec-encoded Stream
        # per_zfc = % Error between Theoretical (true value) and Actual (estimated value)
        TF_rec, plt_zfc, pla_zfc, per_zfc = zfec_BEC(nb_TF_k_zfc, nb_TF_n_zfc, eps, TM_frame_k, bool_txt)

        # TFs_rec_BEC_list_list.append(TF_rec_list) 
        pla_zfc_BEC_TF[ind_sim] = pla_zfc

    pla_zfc_BEC_TF_stat[0][ind_eps] = np.mean(pla_zfc_BEC_TF)


# %%
"""
------------------------------------------------------------------------------
                                PLOTTING
------------------------------------------------------------------------------
"""

# Plotting Parameters
plt.close('all')
nb_fig = 0
tit_fontsize = 20
txt_fontsize = 16
figsize_X = 8
figsize_Y = 6
dpi_ind = 100

'----------------------------- 2D Plots --------------------------------------'
'''
nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X + 2, figsize_Y), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

plt.plot(sim_time, ant_int_fad, color = 'maroon', label = r'$A_{int}$')
plt.plot(sim_time, ant_clr_vis, color = 'green', label = r'$clr_{vis}$')
plt.grid()

'XY Labels'
plt.title(r'Antenna Intermittence', size = tit_fontsize, fontweight = "bold")
plt.xlabel('Time [s]', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)
plt.ylabel('State [-1 / 0 / 1]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)

axis_R = axis.twinx()
plt.plot(sim_time[ind_fading], time_fading, color = 'goldenrod', label = r'$t_{fad}$')

plt.plot(sim_time, (ant_clr_dwl_PRX > 1), color = 'goldenrod')
plt.ylabel('Fading Time [s]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)

axis.legend(loc = 'upper left', fontsize = txt_fontsize)
axis.yaxis.get_offset_text().set_size(txt_fontsize)
axis_R.legend(loc = 'upper right')
axis_R.yaxis.get_offset_text().set_size(txt_fontsize)

ind_xlim_ini = 5 * ind_1op
ind_xlim_end = 6 * ind_1op

#xlim_ini =  sim_time[0]
# xlim_ini = sim_time[ind_xlim_ini]
xlim_ini = 30800
# xlim_end =  sim_time[-1]
# xlim_end =  sim_time[ind_xlim_end]
xlim_end = 31500

plt.xlim(xlim_ini, xlim_end)
plt.ylim(-np.max(time_fading) - 5, np.max(time_fading) + 5)

fig_label = ('plots/erasure correcting/AntInt_LagTimeExp-t_end=' + str(round(xlim_end)) + 's_pres.png')
fig0.savefig(fig_label)
'''

nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X + 2, figsize_Y + 2), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

axis.yaxis.get_offset_text().set_size(txt_fontsize)

plt.plot(sim_eps_SP, sim_eps_SP, '+-', color = 'black', label = 'Uncoded')

colors_exp = ['deepskyblue', 'dodgerblue', 'gold', 'crimson']
symbols_exp = ['*-', 'o-', 'v-', 's-']
pla_labels = []
for ind_R, R_real in enumerate(R_real_ran):
    pla_labels.append('zfec (R = ' + str(round(R_real, 2)) + ')')
    plt.plot(sim_eps_SP, pla_zfc_BEC_SP_stat[ind_R], symbols_exp[ind_R], color = colors_exp[ind_R], label = pla_labels[ind_R])

axis.set_yscale('log')

# legend1 = axis.legend(loc = 'upper left', fontsize = txt_fontsize)    
axis.yaxis.get_offset_text().set_size(txt_fontsize)
    
colors_exp = ['blue', 'orange', 'grey', 'red']
line0 = plt.axvline(x = eps_real_min, color = colors_exp[0], linestyle = '--', label = r'$\epsilon_{min}$')
line1 = plt.axvline(x = eps_vis_mean, color = colors_exp[2], linestyle = '--', label = r'$\epsilon_{vis}$')
line2 = plt.axvline(x = eps_real_max, color = colors_exp[3], linestyle = '--', label = r'$\epsilon_{max}$')

line_list = [line0, line1, line2]
# legend2 = axis.legend(handles = line_list, loc = 'lower right', fontsize = txt_fontsize)
axis.yaxis.get_offset_text().set_size(txt_fontsize)

axis.legend(loc = 'lower right', fontsize = txt_fontsize)    

plt.grid()

# axis.add_artist(legend1)
# axis.add_artist(legend2)

'XY Labels'
# plt.title('BLock (Space Packet = ' + str(len_SP) + ' bytes) Erasure Rate \n vs. \n Erasure Probability', size = tit_fontsize, fontweight = "bold")
plt.title('Code Performance (SP Encoding)', size = tit_fontsize, fontweight = "bold")
plt.xlabel(r'Erasure Probability ($\epsilon$) []', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)
plt.ylabel('BLock (SP) Erasure Rate (BLER) []', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)
'''
fig_label = ('plots/erasure correcting/zfec_SP_Performance_pres')
fig0.savefig(fig_label)
'''
pdf_name = ('plots/erasure correcting/zfec_SP_Performance_pres.pdf')
with PdfPages(pdf_name) as pdf:
    pdf.savefig(fig0)

nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X + 2, figsize_Y + 2), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

plt.plot(sim_eps_TF, sim_eps_TF, '+-', color = 'black', label = 'Uncoded')

colors_exp = ['deepskyblue', 'dodgerblue', 'gold', 'crimson']
symbols_exp = ['*-', 'o-', 'v-', 's-']
pla_labels = []

pla_labels.append('zfec (R = ' + str(round(R_real_ran[0], 2)) + ')')
plt.plot(sim_eps_TF, pla_zfc_BEC_TF_stat[0], symbols_exp[0], color = colors_exp[0], label = pla_labels[0])

# legend1 = axis.legend(loc = 'upper left', fontsize = txt_fontsize)    
axis.yaxis.get_offset_text().set_size(txt_fontsize)
axis.set_yscale('log')
    
colors_exp = ['blue', 'orange', 'grey', 'red']
line0 = plt.axvline(x = eps_real_min, color = colors_exp[0], linestyle = '--', label = r'$\epsilon_{min}$')
line1 = plt.axvline(x = eps_vis_mean, color = colors_exp[2], linestyle = '--', label = r'$\epsilon_{vis}$')
line2 = plt.axvline(x = eps_real_max, color = colors_exp[3], linestyle = '--', label = r'$\epsilon_{max}$')

line_list = [line0, line1, line2]
# legend2 = axis.legend(handles = line_list, loc = 'lower right', fontsize = txt_fontsize)
axis.yaxis.get_offset_text().set_size(txt_fontsize)

axis.legend(loc = 'lower right', fontsize = txt_fontsize)    

plt.grid()

# axis.add_artist(legend1)
# axis.add_artist(legend2)

'XY Labels'
plt.title('Code Performance (TF Encoding)', size = tit_fontsize, fontweight = "bold")
plt.xlabel(r'Erasure Probability ($\epsilon$) []', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)
plt.ylabel('BLock (TF) Erasure Rate (BLER) []', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)
'''
fig_label = ('plots/erasure correcting/zfec_TF_Performance_pres')
fig0.savefig(fig_label)
'''
pdf_name = ('plots/erasure correcting/zfec_TF_Performance_pres.pdf')
with PdfPages(pdf_name) as pdf:
    pdf.savefig(fig0)
'''
nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X + 2, figsize_Y + 2), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

plt.plot(len_SP_ran, b_SP_daily / 10**6, '*-')
b_max = b_daily * (1 / R_vis_mean)
# plt.axhline(y = b_max / 10**6, linestyle = '--', label = r'$b_{max}$')
plt.axvline(x = 64, linestyle = '--')

plt.grid()

'XY Labels'
plt.title('Daily Data vs. SP size', size = tit_fontsize, fontweight = "bold")
plt.xlabel('SP size [bytes]', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)
plt.ylabel('Daily Data [MB]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)

axis.legend(['R=2/3', 'R=3/4', 'R=5/6', 'R=7/8', 'R=1'])

fig_label = ('plots/erasure correcting/b_SP_sizes')
fig0.savefig(fig_label)
'''
