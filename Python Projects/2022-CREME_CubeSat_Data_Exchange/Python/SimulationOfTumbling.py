#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 14:48:32 2022

@author: d.opoka
"""

"""
##############################################################################
                            SIMULATION OF TUMBLING
##############################################################################
"""

import numpy as np  
import matplotlib.pyplot as plt  
from matplotlib.pyplot import figure
import yaml
from yaml.loader import SafeLoader
import scipy.spatial as ss
from matplotlib.backends.backend_pdf import PdfPages

# %%
"""
------------------------------------------------------------------------------
                             ASTRODYNAMICS OF TUMBLING
------------------------------------------------------------------------------
"""

"""
Trigonometric functions are in degrees.

The Earth-Centered Inertial (ECI) frame is used.

The simulation begins on the most recent Vernal/Spring Equinox: 
    20/03/2022 - 15:33:00
"""

'-------------------------------- FUNCTIONS ----------------------------------'
# Reference frame rotation about its X-axis
def R1_RX(angle_X):
    return np.array([[1, 0, 0],
                     [0, np.cos(angle_X), np.sin(angle_X)],
                     [0, -np.sin(angle_X), np.cos(angle_X)]])

# Reference frame rotation about its Y-axis
def R1_RY(angle_Y):
    return np.array([[np.cos(angle_Y), 0, -np.sin(angle_Y)],
                     [0, 1, 0],
                     [np.sin(angle_Y), 0, np.cos(angle_Y)]])

# Reference frame rotation about its Z-axis
def R1_RZ(angle_Z):
    return np.array([[np.cos(angle_Z), np.sin(angle_Z), 0],
                     [-np.sin(angle_Z), np.cos(angle_Z), 0],
                     [0, 0, 1]])

# Reference Frame 2 in terms of Reference Frame 1 converter
def R2_R1(DTM, coord):
    # coord - Input coordinates in terms of array([[X], [Y], [Z]])
    # DTM - Direct Transformation Matrix from R2 to R1
    coord_conv = np.matmul(DTM, coord)     # Converted Coordinates
    return coord_conv

# Calculates Perpendicular Line Point (3D)
# The perpendicular line with respect to a line passing through A & C 
# is defined as a line passing through B and D.
# D is the result of this function.
'Source:' 
'https://math.stackexchange.com/questions/4347497/find-a-point-on-a-line-that-creates-a-perpendicular-in-3d-space'
def PerpLinePoint(A, B, C):
    t = (np.dot((B - A), (C - A))) / (np.dot((C - A), (C - A)))
    D = A + t * (C - A)
    return D

# Magnitude of a 3D Vector (AB)
def magVec(A, B):
    mag = np.sqrt((B[0] - A[0]) ** 2 + (B[1] - A[1]) ** 2 + (B[2] - A[2]) ** 2)
    return mag

def unit2dB(input_unit):
    output_dB = 10 * np.log10(input_unit)
    return output_dB

def dB2unit(input_dB):
    output_unit = 10 ** (input_dB / 10)
    return output_unit

def read_antenna_rad_pat_file(filepath): 
    f=open(filepath,"r")
    lines=f.readlines()
    theta = []
    phi = []
    PCG_gain=[]
    for x in lines:
        split_line = x.split('\t')
        theta.append(float(split_line[0]))
        phi.append(float(split_line[1]))
        PCG_gain.append(float(split_line[4]))
    f.close()
    
    return [np.array(theta), np.array(phi), np.array(PCG_gain)]

def radiation_pattern(theta, phi, PCG_gain):
    PCG_gain_scal = np.power(10, PCG_gain/10)
    x=[]
    y=[]
    z=[]
    for i in range(len(theta)):
        t = np.deg2rad(theta[i])
        p = np.deg2rad(phi[i])
        x.append(PCG_gain_scal[i]*np.sin(t)*np.cos(p))
        y.append(PCG_gain_scal[i]*np.sin(t)*np.sin(p))
        z.append(PCG_gain_scal[i]*np.cos(t))
        
    return np.array(x), np.array(y), np.array(z)

def compute_angle_vector(v1, v2): 
    dot_product = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    norm_v1 = np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    norm_v2 = np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    angle =  np.arccos(dot_product/(norm_v1*norm_v2))
    
    return angle

def compute_depointing_angle(theta, phi): 
    depointing_angle_mat = []
    for i in range(len(theta)):
        t = theta[i]*np.pi/180
        p = phi[i] *np.pi/180
        pointing_vector = [np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)]
        if t >= 0: 
            antenna_normal_direction = [0,1,0]
        else:
            antenna_normal_direction = [0,-1,0]

        depointing_angle = compute_angle_vector(pointing_vector, antenna_normal_direction)*180/np.pi
        depointing_angle_mat.append(depointing_angle)

    mean_depointing_loss = []
    percentile75_depointing_loss = []
    percentile95_depointing_loss = []
    step_lissage = 3
    angle_step = np.arange(0,90,step_lissage)
    for k in angle_step: 
        gain_per_degree = []
        for i in range(len(depointing_angle_mat)):
            if k < depointing_angle_mat[i] and k+step_lissage > depointing_angle_mat[i]:
                gain_per_degree.append(PCG_gain[i])
        gain_per_degree = np.array(gain_per_degree)
        mean_depointing_loss.append(np.mean(gain_per_degree))
        percentile75_depointing_loss.append(np.percentile(gain_per_degree, 10))
        percentile95_depointing_loss.append(np.percentile(gain_per_degree, 5))

    return [mean_depointing_loss, angle_step]

def read_antenna_range_file(filepath): 
    f=open(filepath,"r")
    lines=f.readlines()
    ant_rpt_range = []
    for x in lines:
        split_line = x.split('       ')
        if len(split_line)>= 6:
            ant_rpt_range.append(float(split_line[2]))
    f.close()
    
    return np.array(ant_rpt_range)

# Determine index/position of minimum value
def min_ind(est, real):
    ind = np.abs(est - real)
    ind = np.where(ind == np.amin(ind))
    
    return ind


# %%
'--------------------------------- YAML INPUT FILE ---------------------------'
with open('data/input.yaml') as f:
    data = yaml.load(f, Loader=SafeLoader)

'--------------------- REFERENCE FRAME FOR BODY OF SATELLITE -----------------'
# General Parameters
ear_rad = data['EARTH']['rad'] * 1000              # Radius of the Earth [m]
ear_day = data['EARTH']['day']                     # Day on Earth [s]
ear_ome = np.deg2rad(360.00) / ear_day             # Angular velocity of the Earth [rad / s]
spd_lgt = data['spd_lgt']                          # Speed of Light [m / s] 
stf_bzm_dB = 228.5991672                           # Stefan-Boltzmann Const. [dB]
stf_bzm_JK = 1.380649 * (10 ** -23)                # Stefan-Boltzmann Const. [J / K]

# Orbital Elements (from orbit 2 - alt = 600 km)
orb_OmA = np.deg2rad(data['ORBIT']['asc_node'])    # Longitude of ascending node [rad]
orb_inc = np.deg2rad(data['ORBIT']['inc'])         # Inclination [rad]            
orb_omP = np.deg2rad(data['ORBIT']['arg_per'])     # Argument of periapsis/perigee [rad]
orb_ecc = data['ORBIT']['ecc']                     # Eccentricity []
orb_sma = data['ORBIT']['sma'] * 1000              # Semi-major axis (a) [m]
orb_Tnu = np.deg2rad(data['ORBIT']['tru_ano'])     # True anomaly [rad] (at time 0 s)

orb_per = data['ORBIT']['per'] * 60                # Orbit Period/Duration [s]

# Semi-minor axis0is (b) [m]
orb_smb = np.sqrt((1 - orb_ecc ** 2) * orb_sma ** 2)
# Orbit velocity [m / s]
orb_vel = (np.pi * (orb_sma + orb_smb)) / orb_per

# Direct Transformation Matrix from Orbital Plane to the ECI
# Rot with OmA (Z) -> Rot with inc (X) -> Rot with omP (Z)
RO_RE = np.matmul(R1_RZ(orb_omP), 
                  np.matmul(R1_RX(orb_inc), 
                            R1_RZ(orb_OmA)))

# Translation Matrix from Satellite Reference Frame to the Orbital Plane [m]
'Vector in Orbital Plane of the current position of satellite'
curr_pos = np.array([orb_sma * np.cos(orb_Tnu), orb_smb * np.sin(orb_Tnu), 0])
# Direct Transformation Matrix from Satellite Reference Frame to the ECI    
RR_RE = np.matmul(R1_RZ(np.deg2rad(90) + orb_Tnu), RO_RE)
curr_pos = np.matmul(curr_pos, RO_RE)

# Angles from Satellite Reference Frame to Satellite Body Frame
'The satellite must be perpendicular to the Earth\'s Magnetic Field'
sat_alf_int = np.deg2rad(+00.00)                   # Initial alfa_body (X) [rad]
sat_bet_int = np.deg2rad(+90.00)                   # Initial beta_body (Y) [rad]
sat_gam_int = np.deg2rad(+00.00)                   # Initial gamma_body (Z) [rad]

# Satellite Angular Velocity according to the Satellite Body Frame
'These are a function of data from the on-board sensors'
sat_alf_dot = np.deg2rad(00.00) / 1                # alfa^dot_body (X) [rad / s]
sat_bet_dot = np.deg2rad(00.00) / 1                # beta^dot_body (Y) [rad / s]
sat_gam_dot = np.deg2rad(00.00) / 1                # gamma^dot_body (Z) [rad / s]

# Translation Matrix from Satellite Reference Frame to the Satellite Body Frame [m]
'Vector from the CoG to the origin of the Body'
CgOb = np.array([data['SAT']['CgOb_X'], 
                 data['SAT']['CgOb_Y'], 
                 data['SAT']['CgOb_Z']])
# Direct Transformation Matrix from Satellite Reference Frame to Satellite Body Frame
# -Z_R is aligned with X_O, pointing towards the sun
RB_RE = R1_RY(np.deg2rad(-90.00))

# Current position of the Origin of the Satellite 
curr_Ob = np.matmul(CgOb, RB_RE)

'-------------------- REFERENCE FRAME FOR GROUND SEGMENT ---------------------'
gse_alf = np.deg2rad(data['GSE']['lon'])           # Longitude (E) [rad]
gse_bet = np.deg2rad(data['GSE']['lat'])           # Latitude (N) [rad]
gse_alt = ear_rad + data['GSE']['alt']             # Altitude [m]

# Angle of the GS at time 0 s (Vernal Equinox, Z) [rad]
# Earth's X-axis pointing to the Sun at 0°00' North, 51°08' West at time 0s
# Source: https://www.timeanddate.com/worldclock/sunearth.html?iso=20220320T1532
gse_tht = np.deg2rad(51 + (8 / 60)) + gse_alf

# Direct Transformation Matrix from Ground Segment RF to ECI
# Rot with beta_gse (Z) -> Rot with tht_gse (Y) -> Rot with 90° (Z)
# X -> East, Y -> North, and Z points towards Space
RGS_RE = np.matmul(R1_RZ(np.deg2rad(90)),
                   np.matmul(R1_RY(gse_bet), 
                             R1_RZ(gse_tht)))
# Direct Transformation Matrix from ECI TO Ground Segment RF
RE_RGS = np.linalg.inv(RGS_RE)

# The Ground Segment has a visibility up to the Clearance/Horizon Angle. 
gse_elv_ang = np.deg2rad(data['GSE']['ele'])       # Ground segment elevation Angle [rad]
gse_clr_ang = np.deg2rad(90) - gse_elv_ang         # Ground segment horizon Angle [rad]

# Current position of GS
curr_gse = np.array([np.cos(gse_bet) * np.cos(gse_tht),
                     np.cos(gse_bet) * np.sin(gse_tht),
                     np.sin(gse_bet)]) * gse_alt

'------------------------ SATELLITE ANTENA MODEL -----------------------------'                 
# Distance between Ground Segment and Satellite for Downlink (DGS) [m]
DGS_dwl = 1932.26 * 1000

# Antenna Data for Range
filepath = 'data/ant_range.txt'

# List of antenna radiation pattern range (Edge-to-edge distance) [km]    
ant_rpt_ran_list = read_antenna_range_file(filepath)
# Mean of 75% antenna radiation pattern range (Edge-to-edge distance) [m]    
ant_rpt_ran = np.percentile(ant_rpt_ran_list, 75) * 1000
ant_vol_real = (2 * ant_rpt_ran) ** 3              # Real volume of a Cube [m^3]

# Antenna Data for Radiation Pattern
filepath = 'data/BandeS_PCD-PCG-2218M-norme.txt'

# Radiation Pattern Angles & Gain
[theta, phi, PCG_gain] = read_antenna_rad_pat_file(filepath)

# Radiation XYZ Points
X_rad, Y_rad, Z_rad = radiation_pattern(theta, phi, PCG_gain)
# Normalizing Points 
XYZ_points = np.zeros((len(X_rad), 3))
for ind in range(len(X_rad)):
    XYZ_points[ind][0] = X_rad[ind] / np.max(np.abs(X_rad))
    XYZ_points[ind][1] = Y_rad[ind] / np.max(np.abs(Y_rad))
    XYZ_points[ind][2] = Z_rad[ind] / np.max(np.abs(Z_rad))

ant_rpt_hull = ss.ConvexHull(XYZ_points)
ant_vol_rpt_norm = ant_rpt_hull.volume             # Antenna normalized radiation pattern Volume [unit^3]
ant_vol_norm = 2 ** 3                              # Normalized volume of a cube [unit^3]
ant_vol_real = (2 * ant_rpt_ran) ** 3              # Real volume of a cube [m^3]
# Antenna real radiation pattern volume [m^3]
ant_rpt_vol = (ant_vol_rpt_norm / ant_vol_norm) * ant_vol_real

# The Radiation Pattern is recalculated to be a cone.
# The cone is composered of an height, radius, and angle alfa.
# ant_rpt_rad -> radius, ant_rpt_ran -> height, ant_rpt_alf -> alfa

# Antenna radiation pattern radius [m]
ant_rpt_rad = np.sqrt((3 * ant_rpt_vol) / (2 * np.pi * ant_rpt_ran))
# Antenna radiation pattern cone angle [rad]
ant_rpt_alf = np.arctan((ant_rpt_rad) / ant_rpt_ran)

# Inputting the Mean Depointing Loss of S-band Antennas
[mean_depointing_loss, angle_step] = compute_depointing_angle(theta, phi)

# Curve Fitting the Mean Depointing Loss to a 5% Error (func of degrees)
error = 100
n_poly = 0
while error >= 5:
    n_poly += 1
    
    xdata = angle_step+(angle_step[1]-angle_step[0])/2
    ydata = mean_depointing_loss
    polyfit_data = np.polyfit(xdata, ydata, n_poly)
    polyfit_data = np.poly1d(polyfit_data)
    
    error = np.abs((polyfit_data(xdata) - ydata) / ydata) * 100
    error = np.mean(error)
ant_dpt_loss = polyfit_data

# Antenna FoV Line Coordinates [m]
ant_FoV = np.array([[+00.00, +ant_rpt_ran, +00.20],
                    [+00.00, -ant_rpt_ran, +00.20]])    

# For Antenna Received Power:
# Source from Bilan_Liaison_Creme.xlsx found in the data folder. 
ant_lpl = 3                                        # Polarisation losses [dB]
ant_lat = 1.1                                      # Atmospheric losses [dB] 

# Downlink parameters
ant_dwl_Eb_N0_req = 2.50                           # Required normalized SNR [dB]
ant_dwl_PTX = 3.01                                 # Transmitter output power [dBW]
ant_dwl_LTX = 34.80                                # Transmitter losses [dB]
ant_dwl_GTX = 5.00                                 # Transmitter antenna gain [dB]

ant_frq = 2263.5 * (10 ** 6)                       # Antenna frequency [Hz]
ant_bwd = 265 * (10 ** 6)                          # Antenna bandwidth [Hz]
# Free apace loss [dB] as a function of DGS
ant_dwl_FSL = (unit2dB(((4 * np.pi * DGS_dwl * ant_frq) / (spd_lgt)) ** 2) +
           ant_lpl + ant_lat)                           

ant_dwl_LGP = 3.00                                 # Ground pointing losses [dB]
ant_dwl_G_T = 6.26                                 # Gain-to-noise-temperature [dB/K]
ant_k = stf_bzm_dB                                 # Stefan-Boltzmann const. [dB]
ant_dwl_b_r_spd = 85 * 1000                        # Bit-rate [bps]
ant_dwl_b_r_dB = unit2dB(ant_dwl_b_r_spd)          # Bit-rate [dB]
ant_dwl_smr = 3.00                                 # System margin [dB]

# Received power [dB]
ant_dwl_Eb_N0_real = (ant_dwl_PTX - ant_dwl_LTX + ant_dwl_GTX - ant_dwl_FSL - 
                      ant_dwl_LGP + ant_dwl_G_T + ant_k - ant_dwl_b_r_dB - ant_dwl_smr)
# System Margin [dB] (any received power above the margin is communicated)
ant_dwl_mar = ant_dwl_Eb_N0_real - ant_dwl_Eb_N0_req 

# Uplink parameters
# Not required as the COP-1 is already a solution for the Uplink communication.


# %%
'--------------------------------- SIMULATION --------------------------------'
# Parameters
sim_tini = data['SIM']['int_time']                 # Beginning of time [s]
sim_tend = data['SIM']['end_time']                 # Ending of time [s]
sim_tstp = data['SIM']['time_stp']                 # Time Step [s]
# Number of elements for time [#]
sim_tele = round((sim_tend - sim_tini) / sim_tstp)              
# Time [s]
sim_time = np.linspace(sim_tini, sim_tend, sim_tele)

'For variables that DO NOT CHANGE with Sat. Angular Velocities:'
sat_pos = np.zeros([sim_tele, 3])
gse_pos = np.zeros([sim_tele, 3])
sat_gse = np.zeros([sim_tele, 3])

# Satellite Orbit Dynamics
if (orb_ecc == 00.00):
    # For circular orbit
    orb_Tnu_dot = orb_vel / orb_sma                # Orbit angular velocity [rad / s]
else:
    # For elliptic orbit
    # Orbit distance from orbitting Focus to Satellite [m]
    orb_dis = (orb_sma * (1 - orb_ecc ** 2)) / (1 + orb_ecc * np.cos(orb_Tnu))
    # Orbit angular velocity [rad / s] 
    orb_Tnu_dot = (2 * np.pi * orb_sma * orb_smb) / (orb_per * orb_dis)

for sim_iter, curr_time in enumerate(sim_time): 
    # For elliptic orbit
    if (orb_ecc != 00.00):
        # Orbit distance from orbitting Focus to Satellite [m]
        orb_dis = (orb_sma * (1 - orb_ecc ** 2)) / (1 + orb_ecc * np.cos(orb_Tnu))
        # Orbit angular velocity [rad / s] 
        orb_Tnu_dot = (2 * np.pi * orb_sma * orb_smb) / (orb_per * orb_dis ** 2)

    # Updating True Anomaly
    orb_Tnu = orb_Tnu_dot * curr_time
    
    'Satellite Position in Orbital Plane' 
    curr_sat_pos = np.array([orb_sma * np.cos(orb_Tnu), orb_smb * np.sin(orb_Tnu), 0])
    curr_sat_pos = np.matmul(curr_sat_pos, RO_RE)
    sat_pos[sim_iter] = curr_sat_pos
    
    # Updating the Theta angle of GS
    gse_tht += ear_ome * sim_tstp
    
    'Ground Segment Position in ECI'
    curr_gse_pos = np.array([np.cos(gse_bet) * np.cos(gse_tht),
                             np.cos(gse_bet) * np.sin(gse_tht),
                             np.sin(gse_bet)]) * gse_alt
    gse_pos[sim_iter] = curr_gse_pos
    
    # Updating DTM from Ground Segment RF to ECI
    RGS_RE = np.matmul(R1_RZ(np.deg2rad(90)),
                   np.matmul(R1_RY(gse_bet), 
                             R1_RZ(gse_tht)))
    # Updating DTM from ECI to Ground Segment RF
    RE_RGS = np.linalg.inv(RGS_RE)
    # Current Position of Satellite according to Ground Segment RF
    curr_sat_gse = np.matmul(curr_sat_pos, RE_RGS) - np.matmul(curr_gse_pos, RE_RGS)
    sat_gse[sim_iter] = curr_sat_gse
    

# Parameters
sim_aini = 0                                      # Beginning of angular velocity [rad / s]
max_angle = 15.00                                 # Maximum Angular Velocity [deg]
sim_aend = np.deg2rad(max_angle) / 1              # Ending of angular velocity [rad / s]
sim_aele = round(max_angle) + 1                   # Number of Elements for angular velocity
# Range of Z-axis angular velocities [rad / s]
sim_gam_dot = np.linspace(sim_aini, sim_aend, sim_aele)

DGA_stat = np.zeros([sim_aele, sim_tele])
DAS_stat = np.zeros([sim_aele, sim_tele])
ant_alf_stat = np.zeros([sim_aele, sim_tele])

'For variables that CHANGE with Sat. Angular Velocities:'
# Satellite Tumbling about the Z-axis [0 rad /s, 15 rad/s] 
for sim_iter_ang, sat_gam_dot in enumerate(sim_gam_dot):
    DGS = np.zeros(sim_tele)
    DGA = np.zeros(sim_tele)
    DAS = np.zeros(sim_tele)
    ant_alf = np.zeros(sim_tele)
    
    for sim_iter, curr_time in enumerate(sim_time): 

        # Updating Satellite Body Angles
        sat_alf = sat_alf_dot * sim_tstp
        sat_bet = sat_bet_dot * sim_tstp
        sat_gam = sat_gam_dot * sim_tstp
        # Updating DTM from Satellite Body RF to the ECI
        RB_RE = np.matmul(R1_RX(sat_alf), 
                          np.matmul(R1_RY(sat_bet), 
                                    np.matmul(R1_RZ(sat_gam), RB_RE)))
        # Updating Origin of the Satellite Body Reference Frame
        curr_Ob = np.matmul(CgOb, RB_RE)
        
        'Antenna Attitude'
        # Updating the Field of View direction
        curr_ant = np.zeros([2, 3])
        curr_ant[0] = np.matmul(ant_FoV[0], RB_RE) + sat_pos[sim_iter]+ curr_Ob
        curr_ant[1] = np.matmul(ant_FoV[1], RB_RE) + sat_pos[sim_iter] + curr_Ob
        
        # Updating the current Point of Perpecndicular between Antenna FoV and GS
        gse_ant_point = PerpLinePoint(curr_ant[0], gse_pos[sim_iter], curr_ant[1])
        
        # Current Distance between GS and Satellite (curr_DGS) [m]
        DGS[sim_iter] = magVec(gse_pos[sim_iter], sat_pos[sim_iter])
        # Distance between GS and Antenna FoV Point (DGA) [m]
        DGA[sim_iter] = magVec(gse_pos[sim_iter], gse_ant_point)
        # Distance between Antenna FoV Point and Satellite (DAS) [m]
        DAS[sim_iter] = magVec(gse_ant_point, sat_pos[sim_iter])
        # Antenna angle between DAS and DGS [deg]
        ant_alf[sim_iter] = np.rad2deg(np.arctan(DGA[sim_iter] / DAS[sim_iter]))
       
    DGA_stat[sim_iter_ang] = DGA
    DAS_stat[sim_iter_ang] = DAS
    ant_alf_stat[sim_iter_ang] = ant_alf

ant_alf_avg = np.zeros(sim_tele)
ant_alf_std = np.zeros(sim_tele)

for ind_tele in range(sim_tele):
    ant_alf_aele = np.zeros(sim_aele)
    
    for ind_aele in range(sim_aele):
        ant_alf_aele[ind_aele] = ant_alf_stat[ind_aele][ind_tele]
        
    ant_alf_avg[ind_tele] = np.mean(ant_alf_aele)
    ant_alf_std[ind_tele] = np.std(ant_alf_aele)

ant_alf_inp = ant_alf_avg
    
clr_DGS = np.zeros(sim_tele)                   # Clearance distance between GS and Satellite [m]
clr_tht = np.zeros(sim_tele)                   # Clearance angle between GS and Satellite [°]
clr_vst = np.zeros(sim_tele)                   # Counting for visibility [0/1]

ant_clr_alf = np.zeros(sim_tele)               # Visible/Cleared antenna angle to GS [°]
ant_clr_int = np.zeros(sim_tele)               # Antenna intermittency for communication to GS [0/1]
ant_clr_dwl_PRX = np.zeros(sim_tele)           # Visible/Cleared antenna downlink received power [dB]
ant_clr_dwl_C = np.zeros(sim_tele)             # Visible/Cleared antenna downlink channel capacity [bps]
  
for sim_iter in range(sim_tele): 
    # Antenna Ground Pointing Losses from Mean Depointing Loss [dBi]
    ant_dwl_LGP = ant_dpt_loss(ant_alf_inp[sim_iter]) * -1
    # Antenna Free Space loss[dB]: as a function of DGS
    ant_dwl_FSL = (unit2dB(((4 * np.pi * DGS[sim_iter] * ant_frq) / (spd_lgt)) ** 2) +
                   ant_lpl + ant_lat)     
    # Real Signal-to-Noise from Satellite to GS 
    ant_dwl_Eb_N0_real = (ant_dwl_PTX - ant_dwl_LTX + ant_dwl_GTX - ant_dwl_FSL - 
                          ant_dwl_LGP + ant_dwl_G_T + ant_k - ant_dwl_b_r_dB - ant_dwl_smr)
    # Received Power from Satellite to GS
    ant_dwl_PRX = (ant_dwl_Eb_N0_real - ant_dwl_Eb_N0_req) - ant_dwl_mar
    # Signal-to-Noise Ratio [dB]
    ant_SNR = (ant_dwl_PRX + ant_dwl_mar) + unit2dB(ant_dwl_b_r_spd / ant_bwd)
    # Channel Capacity of the Antenna Downlink [bps]
    ant_dwl_C = ant_bwd * np.log2(1 + dB2unit(ant_SNR))

    # Determining the Clearance/Visibility between Satellite and GS
    if (sat_gse[sim_iter][2] > 0):
        clr_alt = sat_gse[sim_iter][2]
        clr_rad = np.sqrt(sat_gse[sim_iter][0] ** 2 + sat_gse[sim_iter][1] ** 2)
        curr_clr_tht = np.arctan(clr_rad / clr_alt)
        
        # If Visibility Angle is less than the Horizon/GS Clearance Angle 
        if (curr_clr_tht < gse_clr_ang):
            clr_DGS[sim_iter] = DGS[sim_iter]
            clr_tht[sim_iter] = np.rad2deg(curr_clr_tht)
            clr_vst[sim_iter] += 1
            
            ant_clr_alf[sim_iter] = ant_alf_inp[sim_iter]

            if (ant_dwl_PRX > 0): 
                ant_clr_dwl_PRX[sim_iter] = ant_dwl_PRX
                ant_clr_dwl_C[sim_iter] = ant_dwl_C 

data = [clr_vst, ant_clr_dwl_PRX]
data_str = 'data/ant_clr_avg-t_end='+ str(sim_time[-1]) + 's-t_stp=' + str(sim_tstp) + 's.dat'
np.savetxt(data_str, data)


# %%
'''
'--------------------------------- ANALYSIS ----------------------------------'
'''
ind_1op = min_ind(sim_time, orb_per)
nb_orb = int(np.floor(sim_tend / orb_per))
time_clr_dwl_avg = (np.count_nonzero(ant_clr_dwl_PRX) * sim_tstp) / 60
print('Average Visibility Time (dwl_PRX) = ' + str(time_clr_dwl_avg) + ' min \n')

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

ind_xlim_ini = 5 * ind_1op[0][0]
ind_xlim_end = 6 * ind_1op[0][0]

ind_xlim_ini = 0
ind_xlim_end = -1

'----------------------------- 2D Plots --------------------------------------'

nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

plt.plot(sim_time, ant_alf_inp.transpose(), color = 'salmon', label = 'Look Angle')
plt.ylabel(r'Look Angle [°]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)
'''
plt.axvline(x = orb_per, color = 'r', linestyle = '--', label = 'nb_orb = ' + str(nb_orb))
for ind in range(nb_orb - 1):
    plt.axvline(x = orb_per * (ind + 2), color = 'r', linestyle = '--')  
'''
plt.grid()

# xlim_ini = sim_time[ind_xlim_ini]
xlim_ini = 30800
# xlim_end =  sim_time[ind_xlim_end]
xlim_end = 31600
plt.xlim(xlim_ini, xlim_end)


'XY Labels'
plt.title(r'Satellite Attitude', size = tit_fontsize, fontweight = "bold")
plt.xlabel('Time [s]', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)

axis_R = axis.twinx()
axis_R.plot(sim_time, DGS / 10**6, color = 'navy', label = 'Slant Range')
plt.ylabel('Slant Range [Mm]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)

axis.legend(loc = 'upper left', fontsize = txt_fontsize)
axis.yaxis.get_offset_text().set_size(txt_fontsize)
axis_R.legend(loc = 'upper right', fontsize = txt_fontsize)
axis_R.yaxis.get_offset_text().set_size(txt_fontsize)
'''
fig_label = ('plots/orbital mechanics/DGS_alpha/DGS_alpha-t_end='+ str(round(xlim_end)) + 's-t_stp=' + str(sim_tstp) + 's.png')
fig0.savefig(fig_label)
'''
pdf_name = 'plots/orbital mechanics/DGS_alpha/DGS_alpha-t_end='+ str(round(xlim_end)) + 's-t_stp=' + str(sim_tstp) + 's.pdf'
with PdfPages(pdf_name) as pdf:
    pdf.savefig(fig0)

nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

plt.plot(sim_time, clr_vst.transpose(), color = 'green', label = '$clr_{vis}$')
'''
plt.axvline(x = orb_per, color = 'r', linestyle = '--', label = 'nb_orb = ' + str(nb_orb))
for ind in range(nb_orb - 1):
    plt.axvline(x = orb_per * (ind + 2), color = 'r', linestyle = '--')
'''
plt.grid()

plt.xlim(xlim_ini, xlim_end)

'XY Labels'
plt.title('Antenna Clearance', size = tit_fontsize, fontweight = "bold")
plt.xlabel('Time [s]', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)
plt.ylabel('Antenna Visibility [0/1]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)

axis_R = axis.twinx()
axis_R.plot(sim_time, ant_clr_dwl_PRX.transpose(), color = 'goldenrod', label = '$M_{S}$')
plt.ylabel(r'Margin Power ($M_{S}$) [dB]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)

axis.legend(loc = 'upper left', fontsize = txt_fontsize)
axis.yaxis.get_offset_text().set_size(txt_fontsize)
axis_R.legend(loc = 'upper right', fontsize = txt_fontsize)
axis_R.yaxis.get_offset_text().set_size(txt_fontsize)
'''
fig_label = ('plots/orbital mechanics/ant_clr/ant_clr-t_end='+ str(round(xlim_end)) + 
             's-t_stp=' + str(sim_tstp) + 's.png')
fig0.savefig(fig_label)
'''
pdf_name = ('plots/orbital mechanics/ant_clr/ant_clr-t_end='+ str(round(xlim_end)) + 
            's-t_stp=' + str(sim_tstp) + 's.pdf')
with PdfPages(pdf_name) as pdf:
    pdf.savefig(fig0)

# %%

ant_alf_hist = [ant_alf_stat[0], ant_alf_stat[2], ant_alf_stat[-1], ant_alf_avg]
ant_alf_hist_str = ['$\dot{\gamma} = 0 \: [°/s]$', 
                    '$\dot{\gamma} = 2 \: [°/s]$', 
                    '$\dot{\gamma} = 15 \: [°/s]$', 
                    'avg_aele']
ant_alf_hist_str_title = ['gd=0', 'gd=2', 'gd=15', 'avg_aele']
ant_alf_hist_col = ['red', 'green', 'blue', 'salmon']

nb_bins = int(np.ceil(np.sqrt(sim_tele)))

xlim_ini =  sim_time[0]
xlim_ini = sim_time[ind_xlim_ini]
xlim_ini = 30800
xlim_end =  sim_time[-1]
xlim_end =  sim_time[ind_xlim_end]
xlim_end = 31600

# Plotting Parameters
plt.close('all')
nb_fig = 0
tit_fontsize = 18
txt_fontsize = 14
figsize_X = 8
figsize_Y = 6
dpi_ind = 100

nb_fig += 1
fig0 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi = dpi_ind)
plt.ion()
axis = plt.axes()

plt.plot(sim_time, ant_alf_stat[0], color = ant_alf_hist_col[0], label = ant_alf_hist_str[0])
plt.grid()

plt.xlim(xlim_ini, xlim_end)

'XY Labels'
plt.title(r'$\alpha_{A}$ (' + ant_alf_hist_str[0]
          + ') vs. Time', size = tit_fontsize, fontweight = "bold")
plt.xlabel('Time [s]', size = txt_fontsize)
plt.xticks(fontsize = txt_fontsize)
plt.ylabel(r'Angle between GS and Satellite ($\alpha_{A}$) [°]', size = txt_fontsize)
plt.yticks(fontsize = txt_fontsize)
'''
axis.legend(loc = 'upper left')
axis.yaxis.get_offset_text().set_size(txt_fontsize)
'''
fig_label = ('plots/orbital mechanics/alpha/alpha_' + ant_alf_hist_str_title[0] + 
             '-t_end='+ str(round(xlim_end)) + 's-t_stp=' + str(sim_tstp) + 's.png')
fig0.savefig(fig_label)


# %%

for ind_hist in range(3):
    nb_fig += 1
    fig0 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi = dpi_ind)
    plt.ion()
    axis = plt.axes()

    plt.plot(sim_time, ant_alf_stat[0], color = ant_alf_hist_col[0], label = ant_alf_hist_str[0])
    plt.plot(sim_time, ant_alf_stat[ind_hist + 1], color = ant_alf_hist_col[ind_hist + 1], label = ant_alf_hist_str[ind_hist + 1])
    
    plt.grid()
    plt.xlim(xlim_ini, xlim_end)
    
    axis.legend(loc = 'lower right')
    axis.yaxis.get_offset_text().set_size(txt_fontsize)
        
    plt.xlim(xlim_ini, xlim_end)
    
    'XY Labels'
    plt.title(r'$\alpha_{A}$ vs. Time', size = tit_fontsize, fontweight = "bold")
    plt.xlabel('Time [s]', size = txt_fontsize)
    plt.xticks(fontsize = txt_fontsize)
    plt.ylabel(r'Angle between GS and Satellite ($\alpha_{A}$) [°]', size = txt_fontsize)
    plt.yticks(fontsize = txt_fontsize)

    fig_label = ('plots/orbital mechanics/alpha/alpha_' + ant_alf_hist_str_title[ind_hist] + 
                 '-t_end=' + str(round(xlim_end)) + 's-t_stp=' + str(sim_tstp) + 's.png')
    fig0.savefig(fig_label)


# %%

# Plotting Parameters
plt.close('all')
nb_fig = 0
tit_fontsize = 18
txt_fontsize = 14
figsize_X = 16
figsize_Y = 10
dpi_ind = 100

xlim_ini = sim_time[0]
xlim_end = sim_time[-1]

for ind_hist, elem_hist in enumerate(ant_alf_hist):
    nb_fig += 1
    
    fig0 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi = dpi_ind)
    ax1 = fig0.add_subplot(121)
    ax2 = fig0.add_subplot(122)

    plt.ioff()
    
    ax1.plot(sim_time, elem_hist, color = ant_alf_hist_col[ind_hist])
    ax1.grid()
    
    ax1.set_xlim(xlim_ini, xlim_end)
    
    'XY Labels'
    ax1.set_title(r'$\alpha_{A}$ (' + ant_alf_hist_str[ind_hist]
                        + ') vs. Time', size = tit_fontsize, fontweight = "bold")
    ax1.set_xlabel('Time [s]', size = txt_fontsize)
    ax1.tick_params(axis = 'x', labelsize = txt_fontsize)
    ax1.set_ylabel(r'Angle between GS and Satellite ($\alpha_{A}$) [°]', size = txt_fontsize)
    ax1.tick_params(axis = 'y', labelsize = txt_fontsize)
    
    ax2.hist(x = elem_hist, bins = int(np.ceil(np.sqrt(sim_tele))), density = True, color = ant_alf_hist_col[ind_hist])
    # Even though it says density, as we are with discrete values, this is a probability mass function
    ax2.grid()
    
    'XY Labels'
    ax2.set_title(r'PMF (' + ant_alf_hist_str[ind_hist] + ')', 
                  size = tit_fontsize, fontweight = "bold")
    ax2.set_xlabel(r'$\alpha_{A}$ [°]', size = txt_fontsize)
    ax2.tick_params(axis = 'x', labelsize = txt_fontsize)
    ax2.set_ylabel(r'Probability (PMF) []', size = txt_fontsize)
    ax2.tick_params(axis = 'y', labelsize = txt_fontsize)

    fig_label = ('plots/orbital mechanics/alpha/alpha_PMF-' + ant_alf_hist_str_title[ind_hist] + '.png')
    fig0.savefig(fig_label)

