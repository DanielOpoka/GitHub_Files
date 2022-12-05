# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:18:27 2022

@author: d.opoka
"""

"""
##############################################################################
                            CLIP MAKER OF TUMBLING
##############################################################################
"""

'Code used to create 3D renderings of the CREME Tumbling'

'Install OpenCV in the Conda navigator/Linux Terminal with the following line:'
'conda install opencv'

import numpy as np  
import matplotlib.pyplot as plt  
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.pyplot import figure
import glob
import os
import cv2

# %%
"""
------------------------------------------------------------------------------
                           v2 - ASTRODYNAMICS OF TUMBLING
------------------------------------------------------------------------------
"""

"""
Trigonometric functions are in degrees.

The Earth-Centered Inertial (ECI) frame is used.

The simulation begins on the most recent Vernal Equinox: 
    20/03/2022 - 15:33:00
"""

'-------------------------------- FUNCTIONS ----------------------------------'
# Reference frame rotation about its X-axis0is
def R1_RX(angle_X):
    return np.array([[1, 0, 0],
                     [0, np.cos(angle_X), np.sin(angle_X)],
                     [0, -np.sin(angle_X), np.cos(angle_X)]])

# Reference frame rotation about its Y-axis0is
def R1_RY(angle_Y):
    return np.array([[np.cos(angle_Y), 0, -np.sin(angle_Y)],
                     [0, 1, 0],
                     [np.sin(angle_Y), 0, np.cos(angle_Y)]])

# Reference frame rotation about its Z-axis0is
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

# Converts degrees to radians
def d2r(degrees):
    radians = degrees * (np.pi / 180)
    return radians

# Converts radians to degrees
def r2d(radians):
    degrees = radians * (180 / np.pi)
    return degrees

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

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

def drawLine(axis, ind_int, ind_end, ind_color, ind_zorder):
    arrow_prop_dict = dict(mutation_scale=20, arrowstyle='-', shrinkA=0, shrinkB=0)
    a = Arrow3D([ind_int[0], ind_end[0]],
                [ind_int[1], ind_end[1]],
                [ind_int[2], ind_end[2]],
                **arrow_prop_dict, color = ind_color, zorder = ind_zorder)
    axis.add_artist(a)   

# For drawing the axes arrows of a Reference Frame
def drawAxes(axis, ind_int, ind_end, axes_label, ind_color, ind_zorder):
    arrow_prop_dict = dict(mutation_scale=20, arrowstyle='-|>', shrinkA=0, shrinkB=0)
    
    for ind in [0, 1, 2]:
        a = Arrow3D([ind_int[0], ind_end[ind][0]], 
                    [ind_int[1], ind_end[ind][1]], 
                    [ind_int[2], ind_end[ind][2]], 
                    **arrow_prop_dict, color = ind_color, zorder = ind_zorder)
        axis.add_artist(a)                          # X-axis arrow
        axis.text(ind_end[ind][0], ind_end[ind][1], ind_end[ind][2], 
                  '$\hat{' + str(axes_label[ind]) + '}$', color = ind_color)

def create_dir(dir):
  if not os.path.exists(dir):
    os.makedirs(dir)
    print("Created Directory : ", dir)
  else:
    print("Directory already existed : ", dir)
  return dir

'--------------------- REFERENCE FRAME FOR BODY OF SATELLITE -----------------'
# Earth Parameters
ear_rad = 6378.14 * 1000                    # Radius of the Earth [m]
ear_day = 24 * 60 * 60                      # Day on Earth [s]
ear_ome = d2r(360.00) / ear_day             # Angular velocity of the Earth [rad / s]

# Orbital Elements (from orbit 2 - alt = 600 km)
orb_OmA = d2r(00.00)                   # Longitude of ascending node [rad]
orb_inc = d2r(97.79)                   # Inclination [rad]            
orb_omP = d2r(00.00)                   # Argument of periapsis/perigee [rad]
orb_ecc = 00.00                             # Eccentricity []
orb_sma = ear_rad + 600 * 1000              # Semi-major axis0is (a) [m]
orb_Tnu = d2r(00.00)                   # True Anomaly [rad] (at time 0 s)

# Orbit Parameters
# Semi-minor axis0is (b) [m]
orb_smb = np.sqrt((1 - orb_ecc ** 2) * orb_sma ** 2)
orb_per = 96.69 * 60                        # Orbit Period/Duration [s]
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
RR_RE = np.matmul(R1_RZ(d2r(90) + orb_Tnu), RO_RE)

# Angles from Satellite Reference Frame to Satellite Body Frame
'These are a function of data from the on-board sensors'
sat_alf_int = d2r(+00.00)                        # Initial alfa_body (X) [rad]
sat_bet_int = d2r(+90.00)                        # Initial beta_body (Y) [rad]
sat_gam_int = d2r(+00.00)                        # Initial gamma_body (Z) [rad]
'''
# Initial Angles according to the Satellite Body Frame
sat_alf = d2r(00.00)
sat_bet = d2r(00.00)
sat_gam = d2r(00.00)
'''
# Satellite Angular Velocity according to the Satellite Body Frame
'These are a function of data from the on-board sensors'
sat_alf_dot = d2r(00.00) / 1                    # delta^dot_body (X) [rad / s]
sat_bet_dot = d2r(00.00) / 1                    # epsilon^dot_body (Y) [rad / s]
sat_gam_dot = d2r(02.00) / 1                    # zeta^dot_body (Z) [rad / s]

# Translation Matrix from Satellite Reference Frame to the Satellite Body Frame [m]
'Vector from the CoG to the origin of the Body'
CgOb = np.array([0, 0, -0.15])
# Direct Transformation Matrix from Satellite Reference Frame to Satellite Body Frame
# -Z_R is aligned with X_O, pointing towards the sun
RB_RE = R1_RY(np.deg2rad(-90.00))

# Current position of the Origin of the Satellite 
curr_Ob = np.matmul(CgOb, RB_RE)

'-------------------- REFERENCE FRAME FOR GROUND SEGMENT ---------------------'
gse_alf = d2r(1.444209)                     # Longitude (E) [rad]
gse_bet = d2r(43.604652)                    # Latitude (N) [rad]
gse_alt = ear_rad + 130                     # Altitude [m]

# Angle of the GS at time 0 s (Vernal Equinox, Z) [rad]
# Earth's X-axis pointing to the Sun at 0°00' North, 51°08' West at time 0s
# Source: https://www.timeanddate.com/worldclock/sunearth.html?iso=20220320T1532
gse_tht = d2r(51 + (8 / 60)) + gse_alf     

# Direct Transformation Matrix from Ground Segment RF to ECI
# Rot with beta_gse (Z) -> Rot with tht_gse (Y) -> Rot with 90° (Z)
# X -> East, Y -> North, and Z points towards Space
RGS_RE = np.matmul(R1_RZ(d2r(90)),
                   np.matmul(R1_RY(gse_bet), 
                             R1_RZ(gse_tht)))
# Direct Transformation Matrix from ECI TO Ground Segment RF
RE_RGS = np.linalg.inv(RGS_RE)

# The Ground Segment has a visibility up to the Clearance/Horizon Angle. 
gse_clr_ang_E = d2r(90 - 15)                 # Ground segment Horizon Angle (East) [rad]
gse_clr_ang_W = d2r(90 - 15)                 # Ground segment Horizon Angle (West) [rad]

'------------------------ SATELLITE ANTENA GAIN MODEL ------------------------'
ant_ran = (1626.24 * 1000) / 2              # Antenna Range (Edge-to-edge distance) [m] 
# Antenna FoV Line Coordinates [m]
ant_FoV = np.array([[+00.00, +ant_ran, +00.20],
                    [+00.00, -ant_ran, +00.20]])                 
# ...

'--------------------------------- SIMULATION --------------------------------'
# Parameters
sim_tini = 0                              # Beginning of time [s]
sim_tend = 60                            # Ending of time [s]
# sim_tend = 600                               # Ending of time [s]
fps_nb = 60                                 # FPS for video [#]
sim_tele = fps_nb * 5                      # Number of elements for time []
# Time [s]
sim_time = np.linspace(sim_tini, sim_tend, sim_tele)
sim_tstp = sim_time[1] - sim_time[0]        # Time Step [s]

# Satellite Orbit Dynamics
if (orb_ecc == 00.00):
    # For circular orbit
    orb_Tnu_dot = orb_vel / orb_sma         # Orbit angular velocity [rad / s]
else:
    # For elliptic orbit
    # Orbit distance from orbitting Focus to Satellite [m]
    orb_dis = (orb_sma * (1 - orb_ecc ** 2)) / (1 + orb_ecc * np.cos(orb_Tnu))
    # Orbit angular velocity [rad / s] 
    orb_Tnu_dot = (2 * np.pi * orb_sma * orb_smb) / (orb_per * orb_dis)

orb_pos = np.zeros([3, sim_tele])           # Satellite's position in orbit [m]
gse_pos = np.zeros([3, sim_tele])           # Ground Segment's position in ECI [m]
gse_ant_pos = np.zeros([3, sim_tele])       # Position of Point between GS and Antenna FoV [m]
gse_ant_dis = np.zeros(sim_tele)            # Distance between GS and Antenna FoV [m]

# For images and clips
if (sim_tend == orb_per):
    time_end = 'orb_per'
elif (sim_tend == ear_day):
    time_end = '24hr'
else:
    time_end = str(sim_tend) + 's'

clip_name = ('ad=' + str(round(r2d(sat_alf_dot))).rjust(2, '0') + '-' + 
             'bd=' + str(round(r2d(sat_bet_dot))).rjust(2, '0') + '-' + 
             'gd=' + str(round(r2d(sat_gam_dot))).rjust(2, '0') + 
             '_time=' + time_end)

create_dir('images4clips/' + clip_name)  

# Clearing Images for Orbit Clip
folder_path0 = 'images4clips/' + clip_name + '/Orbit'
create_dir('images4clips/' + clip_name + '/Orbit') 
list_images = os.listdir(folder_path0)
for image in list_images:
    if image.endswith('.png'):
        os.remove(os.path.join(folder_path0, image))

# Clearing Images for Orbit Clip
folder_path1 = 'images4clips/' + clip_name + '/Satellite'
create_dir('images4clips/' + clip_name + '/Satellite') 
list_images = os.listdir(folder_path0)
for image in list_images:
    if image.endswith('.png'):
        os.remove(os.path.join(folder_path0, image))

# Plotting Parameters
tit_fontsize = 18
txt_fontsize = 14
figsize_X = 8
figsize_Y = 8
dpi_ind = 75

# Determining the number of trailing zeros of sim_tele
N = 0
while (sim_tele / (10 ** N)) >= (10 ** 1): 
    N += 1
nb_zeros = N + 1

DGA = np.zeros(sim_tele)            # Distance between GS and Antenna FoV (DGA) [m]
clr_tht = np.zeros(sim_tele)        # Clearange Angle between GS and Satellite [°]
clr_DGA = np.zeros(sim_tele)        # DGA during clearance [m]
ant_clr_alf = np.zeros(sim_tele)    # Antenna Angle to GS [°]

for sim_iter, curr_time in enumerate(sim_time): 
    print('sim_iter = ' + str(sim_iter + 1) + ' / ' + str(sim_tele))
    
    # For elliptic orbit
    if (orb_ecc != 00.00):
        # Orbit distance from orbitting Focus to Satellite [m]
        orb_dis = (orb_sma * (1 - orb_ecc ** 2)) / (1 + orb_ecc * np.cos(orb_Tnu))
        # Orbit angular velocity [rad / s] 
        orb_Tnu_dot = (2 * np.pi * orb_sma * orb_smb) / (orb_per * orb_dis ** 2)
       
    # Updating True Anomaly
    orb_Tnu = orb_Tnu_dot * curr_time
    
    'Satellite Position in Orbital Plane' 
    # Current Satellite position in the Orbital Plane
    curr_pos = np.array([orb_sma * np.cos(orb_Tnu), orb_smb * np.sin(orb_Tnu), 0])
    curr_pos = np.matmul(curr_pos, RO_RE)
    
    # Updating DTM from Satellite RF to the ECI (in Orbital Plane)
    RR_RE = np.matmul(R1_RZ(d2r(90) + orb_Tnu), RO_RE)
    
    # Updating Satellite Body Angles
    sat_alf = sat_alf_dot * curr_time
    sat_bet = sat_bet_dot * curr_time
    sat_gam = sat_gam_dot * sim_tstp
    # Updating DTM from Satellite Body RF to the ECI
    RB_RE = np.matmul(R1_RX(sat_alf), 
                      np.matmul(R1_RY(sat_bet), 
                                np.matmul(R1_RZ(sat_gam), RB_RE)))
    # Updating Origin of the Satellite Body RF
    curr_Ob = np.matmul(CgOb, RB_RE)
    
    'Antenna Attitude'
    # Updating the Field of View direction
    curr_ant = np.zeros([2, 3])
    curr_ant[0] = np.matmul(ant_FoV[0], RB_RE) + curr_pos + curr_Ob
    curr_ant[1] = np.matmul(ant_FoV[1], RB_RE) + curr_pos + curr_Ob
    
    'Ground Segment Position in ECI'
    curr_gse = np.array([np.cos(gse_bet) * np.cos(gse_tht),
                         np.cos(gse_bet) * np.sin(gse_tht),
                         np.sin(gse_bet)]) * gse_alt
    # Updating the Ground Segment theta angle
    gse_tht += ear_ome * sim_tstp
    # Updating DTM from Ground Segment RF to ECI
    RGS_RE = np.matmul(R1_RZ(d2r(90)),
                   np.matmul(R1_RY(gse_bet), 
                             R1_RZ(gse_tht)))
    # Updating DTM from ECI to Ground Segment RF
    RE_RGS = np.linalg.inv(RGS_RE)
    # Current position of Satellite according to Ground Segment RF
    curr_sat_gse = np.matmul(curr_pos, RE_RGS) - np.matmul(curr_gse, RE_RGS)
            
    # Updating the current Point of Perpecndicular between Antenna FoV and GS
    gse_ant_point = PerpLinePoint(curr_ant[0], curr_gse, curr_ant[1])
    # Distance between GS and Antenna FoV Point (DGA) [m]
    DGA[sim_iter] = magVec(curr_gse, gse_ant_point)
    # Distance between Antenna FoV Point and Satellite (DGA) [m]
    DAS = magVec(curr_pos, gse_ant_point)
    
    # Determining the Clearance/Visibility angle between Satellite and GS
    if (curr_sat_gse[2] > 0):
        clr_alt = curr_sat_gse[2] 
        clr_rad = np.sqrt(curr_sat_gse[0] ** 2 + curr_sat_gse[1] ** 2)
        curr_clr_tht = np.arctan(clr_rad / clr_alt)
        
        # If on the East-Side
        if (((curr_sat_gse[0] > 0) and (curr_sat_gse[1] > 0)) or 
            ((curr_sat_gse[0] > 0) and (curr_sat_gse[1] < 0))):
            if (curr_clr_tht < gse_clr_ang_E):
                clr_tht[sim_iter] = r2d(curr_clr_tht)
                clr_DGA[sim_iter] = DGA[sim_iter]
                ant_clr_alf[sim_iter] = r2d(np.arctan(DGA[sim_iter] / DAS))
                
        # If on the West-Side
        if (((curr_sat_gse[0] < 0) and (curr_sat_gse[1] > 0)) or 
            ((curr_sat_gse[0] < 0) and (curr_sat_gse[1] < 0))):
            if (curr_clr_tht < gse_clr_ang_W):
                clr_tht[sim_iter] = r2d(curr_clr_tht)
                clr_DGA[sim_iter] = DGA[sim_iter]
                ant_clr_alf[sim_iter] = r2d(np.arctan(DGA[sim_iter] / DAS))
    
    # Updating the current Point of Perpecndicular Line between Antenna FoV and GS
    gse_ant_point = PerpLinePoint(curr_ant[0], curr_gse, curr_ant[1])
    
    for ind_XYZ in [0, 1, 2]:
        orb_pos[ind_XYZ][sim_iter] = curr_pos[ind_XYZ]
        gse_pos[ind_XYZ][sim_iter] = curr_gse[ind_XYZ]
        gse_ant_pos[ind_XYZ][sim_iter] = gse_ant_point[ind_XYZ]
    
    gse_ant_dis[sim_iter] = magVec(curr_gse, gse_ant_point)
    

    # When changing current Reference Frame (RF) to ECI RF: 
    # np.matmul(input coord, RF_RE) + origin of new RF
    
    # '%%
    """
    --------------------------------------------------------------------------
                                       3D Plots
    --------------------------------------------------------------------------
    """
    plt.close('all')
    nb_fig = 0
    
    ear_vis = False
    
    '--------------------------------- ORBIT -------------------------------------'
    nb_fig += 1
    fig0 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi=dpi_ind)
    plt.ioff()
    axis0 = plt.axes(projection='3d')
    axis0.set_box_aspect([1,1,1])
    axis0.view_init(45, 45)
    
    mag = ear_rad                               # Magnitude for axis-arrows [m]
    mult = 1.5                                  # Multiplier for magnitude []
    ind_zorder = 10
    
    'Plotting Axes'
    # Plotting Axes-Arrows for ECI (R_{E})
    arr_pos_int = np.zeros(3)
    arr_pos_end = np.identity(3) * mag * mult
    axes_label = ['X_{E}', 'Y_{E}', 'Z_{E}']
    drawAxes(axis0, arr_pos_int, arr_pos_end, axes_label, 'k', ind_zorder)
    
    # Plotting Axes-Arrows for Orbital Plane (R_{O})
    arr_pos_int = np.zeros(3)
    arr_pos_end = np.matmul(RO_RE, np.identity(3) * mag * mult)
    axes_label = ['X_{O}', 'Y_{O}', 'Z_{O}']
    drawAxes(axis0, np.zeros(3), arr_pos_end, axes_label, 'b', ind_zorder)
    
    # Plotting Axes-Arrows for Satellite Reference Frame (R_{R})
    ind = -1
    arr_pos_int = np.array([orb_pos[0][ind], orb_pos[1][ind], orb_pos[2][ind]])
    arr_pos_int = curr_pos
    arr_pos_end = np.matmul(RR_RE, np.identity(3) * mag * (mult / 3)) + curr_pos
    axes_label = ['X_{R}', 'Y_{R}', 'Z_{R}']
    drawAxes(axis0, arr_pos_int, arr_pos_end, axes_label, 'g', ind_zorder)
    
    # Plotting Axes-Arrows for Satellite Body Reference Frame (R_{B})
    arr_pos_int = curr_pos
    arr_pos_end = np.matmul(RB_RE, np.identity(3) * mag * (mult / 3)) + curr_pos + curr_Ob
    axes_label = ['X_{B}', 'Y_{B}', 'Z_{B}']
    drawAxes(axis0, arr_pos_int, arr_pos_end, axes_label, 'c', ind_zorder)
    
    # Plotting Axes-Arrows for Ground Segment Reference Frame (R_{GS})
    ind = -1
    arr_pos_int = curr_gse
    arr_pos_end = np.matmul(RGS_RE, np.identity(3) * mag * (mult / 3)) + curr_gse
    axes_label = ['X_{GS}', 'Y_{GS}', 'Z_{GS}']
    drawAxes(axis0, arr_pos_int, arr_pos_end, axes_label, 'r', ind_zorder)
    
    'Plotting Elements'
    
    # Earth
    if (ear_vis):
        sph_u, sph_v = np.mgrid[0:2*np.pi:36j, 0:np.pi:36j]
        sph_x = ear_rad * np.cos(sph_u)*np.sin(sph_v)
        sph_y = ear_rad * np.sin(sph_u)*np.sin(sph_v)
        sph_z = ear_rad * np.cos(sph_v)
        axis0.plot_wireframe(sph_x, sph_y, sph_z, zorder = 1, color="g")
    
    # Orbit of Satellite
    axis0.plot(orb_pos[0][:], orb_pos[1][:], orb_pos[2][:], 'b', zorder = ind_zorder)
    # Satellite's Starting Position
    ind = 0
    coord3D = [orb_pos[0][ind], orb_pos[1][ind], orb_pos[2][ind]]
    axis0.plot(coord3D[0], coord3D[1], coord3D[2], 'bx', zorder = ind_zorder)
    
    # Satellite's Final Position
    coord3D = [curr_pos[0], curr_pos[1], curr_pos[2]]
    axis0.plot(coord3D[0], coord3D[1], coord3D[2], 'bo', zorder = ind_zorder)
    # Antenna FoV
    lin_pos_int = curr_ant[0]
    lin_pos_end = curr_ant[1]
    drawLine(axis0, lin_pos_int, lin_pos_end, 'crimson', ind_zorder)
    # GS and Antenna FoV Perp Line 2nd Point
    axis0.plot(gse_ant_point[0], gse_ant_point[1], gse_ant_point[2], 'o', 
               color = 'crimson', zorder = ind_zorder)
    
    # Ground Segment's Path
    axis0.plot(gse_pos[0][:], gse_pos[1][:], gse_pos[2][:], 'red', zorder = ind_zorder)
    # Ground Segment's Starting Position
    ind = 0
    coord3D = [gse_pos[0][ind], gse_pos[1][ind], gse_pos[2][ind]]
    axis0.plot(coord3D[0], coord3D[1], coord3D[2], 'rx', zorder = ind_zorder)
    # Ground Segment's Final Position
    coord3D = curr_gse
    axis0.plot(coord3D[0], coord3D[1], coord3D[2], 'ro', zorder = ind_zorder)
    
    'Labels'
    title_label = (r'$\dot{\alpha}$ (X) = ' + str(round(r2d(sat_alf_dot), 2)) + ' [°/s], ' +
                   r'$\dot{\beta}$ (Y) = ' + str(round(r2d(sat_bet_dot), 2)) + ' [°/s], ' +
                   r'$\dot{\gamma}$ (Z) = ' + str(round(r2d(sat_gam_dot), 2)) + ' [°/s] \n ' +
                   'Time = ' + str(round(curr_time, 2)) + ' s')
    plt.title(title_label, fontsize = txt_fontsize)
    axis0.set_xlabel('$X-axis$ [m]', fontsize = txt_fontsize)
    axis0.set_ylabel('$Y-axis$ [m]', fontsize = txt_fontsize)
    axis0.set_zlabel('$Z-axis$ [m]', fontsize = txt_fontsize)
    
    'XYZ Limits'
    axis0.set_xlim(0 - mag * mult, 0 + mag * mult)
    axis0.set_ylim(0 - mag * mult, 0 + mag * mult)
    axis0.set_zlim(0 - mag * mult, 0 + mag * mult)
    
    fig_label = ('images4clips/' + clip_name + '/Orbit/Clip' + 
                 str(sim_iter).rjust(nb_zeros, '0') + '.png')
    fig0.savefig(fig_label)
    
    '------------------------------------ SATELLITE ------------------------------'
    
    nb_fig += 1
    fig1 = figure(nb_fig, figsize=(figsize_X, figsize_Y), dpi=dpi_ind)
    plt.ioff()
    axis1 = plt.axes(projection='3d')
    axis1.set_box_aspect([1,1,1])
    axis1.view_init(45, 45)
    
    mag = 0.30                               # Magnitude for axis-arrows [m]
    mult = 1.0                               # Multiplier for magnitude []
    ind_zorder = 10
    
    'Plotting Axes'
    # Plotting Axes-Arrows for Satellite Reference Frame (R_{R})
    arr_pos_int = np.zeros(3)
    arr_pos_end = np.matmul(RR_RE, np.identity(3) * mag * mult)
    axes_label = ['X_{R}', 'Y_{R}', 'Z_{R}']
    drawAxes(axis1, arr_pos_int, arr_pos_end, axes_label, 'g', ind_zorder)
    
    # Plotting Axes-Arrows for Satellite Body Reference Frame (R_{B})
    arr_pos_int = curr_Ob
    arr_pos_end = np.matmul(RB_RE, np.identity(3) * mag * mult) + curr_Ob
    axes_label = ['X_{B}', 'Y_{B}', 'Z_{B}']
    drawAxes(axis1, arr_pos_int, arr_pos_end, axes_label, 'c', ind_zorder)
    
    'Plotting Elements'
    # Center of Gravity
    axis1.plot(0, 0, 0, 'go', zorder = ind_zorder)
    axis1.text(0, 0, 0, 'CoG', color = 'g')
    
    axis1.plot(curr_Ob[0], curr_Ob[1], curr_Ob[2], 'co', zorder = ind_zorder)
    axis1.text(curr_Ob[0], curr_Ob[1], curr_Ob[2], 'Ob', color = 'c')
    
    # Satellite Body
    coord_body = np.array([[+00.05, +00.05, +00.00],
                           [-00.05, +00.05, +00.00],
                           [-00.05, -00.05, +00.00],
                           [+00.05, -00.05, +00.00],
                           [+00.05, +00.05, +00.30],
                           [-00.05, +00.05, +00.30],
                           [-00.05, -00.05, +00.30],
                           [+00.05, -00.05, +00.30]])
    for ind0 in range(4):
        if (ind0 == 3):
            ind1 = 0
        else:
            ind1 = ind0 + 1
            
        lin_pos_int = np.matmul(coord_body[ind0], RB_RE) + curr_Ob
        lin_pos_end = np.matmul(coord_body[ind1], RB_RE) + curr_Ob
        drawLine(axis1, lin_pos_int, lin_pos_end, 'k', ind_zorder)
        
        lin_pos_int = np.matmul(coord_body[ind0 + 4], RB_RE) + curr_Ob
        lin_pos_end = np.matmul(coord_body[ind1 + 4], RB_RE) + curr_Ob
        drawLine(axis1, lin_pos_int, lin_pos_end, 'k', ind_zorder)
        
        if (ind0 == 3):
            ind1 = ind0 + 1
        
        lin_pos_int = np.matmul(coord_body[ind0], RB_RE) + curr_Ob
        lin_pos_end = np.matmul(coord_body[ind1 + 3], RB_RE) + curr_Ob
        drawLine(axis1, lin_pos_int, lin_pos_end, 'k', ind_zorder)
    
    # Antenna FoV
    lin_pos_int = np.matmul(np.array([+00.00, +01.00, +00.20]), RB_RE) + curr_Ob
    lin_pos_end = np.matmul(np.array([+00.00, -01.00, +00.20]), RB_RE) + curr_Ob
    drawLine(axis1, lin_pos_int, lin_pos_end, 'crimson', ind_zorder)
    # Antenna Surfaces
    ant_surf1 = np.zeros([3, 100])
    ant_surf2 = np.zeros([3, 100])
    ant_surf_tht = np.linspace(0, d2r(360), 100)
    ant_surf_rad = 0.05
    ant_surf_01 = np.matmul(np.array([+00.00, +00.05, +00.20]), RB_RE)
    ant_surf_02 = np.matmul(np.array([+00.00, -00.05, +00.20]), RB_RE)
    
    for ind, elem in enumerate(ant_surf_tht):
        curr_surf0 = np.array([ant_surf_rad * np.cos(elem), 0, ant_surf_rad * np.sin(elem)])
        
        curr_surf1 = np.matmul(curr_surf0, RB_RE) + curr_Ob + ant_surf_01
        curr_surf2 = np.matmul(curr_surf0, RB_RE) + curr_Ob + ant_surf_02
        
        for ind_XYZ in [0, 1, 2]:
            ant_surf1[ind_XYZ][ind] = curr_surf1[ind_XYZ]
            ant_surf2[ind_XYZ][ind] = curr_surf2[ind_XYZ]
      
    axis1.plot(ant_surf1[0][:], ant_surf1[1][:], ant_surf1[2][:], 'crimson', zorder = ind_zorder)
    axis1.plot(ant_surf2[0][:], ant_surf2[1][:], ant_surf2[2][:], 'crimson', zorder = ind_zorder)      
        
    'Labels'
    title_label = (r'$\dot{\alpha}$ (X) = ' + str(round(r2d(sat_alf_dot), 2)) + ' [°/s], ' +
                   r'$\dot{\beta}$ (Y) = ' + str(round(r2d(sat_bet_dot), 2)) + ' [°/s], ' +
                   r'$\dot{\gamma}$ (Z) = ' + str(round(r2d(sat_gam_dot), 2)) + ' [°/s] \n ' +
                   'Time = ' + str(round(curr_time, 2)) + ' s')
    plt.title(title_label, fontsize = txt_fontsize)
    axis1.set_xlabel('$X-axis$ [m]', fontsize = txt_fontsize)
    axis1.set_ylabel('$Y-axis$ [m]', fontsize = txt_fontsize)
    axis1.set_zlabel('$Z-axis$ [m]', fontsize = txt_fontsize)
    
    'XYZ Limits'
    axis1.set_xlim(0 - mag * mult, 0 + mag * mult)
    axis1.set_ylim(0 - mag * mult, 0 + mag * mult)
    axis1.set_zlim(0 - mag * mult, 0 + mag * mult)
    
    
    fig_label = ('images4clips/' + clip_name + '/Satellite/Clip' + 
                 str(sim_iter).rjust(nb_zeros, '0') + '.png')
    fig1.savefig(fig_label)
    
# %%   
"""
------------------------------------------------------------------------------
                             SAVING VIDEOS
------------------------------------------------------------------------------
"""

def makeMP4(input_folder, fps_nb, output_folder):
    img_array = []
    
    print('Making Video of ' + input_folder + '...')
    for filename in sorted(glob.glob(input_folder)):
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width, height)
        img_array.append(img)
        
    fourcc = cv2.VideoWriter_fourcc('m','p','4','v')
    out = cv2.VideoWriter(output_folder, fourcc, fps_nb, size)
    
    for ind in range(len(img_array)):
        out.write(img_array[ind])
    
    out.release()
    print('Video Ready at ' + output_folder + '!')

input_folder = 'images4clips/' + clip_name + '/Orbit/*.png'
output_folder = 'clips/Orbit/Orbit_' + clip_name + '.mp4'
makeMP4(input_folder, fps_nb, output_folder)

input_folder = 'images4clips/' + clip_name + '/Satellite/*.png'
output_folder = 'clips/Satellite/Satellite_' + clip_name + '.mp4'
makeMP4(input_folder, fps_nb, output_folder)