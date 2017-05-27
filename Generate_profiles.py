""""
This Python script makes several synthetic profiles of salt marsh and tidal flat environments.
It then converts them to .csv files.
These .csv files need to be opened and saved in Excel


UNLESS IT CAN DIRECTLY MAKE THEM INTO .xyz files!!!!"""


#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')


#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from matplotlib import cm
from pylab import *
import functools
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools as itt
import numpy as np
from osgeo import gdal, osr
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab

#------------------------------------------------------------------
#2. Generate variables

# The spring tidal range [m]
TR_spring = np.array([1., 3., 7.]) # 7
# The neap tidal range [m]. We will probably not use this
TR_neap = TR_spring/2

# The elevation budget: relative elevation of the lowest proint of the ideal marsh platform
# It is expressed as a fraction of TR_spring. This idea comes from D. Cahoon
Z_budget = np.array([-0.4, 0.6]) # 5

# The scarp height [m]
Scarp_height = np.array([0, 0.25, 1, 2.5]) # 6
#Scarp_height = np.array([0, 0.25, 1, 2.5]) # 6

# The marsh platform slope [m/m]
SM_slope = 0.00005
# The tidal flat slope [m/m]
TF_slope = np.array([0.5, 2.]) *10**(-4) # 5

# Length of the marsh platform [m]
SM_length = 350
# Length of the Transition [m]
Tr_length = 100
# Length of the Shallow tidal flat [m]
STF_length = 350
# Length of the tidal flat [m]
TF_length = 19250

Frame_length = SM_length + Tr_length + STF_length + TF_length - 50

# X,Y resolutions [m]
step_large = 50
step_medium = 5
step_small = 2

# Y size
Frame_width = 5000

#------------------------------------------------------------------
#------------------------------------------------------------------
# Functions

# Generate 1D profile frame
def Profile_frames_1D (SM_length, Tr_length, STF_length, TF_length, step_large, step_medium, step_small):
    SM_1D = np.zeros((2,SM_length/step_medium),dtype=np.float)
    Tr_1D = np.zeros((2,Tr_length/step_small),dtype=np.float)
    STF_1D = np.zeros((2,STF_length/step_medium),dtype=np.float)
    TF_1D = np.zeros((2,TF_length/step_large),dtype=np.float)

    SM_1D [0, :] = np.arange(0, SM_length, step_medium)
    Tr_1D [0, :] = np.arange(SM_length, SM_length+Tr_length, step_small)
    STF_1D [0, :] = np.arange(SM_length+Tr_length, SM_length+Tr_length+STF_length, step_medium)
    TF_1D [0, :] = np.arange(SM_length+Tr_length+STF_length, SM_length+Tr_length+STF_length+TF_length, step_large)

    return SM_1D, Tr_1D, STF_1D, TF_1D


# Calculate the platform profile
def SM_profile (Profile_1D, TR_spring, Z_budget, SM_slope, SM_length, step_medium):
    Start = 0
    Stop = SM_length/step_medium
    for x in range(Start, Stop):
        Profile_1D [1, x] = SM_slope*(SM_length-x*step_medium) + (TR_spring*Z_budget)/2

    return Profile_1D


# Calculate the Transition profile
def Tr_profile (Profile_1D_prolongated, Profile_1D, Scarp_height, SM_length, Tr_length, step_small):
    Start = 1
    Stop = int(Start + Tr_length/step_small) - 1
    Profile_1D_prolongated [1, 0] = Profile_1D [1, -1]
    x0 = (Start + Stop)/2
    L = Scarp_height
    if L<=0.5:
        k = 0.1
    else:
        k = 10

    for x in range(Start, Stop):
        Profile_1D_prolongated [1, x] = L / (1+np.exp(-k*(x0-x))) - (L-Profile_1D_prolongated [1, 0])

    return Profile_1D_prolongated


# Calculate the Shallow tidal flat profile
def STF_profile (Profile_1D_prolongated, Profile_1D, TF_slope, SM_length, Tr_length, STF_length, TF_length, step_medium):
    Start = 1
    Stop = int(Start + STF_length/step_medium) - 1
    Profile_1D_prolongated [1, 0] = Profile_1D [1, -1]
    for x in range(Start, Stop):
        Profile_1D_prolongated [1, x] = Profile_1D_prolongated [1, x-1] - TF_slope*step_medium

    return Profile_1D_prolongated


# Calculate the Tidal flat profile
def TF_profile (Profile_1D_prolongated, Profile_1D, TF_slope, TF_length, step_large):
    Start = 1
    Stop = int(Start + TF_length/step_large) - 1
    Profile_1D_prolongated [1, 0] = Profile_1D [1, -1]
    for x in range(Start, Stop):
        Profile_1D_prolongated [1, x] = Profile_1D_prolongated [1,x-1] - TF_slope*step_large

    return Profile_1D_prolongated



def Add_Y_values (Profile_1D):
    #Add the values
    Y_values = np.array([0*np.ones(len(Profile_1D[0,:]),dtype=np.float)])
    Profile_1D = np.concatenate((Y_values, Profile_1D),axis = 0)

    # Rearrange the values
    Profile_1D[0,:] = Profile_1D[1,:]
    Profile_1D[1,:] = 0

    return Profile_1D


def Twodimensionize (Profile_1D, Frame_width, step):
    Profile_1D_copy = np.copy(Profile_1D)
    for i in range(1,int(Frame_width/step)+1):
        Profile_1D_copy [1,:]=i*step
        Profile_1D = np.concatenate((Profile_1D, Profile_1D_copy),axis = 1)

    return Profile_1D



#------------------------------------------------------------------
#------------------------------------------------------------------
# Initiate profiles
SM_1D, Tr_1D, STF_1D, TF_1D = Profile_frames_1D (SM_length, Tr_length, STF_length, TF_length, step_large, step_medium, step_small)


#------------------------------------------------------
#------------------------------------------------------
# Setup the figure



#----------------------------------------------------------
#----------------------------------------------------------

for TR in range(len(TR_spring)):
    for z0 in range(len(Z_budget)):
        for h in range(len(Scarp_height)):
            for s in range(len(TF_slope)):
                # Initiate profiles
                SM_1D, Tr_1D, STF_1D, TF_1D = Profile_frames_1D (SM_length, Tr_length, STF_length, TF_length, step_large, step_medium, step_small)

                # Marsh elevations
                SM_1D = SM_profile (SM_1D, TR_spring[TR], Z_budget[z0], SM_slope, SM_length, step_medium)
                # Transition elevations
                Tr_1D = Tr_profile (Tr_1D, SM_1D, Scarp_height[h], SM_length, Tr_length, step_small)
                # Shallow tidal flat elevations
                STF_1D = STF_profile (STF_1D, Tr_1D, TF_slope[s], SM_length, Tr_length, STF_length, TF_length, step_medium)
                # Tidal flat elevations
                TF_1D = TF_profile (TF_1D, STF_1D, TF_slope[s], TF_length, step_large)

                # Make profiles
                Profile_1D = np.concatenate((SM_1D, Tr_1D),axis = 1)
                Profile_1D = np.concatenate((Profile_1D, STF_1D),axis = 1)
                Profile_1D = np.concatenate((Profile_1D, TF_1D),axis = 1)

                # Plot Profiles

                TR_value = TR_spring[TR]
                z0_value = Z_budget[z0]
                h_value = Scarp_height[h]
                s_value = TF_slope[s]


                fig=plt.figure(1, facecolor='White',figsize=[30,15])
                ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)
                ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
                ax.set_ylabel('Z - Elevation (m.a.m.s.l.)', fontsize = 18)
                ax2 = plt.axes([.47, .47, .4, .4])
                ax2.set_title('This is a zoom!', fontsize = 22)
                ax2.set_xlabel('X - Distance from origin (m)', fontsize = 12)
                ax2.set_ylabel('Z - Elevation (m.a.m.s.l.)', fontsize = 12)
                ax.grid(True)
                ax2.grid(True)

                ax.set_xlim (xmin = 0, xmax = 20100)
                ax.set_ylim (ymin = -6, ymax = TR_spring[TR]/2 + 0.5)
                ax2.set_xlim(xmin = SM_length-75, xmax = SM_length+Tr_length+STF_length+75)
                ax2.set_ylim (ymin = -TR_spring[TR]/2 - 0.5, ymax = TR_spring[TR]/2 + 0.5)

                #Filling the tide levels in the plot
                ax.fill_between (np.arange(0,Frame_length), -20, -TR_spring[TR]/2, alpha = 0.1)
                ax.fill_between (np.arange(0,Frame_length), -20, -TR_neap[TR]/2, alpha = 0.1)
                ax.fill_between (np.arange(0,Frame_length), -20, 0, alpha = 0.1)
                ax.fill_between (np.arange(0,Frame_length), -20, TR_neap[TR]/2, alpha = 0.1)
                ax.fill_between (np.arange(0,Frame_length), -20,TR_spring[TR]/2, alpha = 0.1)

                ax.set_title('Profile for TR=%g, Zo=%g, H=%g, S=%g' % (TR_value, z0_value, h_value, s_value), fontsize = 22)
                Map = ax.plot(Profile_1D[0,:], Profile_1D[1,:], color=plt.cm.gist_heat(s*30))
                Map2 = ax2.plot(Profile_1D[0,:], Profile_1D[1,:], color=plt.cm.gist_heat(s*30))
                plt.savefig('Figures/Profile_TR%g_Zo%g_H%g_S%g.png' % (TR_value, z0_value, h_value, s_value))


                #Profile_1D = Profile_1D.T

                np.savetxt('Profile_TR%g_Zo%g_H%g_S%g.txt' % (TR_value, z0_value, h_value, s_value), Profile_1D)

                #----------------------------------
                # IF YOU ONLY NEED PROFILES, STOP THERE


                """# Add y values for the y=0 line
                SM_1D = Add_Y_values(SM_1D)
                Tr_1D = Add_Y_values(Tr_1D)
                STF_1D = Add_Y_values(STF_1D)
                TF_1D = Add_Y_values(TF_1D)

                # Twodimensionize
                SM_1D = Twodimensionize (SM_1D, Frame_width, step_medium)
                Tr_1D = Twodimensionize (Tr_1D, Frame_width, step_small)
                STF_1D = Twodimensionize (STF_1D, Frame_width, step_medium)
                TF_1D = Twodimensionize (TF_1D, Frame_width, step_large)

                #Now aggregate the high and low resolution arrays
                Profile_1D = np.concatenate((SM_1D, Tr_1D),axis = 1)
                Profile_1D = np.concatenate((Profile_1D, STF_1D),axis = 1)
                Profile_1D = np.concatenate((Profile_1D, TF_1D),axis = 1)

                #fig=plt.figure(2, facecolor='White',figsize=[30,30])
                #ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1)
                #ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
                #ax.set_ylabel('Y - Distance from origin (m)', fontsize = 18)
                #ax.scatter(Profile_1D[0,:], Profile_1D[1,:], c = Profile_1D[2,:], linewidth = 0, alpha = 0.5)
                #ax.set_xlim (xmin = 0, xmax = 1000)
                #ax.set_ylim (ymin = 0, ymax = 1000)
                #plt.savefig('Figures/Profilemap_5km.png')


                Profile_1D = Profile_1D.T
                # Name the files correctly
                TR_value = TR_spring[TR]
                z0_value = Z_budget[z0]
                h_value = Scarp_height[h]
                s_value = TF_slope[s]

                np.savetxt('Points_TR%g_Zo%g_H%g_S%g.xyz' % (TR_value, z0_value, h_value, s_value), Profile_1D)"""
