"""
This Python script does two things:
1. It runs a batch of SWAN runs
2. It reads the SWAN outputs and writes them in a text file
"""

#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt


#----------------------------
# These are a bunch of functions

# This one reads the geometry of the nodes
def Read_geo (Input, Header_end, Data_end):
    f = open(Input + '.geo', 'r')
    lines = f.readlines()
    Data = np.zeros((3,Data_end-Header_end), dtype = np.float)

    i = 0
    for line in lines[Header_end:Data_end]:
        Data[0,i] = float(line [0:15])
        Data[1,i] = float(line [16:31])
        Data[2,i] = float(line [31:42])
        i = i + 1

    return Data


#----------------------------------------
# This one imports the outputs from the .mat file
def Read_mat (Input, Geometry):
    Data = io.loadmat(Input + '.mat')

    Hs = Data ['Hsig']
    Tb = Data ['TmBot']
    Urms = Data ['Urms']
    Tp = Data ['TPsmoo']
    Ub = Data ['Ubot']

    Output = np.concatenate((Geometry, Hs),axis=0)
    Output = np.concatenate((Output, Tb),axis=0)
    Output = np.concatenate((Output, Urms),axis=0)
    Output = np.concatenate((Output, Tp),axis=0)
    Output = np.concatenate((Output, Ub),axis=0)


    return Output

#----------------------------------------
# This one makes a profile out of the outputs, to we can pretent to have 1D data
def make_profile(Data, SM_length, Tr_length, STF_length, TF_length, Frame_width, step_large, step_medium, step_small):

    # Sort and clip the data to only have the central band of the domain
    Data= Data.T
    Sorted_Data = sorted(Data.tolist(), key = lambda a:a[0])
    Sorted_Data=np.asarray(Sorted_Data).T
    Clip = np.where(np.logical_and(Sorted_Data[1]<(Frame_width/2+200),Sorted_Data[1]>(Frame_width/2-200)))
    Profile = Sorted_Data[:,Clip[0]]

    """Now average these values in irregular bins. But not today"""
    #Start = 0
    #Stop = Start + SM_length
    #x_bins_SM = np.arange(Start, Stop, step_medium)

    #Start = Stop
    #Stop = Start + Tr_length
    #x_bins_Tr = np.arange(Start, Stop, step_small)

    #Start = Stop
    #Stop = Start + STF_length
    #x_bins_STF = np.arange(Start, Stop, step_medium)

    #Start = Stop
    #Stop = Start + TF_length
    #x_bins_TF = np.arange(Start, Stop, step_large)


    return Profile


#----------------------------
#The script starts here

#These are the dimensions of each section of the mesh. We don't use them yet
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


#------
#This is an example
#INPUT = "Input_TR1_Zo-0.4_H0.25_S0.0002"
#GEOMETRY = "Points_TR1_Zo-0.4_H0.25_S0.0002"
#OUTPUT = "Sim_TR1_Zo-0.4_H0.25_S0.0002"
#------


# The spring tidal range [m]
#TR_spring = ["1", "3", "7"]
TR_spring = ["7"]
# The elevation budget: relative elevation of the lowest proint of the ideal marsh platform
#Z_budget = ["-0.4", "0.6"]
Z_budget = ["0.6"]
# The scarp height [m]
#Scarp_height = ["0", "0.25", "1", "2.5"]
Scarp_height = ["2.5"]
# The tidal flat slope [m/m]
TF_slope = ["5e-05", "0.0002"]

for TR in TR_spring:
    for z0 in Z_budget:
        for h in Scarp_height:
            for s in TF_slope:
                INPUT = "Input_TR%s_Zo%s_H%s_S%s" % (TR, z0, h, s)
                GEOMETRY = "Points_TR%s_Zo%s_H%s_S%s" % (TR, z0, h, s)
                OUTPUT = "Sim_TR%s_Zo%s_H%s_S%s" % (TR, z0, h, s)

                # 1. Run SWAN
                os.system( "./swanrun -input " + INPUT)

                # 2.1. Read the .GEO file
                Data = Read_geo (GEOMETRY, 7, 14393)

                # 2.2. Add the results to the array
                Data = Read_mat (OUTPUT, Data)
                # This is to save all the data
                np.savetxt("Output_TR%s_Zo%s_H%s_S%s.txt" % (TR, z0, h, s),Data)

                # 3. Make a profile in the data array
                Profile = make_profile(Data, SM_length, Tr_length, STF_length, TF_length, Frame_width, step_large, step_medium, step_small)
                # This is to save a profile
                np.savetxt("Profile_TR%s_Zo%s_H%s_S%s.txt" % (TR, z0, h, s), Profile)




"""#-------------------------------
# Plot a profile
fig=plt.figure(1, facecolor='White',figsize=[30,30])
ax = plt.subplot2grid((2,2),(0,0),colspan=1, rowspan=1)
ax.set_title('Z', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Z (m)', fontsize = 18)
ax.plot(Profile[0,:], Profile[2,:])

ax = plt.subplot2grid((2,2),(0,1),colspan=1, rowspan=1)
ax.set_title('Hs', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Hs (m)', fontsize = 18)
ax.plot(Profile[0,:], Profile[3,:])

ax = plt.subplot2grid((2,2),(1,0),colspan=1, rowspan=1)
ax.set_title('Tb', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Tb (s)', fontsize = 18)
ax.plot(Profile[0,:], Profile[4,:])

ax = plt.subplot2grid((2,2),(1,1),colspan=1, rowspan=1)
ax.set_title('Urms', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Urms (m/s)', fontsize = 18)
ax.plot(Profile[0,:], Profile[5,:])

plt.savefig('BBB.png')





fig=plt.figure(2, facecolor='White',figsize=[30,30])
ax = plt.subplot2grid((2,2),(0,0),colspan=1, rowspan=1)
ax.set_title('Z', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y - Distance from origin (m)', fontsize = 18)
ax.scatter(Data[0,:], Data[1,:], c = Data[2,:], linewidth = 0, alpha = 0.5)

ax = plt.subplot2grid((2,2),(0,1),colspan=1, rowspan=1)
ax.set_title('Hs', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y - Distance from origin (m)', fontsize = 18)
ax.scatter(Data[0,:], Data[1,:], c = Data[3,:], linewidth = 0, alpha = 0.5)

ax = plt.subplot2grid((2,2),(1,0),colspan=1, rowspan=1)
ax.set_title('Tb', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y - Distance from origin (m)', fontsize = 18)
ax.scatter(Data[0,:], Data[1,:], c = Data[4,:], linewidth = 0, alpha = 0.5)

ax = plt.subplot2grid((2,2),(1,1),colspan=1, rowspan=1)
ax.set_title('Urms', fontsize = 18)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y - Distance from origin (m)', fontsize = 18)
ax.scatter(Data[0,:], Data[1,:], c = Data[5,:], linewidth = 0, alpha = 0.5)


plt.savefig('AAA.png')"""
