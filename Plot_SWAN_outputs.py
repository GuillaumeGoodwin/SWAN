"""
This Python script plots the results of Swan runs
"""

#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import math as mt
import scipy as sp
from datetime import datetime
from matplotlib import cm
from pylab import *
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools as itt
from osgeo import gdal, osr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import cPickle

#----------------------------
# These are a bunch of functions


#-----------------------
# Calculates the bottom shear stress
def bottom_shear_stress(Ubot,Tbot, Urms, D50):
    Tau = np.zeros(len(Ubot), dtype  = np.float)

    RHOw = 1027
    ni = 0.00000101

    for i in range(len(Ubot)):
        if Ubot[i] != 0:
            Xi = Tbot[i]*Urms[i]/2*np.pi
            ks = 2.5*D50
            fw_rough = 1.39*(Xi/(ks/30))**(-0.52)

            Rewaves=Ubot[i]**2*Tbot[i]/(2*np.pi*ni)
            if (Rewaves<=500000):
                fw_smooth=2*Rewaves**(-0.5)
            else:
                fw_smooth=0.0521*Rewaves**(-0.187)

            Tau[i] = 0.5 * RHOw * max(fw_rough, fw_smooth) * Ubot[i]**2

    return Tau


#-----------------------
# Calculates the wavenumber
def Bisection_k (depth, Tp):
    k = np.zeros(len(Tp), dtype  = np.float)

    g = 9.806
    omega = 2*np.pi/Tp

    for i in range(len(k)):
        kmin = 0.
        kmax = 1000.
        deltak_min = 0.001

        iteration=0
        knext=0.

        if depth[i] > 0:
            while abs(kmin-kmax)>deltak_min and iteration<1000:
                iteration=iteration+1
                Fk_0=np.sqrt(g*kmin*np.tanh(kmin*depth[i]))-omega[i] #Lamb's equation expressed as f(k)
                Fk_1=np.sqrt(g*kmax*np.tanh(kmax*depth[i]))-omega[i] #Lamb's equation expressed as f(k)
                if Fk_0*Fk_1<0:
                    knext=(kmin+kmax)/2
                    Fk_0=np.sqrt(g*knext*np.tanh(knext*depth[i]))-omega[i]
                    if Fk_0*Fk_1<0:
                        kmin=knext
                    else:
                        kmax=knext

                k[i]=knext

    return k






#----------------------------
# Calculates the wave power density
def Wave_power(depth, Hs, Tp):
    RHOw=1027
    g=9.806
    omega = 2*np.pi/Tp

    k = Bisection_k(depth, Tp)


    Cg = g/(2*omega) * (k/np.cosh(k*depth)**2 + np.tanh(k*depth))


    Pw =  Cg * Hs**2 * RHOw*g/8 /1000

    return Pw




#-----------------------
# See Priestas2015
def Wave_thrust(Z0, TR, H, Hs):
    Sea_level = -Z0*TR/2
    Thrust = np.zeros(len(Hs), dtype = np.float)

    if H > 0:
        if Sea_level > 0.4:
            Thrust = 4.4 * Hs
        elif Sea_level > 0 and Sea_level <= 0.4:
            Thrust = 8.8*Hs - 11*(Sea_level)
        elif Sea_level > -H and Sea_level <= 0:
            Thrust = 8.8*Hs

    Thrust = Thrust * 1

    return Thrust



#----------------------------
# Here begins the scripture

TR_spring = ["1", "3", "7"]
Z_budget = ["-0.4", "0.6"]
Scarp_height = ["0", "0.25", "1", "2.5"]
TF_slope = ["5e-05", "0.0002"]


#TR_spring = ["1"]
#Z_budget = ["-0.4", "0.6"]
#Scarp_height = ["2.5"]
#TF_slope = ["0.0002"]

TR_spring_num = np.asarray(TR_spring).astype(np.float)
Z_budget_num = np.asarray(Z_budget).astype(np.float)
Scarp_height_num = np.asarray(Scarp_height).astype(np.float)
TF_slope_num = np.asarray(TF_slope).astype(np.float)



# The median grain size [m]
D50 = np.array([20, 50, 100, 200, 500, 1000, 2000])*10**(-6)




#-----------------------------------
for TR in range(len(TR_spring_num)):
    for z0 in range(len(Z_budget_num)):
        for h in range(len(Scarp_height_num)):
            for s in range(len(TF_slope_num)):
                INPUT = "Swan_15ms/Out/Profile_TR%s_Zo%s_H%s_S%s.txt" % (TR_spring[TR], Z_budget[z0], Scarp_height[h], TF_slope[s])
                Data = np.loadtxt(INPUT, comments="#", delimiter=" ", unpack=False)
                # 0 -> X ; 1 -> Y ; 2 -> Z ; 3 -> Hs ; 4 -> Tbot ; 5 -> Urms ; 6 -> Tp ; 7 -> Ubot

                Data[isnan(Data)] = 0

                Sea_level = - Z_budget_num[z0]*TR_spring_num[TR]/2
                Ground = Data[2] + Sea_level
                Waves = Data[3] + Sea_level; Waves[isnan(Waves)] = Sea_level


                Tau = np.zeros((len(D50), len(Data[0])), dtype = np.float)
                for d in range (len(D50)):
                    Tau[d] = bottom_shear_stress(Data[7],Data[4], Data[5], D50[d])

                Pw = Wave_power (Sea_level - Ground, Data[3], Data[4]); Pw[isnan(Pw)] = 0
                Thrust = Wave_thrust (Z_budget_num[z0], TR_spring_num[TR], Scarp_height_num[h], Data[3]); Thrust[Thrust<0] = 0; Thrust[isnan(Thrust)] = 0


                #cPickle.dump(k,open("0_WWTM_k.pkl", "wb"))








                #######################################################
                # Plot Data
                fig=plt.figure(1, facecolor='White',figsize=[30,15])


                #------------------------------------------------------
                ax1 = plt.subplot2grid((4,4),(2,0),colspan=3, rowspan=2)
                ax1.set_ylabel('Hs (m)', fontsize = 18)
                ax1.set_xlabel('Distance (m)', fontsize = 18)
                ax1.grid(True)
                if max(Ground)>Sea_level:
                    ax1.set_ylim(ymin = -3.7, ymax = max(Ground)+max(Data[3]) +0.2)
                else:
                    ax1.set_ylim(ymin = -3.7, ymax = max(Sea_level+Waves)+0.2)
                ax1.set_xlim(xmin = 0, xmax = 20000)


                Scatt = ax1.scatter(Data[0],  Waves, c = Data[6], linewidth = 0)
                divider = make_axes_locatable(ax1)
                cax = divider.append_axes("bottom", size="5%", pad=0.60)
                cbar = plt.colorbar(Scatt, orientation = 'horizontal', cax=cax)
                cbar.set_label('Wave period (s)', fontsize = 20)

                ax1.plot(Data[0,:], Ground,'k')
                ax1.fill_between(Data[0],Ground, Sea_level, facecolor ='blue' , alpha = 0.5, linewidth = 0)
                ax1.fill_between(Data[0],Sea_level, Waves, facecolor ='blue' , alpha = 0.2)
                ax1.fill_between(Data[0],-10, Ground, facecolor ='grey' , alpha = 1)






                #---------------------------------------------------------
                ax2 = plt.subplot2grid((4,4),(1,0),colspan=3, rowspan=1)
                ax2.set_ylabel('Tau (Pa)', fontsize = 18)
                ax2.set_xticklabels([])
                ax2.grid(True)
                ax2.set_ylim(ymin = -0.15, ymax = np.amax(Tau)+0.2)
                ax2.set_xlim(xmin = 0, xmax = 20000)

                for d in range (len(D50)):
                    ax2.plot(Data[0,:], Tau[d,:], color=plt.cm.jet(d*40))

                #---------------------------------------------------------
                ax3 = plt.subplot2grid((4,4),(0,0),colspan=3, rowspan=1)
                ax3.set_ylabel('Power density (kW/m)', fontsize = 18)
                ax3.set_xticklabels([])
                ax3.grid(True)
                ax3.set_ylim(ymin = -0.15, ymax = max(max(Pw),max(Thrust))+0.2)
                ax3.set_xlim(xmin = 0, xmax = 20000)

                #ax3.plot(Data[0,:], Pw, 'k-')
                ax3.plot(Data[0,:], Data[4,:], 'k-')
                ax3.plot(Data[0,:], Data[6,:], 'r-')







                ###############----------########################
                #Zoom in
                #------------------------------------------------------
                ax11 = plt.subplot2grid((4,4),(2,3),colspan=1, rowspan=2)
                ax11.set_yticklabels([])
                ax11.set_xlabel('Distance (m)', fontsize = 18)
                ax11.grid(True)
                if max(Ground)>Sea_level:
                    ax11.set_ylim(ymin = -3.7, ymax = max(Ground)+max(Data[3]) +0.2)
                else:
                    ax11.set_ylim(ymin = -3.7, ymax = max(Sea_level+Waves)+0.2)
                ax11.set_xlim(xmin = 390, xmax = 410)

                Scatt = ax11.scatter(Data[0],  Waves, c = Data[6], linewidth = 0)
                divider = make_axes_locatable(ax11)
                cax = divider.append_axes("bottom", size="5%", pad=0.60)
                cbar = plt.colorbar(Scatt, orientation = 'horizontal', cax=cax)
                cbar.set_label('Wave period (s)', fontsize = 20)

                ax11.plot(Data[0,:], Ground,'k')
                ax11.fill_between(Data[0],Ground, Sea_level, facecolor ='blue' , alpha = 0.5, linewidth = 0)
                ax11.fill_between(Data[0],Sea_level, Waves, facecolor ='blue' , alpha = 0.2)
                ax11.fill_between(Data[0],-10, Ground, facecolor ='grey' , alpha = 1)



                #---------------------------------------------------------
                ax21 = plt.subplot2grid((4,4),(1,3),colspan=1, rowspan=1)
                ax21.set_yticklabels([])
                ax21.set_xticklabels([])
                ax21.grid(True)
                ax21.set_ylim(ymin = -0.2, ymax = np.amax(Tau)+0.2)
                ax21.set_xlim(xmin = 201, xmax = 600)

                for d in range (len(D50)):
                    ax21.plot(Data[0,:], Tau[d,:], color=plt.cm.jet(d*40))

                #---------------------------------------------------------
                ax31 = plt.subplot2grid((4,4),(0,3),colspan=1, rowspan=1)
                ax31.set_yticklabels([])
                ax31.set_xticklabels([])
                ax31.grid(True)
                ax31.set_ylim(ymin = -0.2, ymax = max(max(Pw),max(Thrust))+0.2)
                ax31.set_xlim(xmin = 201, xmax = 600)
                #ax31.plot(Data[0,:], Pw, 'k-')
                #ax31.plot(Data[0,:], Thrust, 'k--')
                ax31.plot(Data[0,:], Data[4,:], 'k-')
                ax31.plot(Data[0,:], Data[6,:], 'r-')


                #-----------------------------
                fig.subplots_adjust(hspace=0)




                plt.savefig('Swan_15ms/Fig/Profile_TR%s_Zo%s_H%s_S%s.png' % (TR_spring[TR], Z_budget[z0], Scarp_height[h], TF_slope[s]))




"""#----------------------------------------------------------------
ax = plt.subplot2grid((4,4),(0,3),colspan=1, rowspan=2)
ax.set_title('Zoom on the marshes', fontsize = 18)
ax.set_xticklabels([])
ax.grid(True)
ax.set_ylim(ymin = -1, ymax = 1)
ax.set_xlim(xmin = 0, xmax = 1000)
ax.fill_between(Data[0],-10, Data[2], facecolor ='grey' , alpha = 0.5)
ax.fill_between(Data[0],Data[2], Sea_level, facecolor ='blue' , alpha = 0.5, linewidth = 0)
ax.fill_between(Data[0],Sea_level, Data[3], facecolor ='blue' , alpha = 0.2)

ax.plot(Data[0,:], Data[2,:],'k')
ax.plot(Data[0,:], Data[3,:], '--')



ax = plt.subplot2grid((4,4),(2,3),colspan=1, rowspan=1)
ax.set_xticklabels([])
ax.grid(True)
ax.set_xlim(xmin = 0, xmax = 1000)

ax.plot(Data[0,:], Data[6,:])



ax = plt.subplot2grid((4,4),(3,3),colspan=1, rowspan=1)
ax.set_xlabel('X - Distance from origin (m)', fontsize = 18)
ax.grid(True)
ax.set_xlim(xmin = 0, xmax = 1000)

#for d in range (len(D50)):
    #ax.plot(Data[0,:], Tau[d,:], color=plt.cm.jet(d*40))"""


#plt.savefig('Swan_15ms/Fig/Profile_TR%s_Zo%s_H%s_S%s.png' % (TR_spring[TR], Z_budget[z0], h, s))
