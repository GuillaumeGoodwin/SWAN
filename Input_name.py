"""This file generates an input file for every bathymetry file"""

import os
import numpy as np


# The spring tidal range [m]
TR_spring = ["1", "3", "7"]
# The elevation budget: relative elevation of the lowest proint of the ideal marsh platform
Z_budget = ["-0.4", "0.6"] # 5
# The scarp height [m]
Scarp_height = ["0", "0.25", "1", "2.5"] # 6

# The tidal flat slope [m/m]
TF_slope = ["5e-05", "0.0002"]# *10**(-4) # 5


for TR in TR_spring:
    for z0 in Z_budget:
        for h in Scarp_height:
            for s in TF_slope:
                src = open("Template.swn", "r")
                dst = open("Input_TR%s_Zo%s_H%s_S%s.swn" % (TR, z0, h, s), "w")

                GEOFILE = "Points_TR%s_Zo%s_H%s_S%s" % (TR, z0, h, s)
                BOTFILE = "Points_TR%s_Zo%s_H%s_S%s.bot" % (TR, z0, h, s)
                OUTPUTFILE = "Sim_TR%s_Zo%s_H%s_S%s.mat" % (TR, z0, h, s)

                counter = 0
                for columns in (raw.strip().split() for raw in src):
                    line = ''
                    for i in range(len(columns)):
                        line = line + ' ' + columns[i]
                    A = line + '\n'

                    if counter == 12:
                        A = " READ UNSTRUCTURED triangle '" + GEOFILE + "'"+ '\n'
                    elif counter == 15:
                        A = " READINP BOTTOM -1. '" + BOTFILE +  "' 1 0 FREE" + '\n'
                    elif counter == 31:
                        A = " BLOCK 'COMPGRID' NOHEAD '" + OUTPUTFILE + "' LAY 3 HS TPS UBOT URMS TMBOT" + '\n'

                    dst.write(A)
                    counter = counter + 1

                dst.close
