#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 14:42:51 2018

@author: koenig.g
"""

#################################################
# Script to analyze tidal components on         #
# An entire CROCO model for the Ouano lagoon    #
# Using utide and compare them with SPM outputs #
#################################################

#**********Packages import********************#
import matplotlib.pyplot as plt
import numpy as np
import utide
import netCDF4

#**********Data import************************#
nc_data=netCDF4.Dataset('his_psource_Nowave.nc','r')
zeta=nc_data.variables['zeta'][:,:,:].data
land=nc_data.variables['mask_rho'][:,:].data
time=nc_data.variables['time_step'][:,0].data
lat=nc_data.variables['lat_rho'][:,:].data
nc_data.close()

#*********Creating data matrices for storing#
#*********Tidal coefficients****************#
Amp_M2_CROCO=np.zeros(land.shape)
Phase_M2_CROCO=np.zeros(land.shape) # I should think about using dictionnaries
# For flexibility
#***Analyzing tidal components**************#
# First I convert time in days
time=time/(1440) # I need to know the time step

for i in range(Amp_M2.shape[0]):
    for j in range(Amp_M2.shape[1]):
        # To reduce computation time I do not compute on land points
        if land[i,j] :
            # Now I solve for the tidal coefficients
            Coeff=utide.solve(time,zeta[:,i,j],lat=lat[i,j],nodal=True,
                          trend=True,method='robust',conf_int='linear',
                          Rayleigh_min=.95)
            # And I extract M2 coefficients for example
            Amp_M2_CROCO[i,j]=Coeff['A'][0]
            Phase_M2_CROCO[i,j]=Coeff['g'][0]
            

#*********Now we have to read SPM data***#
# Importing data
nc_data=netCDF4.Dataset('Serpent_de_Mer_Ouano.nc','r')
Amp_M2_SPM=nc_data.variables['Amplitude'][0,0,:,:].data
Phase_M2_SPM=nc_data.variables['phase'][0,0,:,:].data
nc_data.close()
   
#*******Plotting the results*************#
fig=plt.figure()
# First ax for amplitude
ax=fig.add_subplot(1,1,1)
cmap_1=ax.imshow(Amp_M2_CROCO-Amp_M2_SPM)
cbar_1=plt.colorbar(cmap_1,cmap='Reds')
cbar_1.set_label('Amplitude (m)')
# Then for the phase
ax=fig.add_subplot(1,2,1)
cmap_2=ax.imshow(Phase_M2_CROCO-Phase_M2_SPM*180/np.pi)
cbar_2=plt.colorbar(cmap_2,cmap='Blues')
cbar_2.set_label('Phase (degrees)')