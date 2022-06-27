#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:48:37 2020
--------------------------------------------
            GRAFICAS
--------------------------------------------
@author: soporte
"""
import matplotlib.pyplot as plt
from matplotlib import colors,ticker 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import numpy as np
import random as rd

days = 30*9
n_stations = 4

data = np.random.randint(4,size=(n_stations,days))
#print(data)
font_value = 19

def add_bounds(data_matrix):
  lenx, leny = len(data_matrix[0]), len(data_matrix)
  nData = np.empty((2*leny,lenx))
  nData[:] = font_value #np.nan
  for i in range(leny):
    nData[2*i+1,:] = data_matrix[i,:]

  return nData

x_min, x_max = 1, days # Dias
y_min, y_max = 1, n_stations # Estaciones

# Create array axis values
dias = list(range(1,days+1))
estaciones = ["","Jicamarca","","Puerto","","Tupiza","","Tucuman"]

# Get data
data_plot = add_bounds(data)

# Make a color map
color_map = {"Ok": "#00B1FF",
             "Lost": "#0019FF",
             "Zero": "#53FFA4",
             "Same": "#F8F500",
             "Font": "#000080"}

bounds_map = {"Ok": 0,
             "Lost": 1,
             "Zero": 2,
             "Same": 3,
             "Font": font_value}

colores = [color_map["Ok"],color_map["Lost"],color_map["Zero"],color_map["Same"],color_map["Font"]]
cmap2 = colors.ListedColormap(colores)

bounds = [bounds_map["Ok"]-0.5,bounds_map["Lost"]-0.5,bounds_map["Zero"]-0.5,bounds_map["Same"]-0.5,bounds_map["Same"]+0.5,bounds_map["Font"]+0.5]
norm2 = colors.BoundaryNorm(bounds,cmap2.N)

#cmap = colors.ListedColormap(['g','r','b','y'])
#bounds = [0,1,2,3,4]
#norm = colors.BoundaryNorm(bounds,cmap.N)

# Create legend
class1_box = mpatches.Patch(color=colores[0], label="Correcto")
class2_box = mpatches.Patch(color=colores[1], label="No recibido")
class3_box = mpatches.Patch(color=colores[2], label="Tamaño cero")
class4_box = mpatches.Patch(color=colores[3], label="Mismo tamaño")

# Plotting
fig, ax = plt.subplots(figsize=(20,4))
im = ax.imshow(data_plot, cmap=cmap2, norm=norm2,origin='upper',aspect="auto")#pcolormesh(data_plot) #,extent=[x_min,x_max,y_min,y_max])
ax.set(xlabel="Dias",ylabel="Estaciones",title="Monitoreo de datos de Ionosondas")
ax.grid(axis='x',color='white')
ax.grid(axis='y',which='minor',color='white')

ax.legend(handles=[class1_box,class2_box,class3_box,class4_box],handlelength=0.9,bbox_to_anchor=(1.02,0.7),loc='lower left',borderaxespad=0.)

#ax.minorticks_on()

# Colorbar
"""
diff = norm_bins[1:] - norm_bins[:-1]
tickz = norm_bins[:-1] + diff/2

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.5)

plt.colorbar(im, format = fmt, ticks = tickz, cax=cax)#(im, cmap=cmap, norm=norm, boundaries=bounds,ticks=bounds)
"""
# Add personalized labels 
ax.set_xticks(np.arange(-0.5,len(dias),10))
ax.set_yticks(np.arange(0,len(estaciones),1))
ax.set_yticks(np.arange(-0.5,len(estaciones),1),minor=True)

ax.set_xticklabels(dias)
ax.set_yticklabels(estaciones)

ax.tick_params(axis='both', which='both',length=0) # Erase little tick lines due to grid lines in both axis 

#fig.savefig("lisn1.png")


