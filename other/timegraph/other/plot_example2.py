#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:48:37 2020
--------------------------------------------
            GRAFICAS
    (Estructurado en funciones)
--------------------------------------------
@author: soporte
"""
import matplotlib.pyplot as plt
from matplotlib import colors,ticker 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import numpy as np
import random as rd

T = 6 # (hours), this is the period data is sent 
N_datos = 19

bounds_map = {"Ok": 0,
             "Lost": 1,
             "Zero": 2,
             "Same": 3,
             "Font": 19}

stations = ["San Mateo","Jicamarca","Puerto","Tupiza","Tucuman","Cuenca"] # The first name will appear on top of Image 

# Create font rows(zonas intermedias entre filas de datos) 
def add_bounds(data_matrix):
  lenx, leny = len(data_matrix[0]), len(data_matrix)
  nData = np.empty((2*leny,lenx))
  nData[:] = bounds_map['Font'] #np.nan
  for i in range(leny):
    nData[2*i+1,:] = data_matrix[i,:]

  return nData

# Plot ionogram
def plot_ionogram(data_plot): 
    n_stations = len(stations)
    n_datos = len(data_plot[0])
    
    # Make a personalized color map
    color_map = {"Ok": "#00B1FF", # Sky blue
                 "Lost": "#0019FF", # Blue
                 "Zero": "#53FFA4", # Green
                 "Same": "#F8F500", # Yellow
                 "Font": "#000080"} # Dark blue

    colores = [color_map["Ok"],color_map["Lost"],color_map["Zero"],color_map["Same"],color_map["Font"]]
    cmap2 = colors.ListedColormap(colores)

    bounds = [bounds_map["Ok"]-0.5,bounds_map["Lost"]-0.5,bounds_map["Zero"]-0.5,bounds_map["Same"]-0.5,bounds_map["Same"]+0.5,bounds_map["Font"]+0.5]
    norm2 = colors.BoundaryNorm(bounds,cmap2.N)

    # Plotting
    fig, ax = plt.subplots(figsize=(20,4))
    im = ax.imshow(data_plot, cmap=cmap2, norm=norm2,origin='upper',aspect="auto")#pcolormesh(data_plot) #,extent=[x_min,x_max,y_min,y_max])
    ax.set(xlabel="Dias",ylabel="Estaciones",title="Monitoreo de datos de Ionosondas")
    ax.grid(axis='x',color='white')
    ax.grid(axis='y',which='minor',color='white')
        
    # Add personalized labels 
        # Create array axis values
    dias = list(range(1,n_datos+1))
    
    estaciones = []
    for i in range(n_stations):
        estaciones.append("")
        estaciones.append(stations[i])

    ax.set_xticks(np.arange(-0.5,len(dias),int(24/T)))
    ax.set_yticks(np.arange(0,len(estaciones),1))
    ax.set_yticks(np.arange(-0.5,len(estaciones),1),minor=True)

    ax.set_xticklabels(dias)
    ax.set_yticklabels(estaciones)

    ax.tick_params(axis='both', which='both',length=0) # Erase little tick lines due to grid lines in both axis 

    # Create legend
    class1_box = mpatches.Patch(color=colores[0], label="Correcto")
    class2_box = mpatches.Patch(color=colores[1], label="No recibido")
    class3_box = mpatches.Patch(color=colores[2], label="Tamaño cero")
    class4_box = mpatches.Patch(color=colores[3], label="Mismo tamaño")

    ax.legend(handles=[class1_box,class2_box,class3_box,class4_box],handlelength=0.9,bbox_to_anchor=(1.02,0.7),loc='lower left',borderaxespad=0.)
    
    return 'Ok'

def main():
    data = np.random.randint(4,size=(len(stations),N_datos)) # Test data

    data_n = add_bounds(data)
    plot_ionogram(data_n)
    

if __name__ == '__main__':
    main()
    