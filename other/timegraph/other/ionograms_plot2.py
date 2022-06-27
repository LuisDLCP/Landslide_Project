#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:59:28 2020
---------------------------------------
    PROGRAMA DE CLASIFICACION DE
    ARCHIVOS Y CAMBIO DE DIRECTORIO
---------------------------------------

@author:LuisDLCP
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
import numpy as np
import glob
import os

## --------------- Defining files's path ----------------------##
# Script's path
current_path = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN/" 
# Folder 1's path
source_path = current_path + "Montecillo/"
# Folder 2's path
destination_path = current_path + "Raid/"

## --------------- Defining variables ----------------------##
states = {'Ok':0, 
          'Lost':1, 
          'Zero':2, 
          'Same':3,
          'Font':19}

zeroFile_size_limit = 5 # (Bytes) All files whose size are lower than this value belongs to Zero file     

T = 6 # (hours), this is the period data is sent 
stations = ["Jicamarca","Puerto","Tupiza","Tucuman"] # The first name will appear on top of Image 


## --------------------------------------------------------- ##

def clasify_files():
    source_files = glob.glob(source_path + '*.txt') # List all the source folder's files 
    destination_files = glob.glob(destination_path + '*.txt') # List all destination folder's files

    state_list = np.zeros(0) # Array of current states
    state_value = 0 # Value of a state 

    if len(source_files) > 0:
        for i in range(len(source_files)):
            state_value = states['Ok'] # Ok data
        
            file_name = source_files[i][len(source_path):] # Get the file's name
            file_size = os.path.getsize(source_path+file_name) # Get the file size in bytes
        
            file_before_name = max(destination_files,key=os.path.getctime) # Get the previous file's name
            file_before_size = os.path.getsize(file_before_name) # Get the previous file's size in bytes
        
            if file_size == file_before_size:
                state_value = states['Same'] # Same data 

            if file_size < zeroFile_size_limit:
                state_value = states['Zero'] # Zero data

            state_list = np.append(state_list,state_value)
            
            move_file(source_path+file_name,destination_path+file_name) # Move file 
            
        print("All the files were transfered!")
    
    else:
        state_value = states['Lost'] # Lost data
        state_list = np.append(state_list,state_value)
        print("There are no files yet!")

    return state_list

def move_file(source_path,dest_path):
    os.rename(source_path, dest_path)
    return 'Ok'
    
def save_states_file(state_list):
    state_file = open(current_path + "states_log.txt","a+")
    
    for element in state_list:
        state_file.write(str(int(element)))
    
    state_file.close()
    
    return 'Ok'

def read_data():
    datan = open(current_path + "states_log.txt","r")
    data_list = datan.read().split(',')
    data_list = data_list[:-1]
    datan.close()
    
    return data_list

# Create font rows for plotting (zonas intermedias entre filas de datos) 
def add_bounds(data_matrix):
  lenx, leny = len(data_matrix[0]), len(data_matrix)
  nData = np.empty((2*leny,lenx))
  nData[:] = states['Font'] #np.nan
  for i in range(leny):
    nData[2*i+1,:] = data_matrix[i,:]

  return nData

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

    bounds = [states["Ok"]-0.5,states["Lost"]-0.5,states["Zero"]-0.5,states["Same"]-0.5,states["Same"]+0.5,states["Font"]+0.5]
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
    st = clasify_files() # Clasifiy and move files
    save_states_file(st) # Save log state file
    
    datos = read_data()
   
    
if __name__ == '__main__':
    main()












































































