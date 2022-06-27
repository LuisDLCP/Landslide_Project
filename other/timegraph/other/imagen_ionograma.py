#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 16:59:28 2020
---------------------------------------
    PROGRAMA DE CLASIFICACION DE
    ARCHIVOS, CAMBIO DE DIRECTORIO Y
                GRAFICAS
---------------------------------------

@author:LuisDLCP
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors
import numpy as np
import datetime
import glob
import os

## --------------- Defining files's paths ----------------------##
# Script's path
current_path = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN/" 
# Folder 1's path
source_path = current_path + "Montecillo/"
source_path1 = source_path + "Jicamarca/"
source_path2 = source_path + "Puerto/"
source_path3 = source_path + "Tucuman/"
source_path4 = source_path + "Tupiza/"
src_path = {0:source_path1, 1:source_path2, 2:source_path3, 3:source_path4}

# Folder 2's path
destination_path = current_path + "Raid/"
destination_path1 = destination_path + "Jicamarca/"
destination_path2 = destination_path + "Puerto/"
destination_path3 = destination_path + "Tucuman/"
destination_path4 = destination_path + "Tupiza/"
dest_path = {0:destination_path1, 1:destination_path2, 2:destination_path3, 3:destination_path4}

## --------------- Defining variables ----------------------##
states = {'Ok':0, 
          'Lost':1, 
          'Zero':2, 
          'Same':3,
          'Font':19}

zeroFile_size_limit = 5 # (Bytes) All files whose size are lower than this value belongs to Zero file     

T = 4 # (minutes), this is the period data is sent 
stations = ["Jicamarca","Puerto","Tucuman","Tupiza"] # The first name will appear on top of Image 

## --------------------------------------------------------- ##

def clasify_files(source_path = None, destination_path = None):
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
    
    # Save log states file 
    save_states_file(data=state_list, path=destination_path)
    
    return 'Ok'

def move_file(source_path,dest_path):
    os.rename(source_path, dest_path)
    return 'Ok'
    
def save_states_file(data=None, path=None):
    state_file = open(path + "states_log.txt","a+")
    
    for element in data:
        state_file.write(str(int(element))+',')
    
    state_file.close()
    
    return 'Ok'

def read_data(path=None):
    datan = open(path + "states_log.txt","r")
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
    n_stations = len(data_plot)
    n_datos = len(data_plot[0])
    
    nData = add_bounds(data_plot)
    
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
    
    # Get current date & time
    #current_time = datetime.datetime.now().strftime("%d %b %Y %H:%M") # Local time
    current_time = datetime.datetime.utcnow().strftime("%d %b %Y %H:%M") # UTC time
    
    # Plotting
    fig, ax = plt.subplots(figsize=(20,4))
    im = ax.imshow(nData, cmap=cmap2, norm=norm2,origin='upper',aspect="auto")#pcolormesh(data_plot) #,extent=[x_min,x_max,y_min,y_max])
    ax.set(xlabel="Horas",ylabel="Estaciones")
    ax.set_title("Perú",fontsize=20)
    ax.set_title(current_time+' (UTC)',loc='left',fontsize=10)
    ax.grid(axis='x',color='white')
    ax.grid(axis='y',which='minor',color='white')
        
    # Add personalized labels 
        # Create array axis values
    dias = list(range(1,n_datos+1))
    
    estaciones = []
    for i in range(n_stations):
        estaciones.append("")
        estaciones.append(stations[i])

    ax.set_xticks(np.arange(-0.5,len(dias),int(60/T)))#int(24/T)))
    ax.set_yticks(np.arange(0,len(estaciones),1))
    ax.set_yticks(np.arange(-0.5,len(estaciones),1),minor=True)

    ax.set_xticklabels(dias)
    ax.set_yticklabels(estaciones)
    
    # Set color to text ytickLabel
    last_state = nData.T[-1]
    def color_state(x):
        return{
                0:'Black',
                1:color_map["Lost"],
                2:color_map["Zero"],
                3:color_map["Same"]                    
                }.get(x,"Red")
    colorText = []
    for i in last_state:
        colorText.append(color_state(i))
    for color,tick in zip(colorText,ax.yaxis.get_major_ticks()):
        tick.label1.set_color(color)
    
    ax.tick_params(axis='both', which='both',length=0) # Erase little tick lines due to grid lines in both axis 

    # Create legend
    class1_box = mpatches.Patch(color=colores[0], label="Correcto")
    class2_box = mpatches.Patch(color=colores[1], label="No recibido")
    class3_box = mpatches.Patch(color=colores[2], label="Tamaño cero")
    class4_box = mpatches.Patch(color=colores[3], label="Mismo tamaño")

    ax.legend(handles=[class1_box,class2_box,class3_box,class4_box],handlelength=0.9,bbox_to_anchor=(1.02,0.7),loc='lower left',borderaxespad=0.)
    
    # Save figure
    current_time2 = datetime.datetime.now().strftime("%d%B_%H:%M")
    plt.savefig(current_path+"Figures/Ionogram_"+current_time2+".png")
    
    return 'Ok'

def main():
    #%
    ex = open(current_path+"hola.txt","a+")
    ex.write("holaaaa")
    ex.close()
    #%    
    # Clasify and move files
    for i in range(len(stations)):
        clasify_files(source_path=src_path[i],destination_path=dest_path[i])
    
    # Read data
    data_set = {}
    for i in range(len(stations)):
        data_set[i] = read_data(path=dest_path[i]) 
        
    raw_data = np.zeros((len(stations),len(data_set[0]))) # Fix len of each station data 
    for i in range(len(raw_data)):
        raw_data[i,:] = data_set[i]
    
    # Plot
    plot_ionogram(raw_data)
  
ex = open(current_path+"hola.txt","a+")
current_time3 = datetime.datetime.now().strftime("%d%B_%H:%M")
ex.write(current_time3 + "start\n")
ex.close()    

main()


"""
if __name__ == '__main__':
    main()
    
"""












































































