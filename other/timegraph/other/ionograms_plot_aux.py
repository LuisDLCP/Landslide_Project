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
import numpy as np
import glob
import os

## --------------- Defining files's path ----------------------##
# Script's path
current_path = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN/" 
# Folder 1's path
source_path = current_path + "Montecillo/"
source_path1 = source_path + "Jicamarca/"
source_path2 = source_path + "Puerto/"

# Folder 2's path
destination_path = current_path + "Raid/"
destination_path1 = destination_path + "Jicamarca/"
destination_path2 = destination_path + "Puerto/"

## --------------- Defining variables ----------------------##
states = {'Ok':0, 
          'Lost':1, 
          'Zero':2, 
          'Same':3}

zeroFile_size_limit = 5 # (Bytes) All files whose size are lower than this value belongs to Zero file     

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

def plot_ionogram():
    pass

def main():
    clasify_files(source_path=source_path1,destination_path=destination_path1) # Clasifiy and move files
    clasify_files(source_path=source_path2,destination_path=destination_path2) # Clasifiy and move files
    
    #save_states_file(st1) # Save log state file
    m = read_data(path=destination_path1)
    print(len(m))
    
if __name__ == '__main__':
    main()












































































