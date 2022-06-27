#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:26:32 2020
---------------------------------------
    PROGRAMA DE CLASIFICACION DE
    ARCHIVOS Y CAMBIO DE DIRECTORIO
---------------------------------------

@author: soporte
"""
import numpy as np
import glob
import os
import shutil

current_path = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN/" 
source_path = current_path + "test/"
destination_path = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN/test2/"

m = glob.glob(source_path + '*.txt') # List all the source folder's files 
n = glob.glob(destination_path + '*.txt')

state = np.zeros(0) # Array of states: v = {0,1,2,3}
v = 0 # Initial state: 0

file_size_limit = 5 # Bytes 

# Save states in a txt
state_file = open(current_path + "state_file.txt","a+")


if len(m) > 0:
    for i in range(len(m)):
        v = 0 # Ok data
        
        file_name = m[i][len(source_path):] # Get the file's name

        file_size = os.path.getsize(source_path+file_name) # Get the file size in bytes
        
        file_before_name = max(n,key=os.path.getctime)
        file_before_size = os.path.getsize(file_before_name)
        
        if file_size == file_before_size:
            v = 3 # Same data 

        if file_size < file_size_limit:
            v = 2 # Zero data

        os.rename(source_path+file_name,destination_path+file_name) # Move files
        #shutil.move(source_path+file_name,destination_path+file_name)
        #os.replace(source_path+file_name,destination_path+file_name)
        state_file.write(str(v))
        
    print("All the files were transfered!")
else:
    v = 1 # Lost data
    state_file.write(str(v))
    
    print("There are no files yet!")

state_file.close()


# Clasify data
