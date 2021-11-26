#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:57:25 2021

@author: soporte
"""
import os 

root_path = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/Landslide_Project/Desarrollo_v2/Software/Procesamiento/Algoritmo_procesamiento/Results/Cuenca_Apr-Oct2021_test2/Interferograms/Interferograms_complete/"
files_input_path = root_path + "Images/"
files_output_path = root_path + "Videos/"
if not os.path.exists(files_output_path):
    os.makedirs(files_output_path)
    
cmd = "ffmpeg -framerate 50 -start_number 0 -i "+ files_input_path + "Itf_BP_dset%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p "+ files_output_path + "SARinterferogram.mp4"

print("Getting started ...")
os.system(cmd)
print("Done!")