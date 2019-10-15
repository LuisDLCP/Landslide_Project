#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:26:37 2019

@author: soporte
"""

import os
#import sys
#sys.path.insert(1,os.getcwd()+'/Results/Imagenes_finales_reconstruidas_BP/')
#import TEXT as tx

#os.chdir(os.path.dirname(__file__))
#print("La ruta es: "+os.getcwd())
i=0
scr_dir = os.getcwd()+"/Results/Imagenes_finales_reconstruidas_BP/"
"""
for filename in os.listdir(scr_dir):
    print(filename)
    dst = "Imagen_"+str(i+10)+".png"#"Imagen_"+str(i).zfill(3)+".png" # .zfill, it zero padded a number /media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/Landslide_Project/Desarrollo/Software/Procesamiento
    src = scr_dir+filename
    dst = scr_dir+dst
    
    os.rename(src,dst)
    i += 1
"""    
os.system("ffmpeg -framerate 5 -start_number 10 -i "+scr_dir+"Imagen_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p "+scr_dir+"Rec_img4.mp4")
#tx.texto()

    
#%% Using terminal instructions in python

import subprocess as sb

sb.check_call(['ls','-la'])
rs = sb.Popen(['ffmpeg',
               '-framerate','5',
               'i','image_name_%d.png',
               '-c:v','libx264',
               '-profile:v','high',
               '-crf','20'
               ])

#%% probes
j = 0    
while True:
    for i in range(10):
        if i==8:
            print("An interruption ocurred!")
            break
        print("Normal sequency!")
        
    j += 1
    
    if j==3:
        print("The second interruption ocurred!")
        break
    
    print("Continue ....")

print("Finishing ...")

#%% creating a txt file
f=open("myfirst_file.txt","w+")
f.write("Hello word!\n")
f.write("How are you?\n")

f.close()

#%% media

a=np.random.randint(10,size=(5))
a[3]=np.nan

































