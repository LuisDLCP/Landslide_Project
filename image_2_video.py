#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:02:33 2019
**************************************************
       CONVERT IMAGES DATA_SET TO VIDEO
**************************************************    
@author: Luis DLCP
"""
import os 

scr_dir = os.getcwd()+"/Results/Imaging_BP/Images/" # Source Path
dst_dir = os.getcwd()+"/Results/Imaging_BP/Videos/" # Destination Path

def im2vid(vd_name):
    os.system("ffmpeg -framerate 5 -start_number 10 -i "+scr_dir+"ImageBP_dset_%d.png -c:v libx264 -profile:v\
              high -crf 20 -pix_fmt yuv420p "+dst_dir+vd_name)
    return "Ok"

def main():
    video_name = "Reflectivity_Images.mp4" # Video name in .mp4 extension
    ans = im2vid(video_name)
    return ans

if __name__ == '__main__':
    main()


