#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:02:33 2019
**************************************************
      CONVERT IMAGES DATA_SET TO VIDEO & GIF
**************************************************    
@author: Luis DLCP
"""
import os 

scr_dir = os.getcwd()+"/Results/Interferograms_BP/Interferograms_complete/Images/" #"/Results/Imaging_BP/Images/" # Source Path
dst_dir = os.getcwd()+"/Results/Interferograms_BP/Interferograms_complete/Videos/" #"/Results/Imaging_BP/Videos/" # Destination Path

def im2vid(vd_src,vd_dst,vd_name):
    os.system("ffmpeg -framerate 20 -start_number 0 -i "+vd_src+"Itf_BP_dset%d.png -c:v libx264 -profile:v\
              high -crf 20 -pix_fmt yuv420p "+vd_dst+vd_name)
    return "Ok"

def vid2gif(gif_src,gif_dst,gif_name):
    os.system("ffmpeg -framerate 40 -start_number 0 -i "+gif_src+"Itf_BP_dset%d.png -r 15 -vf scale=600:-1 "+\
              gif_dst+gif_name)
    return "Ok"
    
def main():
    video_name = "Set_interferograms1.mp4" # Video name in .mp4 extension
    gif_name = "Set_interferograms1.gif" # Gif name in .gif extension
    im2vid(scr_dir,dst_dir,video_name)
    vid2gif(scr_dir,dst_dir,gif_name)
    return "Ok"

if __name__ == '__main__':
    main()


