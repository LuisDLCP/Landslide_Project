#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:02:33 2019
**************************************************
      CONVERT IMAGES DATA_SET TO VIDEO & GIF
************************************************** 
change:
    - scr_dir
    - dst_dir
@author: Luis DLCP
"""
import os 

output_file_name = "RawData_151"
root_dir = os.getcwd()+f"/Results/{output_file_name}/Interferograms/Interferograms_complete/"
#root_dir = os.getcwd()+f"/Results/{output_file_name}/SAR_Images/"
scr_dir = root_dir+"Images/" # Source Path
#scr_dir = os.getcwd() + "/Results/Interferograms_BP/Interferograms_complete/Images/" # Source Path
dst_dir = root_dir+"Videos/" #"/Results/Imaging_BP/Videos/" # Destination Path
#dst_dir = os.getcwd() + "/Results/Interferograms_BP/Interferograms_complete/Videos/" # Destination Path

def im2vid(vd_src,vd_dst,vd_name):
    #os.system("ffmpeg -framerate 20 -start_number 1 -i "+vd_src+"ImageBP_dset_%d.png -c:v libx264 -profile:v\
    #          high -crf 20 -pix_fmt yuv420p "+vd_dst+vd_name)
    #os.system("ffmpeg -framerate 20 -start_number 1 -i "+vd_src+"ImageBP_dset_%d.png -vcodec libx264 -y\
    #           -an " + vd_dst + vd_name)
    os.system("ffmpeg -framerate 20 -start_number 0 -i "+vd_src+"Itf_BP_dset%d.png -c:v libx264 -profile:v\
              high -crf 20 -vframes 50 -pix_fmt yuv420p "+vd_dst+vd_name)
    
    return "Ok"

def vid2gif(gif_src,gif_dst,gif_name):
    #os.system("ffmpeg -framerate 40 -start_number 10 -i "+gif_src+"ImageBP_dset_%d.png -vframes 40 -vf scale=600:-1 "+\
    #          gif_dst+gif_name)
    os.system("ffmpeg -framerate 40 -start_number 0 -i "+gif_src+"Itf_BP_dset%d.png -vframes 50 -vf scale=600:-1 "+\
              gif_dst+gif_name)
    
    return "Ok"
    
def main():
    name = "Interferograms_" + output_file_name
    #name = "Itf_dset_10-1170"
    video_name = name + ".mp4"#"Set_interferograms1.mp4" # Video name in .mp4 extension
    gif_name = name + ".gif" # Gif name in .gif extension
    
    im2vid(scr_dir,dst_dir,video_name)
    vid2gif(scr_dir,dst_dir,gif_name)
    
    return "Ok"

if __name__ == '__main__':
    main()

