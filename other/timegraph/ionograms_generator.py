#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 18:52:54 2020
-------------------------------------------------
        PROGRAMA DE CREACION DE ARCHIVOS 
                 ALEATORIOS
-------------------------------------------------
@author: LuisDLCP
"""
import datetime

drc = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN/Montecillo/"
drc1 = drc + "/Jicamarca/"
drc2 = drc + "/Puerto/"
drc3 = drc + "/Tupiza/"
drc4 = drc + "/Tucuman/"

name = datetime.datetime.now().strftime("%d%B_%H:%M")

file = open(drc1 + "ionogram_" + name + ".txt","w")
file.write("Hi!")
file.close()

file2 = open(drc2 + "ionogram_" + name + ".txt","w")
file2.write("Welcome!")
file2.close()

file3 = open(drc3 + "ionogram_" + name + ".txt","w")
file3.write("How are you doing?")
file3.close()

file4 = open(drc4 + "ionogram_" + name + ".txt","w")
file4.write("Yes!")
file4.close()
