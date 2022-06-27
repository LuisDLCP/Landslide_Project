#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 17:55:42 2020
-------------------------------------------------
        PROGRAMA DE CREACION DE ARCHIVOS 
                 ALEATORIOS
-------------------------------------------------
@author: soporte
"""
import datetime
import os

#drc = os.getcwd()
drc = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/LISN"
drc += "/test/"

name = datetime.datetime.now().strftime("%d%B_%H:%M")

file = open(drc + "test_" + name + ".txt","w")
file.write("Hi!")
file.close()
