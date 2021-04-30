#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:09:44 2019
---------------------------------------------------
                INTERFEROMETRIA
                (Main Program)
---------------------------------------------------
* Se obtendran los mapas de desplazamientos 
* Se obtendran las curvas de desplazamientos 
* Se usara el algoritmo FDBP para la formacion de imagenes

@author: LuisDLCP
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplt
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from glob import glob 
#import RMA_real_main as RMA
import BP_real_main as BP
import sarPrm as sp
import time
import timeit
import datetime
import random
import os

# specify the input files directory 
dir_input_files = "/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/Landslide_Project/Desarrollo_v2/Software/Procesamiento/Data_set/CUENCA_Apr2021_22-04-21_14:47:51/"
# get the current path
os.chdir(os.path.dirname(__file__)) 

show = True # esta variable muestra ciertas figuras 
n_im = len(glob(dir_input_files+"*.hdf5")) #2 #1170 #1805 #2000 #5 #1170 #400 #1170 #4656 # Numero de imagenes a considerar
i_o = 1 #10 #100 #10 # Numero de imagen inicial(10)

directory0 = os.getcwd()+"/Results/RawData_"+str(n_im) # root directory to save all the results
if not os.path.exists(directory0):
    os.makedirs(directory0)

directory1 = directory0+"/Output_Imaging" # directory to save data, Dates & Parameters
if not os.path.exists(directory1):
    os.makedirs(directory1)
            
directory2 = directory1+"/data" # directory to save data
if not os.path.exists(directory2):
    os.makedirs(directory2)

def get_images(algorithm=None):
    """ Esta funcion se encarga de graficar las imagenes SAR usando un algoritmo dado por 
        "algorithm", luego guarda la siguiente informacion en la siguiente ruta:
        -> Imagenes: 
                dir: /Results/Output_.../Im_...
        -> Parametros: Parametros como L,BW,fc que se usaron para obtener el Raw Data
                dir: Parameters_...npy
        -> Fechas: La fecha que se tomo el Raw Data asociado a cada imagen
                dir: Dates_...npy
    """
    if algorithm == "RMA":
        #Ims = {}
        dates = []
        #Ims = np.load("Set_images_RMA.npy").item()
        #dates = np.load("Dates_RMA.npy")
        for i in range(n_im):
            i += i_o # Empieza en la posicion 10
            data = RMA.main("dset_"+str(i)+".hdf5")
            #Ims[10+i] = data['Sf_n']
            dates.append(data['date'])
            np.save(os.getcwd()+"/Results/Output_RMA/Im_"+str(i)+".npy",data['Sf_n'])
        #np.save("Set_images_RMA",Ims) # Para guardar el set de imagenes
        np.save("Parameters_RMA",data)
        np.save("Dates_RMA",np.array(dates))

    elif algorithm == "BP":
        #Ims = {}
        dates = []            
        #Ims = np.load("Set_images_BP.npy").item()
        #dates = np.load("Dates_BP.npy")
        for i in range(n_im): #(4991):
            i += i_o # Empieza en la posicion 10
            data = BP.main("dset_"+str(i)+".hdf5",i-i_o, dir_input= dir_input_files, directory=directory0) 
            #Ims[i] = data['Im']
            dates.append(data['date'])
            np.save(directory2 + "/Im_"+str(i)+".npy",data['Im']) # Imagenes de todo el dataset
        np.save(directory1 + "/Parameters",data) # Parametros geometricos como dimensiones y grilla de la imagen
        np.save(directory1 + "/Dates",np.array(dates)) # Fechas de las iamgenes tomadas de todo el dset
        
        # Saving used raw data dsets
        g = open(directory0 + "/dir_rawdata.txt","a+")
        g.write("\n\nDATA SETs:\n")
        g.write("\ndset_"+str(i_o)+" - dset_"+str(i_o+n_im-1))
        g.close()

    return 'Ok'

# This function calibrates each image after reconstruction
def calibration(msk,Im):
    """    
    # 1) Create a mask with the highest value of coherence in 1 
    M = (coh>=coh.max()) # bool type
    M = M*1 # int type
        # ...This part helps to choose randonly only a value if there's many maximums  
    M_aux = np.linspace(1,100,len(M)*len(M[0]))
    random.shuffle(M_aux)
    M_aux = M_aux.reshape(len(M),len(M[0]))
    M = M*M_aux
    M = (M==M.max())
    M = M*1
    """    
    # 2) Multiply the image by the mask
    M = Im*msk
        # 3) Get the phase factor 
    pf = np.angle(np.sum(M))
    #print("Phase without compensation: ",np.angle(Im[0,0]))
        # 4) Finally calculate the difference of phases 
    M2 = Im*np.exp(-1j*pf)
    #print("Phase with compensation: ",np.angle(M3[0,0]))#np.angle(np.sum(M3*M)))

    return M2

def make_interferometry(data,algorithm=None):
    """ Esta funcion se encarga de obtener los  mapas de desplazamientos, o los 
        interferogramas, esto lo realiza usando 2 imagenes consecutivas para obtener
        una imagen de mapa de desplazamientos.
        Para ello pasa por las siguientes etapas:
        1) Obtencion de mapa de coherencia: Esto se realiza para todo el set de imagenes,
                                            con ello se obtiene la mascara.
        2) Obtencion de los interferogramas usando la mascara
        3) Obtencion de los mapas de desplazamientos
        4) Se divide el mapa en zonas 
        5) Se obtiene el promedio de cada zona por interferograma 
        6) Finalmente se guarda los promedios por cada zona y por cada interferograma 
           en una matriz.
    """
    # Obteniendo el mapa de coherencia
    Im = np.load(directory2+"/Im_"+str(i_o)+".npy") # Se carga la imagen 10
    f1 = np.zeros((len(Im),len(Im.T)),dtype=complex)
    f2 = np.zeros((len(Im),len(Im.T)),dtype=complex)
    coh = np.zeros((len(Im),len(Im.T)),dtype=complex)

    for i in range(n_im):
        i += i_o # Valor inicial
        Im1 = np.load(directory2+"/Im_"+str(i)+".npy")
        if i<n_im+i_o-1:
            Im2 = np.load(directory2+"/Im_"+str(i+1)+".npy")
            f1 += Im1*Im2.conjugate()
        f2 += abs(Im1)**2 #np.sqrt((abs(Im1)**2)*(abs(Im2)**2))

    coh = (n_im/(n_im-1))*abs(f1/f2)
    
    # Graficando la coherencia
    if show:
        plt.close('all')
        cmap ="plasma"
        if algorithm == "BP":
            title_name ='Mapa de coherencia(BP)\n[dset'+str(i_o)+'-'+str(i_o+n_im-1)+']'
            direction ='Coherencia_BP_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'.png'

        elif algorithm == "RMA":
            title_name ='Mapa de coherencia(RMA)\n[dset'+str(i_o)+'-'+str(i_o+n_im-1)+']'
            direction ='Coherencia_RMA_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'].png'
        #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
        #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20

        fig, ax = plt.subplots()
        im=ax.imshow(coh,cmap,origin='lower',aspect='equal', extent=[data['x_min'],data['x_max'],data['y_min'],data['y_max']]) #, vmin=vmin, vmax=vmax)
        ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
        ax.grid()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
        plt.colorbar(im,cax=cax,label='',extend='both')
        
        directory_coh = directory0 + "/Interferograms/Coherence_maps/"
        if not os.path.exists(directory_coh):
            os.makedirs(directory_coh)
        
        fig.savefig(directory_coh+direction, orientation='landscape')

    # Sacando una mascara de [0.7,1] para la coherencia
    mask = coh<=0.7
    coh2 = np.copy(coh)
    coh[mask] = np.nan

    # Graficando con la mascara
    if show:
        cmap ="plasma"
        if algorithm == "BP":
            title_name ='Mapa de coherencia recortada(BP)\n[dset'+str(i_o)+'-'+str(i_o+n_im-1)+']'
            direction ='CoherenciaCut_BP_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'.png'

        elif algorithm == "RMA":
            title_name ='Mapa de coherencia recortada(RMA)\n[dset'+str(i_o)+'-'+str(i_o+n_im-1)+']'
            direction ='CoherenciaCut_RMA_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'.png'
        #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
        #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20
        fig, ax = plt.subplots()
        im=ax.imshow(coh,cmap,origin='lower',aspect='equal', extent=[data['x_min'],data['x_max'],data['y_min'],data['y_max']]) #, vmin=vmin, vmax=vmax)
        ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
        ax.grid()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
        plt.colorbar(im,cax=cax,label='',extend='both')
        fig.savefig(directory_coh+direction, orientation='landscape')

    # Hallando el factor de calibracion
        # 1) Create a mask with the highest value of coherence in 1 
    M = (coh2>=coh2.max()) # bool type
    M = M*1 # int type
        # ...This part helps to choose randonly only a value if there's many maximums  
    M_aux = np.linspace(1,100,len(M)*len(M[0]))
    random.shuffle(M_aux)
    M_aux = M_aux.reshape(len(M),len(M[0]))
    M = M*M_aux
    M = (M==M.max())
    M = M*1
    
    # Obteniendo el mapa de desplazamientos, aka, Interferograma
    par = sp.get_parameters()
    global c,fc
    c,fc = par['c'],par['fc']

    if show:
        for i in range(n_im-1):
            # Hallando el interferograma 'i'-esimo
            i1 = i+i_o
            Im1 = np.load(directory2+"/Im_"+str(i1)+".npy")
            Im1 = calibration(M,Im1)
            
            Im2 = np.load(directory2+"/Im_"+str(i1+1)+".npy")
            Im2 = calibration(M,Im2)
            
            disp = np.angle(Im1*Im2.conjugate())*1000*c/(4*np.pi*fc) # Distancias en mm
            disp[mask] = np.nan
            
            # Graficando y guardando
            if algorithm == "BP":
                title_name ='Interferograma (BP) \ndset'+'['+str(i1)+']-['+str(i1+1)+']'
                #direction ='Interferogramas_BP/'+'Itf_BP_dset'+'['+str(i+10)+']-['+str(i+11)+'].png'
                direction ='Itf_BP_dset'+str(i)+'.png'

            elif algorithm == "RMA":
                title_name ='Interferograma (RMA) \ndset'+'['+str(i1)+']-['+str(i1+1)+']'
                #direction ='Interferogramas_RMA/'+'Itf_RMA_dset'+'['+str(i+10)+']-['+str(i+11)+'].png'
                direction ='Itf_RMA_dset'+str(i)+'.png'
            cmap = plt.cm.plasma #brg
            cmap.set_bad('black',1.)
            #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
            #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20
            fig, ax = plt.subplots()
            im=ax.imshow(disp,cmap,origin='lower',aspect='equal', extent=[data['x_min'],data['x_max'],data['y_min'],data['y_max']]) #, vmin=vmin, vmax=vmax)
            ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
            ax.grid()
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
            plt.colorbar(im,cax=cax,label='desplazamiento LOS(mm)',extend='both')
            
            directory_itc = directory0 + "/Interferograms/Interferograms_complete/Images/"
            if not os.path.exists(directory_itc):
                os.makedirs(directory_itc)
            
            fig.savefig(directory_itc+direction,orientation='landscape')
            plt.close()

    # Hallando la curva de distancias vs tiempo
    #-------Definiendo las zonas--------------
    zone0 = np.array([0,len(Im2),0,len(Im2[0])])# ([0,700,0,600]) Toda la imagen
    zone1 = np.array([0,200,100,300]) # Parte de la imagen
    zone2 = np.array([0,200,300,500])
    zone3 = np.array([200,400,100,300])
    zone4 = np.array([200,400,300,500])
    zone5 = np.array([400,600,100,300])
    zone6 = np.array([400,600,300,500])
    zone_indexes = {0:zone0,1:zone1,2:zone2,3:zone3,4:zone4,5:zone5,6:zone6}

    desp = np.zeros((len(zone_indexes),n_im-1)) # Variable y: desplazamiento
    desp_acc = np.zeros((len(zone_indexes),n_im-1)) # Variable y: desplazamiento acumulado
    desv_std = np.zeros((len(zone_indexes),n_im-1)) # Standar Deviation
    #mag = np.zeros((len(zone_indexes),n_im-1)) # Image magnitud
    snr = np.zeros((len(zone_indexes),n_im-1)) # SNR
    
    for z in range(len(zone_indexes)): # Hallando los desplazamientos promedios acumulados por zona
        idc = zone_indexes[z]
        for i in range(n_im-1):
            i1 = i+i_o # 10 es el valor inicial de las imagenes
            Im1 = np.load(directory2+"/Im_"+str(i1)+".npy")
            Im1 = calibration(M,Im1)
                        
            Im2 = np.load(directory2+"/Im_"+str(i1+1)+".npy")
            Im2 = calibration(M,Im2)
            
            d_i = np.angle(Im1*Im2.conjugate())*1000*c/(4*np.pi*fc) # dist(mm)
            
            # ----- Grafica ------
            #d_i2 = d_i.copy()
            #d_i2[mask] = np.nan
            #d_i2 = d_i2[idc[0]:idc[1],idc[2]:idc[3]]
            # -----   end --------
            d_i[mask] = np.nan
            d_i = d_i[idc[0]:idc[1],idc[2]:idc[3]]
            
            factor = 20*np.log10(np.mean(abs(Im1[550:650,0:200]))) # mean value of error 
            aux_snr = 20*np.log10(abs(Im1)) # snr of the image
            aux_snr = aux_snr - factor
            
            aux_snr[mask] = np.nan
            aux_snr = aux_snr[idc[0]:idc[1],idc[2]:idc[3]]
            
            if i==0: 
                desp[z,i] = np.nanmean(d_i) # mean ignoring nan values
                desp_acc[z,i] = np.nanmean(d_i)
                desv_std[z,i] = np.nanstd(d_i)
                #mag[z,i] = np.mean(abs(Im1))
                snr[z,i] = np.nanmean(aux_snr)
            else: 
                desp[z,i] = np.nanmean(d_i)
                desp_acc[z,i] = np.nanmean(d_i)+desp_acc[z,i-1]
                desv_std[z,i] = np.nanstd(d_i)
                #mag[z,i] = np.mean(abs(Im1))
                snr[z,i] = np.nanmean(aux_snr)
                
            # Graficando la zona z-esima
            if show and i == 0: #
                if algorithm == "BP":
                    title_name ='Interferograma recortado(BP) \nZona '+str(z)+' - dset'+'['+str(i1)+']-['+str(i1+1)+']'
                    direction ='Itf_BP_rec_zona_'+str(z)+'_dset'+'['+str(i1)+']-['+str(i1+1)+'].png'
                    #direction ='Interferogramas_BP/'+'Itf_BP_dset'+str(i)+'.png'

                elif algorithm == "RMA":
                    title_name ='Interferograma recortado(RMA) \nZona '+str(z)+' - dset'+'['+str(i1)+']-['+str(i1+1)+']'
                    direction ='Itf_RMA_rec_zona_'+str(z)+'_dset'+'['+str(i1)+']-['+str(i1+1)+'].png'
                    #direction ='Interferogramas_RMA/'+'Itf_RMA_dset'+str(i)+'.png'
                cmap = plt.cm.plasma #brg
                cmap.set_bad('black',1.)
                #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
                #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20
                fig, ax = plt.subplots()
                im=ax.imshow(d_i,cmap,origin='lower',aspect='equal', extent=[idc[2]-300,idc[3]-300,idc[0]+100,idc[1]+100]) #, vmin=vmin, vmax=vmax)
                ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
                ax.grid()
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
                plt.colorbar(im,cax=cax,label='desplazamiento LOS(mm)',extend='both')
                
                directory_its = directory0 + "/Interferograms/Interferograms_section/Images/"
                if not os.path.exists(directory_its):
                    os.makedirs(directory_its)
            
                fig.savefig(directory_its+direction,orientation='landscape')

    #t = np.arange(len(Ims)-1)+10 # Variable x, tiempo
    time_dset = np.load(directory1 + "/Dates.npy")
    time_dset = time_dset[:-1]
    time_dset2 = np.array([datetime.datetime.strptime(idx,"%d-%m-%y %H:%M:%S") for idx in time_dset])

    # Grafica de desplazamiento vs tiempo por cada algoritmo
    """if show:
        if algorithm == "BP":
            title_name ='Desplazamiento promedio(BP)'
            direction ='desplazamientoPromedio_BP_dset10-100.png'

            fig, ax= plt.subplots()
            ax.plot(t,desp,'b',marker='o',markerfacecolor='b',markeredgecolor='b')
            ax.set(xlabel='Tiempo(unidades)',ylabel='Desplazamiento(mm)', title=title_name)
            #ax.set_xlim([R.min(),R.max()])
            ax.set_ylim([-c*100/(4*fc),c*100/(4*fc)])
            ax.grid()
            plt.show()
            fig.savefig(os.getcwd()+"/Results/Desplazamientos/"+direction,orientation='landscape')

        elif algorithm == "RMA":
            title_name ='Desplazamiento promedio(RMA)'
            direction ='desplazamientoPromedio_RMA_dset10-100.png'

            fig, ax= plt.subplots()
            ax.plot(t,desp,'r',marker='o',markerfacecolor='r',markeredgecolor='r')
            ax.set(xlabel='Tiempo(unidades)',ylabel='Desplazamiento(mm)', title=title_name)
            #ax.set_xlim([R.min(),R.max()])
            ax.set_ylim([-c*50/(4*fc),c*50/(4*fc)])
            ax.grid()
            plt.show()
            fig.savefig(os.getcwd()+"/Results/Desplazamientos/"+direction,orientation='landscape')
      """
    return {'disp':desp, 'disp_acc':desp_acc, 'std':desv_std,'snr':snr,'t':time_dset2}

def get_displacements(ds):
    """ En esta funcion se grafican las curvas de desplazamientos promedios.
    """
    # Se obtienen una matriz de datos con los desplazamientos promedios de cada imagen
    t = ds['t']
    t = t[:n_im-1]
    t = mplt.dates.date2num(t)
    d = ds['disp']
    d_acc = ds['disp_acc']
    d_std = ds['std']
    snr = ds['snr']
    
    directory_stat = directory0 + "/Statistical_graphs/"
    if not os.path.exists(directory_stat):
        os.makedirs(directory_stat)
    
    f = open(directory_stat+"Log_statistic.txt","w+")
    f.write("-----------------------------------------------------------------\n")
    f.write("                 LANDSLIDE STATISTICS PARAMETERS\n")
    f.write("-----------------------------------------------------------------\n")
    f.write("Current time: "+time.strftime("%H:%M:%S - ")+time.strftime("%d/%m/%y")+"\n")
    f.write("-----------------------------------------------------------------\n")
    #f.write("\n")
    
    # Se grafica la curva Desplazamientos promedios vs Tiempo
    formatter = DateFormatter("%d/%m - %H:%M")
    for i in range(len(d_acc)):
        # Hallando el valor promedio final x zona
        mean_bp = d[i].mean()
        mean_acc_bp = d_acc[i].mean()
        mean_std = d_std[i].mean()
        mean_snr = snr[i].mean()
        f.write("Desplazamiento promedio BP_zona"+str(i)+"(mm): "+str(mean_bp)+"\n")
        f.write("Desplazamiento acumulado promedio BP_zona"+str(i)+"(mm): "+str(mean_acc_bp)+"\n")
        f.write("Desviacion estandar promedio BP_zona"+str(i)+"(mm): "+str(mean_std)+"\n")
        f.write("SNR promedio BP_zona"+str(i)+": "+str(mean_snr)+"\n")
        f.write("\n")
        
        # Graficando los desplazamientos 
        direction = 'Stat_graphs_zone'+str(i)+'.png'
        title_name = 'Statistical Results - Section '+str(i)

        fig, ax= plt.subplots(2,2,figsize=(11,8))
        ax[0,0].plot_date(t,d[i],'b',marker='',markerfacecolor='b',markeredgecolor='b') #,label='Back Projection')
        ax[0,0].set(xlabel = r"$time$",ylabel = r"$\overline{\Delta r} \ - \ LOS (mm)$",title= r"$\overline{\Delta r} \ vs \ time $") # \n(Zona "+str(i)+')')
        ax[0,0].xaxis.set_major_formatter(formatter)
        ax[0,0].xaxis.set_tick_params(rotation=20)
        #ax.set_xlim([R.min(),R.max()])
        ax[0,0].set_ylim([-c*1000/(4*fc*2),c*1000/(4*fc*2)]) # En (mm)
        ax[0,0].grid(linestyle='dashed')
        #ax[0,0].legend()
        #plt.show()
        #fig.savefig(os.getcwd()+"/Results/Displacement_BP/"+direction1,orientation='landscape')

        # Graficando los desplazamientos acumulados
        #direction2 = 'desplazamientosPromedios_acc_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'_zona'+str(i)+'.png'

        #fig, ax= plt.subplots(figsize=(10,7))
        ax[0,1].plot_date(t,d_acc[i],'b',marker='',markerfacecolor='b',markeredgecolor='b') #,label='Back Projection')
        ax[0,1].set(xlabel = r"$time$",ylabel = r"$\overline{\Delta r_{acc}} \ - \ LOS (mm)$",title= r"$\overline{\Delta r_{acc}} \ vs \ time $") # \n(Zona "+str(i)+')')
        ax[0,1].xaxis.set_major_formatter(formatter)
        ax[0,1].xaxis.set_tick_params(rotation=20)
        #ax.set_xlim([R.min(),R.max()])
        ax[0,1].set_ylim([-c*1000*2/(4*fc),c*1000*2/(4*fc)]) # En (mm)
        ax[0,1].grid(linestyle='dashed')
        #ax[0,1].legend()
        #plt.show()
        #fig.savefig(os.getcwd()+"/Results/Displacement_BP/"+direction2,orientation='landscape')
        
        # Graficando la desviacion estandar
        #direction3 = 'desviacionEstandar_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'_zona'+str(i)+'.png'

        #fig, ax= plt.subplots(figsize=(10,7))
        ax[1,0].plot_date(t,d_std[i],'b',marker='',markerfacecolor='b',markeredgecolor='b')#,label='Back Projection')
        ax[1,0].set(xlabel = r"$time$",ylabel = r"$\sigma_{\Delta r}(mm)$",title= r"$\sigma_{\Delta r} \ vs \ time $") # \n(Zona "+str(i)+')')
        ax[1,0].xaxis.set_major_formatter(formatter)
        ax[1,0].xaxis.set_tick_params(rotation=20)
        #ax.set_xlim([R.min(),R.max()])
        ax[1,0].set_ylim([0,c*1000/(4*fc)]) # En (mm)
        ax[1,0].grid(linestyle='dashed')
        #ax[1,0].legend()
        #plt.show()
        
        # Graficando el SNR
        ax[1,1].plot_date(t,snr[i],'b',marker='',markerfacecolor='b',markeredgecolor='b')#,label='Back Projection')
        ax[1,1].set(xlabel = r"$time$",ylabel = r"$\overline{snr}(dB)$",title= r"$\overline{snr} \ vs \ time $") # \n(Zona "+str(i)+')')
        ax[1,1].xaxis.set_major_formatter(formatter)
        ax[1,1].xaxis.set_tick_params(rotation=20)
        ax[1,1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        ax[1,1].yaxis.major.formatter._useMathText = True
        #ax.set_xlim([R.min(),R.max()])
        #ax[1,1].set_ylim([0,c*1000/(4*fc)]) # En (mm)
        ax[1,1].grid(linestyle='dashed')
        #ax[1,1].legend()
        
        #fig.subplots_adjust(left=0.065, right=0.95, wspace=0.3)
        fig.suptitle(title_name, y=0.97)
        fig.tight_layout(w_pad=2,h_pad=2,rect=[0,0,1,0.95])
            
        fig.savefig(directory_stat+direction,orientation='landscape')
        
        plt.close('all')
        
    f.close()
        
    return 'Ok'

def compare_displacements(ds1,ds2):
    """ En esta funcion se grafican las curvas de desplazamientos promedios para cada algoritmo y 
        por cada zona. Luego de ello se comparan ambas curvas graficandolos en una sola grafica. 
    """
    # Obteniendo los datos para BP
    t1 = ds1['t']
    t1 = t1[:n_im-1]
    t1 = mplt.dates.date2num(t1)
    d1 = ds1['d_t']
    # Obteniendo los datos para RMA
    t2 = ds2['t']
    t2 = t2[:n_im-1]
    t2 = mplt.dates.date2num(t2)
    d2 = ds2['d_t']

    # Graficando las 2 curvas juntas
    formatter = DateFormatter("%d/%m - %H:%M")
    for i in range(len(d1)):
        # Hallando el valor promedio final x zona
        mean_bp = d1[i].mean()
        mean_rma = d2[i].mean()
        print("Valor promedio BP_zona"+str(i)+": ",mean_bp)
        print("Valor promedio RMA_zona"+str(i)+": ",mean_rma)
        print("")
        # Graficando
        direction = 'desplazamientosPromedios_dset'+str(i_o)+'-'+str(i_o+n_im-1)+'_zona'+str(i)

        fig, ax= plt.subplots(figsize=(10,7))
        ax.plot_date(t1,d1[i],'b',marker='',markerfacecolor='b',markeredgecolor='b',label='Back Projection')
        ax.plot_date(t2,d2[i],'r',marker='',markerfacecolor='r',markeredgecolor='r',label='RMA')
        ax.set(xlabel='Tiempo',ylabel='Desplazamiento(mm)',title="Desplazamientos promedios\n(Zona "+str(i)+')')
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_tick_params(rotation=20)
        #ax.set_xlim([R.min(),R.max()])
        ax.set_ylim([-c*1000*4/(4*fc),c*1000*4/(4*fc)])
        ax.grid(linestyle='dashed')
        ax.legend()
        plt.show()
        fig.savefig(os.getcwd()+"/Results/Desplazamientos/"+direction,orientation='landscape')

    return 'Ok'

def main():
    plt.close('all')
    start_time = timeit.default_timer() 

    # Se hallan las imagenes SAR
    get_images(algorithm='BP')
    
    # Se cargan los parametros del Imaging-SAR
    data1 = np.load(directory1+"/Parameters.npy").item()

    # Se calculan los interferogramas
    disp1 = make_interferometry(data1,algorithm='BP')
    
    # Se obtienen las curvas de analisis estadistico
    get_displacements(disp1)
    
    h = open(directory0 + "/Statistical_graphs/Log_statistic.txt","a+")
    h.write("-----------------------------------------------------------------\n")
    h.write("Tiempo total de procesamiento(s): " + str(timeit.default_timer() - start_time))
    h.close()
        
    return 'Ok'

if __name__ == '__main__':
    main()
