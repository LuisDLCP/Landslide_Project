"""
Created on Mon Jun  3 12:23:13 2019
---------------------------------------------------
            INTERFEROMETRIA GENERAL
---------------------------------------------------
* Se obtendra el mapa de coherencia
* Se obtendra un grafico de desplazamientos a lo largo del tiempo
* Datos reales RMA y BP
* Se compararan los 2 graficos de desplazamientos vs tiempo, para RMA y BP
@author: LuisDLCP
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplt
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import RMA_real_main as RMA
import BP_real_main as BP
import sarPrm as sp
import timeit
import datetime
import os

os.chdir(os.path.dirname(__file__)) # get the current path

show = False
n_im = 3 #4656 # Numero de imagenes a considerar

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
            i += 10 # Empieza en la posicion 10
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
            i += 10 # Empieza en la posicion 10
            data = BP.main("dset_"+str(i)+".hdf5") 
            #Ims[i] = data['Im']
            dates.append(data['date'])
            np.save(os.getcwd()+"/Results/Output_BP/Im_"+str(i)+".npy",data['Im']) # Imagenes de todo el dataset
        np.save("Parameters_BP",data) # Parametros geometricos como dimensiones y grilla de la imagen
        np.save("Dates_BP",np.array(dates)) # Fechas de las iamgenes tomadas de todo el dset

    return 'Ok'

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
    Im = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_10.npy") # Se carga la imagen 10
    f1 = np.zeros((len(Im),len(Im.T)),dtype=complex)
    f2 = np.zeros((len(Im),len(Im.T)),dtype=complex)
    coh = np.zeros((len(Im),len(Im.T)),dtype=complex)

    for i in range(n_im-1):
        i += 10 # Valor inicial
        Im1 = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_"+str(i)+".npy")
        Im2 = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_"+str(i+1)+".npy")
        f1 += Im1*Im2.conjugate()
        f2 += np.sqrt((abs(Im1)**2)*(abs(Im2)**2))

    coh = abs(f1/f2)

    # Graficando la coherencia
    if show:
        plt.close('all')
        cmap ="plasma"
        if algorithm == "BP":
            title_name ='Mapa de coherencia(BP)\n[dset10-5000]'
            direction ='Interferogramas_BP/coherencia_BP_dset10-5000.png'

        elif algorithm == "RMA":
            title_name ='Mapa de coherencia(RMA)\n[dset10-5000]'
            direction ='Interferogramas_RMA/coherencia_RMA_dset10-5000.png'
        #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
        #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20

        fig, ax = plt.subplots()
        im=ax.imshow(coh,cmap,origin='lower',aspect='equal', extent=[data['x_min'],data['x_max'],data['y_min'],data['y_max']]) #, vmin=vmin, vmax=vmax)
        ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
        ax.grid()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
        plt.colorbar(im,cax=cax,label='',extend='both')
        fig.savefig(os.getcwd()+"/Results/"+direction, orientation='landscape')

    # Sacando una mascara de [0.7,1] para la coherencia
    mask = coh<=0.7
    coh[mask] = np.nan

    # Graficando con la mascara
    if show:
        cmap ="plasma"
        if algorithm == "BP":
            title_name ='Mapa de coherencia recortada(BP)\n[dset10-5000]'
            direction ='Interferogramas_BP/coherenciaCut_BP_dset10-5000.png'

        elif algorithm == "RMA":
            title_name ='Mapa de coherencia recortada(RMA)\n[dset10-5000]'
            direction ='Interferogramas_RMA/coherenciaCut_RMA_dset10-5000.png'
        #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
        #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20
        fig, ax = plt.subplots()
        im=ax.imshow(coh,cmap,origin='lower',aspect='equal', extent=[data['x_min'],data['x_max'],data['y_min'],data['y_max']]) #, vmin=vmin, vmax=vmax)
        ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
        ax.grid()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
        plt.colorbar(im,cax=cax,label='',extend='both')
        fig.savefig(os.getcwd()+"/Results/"+direction, orientation='landscape')

    # Obteniendo el mapa de distancias, aka, Interferograma
    par = sp.get_parameters()
    global c,fc
    c,fc = par['c'],par['fc']

    if show:
        for i in range(n_im-1):
            # Hallando el interferograma 'i'-esimo
            i1 = i+10
            Im1 = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_"+str(i1)+".npy")
            Im2 = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_"+str(i1+1)+".npy")
            disp = np.angle(Im1*Im2.conjugate())*1000*c/(4*np.pi*fc) # Distancias en mm
            disp[mask] = np.nan
            # Graficando y guardando
            if algorithm == "BP":
                title_name ='Interferograma (BP) \ndset'+'['+str(i1)+']-['+str(i1+1)+']'
                #direction ='Interferogramas_BP/'+'Itf_BP_dset'+'['+str(i+10)+']-['+str(i+11)+'].png'
                direction ='Interferogramas_BP/'+'Itf_BP_dset'+str(i)+'.png'

            elif algorithm == "RMA":
                title_name ='Interferograma (RMA) \ndset'+'['+str(i1)+']-['+str(i1+1)+']'
                #direction ='Interferogramas_RMA/'+'Itf_RMA_dset'+'['+str(i+10)+']-['+str(i+11)+'].png'
                direction ='Interferogramas_RMA/'+'Itf_RMA_dset'+str(i)+'.png'
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
            plt.colorbar(im,cax=cax,label='desplazamiento(mm)',extend='both')
            fig.savefig(os.getcwd()+"/Results/"+direction,orientation='landscape')
            plt.close()

    # Hallando la curva de distancias vs tiempo
    #-------Definiendo las zonas--------------
    zone0 = np.array([0,650,0,400])
    zone1 = np.array([0,100,100,200])
    zone2 = np.array([0,100,200,300])
    zone3 = np.array([100,200,100,200])
    zone4 = np.array([100,200,200,300])
    zone5 = np.array([200,300,100,200])
    zone6 = np.array([200,300,200,300])
    zone_indexes = {0:zone0,1:zone1,2:zone2,3:zone3,4:zone4,5:zone5,6:zone6}

    desp = np.zeros((len(zone_indexes),n_im-1)) # Variable y: desplazamiento

    for z in range(len(zone_indexes)): # Hallando los desplazamientos promedios por zona
        idc = zone_indexes[z]
        for i in range(n_im-1):
            i1 = i+10 # 10 es el valor inicial de las imagenes
            Im1 = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_"+str(i1)+".npy")
            Im2 = np.load(os.getcwd()+"/Results/Output_"+algorithm+"/Im_"+str(i1+1)+".npy")
            d_i = np.angle(Im1*Im2.conjugate())*1000*c/(4*np.pi*fc)
            # ----- Grafica ------
            d_i2 = d_i.copy()
            d_i2[mask] = np.nan
            d_i2 = d_i2[idc[0]:idc[1],idc[2]:idc[3]]
            # -----   end --------
            d_i[mask] = 0
            d_i = d_i[idc[0]:idc[1],idc[2]:idc[3]]
            if i==0: desp[z,i] = d_i.mean()
            else: desp[z,i] = d_i.mean()+desp[z,i-1]

            # Graficando la zona z-esima
            if show and i == 0: #
                if algorithm == "BP":
                    title_name ='Interferograma recortado(BP) \nZona '+str(z)+' - dset'+'['+str(i1)+']-['+str(i1+1)+']'
                    direction ='Interferogramas_BP/'+'Itf_BP_rec_zona_'+str(z)+'_dset'+'['+str(i1)+']-['+str(i1+1)+'].png'
                    #direction ='Interferogramas_BP/'+'Itf_BP_dset'+str(i)+'.png'

                elif algorithm == "RMA":
                    title_name ='Interferograma recortado(RMA) \nZona '+str(z)+' - dset'+'['+str(i1)+']-['+str(i1+1)+']'
                    direction ='Interferogramas_RMA/'+'Itf_RMA_rec_zona_'+str(z)+'_dset'+'['+str(i1)+']-['+str(i1+1)+'].png'
                    #direction ='Interferogramas_RMA/'+'Itf_RMA_dset'+str(i)+'.png'
                cmap = plt.cm.plasma #brg
                cmap.set_bad('black',1.)
                #vmin = np.amin(20*np.log10(abs(Sf_n)))+55 #dB
                #vmax = np.amax(20*np.log10(abs(Sf_n)))#-20
                fig, ax = plt.subplots()
                im=ax.imshow(d_i2,cmap,origin='lower',aspect='equal', extent=[idc[2]-200,idc[3]-200,idc[0]+300,idc[1]+300]) #, vmin=vmin, vmax=vmax)
                ax.set(xlabel='Azimut(m)',ylabel='Rango(m)', title=title_name)
                ax.grid()
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.1) # pad es el espaciado con la grafica principal
                plt.colorbar(im,cax=cax,label='desplazamiento(mm)',extend='both')
                fig.savefig(os.getcwd()+"/Results/"+direction,orientation='landscape')

    #t = np.arange(len(Ims)-1)+10 # Variable x, tiempo
    time_dset = np.load("Dates_BP.npy")
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
    return {'d_t':desp, 't':time_dset2}

def compare_displacement(ds1,ds2):
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
        direction = 'desplazamientosPromedios_dset10-5000_zona'+str(i)

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

    # Se calcula la data del algoritmo Back Projection
    get_images(algorithm='BP')

    # Se calcula la data de algoritmo RMA
    get_images(algorithm='RMA')

    # Se carga la data del algoritmo Back Projection
    #Ims1 = np.load("Set_images_BP.npy").item()
    data1 = np.load("Parameters_BP.npy").item()

    # Se carga la data de algoritmo RMA
    #Ims2 = np.load("Set_images_RMA.npy").item()
    data2 = np.load("Parameters_RMA.npy").item()

    # Se calculan los interferogramas
    disp1 = make_interferometry(data1,algorithm='BP')
    disp2 = make_interferometry(data2,algorithm='RMA')

    # Se comaparan las curvas de desplazamientos de los 2 algoritmos
    compare_displacement(disp1,disp2)
    print("\nTiempo total: ",timeit.default_timer() - start_time,"s")

    return 'Ok'

if __name__ == '__main__':
    main()
