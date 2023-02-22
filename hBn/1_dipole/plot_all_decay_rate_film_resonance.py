
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo directo con la convencion de z hacia abajo
en z = 0
graficar mapa de color x,y
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
from scipy import special
from scipy.signal import find_peaks
#import seaborn as sns
#sns.set()

#%%
plot_vs_E = 0
plot_vs_c = 0
plot_vs_zp = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/1_dipole','')
path_save = path_basic + '/' + 'decay_rate_film'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_film_resonance.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film_resonance import EELS_film_ana_f_div_gamma0, EELS_film_num_f_div_gamma0, EELS_film_pole_aprox_f_div_gamma0
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

epsi1,epsi3 = 1,1

zp = 0.05
b = -0.01

d_nano = 1
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title4 = r'b = %i nm' %(b*1e3)
labelp = r'_res' 

N = 100



def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta1 = EELS_film_ana_f_div_gamma0(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)

    return rta1



def function_imag_num(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta1 = EELS_film_num_f_div_gamma0(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)

    return rta1


def function_imag_pole_aprox(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta1 = EELS_film_pole_aprox_f_div_gamma0(omegac0,epsi1,epsi3,d_nano,int_v,b,zp)

    return rta1



if plot_vs_c == 1 :
    E0 = 44 # meV
    # z0 = 0.06*1e3
    zp0 = 0.05*1e3 
    
    labelx = r'v/c'   
    title4 = title4 + ', ' + r'$\hbar\omega$ = %i meV, $z_p$ = %i nm' %(E0,zp0)
    label1 = 'vs_v' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(0.005,0.12,N)

if plot_vs_E ==1 :
    # z0 = 0.06*1e3
    int_v0 = 10
    zp0 = 0.05*1e3 
    
    labelx = r'$\hbar\omega$ [meV]'   
    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
    label1 = 'vs_E' + labelp
#    listx = np.linspace(0.0001,2,N)
    listx = np.linspace(15,65,N)

if plot_vs_zp == 1 : 
    int_v0 = 10 ## deberia ser 150 (disp relation) pero funciona con 10 <--- problema con la relacion de dispersion
    E0 = 0.175 # eV
#    E0 = 1.5

    labelx = r'Surface-dipole distance, $z_{\rm 0}$/$\lambda_{\rm p}$'   
    title4 = title4 + ', ' + r'v = c/%i, $\hbar\omega$ = %i eV' %(int_v0,E0)
    label1 = 'vs_zp' + labelp + '_E%imeV_d%inm' %(E0*1e3,d_nano)
#    listx = np.linspace(0.0001,2,N)
    if d_nano == 1:
        if E0 <= 0.187:
            listx = np.linspace(10,500,N)
        else:
            listx = np.linspace(25,400,N)

    else:
        listx = np.linspace(5,1050,N)
    
#    print(minimum_function(E0,int_v0)*1e3)
#    print(np.abs(minimum_function(E0,int_v0))*2*1e3)
    
title =  title4 


#%%
    
    
tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 8
tamnum = 7
tamlegend = 7
labelpady = 2
labelpadx = 3
pad = 3
mk = 1
ms = 2
hp = 0.3
length_marker = 1.5
dpi = 500

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.tick_params(labelsize = tamnum, length = 4 , width=1, direction="in", pad = pad) ## para el congreso de CLEO
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return   


#%%


listy_im_ana = []
listy_im_num = []
listy_im_pole_aprox = []


if plot_vs_E == 1: 

    for value in listx: 

        y_im_ana = function_imag_ana(value,int_v0,zp0)        
        listy_im_ana.append(y_im_ana)

        y_im_num = function_imag_num(value,int_v0,zp0)        
        listy_im_num.append(y_im_num)      

        y_im_pole_aprox = function_imag_pole_aprox(value,int_v0,zp0)        
        listy_im_pole_aprox.append(y_im_pole_aprox)      


elif plot_vs_c == 1:       

    for value in listx: 
        
        value2 = 1/value
#        print(value2)

        y_im_ana = function_imag_ana(E0,value2,zp0)        
        listy_im_ana.append(y_im_ana)

        y_im_num = function_imag_num(E0,value2,zp0)        
        listy_im_num.append(y_im_num)      

        y_im_pole_aprox = function_imag_pole_aprox(E0,value2,zp0)        
        listy_im_pole_aprox.append(y_im_pole_aprox)      

elif plot_vs_zp == 1:
    
    for value in listx: 

        y_im_ana = function_imag_ana(E0,int_v0,value)        
        listy_im_ana.append(y_im_ana)


        y_im_num = function_imag_num(E0,int_v0,value)        
        listy_im_num.append(y_im_num)      

        y_im_pole_aprox = function_imag_pole_aprox(E0,int_v0,value)        
        listy_im_pole_aprox.append(y_im_pole_aprox)      

        
#%%
peaks, _ = find_peaks(listy_im_ana, height=0)
if len(peaks) > 1:
    maxis = []
    for ind in peaks:
        maxis.append(listy_im_ana[ind])
#    print(maxis)
    maxi = listx[peaks[int(np.argmax(maxis))]]
    
else:
    maxi = listx[peaks]
listy_aux  = np.linspace(np.min(listy_im_ana), np.max(listy_im_ana), 10)
labell2 =  r'PP num $k_\parallel/(k_\parallel - k_{\rm p})$'

graph(title,labelx,r'$\Gamma_{SP}/\Gamma_{\rm EELS}$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_ana,'.-',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_im_num,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'darkred',label = 'PP numerical')
#plt.plot(listx,list_ana_parallel,'.-',ms = ms,color = 'darkred',label = r'$\Gamma_{\parallel}$')
plt.plot(np.array(np.ones(10))*maxi, listy_aux,'-k')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.2,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_film_' + label1 + '.png', bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)    


