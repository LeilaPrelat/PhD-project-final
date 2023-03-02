
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
from scipy.signal import find_peaks
import seaborn as sns
sns.set()

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

err = 'decay_rate_film.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_film import EELS_film_ana_f_div_gamma0, EELS_film_num_f_div_gamma0, EELS_film_pole_aprox_f_div_gamma0
except ModuleNotFoundError:
    print(err)



try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x,hBn_lambda_p
except ModuleNotFoundError:
    print('epsilon_z.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
except ModuleNotFoundError:
    print('Silica_epsilon.py no se encuentra en ' + path_basic)
    

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


x1 = 0.09260651629072682 
x2 = 0.10112781954887218

x3 = 0.17030075187969923 # no se puede hacer el 2do intervalo porque Im{Gself} es negativo 
x4 = 0.19937343358395992

x0 = 0.08684904004367106 
xf =  0.8064559363534132

x0 = x1 + 1e-3
xf = x4 - 1e-3

b = -0.01

d_nano_film = 1

d_thickness_disk_nano = 1
D_disk_nano = 100



#### la parte real de lambda_p para d = 0.4 nm es positiva si #####

# 0.08684904004367106 intervalo 1   # 0.093 - 0.1 esta dentro de x1-x2 pero no funciona porque Im{Gself} es negativo ---> el decay rate es negativo 
# 0.13150932790273412

# 0.15239881738519911 intervalo 2   # 0.153 - 0.169 esta dentro de x3-x4 . funciona desde 0.141 eV hasta 0.162 eV (entre 0.163 y 0.169 Im{Gself} es negativo ---> el decay rate es negativo )
# 0.16968667074999771

# 0.20714368637372804 intervalo 3 # es mayor a x4 
# 0.8064559363534132


#########################

def lambda_p(energy0):
    
    epsi_x = epsilon_x(energy0)
    epsi_HBN_par = epsi_x     
    
    omegac = energy0/aux
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsilon_Silica(energy0)*2/(omegac*d_micro*(1 - epsi_HBN_par))
    kp = alfa_p*omegac
      

    return (2*np.pi/kp)*1e3 ## en micro


def lambda_p_v2(energy0):
    
#    d_micros = d_nano*1e-3
    lambda_p_v = hBn_lambda_p(energy0,epsilon_Silica(energy0),epsilon_Silica(energy0))*d_nano_film

    return lambda_p_v ## en nano

#
#E0 = 0.093 # eV
##    E0 = xf # eV
##    E0 = 0.143
#E0 = 0.1625

E0 = 0.171
E0 = 0.195
int_v0 = 10
lambbda_p = np.real(lambda_p(E0))

labelx = r'$z_{\rm 0}$ (nm)'   


 
#title1 = r'$\kappa$ = %.2f$\omega_0$, $\kappa_r$ = %.2f$\kappa$, $E_0$=%i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)     
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'$v$ = c/%i, $b$ = %i nm, $d_{\rm film}$ = %.2f nm' %(int_v0,b*1e3,d_nano_film)
title2 = r'$d_{\rm disk}$ = %.2f nm, $D$ = %.2f nm, $\hbar\omega$ = %.3f eV' %(d_thickness_disk_nano,D_disk_nano,E0)
labelp = r'_dfilm%.2fnm_ddisk%.2fnm_D%inm' %(d_nano_film,d_thickness_disk_nano,D_disk_nano) 
label1 = 'vs_zp_E%imeV' %(E0*1e3) + labelp
title =  title1 + '\n' + title2

N = 30

#%%

def function_imag_ana(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_ana_f_div_gamma0(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    return rta


def function_imag_num(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_num_f_div_gamma0(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    return rta


def function_imag_pole_aprox(energy0,int_v,zp_nano):
    omegac0 = energy0/aux 
    zp = zp_nano*1e-3

    rta = EELS_film_pole_aprox_f_div_gamma0(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    return rta



#%%
    
#
#if plot_vs_c == 1 :
#    E0 = 44 # meV
#    # z0 = 0.06*1e3
#    zp0 = 0.05*1e3 
#    
#    labelx = r'v/c'   
#    title4 = title4 + ', ' + r'$\hbar\omega$ = %.2f eV, $z_p$ = %i nm' %(E0,zp0)
#    label1 = 'vs_v' + labelp
##    listx = np.linspace(0.0001,2,N)
#    listx = np.linspace(0.005,0.12,N)
#
#if plot_vs_E ==1 :
#    # z0 = 0.06*1e3
#    int_v0 = 10
#    zp0 = 0.05*1e3 
#    
#    labelx = r'$\hbar\omega$ [eV]'   
#    title4 = title4 + ', ' + r'v = c/%i, $z_p$ = %i nm' %(int_v0,zp0)
#    label1 = 'vs_E' + labelp
##    listx = np.linspace(0.0001,2,N)
#    listx = np.linspace(15,65,N)
#
#
#

#
#if plot_vs_zp == 1 : 

#    listx = np.linspace(0.0001,2,N)
if d_nano_film == 1:
    if E0 <= 0.187:
        listx = np.linspace(1,150,N)
    else:
        listx = np.linspace(0.1,20,N)

else:
    listx = np.linspace(1,800,N)
    
#if D_disk_nano == 50:
#    
#    listx = np.linspace(0.001,3.5,N)
#else:
#        listx = np.linspace(1,10,N)
#    listx = np.linspace(400,1600,N)
#    listx = np.linspace(250,750,N)
#    
    
    

#%%
    
tamfig = [4.5, 4]
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
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   

#%%


listy_im_ana = []
listy_im_num = []
listy_im_pole_aprox = []

#if plot_vs_E == 1: 
#
#    for value in listx: 
#
#        y_im_ana = function_imag_ana(value,int_v0,zp0)        
#        listy_im_ana.append(y_im_ana)
#        
#        y_im_num = function_imag_num(value,int_v0,zp0)        
#        listy_im_num.append(y_im_num)
#        
#        y_im_pole_aprox = function_imag_pole_aprox(value,int_v0,zp0)        
#        listy_im_pole_aprox.append(y_im_pole_aprox)
#        
#               
#elif plot_vs_c == 1:       
#
#    
#    for value in listx: 
#        
#        value2 = 1/value
##        print(value2)
#
#        y_im_ana = function_imag_ana(E0,value2,zp0)        
#        listy_im_ana.append(y_im_ana)
#
#        y_im_num = function_imag_num(E0,value2,zp0)        
#        listy_im_num.append(y_im_num)
#        
#        y_im_pole_aprox = function_imag_pole_aprox(E0,value2,zp0)        
#        listy_im_pole_aprox.append(y_im_pole_aprox)
#
#        
#elif plot_vs_zp == 1:

    
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
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'darkred',label = labell2)
#plt.plot(listx,list_ana_parallel,'.-',ms = ms,color = 'darkred',label = r'$\Gamma_{\parallel}$')
plt.plot(np.array(np.ones(10))*maxi, listy_aux,'-k')
plt.legend(loc = 'best',markerscale=1.5,fontsize=tamlegend,frameon=0.2,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#if plot_vs_c == 1:
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'EELS_film_' + label1 + '.png', bbox_inches='tight',pad_inches = 0.01, format='png', dpi=dpi)    



#%%


