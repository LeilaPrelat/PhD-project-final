
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
import seaborn as sns
sns.set()

plot_vs_E = 1
plot_vs_zp = 0
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'Green_self_image'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'green_self_image.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num, green_self_ana_residuos,green_self_pole_aprox, green_self_ana_exponential_function
except ModuleNotFoundError:
    print(err)

try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb
aux2 = 1e12/c

#%%

print('Definir parametros del problema')

x1 = 0.09260651629072682 
x2 = 0.10112781954887218


x3 = 0.17030075187969923
x4 = 0.19937343358395992

N = 50

epsi1,epsi3 = 1,1
d_nano = 1

if plot_vs_E == 1:
    zp0_nano = 0.05*1e3
    listx =  np.linspace(0.093,0.195,N)
    title = r'$z_0$ = %i nm, d = %.2f nm' %(zp0_nano,d_nano)
    labelx = r'$\hbar\omega$ (eV)'   
    label1 = '_vs_E_d%.2fnm' %(d_nano) 


elif plot_vs_zp == 1:
    E0 = 0.175
    listx = np.linspace(10,400,N)
    title = r'$\hbar\omega$ = %i eV, d = %.2f nm' %(E0,d_nano)
    labelx = r'$z_0$ (nm)'   
    label1 = '_vs_zp_d%.2fnm' %(d_nano) 



def function_num_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp)
    # print(rta    
    return np.real(rtaself_x)

def function_num_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z = green_self_num(omegac0,epsi1,epsi3,d_nano,zp)
    # print(rta    
    return np.imag(rtaself_x)


#%%

def function_ana_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_residuos(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)

def function_ana_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_residuos(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)

#%%

def function_ana2_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 
    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_exponential_function(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_ana2_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 
    rtaself_x, rtaself_y, rtaself_z  = green_self_ana_exponential_function(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)



#%%
    
def function_pole_aprox_xx_re(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.real(rtaself_x)



def function_pole_aprox_xx_im(energy0,zp_nano):
    zp = zp_nano*1e-3

    omegac0 = energy0/aux 

    rtaself_x, rtaself_y, rtaself_z  = green_self_pole_aprox(omegac0,epsi1,epsi3,d_nano,zp)
    
    return np.imag(rtaself_x)



#%%
    
tamfig = (4.5,3.5)
tamlegend = 12
tamletra = 12
tamtitle = 11
tamnum = 10
labelpady = -1.5
labelpadx = 2
pad = 0
mk = 2
ms = 4
hp = 0.3
length_marker = 1

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))

    return   


#%%
    
listy_re_ana = []
listy_re_ana2 = []
listy_re_num = []
listy_re_pole_aprox = []


listy_im_ana = []
listy_im_ana2 = []
listy_im_num = []
listy_im_pole_aprox = []

if plot_vs_E == 1:

    for value in listx: 
    
        y_re_ana = function_ana_xx_re(value,zp0_nano)       
        y_re_ana2 = function_ana2_xx_re(value,zp0_nano)        
        y_re_num = function_num_xx_re(value,zp0_nano)
        y_re_pole_aprox = function_pole_aprox_xx_re(value,zp0_nano)
        
        listy_re_ana.append(y_re_ana)
        listy_re_ana2.append(y_re_ana2)
        listy_re_num.append(y_re_num)
        listy_re_pole_aprox.append(y_re_pole_aprox)
    
        y_im_ana = function_ana_xx_im(value,zp0_nano)       
        y_im_ana2 = function_ana2_xx_im(value,zp0_nano)   
        y_im_num = function_num_xx_im(value,zp0_nano)
        y_im_pole_aprox = function_pole_aprox_xx_im(value,zp0_nano)
    
        listy_im_ana.append(y_im_ana)
        listy_im_ana2.append(y_im_ana2)
        listy_im_num.append(y_im_num)
        listy_im_pole_aprox.append(y_im_pole_aprox)
    
    
#%%
    
        
labell2 =  r'PP num $R_{\rm p}k_\parallel/(k_\parallel - k_{\rm p})$'

graph(title,labelx,r'Re{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_re_ana3,'.-',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'PP ana residuos')
plt.plot(listx,listy_re_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole_aprox,'.-',ms = 3,color = 'darkred',label = labell2)
plt.plot(listx,listy_re_ana2,'.-',ms = ms,color = 'darkorange',label = 'PP ana exp func')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_Gself' + label1 + '.png', format='png')   


graph(title,labelx,r'Im{G$_{self}$} ($\mu$m$^{-3}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#plt.plot(listx,listy_im_ana3,'.-',ms = ms,color = 'blue',label = 'PP analytical 3')
#plt.plot(listx,listy_im_ana,'.',ms = ms,color = 'purple',label = 'PP ana residuos')
plt.plot(listx,listy_im_num,'.-',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_im_pole_aprox,'.-',ms = 3,color = 'darkred',label = labell2)
plt.plot(listx,listy_im_ana2,'.-',ms = ms,color = 'darkorange',label = 'PP ana exp func')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#    plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_Gself' + label1 + '.png', format='png')   


#%%

#%%