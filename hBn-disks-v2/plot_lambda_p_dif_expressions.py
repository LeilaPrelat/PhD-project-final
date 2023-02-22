
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

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'lambda_p_dif_expressions'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
except ModuleNotFoundError:
    print('Silica_epsilon.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from hBn_PP import hBn_lambda_p, epsilon_x
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_basic)

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


d_nano = 0.4
 
title = r'$d$ = %.2f nm' %(d_nano)
labelp = r'd%.1fnm' %(d_nano)

N = 100

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ (eV)'   
label1 = 'vs_energy_' + labelp
listx = np.linspace(x1 + 1e-3,0.195,N)

def function_kp_Leila(energy0):

    epsi_x = epsilon_x(energy0) ## paralelo 
#    epsi_z = epsilon_z(E)
    
    epsi_HBN_par = epsi_x
    d_micros = d_nano*1e-3
    sigma_2D = d_micros*( 1 - epsi_HBN_par )/(4*pi)  ## se cancela el i*omega del sigma 
    kp = epsilon_Silica(energy0)/(2*np.pi*sigma_2D)

    return kp


def function_kp_Eduardo(energy0):
    epsi1, epsi3 = 1,1
    d_micros = d_nano*1e-3
    lambda_p_v = hBn_lambda_p(energy0,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v

    return kp

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


listy_re_Edu = []
listy_re_Lei = []

listy_im_Edu = []
listy_im_Lei = []

for value in listx: 
    
    print(value)
 
    y_re_Edu = function_kp_Eduardo(value)
    y_re_Lei = function_kp_Leila(value)


    listy_im_Edu.append(np.imag(y_re_Edu))
    listy_im_Lei.append(np.imag(y_re_Lei))


    listy_re_Edu.append(np.real(y_re_Edu))
    listy_re_Lei.append(np.real(y_re_Lei))


#%%
    
graph(title,labelx,r'Im{$k_p$} ($\mu$m$^{-1}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_im_Edu,'.-',ms = ms,color = 'purple',label = r'$r_{\rm p}$ from paper 331')
plt.plot(listx,listy_im_Lei,'.',ms = ms+1,color = 'blue',label = r'$k_{\rm p} = i\omega\epsilon/(2\pi\sigma_{\rm 2D})$')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Im_rp_' + label1 + '.png', format='png')   


graph(title,labelx,r'Re{$k_p$} ($\mu$m$^{-1}$)',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
plt.plot(listx,listy_re_Edu,'.-',ms = ms,color = 'purple',label =  r'$r_{\rm p}$ from paper 331')
plt.plot(listx,listy_re_Lei,'.',ms = ms+1,color = 'blue',label =  r'$r_{\rm p} = i\omega\epsilon/(2\pi\sigma_{\rm 2D})$')
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0.05,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
#plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'Re_rp_' + label1 + '.png', format='png')   
