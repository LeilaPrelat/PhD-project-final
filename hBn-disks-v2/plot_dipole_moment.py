
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


comparar_otro_hBN = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'dipole_moment'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_ana, dipole_moment_pole_aprox, dipole_moment_num
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_basic)
    from Silica_epsilon import epsilon_Silica
except ModuleNotFoundError:
    print('Silica_epsilon.py no se encuentra en ' + path_basic)


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

int_v = 10
#v = c/int_v
#omega = 0.7*1e12

zp = 0.05
b = -0.01

d_thickness_disk_nano = 1
D_disk_nano = 100 

d_nano_film = 1

#omega0THz = 65
#omega0 = omega0THz*1e12 

 
    
#title2 = r'$\hbar\mu$ = %.2feV, $\hbar\gamma$ = %.4feV' %(hbmu,hbgama) 
#title3 = r'$z_p$=%inm, px=%i, py=%i, pz=%i' %(zp*1e3,px,py,pz)
title1 = r'$v$ = c/%i, $z_{\rm p}$ = %i nm, $b$ = %i nm, $d_{\rm film}$ = %.2f nm' %(int_v, zp*1e3,b*1e3,d_nano_film)
title2 = r'$d_{\rm disk}$ = %.2f nm, $D$ = %.2f nm' %(d_thickness_disk_nano,D_disk_nano)
labelp = r'_res_d%.2fnm' %(d_nano_film)

N = 75

    # z0 = 0.06*1e3
labelx = r'$\hbar\omega$ (eV)'   
labelp = 'vs_E' + labelp

x1 = 0.09260651629072682 
x2 = 0.10112781954887218
x3 = 0.17030075187969923
x4 = 0.19937343358395992

x0 = 0.08684904004367106 
xf =  0.8064559363534132

#listx = np.linspace(x1,x2,N)
listx = np.linspace(x3 ,x4,N)
#listx = np.linspace(x0+1e-3,xf - 1e-3,N)
#listx = np.linspace(x1,x2,N)

## entre x1 y x2, y entre x3 y x4 ##

title =  title1 + '\n' + title2  



#%%

def function_num(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f = dipole_moment_num(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    # print(rta)
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 

    return np.sqrt(rta)


#%%

def function_ana(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f  = dipole_moment_ana(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)

#%%

def function_pole_approx(energy0):
    omegac0 = energy0/aux 

    px_f,py_f,pz_f  = dipole_moment_pole_aprox(omegac0,epsilon_Silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)
    
    rta = np.abs(px_f)**2 + np.abs(py_f)**2 + np.abs(pz_f)**2 
    
    return np.sqrt(rta)


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
listy_re_anav2 = []
listy_re_num = []
listy_re_pole = []

for value in listx: 

#    y_re_ana = function_ana(value)       
    y_re_anav2 = function_ana(value)  
    
    y_re_num = function_num(value)
    y_re_pole = function_pole_approx(value)
    
#    listy_re_ana.append(y_re_ana)
    
    listy_re_anav2.append(y_re_anav2)

    listy_re_num.append(y_re_num)
    listy_re_pole.append(y_re_pole)

#Emax = listx[np.argmax(listy_re_num)]
#print(Emax)

#%%
labell2 =  r'PP num $R_{\rm p}k_\parallel/(k_\parallel - k_{\rm p})$'

graph(title,labelx,r'$|p|^2/(e/(2 \pi v))$',tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
#    plt.plot(listx,listy_re_ana,'.',ms = ms,color = 'purple',label = 'analytical')
plt.plot(listx,listy_re_anav2,'.',ms = ms,color = 'purple',label = 'PP analytical')
plt.plot(listx,listy_re_num,'.',ms = ms,color = 'lightseagreen',label = 'full numerical')
plt.plot(listx,listy_re_pole,'.-',ms = 3,color = 'darkred',label = labell2)
ejey_aux = np.linspace(np.min(listy_re_num), np.max(listy_re_num) , 10)
for x in [x3,x4]:
    plt.plot(x*np.ones(10), ejey_aux,'--',color = 'grey' )
plt.legend(loc = 'best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.2, handlelength=length_marker)
plt.tight_layout()
plt.yscale('log')
os.chdir(path_save)
plt.savefig( 'p_tot' + labelp + '.png', format='png')   



#%%
