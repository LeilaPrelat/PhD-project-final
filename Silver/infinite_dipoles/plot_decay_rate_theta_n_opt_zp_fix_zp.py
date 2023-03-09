
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
from scipy.interpolate import interp1d
#import seaborn as sns
#sns.set()
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/infinite_dipoles','')
path_save = path_basic + '/' + 'decay_rate_theta_n'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)

err = 'decay_rate_theta_n.py no se encuentra en ' + path_basic
try:
    sys.path.insert(1, path_basic)
    from decay_rate_theta_n import decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v2, decay_rate_theta_inf_dipoles_ana_res_div_gamma0
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

epsi1, epsi3 = 1,1

print('Definir parametros del problema')

list_n = [0,1,2,3,4]   
b = - 0.01

d_nano = 1
int_v = 10
Nmax = list_n[-1]

labely = r'$\Gamma_{\rm SP, n}/\Gamma_{\rm EELS}$'
#labely = r'Emission probability (eV$^{-1}$)'

tabla = np.loadtxt('zp_optimum_for_decay_rate_resonance_Silver_d%inm_v%i.txt' %(d_nano,int_v), delimiter='\t', skiprows=1)
tabla = np.transpose(tabla)
[listx,listy,listz] = tabla

zp_nano = listy[10]

#zp_nano = listy[-1]
zp_nano = 35
#zp_nano = 50
omegac0_1 = np.max(listx)/(c*hb)
lambda_SP_1 = 2*np.pi/omegac0_1

omegac0_2 = np.min(listx)/(c*hb)
lambda_SP_2 = 2*np.pi/omegac0_2


a_min = np.real(lambda_SP_1)*Nmax/(int_v - 1)
a_max = np.real(lambda_SP_2)*Nmax/(int_v + 1)

a = np.mean([a_min,a_max])
a = a_min
#a = 5031*1e-3
#a = 5500*1e-3
a = 850*1e-3

a_nm = a*1e3


labelx = r'Surface-dipole distance, $z_{\rm 0}/\lambda_{\rm p}$'  
    
title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$' 

labelp = r'_a%inm_zp%inm_d%inm' %(a*1e3,zp_nano,d_nano)
#label1 = 'vs_zp_lambda_p'  + labelp


#elif theta_degree == 60:
#    lim1, lim2 = 2.6,2.9
    
#def find_nearest(array, value):
#    array = np.asarray(array)
#    idx = (np.abs(array - value)).argmin()
#    return array[idx],idx    
#   
#num1,ind1 = find_nearest(np.array(listx), value=lim1)
#num2,ind2 = find_nearest(np.array(listx), value=lim2)

f1 = interp1d(listx, listy)
f2 = interp1d(listx, listz)

N = 225
lim1,lim2 = 18,-60
lim1,lim2 = 28,-52
#lim1,lim2 = 0,-1
listx_2 = np.linspace(listx[lim1], listx[lim2], N)


listy_2 = f1(listx_2)
listz_2 = f2(listx_2)  



title2 = r'v = c/%i, b = %i nm' %(int_v,b*1e3) 
title3 = r'a = %i nm' %(a*1e3)
title4 = r', $z_p$ = $z^{opt}_p$'


#labelp = r'_a%inm_zp_opt' %(a*1e3)
#label1 = 'vs_zp'  + labelp

    

#title1 = r'$\kappa$ = %.2f$\omega_o$, $\kappa_r$ = %.2f$\kappa$, $\hbar\omega_o$ = %i meV' %(kappa_factor_omega0, kappa_r_factor, energy0_pol)   
 
title =  title2 + '\n' +  title3  + title4



#%%


def function_real_ana(energy0_meV,Nmax):
    
#    a = a_nano*1e-3
    omegac0 = energy0_meV/(c*hb)  
    zp = zp_nano*1e-3
         
    rta = decay_rate_theta_inf_dipoles_ana_res_div_gamma0(omegac0,epsi1,epsi3,d_nano,int_v,zp,a,b,Nmax)

    return rta


tamfig = [2.5, 2]
tamletra = 9
tamtitle  = 9
tamnum = 7
tamlegend = 7
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500


#%%
    


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, length = 2 , width=1, direction="in", pad = pad)
#    plt.title(title,fontsize=int(tamtitle*0.9))

    return  
 
#%%
    
from scipy.signal import find_peaks
from scipy.signal import savgol_filter

maxis_zp_lambda_p = []
values_maxis = []  
values_maxis_tot = []
#    if theta_degree != 0:     
#        listx_2 = listx
#        listy_2 = listy
#        listz_2 = listz
list_y_re_tot = []
for n in list_n:
    
    list_y_re = []


    for ind in range(len(listx_2)):
#        zp = listy_2[ind]
        x =  listx_2[ind]
#        x = 43 #meV
        list_y_re.append(function_real_ana(x,n))
        
    maxi, _ = find_peaks(list_y_re, height=0)
    if len(maxi)>1:
        list_maximos = []
        for maxis in maxi: 
            list_maximos.append(list_y_re[maxis])
        index_max = np.argmax(list_maximos)
        maxis_zp_lambda_p.append(listx_2[maxi[int(index_max)]])
        values_maxis.append(np.max(list_maximos))
        values_maxis_tot.append(np.sum(list_maximos))
    else:
#        print(int(maxi))
        maxis_zp_lambda_p.append(listx_2[int(maxi)])
        values_maxis.append(list_y_re[int(maxi)])
        values_maxis_tot.append(list_y_re[int(maxi)])
    

    list_y_re_tot.append(list_y_re)    
    
#    print(n,maxis_zp_lambda_p,values_maxis)
#    list_y_re = np.array(list_y_re)/maxi
   
tot = np.sum(values_maxis)

print('a', a_nm, 'nm', 'zp:',zp_nano,'nm', 'maximo:', np.abs(tot))
#    
#%%

listx_3 = []
for ind in range(len(listy_2)):
    listx_3.append(listy_2[ind]/listz_2[ind])
    
listx_4 = np.linspace(np.min(listx_3),np.max(listx_3),len(list_y_re_tot[0]))

graph(title,labelx,labely ,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
k = 0
for n in list_n:   
    list_y_re = np.array(list_y_re_tot[k]) 
#    list_y_re = savgol_filter(list_y_re, 49, 3)  
    plt.plot(listx_4,list_y_re,'-',ms = ms, label = 'n = %i'%(n))
    k = k + 1
    
plt.legend(loc = 'best',markerscale=mk,fontsize=tamlegend,frameon=False,handletextpad=hp, handlelength=1)
#    plt.grid(1)
plt.tight_layout()

#plt.yscale('log')
os.chdir(path_save)
plt.savefig('decay_rate_fix_zp_' + labelp + '.png', format='png',bbox_inches='tight',pad_inches = 0.008,dpi = dpi)


#%%

