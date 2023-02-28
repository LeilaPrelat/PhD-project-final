#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
"""
import numpy as np
import sys
import os 
from scipy import special

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_constants =  path_basic.replace('/1_dipole','')
#print('Importar modulos necesarios para este codigo')


try:
    sys.path.insert(1, path_constants)
    from dipole_moment_resonance import dipole_moment_anav2_for_decay_rate_resonance, dipole_moment_num_for_decay_rate_resonance, dipole_moment_pole_aprox_for_decay_rate_resonance
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)



try:
    sys.path.insert(1, path_constants)
    from green_self_image import green_self_ana_exponential_function, green_self_num_integral_inside_light_cone
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from green_self_image import green_self_num, green_self_pole_aprox 
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)

try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%

def EELS_film_ana_f_div_gamma0(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     ## normalizando con el paper 149
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """


    px_v,py_v,pz_v = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp) # multiplicar por e/(2*pi*v)
    

#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_exponential_function(omegac,epsi1,epsi3,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi3,d_nano,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
   

    
    factor_K = K0**2 + K1**2  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2

    
#    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
#    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

    k_prima = omegac*np.sqrt(epsi1)
    
    factor_final = k_prima/(12*np.pi*(int_v**2))

    rta = np.imag(Green_self)*factor_final/factor_K    
    

    return 3*rta 

#%%
    
def EELS_film_num_f_div_gamma0(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     ## normalizando con el paper 149
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """

    px_v,py_v,pz_v = dipole_moment_num_for_decay_rate_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp) # multiplicar por e/(2*pi*v)
    

#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_num(omegac,epsi1,epsi3,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi3,d_nano,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
   

    
    factor_K = K0**2 + K1**2  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2

    
#    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
#    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

    k_prima = omegac*np.sqrt(epsi1)
    
    factor_final = k_prima/(12*np.pi*(int_v**2))

    rta = np.imag(Green_self)*factor_final/factor_K    
    

    return 3*rta 


def EELS_film_pole_aprox_f_div_gamma0(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     ## normalizando con el paper 149
    """    
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    z : coordenada z
    xD : coordenada x del dipolo 
    yD : coordenada y del dipolo
    zD : coordenada z del dipolo 
    zp : posicion del plano (>0)
    px : coordenada x del dipolo 
    py : coordenada y del dipolo
    pz : coordenada z del dipolo
    Returns
    -------
    formula del potencial electric con QE approximation, rp con 
    aproximacion del polo y con aprox de principal value para las integrales
    con rp
    """
    px_v,py_v,pz_v = dipole_moment_pole_aprox_for_decay_rate_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp) # multiplicar por e/(2*pi*v)
    

#    print(px_v,py_v,pz_v)
#    px_tot_2 = np.abs(px_v)**2 + np.abs(py_v)**2 + np.abs(pz_v)**2 
    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi3,d_nano,zp)
    
    rtaself_x, rtaself_y, rtaself_z  = rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2
    
    Green_self = rtaself_x*(np.abs(px_v)**2) + rtaself_y*(np.abs(py_v)**2)  + rtaself_z*(np.abs(pz_v)**2)

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
   

    
    factor_K = K0**2 + K1**2  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2

    
#    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
#    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

    k_prima = omegac*np.sqrt(epsi1)
    
    factor_final = k_prima/(12*np.pi*(int_v**2))

    rta = np.imag(Green_self)*factor_final/factor_K    
    

    return 3*rta 
