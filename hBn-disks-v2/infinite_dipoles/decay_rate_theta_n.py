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
path_constants =  path_basic.replace('/infinite_dipoles','')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_constants)
    from hBn_PP import epsilon_x,hBn_lambda_p,hBn_Rp
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_constants)


try:
    sys.path.insert(1, path_basic)
    from dipole_moment import dipole_moment_ana_for_decay_rate,dipole_moment_pole_aprox_for_decay_rate,dipole_moment_num_for_decay_rate
except ModuleNotFoundError:
    print('potential.py no se encuentra en ' + path_basic)


#
#try:
#    sys.path.insert(1, path_basic)
#    from potential_induced import potential_ana_resonance_v1,potential_num_resonance
#except ModuleNotFoundError:
#    print('potential_induced.py no se encuentra en ' + path_basic)
    
try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def decay_rate_theta_inf_dipoles_ana_div_gamma0(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,zp,a,b,n):     
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


    E = omegac*aux

    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x
    
    d_micro = d_nano_film*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(1 - epsi_HBN_par))
    kp = alfa_p*omegac

    Rp = 1

    px,py,pz  = dipole_moment_ana_for_decay_rate(omegac,epsi_silica,d_nano_film,d_thickness_disk_nano,D_disk_nano,int_v,b,zp)  

    kx = omegac*int_v + 2*np.pi*n/a     

    den = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    term_kp_2 = kp_2 + kp
    phi_n = -np.exp(-2*kp_2*zp)*Rp*kp*(px*kx*term_kp/den + py*term_kp + 1j*pz*term_kp_2/den )/(4*np.pi*a)
    

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    seno_theta_n = den/kp     
    k_prima = omegac*np.sqrt(epsi_silica(E))    
    v = c/int_v
    
    imag_alpha = 3*epsi_silica(E)/(2*k_prima**3) ## = np.imag(-3*epsi1/(2*1j*k_prima**3))

    factor_gamma0 = (2*omegac*int_v/v)**2
    Gamma_EELS = factor_gamma0*(K0**2 + K1**2)*imag_alpha/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
       
    Gamma_SPn = a*epsi_silica(E)*np.abs(phi_n)**2/(np.pi*np.abs(Rp)*seno_theta_n*4*np.pi**2*v**2)
    
    return Gamma_SPn*1e9*1e3/Gamma_EELS



