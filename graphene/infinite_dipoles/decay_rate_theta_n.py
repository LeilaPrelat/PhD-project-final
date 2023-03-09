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
    from dipole_moment_resonance import dipole_moment_anav2_for_decay_rate_resonance, dipole_moment_num_for_decay_rate_resonance, dipole_moment_pole_aprox_for_decay_rate_resonance
except ModuleNotFoundError:
    print('dipole_moment.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_constants)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_constants)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%% 

# normalizado con el paper 149 
def decay_rate_theta_inf_dipoles_ana_res_div_gamma0(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a
    
#    print(delta_n/kx)

    den = np.sqrt(kp**2 - kx**2)
    
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    term_kp_2 = kp_2 + kp
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = -np.exp(-2*kp_2*zp)*Rp*kp*(px*kx*term_kp/den + py*term_kp + 1j*pz*term_kp_2/den )/(4*np.pi*a)
    
#    cte_formula = a/(48*(np.pi**2)*Rp)
#    cte_formula = 12*np.pi*a/(Rp) ## hay un extra 1/(2pi) en la formula de phi. necesario para silver 
#      
##    cte_formula = a/(12*Rp) ## hay un extra 1/(2pi) en la formula de phi
#    
##    cte_formula = a*np.pi/Rp  ## hay un extra 1/(2pi) en la formula de phi. necesario para grafeno  
##
#    
#    factor_K = K0**2 + K1**2  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
#  
##    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
##    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2
#
##    cte_extra = (a*omegac)**4
#    extra_cte_adimensional = a*omegac ## para poder comparar con diferentes materiales y que no dependa del periodo "a" 

#    rta = 2*epsi1*(np.abs(phi_n)**2)*cte_formula*k_prima*(int_v**(-2))/(factor_K*seno_theta_n)

    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    seno_theta_n = den/kp     
    k_prima = omegac*np.sqrt(epsi1)    
    v = c/int_v
    
    imag_alpha = 3*epsi1/(2*k_prima**3) ## = np.imag(-3*epsi1/(2*1j*k_prima**3))

    factor_gamma0 = (2*omegac*int_v/v)**2
    Gamma_EELS = factor_gamma0*(K0**2 + K1**2)*imag_alpha/np.pi  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
       
    Gamma_SPn = a*epsi1*np.abs(phi_n)**2/(np.pi*np.abs(Rp)*seno_theta_n*4*np.pi**2*v**2)
    
    return Gamma_SPn*5*1e4/Gamma_EELS



def decay_rate_theta_inf_dipoles_ana_res_div_gamma0_v2(omegac,epsi1,epsi2,hbmu,hbgama,int_v,zp,a,b,n):     
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

#    x, y, z = 0,0,0
    E = omegac*aux
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
   
    cond = 4*np.pi*alfac*sigma_DL(E,hbmu,hbgama)
    Rp = 2*epsi1/(epsi1 + epsi2)
    alfa_p = 1j*(epsi1 + epsi2)/(cond*cte1)
    kp = alfa_p*k1
    
    
    px,py,pz  = dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi2,hbmu,hbgama,int_v,b,zp)  
#    list_dipoles = np.linspace(-Nmax,Nmax,2*Nmax + 1)
#            
    kx = omegac*int_v + 2*np.pi*n/a
    
    delta_n = 2*np.pi/a
#    print(delta_n/kx)

    den = np.sqrt(kp**2 - kx**2)
    
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    term_kp_2 = kp_2 + kp
   # return np.imag(final_2*cte*kp*np.cos(theta))
    phi_n = -np.exp(-2*kp_2*zp)*Rp*kp*(px*kx*term_kp/den + py*term_kp + 1j*pz*term_kp_2/den )/(4*np.pi*a)
    
    cte_formula = 4*np.pi*a/Rp
#    cte_formula = 12*np.pi*a/(Rp) ## hay un extra 1/(2pi) en la formula de phi. necesario para silver 
#      
#    cte_formula = a/(12*Rp) ## hay un extra 1/(2pi) en la formula de phi
    
#    cte_formula = a*np.pi/Rp  ## hay un extra 1/(2pi) en la formula de phi. necesario para grafeno  
#

    gamma = (1 - (int_v**(-2)) )**(-1/2)
    arg = np.abs(b)*omegac*int_v/gamma
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)

    
    factor_K = ( (K0/gamma)**2 + K1**2)/(gamma**2)  ## decay rate de 1 dipolo # pero sin el "e/hbar" se cancela con el momento dipolar^2
    seno_theta_n = den/kp
#    print(seno_theta_n)
    
#    px_dir,py_dir,pz_dir = dipole_moment_anav2_for_decay_rate_resonance_dir(omegac,int_v,b,zp)        
#    denominador = np.abs(px_dir)**2 +  np.abs(py_dir)**2 +  np.abs(pz_dir)**2

#    cte_extra = (a*omegac)**4
    extra_cte_adimensional = a*omegac ## para poder comparar con diferentes materiales y que no dependa del periodo "a" 
   
    k_prima = omegac*np.sqrt(epsi1)
        
    rta = 2*epsi1*(np.abs(phi_n)**2)*cte_formula*k_prima*(int_v**(-2))/(factor_K*seno_theta_n)
        
    return 12*rta

