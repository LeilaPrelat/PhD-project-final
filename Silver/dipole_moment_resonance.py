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
from scipy import integrate

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Silver_PP import Silver_lambda_p, Silver_Rp, epsilon_m
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_ana_exponential_function, green_self_num_integral_inside_light_cone
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from green_self_image import green_self_num, green_self_pole_aprox 
except ModuleNotFoundError:
    print('green_self_image.py no se encuentra en ' + path_basic)




try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb


#%%

def dipole_moment_anav2_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    rtaself_x, rtaself_y, rtaself_z  =  green_self_ana_exponential_function(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
    d_micros = d_nano*1e-3
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    px = alffa_eff_x*1j*omegac*int_v*(2*K0 - 2*np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*( - np.sqrt(2)*omegac*int_v*K1 + 2*np.pi*Rp*kp*expo)
    
    pz = alffa_eff_z*( - np.sqrt(2)*omegac*int_v*K1 + 2*np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz

def dipole_moment_anav2_for_decay_rate_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    px,py,pz en unidades de k*alfa_eff
    """

    E = omegac*aux
#    k0 = omegac #=omega/c
    # x_y = ky/k0
    n1 = epsi1*mu1
    cte1 = np.sqrt(n1)
    k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1

    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_ana_exponential_function(omegac,epsi1,epsi3,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi3,d_nano,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2

    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)


#    charge_electron = 4.806e-10/c
#    cte_uni = int_v/(2*np.pi*c)
    
    d_micros = d_nano*1e-3
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

      
    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)
    
    kx = omegac*int_v
#    expo = np.exp(-np.sqrt(kx**2 + kp**2)*(np.abs(b) + 2*zp))
    expo = np.exp(-kp*(np.abs(b) + 2*zp))
    
#    den = np.sqrt(kx**2 + kp**2) - kp
    ky = np.sqrt(kp**2 - kx**2)
    kp_2 = np.sqrt(kp**2)
    term_kp = 1 + kp/kp_2
    
#    term_extra = 2*np.pi*1j*Rp*kp*np.abs(kp)*expo/ky
    
    
    px = alffa_eff_x*1j*omegac*int_v*(2*K0 - 2*np.pi*1j*Rp*kp*expo/ky)
    
    py = alffa_eff_y*( - np.sqrt(2)*omegac*int_v*K1 + 2*np.pi*Rp*kp*expo)
    
    pz = alffa_eff_z*( - np.sqrt(2)*omegac*int_v*K1 + 2*np.pi*1j*Rp*(kp**2)*expo/ky )
    
    return px, py, pz


#%%



def dipole_moment_num_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    

    rtaself_x, rtaself_y, rtaself_z  =  green_self_num(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)

#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)
    
    
#    rp_num = lambda u: epsi2*1j*kparallel(u) - epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp_den = lambda u: epsi2*1j*kparallel(u) + epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp = lambda u: rp_num(u)/rp_den(u)

    d_micros = d_nano*1e-3
    epsi_2 = epsilon_m(E)

    kz1 =  lambda u: np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 =  lambda u: np.sqrt(epsi_2 - u**2)
    kz3 =  lambda u: np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)



#    kz1_v2 = lambda u : kz1(u) if np.imag(kz1(u)) > 0 else -kz1(u) ### convergencia 
#    kz2_v2 = lambda u : kz2(u) if np.imag(kz2(u)) > 0 else -kz2(u)
#    kz3_v2 = lambda u : kz3(u) if np.imag(kz3(u)) > 0 else -kz3(u)


    r12 =  lambda u: (kz1(u)*epsi_2 - kz2(u)*epsi1)/(kz1(u)*epsi_2 + kz2(u)*epsi1)
    r21 = lambda u: (kz2(u)*epsi1 - kz1(u)*epsi_2)/(kz2(u)*epsi1 + kz1(u)*epsi_2)
    r23 =  lambda u: (kz2(u)*epsi3 - kz3(u)*epsi_2)/(kz2(u)*epsi3 + kz3(u)*epsi_2)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_2)
    t12 = lambda u: 2*kz1(u)*cte_t/(kz1(u)*epsi_2 + kz2(u)*epsi1)
    t21 = lambda u: 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_2)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1 - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)  
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(alpha_parallel(u))*expo(u))
    int_f_im_z = lambda u: np.imag(rp(alpha_parallel(u))*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(2*K0 - INT_x)
    
    py = alffa_eff_y*(  - np.sqrt(2)*omegac*int_v*K1 - 1j*k0*INT_y)
    
    pz = alffa_eff_z*(  - np.sqrt(2)*omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz





def dipole_moment_num_for_decay_rate_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    

    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_num(omegac,epsi1,epsi3,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi3,d_nano,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2


    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)

#    charge_electron = alffa_eff*4.806e-10*int_v/(2*np.pi)
    charge_electron = 4.806e-10/c
    cte_uni = int_v/(2*np.pi*c)
    
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)
    
    
#    rp_num = lambda u: epsi2*1j*kparallel(u) - epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp_den = lambda u: epsi2*1j*kparallel(u) + epsi1*1j*kparallel(u) - cond*kparallel(u)**2/k0
#    rp = lambda u: rp_num(u)/rp_den(u)

    d_micros = d_nano*1e-3
    epsi_2 = epsilon_m(E)

    kz1 =  lambda u: np.sqrt(epsi1 - u**2) if (u**2 <= epsi1) else 1j*np.sqrt(u**2 - epsi1)
    kz2 =  lambda u: np.sqrt(epsi_2 - u**2)
    kz3 =  lambda u: np.sqrt(epsi3 - u**2) if (u**2 <= epsi3) else 1j*np.sqrt(u**2 - epsi3)



#    kz1_v2 = lambda u : kz1(u) if np.imag(kz1(u)) > 0 else -kz1(u) ### convergencia 
#    kz2_v2 = lambda u : kz2(u) if np.imag(kz2(u)) > 0 else -kz2(u)
#    kz3_v2 = lambda u : kz3(u) if np.imag(kz3(u)) > 0 else -kz3(u)


    r12 =  lambda u: (kz1(u)*epsi_2 - kz2(u)*epsi1)/(kz1(u)*epsi_2 + kz2(u)*epsi1)
    r21 = lambda u: (kz2(u)*epsi1 - kz1(u)*epsi_2)/(kz2(u)*epsi1 + kz1(u)*epsi_2)
    r23 =  lambda u: (kz2(u)*epsi3 - kz3(u)*epsi_2)/(kz2(u)*epsi3 + kz3(u)*epsi_2)
    

    exp_fresnel = lambda u: np.exp(1j*2*kz2(u)*omegac*d_micros)
    
    cte_t = np.sqrt(epsi1*epsi_2)
    t12 = lambda u: 2*kz1(u)*cte_t/(kz1(u)*epsi_2 + kz2(u)*epsi1)
    t21 = lambda u: 2*kz2(u)*cte_t/(kz2(u)*epsi1 + kz1(u)*epsi_2)

    rp_num = lambda u: t12(u)*t21(u)*r23(u)*exp_fresnel(u)
    rp_den = lambda u: 1 - r21(u)*r23(u)*exp_fresnel(u)
    rp = lambda u: r12(u) +  rp_num(u)/rp_den(u)  
 
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(alpha_parallel(u))*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(alpha_parallel(u))*expo(u))
    int_f_im_z = lambda u: np.imag(rp(alpha_parallel(u))*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(2*K0 - INT_x)
    
    py = alffa_eff_y*(  - np.sqrt(2)*omegac*int_v*K1 - 1j*k0*INT_y)
    
    pz = alffa_eff_z*(  - np.sqrt(2)*omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz






#%%



def dipole_moment_pole_aprox_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    d_micros = d_nano*1e-3
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

    rtaself_x, rtaself_y, rtaself_z  =  green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp)
    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
 
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: Rp*alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(2*K0 - INT_x)
    
    py = alffa_eff_y*(  - np.sqrt(2)*omegac*int_v*K1 - 1j*k0*INT_y)
    
    pz = alffa_eff_z*(  - np.sqrt(2)*omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz

#%%


def dipole_moment_pole_aprox_for_decay_rate_resonance(omegac,epsi1,epsi3,d_nano,int_v,b,zp):     
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
    k0 = omegac #=omega/c
    # x_y = ky/k0
#    n1 = epsi1*mu1
#    cte1 = np.sqrt(n1)
    # k1 = omegac*cte1
#    k1_2 = (k0*cte1)**2
 #   n_v1 = int_v/cte1
    
    d_micros = d_nano*1e-3
    Rp = Silver_Rp(E,epsi1,epsi3)
    lambda_p_v = Silver_lambda_p(E,epsi1,epsi3)*d_micros
    kp = 2*np.pi/lambda_p_v
    alfa_p = kp/omegac

    rtaself_x1, rtaself_y1, rtaself_z1  =  green_self_pole_aprox(omegac,epsi1,epsi3,d_nano,zp)
    rtaself_x2, rtaself_y2, rtaself_z2  =  green_self_num_integral_inside_light_cone(omegac,epsi1,epsi3,d_nano,zp)

    rtaself_x, rtaself_y, rtaself_z  =  rtaself_x1 - rtaself_x2, rtaself_y1 - rtaself_y2, rtaself_z1 - rtaself_z2


    alffa_eff_x = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_x))**(-1)
    alffa_eff_y = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_y))**(-1)
    alffa_eff_z = 1j*(2*omegac**3/(3*epsi1) +  np.imag(rtaself_z))**(-1)
 
    cota_inf = 0.01/k0
    cota_sup = 600/k0   
    
    alpha_x = int_v    
    
#    term6 = np.sign(z)*pz*Rp*kp_2*J0*exp_electron
    alpha_parallel = lambda u: np.sqrt(alpha_x**2 + u**2)

    rp = lambda u: Rp*alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
      
    expo = lambda u: np.exp(-np.sqrt(alpha_x**2 + u**2)*k0*(2*zp + np.abs(b)))
    
    int_f_re_x = lambda u: np.real(rp(u)*expo(u)/alpha_parallel(u))
    int_f_im_x = lambda u: np.imag(rp(u)*expo(u)/alpha_parallel(u))
    
    INT_re_x,err = integrate.quad(int_f_re_x, cota_inf, cota_sup) 
    INT_im_x,err = integrate.quad(int_f_im_x, cota_inf, cota_sup) 
    
    INT_x = INT_re_x + 1j*INT_im_x
    
    
    int_f_re_y = lambda u: np.real(rp(u)*expo(u)*u/alpha_parallel(u))
    int_f_im_y = lambda u: np.imag(rp(u)*expo(u)*u/alpha_parallel(u))
    
    INT_re_y,err = integrate.quad(int_f_re_y, cota_inf, cota_sup) 
    INT_im_y,err = integrate.quad(int_f_im_y, cota_inf, cota_sup) 
    
    INT_y = INT_re_y + 1j*INT_im_y
        


    int_f_re_z = lambda u: np.real(rp(u)*expo(u))
    int_f_im_z = lambda u: np.imag(rp(u)*expo(u))
    
    INT_re_z,err = integrate.quad(int_f_re_z, cota_inf, cota_sup) 
    INT_im_z,err = integrate.quad(int_f_im_z, cota_inf, cota_sup) 
    
    INT_z = INT_re_z + 1j*INT_im_z


    arg = np.abs(b)*omegac*int_v
    K1 = special.kn(1,arg)
    K0 = special.kn(0,arg)    
    
    
    px = alffa_eff_x*1j*omegac*int_v*(2*K0 - INT_x)
    
    py = alffa_eff_y*(  - np.sqrt(2)*omegac*int_v*K1 - 1j*k0*INT_y)
    
    pz = alffa_eff_z*(  - np.sqrt(2)*omegac*int_v*K1 + k0*INT_z )
    
    return px, py, pz












