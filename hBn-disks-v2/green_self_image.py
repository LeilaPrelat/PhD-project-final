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
    from rp_coefficient import epsilon_x, epsilon_z
except ModuleNotFoundError:
    print('hBn_PP.py no se encuentra en ' + path_basic)


try:
    sys.path.insert(1, path_basic)
    from constants import constantes
except ModuleNotFoundError:
    print('constants.py no se encuentra en ' + path_basic)

pi,hb,c,alfac,mu1,mu2 = constantes()
aux = c*hb

#%%


def green_self_num(omegac,epsi_silica,d_nano,zp_micro):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux  
    k1 = omegac
    k1_3 = k1**3


    epsi_x = epsilon_x(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x

    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))


    z_dip_barra_self = k1*2*zp_micro  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(1 - epsi_HBN_par))
    
    
#    d_micros = d_nano*1e-3
#    lambda_p_v = hBn_lambda_p(E,epsi_silica(E),epsi_silica(E))*d_micros
#    kp = 2*np.pi/lambda_p_v  
#    alfa_p = kp/omegac



    rp = lambda u: u/(u*(1-r*expB_self(u)) - alfa_p)


    cota_sup1 = 400/omegac
    cota_inf1 = 0.01/omegac
    
####       
    cte_x = k1_3*0.5/epsi_silica(E) # antes epsi1 era 1

    
  


    IntselfB_function_re_xx = lambda u: np.real((u**2)*rp(u)*expB_self(u))
    IntselfB_function_im_xx = lambda u: np.imag((u**2)*rp(u)*expB_self(u))

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_xx, cota_inf1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_xx, cota_inf1, cota_sup1)
#    
#    print(intselfB_re_x)
    
    
    rtaself_x = (intselfB_re_x + 1j*intselfB_im_x)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z




#%%


def green_self_num_integral_inside_light_cone(omegac,epsi_silica,d_nano,zp_micro):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux  
    k1 = omegac
    k1_3 = k1**3


    epsi_x = epsilon_x(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x

    r = (1 - epsi_silica(E))/(1 + epsi_silica(E))


    z_dip_barra_self = k1*2*zp_micro  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 
    
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(1 - epsi_HBN_par))
    



    rp = lambda u: u/(u*(1-r*expB_self(u)) - alfa_p)


    cota_sup1 = np.sqrt(epsi_silica(E))
    cota_inf1 = 0.01
    
####       
    cte_x = k1_3*0.5/epsi_silica(E) # antes epsi1 era 1

    



    IntselfB_function_re_xx = lambda u: np.real((u**2)*rp(u)*expB_self(u))
    IntselfB_function_im_xx = lambda u: np.imag((u**2)*rp(u)*expB_self(u))

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_xx, cota_inf1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_xx, cota_inf1, cota_sup1)
#    
#    print(intselfB_re_x)
    
    
    rtaself_x = (intselfB_re_x + 1j*intselfB_im_x)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z



#%%
    #     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
def green_self_pole_aprox(omegac,epsi_silica,d_nano,zp_micro):
    """
    Parameters
    ----------
    omegac : omega/c = k0 en 1/micrometros    
    epsi1 : epsilon del medio de arriba del plano
    epsi2 : epsilon del medio de abajo del plano
    hbmu : chemical potential in eV  
    hbgama : collision frequency in eV
    int_v : fraccion de la vel de la c que vale la velocidad uniforme del electron en x, si v = c/int
    b: electron position en z, b < 0 (dipolo en 0 y eje de z apunta hacia abajo) 
    zp : coordenada zp del plano, zp > 0  (dipolo en 0 y eje de z apunta hacia abajo) 
    x : punto x donde miramos el campo en micrones
    omega0 : resonant frequency 
    kappa_factor_omega0 : kappa (collision frequency) = kappa_factor*omega0
    kappa_r_factor : kappa_r = kappa_r_factor*kappa with kappa_r_factor<1    
    Returns
    -------
    alfa effectivo en QE approx
    """

    E = omegac*aux  
    k1 = omegac
    k1_3 = k1**3

    epsi_x = epsilon_x(E)
    
#    if np.imag(np.sqrt(epsi_x*epsi_z))>0:
#        epsi_m = np.sqrt(epsi_x*epsi_z)
#    else:
#        epsi_m = -np.sqrt(epsi_x*epsi_z)
#    

    epsi_HBN_par = epsi_x

#    sigma_2D = d_nano*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(1 - epsi_HBN_par ))
    
    rp = lambda u: u/(u - alfa_p)   
#    print(np.real(rp(10)))
    
#    rs_num = lambda u: 1j*u - 1j*u - cond/cte1
#    rs_den = lambda u: 1j*u + 1j*u + cond/cte1
#    rs = lambda u: rs_num(u)/rs_den(u)
    cota_sup1 = 400/omegac
    cota_inf1 = 0.01/omegac
    
    
####       
    cte_x = k1_3*0.5/epsi_silica(E) # antes epsi1 era 1

    

    
    z_dip_barra_self = k1*2*zp_micro  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = lambda u: np.exp(-u*z_dip_barra_self) 

    IntselfB_function_re_xx = lambda u: u**2*np.real(rp(u))*expB_self(u)
    IntselfB_function_im_xx = lambda u: u**2*np.imag(rp(u))*expB_self(u)

    intselfB_re_x,err = integrate.quad(IntselfB_function_re_xx, cota_inf1, cota_sup1)
    intselfB_im_x,err = integrate.quad(IntselfB_function_im_xx, cota_inf1, cota_sup1)
    
#    print(intselfB_re_x)
    
    
    rtaself_x = (intselfB_re_x + 1j*intselfB_im_x)*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z


#%%
    

def green_self_ana_residuos(omegac,epsi_silica,d_nano,zp_micro):     
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


    Rp = 1
    epsi_x = epsilon_x(E)
    epsi_HBN_par = epsi_x

#    sigma_2D = d_nano*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(1 - epsi_HBN_par))
    kp = alfa_p*omegac

    kp_3 = kp**3
    
#    kp = np.real(kp)   ### neglecting the imaginary part of kp para la formula de EELS (la del campo)
#    alfa_p = np.real(alfa_p)
      

    cte_x = 1j*np.pi*Rp*kp_3*0.5/epsi_silica(E) # antes epsi1 era 1

    z_dip_barra_self = omegac*2*zp_micro  #NO ESTOY SEGURA DE ESTA EXPONENCIAL EH 
    expB_self = np.exp(-alfa_p*z_dip_barra_self) 

    
    rtaself_x = expB_self*cte_x
    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    

    return rtaself_x, rtaself_y, rtaself_z

#     rp = lambda u: alpha_parallel(u)/(alpha_parallel(u) - alfa_p)
def green_self_ana_exponential_function(omegac,epsi_silica,d_nano,zp_micro):     
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

#    sigma_2D = d_nano*(epsi_HBN_par - 1 )/(4*pi)  ## se cancela el i*omega del sigma 
    d_micro = d_nano*1e-3
    alfa_p = epsi_silica(E)*2/(omegac*d_micro*(1 - epsi_HBN_par))
    kp = alfa_p*omegac



    arg = -2*kp*zp_micro  
    
    try:

        dif_term = kp*special.exp1(arg)*np.exp(arg) + (2*zp_micro)**(-1)
    except RuntimeWarning:
        dif_term = kp*np.pi*1j*np.exp(arg)
        
    
    rtaself_x = 0.5*( 2*(2*zp_micro)**(-3) + kp*(2*zp_micro)**(-2) + kp**2*dif_term   )/epsi_silica(E)
    
    

    rtaself_y = rtaself_x
    rtaself_z = 2*rtaself_x
    
    return rtaself_x, rtaself_y, rtaself_z



