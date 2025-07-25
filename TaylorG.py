# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 09:26:04 2025

Class for the resolution of the Taylor-Goldstein equation for a squared 
periodic lattice with an incident wave as given by :
    Ghaemsaidi et al. JFM 789 - 2016

@author: paulb (Brancher Paul)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.interpolate import interp1d


class TaylorGoldstein(object):

    def __init__(self, y_lat, N_lat):
        self.z_tot, self.N_tot = y_lat, N_lat
        self.phi_0 = 1
        self.Mpts = None
        self.phi, self.phi_1, self.phi_2, self.phi_3 = None, None, None, None
        self.module = None
        self.lambda_z = None
        self.N_i, self.N_o = None, None
        self.L=None


    def __str__(self):
        """ object description """
        m0='-'*100
        m1 = "\nCreation of an periodic squared lattice of density"
        m2="\nResolution of the TaylorGoldstein equation on it"
        m3="for an incident wave as given by :\n"
        m4="\n"
        message=m0+m1+m2+m3+m4+m0
        return message


    def resolution(self, k, omega):
        """ From now, TaylorGoldstein.resolution(k, omega) always return y and phi """

        self.L = self.z_tot[-1]
        
        self.N_i = self.N_tot[0]
        self.N_o = self.N_tot[-1]

        self.Mpts = len(self.z_tot)

        self.solve_wave(k, omega)
        

        return self.z_tot, self.phi

   
    def N_z(self, z) :
        """ for a given z, retrun the value N(z) """
        fct = interp1d(self.z_tot, self.N_tot)
        # print(fct(z))
        return fct(z)
    
    
    @staticmethod
    def WaveNB(k, omega, nu, N1, N2):
        """ set the vertical wavenumbers, m1 and n1 """
        gamma1 = (4*k**2*N1**2*nu)/(omega**3)
        gamma2 = (4*k**2*N2**2*nu)/(omega**3)
        
        Aneg1  = k**2 - (omega*(1+gamma1**2)**0.25*np.sin(np.arctan2(gamma1,1)/2))/(2*nu)
        Aneg2  = k**2 - (omega*(1+gamma2**2)**0.25*np.sin(np.arctan2(gamma2,1)/2))/(2*nu)
        
        Apos1  = k**2 + (omega*(1+gamma1**2)**0.25*np.sin(np.arctan2(gamma1,1)/2))/(2*nu)
        Apos2  = k**2 + (omega*(1+gamma2**2)**0.25*np.sin(np.arctan2(gamma2,1)/2))/(2*nu)
        
            
        Bpos1 = (-omega + (omega*(1+gamma1**2)**0.25*np.cos(np.arctan2(gamma1,1)/2)))/(2*nu)
        Bpos2 = (-omega + (omega*(1+gamma2**2)**0.25*np.cos(np.arctan2(gamma2,1)/2)))/(2*nu)
        
        Bneg1  = (-omega - (omega*(1+gamma1**2)**0.25*np.cos(np.arctan2(gamma1,1)/2)))/(2*nu)
        Bneg2  = (-omega - (omega*(1+gamma2**2)**0.25*np.cos(np.arctan2(gamma2,1)/2)))/(2*nu)
            
        m1 = (Aneg1**2+Bpos1**2)**0.25*(np.cos(np.arctan2(Bpos1,Aneg1)/2) + 1j*np.sin(np.arctan2(Bpos1,Aneg1)/2))
        m2 = -m1
        m3 = (Apos1**2+Bneg1**2)**0.25*(np.cos(np.arctan2(Bneg1,Apos1)/2) + 1j*np.sin(np.arctan2(Bneg1,Apos1)/2))
        m4 = -m3
        
        n1 = (Aneg2**2+Bpos2**2)**0.25*(np.cos(np.arctan2(Bpos2,Aneg2)/2) + 1j*np.sin(np.arctan2(Bpos2,Aneg2)/2))
        n2 = -n1
        n3 = (Apos2**2+Bneg2**2)**0.25*(np.cos(np.arctan2(Bneg2,Apos2)/2) + 1j*np.sin(np.arctan2(Bneg2,Apos2)/2))
        n4 = -n3
        
        return -m1, m2, m3, m4, -n1, n2, n3, n4
    
    
    
    def solve_wave(self, k, omega):
        """ solve the Taylor-Goldstein equation for given omega and k_x """
        nu=1e-6

        m1, m2, m3, m4, n1, n2, n3, n4 = TaylorGoldstein.WaveNB(k, omega, nu, self.N_i, self.N_o)
        z_init=np.linspace(0, self.L, self.Mpts)
        phi_guess = np.zeros((8, z_init.size))

        self.lambda_z = 2*np.pi / np.imag(m1)

        def fun(z, y):
            A = (1j * omega / nu) - 2 * k**2
            B = k**2 * (k**2 + (1j * omega/nu) * (self.N_z(z)**2 / (omega**2) - 1))

            a=np.real(A)
            b=np.imag(A)
            c=np.real(B)
            d=np.imag(B)
            
            dy0 = y[1]
            dy1 = y[2]
            dy2 = y[3]
            dy3 = -a * y[2] + b * y[6] - c * y[0] + d * y[4]
            
            dy4 = y[5]
            dy5 = y[6]
            dy6 = y[7]
            dy7 = -a * y[6] - b * y[2] - c * y[4] - d * y[0]
            
            equation = np.vstack((dy0, dy1, dy2, dy3, dy4, dy5, dy6, dy7))
            return equation
        
        def bc(ya, yb):
            alpha = np.real(m1**2)
            beta = np.imag(m1**2)
            gamma = np.real(n1)
            delta = np.imag(n1)
            epsilon = np.real(n1**2)
            zeta = np.imag(n1**2)
            
            c1 = ya[0]-self.phi_0
            c2 = ya[4]
            c3 = ya[2]- alpha * ya[0] + (beta * ya[4])
            c4 = ya[6] - beta * ya[0] - (alpha * ya[4])

            c5 = yb[1]- gamma * yb[0] + delta * yb[4] 
            c6 = yb[5]- delta * yb[0] - gamma * yb[4]
            c7 = yb[2] - epsilon * yb[0] + zeta * yb[4]
            c8 = yb[6] - zeta * yb[0] - epsilon * yb[4]
    
            bc = np.array([c1, c3, c2, c4, c5, c6, c7, c8])

            return bc
        
        sol = solve_bvp(fun, bc, z_init, phi_guess)
        r_phi=sol.sol(self.z_tot)[0]            # Re phi
        i_phi = sol.sol(self.z_tot)[4]          # Im phi
        
        self.module = np.sqrt(r_phi**2 + i_phi**2)
        self.phi = r_phi + 1j * i_phi




    def get_results(self, indeD, indeA, indeB, indeC) :
        """ 
        Return the transmision with 3 different methods
        
        indeA is the index of the start of L_bot
        
        """

        phi_edge_tot = np.mean(abs(self.phi[indeA:]))
        T_tot = phi_edge_tot**2
        
        phi_edge_2 = np.mean(abs(self.phi[indeA:indeD]))
        T_2 = phi_edge_2**2
        
        phi_edge_1 = np.mean(abs(self.phi[indeC:indeB]))
        T_1 = phi_edge_1**2
        

        return T_1, T_2, T_tot
        
