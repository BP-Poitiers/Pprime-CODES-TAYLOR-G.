# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 14:52:37 2025

@author: Paul Brancher

Tutorial for a resolution from the TaylorGoldstein class, show for a given
1/3 ; 1/3 ; 1/3, ratio lattice with 4 steps

"""

import numpy as np
from Reader import XDMF_reader
from Affichage import affichage
import matplotlib.pyplot as plt
from CMAP import cmpw
from TaylorG import TaylorGoldstein


# =============================================================================
# Lattice for the example
# =============================================================================

dy=0.007
k=27.6
Ly = 2*np.pi / k * 21 / 2

def lat_cte(taille, dy, N):
    pts = taille / dy + 1
    lat = np.linspace(N, N, int(pts))

    return lat


def MARCHE_1():
    haut = lat_cte(1/3*0.1, dy, 1)
    bas = lat_cte(1/3*0.1, dy, 0)
    
    marche_1 = np.concatenate((haut, bas, haut))
    return marche_1


def RES_1():
    
    L_HAUT = 0.4
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_1(),MARCHE_1(),MARCHE_1(),MARCHE_1(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT


# =============================================================================
# Creation of y and N(y), and other parameters
# =============================================================================


y_lat, N_lat, lbot = RES_1()

omega=0.4


# =============================================================================
# Use of the Taylor-Goldstein class
# =============================================================================

simu = TaylorGoldstein(y_lat, N_lat)
y, phi = simu.resolution(k, omega)


# exemple de plot
affichage.simple_plot(np.real(phi)[::-1], y, 'red', [r'$\phi$', 'y'], r'$\phi(y)$', shape=(8,10))



# =============================================================================
# Other usefull functions
# =============================================================================

def plot_res(tup, percent=0.75):
    """
    tup : tuple of 3 elements - y, N(y) and Lbot
    Plot a lattice for a window (percent*Ly, Ly)
    """
    y, res, L_BOT = tup
    affichage.simple_plot(res[::-1], y, 'red', ['N', None], 'lattice', shape=(7, 10), nbplot=2, lw=[3.5])
    plt.ylim(percent*Ly, Ly)
    plt.axhspan(0, L_BOT, facecolor='none', alpha=1, edgecolor='lightseagreen',lw=0, hatch='//', label="Mean Zone")
    plt.legend(fontsize=20, loc='lower left')
    


def get_transmi(tup, npt):
    """
    Get T=f(omega)
    Parameters
    ----------
    tup : tuple of 3 elements - y, N(y) and Lbot
    npt : number of points on the transmission function
    Returns
    -------
    omegaL : 1D array of omega
    np.array(T) : 1D array of the corresponding transmission
    """
    y, res, L_BOT = tup
    omegaL = np.linspace(0.1, 1, npt)
    T = []

    L_A = Ly-L_BOT
    indeA = np.abs(y - L_A).argmin()

    for i in range(len(omegaL)) :
        if i%100 == 0 :
            print('ite :', i, 'sur', len(omegaL))
        simu = TaylorGoldstein(y, res)
        useless1, useless2 = simu.resolution(k, omegaL[i])
        T_1, T_2, T_tot = simu.get_results(5, indeA, 5, 5)
        T.append(T_tot)

    return omegaL, np.array(T)


def plot_transmi(tup, npt, color='purple'):
    """
    Plot the transmission function for a given lattice OR a given T=f(omega)
    Parameters
    ----------
    tup : tuple of 3 elements - y, N(y) and Lbot OR omega and T(omega)
    npt : number of points on the transmission function
    color : str, optional - The default is 'purple'. (color of the plot)
    """
    if len(tup)==3 :
        o4, t4 = get_transmi(tup, npt)
    else :
        o4, t4 = tup
    affichage.simple_plot(o4, t4, color, [r'$\omega$', 'T'], None, nbplot=3, cs='%.0f', csx='%.1f')
    plt.ylim(-0.05, 3.5)
    plt.xlim(0.1, 1)
    plt.tight_layout()



def get_SurfState(tup, npt):
    """
    Get all |phi|=f(omega, y)
    Parameters
    ----------
    tup : tuple of 3 elements - y, N(y) and Lbot OR omega and T(omega)
    npt : number of points on the transmission function
    Returns
    -------
    All parameters used in plot_SURF_STATE
    """
    y_lat, N_lat, LB = tup

    omegaL = np.linspace(0.1, 1, npt)
    L_phi = []

    for i in range(len(omegaL)) :
        if i%50 == 0 :
            print('ite :', i, 'sur', len(omegaL))
        simu = TaylorGoldstein(y_lat, N_lat)
        zzz, phi = simu.resolution(k, omegaL[i], retuP=True)
        L_phi.append(np.real(phi))

    OMEG, ZZZ = np.meshgrid(np.array(omegaL), zzz)

    L_phi = np.array(L_phi).T
    L_phi = abs(L_phi)
    
    return OMEG, ZZZ, L_phi, N_lat, zzz



def plot_SURF_STATE(tup, npt, maxx=5, percent=0.75):
    """ Plot |phi|=f(omega, y) from a tup-lattice OR from get_SurfState """
    if len(tup)==5 :
        OMEG, ZZZ, L_phi, N_lat, zzz = tup
    else :
        OMEG, ZZZ, L_phi, N_lat, zzz = get_SurfState(tup, npt)

    ax = [r'$\omega$', 'y']
    cs='%.1f'
    cb = [r'|$\Phi$|', 0, maxx]
    cp='inferno_r'

    affichage.graph2D(OMEG, ZZZ[::-1], L_phi, 0, cbar=cb, axes=ax, nbplot=4, cs=cs, cmap=cp, shape=(17, 10))
    plt.xlim(OMEG.min(), OMEG.max())
    plt.ylim(percent*ZZZ.max(), ZZZ.max())
    plt.plot(N_lat[::-1]*100-50, zzz, '--', color='k')
    plt.tight_layout()

