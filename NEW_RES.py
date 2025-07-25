# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 09:06:09 2025

@author: paulb


# o8, t8 = get_transmi(RES_81(), 500)
# XDMF_reader.creaFichier('T_res81.csv', o8, t8)
# freq, trans = XDMF_reader.lecture('T_res81.csv')
# plot_transmi((freq, trans), 3, color='red')

# yyy, nnn, lb = RES_21()



# XDMF_reader.creaFichier('data_surf21.csv', yyy/0.003, nnn[::-1])


"""

import numpy as np
# from parameters import Simulation
from Reader import XDMF_reader
from Affichage import affichage
import matplotlib.pyplot as plt
from CMAP import cmpw
from TaylorG import TaylorGoldstein





k=27.6
Ly = 2*np.pi / k * 21 / 2
# dy = 0.0007
dy = 0.0007

def lat_cte(taille, dy, N):
    pts = taille / dy + 1
    lat = np.linspace(N, N, int(pts))

    return lat

def periode(taille, ratio, dy, N_1, N_2):
    pts = taille / dy + 1
    creux = np.linspace(N_2, N_2, int(pts*ratio))
    marche = np.linspace(N_1, N_1, int(pts*(1-ratio)))
    periode = np.concatenate((creux, marche))

    return periode

def res_I(n_marche, l_marche, ratio, dy, N_1, N_2):
    marches = []

    for i in range(n_marche) :
        perio = periode(l_marche[i], ratio[i], dy, N_1, N_2)
        marches.append(perio)
    
    return np.concatenate(tuple(marches))

def MARCHE_1():
    haut = lat_cte(1/3*0.1, dy, 1)
    bas = lat_cte(1/3*0.1, dy, 0)
    
    marche_1 = np.concatenate((haut, bas, haut))
    return marche_1

def MARCHE_2():
    haut = lat_cte(2/3*0.1, dy, 1)
    bas = lat_cte(1/6*0.1, dy, 0)
    
    
    marche_2 = np.concatenate((bas, haut, bas))

    return marche_2



def RES_1():
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.4
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_1(),MARCHE_1(),MARCHE_1(),MARCHE_1(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT

def RES_2():
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.4
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_2(),MARCHE_2(),MARCHE_2(),MARCHE_2(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT


def plot_res(tup, percent=0.75):
    y, res, L_BOT = tup
    affichage.simple_plot(res[::-1], y, 'red', ['N', None], 'lattice', shape=(7, 10), nbplot=2, lw=[3.5])
    plt.ylim(percent*Ly, Ly)
    plt.axhspan(0, L_BOT, facecolor='none', alpha=1, edgecolor='lightseagreen',lw=0, hatch='//', label="Mean Zone")
    plt.legend(fontsize=20, loc='lower left')
    




def get_transmi(tup, npt):
    y, res, L_BOT = tup

    omegaL = np.linspace(0.1, 1, npt)
    T = []

    L_A = Ly-L_BOT
    indeA = np.abs(y - L_A).argmin()

    for i in range(len(omegaL)) :
        if i%100 == 0 :
            print('ite :', i, 'sur', len(omegaL))
        simu = TaylorGoldstein(y, res)
        simu.resolution(k, omegaL[i])
        T_1, T_2, T_tot = simu.get_results(5, indeA, 5, 5)
        T.append(T_tot)


    return omegaL, np.array(T)







def plot_transmi(tup, npt, color='purple'):
    if len(tup)==3 :
        o4, t4 = get_transmi(tup, npt)
    else :
        o4, t4 = tup
    affichage.simple_plot(o4, t4, color, [r'$\omega$', 'T'], None, nbplot=3, cs='%.0f', csx='%.1f')
    plt.ylim(-0.05, 3.5)
    plt.xlim(0.1, 1)
    plt.tight_layout()



def get_SurfState(tup, npt):
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
    if len(tup)==5 :
        OMEG, ZZZ, L_phi, N_lat, zzz = tup
    else :
        OMEG, ZZZ, L_phi, N_lat, zzz = get_SurfState(tup, npt)
    
    print('omega_max', get_ValueO_Surf(OMEG, ZZZ, L_phi))

    ax = [r'$\omega$', 'y']
    cs='%.1f'
    cb = [r'|$\Phi$|', 0, maxx]
    cp='inferno_r'

    affichage.graph2D(OMEG, ZZZ[::-1], L_phi, 0, cbar=cb, axes=ax, nbplot=4, cs=cs, cmap=cp, shape=(17, 10))
    plt.xlim(OMEG.min(), OMEG.max())
    plt.ylim(percent*ZZZ.max(), ZZZ.max())
    plt.plot(N_lat[::-1]*100-50, zzz, '--', color='k')
    plt.tight_layout()



def get_bande(freq, trans, retuO=None, plott=True) :
    depart=0.5
    indePart = np.abs(freq - depart).argmin()
    while trans[indePart]<0.05 :
        indePart=indePart+1
    o1m = freq[indePart]
    
    depart=0.5
    indePart = np.abs(freq - depart).argmin()
    while trans[indePart]<0.025 :
        indePart=indePart-1
    o2p = freq[indePart]
    
    depart=0.9
    indePart = np.abs(freq - depart).argmin()
    while trans[indePart]<0.05 :
        indePart=indePart-1
    o1p = freq[indePart]
    
    depart=0.2
    indePart = np.abs(freq - depart).argmin()
    while trans[indePart]<0.025 :
        indePart=indePart+1
    o2m = freq[indePart]
    
    if plott :
        plt.axvline(o2m, color='purple', linestyle='--')
        plt.axvline(o2p, color='purple', linestyle='--')
        plt.axvspan(o2m, o2p, facecolor='none', alpha=1, edgecolor='magenta',lw=0, hatch='/', label="Mean Zone")

        
        plt.axvline(o1m, color='darkviolet', linestyle='--')
        plt.axvline(o1p, color='darkviolet', linestyle='--')
        plt.axvspan(o1m, o1p, facecolor='none', alpha=1, edgecolor='violet',lw=0, hatch='\\', label="Mean Zone")

        print(f"zone 1 : {o1m:.2f} à {o1p:.2f}")
        print(f"zone 2 : {o2m:.2f} à {o2p:.2f}")


    if retuO is not None :
        return o1m, o1p, o2m, o2p



def plot_phi_bande(tup, freq, trans, percent=0.75):
    y_lat, N_lat, bottt = tup
    o1m, o1p, o2m, o2p = get_bande(freq, trans, retuO=True, plott=False)


    SIM1m = TaylorGoldstein(y_lat, N_lat)
    zzz, P1m = SIM1m.resolution(k, o1m, retuP=True)
    
    SIM2m = TaylorGoldstein(y_lat, N_lat)
    zzz, P2m = SIM2m.resolution(k, o2m, retuP=True)
    
    SIM1p = TaylorGoldstein(y_lat, N_lat)
    zzz, P1p = SIM1p.resolution(k, o1p, retuP=True)
    
    SIM2p = TaylorGoldstein(y_lat, N_lat)
    zzz, P2p = SIM2p.resolution(k, o2p, retuP=True)
    
    fig, axs = plt.subplots(1, 4, figsize=(28, 10))
    axs = axs.flatten()

    axs[0].plot(P1m[::-1], y_lat, '-', color='darkviolet', label='t 1', lw=3, markersize=15)
    axs[0].set_title(r'$\phi_1^-$', fontsize='35')
    moxx = max(abs(P1m.max()), abs(P1m.min()))
    axs[0].set_xlim(-1.1*moxx, 1.1*moxx)
    
    axs[1].plot(P1p[::-1], y_lat, '-', color='darkviolet', label='t 2', lw=3, markersize=15)
    axs[1].set_title(r'$\phi_1^+$', fontsize='35')
    moxx = max(abs(P1p.max()), abs(P1p.min()))
    axs[1].set_xlim(-1.1*moxx, 1.1*moxx)
    
    axs[2].plot(P2m[::-1], y_lat, '-', color='purple', label='T tot', lw=3, markersize=15)
    axs[2].set_title(r'$\phi_2^-$', fontsize='35')
    moxx = max(abs(P2m.max()), abs(P2m.min()))
    axs[2].set_xlim(-1.1*moxx, 1.1*moxx)
    
    axs[3].plot(P2p[::-1], y_lat, '-', color='purple', label='T tot', lw=3, markersize=15)
    axs[3].set_title(r'$\phi_2^+$', fontsize='35')
    moxx = max(abs(P2p.max()), abs(P2p.min()))
    axs[3].set_xlim(-1.1*moxx, 1.1*moxx)
    
    axs[0].plot(-P1m[::-1], y_lat, ':', color='blue', label='t 1', lw=1.5, markersize=15)
    axs[1].plot(-P1p[::-1], y_lat, ':', color='blue', label='t 2', lw=1.5, markersize=15)
    axs[3].plot(-P2p[::-1], y_lat, ':', color='blue', label='t 2', lw=1.5, markersize=15)
    axs[2].plot(-P2m[::-1], y_lat, ':', color='blue', label='t 2', lw=1.5, markersize=15)
        
    for iii in range(len(axs)) :
        axs[iii].plot(N_lat[::-1]*100-50, zzz, '--', color='k', lw=1)
        axs[iii].axhline(y_lat[-1]-0.1, color='red', lw=3)
        axs[iii].axhline(y_lat[-1]-0.2, color='red', lw=3)
        axs[iii].axhline(y_lat[-1]-0.3, color='red', lw=3)
        axs[iii].axhline(y_lat[-1]-0.4, color='red', lw=3)
        axs[iii].set_ylim(percent*y_lat[-1], y_lat[-1])
        axs[iii].grid()
        axs[iii].tick_params(axis='x', which='major', labelsize=28)
        axs[iii].tick_params(axis='y', which='major', labelsize=0)
        if percent != 0.75 :
            axs[iii].axhline(zzz[-1]-0.5, color='red', lw=3)
            axs[iii].axhline(zzz[-1]-0.6, color='red', lw=3)
            axs[iii].axhline(zzz[-1]-0.7, color='red', lw=3)
            axs[iii].axhline(zzz[-1]-0.8, color='red', lw=3)
        


def get_phi_noabs(tup, npt):
    y_lat, N_lat, LB = tup

    omegaL = np.linspace(0.1, 1, npt)
    L_phi = []

    for omega in omegaL :
        simu = TaylorGoldstein(y_lat, N_lat)
        zzz, phi = simu.resolution(k, omega, retuP=True)
        L_phi.append(np.real(phi))

    OMEG, ZZZ = np.meshgrid(np.array(omegaL), zzz)

    L_phi = np.array(L_phi).T
    L_phi = L_phi
    
    return OMEG, ZZZ, L_phi, N_lat, zzz

def plot_phi_noabs(tup, npt, maxx=5, percent=0.75):
    if len(tup)==5 :
        OMEG, ZZZ, L_phi, N_lat, zzz = tup
    else :
        OMEG, ZZZ, L_phi, N_lat, zzz = get_phi_noabs(tup, npt)

    ax = [r'$\omega$', 'y']
    cs='%.1f'
    cb = [r'|$\Phi$|', -maxx, maxx]
    cp='jet'
    cp=cmpw.neonCUSTOM('green', 'blue')

    affichage.graph2D(OMEG, ZZZ[::-1], L_phi, 0, cbar=cb, axes=ax, nbplot=4, cs=cs, cmap=cp, shape=(17, 10))
    plt.xlim(OMEG.min(), OMEG.max())
    plt.ylim(percent*ZZZ.max(), ZZZ.max())
    plt.axhline(zzz[-1]-0.1, color='red', lw=3)
    plt.axhline(zzz[-1]-0.2, color='red', lw=3)
    plt.axhline(zzz[-1]-0.3, color='red', lw=3)
    plt.axhline(zzz[-1]-0.4, color='red', lw=3)
    if percent != 0.75 :
        plt.axhline(zzz[-1]-0.5, color='red', lw=3)
        plt.axhline(zzz[-1]-0.6, color='red', lw=3)
        plt.axhline(zzz[-1]-0.7, color='red', lw=3)
        plt.axhline(zzz[-1]-0.8, color='red', lw=3)
    plt.tight_layout()





def RES_12():
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_1(),MARCHE_1(),MARCHE_1(),MARCHE_1(), MARCHE_2(),MARCHE_2(),MARCHE_2(),MARCHE_2(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT

def RES_21():
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_2(),MARCHE_2(),MARCHE_2(),MARCHE_2(), MARCHE_1(),MARCHE_1(),MARCHE_1(),MARCHE_1(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    print("- points :", len(reseau))
    
    return y, reseau, L_BOT



def RES_81():
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_1(),MARCHE_1(),MARCHE_1(),MARCHE_1(), MARCHE_1(),MARCHE_1(),MARCHE_1(),MARCHE_1(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT

def RES_82():
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_2(),MARCHE_2(),MARCHE_2(),MARCHE_2(), MARCHE_2(),MARCHE_2(),MARCHE_2(),MARCHE_2(), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT




def RES_test(l_1):
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    Layer_1 = lat_cte(l_1, dy, 0)
    Layer_2 = lat_cte(0.03333333333333, dy, 1)
      
    L_HAUT = l_1+0.033
    L_BOT = Ly - L_HAUT
    


    
    lbot = lat_cte(L_BOT, dy, 0)
    reseau = np.concatenate((Layer_1, Layer_2, lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    return y, reseau, L_BOT


def get_ValueO_Surf(OMEG, ZZZ, L_phi):
    indexZ = np.where(L_phi==L_phi.max())[0][0]
    indexO = np.where(L_phi==L_phi.max())[1][0]
    return OMEG[indexZ, indexO]



def SurfaceLat(ln2):
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    n = 4
    L = [0.1]*n
    R = [1/3]*n
    L_N2 = ln2
    L_TOP = 0.033
    L_HAUT = L_N2 + L_TOP + np.sum(L)
    L_BOT = Ly - L_HAUT
    res = res_I(n, L, R, dy, 1, 0)
    ln2 = lat_cte(L_N2, dy, 0)
    ltop = lat_cte(L_TOP, dy, 1)
    lbot = lat_cte(L_BOT, dy, 1)
    LATTICE = np.concatenate((ln2, ltop, res, lbot))[::-1]
    pt = len(LATTICE)
    y = np.linspace(0, Ly, pt)
    
    
    return y, LATTICE[::-1], L_BOT


def SurfaceIN(lin):
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    n = 4
    L = [0.1]*n
    R = [1/3]*n
    L_IN = lin
    L_N2 = 0.03
    L_TOP = 0.033
    L_HAUT = L_N2 + L_TOP + np.sum(L) + L_IN
    L_BOT = Ly - L_HAUT
    res = res_I(n, L, R, dy, 1, 0)
    l_in = lat_cte(L_IN, dy, 1)
    ln2 = lat_cte(L_N2, dy, 0)
    ltop = lat_cte(L_TOP, dy, 1)
    lbot = lat_cte(L_BOT, dy, 1)
    LATTICE = np.concatenate((l_in, ln2, ltop, res, lbot))[::-1]
    pt = len(LATTICE)
    y = np.linspace(0, Ly, pt)
    
    
    return y, LATTICE[::-1], L_BOT




def plot_marche(marche):
    absc=np.linspace(0, 0.1, len(marche))
    affichage.simple_plot(absc, marche, 'red', [None, None], None)


# %%

# =============================================================================
# Variation simultanée de 2 cristaux différents
# =============================================================================


dy =0.0007

def verif(ratio1, ratio2):
    testtte = 2*ratio1+ratio2
    if testtte !=1 :
        raise ValueError("pas bon ratio", testtte)


def MARCHE_1_bis(ratio1, ratio2):
    haut = lat_cte(ratio1*0.1, dy, 1)
    bas = lat_cte(ratio2*0.1, dy, 0)
    
    marche_1 = np.concatenate((haut, bas, haut))
    return marche_1

def MARCHE_2_bis(ratioA, ratioB):
    haut = lat_cte(ratioB*0.1, dy, 1)
    bas = lat_cte(ratioA*0.1, dy, 0)
    
    
    marche_2 = np.concatenate((bas, haut, bas))

    return marche_2

def RES_21_RATIO(ratio1, ratio2, ratioA, ratioB):
    verif(ratio1, ratio2)
    verif(ratioA, ratioB)
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_2_bis(ratioA, ratioB),MARCHE_2_bis(ratioA, ratioB),MARCHE_2_bis(ratioA, ratioB),MARCHE_2_bis(ratioA, ratioB),MARCHE_1_bis(ratio1, ratio2),MARCHE_1_bis(ratio1, ratio2),MARCHE_1_bis(ratio1, ratio2),MARCHE_1_bis(ratio1, ratio2), lbot))
    
    y = np.linspace(0, Ly, len(reseau))
    
    print("- points :", len(reseau))
    
    return y, reseau, L_BOT


plot_res(RES_21_RATIO(16/48, 16/48, 8/48, 32/48), percent=0.6)
plot_res(RES_21(), percent=0.6)



ratio1, ratio2 = 16/48, 16/48
ratioA, ratioB = 8/48, 32/48
LAT = RES_21_RATIO(ratio1, ratio2, ratioA, ratioB)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6, maxx=1)

ratio1, ratio2 = 15/48, 18/48
ratioA, ratioB = 9/48, 30/48
LAT = RES_21_RATIO(ratio1, ratio2, ratioA, ratioB)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6, maxx=1)

ratio1, ratio2 = 14/48, 20/48
ratioA, ratioB = 10/48, 28/48
LAT = RES_21_RATIO(ratio1, ratio2, ratioA, ratioB)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6, maxx=1)

ratio1, ratio2 = 13/48, 22/48
ratioA, ratioB = 11/48, 26/48
LAT = RES_21_RATIO(ratio1, ratio2, ratioA, ratioB)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6, maxx=1)

ratio1, ratio2 = 12/48, 24/48
ratioA, ratioB = 12/48, 24/48
LAT = RES_21_RATIO(ratio1, ratio2, ratioA, ratioB)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6, maxx=1)



# %%

# =============================================================================
# Variation simultanée de 2 cristaux avec une origine commune
# =============================================================================

dy=0.0007

def MARCHE_1_R(ratio1):
    
    ratio2 = 1 - ratio1 - ratio1
    
    haut = lat_cte(ratio1*0.1, dy, 1)
    bas = lat_cte(ratio2*0.1, dy, 0)
    
    marche_1 = np.concatenate((haut, bas, haut))
    # print(len(marche_1)*dy)
    return marche_1




def RES_11_Difer(ratio1, ratio2):
 
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_1_R(ratio1),MARCHE_1_R(ratio1),MARCHE_1_R(ratio1),MARCHE_1_R(ratio1),MARCHE_1_R(ratio2),MARCHE_1_R(ratio2),MARCHE_1_R(ratio2),MARCHE_1_R(ratio2),  lbot))
    
    y = np.linspace(0, Ly, len(reseau))

    
    return y, reseau, L_BOT




ratio1, ratio2 = 16/48, 16/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm1.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 17/48, 14/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm2.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 18/48, 12/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm3.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 19/48, 10/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm4.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 20/48, 8/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm5.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 21/48, 6/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm6.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 22/48, 4/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm7.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 23/48, 2/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm8.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 24/48, 0/48
LAT = RES_11_Difer(ratio1, ratio2)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/norm9.png', bbox_inches='tight', dpi=100)

plt.show()


plt.figure()
plt.plot([], [])
plt.title("Cas inverse")
plt.show()

ratio1, ratio2 = 16/48, 16/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv1.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 17/48, 14/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv2.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 18/48, 12/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv3.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 19/48, 10/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv4.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 20/48, 8/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv5.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 21/48, 6/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv6.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 22/48, 4/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv7.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 23/48, 2/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv8.png', bbox_inches='tight', dpi=100)


ratio1, ratio2 = 24/48, 0/48
LAT = RES_11_Difer(ratio2, ratio1)
# plot_res(LAT, percent=0.6)
plot_SURF_STATE(LAT, 200, percent=0.6)
plt.savefig('IMag/inv9.png', bbox_inches='tight', dpi=100)

plt.show()



# %%

# =============================================================================
# Déphasage d'un même cristal
# =============================================================================

dy=0.0007

def confirm_rat(marche):
    pt0 = 0
    pt1 = 1
    pt_tot = len(marche)
    
    for i in range(len(marche)):
        if marche[i]==0 :
            pt0 += 1
        else :
            pt1 += 1
    if pt0/pt_tot != 1/3 :
        raise ValueError("probleme de resolut", pt0/pt_tot)


def MARCHE_1_Deph(phase):
    L_total = 0.1
    L = L_total / 3
    haut1 = lat_cte(L, dy, 1)
    bas = lat_cte(L, dy, 0)
    haut2 = lat_cte(L, dy, 1)

    y_period = np.concatenate([haut1, bas, haut2])
    n_pts = len(y_period)
    
    x_period = np.linspace(0, L_total, n_pts, endpoint=False)

    x_shifted = (x_period + phase * L_total) % L_total
    
    sort_idx = np.argsort(x_shifted)
    x_sorted = x_shifted[sort_idx]
    y_sorted = y_period[sort_idx]

    f_shifted = np.interp(x_period, x_sorted, y_sorted)
    
    for i in range(len(f_shifted)):
        if f_shifted[i]!=0 and f_shifted[i]!=1 :
            if f_shifted[i] < 0.5 :
                f_shifted[i] = 0
            else : 
                f_shifted[i] = 1
    
    return f_shifted


def RES_11_dephase(phase1, phase2):
 
    k=27.6
    Ly = 2*np.pi / k * 21 / 2
    
    L_HAUT = 0.8
    L_BOT = Ly - L_HAUT
    
    lbot = lat_cte(L_BOT, dy, 1)
    reseau = np.concatenate((MARCHE_1_Deph(phase1),MARCHE_1_Deph(phase1),MARCHE_1_Deph(phase1),MARCHE_1_Deph(phase1),MARCHE_1_Deph(phase2),MARCHE_1_Deph(phase2),MARCHE_1_Deph(phase2),MARCHE_1_Deph(phase2),  lbot))
    
    y = np.linspace(0, Ly, len(reseau))

    confirm_rat(MARCHE_1_Deph(phase1))
    confirm_rat(MARCHE_1_Deph(phase2))
    
    return y, reseau, L_BOT


# plot_res(RES_11_dephase(0, 0.5), percent=0.6)
# plot_res(RES_12(), percent=0.6)

# phase1=0
# LAT=RES_11_dephase(0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p0.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.05
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p05.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.1
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p1.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.15
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p15.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.2
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p2.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.25
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p25.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.3
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p3.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.35
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p35.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.4
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p4.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.45
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p45.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.5
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p5.png', bbox_inches='tight', dpi=100)
# plt.show()



# phase1=0.55
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p55.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.6
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p6.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.65
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p65.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.7
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p7.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.75
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p75.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.8
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p8.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.85
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p85.png', bbox_inches='tight', dpi=100)
# plt.show()

# phase1=0.9
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p9.png', bbox_inches='tight', dpi=100)
# plt.show()


# phase1=0.95
# LAT=RES_11_dephase(-0, phase1)
# plot_SURF_STATE(LAT, 200, percent=0.6, maxx=2.5)
# plt.savefig('PHA_REL/T12_ph0p95.png', bbox_inches='tight', dpi=100)
# plt.show()



# %%

plot_marche(MARCHE_1_Deph(0))
plot_marche(MARCHE_1_Deph(0.1))
plot_marche(MARCHE_1_Deph(0.2))
plot_marche(MARCHE_1_Deph(0.3))
plot_marche(MARCHE_1_Deph(0.4))
plot_marche(MARCHE_1_Deph(0.5))



