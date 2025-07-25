# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 16:22:27 2025

@author: Brancher Paul

Tutorial for using the two functions of the XDMF_reader class, giving access
to the ux, uy, pp, and phi (scalar 1) fields as 2D numpy arrays.

"""

import numpy as np
from Reader import XDMF_reader
from Affichage import affichage
import matplotlib.pyplot as plt
from CMAP import cmpw


# =============================================================================
# Champs à un temp donné
# =============================================================================


ux, uy, pp, phi = XDMF_reader.load_data('snapshot-100.xdmf', basepath='data_0p7/', verbose=True)

i_y = np.linspace(1, np.shape(ux)[0], np.shape(ux)[0])
i_x = np.linspace(1, np.shape(ux)[1], np.shape(ux)[1])

k=27.6
A=0.003
Lx = 2*np.pi / k * 21 / A
Ly = 2*np.pi / k * 21 / A / 2
dx = Lx / (len(i_x) - 1)
dy = Ly / (len(i_y) - 1)
x = (i_x-1) * dx
y = (i_y-1) * dy
lambX = 2*np.pi / k / A

# exemple de plot
affichage.graph2D(x, y, uy, 100, cmap=cmpw.neonCUSTOM('red', 'blue'), axes='zero', shape=(17,10), nbplot=4, cs='%.2f', cbar='ok')



# =============================================================================
# Listes 1D des champs à chaque pas de temps
# =============================================================================


UX, UY, PP, PHI = XDMF_reader.timed_data('data_0p7/', snapMAX=100)

# exemple de plot pour la snapshot 33
affichage.graph2D(x, y, PHI[33], 100, cmap=cmpw.neonCUSTOM('orange', 'blue'), axes='zero', shape=(17,10), nbplot=4, cs='%.2f', cbar='ok')
