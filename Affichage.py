# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:56:36 2024

@author: Paul BRANCHER

Code développé par Paul BRANCHER depuis octobre 2024
Classe permettant l'automatisation du processus de génération de belles figures

VERSION DU 25/07/25 :
==> le traitement de |phi|(y, omega) est possible mais pas automatisé
==> i.e UN ENORME IF PAS TRES JOLI DANS LE CODE
==> Mais ça marche bien comme ça...

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.integrate import solve_ivp
from matplotlib.tri import Triangulation
from scipy.interpolate import RectBivariateSpline
import matplotlib.animation as animation
import matplotlib.ticker as ticker



class affichage(object):
    """ Classe comprenant les macros de tout les tracés"""
    
    @staticmethod
    def configTeX(number) :
        if number==2 :
            gradsize = 27           
            lablsize = gradsize+3
            legendsize=28
            titlesize=30
        
        elif number==3:
            gradsize = 40           
            lablsize = gradsize+3
            legendsize=38
            titlesize=41
        
        elif number==4 :
            gradsize = 50
            lablsize = gradsize+3
            lablsize=0
            legendsize=55
            titlesize=52
            
        return gradsize, lablsize, legendsize, titlesize
    
    @staticmethod
    def getSize(cbar=None, shape=None, titre=None, nbplot=None):
        """ pose la size de la figure proprement"""
        if nbplot is None :
            nbplot=2
        gradsize, lablsize, legendsize, titlesize = affichage.configTeX(nbplot)
        orient='vertical'

        if shape is None :
            shape=(10, 10)
            if cbar is None :
                fig = plt.figure(figsize=shape)

            elif cbar is not None :
                shape=(12, 10)
                fig = plt.figure(figsize=shape)
 
        elif shape=='pure' :
            shape=(10, 10)
            gradsize = 0
            lablsize = 0
            if cbar is None :
                fig = plt.figure(figsize=shape)
            
            elif cbar is not None :
                shape=(12,10)
                fig = plt.figure(figsize=shape)

        else :
            fig = plt.figure(figsize=shape)

        if titre is not None :
            fig = plt.title(titre, fontsize=titlesize)

        return fig, gradsize, lablsize, legendsize, orient

    
    @staticmethod
    def getcbar(cbar, form, orient, contour, shape, Z, cs, gradsize, lablsize, nbplot):
        
        co_min = Z.min()
        co_max = Z.max()
        
        if type(cbar) == list :
            co_min = cbar[1]
            co_max=cbar[2]
            cbar = cbar[0]
        
        if type(cbar) == str and cbar[0]!='?':
            if type(form)==dict :
                colorbar = plt.colorbar(orientation=orient, label=cbar)
            else :
                colorbar = plt.colorbar(contour, orientation=orient, label=cbar)      

            if shape=='pure':
                gradsize = 50
                lablsize = gradsize+3
                lablsize=0
                legendsize=55
                titlesize=52
            
            contour.set_clim(co_min, co_max)
            colorbar.ax.tick_params(labelsize=gradsize)
            colorbar.ax.yaxis.label.set_size(lablsize)
            colorbar.set_ticks(np.linspace(co_min, co_max, 5))
            colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter(cs))
            
            
            if cbar==r'|$\Phi$|' :
                contour.set_clim(co_min, co_max)
                colorbar.ax.tick_params(labelsize=gradsize)
                colorbar.ax.yaxis.label.set_size(gradsize)
                colorbar.set_ticks(np.linspace(co_min, co_max, 4))
                colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter(cs))
                colorbar.set_label(r'|$\Phi$|', labelpad=-50)
                

        elif type(cbar) == str and cbar[0]=='?' :
            lettre = cbar[1]
            var=lettre[2]
            cbar=cbar[3:]
            print(cbar)
            colorbar = plt.colorbar(contour)

            Mini = Z.min()
            Maxi = Z.max()
            Moyen = (Mini + Maxi) / 2
            positions = [Mini, Moyen, Maxi]

            graduation = [cs % Z.min(), cbar, r"$\lim_{y \to 0} %s$" % lettre]
            colorbar.set_ticks(positions)
            colorbar.set_ticklabels(graduation)

            titre = colorbar.ax.get_yticklabels()[1]
            titre.set_rotation(90)
            titre.set_horizontalalignment('left')  # Centrer horizontalement les étiquettes
            titre.set_verticalalignment('center')
            titre.set_position((titre.get_position()[0] + 0.5,
                                titre.get_position()[1]))

            for grad, size in zip(colorbar.ax.get_yticklabels(),
                                  [gradsize, lablsize, lablsize+2]):
                grad.set_fontsize(size)
        
        if Z.max() < 2*10**(-3):
            exponent = -int(cs[2])
            scale = 10 ** exponent
            def sci_tick_formatter(val, pos):
                return f"{val / scale:.2f}"

            colorbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(sci_tick_formatter))

            sci_label = r"$\times 10^{{{}}}$".format(exponent)
            if shape==None and nbplot==4 :
                colorbar.ax.text(-0.975, 0.5, sci_label, transform=colorbar.ax.transAxes,
                                 fontsize=gradsize/1.35, va='center', ha='left', rotation=90)  
            elif shape==None and nbplot==3 :
                colorbar.ax.text(-0.9, 0.5, sci_label, transform=colorbar.ax.transAxes,
                                 fontsize=gradsize/1.1, va='center', ha='left', rotation=90)
            elif shape==None and nbplot==2 :
                colorbar.ax.text(-0.8, 0.5, sci_label, transform=colorbar.ax.transAxes,
                                 fontsize=gradsize, va='center', ha='left', rotation=90)
            else :
                colorbar.ax.text(-1.4, 0.5, sci_label, transform=colorbar.ax.transAxes,
                                 fontsize=gradsize/1.2, va='center', ha='left', rotation=90)


    @staticmethod
    def getplot(X, Y, Z, form, triang=None, cmap=None):
        """ choisis le mode de plot : sur mesh triang ou sur mesh equireparti """

        # norm = TwoSlopeNorm(vmin=Z.min(), vcenter=0, vmax=Z.max())
        norm = None

        if triang is None:
            if form == 0:
                contour = plt.pcolor(X, Y, Z, cmap=cmap, shading='auto', norm=norm)
            elif type(form)==str and form[0]!='b' and form[0]!='p' :
                form=int(form)
                contour = plt.contour(X, Y, Z, levels=form, cmap=cmap)
            elif type(form)==dict :
                npdf = form["bins"]
                density = form["density"]
                contour = plt.hist2d(X, Y, bins=npdf, cmap = cmap, density=density)
            elif type(form)==str and form[0]=='b' :
                form=form[1:]
                form=int(form)
                contour = plt.contourf(X, Y, Z, levels=form, cmap=cmap)
                plt.contour(X, Y, Z, levels=form, colors='k')
            elif type(form)==str and form[0]=='p' :
                form=form[1:]
                form=int(form)
                contour = plt.pcolor(X, Y, Z, cmap=cmap, shading='auto', norm=norm)
                plt.contour(X, Y, Z, levels=form, colors='k', linewidths=1)
                # plt.contour(X, Y, Z, levels=form, cmap='seismic_r', linewidths=1)
            else:
                contour = plt.contourf(X, Y, Z, levels=form, cmap=cmap)
            
        
        if triang is not None:
            triangle=Triangulation(X, Y)
            if form == 0:
                contour = plt.tripcolor(triangle, Z, cmap=cmap, shading='auto')
            elif type(form)==str and form[0]!='b' :
                form=int(form)
                contour = plt.tricontour(triangle, Z, levels=form, cmap=cmap)
            else:
                contour = plt.tricontourf(triangle, Z, levels=form, cmap=cmap)
            
        return contour


    @staticmethod
    def graph2D(X, Y, Z, form, cmap=None, triang=None, cs=None, cbar=None, titre=None, shape=None, axes=None, nbplot=None):
        if cs is None:
            cs = '%.3f'
            
        if cmap is None :
            cmap='jet'
            
        fig, gradsize, lablsize, legendsize, orient = affichage.getSize(cbar, shape, titre, nbplot)
        
        contour = affichage.getplot(X, Y, Z, form, triang, cmap)
        
        gr = gradsize
        
        if axes=='zero':
            gr = 0
            axes = False
            
        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=gr)
        ax.tick_params(axis='y', labelsize=gr)

        if axes is None :
            if shape!='pure':
                plt.xlabel('x', fontsize='%i' % lablsize, labelpad=-10)
                plt.ylabel('y', fontsize='%i' % lablsize, labelpad=-10)
                
        
        elif axes is not None and axes!=False:
            if shape!='pure':
                plt.xlabel(axes[0], fontsize='%i' % lablsize, labelpad=-10)
                # plt.xlabel(axes[0], fontsize='%i' % lablsize)
                plt.ylabel(axes[1], fontsize='%i' % lablsize, labelpad=-10)
                if axes[0]==r'$\omega$' :
                    plt.xlabel(axes[0], fontsize='%i' % gradsize, labelpad=-15, loc='center')
                    plt.ylabel(axes[0], fontsize='%i' % 0, labelpad=-15, loc='center')
                    # ax.xaxis.set_label_coords(0.99, 0.075)
                    ax.tick_params(axis='y', labelsize=0)
                             

        if nbplot==3 and shape==None :
            positions = ax.get_xticks()
            labels = np.copy(positions).tolist()
            labels[1], positions[1] = 0.2, 0.2
            labels[2], positions[2] = 0.4, 0.4
            labels[0], positions[0] = 0.6, 0.6
            labels[-1], positions[-1] = 0.8, 0.8
            if cbar is not None :
                labels[3], positions[3] = -1, -1
                labels[4], positions[4] = -1, -1
            ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            
            posY = ax.get_yticks()
            llbY=[]
            for e in posY :
                llbY.append(f"{e:.1f}")
            llbY[1]=''
            ax.yaxis.set_major_locator(ticker.FixedLocator(posY))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(llbY))


        affichage.getcbar(cbar, form, orient, contour, shape, Z, cs, gradsize, lablsize, nbplot)

        plt.tight_layout()


    @staticmethod
    def simple_plot(x, y, couleur, axes, labels, style=None, nogrid=None, lw=None, titre=None, nbplot=None, shape=None, cs=None, csx=None):
        cbar=None
        if shape is None :
            shape=(12, 10)

        if lw is None:
            lw = [3]*len(y)
            
        if cs is None:
            cs = False
            
        if csx is None:
            csx = False
        
        nb=10

        if type(x) != list and type(y) == list:
            nb=len(y)
            x = [x]*len(y)
            if style is None :
                style=['xy']*nb

                
        if type(y) != list and type(x) == list:
             y = [y]*len(x)
             nb=len(x)
             if style is None :
                 style=['xy']*len(x)
                

        if style is None:
            if type(x)==list:
                style=['xy']*nb
            else:
                style = 'xy'

        if type(y) != list:
            y = [y]
            x = [x]
            couleur = [couleur]
            labels = [labels]
            style = [style]

        fig, gradsize, lablsize, legendsize, orient = affichage.getSize(cbar, shape, titre, nbplot)
        
        if axes !=False and axes != 'zero' :
            plt.ylabel(axes[1], fontsize=lablsize)
            plt.xlabel(axes[0], fontsize=lablsize)
        
        plt.tick_params(axis='both', which='major', labelsize=gradsize)
        if axes == 'zero' :
            plt.tick_params(axis='y', which='major', labelsize=0)
            plt.xlabel('N', fontsize=lablsize)

        ax = plt.gca()
        if cs !=False :
            ax.yaxis.set_major_formatter(FormatStrFormatter(cs))
        
        if csx !=False :
            ax.xaxis.set_major_formatter(FormatStrFormatter(csx))

        lst_p1=['d', 's', 'o', 'd', '.', '1', '2', 'p', 'h', 'H', 'x', 'X', '*']
        lst_p2=[]
        for elem in lst_p1:
            lst_p2.append('-%s'%elem)
            lst_p2.append('--%s'%elem)
        lst_p=lst_p1+lst_p2

        for i in range(len(labels)):
            if style[i] == "xy":
                plt.plot(x[i], y[i], color=couleur[i], label=labels[i],
                         linewidth=lw[i])
            elif style[i] in lst_p:
                if style[i][-1]=='*':
                    mksize=lw[i]*7
                else :
                    mksize=lw[i]*5
                plt.plot(x[i], y[i], style[i], color=couleur[i],
                         label=labels[i], markersize=mksize, linewidth=lw[i])
            elif type(style[i]) == dict :
                dico = style[i]
                density= dico['density']
                bins = dico['bins']
                plt.hist(x[i], color=couleur[i], label=labels[i], density=density, bins=bins)
            else:
                plt.plot(x[i], y[i], linestyle=style[i], color=couleur[i],
                         label=labels[i], linewidth=lw[i])
        
        if nogrid is None :
            plt.grid()
        if labels[0] is not None :
            plt.legend(fontsize=legendsize, loc='upper left')
        if len(labels)>=4 :
            plt.legend(fontsize=legendsize, ncol=2)
        plt.tight_layout()


    @staticmethod
    def graph3D(X, Y, Z, labl, cmap, cs=None, cbar=None, titre=None):
        
        print("3D : ", np.shape(X), len(np.shape(X)))
        if len(np.shape(X))==1 :
            X, Y = np.meshgrid(X, Y)
            
        if cs is None:
            cs = '%.3f'

        fig = plt.figure(figsize=(10, 10))

        ax = fig.add_subplot(111, projection='3d')
        if titre is not None:
            # ax.title.set_text(titre, fontsize='24')
            ax.set_title(titre, fontsize='23')
            # plt.title(titre, fontsize=28)

        gradsize = 23
        lablsize = gradsize+1

        surf = ax.plot_surface(X, Y, Z, cmap=cmap, edgecolor='none')

        ax2 = plt.gca()
        ax2.tick_params(axis='x', labelsize=gradsize-3)
        ax2.tick_params(axis='y', labelsize=gradsize-3)
        ax2.tick_params(axis='z', labelsize=gradsize-3)

        if cbar is True :
            ax.set_xlabel('X', labelpad=19, fontsize=lablsize)
            ax.set_ylabel('Y', labelpad=15, fontsize=lablsize)
            colorbar = fig.colorbar(surf, ax=ax, shrink=0.7)
            colorbar.set_label(labl, fontsize=lablsize)
            colorbar.ax.tick_params(labelsize=gradsize-3)
            colorbar.set_ticks(np.linspace(Z.min(), Z.max(), 10))
            colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter(cs))
            ax.set_box_aspect([15, 15, 8])
            ax.view_init(elev=25, azim=-70)
        elif type(cbar) == str:
            ax.set_xlabel('X', labelpad=19, fontsize=lablsize)
            ax.set_ylabel('Y', labelpad=9, fontsize=lablsize)
            colorbar = fig.colorbar(surf, ax=ax, shrink=0.7)
            Mini = Z.min()
            Maxi = Z.max()
            Moyen = (Mini + Maxi) / 2
            positions = [Mini, Moyen, Maxi]

            graduation = [cs % Z.min(), labl, r"$\lim_{y \to 0} %s$" % cbar]
            colorbar.set_ticks(positions)
            colorbar.set_ticklabels(graduation)

            titre = colorbar.ax.get_yticklabels()[1]
            titre.set_rotation(90)
            titre.set_horizontalalignment('left')  # Centrer horizontalement les étiquettes
            titre.set_verticalalignment('center')
            titre.set_position((titre.get_position()[0] + 0.5,
                                titre.get_position()[1]))

            for grad, size in zip(colorbar.ax.get_yticklabels(),
                                  [gradsize, lablsize, lablsize+2]):

                grad.set_fontsize(size)

            ax.set_box_aspect([25, 13, 8])
            ax.view_init(elev=30, azim=-75)

        plt.tight_layout()


    @staticmethod
    def preprocess_for_quiver(x, y, u, v):
        """
        Prépare les données pour quiver, en détectant automatiquement si elles sont sur une grille régulière.
        Parameters:
            x (np.ndarray): Coordonnées x (1D)
            y (np.ndarray): Coordonnées y (1D)
            u (np.ndarray): Composantes du vecteur en x (1D) ou (2D)
            v (np.ndarray): Composantes du vecteur en y (1D) ou (2D)
        Returns:
            tuple: (X, Y, U, V) si les données sont régulières (2D),
                   (x, y, u, v) sinon (1D).
        """

        # Identifier les valeurs uniques de x et y
        x_unique = np.unique(x)
        y_unique = np.unique(y)

        # Vérifier si les valeurs uniques sont équiréparties
        x_diff = np.diff(x_unique)
        y_diff = np.diff(y_unique)

        if np.allclose(x_diff, x_diff[0]) and np.allclose(y_diff, y_diff[0]):
            # Grille régulière détectée
            X, Y = np.meshgrid(x_unique, y_unique)
            try:
                U = u.reshape(len(y_unique), len(x_unique))
                V = v.reshape(len(y_unique), len(x_unique))
                print('Données modif pour correspondre')
                return X, Y, U, V
            except ValueError:
                # En cas d'erreur de reshape, retourner les données telles quelles
                print("Erreur : Les dimensions des vecteurs u/v ne correspondent pas à une grille.\
                      \n il faut que u et v soient  2D !!!")
                return x, y, u, v
        else:
            # Données non régulières
            print('Données pas modif')
            return x, y, u, v


    @staticmethod
    def vectorize(x, y, u, v, step, cmap=None, scale=None, cs=None, cbar=None, titre=None, shape=None, nbplot=None, axes=None):
        """ prends des données, les rends exploitables, les plots en vect """

        x,y,u,v = affichage.preprocess_for_quiver(x, y, u, v)   # exploitabilité

        if cmap is None :
            cmap = 'jet'

        mag = np.sqrt(u**2 + v**2)
        norm = plt.Normalize(mag.min(), mag.max())
        cmap = plt.get_cmap(cmap)

        fig, gradsize, lablsize, legendsize, orient = affichage.getSize(cbar, shape, titre, nbplot)
        
        gr=gradsize
        
        if axes=='zero':
            gr = 0
            axes = False

        
        if scale is None :
            scale = 7

        if cs is None:
            cs = '%.3f'

            
        if len(x.shape) == 1:  # Cas des données 1D (non reshaped)
            colors = cmap(norm(mag))[::step]
            x, y, u, v = x[::step], y[::step], u[::step], v[::step]
        else:  # Cas des données 2D (reshaped)
            colors=cmap(norm(mag))
            print('modif des couleurs')
            colors = colors[::step, ::step]
            colors = colors.reshape(-1, 4)
            x = x[::step, ::step]
            y = y[::step, ::step]
            u = u[::step, ::step]
            v = v[::step, ::step]

        contour = plt.quiver(x, y, u, v, scale=scale, color=colors, headwidth=8, headlength=3, headaxislength=2.5)
        m=plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        m.set_array([])
        
        if axes !=False :
            plt.xlabel('x', fontsize='%i' % lablsize)
            plt.ylabel('y', fontsize='%i' % lablsize)

        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=gr)
        ax.tick_params(axis='y', labelsize=gr)
        
        if type(cbar)==str :
            colorbar = plt.colorbar(m, ax=plt.gca(), label=cbar)
            colorbar.ax.tick_params(labelsize=gradsize)
            colorbar.ax.yaxis.label.set_size(lablsize)
            colorbar.set_ticks(np.linspace(mag.min(), mag.max(), 5))
            colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter(cs))
        
        plt.xlim(x.min(), x.max())
        plt.ylim(y.min(), y.max())
        plt.tight_layout()


    @staticmethod
    def ldc(x, y, ux, uy, lenght=None, start_points=None, nb=None, lw=None, shape=None, titre=None, nbplot=None, axes=None):
        if shape is None :
            shape=(10, 10)

        if lw is None :
            lw=1

        if lenght is None :
            lenght=10
        
        if start_points is None :
            if nb is None :
                nb=5
            x_start = np.linspace(x[3], x[-3], nb)
            y_start = np.linspace(y[3], y[-3], nb)

            start_points = np.array([(xi, yi) for xi in x_start for yi in y_start])
        
        interp_U = RectBivariateSpline(y, x, ux)
        interp_V = RectBivariateSpline(y, x, uy)
        
        def velocity_field(t, pos):
            x, y = pos
            u = interp_U(y, x)[0, 0]
            v = interp_V(y, x)[0, 0]
            return [u, v]

        ABS, ORD = [], []
        COUL, LW = [], []
        if axes is None :
            AX = ['x', 'y']
        else :
            AX=axes
        
        i=0
        for start in start_points:
            print("lcd plot :", i, 'sur', len(start_points))
            i+=1
            sol = solve_ivp(velocity_field, [0, lenght], start, max_step=0.1, t_eval=np.linspace(0, lenght, 200))
            ABS.append(sol.y[0])
            ORD.append(sol.y[1])
            COUL.append('k')
            LW.append(lw)

        LLB=[None]*len(ABS)
        ST=["xy"]*len(ABS)
        affichage.simple_plot(ABS, ORD, COUL, AX, LLB, shape=shape, style=ST, lw=LW, titre=titre, nbplot=nbplot)
        plt.xlim(x.min(), x.max())
        plt.ylim(y.min(), y.max())
        plt.tight_layout()
        

    @staticmethod
    def animate2D(X, Y, data, form, cs=None, cmap=None, cbar=None, shape=None, titre=None, nbplot=None, inter=None, axes=None):
        if cs is None:
            cs = '%.3f'

        if cmap is None :
            cmap='jet'
        
        if inter is None :
            inter=10
        
        fig, gradsize, lablsize, legendsize, orient = affichage.getSize(cbar, shape, titre, nbplot)

        timesteps = len(data)
        # if axes!=False :
        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=gradsize)
        ax.tick_params(axis='y', labelsize=gradsize)
        if form == 0 :
            contour = plt.pcolor(X, Y, data[0], cmap=cmap, shading="auto")
            def update(frame):
                contour.set_array(data[frame].ravel())
                return contour,
            ani = animation.FuncAnimation(fig, update, frames=timesteps, interval=inter)
            affichage.getcbar(cbar, form, orient, contour, shape, data[0], cs, gradsize, lablsize, nbplot)
            return ani
        elif type(form)==str and form[0]!='b' :
            form=int(form)
            contour = plt.contour(X, Y, data[0], levels=form, cmap=cmap)
            def update(frame):
                nonlocal contour
                for coll in contour.collections:
                    coll.remove()
                contour = plt.contour(X, Y, data[frame], cmap=cmap)
                return contour.collections
            affichage.getcbar(cbar, form, orient, contour, shape, data[0], cs, gradsize, lablsize, nbplot)
            ani = animation.FuncAnimation(fig, update, frames=timesteps, interval=inter)
            return ani
        elif type(form)==str and form[0]=='b' :
            form=form[1:]
            form=int(form)
            contour = plt.contourf(X, Y, data[0], levels=form, cmap=cmap)
            contour2 = plt.contour(X, Y, data[0], levels=form, colors='k')
            def update(frame):
                nonlocal contour
                nonlocal contour2
                for coll in contour.collections:
                    coll.remove()
                for coll in contour2.collections:
                    coll.remove()
                contour = plt.contourf(X, Y, data[frame], levels=form, cmap=cmap)
                contour2 = plt.contour(X, Y, data[frame], levels=form, colors='k')
                return contour.collections
            affichage.getcbar(cbar, form, orient, contour, shape, data[0], cs, gradsize, lablsize, nbplot)
            ani = animation.FuncAnimation(fig, update, frames=timesteps, interval=inter)
            return ani
        else :
            form=int(form)
            contour = plt.contourf(X, Y, data[0], levels=form, cmap=cmap)
            def update(frame):
                nonlocal contour
                for coll in contour.collections:
                    coll.remove()
                contour = plt.contourf(X, Y, data[frame], levels=form, cmap=cmap)
                return contour.collections
            affichage.getcbar(cbar, form, orient, contour, shape, data[0], cs, gradsize, lablsize, nbplot)
            ani = animation.FuncAnimation(fig, update, frames=timesteps, interval=inter)
            return ani

    @staticmethod
    def plot_TRI_mesh(x, y, lw=None, titre=None, shape=None, nbplot=None):
        cbar=None
        fig, gradsize, lablsize, legendsize, orient = affichage.getSize(shape, titre, nbplot=None)
        if lw is None:
            lw=0.15
        
        
        triang = Triangulation(x, y)
        plt.triplot(triang, color='k', lw=lw)
        plt.xlim(x.min(), x.max())
        plt.ylim(y.min(), y.max())

        plt.xlabel('x', fontsize='%i' % lablsize)
        plt.ylabel('y', fontsize='%i' % lablsize)
        ax = plt.gca()
        ax.tick_params(axis='x', labelsize=gradsize)
        ax.tick_params(axis='y', labelsize=gradsize)
        
        plt.tight_layout()


    @staticmethod
    def maillage(x, y, win=None, shape=None) :
        if shape is None :
            shape=(10, 10)
        X, Y = np.meshgrid(x, y)
        ABS, ORD = [X, X.T], [Y, Y.T]
        COUL, LLB = ['b']*2, [None]*2
        AX, LW = ['x', "y"], [0.8]*2
        t = "Mesh"
        affichage.simple_plot(ABS, ORD, COUL, AX, LLB, lw=LW, shape=shape, titre=t)
        if win is None :
            plt.xlim(X.min(), X.max())
            plt.ylim(Y.min(), Y.max())
        else :
            plt.xlim(win[0], win[1])
            plt.ylim(win[0], win[1])
        plt.grid()
        plt.tight_layout()

    @staticmethod
    def XYlim(x, y) :
        plt.xlim(x.min(), x.max())
        plt.ylim(y.min(), y.max())
        plt.tight_layout()