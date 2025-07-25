# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 15:15:27 2025

@author: Paul Brancher

New colormaps !
 - 2 beautiful coolwarm with orange
 - More than 100 option diverging from black 
 - The same but from white (Ugly)
 - Some variations from black to another color
 - Some diverging as black / color light / color dark
 - Some macros about cmap
 - Some white diverging cmap
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, to_rgb
# import seaborn as sns
from Affichage import affichage


class cmpw(object):
    
    @staticmethod
    def printCUST():
        """ display of all diverging black cmap """
        x = np.linspace(0, 4*np.pi, 200)
        X, Y = np.meshgrid(x, x)
        Z = np.cos(X) + np.sin(Y) + 0.5*X - 0.5*Y

        color=['red', 'violet', 'purple', 'blue', 'green', 'orange', 'steel', 'mint', 'cream', 'salmon']
        fig, axs = plt.subplots(len(color), len(color), figsize=(40, 40))

        for i in range(len(color)):
            for j in range(len(color)):
                axs[i, j].pcolor(X, Y, Z, shading='auto', cmap=cmpw.neonCUSTOM(color[i], color[j]))
                axs[i, j].tick_params(axis='both', which='major', labelsize=0)
        plt.tight_layout()
        
                

    @staticmethod
    def cmptest(cmap, titre=None):
        """ display of a chosen cmap """
        x = np.linspace(0, 4*np.pi, 200)
        X, Y = np.meshgrid(x, x)
        Z = np.cos(X) + np.sin(Y) + 0.5*X - 0.5*Y
        affichage.graph2D(X, Y, Z, 0, cmap=cmap, cbar='ok', axes='zero', nbplot=4, cs='%.2f', titre=titre)
        plt.tight_layout()
    
    # @staticmethod
    # def icefire():
    #     return sns.color_palette("icefire", as_cmap=True)
    
    # @staticmethod
    # def vlag():
    #     return sns.color_palette("vlag", as_cmap=True)
    
    @staticmethod
    def orangeblue():
        points = [
                (0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275),
                (0.030334, 0.113725, 0.0235294, 0.45098),
                (0.055527, 0.105882, 0.0509804, 0.509804),
                (0.073008, 0.0392157, 0.0392157, 0.560784),
                (0.089974, 0.0313725, 0.0980392, 0.6),
                (0.106427, 0.0431373, 0.164706, 0.639216),
                (0.130077, 0.054902, 0.243137, 0.678431),
                (0.16144, 0.054902, 0.317647, 0.709804),
                (0.2, 0.0509804, 0.396078, 0.741176),
                (0.225, 0.0392157, 0.466667, 0.768627),
                (0.25, 0.0313725, 0.537255, 0.788235),
                (0.276093, 0.0313725, 0.615686, 0.811765),
                (0.302828, 0.0235294, 0.709804, 0.831373),
                (0.329563, 0.0509804, 0.8, 0.85098),
                (0.351671, 0.0705882, 0.854902, 0.870588),
                (0.372237, 0.262745, 0.901961, 0.862745),
                (0.390231, 0.423529, 0.941176, 0.87451),
                (0.417995, 0.572549, 0.964706, 0.835294),
                (0.436504, 0.658824, 0.980392, 0.843137),
                (0.456041, 0.764706, 0.980392, 0.866667),
                (0.468895, 0.827451, 0.980392, 0.886275),
                (0.482262, 0.890196, 0.988235, 0.92549),
                (0.492545, 0.913725, 0.988235, 0.937255),
                (0.501285, 1, 1, 0.972549),
                (0.510026, 0.988235, 0.988235, 0.905882),
                (0.526478, 0.992157, 0.972549, 0.803921),
                (0.539846, 0.992157, 0.964706, 0.713725),
                (0.554756, 0.988235, 0.956863, 0.643137),
                (0.576864, 0.980392, 0.917647, 0.509804),
                (0.599486, 0.968627, 0.87451, 0.407843),
                (0.620051, 0.94902, 0.823529, 0.321569),
                (0.636504, 0.929412, 0.776471, 0.278431),
                (0.660668, 0.909804, 0.717647, 0.235294),
                (0.682262, 0.890196, 0.658824, 0.196078),
                (0.7, 0.878431, 0.619608, 0.168627),
                (0.725, 0.870588, 0.54902, 0.156863),
                (0.75, 0.85098, 0.47451, 0.145098),
                (0.775, 0.831373, 0.411765, 0.133333),
                (0.8, 0.811765, 0.345098, 0.113725),
                (0.825, 0.788235, 0.266667, 0.0941176),
                (0.85, 0.741176, 0.184314, 0.0745098),
                (0.875, 0.690196, 0.12549, 0.0627451),
                (0.9, 0.619608, 0.0627451, 0.0431373),
                (0.923393, 0.54902, 0.027451, 0.0705882),
                (0.943959, 0.470588, 0.0156863, 0.0901961),
                (0.967095, 0.4, 0.00392157, 0.101961),
                (1, 0.188235, 0, 0.0705882)
                ]
        
        positions = [p[0] for p in points]
        colors = [p[1:] for p in points]

        custom_cmap = LinearSegmentedColormap.from_list("Divergent_1", list(zip(positions, colors)))
        return custom_cmap
    
    @staticmethod
    def orangeblue_light():
        points = [
            (0, 0.03010016365396975, 0.08341997087738243, 0.43000000542450717),
            (0.16144, 0.14905897070127497, 0.5322347084650014, 0.7098040101429154),
            (0.351671, 0.4140000971097244, 0.9000000186376352, 0.8351999703415328),
            (0.501285, 0.9874999531995236, 1, 0.8499999826260306),
            (0.620051, 0.9490199368674209, 0.8235290793622709, 0.32156893407857895),
            (0.8, 0.8117648716952724, 0.36367082948247303, 0.06494106256058779),
            (1, 0.47999995904609954, 0.06096038788114378, 0.014399937913024172)
            ]
        
        positions = [p[0] for p in points]
        colors = [p[1:] for p in points]

        custom_cmap = LinearSegmentedColormap.from_list("Divergent_1", list(zip(positions, colors)))
        return custom_cmap


    @staticmethod
    def neonCUSTOM(chaud, froid):
        
        dico = dict(red=['red', 'lightcoral'],
                    violet=['violet', 'pink'],
                    purple=['blueviolet', 'pink'],
                    blue=['cyan', 'lavender'],
                    green=['lime', 'yellow'],
                    orange=['orange', 'yellow'],
                    steel=['steelblue', 'powderblue'],
                    cream=["burlywood", 'papayawhip'],
                    mint=["springgreen", "mintcream"],
                    salmon=['lightsalmon', 'antiquewhite'])
        
        chaud, froid = dico[chaud], dico[froid]
        custom_colors = froid[::-1] + ['k'] + chaud
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap
    
    @staticmethod
    def divergeCUSTOM(chaud, froid):
        
        dico = dict(red=['red', 'lightcoral'],
                    violet=['violet', 'pink'],
                    purple=['blueviolet', 'pink'],
                    blue=['cyan', 'lavender'],
                    green=['lime', 'yellow'],
                    orange=['orange', 'yellow'],
                    steel=['steelblue', 'powderblue'],
                    cream=["burlywood", 'papayawhip'],
                    mint=["springgreen", "mintcream"],
                    salmon=['lightsalmon', 'antiquewhite'])
        
        chaud=dico[chaud]
        froid=dico[froid]
        
        custom_colors = froid + ['white'] + chaud[::-1]
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap



    @staticmethod
    def DARK2Blue():
        custom_colors = ['black', 'navy', 'blue', 'deepskyblue', 'cyan', 'aqua']
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap
    
    @staticmethod
    def DARK2Red():
        custom_colors = ['black', 'maroon', 'red', 'orangered', 'gold']
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap
    
    @staticmethod
    def DARK2Green():
        custom_colors = ['black', 'forestgreen', 'lime', 'gold']
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap
    
    @staticmethod
    def DARK2Blue2():
        custom_colors = ['black', 'blue', 'cyan', 'gold', 'yellow']
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap
    
    
    @staticmethod
    def lineDARK():
        custom_colors = ['navy','blue', 'deepskyblue', 'black', 'lightcoral', 'red', 'darkred']
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap
    
    @staticmethod
    def lineDARK2():
        chaud = ['blue', 'cyan']
        froid = ['forestgreen', 'lime']
        custom_colors = froid + ['k'] + chaud[::-1]
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        return custom_cmap

    
    @staticmethod
    def divergeW():
        chaud = [ 'darkred', 'red', 'tomato']
        froid = ['dodgerblue', 'blue',  'navy']
        custom_colors = froid[::-1] + ['seashell']+ chaud[::-1]
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        
        return custom_cmap
    
    @staticmethod
    def adjust_brightness(color, factor):
        r, g, b = to_rgb(color)
        r = min(r * factor, 1.0)
        g = min(g * factor, 1.0)
        b = min(b * factor, 1.0)
        return (r, g, b)

    @staticmethod
    def stretched_colors(cmap_name='coolwarm', n=256, stretch=0.7):
        base = plt.get_cmap(cmap_name)
        stretched = np.linspace(0, 1, n)
        stretched = 0.5 * (1 + np.tanh((stretched - 0.5) / stretch))  # sigmoïde aplatie
        colors = [base(val) for val in stretched]
        return ListedColormap(colors, name=f"stretched_{cmap_name}")
    
    @staticmethod
    def divergeW2():
        chaud = ['red', 'orange']
        froid = ['blue', 'lightskyblue']
        custom_colors = froid[::-1] + ['#ffebff'] + chaud
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        
        return custom_cmap
    
    @staticmethod
    def divergeW3():
        chaud = ['red','#e70000', '#ff8c00']
        froid = ['lightskyblue', '#0044ff', 'blue']
        
        chaud = ['firebrick','red','#e70000', '#ff8c00']
        froid = ['lightskyblue', '#0044ff', 'blue', 'mediumblue']
        custom_colors = froid[::-1] + ['w'] + chaud[::-1]
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        
        return custom_cmap
    
    @staticmethod
    def divergeW4():
        chaud = ['gold', 'orange', 'darkorange', 'red', 'firebrick', 'darkred']
        froid = ['aqua', 'deepskyblue', 'dodgerblue', 'blue', 'mediumblue', 'navy']
        custom_colors = froid[::-1] + ['w'] + chaud
        custom_cmap = LinearSegmentedColormap.from_list("custom_gradient", custom_colors)
        
        return custom_cmap
