# -*- coding: utf-8 -*-
"""
Created on Wed May 28 16:35:15 2025

@author: un génie du code
"""

import numpy as np
from lxml import etree

class XDMF_reader(object):

    @staticmethod
    def read_binary_dataitems(xdmf_file):
        from lxml import etree
        import numpy as np

        with open(xdmf_file, 'rb') as f:
            tree = etree.parse(f)

        binary_items = []
        for data in tree.xpath('//DataItem'):
            format_ = data.attrib.get("Format", "").lower()
            if format_ == "binary":
                dims = tuple(map(int, data.attrib.get("Dimensions").split()))
                dtype = data.attrib.get("NumberType", "Float").lower()
                precision = int(data.attrib.get("Precision", "4"))
                file_name = data.text.strip()

                if dtype == "float":
                    np_dtype = np.float32 if precision == 4 else np.float64
                elif dtype == "int":
                    np_dtype = np.int32 if precision == 4 else np.int64
                else:
                    continue  # skip unsupported

                binary_items.append({
                    'file': file_name,
                    'shape': dims,
                    'dtype': np_dtype
                })

        return binary_items

    @staticmethod
    def load_field(info, base_path=""):
        path = base_path + info['file']
        data = np.fromfile(path, dtype=info['dtype'])
        return data.reshape(info['shape'])

    @staticmethod
    def load_data(path, basepath="", verbose=False):
        path = basepath + path
        items = XDMF_reader.read_binary_dataitems(path)


        ux = XDMF_reader.load_field(items[0], base_path=basepath)
        ux=ux[0,:]
        
        uy = XDMF_reader.load_field(items[1], base_path=basepath)
        uy=uy[0,:]
        
        pp = XDMF_reader.load_field(items[3], base_path=basepath)
        pp=pp[0,:]

        phi = XDMF_reader.load_field(items[4], base_path=basepath)
        phi=phi[0,:]
        
        if verbose :
            nbLINE=45
            print('='*(nbLINE+1))
            print("- Lecture de :", path)
            print('-'*nbLINE)
            for i, item in enumerate(items):
                print(f"- {i}: {item['file']} - shape={item['shape']}")
            print('-'*nbLINE)
            if np.size(pp)==np.size(ux) :
                print("- Taille d'un champ", np.shape(pp))
                print("- Nombre de points", np.size(pp))
            print('='*(nbLINE+1))

        return ux, uy, pp, phi
    
    @staticmethod
    def timed_data(repertoire, snapMAX=100, snapMIN=0) :
        PP=[]
        UX=[]
        UY=[]
        PHI=[]
        
        for time in range(snapMIN, snapMAX):
            time+=1
            path=f'snapshot-{time}.xdmf'
            ux, uy, pp, phi = XDMF_reader.load_data(path, basepath=repertoire)
            UX.append(ux)
            UY.append(uy)
            PP.append(pp)
            PHI.append(phi)
            
            if time==snapMIN+1 :
                print(repertoire+" ----------------------")
                print("- Premier lu :", path)
            elif time==snapMAX :
                print("- Dernier lu :", path)
        
        return UX, UY, PP, PHI
    
    @staticmethod
    def lecture(path, delimiter = None):
        if delimiter is None :
            delimiter = ","
        """ Read a file - retrun its columns as lists. """
        tableau = np.loadtxt(path, delimiter=delimiter, skiprows=0)
        liste_colonnes = []  # mets les colonnes dans une liste
        for i in range(tableau.shape[1]):
            liste_colonnes.append(tableau[:, i])
        return tuple(liste_colonnes)
    
    @staticmethod
    def creaFichier(path, *colonnes, header=None):
        data_table = np.column_stack(colonnes)
        if  header != None:
            np.savetxt(path, data_table, delimiter=',', fmt="%.6f",
                       header=header)
        else:
            np.savetxt(path, data_table, delimiter=',', fmt="%.6f")