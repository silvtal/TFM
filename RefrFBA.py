#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 12:24:06 2020

@author: urihs
"""

from reframed import load_cbmodel
import os
from reframed import FBA, Environment
from carveme.reconstruction.utils import load_media_db
import statistics as stats


wd = "/home/urihs/Desktop/TFM_private/08_reframed/sin_gapfill/"


models = {}
for medium in os.listdir(wd):
    models[medium] = {}
    for node in os.listdir(wd+medium):
        if "Node" in node:
            models[medium][node] = {}
        else:
            continue
        for sbml in os.listdir(wd+medium+"/"+node):
            if ".xml" in sbml:
                models[medium][node][sbml]=wd+medium+"/"+node+"/"+sbml

media_db = load_media_db("/home/urihs/Desktop/TFM_private/08_reframed/test_media.txt", compound_col="compound")

for medium in models.keys():
    medio = "M9["+medium+"]"
    # medio = "LB"
    init_env = Environment.from_compounds(media_db[medio])
    print("\nMedio: "+medio+"\nNodos de: "+medium+"\n==============")
    media = 0 # !!!
    listofkeys = list(models[medium].keys())
    listofkeys.sort()
    for node in listofkeys:
        media_prev = 0+media
        print("Nodo: "+node+"\n===============")
        total=[]
        n=0# para hacer la media de solution.fobj
        for f in models[medium][node]:
            model = load_cbmodel(models[medium][node][f],flavor="fbc2")
            init_env.apply(model)
            solution = FBA(model,   objective="Growth", get_values=False)    
            print(f)
            print(solution)
            if solution.fobj!=0:
                n+=1
                total.append(solution.fobj)
        try:
            media = sum(total)/n
            print("media del nodo: ", media)
            print("desviacion típica: ", stats.stdev(total),"\n")
        except:
            pass
    try:
        print("RATIO : ",media_prev/media) # es E/P para leu y P/E para cit y glc
    except:
        pass
