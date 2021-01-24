#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 12:24:06 2020

@author: urihs
"""

from reframed import load_cbmodel
import os, sys
from reframed import FBA, Environment
from carveme.reconstruction.utils import load_media_db
import statistics as stats
import pandas as pd

# INPUT: RefrFBA_E2.py input_folder/wd mediadb medium outputdir outputname

# IMPRIME la media y la desviación típica del crecimiento de todos los 
#         modelos de la carpeta dada como entrada
# GENERA un archivo .csv con todas las tasas de crecimiento

wd = sys.argv[1] 
mediadb = sys.argv[2]
medium = sys.argv[3]
    
if len(sys.argv) > 4:
    outputdir = sys.argv[4].rstrip("/") # si el usuario da el outputdir
    
    if len(sys.argv) > 5:
        outputname = sys.argv[5] # si el usuario da el nombre del output
    else:
        outputname = "report_reframed"+wd.split("/")[-1]+".csv"
else:
    outputdir = wd


models = []
for sbml in os.listdir(wd):
  if ".xml" in sbml:
    models.append(wd+"/"+sbml)

media_db = load_media_db(mediadb, compound_col="compound")

results = [] # lista de listas que finalmente dará un dataframe/csv con todos los datos

init_env = Environment.from_compounds(media_db[medium])
print("\nMedio: "+medium+"\n==============")
total=[]
n=0# para hacer la media de solution.fobj
for f in models:
    model = load_cbmodel(f,flavor="fbc2")
    init_env.apply(model)
    solution = FBA(model,   objective="Growth", get_values=False)    
    print(solution)
    if solution.fobj!=0:
        n+=1
        total.append(solution.fobj)
    results.append([medium,f.split("/")[-1],solution.fobj])
try:
    media = sum(total)/n
    print("media del nodo: ", media)
    print("desviacion típica: ", stats.stdev(total),"\n")
except:
    pass

all_data=pd.DataFrame(results,columns=["Medio","Modelo","Crecimiento"])
all_data.to_csv(outputdir+outputname,index=False)
