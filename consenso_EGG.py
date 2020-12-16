#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 11:18:55 2020

@author: urihs
"""

import os
import sys
from carveme.reconstruction.eggnog import load_eggnog_data
import datetime, time

start = time.time() # para devolver al final el tiempo de ejecución
# =============================================================================
# Leer todos los modelos
# =============================================================================
# INPUT: consenso_EGG.py input_folder output_folder (outputname) (percentage)
# SOLO LEE ARCHIVOS QUE CONTENGAN "annotations" EN EL NOMBRE
# TIENE QUE USARSE EL "REALPATH"
wd = sys.argv[1].rstrip("/") # !!! ponerle un input más guay (con flags)

if len(sys.argv) > 2:
    outputdir = sys.argv[2].rstrip("/") # si el usuario da el outputdir
    
    if len(sys.argv) > 3:
        outputname = sys.argv[3].rstrip("/") # si el usuario da el nombre del output
        if len(sys.argv) > 4: 
            perc = float(sys.argv[4]) # si el usuario da el porcentaje para filtrar.
        else:
            perc = 0.80
    else:
        outputname = wd.split("/")[-1]
else:
    outputdir = wd

models = [] # guardamos aquí los nombres de los archivos
for filename in os.listdir(wd):
    if "annotations" in filename: # ignore seed orthologs files
        models.append(wd+"/"+filename)
    # !!! estaría bien un warning AL ABRIR el modelo, si 
    # no se puede abrir o si no tiene todas las columnas !

if len(models) == 0:
    quit("The input file is empty.")
else:
    num_of_models = len(models)


# =============================================================================
# Crear un almacén de todas las reacciones a partir de un modelo cualquiera.
# =============================================================================
all_reactions = load_eggnog_data(models[0],drop_unused_cols=False)
    # !!! va MUCHO más lento si uso las 17 columnas (drop_unused_cols=False). poniendo 
    # !!! True funciona exactamente igual pero más rapido, la unica diferencia es que el 
    # !!! output no tiene todas las columnas. podria ajustar esto para que trabaje con 4 
    # !!! solo y luego genere una matriz sparse con las 17

# Agrupamos por reacción para evitar repeticiones:
all_reactions = all_reactions.sort_values(by='score', ascending=False) \
                             .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
# Índice para acceder eficientemente a la fila de cada reacción:
reac2idx = {reac:idx for idx,reac in enumerate(all_reactions["BiGG_gene"])}
last_index = len(all_reactions)-1

# =============================================================================
#  Voy abriendo los demás y contando las apariciones de cada reacción
# =============================================================================

# Inauguramos el diccionario de conteo con las reacciones del primer modelo
conteo = {reac:1 for reac in all_reactions["BiGG_gene"]}
for model in models[1:]:
    # abrimos el modelo
    open_model = load_eggnog_data(model,drop_unused_cols=False)
    open_model = open_model.sort_values(by='score', ascending=False) \
                          .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])

    for i, df_row in open_model.iterrows():
        reaction = df_row["BiGG_gene"]
        
        if reaction not in conteo.keys():
            conteo[reaction] = 1         # la contamos
            all_reactions = all_reactions.append(df_row,ignore_index=True) # la añadimos al almacén
            last_index += 1              # y anotamos su índice
            reac2idx[reaction] = last_index
            
        
        else:
            conteo[reaction] += 1        # la contamos
            # puntuación nueva ponderada (el resultado es la media para esa reacción):
            old_score = all_reactions.at[reac2idx[reaction],"score"]
            new_score = ( old_score *(conteo[reaction]-1) + df_row["score"] ) /conteo[reaction]
            all_reactions.at[reac2idx[reaction],"score"] = new_score # actualizar
    # cerramos el modelo (lo quitamos de memoria)
    del(open_model)
    
    
# =============================================================================
# Filtro y elimino las reacciones que están en menos de un 80% de modelos:
# =============================================================================
for r in reac2idx.keys():
    if conteo[r] < int(perc*num_of_models):
        all_reactions = all_reactions.drop(reac2idx[r],axis=0)
# =============================================================================
# Guardo el archivo de anotaciones consenso
# =============================================================================
if not os.path.exists(outputdir):
    os.mkdir(outputdir)
f = open(outputdir+"/"+outputname+".tsv","a+")

f.write("# emapper version: emapper-2.0.1 emapper DB: 2.0\n")
f.write("# consensus annotation file for "+outputname+"\n")
f.write("# time: "+str(datetime.datetime.now())+"\n")
f.write("#query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	best_tax_level	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TCCAZy	BiGG_Reaction	taxonomic scope	eggNOG OGs	best eggNOG OG	COG Functional cat.	eggNOG free text desc.\n")
all_reactions.to_csv(f,sep="\t",header=False, index=False)

f.close()

end = time.time()
print("Running time: ",end - start)

