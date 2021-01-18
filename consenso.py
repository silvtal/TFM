#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:42:54 2020

@author: urihs
"""

import os
import sys
import xml.etree.ElementTree as ET
import time, datetime

start = time.time() # para devolver al final el tiempo de ejecución
# =============================================================================
# Leer todos los modelos
# =============================================================================
# INPUT: consenso.py input_folder output_folder outputname
wd = sys.argv[1] # !!! ponerle un input más guay (con flags)
os.chdir(wd)

if len(sys.argv) > 2:
    outputdir = sys.argv[2].rstrip("/") # si el usuario da el outputdir
    
    if len(sys.argv) > 3:
        nodename = sys.argv[3] # si el usuario da el nombre del output
    else:
        nodename = wd.split("/")[-1]
else:
    outputdir = wd

models = [] # guardamos aquí los nombres de los archivos
for filename in os.listdir(wd):
    if ".xml" in filename:
        models.append(filename)

# !!! necesito mirar campos distintos si es gapfilled. La reacción de crecimiento 
        # también es distinta. Por ahora NO SIRVE SI ES GAPFILLED.
# =============================================================================
# Preparar el formato SBML
# =============================================================================
ET.register_namespace('', "http://www.sbml.org/sbml/level3/version1/core")
ET.register_namespace('fbc',"http://www.sbml.org/sbml/level3/version1/fbc/version2")
# !!! quizá estaría bien generalizar esto

# =============================================================================
# Crear el "esqueleto" de nuestro modelo consenso. Usamos un modelo cualquiera.
# =============================================================================
try:
    first_model = ET.parse(models[0])
    consensus   = first_model
    cons_root   = consensus.getroot()
    cons_root[0].attrib['id'] = nodename
    cons_root[0][0][0][0].text = 'Description: This model is a consensus model from SBML FBC2 CarveMe output'
except:
    print("No se ha podido cargar ningún modelo.")
    raise(SystemExit(0))

# =============================================================================
# for child in root:
#     for grandchild in child:
#         print(grandchild.tag, grandchild.attrib)
#
#
# # {http://www.sbml.org/sbml/level3/version1/core}notes {}
#         # este es root[0][0] 
# # {http://www.sbml.org/sbml/level3/version1/core}listOfCompartments {}
#         #este es root[0][1]
# # {http://www.sbml.org/sbml/level3/version1/core}listOfSpecies {}
# # {http://www.sbml.org/sbml/level3/version1/core}listOfParameters {}
# # {http://www.sbml.org/sbml/level3/version1/core}listOfReactions {}
# # {http://www.sbml.org/sbml/level3/version1/fbc/version2}listOfObjectives {'{http://www.sbml.org/sbml/level3/version1/fbc/version2}activeObjective': 'objective'}
# # {http://www.sbml.org/sbml/level3/version1/fbc/version2}listOfGeneProducts {}
# 
# Nos quedamos con los campos 0 (modificado), 1, 3 y 5 del first_model
# El 6 sí lo cambiamos porque hay algunas reacciones que incluyen estas 
# anotaciones. Así pues, igual que las reacciones nos hacen seleccionar 
# reactivos, también nos hacen seleccionar geneproducts. (https://www.researchgate.net/publication/292228111_SBML_Level_3_Package_Flux_Balance_Constraints_'fbc')
# Después, tenemos que ir añadiendo Y seleccionando 2, 4 y 6
# =============================================================================
# ??? preguntarle a daniel si opina lo mismo del campo 6
    

# =============================================================================
# Hacer un conteo de cuántas veces aparece cada elemento en todos los
# modelos. *Evitamos cargar todos en memoria a la vez.*
# Nos quedaremos con aquellos que aparecen en el 80% o más de los modelos
# (4 de cada 5)
# =============================================================================

# Inauguramos una lista de reacciones de crecimiento (a partir de la cual obtendremos la reacción de biomasa consenso)
growth = [first_model.getroot()[0][4][-2]]
cons_root[0][4].remove (cons_root[0][4][-2]) # y la eliminamos

# Hago una lista de diccionarios donde se hará el conteo
conteo = [{},{},{}]

# Los ID de cada elemento tienen diferente nombre según el campo
ID = ['id', 'id', '{http://www.sbml.org/sbml/level3/version1/fbc/version2}id']


# Inauguramos esta lista con el primer modelo
for n, campo in enumerate([2, 4, 6]):
    conteo[n] = {cons_root[0][campo][e].attrib[ID[n]]:1 for e in range(len(cons_root[0][campo]))}
    # <dicci>   <------------key-------------------><1>
del(first_model)

# Vamos abriendo los demás modelos y contando
# Guardamos los elementos de los campos 2 y 4 que no estaban ya
for m in models[1:]:
    model = ET.parse(m)
    # ??? iterparse? mirar y tal
    model_root = model.getroot()
    # dejamos aparte la función de crecimiento
    # ??? hacer esto más bonito?
    growth.append(model_root[0][4][-2])
    model_root[0][4].remove (model_root[0][4][-2])

    for n, campo in enumerate([2, 4, 6]):
        for e in range(len(model_root[0][campo])):
            ident = model_root[0][campo][e].attrib[ID[n]]
            if ident not in conteo[n].keys():
                conteo[n][ident] = 1
                      #dic  #<-------key------------------------->
                cons_root[0][campo].append(model_root[0][campo][e]) 
            else:
                conteo[n][ident] += 1
                
    del(model) # !!! haciendo esto borro en memoria
    del(model_root)

# !!! probar varios threshold (90%, 80%)
# Filtro y elimino las entradas que están en menos de un 80% de modelos:
total = len(models)
for n, campo in enumerate([2,4,6]):
    removed = 0 # actualizamos el índice para evitar errores de indexación tras borrar
    for e in range(len(cons_root[0][campo])):
        ident = cons_root[0][campo][e-removed].attrib[ID[n]]
        if conteo[n][ident] < 0.80*total:
            cons_root[0][campo].remove (cons_root[0][campo][e-removed])
            removed+=1
    
## Ejemplo:
# =============================================================================
# len(old_cons[0][4])
# Out[312]: 3458
# 
# len(cons_root[0][4])
# Out[313]: 2149
# =============================================================================
# !!! raise: "solo sobreviven X", comparar con la media de reacciones...

# =============================================================================
# Por último añadimos la reacción de crecimiento
# =============================================================================

pure_list = [item for item in growth]
consensus_growth = max(set(pure_list), key=pure_list.count)
# !!! warning: no va a pasar, pero meter alguno para por si growth coge en algún 
    #momento un growth que no es growth (como cuando coge consensus.xml) y que avise que archivo es.
# !!! otro warning chulo: decir cuantos outliers hay, de alguna forma!
cons_root[0][4].append(consensus_growth)

# Guardamos el modelo
consensus.write(outputdir+"/"+nodename+"_consensus.xml", xml_declaration=True)

end = time.time()
print("Running time: ",str(datetime.timedelta(seconds = end - start)))
