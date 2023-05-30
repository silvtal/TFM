# TFM

[![en](https://img.shields.io/badge/lang-en-red.svg)](https://github.com/silvtal/TFM/blob/master/README.md)

Este repositorio recoge los *scripts* desarrollados para automatizar los análisis correspondientes al Trabajo de Fin de Máster "Modelización metabólica de comunidades microbianas estables crecidas con fuentes de carbono y energía únicas y simples" (Silvia Talavera Marcos, Máster en Bioinformática y Biología Computacional, Universidad Autónoma de Madrid). El trabajo puede descargarse [aquí](http://hdl.handle.net/10486/695123).

----------

## Esquema general del _pipeline_ propuesto

![](https://github.com/silvtal/TFM/blob/master/Anexo/Esquema.png)

## Descripción de los _scripts_

**modelado.R** incluye el alineamiento con Nucmer, la creación de modelos con CarveMe y (opcionalmente) el análisis con Smetana para una pareja de nodos dada. Llama al resto de scripts, pero tiene hard-coded los parámetros.

**annotate.R** se encarga de la anotación funcional con eggNOG-mapper y, opcionalmente, la creación de un modelo consenso. También es capaz de llamar a Nucmer para iniciar el proceso desde el principio. **create_compatible_database.R** es un _script_ adicional que capaz de generar una base de datos compatible con la versión de Diamond que sea utilizada automáticamente por eggNOG-mapper desde _annotate.R_.

Las funciones definidas para estos scripts se encuentran en **utils.R**.

**RefrFBA.py** es un _script_ que se encarga de ejecutar análisis de balance de flujo (FBA) para todos los modelos de un grupo filogenético central (PCG) dado. Devuelve las tasas de crecimiento tanto en forma de archivos _.csv_ como por la salida estándar. 

**parser.R** se encarga de generar informes en formato de texto plano que resumen los resultados de Smetana.

Por último, se han desarrollado dos scripts en Python para crear los consensos (tanto anotaciones como modelos consenso): **consenso.py** crea un modelo SBML a partir de otros modelos SBML y **consenso_EGG.py** crea una tabla de anotaciones consenso a partir de múltiples archivos de anotaciones.

## Input

La descripción del formato de los ficheros de entrada está desactualizada en el Anexo: ahora no se necesita un árbol filogenético, y además el fichero con la información de los nodos de interés es diferente (incluye una lista de las hojas de cada nodo). Se incluyen en este repositorio ficheros de ejemplo, además de un script de prueba ("test")

## Troubleshooting

- Si se instalan CarveMe/Python3.7 en un entorno de Conda, se puede activar este entorno desde el propio script. Por ejemplo, la llamada de la función `carve()` puede modificarse a: `conda activate carveme && carve ",db_protein_folder`(...).

- Si Smetana y/o CarveMe no detectan CPLEX o la versión correcta de CPLEX, hay que [definir la variable global `PYTHONPATH`](https://www.ibm.com/docs/en/icos/20.1.0?topic=cplex-setting-up-python-api) correctamente 
