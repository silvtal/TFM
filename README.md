# TFM
**ES**: Este repositorio recoge los *scripts* desarrollados para automatizar los análisis correspondientes al Trabajo de Fin de Máster "Modelización metabólica de comunidades microbianas estables crecidas con fuentes de carbono y energía únicas y simples" (Silvia Talavera Marcos, Máster en Bioinformática y Biología Computacional, Universidad Autónoma de Madrid). El trabajo puede descargarse [aquí](http://hdl.handle.net/10486/695123).

**EN**: This repository includes scripts from Master's Thesis "Metabolic modeling of microbial communities grown in simple energy and carbon sources" (Silvia Talavera Marcos, Master's Degree in Bioinformatics and Computational Biology, Autonomous University of Madrid). The work is downloadable [here](http://hdl.handle.net/10486/695123).

----------

## Esquema general del _pipeline_ propuesto

![](https://github.com/urihs/TFM/blob/master/Anexo/Esquema.png)

## Descripción de los _scripts_

**modelado.R** incluye el alineamiento con Nucmer, la creación de modelos con CarveMe y (opcionalmente) el análisis con Smetana para una pareja de nodos dada.

**annotate.R** se encarga de la anotación funcional con eggNOG-mapper y, opcionalmente, la creación de un modelo consenso. También es capaz de llamar a Nucmer para iniciar el proceso desde el principio. **create_compatible_database.R** es un _script_ adicional que genera una base de datos compatible con la versión de Diamond que sea utilizada automáticamente por eggNOG-mapper desde _annotate.R_.

Las funciones definidas para estos scripts se encuentran en **utils.R**.

**RefrFBA.py** es un _script_ que se encarga de ejecutar simulaciones de FBA para todos los modelos de un PCG dado. Devuelve las tasas de crecimiento en forma de archivos _.csv_ y por la salida estándar. 

**parser.R** se encarga de generar informes en formato de texto plano que resumen los resultados de Smetana.

Por último, se han desarrollado dos scripts en Python para crear los consensos: **consenso.py** crea un modelo SBML a partir de otros modelos SBML y **consenso_EGG.py** crea una tabla de anotaciones consenso a partir de múltiples archivos de anotaciones.

