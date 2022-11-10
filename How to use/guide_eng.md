Description of the creation process and analysis of CNT systems that present different degrees of aggregation inside a homogeneous medium.


![Sistemas de cadenas con distintos grados de dispersión generados mediante esta metodología ](https://i.ibb.co/F56rJ5F/figura-markdown.png)



## Index
* [Systems construction](#systems-construction)
  * [Prole chain](#prole-chain)
  * [Homogeneous chains](#homogeneous-chains)
  * [Arrangement (Ordenador)](#arrangement-ordenador)
  * [Equilibrium state](#equilibrium-state)
  * [TakeConf](#takeconf)
* [Post-processing and analysis](#post-processing-and-analysis)
  * [Make up (Maquillador)](#make-up-maquillador)
  * [Pair correlation function and scattering intensity](#pair-correlation-function-and-scattering-intensity)
  * [External radius and radius of gyration](#external-radius-and-radius-of-gyration)

<div id='id-section1'/>

## Systems construction ## 
<div id='id-section1-1'/>

### Prole chain ###
The code <code>CadenaProle.f90</code> generates an isolated chain by concatenating spheres in a hard contact condition. During the process of creating the chain, the spheres are concatenated successively, checking that there is no overlap between them before accepting each new sphere. The parameters that can be set for string generation are:
- **nPartCadena**: Number of spheres that form the chain(*line 40*).
- **theta**: Angle to place the next sphere in the XY plane (*line 76*).
- **phi**: Angle to place the next sphere with respect to the Z axis (*line 78*).

The parameter **nPartCadena** limits the length of the chain, while **theta** and **phi** define the curvature. By default, these parameters take the values **nPartCadena = 1000**, **theta = pi/4** y **phi = pi/4**.

```fortran
nPartCadena=1000 !Numero de esferas para cada cadena

...

theta=theta*pi/4	!Angulo plano XY
phi=phi*pi/4		!Angulo eje Z
```

If the values are set as **theta = 0** y **phi = 0**, a rigid chain is generated.

The results are stored in the file <code> CadenaProle.dat</code>.

<div id='id-section1-2'/>

### Homogeneous chains ###
The code <code>CadenasHomogeneas_v3.f90</code> generates a random and uniform distribution of strings in a cubic simulation box. The input parameters to generate the system are:
- **Filename.dat**: Name of the file which contains the data of the chain that is going to be used as a model of the CNTs. This file must contain the positions of the spheres in format *OVITO*: a) the first line of the file is an integer that matches the number of spheres that make up the string; b) the second line contains a short description of the system (comment line); c) the following lines contain (by rows) the index of the sphere, the coordinates *(x,y,z)* and the radius.
- **Lcaja**: Edge size of the cubic simulation box .
- **nCad**: Number of chains that are placed in the simulation box.

These parameters are automatically read from the <code>Parametros.inp</code> file when executing the code. By default, they take the values:
```txt
CadenaProle.dat		!Nombre del archivo que contiene la cadena modelo (.dat) formato OVITO
750.0d0			!Arista de la caja total dentro de la cual se generan los CM de las cadenas
500			!Numero de cadenas total del sistema
```

The algorithm builds a cubic simulation box with edge **Lbox** and generates a total of **nCad** random points uniformly distributed inside it. These random points represent the position of the centers of mass of the chains that model the CNTs. Also, each chain is randomly rotated before being added to the system. When all the strings have been placed, the position of all the spheres in the system is readjusted by applying periodic boundary conditions (PBC).

The chain configuration is saved in *OVITO* format in the file <code>SistCadenasHomogeneas_v3.dat</code>.

<div id='id-section1-3'/>

### Arrangement (Ordenador) ###
The algorithm <code>Ordenador_v0.97.f90</code> performs a Monte Carlo process in which the chains progressively aggregate around different centers of attraction to give rise to the formation of aggregates. The input parameters are:
- **Filename.dat**: Name of the file that contains the initial system configuration data. This file must contain the positions of the spheres in *OVITO* format.
- **NEsf_Cad**: Number of spheres that make up each chain.
- **Ncentros**:  Number of attracting centers generated to form the chain aggregates.
- **Niter**: Maximum number of iterations to relax the system.
- **dr**: Size of the random step that the strings move in each iteration.
- **r0**: Radius of the spheres.
- **NPelicula**: Iteration interval to save the system configuration.

These parameters are automatically read from the file <code>ordenador.in</code> when the code is executed. By default, they take the values:
```txt
SistCadenasHomogeneas_v3.dat	!Fichero con la configuración inicial de las cadenas
20				!Numero de esferas en cada cadena
1				!Numero de centros atractores
10000				!Numero de iteraciones maximas para relajar el sistema. 
0.1				!Tamanio del salto aleatorio. Como R=0.5, pueeees.... dr=0.1?
0.5				!Radio de las esferas
100				!Intervalo de iteraciones para grabar la configuracion
```

After loading the initial configuration, the algorithm generates **Ncenters** random points inside the simulation box that act as attracting centers and calculates, for each chain in the system, which is its closest attracting center in the initial configuration. In addition, it calculates the initial energy of the system, defined as the sum of the distances from the centers of mass to their respective centers of attraction. In each iteration, the centers of mass of the chains are randomly displaced by a distance **dr**. The movement is only accepted if the center of mass approaches its corresponding center of attraction. This process is repeated until completing the **Niter** iterations, calculating the energy of the system in each of them.

At the end of the aggregation process, the files are generated:

- **Registro.out**: It contains the information of the parameters used in the simulation and a list that enumerates the attractor center to which each chain approaches.
- **CM_cadenas.dat**: Contains the coordinates of the mass centers of the chains in the initial configuration, in *OVITO* format.
-  **Centros.dat**: Contains the coordinates of the attracting centers, in *OVITO* format.
- **Registro_cad_virtual.dat**: Contains the coordinates of the spheres in the initial configuration, but eliminating the periodic boundary conditions that are applied at the end of the algorithm <code>CadenasHomogeneas_v3.f90</code>, so that the chains on the border of the simulation box are not cut off. This format is useful to visualize the system more clearly.
- **Energia.dat**: Contains the value of the energy of the system in each iteration of the Monte Carlo process.
- **Peli.dat**: Contains the system configuration every **NPelicula** iterations, in *OVITO* format.

<div id='id-section1-4'/>

### Equilibrium state ###
The algorithm <code>Nstep_Max.f95</code>  takes as input the file <code>Energia.dat</code> to calculate the Monte Carlo iteration in which the equilibrium configuration is reached. This is defined as the configuration in which the energy varies less than 1% with respect to the previous configuration. The result is stored in the file <code> EnergiaNstepmax.dat</code>.

Since the system configuration throughout the aggregation process is stored every **NPelicula** iterations, the maximum aggregation configuration is chosen as the one that appears in the **Peli.dat** file and is greater than the iteration defined in <code> EnergiaNstepmax.dat</code>.

<div id='id-section1-5'/>

### TakeConf ###
The file <code>TakeConf.sh</code> allows you to extract a configuration from the movie generated in the aggregation process. To do this, use as input parameters:

- **Configuración**: Number of the configuration that you want to extract. Since the configurations are stored every **NPelicula** iterations, choosing the *N* configuration implies selecting the Monte Carlo iteration *Iter = N x NPelicula*.
- **Esferas_conf**: Number of spheres that make up each configuration.
- **File**: Name of the file that stores the configurations during the aggregation process.

By default, these parameters take the values **Configuración=10**, **Esferas_conf=500000** y **File=Peli.dat**:
```bash
#!/bin/bash
# ----------------------------------
# Sacar una configuracion de la Pelicula de Ovito
# ----------------------------------

# Configuracion y caracteristicas del sistema
Configuracion=10
Esferas_conf=500000
File=Peli.dat

# Eleccion de lineas inicial y final
Start_line=$[($Configuracion)*($Esferas_conf+2)+1]
End_line=$[($Configuracion+1)*($Esferas_conf+2)]

# Creacion de la configuracion en fichero
sed -n -e "$Start_line,$End_line p" $File >> Conf_$Configuracion.dat 
```

The isolated configuration is stored in the file <code>Conf_n.dat</code>, where *n* is the number of the configuration.

<div id='id-section2'/>

## Post-processing and analysis ##
<div id='id-section2-1'/>

### Make up (Maquillador) ###
The algorithm <code>Maquillador2021Cadenas.f95</code> deletes the particles that have been overlapping in the process of aggregation of a configuration. When executing the code, the console asks you to enter by keyboard the name of the file (including the extension) that contains the configuration to make up, which is stored in *OVITO* format.
``` 
Dame el nombre del archivo (formato ovito):
>> Conf_n.dat
```

At the end, the number of spheres before and after the make-up process is shown on the screen. The new configuration is stored in the file <code>Maquillado.dat</code>.

<div id='id-section2-2'/>

### Pair correlation function and scattering intensity ###
The algorithm <code>grIq.f90</code> computes the pair correlation function and scattering intensity curve of a sphere configuration from its coordinates. To do this, it uses the parameters as input:

- **Nombre_fichero.dat**: Name of the file containing the configuration, in *OVITO* format.
- **r0**: Radius of the spheres that form the system.
- **dr**: *''Infinitesimal''* interval to discretize the radial variable in the computation of the pair correlation function.
- **qmin**: Minimum value of *q* for the calculation of the scattering intensity.
- **qmax**: Maximum value of *q* for the calculation of the scattering intensity.
- **dq**: *''Infinitesimal''* interval to discretize the variable *q*.

These parameters are read from the file <code>grIq.in</code>, whose default content is:
```
Conf_0_maq.dat		!Fichero de configuración
0.5			!Radio esferas
0.08			!dr
0.0010d0		!qmin
25.0d0			!qmax
1e-4			!dq
```

The pair correlation function is obtained by computing the number of pairs of particles that are within a certain distance *r*. This value is normalized according to the volume of a spherical shell of radius *r* and thickness *dr* and the number density. Subsequently, the scattering intensity is calculated as the product of the shape factor (in this case, that of a sphere of radius **r0**) and the structure factor, related to *g( r )* through an intergal expression ([Hasmy et al., 1995](https://www.sciencedirect.com/science/article/pii/0022309395000461)).

The code generates the files:

- **grIq.out**:  Contains a record of the parameters entered as input and calculated during the simulation.
- **gr.dat**: Contains the results of the pairwise correlation function. By columns: *r* and *g( r )*.
- **Iq.dat**: Contains the results of the scattering intensity curve. By columns: scattering vector, *q*, shape factor, *P( q )*, structure factor, *S( q )*, and scattering intensity, *I( q )*.


<div id='id-section2-3'/>

### External radius and radius of gyration ###
The algorithm <code>R_ext_giro.f90</code> calculates the external radius and the radius of gyration of the clusters of chains obtained after the aggregation process carried out by the algorithm [Arrangement (Ordenador)](#arrangement-ordenador). The code reads as input the data stored in the file <code>Centros.dat</code> generated (*see*  [Arrangement (Ordenador)](#arrangement-ordenador)). In addition, it asks that the name of the file (including the extension) containing the final configuration be entered by keyboard in *OVITO* format. 
``` 
Dame el nombre del archivo (formato ovito):
>> Conf_n.dat
```

As a result, the external radius and the radius of gyration of each aggregate are obtained, as well as the average of both values. All these values are stored in the file <code>R_exterior_giro.dat</code>. In addition, the file <code>CM_ovillos.dat</code>, which contains the coordinates of the centers of mass of the aggregates in *OVITO* format; and the file <code>Reposicionamiento.dat</code>, which rewrites the coordinates of all the spheres of the system eliminating the condition of periodic boundary conditions, so that the aggregates do not appear cut off. This format is useful to visualize the system more clearly.






