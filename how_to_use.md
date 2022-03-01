# Agregados de CNT #

Descripción del proceso de creación y análisis de sistemas de CNT que presentan distintos grados de agregación en el interior de un medio homogéneo.


![Sistemas de cadenas con distintos grados de dispersión generados mediante esta metodología ](https://i.ibb.co/F56rJ5F/figura-markdown.png)



## Índice
1. [Construcción de los sistemas](#id-section1)
	1.1. [Cadena prole](#id-section1-1)
	1.2. [Cadenas homogéneas](#id-section1-2)
	1.3. [Ordenador](#id-section1-3)
	1.4. [Configuración de equilibrio](#id-section1-4)
	1.5. [TakeConf](#id-section1-5)
2. [Post-procesado y análisis](#id-section2)
	2.1. [Maquillador](#id-section2-1)
	2.2. [Función de correlación de pares e intensidad de scattering ](#id-section2-2)
	2.3. [Radio externo y radio de giro ](#id-section2-3)

<div id='id-section1'/>

## Construcción de los sistemas ## 
<div id='id-section1-1'/>

### Cadena prole ###
El código <code>CadenaProle.f90</code> genera una cadena aislada mediante concatenación de esferas en condición de contacto duro. Durante el proceso de creación de la cadena, las esferas se concatenan sucesivamente, comprobando que no existe solapamiento entre ellas antes de aceptar cada nueva esfera. Los parámetros que pueden ajustarse para la generación de la cadena son:
- **nPartCadena**: Número de esferas que forman la cadena (*línea 40*).
- **theta**: Ángulo para colocar la siguiente esfera en el plano XY (*línea 76*).
- **phi**: Ángulo para colocar la siguiente esfera respecto al eje Z (*línea 78*).

El parámetro **nPartCadena** limita la longitud de la cadena, mientras que **theta** y **phi** definen la curvatura. Por defecto, estos parámetros toman los valores **nPartCadena = 1000**, **theta = pi/4** y **phi = pi/4**.

```fortran
nPartCadena=1000 !Numero de esferas para cada cadena

...

theta=theta*pi/4	!Angulo plano XY
phi=phi*pi/4		!Angulo eje Z
```

Si se fijan los valores **theta = 0** y **phi = 0** se genera una cadena rígida.

Los resultados se almacenan en el fichero <code> CadenaProle.dat</code>.

<div id='id-section1-2'/>

### Cadenas homogéneas ###
El código <code>CadenasHomogeneas_v3.f90</code> genera una distribución aleatoria y uniforme de cadenas en una caja de simulación cúbica. Los parámetros de entrada para generar el sistema son:
- **Nombre_fichero.dat**: Nombre del fichero que contiene los datos de la cadena que se va a utilizar como modelo de los CNT. Este archivo debe contener las posiciones de las esferas en formato *OVITO*: a) la primera línea del archivo es un número entero que coincide con el número de esferas que forma la cadena; b) la segunda línea contiene una descripción breve del sistema (línea de comentarios); c) las siguientes líneas contienen (por filas) el índice de la esfera, las coordenadas *(x,y,z)* y el radio.
- **Lcaja**: Tamaño de la arista de la caja de simulación cúbica.
- **nCad**: Número de cadenas que se colocan en la caja de simulación.

Estos parámetros se leen automáticamente del fichero <code>Parametros.inp</code> al ejecutar el código. Por defecto, toman los valores:
```txt
CadenaProle.dat		!Nombre del archivo que contiene la cadena modelo (.dat) formato OVITO
750.0d0			!Arista de la caja total dentro de la cual se generan los CM de las cadenas
500			!Numero de cadenas total del sistema
```

El algoritmo construye una caja de simulación cúbica de arista **Lcaja** y genera un total de **nCad** puntos aleatorios distribuidos de manera uniforme en su interior. Estos puntos aleatorios representan la posición de los centros de masa de las cadenas que modelan los CNT. Además, cada cadena se gira aleatoriamente antes de añadirse al sistema. Cuando todas las cadenas se han colocado, se reajusta la posición de todas las esferas del sistema aplicando condiciones de contorno periódicas (PBC).

La configuración de las cadenas se guarda en formato *OVITO* en el fichero <code>SistCadenasHomogeneas_v3.dat</code>.

<div id='id-section1-3'/>

### Ordenador ###
El algoritmo <code>Ordenador_v0.97.f90</code> realiza un proceso de Monte Carlo en el que las cadenas se agregan progresivamente en torno a distintos centros de atracción para dar lugar a la formación de agregados. Los parámetros  de entrada son:
- **Nombre_fichero.dat**: Nombre del fichero que contiene los datos de la configuración inicial del sistema. Este archivo debe contener las posiciones de las esferas en formato *OVITO*.
- **NEsf_Cad**: Número de esferas que compone cada cadena.
- **Ncentros**:  Número de centros atractores generados para formar los agregados de cadenas.
- **Niter**: Número de iteraciones máximas para relajar el sistema.
- **dr**: Tamaño del salto aleatorio que se desplazan las cadenas en cada iteración.
- **r0**: Radio de las esferas.
- **NPelicula**: Intervalo de iteraciones para grabar la configuración del sistema.

Estos parámetros se leen automáticamente del fichero <code>ordenador.in</code> al ejecutar el código. Por defecto, toman los valores:
```txt
SistCadenasHomogeneas_v3.dat	!Fichero con la configuración inicial de las cadenas
20				!Numero de esferas en cada cadena
1				!Numero de centros atractores
10000				!Numero de iteraciones maximas para relajar el sistema. 
0.1				!Tamanio del salto aleatorio. Como R=0.5, pueeees.... dr=0.1?
0.5				!Radio de las esferas
100				!Intervalo de iteraciones para grabar la configuracion
```

Tras cargar la configuración inicial, el algoritmo genera **Ncentros** puntos aleatorios en el interior de la caja de simulación que actúan como centros atractores y calcula, para cada cadena del sistema, cuál es su centro atractor más cercano en la configuración inicial. Además, calcula la energía inicial del sistema, definida como la suma de las distancias de los centros de masa a sus respectivos centros de atracción. En cada iteración, los centros de masas de las cadenas se desplazan aleatoriamente una distancia **dr**. El movimiento sólo se acepta si el centro de masa se acerca a su correspondiente centro atractor. Este proceso se repite hasta completar las **Niter** iteraciones, calculando en cada una de ellas las energía del sistema.

Al terminar el proceso de agregación se generan los ficheros:

- **Registro.out**: Contiene la información de los parámetros empleados en la simulación y una lista que enumera el centro atractor al que se acerca cada cadena.
- **CM_cadenas.dat**: Contiene las coordenadas de los centros de masa de las cadenas en la configuración inicial, en formato *OVITO*.
-  **Centros.dat**: Contiene las coordenadas de los centros atractores, en formato *OVITO*.
- **Registro_cad_virtual.dat**: Contiene las coordenadas de las esferas en la configuración inicial, pero eliminando las condiciones de contorno periódicas que se aplican al finalizar el algoritmo <code>CadenasHomogeneas_v3.f90</code>, de forma que las cadenas en la frontera de la caja de simulación no aparecen cortadas. Este formato es útil para visualizar el sistema con mayor claridad.
- **Energia.dat**: Contiene el valor de la energía del sistema en cada iteración del proceso de Monte Carlo.
- **Peli.dat**: Contiene la configuración del sistema cada **NPelicula** iteraciones, en formato *OVITO*.

<div id='id-section1-4'/>

### Configuración de equilibrio ###
El algoritmo <code>Nstep_Max.f95</code>  toma como entrada el fichero <code>Energia.dat</code> para calcular la iteración de Monte Carlo en la que se alcanza la configuración de equilibrio. Esta se define como la configuración en la que la energía varía menos de un 1% respecto a la configuración anterior. El resultado se almacena en el fichero <code> EnergiaNstepmax.dat</code>.

Dado que la configuración del sistema a lo largo del proceso de agregación se almacena cada **NPelicula** iteraciones, se elige como configuración de máxima agregación aquella que aparece en el fichero **Peli.dat** y que es superior a la iteración definida en <code> EnergiaNstepmax.dat</code>.

<div id='id-section1-5'/>

### TakeConf ###
El fichero <code>TakeConf.sh</code> permite extraer una configuración de la película generada en el proceso de agregación. Para ello, usa como parámetros de entrada:

- **Configuración**: Número de la configuración que quiere extraerse. Dado que las configuraciones se almacenan cada **NPelicula** iteraciones, elegir la configuración *N* implica seleccionar la iteración de Monte Carlo *Iter = N x NPelicula*.
- **Esferas_conf**: Número de esferas que conforma cada configuración.
- **File**: Nombre del fichero que almacena las configuraciones durante el proceso de agregación.

Por defecto, estos parámetros toman los valores **Configuración=10**, **Esferas_conf=500000** y **File=Peli.dat**:
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

La configuración aislada se almacena en el fichero <code>Conf_n.dat</code>, donde *n* es el número de la configuración.

<div id='id-section2'/>

## Post-procesado y análisis ##
<div id='id-section2-1'/>

### Maquillador ###
El algoritmo <code>Maquillador2021Cadenas.f95</code> elimina las partículas que han quedado solapando en el proceso de agregación de una configuración. Al ejecutar el código, pide que se introduzca por teclado el nombre del fichero  (incluyendo la extensión) que contiene la configuración a maquillar, que se encuentra almacenado en formato *OVITO*. 
``` 
Dame el nombre del archivo (formato ovito):
>> Conf_n.dat
```

Al finalizar, se muestra por pantalla el número de esferas antes y después del proceso de maquillado. La nueva configuración se almacena en el fichero <code>Maquillado.dat</code>.

<div id='id-section2-2'/>

### Función de correlación de pares e intensidad de scattering ###
El algoritmo <code>grIq.f90</code> calcula la función de correlación de pares y la curva de intesidad de scattering de una configuración de esferas a partir de sus coordenadas. Para ello, utiliza como entrada los parámetros:

- **Nombre_fichero.dat**: Nombre del fichero que contiene la configuración, en formato *OVITO*.
- **r0**: Radio de las esferas que forman el sistema.
- **dr**: Intervalo *''infinitesimal''* para discretizar la variable radial en el cálculo de la función de correlación de pares.
- **qmin**: Mínimo valor de *q* para el cálculo de la intensidad de scattering.
- **qmax**: Máximo valor de *q* para el cálculo de la intensidad de scattering.
- **dq**: Intervalo *''infinitesimal''* para discretizar la variable *q*.

Estos se leen del fichero <code>grIq.in</code>, cuyo contenido por defecto es:
```
Conf_0_maq.dat		!Fichero de configuración
0.5			!Radio esferas
0.08			!dr
0.0010d0		!qmin
25.0d0			!qmax
1e-4			!dq
```

La función de correlación de pares se obtiene computando el número de parejas de partículas que se encuentran a cierta distancia *r*. Este valor se normaliza según el volumen de una corteza esférica de radio *r* y espesor *dr* y la densidad numérica. Posteriormente, la intensidad de scattering se calcula como el producto del factor de forma (en este caso, el de una esfera de radio **r0**) y el factor de estrucutra, realcionado con *g( r )* mediante una expresión intergal ([Hasmy et al., 1995](https://www.sciencedirect.com/science/article/pii/0022309395000461)).

El código genera los ficheros:

- **grIq.out**:  Contiene un registro de los parámetros introducidos como entrada y calculados durante la simulación.
- **gr.dat**: Contiene los resultados de la función de correlación de pares. Por columnas: *r* y *g( r )*.
- **Iq.dat**: Contiene los resultados de la curva de intensidad de scattering. Por columnas: vector de scattering, *q*, factor de forma, *P( q )*, factor de estructura, *S( q )*, e intensidad de scattering, *I( q )*.


<div id='id-section2-3'/>

### Radio externo y radio de giro ###
El algoritmo <code>R_ext_giro.f90</code> calcula el radio externo y el radio de giro de los aglomerados de cadenas obtenidos tras el proceso de agregación llevado a cabo por el algoritmo [Ordenador](#id-section1-3). El código lee como input los datos almacenados en el fichero <code>Centros.dat</code> generado (*ver*  [Ordenador](#id-section1-3)). Además, pide que se introduzca por teclado el nombre del fichero (incluyendo la extensión) que contiene la configuración final cen formato *OVITO*. 
``` 
Dame el nombre del archivo (formato ovito):
>> Conf_n.dat
```

Como resultado se obtiene el radio externo y el radio de giro de cada agregado, así como el promedio de ambos valores. Todos estos valores se almacenan en el fichero <code>R_exterior_giro.dat</code>. Además, se genera el archivo <code>CM_ovillos.dat</code>, que contiene las coordenadas de los centros de masas de los agregados en formato *OVITO*; y el archivo <code>Reposicionamiento.dat</code>, que reescribe las coordenadas de todas las esferas del sistema eliminando la condición de condiciones de contorno periódicas, de forma que los agregados no aparezcan cortados. Este formato es útil para visualizar el sistema con mayor claridad.





