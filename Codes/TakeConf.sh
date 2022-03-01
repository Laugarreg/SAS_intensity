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

