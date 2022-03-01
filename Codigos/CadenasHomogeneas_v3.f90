program SistCadenasJuntas
implicit none
!****************************************
!****************************************
! Sevilla, 2021
! Creacion de un sistema de cadenas homogeneamente dispersas en el 
! interior de una caja de simulacion cubica centrada en el origen de 
! coordenadas. La cadena modelo se carga a partir de un archivo externo.
! laugarreg@alum.us.es
!****************************************
!****************************************

!----------------------------------------
! DECLARACION DE VARIABLES
REAL(8), PARAMETER :: pi=3.1415926535897932384626
CHARACTER(len=35) :: nombfich, comentarios
INTEGER :: nPartCadena, i, k, m, mudaInt, nPartSistema, nCad, countCad
REAL(8), DIMENSION(:), ALLOCATABLE :: xCad, yCad, zCad
REAL(8) :: r0, xCentCadena, yCentCadena, zCentCadena, Lcaja
REAL(8) :: xx, yy, zz
REAL(8), DIMENSION(:,:), ALLOCATABLE :: Sist, At
REAL(8), DIMENSION(:), ALLOCATABLE :: xGiro, CadGir
CHARACTER(len=10) :: fecha,hora !Semilla aletoria
INTEGER :: the_size !Semilla aletoria
INTEGER, DIMENSION(:), ALLOCATABLE :: seed !Semilla aletoria
REAL(8) :: xxx !Semilla aletoria
REAL(8) :: theta, phi	!Angulos para direccion aleatoria de ur
REAL(8), DIMENSION(:), ALLOCATABLE :: ur
REAL(8) :: urmod
REAL(8) :: ang_giro
REAL(8) :: xCentSist, yCentSist, zCentSist
REAL(8)	:: Lmax_ini, Lmax_end

!----------------------------------------
! INICIALIZACION DE LA SEMILLA ALEATORIA  !!Creador2019 de Victor
CALL DATE_AND_TIME(DATE = fecha,TIME = hora)
the_size=1
ALLOCATE(seed(1:the_size))
READ(hora,*) xxx
seed=INT(1000*(xxx-FLOAT(INT(xxx)))) ! seed=1000veces la parte entera de xxx.
CALL RANDOM_SEED
call random_seed(size=the_size)
CALL RANDOM_SEED(PUT=SEED(1:the_size))

!----------------------------------------
! LECTURA DE PARAMETROS DE CREACION DE UN SISTEMA DE CADENAS DISTRIBUIDAS ALEATORIAMENTE EN 3D
OPEN(UNIT=11, FILE='Parametros.inp', STATUS='OLD')
	READ(11,*) nombfich	!Nombre del archivo que contiene la cadena modelo (.dat) formato OVITO
	READ(11,*) Lcaja	!Arista de la caja total dentro de la cual se generan los CM de las cadenas
	READ(11,*) nCad		!Numero de cadenas total del sistema
CLOSE(11)

!----------------------------------------
! CARGAR LOS DATOS DE UNA CADENA MODELO
OPEN(UNIT=12,FILE=nombfich,STATUS='OLD')
READ(12,*) nPartCadena	!Numero de esferas por cadena
READ(12,*) comentarios
ALLOCATE(xCad(nPartCadena))
ALLOCATE(yCad(nPartCadena))
ALLOCATE(zCad(nPartCadena))
DO i=1,nPartCadena
	READ(12,*) mudaInt, xCad(i), yCad(i), zCad(i), r0
ENDDO
CLOSE(12)

nPartSistema = nCad * nPartCadena	!Numero de particulas esfericas
ALLOCATE(Sist(nPartSistema,6))

!----------------------------------------
! CALCULO DEL CENTRO DE MASAS DE LA CADENA MODELO
xCentCadena = 0.0d0
yCentCadena = 0.0d0
zCentCadena = 0.0d0
DO i=1,nPartCadena
	xCentCadena = xCentCadena + xCad(i)
	yCentCadena = yCentCadena + yCad(i)
	zCentCadena = zCentCadena + zCad(i)
ENDDO
xCentCadena = xCentCadena / nPartCadena
yCentCadena = yCentCadena / nPartCadena
zCentCadena = zCentCadena / nPartCadena
WRITE(*,*) 'CM cadena modelo: ', xCentCadena, yCentCadena, zCentCadena

!----------------------------------------
! REESCRIBIR LA POSICION DE CADA ESFERA SEGUN LA POSICION RELATIVA AL CM EN LA CADENA MODELO
DO i=1,nPartCadena
	xCad(i) = xCad(i) - xCentCadena
	yCad(i) = yCad(i) - yCentCadena
	zCad(i) = zCad(i) - zCentCadena
ENDDO

!----------------------------------------
! CREACION DEL SISTEMA COMPLETO
ALLOCATE(At(3,3))		!Matriz de giro
ALLOCATE(xGiro(3))		!Coordenadas de cada esfera antes del giro
ALLOCATE(CadGir(3))		!Coordenadas de cada esfera despues del giro
ALLOCATE(ur(3))			!Vector de tres componentes: eje de giro de cada cadena

mudaInt = 1 	!Cuenta todas las esferas del sistema
countCad = 1  	!Cuenta todas las cadenas del sistema

DO k=1,nCad   !Genera todas las cadenas del sistema en la caja Lcaja
101 mudaInt = (k-1)*nPartCadena + 1

	!CM de cada cadena dentro de la caja Lcaja	
	CALL RANDOM_NUMBER(xx)  
	xx = Lcaja * ( xx - 0.5d0 )
	CALL RANDOM_NUMBER(yy)
	yy = Lcaja * ( yy - 0.5d0 )
	CALL RANDOM_NUMBER(zz)
	zz = Lcaja * ( zz - 0.5d0 )
		
	!Definicion del vector w que pasa por CM de la cadena y en torno al que se rota la cadena
	!Generacion de una direccion aleatoria en la superficie de una corteza esferica     
     CALL RANDOM_NUMBER(theta)			!Numero al azar en [0,1]
     theta = ACOS(2.0d0*theta - 1.0d0)	!Numero al azar en [-1,1]
	 CALL RANDOM_NUMBER(phi)			!Numero al azar en [0,1]
     phi = 2.0d0*pi*phi					!Numero al azar en [0,2pi]
	
	!Componentes y modulo del vector de direccion aleatoria
	ur(1) = 1.0d0 * SIN(theta) * COS(phi)
	ur(2) = 1.0d0 * SIN(theta) * SIN(phi)
	ur(3) = 1.0d0 * COS(theta)
	urmod = DSQRT(SUM(ur*ur))	!Modulo del vector w
	ur = ur / urmod				!Vector unitario en la direccion de w
		
	!Angulo que se gira en torno al eje ur
	CALL RANDOM_NUMBER(ang_giro)  		
	ang_giro = 2.0d0 * pi * ang_giro
	
	!Creacion de la matriz de rotacion	
	At(1,1) = dble( ur(1)**2.0*(1-cos(ang_giro)) + cos(ang_giro) )     
	At(1,2) = dble( ur(1)*ur(2)*(1-cos(ang_giro)) + ur(3)*sin(ang_giro) )
	At(1,3) = dble( ur(1)*ur(3)*(1-cos(ang_giro)) - ur(2)*sin(ang_giro) )
	At(2,1) = dble( ur(1)*ur(2)*(1-cos(ang_giro)) - ur(3)*sin(ang_giro) ) 
	At(2,2) = dble( ur(2)**2.0*(1-cos(ang_giro)) + cos(ang_giro) )
	At(2,3) = dble( ur(2)*ur(3)*(1-cos(ang_giro)) + ur(1)*sin(ang_giro) )
	At(3,1) = dble( ur(1)*ur(3)*(1-cos(ang_giro)) + ur(2)*sin(ang_giro) )
	At(3,2) = dble( ur(2)*ur(3)*(1-cos(ang_giro)) - ur(1)*sin(ang_giro) )
	At(3,3) = dble( ur(3)**2.0*(1-cos(ang_giro)) + cos(ang_giro) )
				
	DO m=1,nPartCadena !Se suma cada esfera de la cadena a la posicion inicial CM_Cadena=(xx,yy,zz)
		xGiro = (/xCad(m),yCad(m),zCad(m)/)	!Posicion de cada esfera de la cadena modelo respecto al CM
		xGiro = matmul(At,xGiro)				!Se gira cada esfera de la cadena segun la matriz de giro At
		CadGir(1) = xGiro(1)	!Componentes de la posicion de cada esfera de la cadena girada respecto al eje ur
		CadGir(2) = xGiro(2)
		CadGir(3) = xGiro(3)
			
		Sist(mudaInt,1) = mudaInt			!Esfera mudaInt
		Sist(mudaInt,2) = xx + CadGir(1)	!Componente x
		Sist(mudaInt,3) = yy + CadGir(2)	!Componente y
		Sist(mudaInt,4) = zz + CadGir(3)	!Componente z
		Sist(mudaInt,5) = r0				!Radio
		Sist(mudaInt,6) = countCad		!Cadena a la que pertenece
		
		mudaInt=mudaInt+1	!Se aumenta en 1 el indice de la esfera
	ENDDO !Genera cada esfera de la cadena k actual
	
	countCad = countCad + 1 !Se aumenta en 1 el indice de la cadena
ENDDO !Genera todas las cadenas

!----------------------------------------
!~ ! APLICACION PBC
Lmax_ini = 0.0d0
Lmax_end = 0.0d0
DO i=1, nPartSistema
	Lmax_ini = MAX( Lmax_ini, ABS(Sist(i,2)), ABS(Sist(i,3)), ABS(Sist(i,4)) )
	Sist(i,2) = Sist(i,2) - Lcaja * DNINT( Sist(i,2) / Lcaja )
	Sist(i,3) = Sist(i,3) - Lcaja * DNINT( Sist(i,3) / Lcaja )
	Sist(i,4) = Sist(i,4) - Lcaja * DNINT( Sist(i,4) / Lcaja )
	Lmax_end = MAX( Lmax_end, ABS(Sist(i,2)), ABS(Sist(i,3)), ABS(Sist(i,4)) )
ENDDO
print*, 'Lmax_ini = ', 2*Lmax_ini
print*, 'Lmax_end = ', 2*Lmax_end

!----------------------------------------
! PARAMETROS GEOMETRICOS FINALES
! Centro de masas del sistema
xCentSist = 0.0d0
yCentSist = 0.0d0
zCentSist = 0.0d0
DO i=1,nPartSistema
	xCentSist = xCentSist + Sist(i,2)
	yCentSist = yCentSist + Sist(i,3)
	zCentSist = zCentSist + Sist(i,4)
ENDDO
xCentSist = xCentSist / nPartSistema
yCentSist = yCentSist / nPartSistema
zCentSist = zCentSist / nPartSistema
WRITE(*,*) 'CM del sistema: ', xCentSist, yCentSist, zCentSist

!----------------------------------------
! SAVE FOR OVITO
! Esferas del sistema de cadenas homogeneamente distribuidas 
OPEN(UNIT=13,FILE='SistCadenasHomogeneas_v3.dat',STATUS='UNKNOWN')
WRITE(13,*) nPartSistema
WRITE(13,*) 'COMMENTS'
DO i=1,nPartSistema
	WRITE(13,*) Sist(i,:)
ENDDO
CLOSE(13)

!----------------------------------------
WRITE(*,*) 'Done!'
!----------------------------------------
end program SistCadenasJuntas
