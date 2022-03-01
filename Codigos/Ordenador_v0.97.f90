! Sevilla, 20/02/2021 * 2 ventoso CCXXIX  * mes 12 del anio COVID
! vmorales@us.es
! laura25041997@gmail.com
! Universidad de Sevilla
! ------------------------------------------------
! Programa para crear agrupaciones de puntos en el espacio, 
! con distintos grados de homogenizacion:
! 	Partiremos de una serie de puntos al azar en el espacio
! 	y los vamos acercando hacia una serie de pozos de potencial 
! 	mediante relajación por Monte Carlo.
! Se trata de conseguir sistemas de puntos aleatorios, pero con distinto grado de homogenización
! ------------------------------------------------
! 28Abril2021: LGR 
! Modificaciones sobre el script inicial de Victor para esferas dispersas
! Se cambia la lectura de datos para iniciar con un sistema de cadenas dispersas y calcular
! las posiciones de sus centros de masas.
!--------------------------------------------------
! 21Mayo2021: LGR
! Se añaden PBC a la hora de calcular el centro atractivo mas cercano a cada cadena
! y tras mover cada una de estas.


PROGRAM ordenador

!--------------------------------------------------
! 1. Declaración de variables
!--------------------------------------------------
REAL(8), ALLOCATABLE 	:: XYZ(:,:)		!Configuracion del sistema (Indice,x,y,z,radio,nCadena)
REAL(8), ALLOCATABLE 	:: CENTROS(:,:)	!Posiciones de los centros atractores
REAL(8), ALLOCATABLE 	:: CM_CAD(:,:)	!Posiciones de los CM de las cadenas
REAL(8), ALLOCATABLE 	:: CAD_VIR(:,:)	!Cadena virtual para calcular los CM
INTEGER, ALLOCATABLE 	:: Cluster(:)	!Vector para guardar la etiqueta del centro atractor más cercano a cada cadena
REAL(8) 				:: Ldummy(1:3)	!Dummy para calcular el tamaño de la caja de simulacion
REAL(8) 				:: Lcaja		!Tamanio de la caja de simulacion
REAL(8) 				:: dist			!Distancia (variable temporal)
REAL(8) 				:: distMIN		!Distancia minima (variable temporal)
REAL(8) 				:: dist0		!Distancia en cada iteracion (variable temporal)
REAL(8) 				:: Etotal		!Energia total del sistema
REAL(8) 				:: dr 			!Tamanio del salto aleatorio
REAL(8) 				:: r0 			!Radio de las esferas
INTEGER					:: NPelicula	!Intervalo de iteraciones para grabar la configuracion
REAL(8)					:: azar(1:3)	!Generacion de los centros aleatorios
INTEGER 				:: Npart		!Numero de particulas (esferas)
INTEGER 				:: Ncentros 	!Numero de centros atractores
INTEGER 				:: Niter		!Numero de iteraciones para relajar el sistema
INTEGER 				:: NEsf_Cad		!Numero de esferas por cadena
INTEGER 				:: NCad			!Numero de cadenas
CHARACTER(len=50)		:: nombrefich	!Nombre del fichero con la configuracion inicial
CHARACTER(len=50)		:: DummyChar	!Muda para variables de caracteres
REAL(8) 				:: xCent		!Dummy para calcular los CM de las cadenas
REAL(8) 				:: yCent		!Dummy para calcular los CM de las cadenas
REAL(8) 				:: zCent		!Dummy para calcular los CM de las cadenas
REAL(8) 				:: xdist		!Dummy para calcular interdistancias
REAL(8) 				:: ydist		!Dummy para calcular interdistancias
REAL(8) 				:: zdist		!Dummy para calcular interdistancias
REAL(8) 				:: x			!Dummy para desplazar los CM de las cadenas
REAL(8) 				:: y			!Dummy para desplazar los CM de las cadenas
REAL(8) 				:: z			!Dummy para desplazar los CM de las cadenas
REAL(8)					:: Xcambio		!Cambio de la posicion del CM de las cadenas
REAL(8)					:: Ycambio		!Cambio de la posicion del CM de las cadenas
REAL(8)					:: Zcambio		!Cambio de la posicion del CM de las cadenas
INTEGER					:: aux			!Auxiliar entero para conteo de CAD_VIR
INTEGER					:: i, j, k		!Indices

! Paranoia del random (semilla aleatoria)
INTEGER the_size
INTEGER, ALLOCATABLE :: seed(:)
! =============================================
! La movida aleatoria: esto es muy ortopédico; 	--> MIRAR SEMILLA ALEATORIA GUILLERMO (MONTECARLO B1)
! hay que buscar otra fórmula más elegante
! open(12,file='seed.in')	
! the_size=4
! write (*,*) 'the_size=',the_size
! call random_seed(size=the_size)
! allocate (seed(1:the_size))
! do i=1,the_size
! read(12,*) seed(i)
! write (*,*) i,the_size !,seed(i)
! enddo
! close(12)
! Call random_seed(put=seed)
! =============================================


!--------------------------------------------------
! 2. INICIALIZACION Y CARGA DE DATOS
!--------------------------------------------------
! Lectura de los parametros iniciales (archivo .IN)
OPEN(UNIT=21,FILE="ordenador.in",STATUS='OLD')
	READ(21,*) nombrefich	 	!Fichero con la configuracion inicial de las cadenas
	READ(21,*) NEsf_Cad			!Numero de esferas en cada cadena
	READ(21,*) Ncentros 		!Numero de centros atractores
	READ(21,*) Niter 			!Numero de iteraciones maximas para relajar el sistema. 
	READ(21,*) dr 				!Tamanio del salto aleatorio. Como R=0.5, pueeees.... dr=0.1?
	READ(21,*) r0 				!Radio de las esferas
	READ(21,*) NPelicula		!Intervalo de iteraciones para grabar la configuracion
CLOSE(21)

! Lectura de la configuracion inicial del sistema
Lcaja = 0.0d0 !Inicializa tamaño de la caja a cero
OPEN(UNIT=22,FILE=nombrefich,STATUS='OLD')
	READ(22,*) Npart	 	!Numero de esferas totales
	READ(22,*) DummyChar	!Comentarios (formato OVITO)
	ALLOCATE(XYZ(Npart,7))
	DO i=1,Npart
		READ(22,*) XYZ(i,1), XYZ(i,2), XYZ(i,3), XYZ(i,4), XYZ(i,5), XYZ(i,6)
		Ldummy(1) =  XYZ(i,2)
		Ldummy(2) =  XYZ(i,3)
		Ldummy(3) =  XYZ(i,4)
!~ 		XYZ(i,6) = 1
		Lcaja = MAX(MAXVAL(Ldummy),Lcaja)
	ENDDO
CLOSE(22)
WRITE(*,*) 'Sistema cargado.'


! Tamanio de la caja de simulacion
Lcaja = 2.0*Lcaja
PRINT*, 'Lcaja = ', Lcaja

! Generacion de los centros atractores en posiciones aleatorias
NCad = INT(XYZ(Npart,6))
ALLOCATE(CENTROS(3,Ncentros))
ALLOCATE(Cluster(NCad))  

OPEN(UNIT=23, FILE="Centros.dat",STATUS='UNKNOWN')
WRITE (23,*) Ncentros
WRITE (23,*) 'COMMENTS'
DO i=1,Ncentros
	CALL RANDOM_NUMBER (azar) 
	CENTROS(1,i) = Lcaja * (azar(1) - 0.5)
	CENTROS(2,i) = Lcaja * (azar(2) - 0.5)
	CENTROS(3,i) = Lcaja * (azar(3) - 0.5)
	WRITE (23,*) i, CENTROS(1,i), CENTROS(2,i), CENTROS(3,i), 0.2d0  !Radio distinto para representacion
ENDDO
WRITE(*,*) 'Centros atractores generados.'

OPEN(UNIT=31, FILE='Registro.out', STATUS='UNKNOWN')

	WRITE(31,*) 'Fichero configuracion inicial: ', nombrefich	
	WRITE(31,*) 'Numero de esferas en cada cadena: ', NEsf_Cad		
	WRITE(31,*) 'Numero de centros atractores: ', Ncentros 	
	WRITE(31,*) 'Numero de iteraciones maximas para relajar el sistema: ', Niter 			
	WRITE(31,*) 'Tamanio del salto aleatorio: ', dr 				
	WRITE(31,*) 'Radio de las esferas: ', r0 				
	WRITE(31,*) 'Intervalo de iteraciones para grabar la configuracion: ', NPelicula		
	WRITE(31,*) 'Numero de cadenas: ', NCad
	WRITE(31,*) 'Lcaja: ', Lcaja
	WRITE(31,*) ' '

!--------------------------------------------------
! 3. CALCULO DEL CENTRO DE MASA DE CADA CADENA
!--------------------------------------------------
OPEN(55,FILE='Registro_cad_virtual.dat',STATUS='UNKNOWN')
WRITE(55,*) Npart
WRITE(55,*) 'COMMENTS'

ALLOCATE(CAD_VIR(NEsf_Cad,3))
ALLOCATE(CM_CAD(NCad,3))
DO i=1,NCad
	CAD_VIR = 0.0d0
	
	! La primera esfera de la cadena virtual coincide con  la primera 
	! esfera de la cadena real
	CAD_VIR(1,1) = XYZ(((i-1)*NEsf_Cad+1),2)
	CAD_VIR(1,2) = XYZ(((i-1)*NEsf_Cad+1),3)
	CAD_VIR(1,3) = XYZ(((i-1)*NEsf_Cad+1),4)
	WRITE(55,*) 0, CAD_VIR(1,1), CAD_VIR(1,2), CAD_VIR(1,3), 0.5
	
	
	! En cada nueva cadena se inicializa el centro de masas con la
	! primera esfera (original = virtual)
	xCent = CAD_VIR(1,1)
	yCent = CAD_VIR(1,2)
	zCent = CAD_VIR(1,3)
	aux = 1
	
	DO j=((i-1)*NEsf_Cad+1),(i*NEsf_Cad-1)
		! Posicion relativa r_i y r_i+1
		xdist = XYZ(j+1,2) - CAD_VIR(aux,1)
		ydist = XYZ(j+1,3) - CAD_VIR(aux,2)
		zdist = XYZ(j+1,4) - CAD_VIR(aux,3)
		
		! Correccion de la nueva esfera si la distancia es mayor de L/2
		CAD_VIR(aux+1,1) = XYZ(j+1,2) - Lcaja * DNINT( xdist / Lcaja )
		CAD_VIR(aux+1,2) = XYZ(j+1,3) - Lcaja * DNINT( ydist / Lcaja )
		CAD_VIR(aux+1,3) = XYZ(j+1,4) - Lcaja * DNINT( zdist / Lcaja )
		
		WRITE(55,*) 0, CAD_VIR(aux+1,1), CAD_VIR(aux+1,2), CAD_VIR(aux+1,3), 0.5
		
		! Calcula el CM con la CAD_VIR
		xCent = xCent + CAD_VIR(aux+1,1)
		yCent = yCent + CAD_VIR(aux+1,2)
		zCent = zCent + CAD_VIR(aux+1,3)
		aux = aux + 1
	ENDDO
	
	! Normalizacion con NEsf_Cad para calcular el CM
	xCent = xCent / NEsf_Cad
	yCent = yCent / NEsf_Cad
	zCent = zCent / NEsf_Cad
	
	! PBC por si el CM calculado con CAD_VIR se sale de la caja
	CM_CAD(i,1) = xCent - Lcaja * DNINT( xCent / Lcaja )
	CM_CAD(i,2) = yCent - Lcaja * DNINT( yCent / Lcaja )
	CM_CAD(i,3) = zCent - Lcaja * DNINT( zCent / Lcaja )
ENDDO
WRITE(*,*) 'Posiciones de los CM de las cadenas calculados.'

OPEN(UNIT=32, FILE='CM_cadenas.dat', STATUS='UNKNOWN')
WRITE(32,*) NCad+Ncentros
WRITE(32,*) 'COMMENTS'
DO i=1,Ncentros
	WRITE(32,*) i, CENTROS(1,i), CENTROS(2,i), CENTROS(3,i), 2
ENDDO
DO i=1,NCad
	WRITE(32,*) i, CM_CAD(i,:), 0.2
ENDDO


!--------------------------------------------------
! 4. ENERGIA INICIAL DEL SISTEMA Y CALCULO DEL CENTRO MAS CERCANO A CADA CADENA
!--------------------------------------------------
! Calculamos a qué centro atractor se acercará cada punto para evitar tener que calcularlo cada vez en el bruten-bucle....
! aprovecho para calcular la "Etotal" inicial:
! En realidad, calculamos las distancias, Ei=dist(particula_i - su centro atractor)
! Etotal=sum Ei = sum [i=1,Npart] (dist(i,centro(i)))
! Creo ques lo más sencillo => lo más rápido

! Inicializacion de la "energía" total. Este valor se minimiza en el apartado 5
Etotal = 0.00000

! Calculo de la "energia" total inicial y del centro atractor al que se acerca cada cadena
DO i=1,NCad
	distMIN = Lcaja*Lcaja  ! Inicializacion con una distancia absurdamente enorme
	DO j=1,Ncentros
	
		xdist = CM_CAD(i,1) - CENTROS(1,j)	!rij CM de la cadena i al centro atractor j
		ydist = CM_CAD(i,2) - CENTROS(2,j)
		zdist = CM_CAD(i,3) - CENTROS(3,j)
		
		xdist = xdist - Lcaja * DNINT( xdist / Lcaja )	!Corrección rij mediante MIC
		ydist = ydist - Lcaja * DNINT( ydist / Lcaja )
		zdist = zdist - Lcaja * DNINT( zdist / Lcaja )
	
		dist = DSQRT( xdist ** 2.0d0 + ydist ** 2.0d0 + zdist ** 2.0d0 )
		IF (dist.LT.distMIN) THEN
			distMIN = dist
			Cluster(i) = j
		ENDIF
	ENDDO
	Etotal=Etotal+distMIN
ENDDO
! Ya sabemos a qué centro atractor se tendrá que acercar cada punto del 
! sistema y ya tenemos la Etotal del sistema, que deberá ir reduciéndose

! Se guarda en el fichero de registro a que cluster se agrega cada cadena
WRITE(31,*) 'CLUSTER AL QUE SE ACERCA CADA CADENA'
DO i=1,NCad
	WRITE(31,*) i, Cluster(i)
ENDDO
WRITE(31,*) ' '


! Registro de la pelicula y la variacion de energia total
OPEN(UNIT=24,FILE='Energia.dat',STATUS='UNKNOWN')
WRITE(24,*) 0, Etotal

OPEN(UNIT=25,FILE='Peli.dat',STATUS='UNKNOWN')
WRITE(25,*) Npart
WRITE(25,*) 'COMMENTS'
DO i=1,Npart
	WRITE(25,*) XYZ(i,1), XYZ(i,2), XYZ(i,3), XYZ(i,4), XYZ(i,5), XYZ(i,6)
ENDDO

WRITE(*,*) 'Energia y configuracion inicial registradas.'

!--------------------------------------------------
! 5. RELAJACION DEL SISTEMA (SE ACERCAN LAS CADENAS ALEATORIAMENTE)
!--------------------------------------------------
! Bucle para el numero de iteraciones para minimizar la energia
DO k=1,Niter 	! Variando Niter y dr, el ancho del salto aleatorios, 
				! se controla la velocidad de acercamiento
WRITE(*,*) 'ITERACION', k, 'de', Niter


IF(MOD(k, NPelicula) == 0) THEN
	WRITE(32,*) NCad+Ncentros
	WRITE(32,*) 'COMMENTS'
	DO i=1,Ncentros
		WRITE(32,*) i, CENTROS(1,i), CENTROS(2,i), CENTROS(3,i), 2
	ENDDO
	DO i=1,NCad
		WRITE(32,*) i, CM_CAD(i,:), 0.2
	ENDDO
ENDIF

! Inicializamos la Etotal a cero en cada iteracion
Etotal = 0.00000

	! Bucle para cada CM de las cadenas
	DO i=1,NCad
		xdist = CM_CAD(i,1) - CENTROS(1,Cluster(i))	!rij inicial del CM de la cadena i al centro atractor Cluster(i)=j
		ydist = CM_CAD(i,2) - CENTROS(2,Cluster(i))
		zdist = CM_CAD(i,3) - CENTROS(3,Cluster(i))
			
		xdist = xdist - Lcaja * DNINT( xdist / Lcaja )	!Corrección rij mediante MIC
		ydist = ydist - Lcaja * DNINT( ydist / Lcaja )
		zdist = zdist - Lcaja * DNINT( zdist / Lcaja )	
	
		dist0 = DSQRT( xdist ** 2.0d0 + ydist ** 2.0d0 + zdist ** 2.0d0 )
		
		IF (dist0.LE.dr) GOTO 40 !Si dist0.LE.dr, es jodido o imposible que se acerque más: Ve a por otra cadena
			
		dist = Lcaja !Inicializo una distancia grande para entrar en el do-while	
		! Movemos al azar hasta que se acerque un poquillo:
		DO WHILE (dist.GE.dist0)
			CALL RANDOM_NUMBER (azar)
			x = CM_CAD(i,1) + dr * (azar(1)-0.5)
			y = CM_CAD(i,2) + dr * (azar(2)-0.5)
			z = CM_CAD(i,3) + dr * (azar(3)-0.5)
			
			xdist = x - CENTROS(1,Cluster(i))	!Nuevo rij CM de la cadena i al centro atractor Cluster(i)=j
			ydist = y - CENTROS(2,Cluster(i))
			zdist = z - CENTROS(3,Cluster(i))
			
			xdist = xdist - Lcaja * DNINT( xdist / Lcaja )	!Corrección rij mediante MIC
			ydist = ydist - Lcaja * DNINT( ydist / Lcaja )
			zdist = zdist - Lcaja * DNINT( zdist / Lcaja )	
			
			dist = DSQRT( xdist ** 2.0d0 + ydist ** 2.0d0 + zdist ** 2.0d0 )
		ENDDO	!Fin bucle DO WHILE
		
		! Si llega aqui, es que el CM de la cadena i se ha acercado a su centro atractor
		Etotal=Etotal+dist
		
		! Variacion del centro de masa de cada cadena
		Xcambio = x - CM_CAD(i,1)	
		Ycambio = y - CM_CAD(i,2)
		Zcambio = z - CM_CAD(i,3)
		
		! Almacenamos las nuevas coordenadas de cada esfera de la cadena i
		CM_CAD(i,1) = x	- Lcaja * DNINT( x / Lcaja )	!Nueva posicion del CM de la cadena i
		CM_CAD(i,2) = y - Lcaja * DNINT( y / Lcaja )
		CM_CAD(i,3) = z - Lcaja * DNINT( z / Lcaja )
		
		DO j=1,NEsf_Cad		!Nueva posicion de todas las esferas de la cadena i
			XYZ(NEsf_Cad*(i-1)+j,2) = XYZ(NEsf_Cad*(i-1)+j,2) + Xcambio
			XYZ(NEsf_Cad*(i-1)+j,3) = XYZ(NEsf_Cad*(i-1)+j,3) + Ycambio
			XYZ(NEsf_Cad*(i-1)+j,4) = XYZ(NEsf_Cad*(i-1)+j,4) + Zcambio
			
			XYZ(NEsf_Cad*(i-1)+j,2) = XYZ(NEsf_Cad*(i-1)+j,2) - Lcaja * DNINT( XYZ(NEsf_Cad*(i-1)+j,2) / Lcaja )	!Aplicacion PBC tras el movimiento de la cadena
			XYZ(NEsf_Cad*(i-1)+j,3) = XYZ(NEsf_Cad*(i-1)+j,3) - Lcaja * DNINT( XYZ(NEsf_Cad*(i-1)+j,3) / Lcaja )
			XYZ(NEsf_Cad*(i-1)+j,4) = XYZ(NEsf_Cad*(i-1)+j,4) - Lcaja * DNINT( XYZ(NEsf_Cad*(i-1)+j,4) / Lcaja )
		ENDDO
		! NOTA: Ignoramos si solapan o no entre sí, porque después, 
		! usaremos estos puntos para construir sistemas de esferas (alineadas, racimos...)
		! Es entonces cuando esos sistemas de esferas deberán garantizar que sus esferas no solapen.
		
40 ENDDO 

	WRITE (24,*) k, Etotal  ! NOTA: Registramos la energia para ir viendo cómo desciende

	! Salvar las configuraciones en formato xyz para Ovito 
	! cada NPelicula iteraciones
	IF(MOD(k, NPelicula) == 0) THEN
		WRITE(25,*) Npart
		WRITE(25,*) 'Lennard-Jones esfera'
		DO i = 1, Npart
			WRITE(25,*) XYZ(i,1), XYZ(i,2), XYZ(i,3), XYZ(i,4), XYZ(i,5), XYZ(i,6)
		ENDDO  
	ENDIF
	
END DO  !Final bucle para Niter

CLOSE(32)
CLOSE(24)
CLOSE(25)
WRITE(*,*) 'Minimizacion de la energia terminada.'
!--------------------------------------------------
WRITE(*,*) 'Done!'
END PROGRAM ordenador


