program CadenaProleBien
implicit none
!****************************************
!****************************************
!
! Cracion de una cadena con ángulo reestringido y sin eje privilegiado
!
!****************************************
!****************************************

!----------------------------------------
! DECLARACION DE VARIABLES
REAL(8), PARAMETER :: pi=3.1415926535897932384626
INTEGER :: nPartCadena, nIntentos, i, j
REAL(8) :: r0, x0, y0, z0, d, aux, theta, phi, x, y, z, dist
REAL(8) :: xCent, yCent, zCent
!CHARACTER(len=35) :: nomfich
REAL(8), DIMENSION(:), ALLOCATABLE :: xPart, yPart, zPart
CHARACTER(len=10) :: fecha,hora !Semilla aletoria
INTEGER :: the_size !Semilla aletoria
INTEGER, DIMENSION(:), ALLOCATABLE :: seed !Semilla aletoria
REAL(8) :: xxx !Semilla aletoria
REAL(8), DIMENSION(:), ALLOCATABLE :: xAxis,yAxis, zAxis, dnew, dcoord
REAL(8), DIMENSION(:,:), ALLOCATABLE :: matAxis, matR

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
! PARAMETROS DE SIMULACION
r0=0.5 ! Radio de la particula elemental (unidades reducidas)
nPartCadena=1000 ! Numero de esferas para cada cadena (sale del ratio longitud/diametro de los CNTS)
x0=0.0 ! Coordenadas de la primera partícula
y0=0.0
z0=0.0
ALLOCATE(xAxis(3))  !Ejes de referencia iniciales
ALLOCATE(yAxis(3))
ALLOCATE(zAxis(3))
xAxis=(/1.0,0.0,0.0/)
yAxis=(/0.0,1.0,0.0/)
zAxis=(/0.0,0.0,1.0/)

!----------------------------------------
! Creacion de la cadena
ALLOCATE(xPart(0:nPartCadena)) ! Se reserva xiPart(0) para el centro de masas
ALLOCATE(yPart(0:nPartCadena))
ALLOCATE(zPart(0:nPartCadena))
xPart(1)=x0 ! Primera partícula
yPart(1)=y0
zPart(1)=z0

ALLOCATE(matAxis(3,3)) !Matriz de ejes iniciales
matAxis(1,:)=xAxis
matAxis(2,:)=yAxis
matAxis(3,:)=zAxis

ALLOCATE(matR(3,3)) !Matriz de rotacion inicial
matR=matAxis

ALLOCATE(dnew(3))
ALLOCATE(dcoord(3))
d=2*r0  !Distancia de contacto duro entre esferas
DO i=2,nPartCadena ! Bucle para generar el resto de esferas de la cadena
	nIntentos=1
	aux=1
	DO WHILE ((aux/=0).AND.(nIntentos<1000))
		CALL RANDOM_NUMBER(theta)  ! Angulos aleatorios para generar la posicion de la nueva particula
		theta=theta*pi/4		!Angulo plano XY
		CALL RANDOM_NUMBER(phi)
		phi=phi*pi/4			!Angulo eje Z		
		
		dcoord(1)=d*sin(theta)*cos(phi)  !Vector d expresado en los nuevos ejes
		dcoord(2)=d*sin(theta)*sin(phi)
		dcoord(3)=d*cos(theta)
		
		matAxis=matmul(matR,matAxis)
		dnew=matmul(dcoord,matAxis)   !Vector d en los ejes originales
		
		x=xPart(i-1) + dnew(1) ! Coordenadas de la nueva particula en los ejes de referencia iniciales
		y=yPart(i-1) + dnew(2)
		z=zPart(i-1) + dnew(3)
		
		dist=0 ! Inicializamos dist=0
		aux=0		
		DO j=1,(i-1) ! Bucle para comprobar que la nueva particula no solapa con ninguna de las anteriores
			dist=sqrt((x-xPart(j))**2 + (y-yPart(j))**2 + (z-zPart(j))**2)
			IF (dist<1) THEN 
				aux=1
				EXIT  !Si alguna solapa, sale del bucle y genera un nuevo intento
			ENDIF
		ENDDO 
		
		IF (aux==0) THEN !Bucle para asignar la nueva particula o aumentar nIntentos
			xPart(i)=x
			yPart(i)=y
			zPart(i)=z
		ELSE 
			nIntentos=nIntentos+1
		ENDIF
		
		matR(1,1)=dble(sin(phi))  !Matriz R de los ejes nuevos (i) respecto a los ejes anteriores (i-1)
		matR(1,2)=dble(-cos(phi))
		matR(1,3)=dble(0.0)
		matR(2,1)=dble(cos(theta)*cos(phi))  
		matR(2,2)=dble(cos(theta)*sin(phi))
		matR(2,3)=dble(-sin(theta))
		matR(3,1)=dble(sin(theta)*cos(phi))  
		matR(3,2)=dble(sin(theta)*sin(phi))
		matR(3,3)=dble(cos(theta))
		
	ENDDO ! End While
ENDDO ! End bucle i para cada particula nueva

!----------------------------------------
! Calculo del centro de masas xiPart(0)
xCent=0
yCent=0
zCent=0
DO i=1,nPartCadena
	xCent=xCent+xPart(i)
	yCent=yCent+yPart(i)
	zCent=zCent+zPart(i)
ENDDO
xPart(0)=xCent/nPartCadena
yPart(0)=yCent/nPartCadena
zPart(0)=zCent/nPartCadena
WRITE(*,*) 'Centro de masas: ', xPart(0), yPart(0), zPart(0)

!----------------------------------------
! SAVE FOR OVITO
OPEN(UNIT=12,FILE='CadenaProle.dat',STATUS='UNKNOWN')
WRITE(12,*) nPartCadena
WRITE(12,*) 'COMMENTS'
DO i=1,nPartCadena
	WRITE(12,*) i,xPart(i),yPart(i),zPart(i),r0
ENDDO
CLOSE(12)

!----------------------------------------
! SAVE incluyendo el centro de masas
!~ OPEN(UNIT=13,FILE='CadenaBurguesaCM.dat',STATUS='UNKNOWN')
!~ WRITE(13,*) 'CM y posiciones de las esferas'
!~ DO i=0,nPartCadena
!~ 	WRITE(13,*) i,xPart(i),yPart(i),zPart(i),r0
!~ ENDDO
!~ CLOSE(13)

!----------------------------------------
end program CadenaProleBien
