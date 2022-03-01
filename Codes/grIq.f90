program grIq
implicit none
!****************************************
!****************************************
!
! Script para calcular la función de correlación de pares, g(r),  y la
! curva de intensidad de scattering, I(q), a partir de las coordenadas 
! de todas las esferas de un sistema.
! 
! INPUTS: 
! 	Fichero: gr.in, que contiene, por filas: el nombre del fichero 
! 			 con las coordenadas de las esferas, el radio de la esfera
!			 elemental, y el dr empleado en la discretización de r.
!
! OUTPUTS:
! 	Fichero: grIq.out, que contiene información de los parámetros 
!			 empleados en la simulación.
! 	Fichero: gr.dat, que contiene, por columnas, los valores del vector 	
!			 r y el vector gr.
! 	Fichero: Iq.dat, que contiene, por columnas, los valores del vector 	
!			 q, Pq, Sq, Iq.
!
!****************************************
!****************************************

!----------------------------------------
! DECLARACION DE VARIABLES
REAL(8), PARAMETER :: pi=3.1415926535897932384626
INTEGER :: npartRacimo, rlong1, rlong2, i, j, k
REAL(8) :: muda, L1, L2, Lcaja, rmax, rho, dist, bucle 
CHARACTER(len=25) :: nombrefich, comentarios
REAL(8), DIMENSION(:), ALLOCATABLE :: xRacimo, yRacimo, zRacimo
REAL(8), DIMENSION(:), ALLOCATABLE :: gr, rSist, rSist2, grMayus! gr del sistema 
REAL(8), DIMENSION(:), ALLOCATABLE :: xrel, yrel, zrel
REAL(8) :: r0, dr ! Radio esfera elemental (r0) e intervalo para recorrer r (dr)
REAL(8) :: Suma
REAL(8) :: muda1, muda2, muda3 !Variables mudas
CHARACTER(len=10) fecha,hora   !Registro
REAL(8) :: qmin, qmax, dq, g0, c, num, den, Acum  !Calculo intensidad
INTEGER :: qlong								  !Calculo intensidad
REAL(8), DIMENSION(:), ALLOCATABLE :: q, Pq, Sq, Iq, numVec, denVec, Int0, Integrando	!Calculo intensidad

!----------------------------------------
! 1. LECTURA DE PARAMETROS DE SIMULACION Y LAS COORDENADAS x,y,z DE LAS ESFERAS DEL SISTEMA
!----------------------------------------
! Fecha y hora
CALL DATE_AND_TIME(DATE = fecha,TIME = hora)

! Cargamos los parámetros de ejecución a partir de un archivo
OPEN(UNIT=11, FILE='grIq.in',STATUS='OLD')
	READ(11,*) nombrefich  !nombre del archivo que contiene las coordenadas de las particulas del sistema
	READ(11,*) r0  		! Radio de la esfera elemental
	READ(11,*) dr  		! dr que queremos utilizar
	READ(11,*) qmin  	! q minimo para I(q)
	READ(11,*) qmax 	! q maximo para I(q)
	READ(11,*) dq  		! dq para discretizar q entre qmin y qmax
CLOSE(11)

! Abrimos el archivo de registro
OPEN(UNIT=12, FILE='grIq.out',STATUS='UNKNOWN')
	WRITE(12,*) 'Ejecutado el día',fecha,'  a las',hora
	WRITE(12,*) '------------------------------------------------------'
	WRITE(12,*) ''
	WRITE(12,*) 'r0=',r0
	WRITE(12,*) 'dr=',dr
	WRITE(12,*) 'Sistema:',nombrefich

! Intentamos primero cargar los datos de un archivo formato .xyz de 3 columnas
OPEN(UNIT=13,FILE=nombrefich,STATUS='OLD')
	npartRacimo=0
	DO i=1,100000000  !max. 100.000.000 partículas
		READ(13,*,end=22,err=23) muda, muda1, muda2, muda3
		npartRacimo=npartRacimo+1
	ENDDO
	
! Si se acaban las líneas del fichero, salta aquí
22  REWIND (13)   !Reseteo de lectura del fichero de entrada
	ALLOCATE(xRacimo(npartRacimo))  !Dimensionar los vectores que contienen las coordenadas de las particulas
	ALLOCATE(yRacimo(npartRacimo))
	ALLOCATE(zRacimo(npartRacimo))

	DO i=1,npartRacimo   !Lectura de las coordenadas del fichero de entrada
		READ(13,*) muda, xRacimo(i), yRacimo(i), zRacimo(i)
	ENDDO
CLOSE(13)
WRITE(12,*) 'Archivo de sistema cargado en formato xyz.dat' !Fichero de registro
GOTO 24		!Si se ha leido en fomato .xyz, no se intenta leer en formato OVITO

! Si el fichero estaba en formato OVITO, error en el READ y salta aquí
23 REWIND (13)	!Reseteo de lectura del fichero de entrada
	READ(13,*) npartRacimo		!Lee el numero de particulas
	READ(13,*) comentarios		!Lee linea de comentarios
	ALLOCATE(xRacimo(npartRacimo))	!Dimensionaliza los arrays que contienen las coordenadas de las particulas
	ALLOCATE(yRacimo(npartRacimo))
	ALLOCATE(zRacimo(npartRacimo))
	DO i=1,npartRacimo	!Lectura de las coordenadas del fichero de entrada
		READ(13,*) muda, xRacimo(i), yRacimo(i), zRacimo(i), muda1
	ENDDO
CLOSE(13)
WRITE(12,*) 'Archivo de sistema cargado en formato OVITO'	!Fichero de registro
GOTO 24

24 CONTINUE


!----------------------------------------
! 2. PARÁMETROS DE SIMULACIÓN
!----------------------------------------

! Dimensiones de la caja de simulacion (sumando 2r0 para que el calculo de la densidad sea mas similar al obtenido por MC)
L1=MAX(MAXVAL(xRacimo),MAXVAL(yRacimo),MAXVAL(zRacimo))
L2=MIN(MINVAL(xRacimo),MINVAL(yRacimo),MINVAL(zRacimo))
	!esto es un poco a lo bruto: habría que afinar más para ahorrar cálculos.
Lcaja= 2.0*MAX(ABS(L1),ABS(L2)) + 2.0*r0
WRITE(12,*) ''
WRITE(12,*) 'SISTEMA'
WRITE(12,*) 'Numero de particulas:',npartRacimo
WRITE(12,*) 'Lmin=',L2,', Lmax=',L1,', Lcaja=',Lcaja

! Densidad numérica
rho=npartRacimo/(Lcaja**dble(3))
WRITE(12,*) 'Densidad=',rho
 
! Valores que toma la variable r
! puede ser desde 0 hasta la diagonal del cubo:
rmax=INT(0.5*Lcaja)		!Caida del plateau g(r)=1 en racimos
rlong1 = INT(Lcaja/dr)	!Dimension de g(r) para ver la caida a cero
rlong2 = INT(rmax/dr)	!Dimension de g(r) cortada en el plateau para S(q)
WRITE(12,*) 'r_max:', rmax
WRITE(12,*) 'Numero de componentes de g(r):', rlong1
WRITE(12,*) 'Numero de componentes de g(r) cortada:', rlong2


!----------------------------------------
! 3. CALCULO DE g(r)
!----------------------------------------
ALLOCATE(xrel(npartRacimo))		!Dimensionalizacion para distancias
ALLOCATE(yrel(npartRacimo))
ALLOCATE(zrel(npartRacimo))


! Inicialización de g(r) a cero
ALLOCATE(gr(rlong1)) 
DO i=1,rlong1
gr(i)=0
ENDDO

! Registro de las particulas que violan el contacto duro tras aplicar las PBC
WRITE(12,*) ''
WRITE(12,*) '---- PARTICULAS QUE VIOLAN EL CONTACTO DURO ----'
WRITE(12,*) ''


! Bucle para el calculo de g(r) 
DO i=1, npartRacimo-1  		!i recorre todas las particulas menos la ultima
	DO j=i+1, npartRacimo   !j recorre desde la particula i+1 hasta la ultima
		xrel(j)=xRacimo(j)-xRacimo(i)  !Coordenadas relativas ij
		yrel(j)=yRacimo(j)-yRacimo(i)
		zrel(j)=zRacimo(j)-zRacimo(i)
		xrel(j)= xrel(j) - Lcaja*DNINT(xrel(j)/Lcaja)   !Correcion de las coordenadas relativas ij mediante PBC (MIC)
		yrel(j)= yrel(j) - Lcaja*DNINT(yrel(j)/Lcaja)
		zrel(j)= zrel(j) - Lcaja*DNINT(zrel(j)/Lcaja)

		!Distancia ij corregida mediante PBC (MIC)
		dist = sqrt((xrel(j))**2.0 + (yrel(j))**2.0 + (zrel(j))**2.0)

		IF (dist.LT.0.98) WRITE(12,*) i,j,dist   !Registro de particulas que violan el contacto duro
		k=INT(dist/dr) + 1   !Calculo de la componente donde debe registrarse esta distancia (r-dr/2, r+dr/2)
		gr(k)=gr(k)+1	  !Sumar la contribucion a g(r)
	ENDDO

	!Contador del bucle i por pantalla
	bucle=dble(i)/dble(npartRacimo)*100 
	WRITE(*,*) 'Particula: ', i, '; Porcentaje calculado: ', bucle, '%'
ENDDO

! Calculo de g(r) y escritura de resultados en fichero
ALLOCATE(rSist(rlong1))
Suma=0.0
OPEN(UNIT=21,FILE="gr.dat")
DO i=1,rlong1
	rSist(i)=(dble(i)-0.5)*dr
	gr(i)=2.0*gr(i)/(4*pi*rSist(i)**2.0*dr*rho*npartRacimo)  !Se multiplica por 2 para considerar las parejas i,j 	!f(r)/rho de HASMY
	WRITE(21,*) rSist(i), gr(i)
	Suma=Suma+gr(i)
ENDDO
CLOSE(21)
WRITE (*,*)'Calculo de g(r) terminado.'

! Comprobación de la normalizacion de g(r) (HASMY) --> Pendiente de comprobacion




!----------------------------------------
! 4. CALCULO DE P(q), S(q), I(q)
!----------------------------------------
! Parametros para definir el rango de la variable q
!~ qmin=0.0010d0
!~ qmax=25.0d0
!~ dq=1e-4
qlong=NINT((qmax-qmin)/dble(dq))

! Generacion de las componentes del vector q
ALLOCATE(q(qlong))
DO i=1,qlong
   q(i)=qmin+(i-1)*dq
ENDDO

! Valor g0
c=(pi/6)*(npartRacimo/Lcaja**dble(3))
ALLOCATE(numVec(rlong2))
ALLOCATE(denVec(rlong2))
DO i=1,rlong2
   numVec(i)=4*pi*dr*rSist(i)**dble(2)*gr(i)
   denVec(i)=4*pi*dr*rSist(i)**dble(2)
ENDDO
num=SUM(numVec)
den=SUM(denVec)
g0=(num+pi/(6*c))/den

! Escritura de parametros en el fichero de salida
WRITE(12,*) '---- PARAMETROS Iq ----'
   WRITE(12,*) 'qmin: ', qmin 
   WRITE(12,*) 'qmax: ', qmax 
   WRITE(12,*) 'dq: ', dq
   WRITE(12,*) 'qlong: ', qlong
   WRITE(12,*) 'c: ', c
   WRITE(12,*) 'g0: ', g0
CLOSE(12)

! Dimensionalización de Pq, Sq, Iq
ALLOCATE(rSist2(rlong2))
ALLOCATE(grMayus(rlong2))
ALLOCATE(Int0(rlong2))
ALLOCATE(Integrando(rlong2))
ALLOCATE(Pq(qlong))
ALLOCATE(Sq(qlong))
ALLOCATE(Iq(qlong))

! Definicion de los vectores rSist2 y grMayus
DO i=1,rlong2
	rSist2(i)=rSist(i)
	grMayus(i)=gr(i)-g0
ENDDO

! Calculo de Pq, Sq, Iq
DO i=1,qlong
!~ ! Calculo directo de Hasmy para S(q)
!~ 	sumaSq=0.0
!~ 	DO j=1,npartRacimo
!~ 		DO k=j,npartRacimo
!~ 			dist=sqrt((xRacimo(j)-xRacimo(k))**2.0+(yRacimo(j)-yRacimo(k))**2.0+(zRacimo(j)-zRacimo(k))**2.0)
!~ 			sumaSq=sumaSq + sin(q(i)*dist)/(q(i)*dist)
!~ 		ENDDO
!~ 	ENDDO
!~ 	Sq(i)=1.0+sumaSq/npartRacimo

	!Cálculo de S(q) habitual a partir de g(r=
	Int0=sin(q(i)*rSist2)/(q(i)*rSist2)
    Integrando=4*pi*rSist2**2.0*dr*Int0*grMayus
    Acum=SUM(Integrando)
   
    Sq(i)=1 + 6.0*c*Acum/(pi)
	Pq(i)=(24.0*(sin(r0*q(i))-r0*q(i)*cos(r0*q(i)))/(q(i)**3.0))**2.0
	Iq(i)=Pq(i)*Sq(i)
	
	bucle=dble(i)/dble(qlong)*100 
	WRITE(*,*) 'q(i): ', i, '; Porcentaje calculado: ', bucle, '%'
ENDDO
WRITE (*,*)'Calculo de I(q) terminado.'

! Escritura de resultados en fichero
OPEN(UNIT=17,FILE='Iq.dat',STATUS='UNKNOWN')
   DO i=1,qlong
      WRITE(17,*) q(i), Pq(i), Sq(i), Iq(i)
   ENDDO
CLOSE(17)


!----------------------------------------
WRITE(*,*) ''
WRITE(*,*) '--FIN DEL PROGRAMA--'
end program grIq
