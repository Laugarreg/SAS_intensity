program R_ext_giro
implicit none

!*****************************************
!*****************************************
! 		DECLARACION DE VARIABLES
!*****************************************
!*****************************************
INTEGER					:: Ncentros, i, Npart, j, dummy, aux_int
REAL(8)					:: dist, dist_0, dummyreal, aux
REAL(8)					:: xrel, yrel, zrel, Lcaja
REAL(8), ALLOCATABLE	:: Centros(:,:), Particulas(:,:), CM_ovillo(:,:)
REAL(8), ALLOCATABLE	:: R_ext(:)
REAL(8), ALLOCATABLE	:: R_giro(:)
INTEGER, ALLOCATABLE	:: Cluster(:)
CHARACTER(len=50)		:: comentarios, nomfich


!*****************************************
!*****************************************
! 			INICIALIZACION
!*****************************************
!*****************************************
! Leer centros atractores
OPEN(12,FILE='Centros.dat',STATUS='OLD')
	READ(12,*) Ncentros
	ALLOCATE(Centros(Ncentros,3))
	READ(12,*) comentarios
	DO i=1,Ncentros
		READ(12,*) dummy, Centros(i,:)
	ENDDO
CLOSE(12)
WRITE(*,*) 'Centros atractores leidos'

! Leer sistema particulas
Lcaja = 0.0d0
WRITE(*,*) 'Dame el nombre del archivo (formato ovito):'
READ(*,*) nomfich
OPEN(13,FILE=nomfich,STATUS='OLD')
	READ(13,*) Npart
	ALLOCATE(Particulas(Npart,3))
	READ(13,*) comentarios
	DO i=1,Npart
		READ(13,*) dummyreal, Particulas(i,:)
		Lcaja = MAX(MAXVAL(Particulas(i,:)),Lcaja)
	ENDDO
CLOSE(13)
Lcaja = 2.0d0*Lcaja
WRITE(*,*) 'Particulas del sistema leidas'


!*****************************************
!*****************************************
! 			CALCULO PRINCIPAL
!*****************************************
!*****************************************
! Centro atractor mas cercano
ALLOCATE(Cluster(Npart))
Cluster = 0
DO i=1,Npart
	dist_0 = 1000000000.0d0
	DO j=1,Ncentros
	
		! Distancia al centro atractor
		xrel = Particulas(i,1) - Centros(j,1)
		yrel = Particulas(i,2) - Centros(j,2)
		zrel = Particulas(i,3) - Centros(j,3)
		
		! Correcion PBC
		xrel = xrel - Lcaja*DNINT(xrel/Lcaja)
		yrel = yrel - Lcaja*DNINT(yrel/Lcaja)
		zrel = zrel - Lcaja*DNINT(zrel/Lcaja)
		
		dist = DSQRT( xrel**2.0d0 + yrel**2.0d0 + zrel**2.0d0 )
		
		! Actualizacion del cluster más cercano
		IF (dist.LT.dist_0) THEN
			Cluster(i) = j
			dist_0 = dist
		ENDIF		
	ENDDO
ENDDO
WRITE(*,*) 'Centro atractor mas cercano calculados'

! Reposicionamiento de las particulas
DO i=1,Ncentros
	DO j=1,Npart
		IF (Cluster(j).EQ.i) THEN
			xrel = Particulas(j,1) - Centros(i,1)
			yrel = Particulas(j,2) - Centros(i,2)
			zrel = Particulas(j,3) - Centros(i,3)
			
			Particulas(j,1) = Particulas(j,1) - Lcaja*DNINT(xrel/Lcaja)
			Particulas(j,2) = Particulas(j,2) - Lcaja*DNINT(yrel/Lcaja)
			Particulas(j,3) = Particulas(j,3) - Lcaja*DNINT(zrel/Lcaja)
		ENDIF	
	ENDDO
ENDDO
WRITE(*,*) 'Reposicionamiento realizado'

! Centro de masa de cada ovillo
ALLOCATE(CM_ovillo(Ncentros,3))
CM_ovillo = 0.0d0
DO  i=1,Ncentros
	aux_int = 0
	DO j=1,Npart
		
		IF (Cluster(j).EQ.i) THEN		
		
			CM_ovillo(i,1) = CM_ovillo(i,1) + Particulas(j,1)
			CM_ovillo(i,2) = CM_ovillo(i,2) + Particulas(j,2)
			CM_ovillo(i,3) = CM_ovillo(i,3) + Particulas(j,3)
			
			aux_int = aux_int + 1
		ENDIF
		
	ENDDO
	CM_ovillo(i,1) = CM_ovillo(i,1) / aux_int
	CM_ovillo(i,2) = CM_ovillo(i,2) / aux_int
	CM_ovillo(i,3) = CM_ovillo(i,3) / aux_int
ENDDO
WRITE(*,*) 'Centro de masa de los ovillos calculados'

! Radio externo y radio de giro de cada cluster
ALLOCATE(R_ext(Ncentros))
ALLOCATE(R_giro(Ncentros))
R_ext = 0.0d0
R_giro = 0.0d0
dist = 0.0d0

DO i=1,Ncentros
	! Inicializar acumuladores a cero
	aux = 0.0d0
	aux_int = 0
	
	DO j=1,Npart
		IF (Cluster(j).EQ.i) THEN
		
			! Distancia relativa por componentes al CM del ovillo
			xrel = Particulas(j,1) - CM_ovillo(i,1)
			yrel = Particulas(j,2) - CM_ovillo(i,2)
			zrel = Particulas(j,3) - CM_ovillo(i,3)
		
			! Correccion PBC
			xrel = xrel - Lcaja*DNINT(xrel/Lcaja)
			yrel = yrel - Lcaja*DNINT(yrel/Lcaja)
			zrel = zrel - Lcaja*DNINT(zrel/Lcaja)
		
			! Distancia
			dist = DSQRT( xrel**2.0d0 + yrel**2.0d0 + zrel**2.0d0 )
			
			! Aceptacion o rechazo del nuevo R_ext
			IF (R_ext(i).LT.dist) R_ext(i) = dist
			
			! Calculo R_giro
			aux = aux + xrel**2.0d0 + yrel**2.0d0 + zrel**2.0d0
			aux_int = aux_int + 1
		ENDIF
	ENDDO ! fin bucle j=1,Npart
	
	R_giro(i) = DSQRT(aux/aux_int)
ENDDO ! Fin bucle i=1,Ncentros
WRITE(*,*) 'Rext y R_giro calculado'


!*****************************************
!*****************************************
! 			ESCRITURA EN FICHERO
!*****************************************
!*****************************************
! Resultados de los CM de los ovillos
OPEN(15,FILE='CM_ovillos.dat')
	WRITE(15,*) Ncentros
	WRITE(15,*) 'Centro de masas de los ovillos'
	DO i=1,Ncentros
		WRITE(15,*) i, CM_ovillo(i,:)
	ENDDO
CLOSE(15)
WRITE(*,*) 'CM de los ovillos salvados'

! Resultados del radio externo y el radio de giro
OPEN(14,FILE='R_exterior_giro.dat')
	WRITE(14,*) 'Radio externo:'
	DO i=1,Ncentros
		WRITE(14,*) R_ext(i)
	ENDDO
	WRITE(14,*) 'R_ext medio = ', SUM(R_ext)/Ncentros
	WRITE(14,*) ''
	WRITE(14,*) 'Radio giro:'
	DO i=1,Ncentros
		WRITE(14,*) R_giro(i)
	ENDDO
	WRITE(14,*) 'R_giro medio = ', SUM(R_giro)/Ncentros
CLOSE(14)
WRITE(*,*) 'R_ext y R_giro salvados'

! Sistema reposicionado
OPEN(15,FILE='Reposicionamiento.dat')
	WRITE(15,*) Npart
	WRITE(15,*) 'Sistema reposicionado'
	DO i=1,Npart
		WRITE(15,*) i, Particulas(i,:), 0.5
	ENDDO
CLOSE(15)
WRITE(*,*) 'Reposicionamiento salvado'


WRITE(*,*) ''
WRITE(*,*) 'Done!'
end program R_ext_giro
