program Nstep_Max
implicit none
! *********************************************
! *********************************************
INTEGER		:: Nstep, i, j
REAL(8) 	:: energia1, energia2
REAL(8)		:: dif_rel
! *********************************************
! *********************************************
OPEN(UNIT=12,FILE='Energia.dat',STATUS='OLD')
	READ(12,*) Nstep, energia1
	DO i=1,1000000000
		READ(12,*) Nstep, energia2
		dif_rel = abs(energia2-energia1)*100.0d0	
		IF (dif_rel.LT.1.0d0) GOTO 100
		energia1 = energia2		
	ENDDO
CLOSE(12)
WRITE(*,*) 'No se ha alcanzado el equilibrio'
GOTO 101

100 CONTINUE
OPEN(UNIT=13,FILE='EnergiaNstepmax.dat',STATUS='UNKNOWN')
	WRITE(13,*) 'El paso en el que se ha alcanzado el equilibrio es: ', Nstep
CLOSE(13)

101 CONTINUE
! *********************************************
! *********************************************
end program Nstep_Max
