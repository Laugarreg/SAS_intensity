! programa para eliminar del sistema
! las esferas que solapan debido al algoritmo de construcción
! 


	program maquillador
	
	implicit none

	REAL(8)	:: Xmin,Xmax, Ymin, Ymax, muda, Zmin, Zmax, Rext, Rp, muda2, muda4, muda5
	REAL(8)	:: x,y,z, dist
	REAL(8)	, ALLOCATABLE :: XYZ(:,:), Racimo(:,:)
	LOGICAL, ALLOCATABLE :: Solapa(:)
	INTEGER i,j,mudaInt,muda3, Npart, Npart2
	CHARACTER*35 nomfich, comentarios

! ----------------------------------------------------------------------
	write (*,*)'Dame el nombre del archivo (formato ovito):'
	read (*,*) nomfich
	
!Abre el fichero y lee el sistema 
	open (21, file=nomfich)
	read (21,*) Npart
	read (21,*) comentarios
	
	ALLOCATE (XYZ(1:Npart,0:5))
	ALLOCATE (Racimo(1:Npart,1:3))
	ALLOCATE (Solapa(1:Npart))
	
	Do i=1,Npart
		read (21,*) XYZ(i,0), XYZ(i,1), XYZ(i,2), XYZ(i,3), XYZ(i,4), XYZ(i,5)
		Solapa(i)=.FALSE.
    enddo

! Sistema cargado
	close (21)

! Etiquetamos las esferas que solapan con un "Solapa(j)=.TRUE.
	Racimo(1,1)=0.00
	Racimo(1,2)=0.00
	Racimo(1,3)=0.00
	
	Do i=1,Npart

		! Contador en pantalla:..................!*
		muda3=100.*i/Npart       	 !*			
		IF (mudaINT.NE.INT(muda3)) then !*		
			mudaINT=INT(muda3)	 !*		
			WRITE(*,*) mudaINT,'%'  !*			
		END IF				 !*		
		! Fin del contador en pantalla...........!*

		if(.NOT.Solapa(i)) then !ignora las que ya sabe que solapan...
		Do j=i+1,Npart			
			dist=sqrt((XYZ(i,1)-XYZ(j,1))**2.+(XYZ(i,2)-XYZ(j,2))**2.+(XYZ(i,3)-XYZ(j,3))**2.)
			IF (dist.LE.0.98) Solapa(j)=.TRUE. !then
			!write (*,*) i,'solapa con',j			
			!Solapa(j)=.TRUE.
			!endif
		enddo
		end if
	enddo

! En Solapa() tenemos almacenadas las etiquetas de las esferas que solapan, al menos, con otra.

! Transferimos a array Racimo solo las que no solapan.

	Npart2=0
	Do i=1, Npart
		If (.NOT.solapa(i)) then
		Npart2=Npart2+1
		Racimo(Npart2,1)=XYZ(i,1)
		Racimo(Npart2,2)=XYZ(i,2)
		Racimo(Npart2,3)=XYZ(i,3)
		end if
	End do
	
	print*, 'N particulas original:',Npart,' y numero que no solapa:',Npart2
	
! guardamos
	
	open(22,file="Maquillado.dat")
	write(22,*) Npart2
	write(22,*) 'COMMENTS'
	Do i=1,Npart2
		write(22,*) XYZ(i,0),Racimo(i,1),Racimo(i,2),Racimo(i,3), XYZ(i,4), XYZ(i,5)
	Enddo
	
	close (22)

stop
end
