program main
implicit none
integer, parameter :: N=50 !Número de espacios
integer i,j,cont
double precision tab((N+1)*(N+1),5), dx, dy, R
double precision, parameter :: Lx=3.0d0, Ly=3.0d0
!Calculo de discretización
	dx = Lx/N
	dy = Ly/N
!Formación de Tabla
	cont = 1
	do j=1,N+1
		do i=1,N+1
			tab(cont,1)= 2.0d0*dx*(i-1) - Lx !coordenada x
			tab(cont,2)= 2.0d0*dy*(j-1) - Ly !coordenada y
			R = sqrt(tab(cont,1)**2.0d0 + tab(cont,2)**2.0d0)
			tab(cont,3)= (2.0d0*tab(cont,1)**2.0d0 - R*tab(cont,2)**2.0d0)/(R**5.0d0) !componente x
			tab(cont,4)= (tab(cont,1)*tab(cont,2)*(2.0d0+R))/(R**5.0d0) !componente y
			tab(cont,5)= sqrt(tab(cont,3)**2.0d0 + tab(cont,4)**2.0d0) !magnitud
			tab(cont,3)= 1.5d0*dx*tab(cont,3)/tab(cont,5)
			tab(cont,4)= 1.5d0*dy*tab(cont,4)/tab(cont,5)
			cont = cont + 1
		enddo
	enddo
!Escritura de Tabla
	open(unit=10,file='dipole.txt')
	cont = 1	
	do j=1,N+1
		do i=1,N+1
			write(10,*) tab(cont,1), tab(cont,2), tab(cont,3), tab(cont,4), tab(cont,5)
			cont = cont + 1
		enddo
	enddo
	close(10)
!Ploteo
	open(unit=40, file="style.gnu")
	write(40,*) "set terminal pngcairo size 1200, 600 font 'Verdana, 12'"
	write(40,*) "set output 'dipole.png'"
	write(40,*) "set xrange [",-Lx,":",Lx,"]"
	write(40,*) "set yrange [",-Ly,":",Ly,"]"
	write(40,*) "unset key"
	write(40,*) 'set grid'
	write(40,*) 'set isosamples 2, 2'
	!write(40,*) 'unset xtics'
	!write(40,*) 'unset ytics'
	write(40,*) 'set cbrange [0:1.0]'
	write(40,*) 'unset colorbox'
	write(40,*) "set palette maxcolors 2000"
	write(40,*) "set palette defined ( 0 'blue',500 'cyan',1000 'green',1500 'yellow', 2000 'red')"
	write(40,*) "plot 'dipole.txt' using 1:2:3:4:5 with vectors filled linecolor palette z"
	close(unit=40)
	call system('gnuplot style.gnu')
	call system('rm style.gnu')
	call system('rm dipole.txt')
end program
