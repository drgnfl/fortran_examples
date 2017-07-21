!Este programa genera un Contour Plot pero con Nodos y vectores del gradiente
!Ideal para tareas donde se utilizen diferencias finitas
program main
implicit none
double precision, parameter :: L = 1.0d0
integer i,j
integer, parameter :: N = 41 !Número de Nodos
double precision dx,y,x, T(N,N), qx, qy, magn
!main
dx = L/(N-1) !Espacio entre nodos
!Escritura de Tabla
!El formato es pos_x pos_y valor
!Notar el espacio entre cada fila (grided data)
open(unit=10,file='heatt.txt')
do i=1,N
	y=dx*(i-1)
	do j=1,N
		x=dx*(j-1)
		T(i,j) = sinh(3.14*(L-x)/L)**0.5*sin(3.14*y/L)*y
		write(10,*) x,y, T(i,j)
	enddo
	write(10,*) 
enddo
close(10)
!Vectores del gradiente
open(unit=20,file='heattvec.txt')
do i=2,N-1
	y = dx*(i-1)
	do j=2,N-1
		x = dx*(j-1)
		qx = (-T(i,j-1)+T(i,j+1))/dx !Diferencias Centradas
		qy = (-T(i-1,j)+T(i+1,j))/dx !Diferencias Centradas
		magn = sqrt(qx**2+qy**2)
		write(20,*) x,y,0.9*dx*qx/18,0.7*dx*qy/18, magn !Normalizados por el valor mas grande (18)
	enddo
enddo
close(20)
!Ploteo llamando a gnuplot con círculos
open(unit=40, file="style.gnu")
write(40,*) "set terminal pngcairo enhanced font 'Verdana,10'"
write(40,*) "set output 'heatt.png'"
write(40,*) "unset key"
write(40,*) 'set grid'
write(40,*) "set palette defined(1 'red',2 'yellow')"
write(40,*) "set palette defined(1 'blue',2 'green',3 'yellow',4 'red')"
write(40,*) 'set style fill solid border lt 16'
write(40,*) 'set style circle radius graph 0.01'
write(40,*) "plot 'heatt.txt' u 1:2:3 w circles palette"
write(40,*) "set output 'heattvec.png'"
write(40,*) "plot 'heattvec.txt' u 1:2:3:4:5 w vectors filled head palette"
close(unit=40)
call system('gnuplot style.gnu')
!Ploteo llamando a gnuplot con círculos
open(unit=50, file="style2.gnu")
write(50,*) "set terminal pngcairo enhanced font 'Verdana,10'"
write(50,*) "set output 'heatt2.png'"
write(50,*) "unset key"
write(50,*) "set xrange [0:1]"
write(50,*) "set yrange [0:1]"
write(50,*) "set palette defined(1 'blue',2 'green',3 'yellow',4 'red')"
write(50,*) "set contour"
write(50,*) "unset surface"
write(50,*) "set cntrparam levels auto 20"
write(50,*) "set view map"
write(50,*) "unset clabel"
write(50,*) "set table 'cont.dat'"
write(50,*) "splot 'heatt.txt' u 1:2:3 w lines"
write(50,*) "unset table"
write(50,*) "unset contour"
write(50,*) "plot 'heatt.txt' u 1:2:3 w image, 'cont.dat' w l lt -1 lw 1.5"
close(unit=50)
call system('gnuplot style2.gnu')
!Se borran los archivos auxiliares
call system('rm style.gnu')
call system('rm style2.gnu')
call system('rm heatt.txt')
call system('rm cont.dat')
call system('rm heattvec.txt')
end program