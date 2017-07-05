program main
implicit none
integer, parameter :: N=1000 !Número de espacios
integer i
double precision x(N+1), f(N+1), dx
double precision, parameter :: Lf=310.0, Li=290.0
dx=(Lf-Li)/N
do i=0,N
	x(i+1)=Li+dx*i
	f(i+1)=x(i+1)-9.5d0*10.0d0**(-9.0d0)*x(i+1)**4.0d0-223.4d0
enddo
!Escritura de Tabla
	open(unit=10,file='plot_eq.txt')
	do i=1,N+1
		write(10,*) x(i), f(i)
	enddo
	close(10)
!Ploteo
	open(unit=40, file="style.gnu")
	write(40,*) "set terminal pngcairo font 'Verdana, 10'"
	write(40,*) "set output 'plot_eq.png'"
	write(40,*) "unset key"
	write(40,*) 'set grid'
	write(40,*) "plot 'plot_eq.txt' with lines ls 6 lw 2"
	close(unit=40)
	call system('gnuplot style.gnu')
	call system('rm style.gnu')
	call system('rm plot_eq.txt')
end program
