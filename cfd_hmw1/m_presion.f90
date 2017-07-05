!Módulo para hacer las operaciones para despejar la presión en n+1. En este módulo:
!Se calcula la divergencia del campo de velocidad v*. 
!Se genera el lado derecho (Q) de la ecuación de Poisson para Phi. 
!Se aplica la transformada de Fourier para este lado derecho.
!Se resuelve el sistema tridiagonal formado por la serie de edos que nacen de la ec. de poisson
!Se aplica la transformada inversa para generar la solución.
!Además se incluye la subrutina para calcular vorticidad
module m_presion
use m_constantes
use m_tridiagonal
use m_iniciar
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Calcula la vorticidad de un campo de velocidades en cada punto de la malla de presión
	subroutine vorticidad(u_aux,v_aux,deltax,deltay,vort)
		integer i,j
		double precision u_aux(Nx+1,Ny+1),v_aux(Nx+1,Ny+1),vort(Nx-1,Ny-1), deltax, deltay
		vort = 0.0d0
		do j = 2,Ny
			do i = 2,Nx
					vort(i-1,j-1) = (v_aux(i+1,j)+v_aux(i+1,j+1)-v_aux(i-1,j)-v_aux(i-1,j+1))*0.25d0*deltax**(-1.0d0)  
					vort(i-1,j-1) = vort(i-1,j-1) - (u_aux(i,j+1)+u_aux(i+1,j+1)-u_aux(i,j-1)-u_aux(i+1,j-1))*0.25d0*deltay**(-1.0d0)  
			enddo
		enddo
	end subroutine vorticidad
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Calcula la divergencia de un campo de velocidades en cada punto de la malla de presión
	subroutine divergencia(u,v,deltax,deltay,div)
		integer i,j
		double precision u(Nx-1,Ny-1), v(Nx-1,Ny-1)
		double precision u_aux(Nx+1,Ny+1),v_aux(Nx+1,Ny+1),div(Nx-1,Ny-1),deltax, deltay
		div = 0.0d0
		call auxiliares(u,v,u_aux,v_aux)
		do j = 2,Ny
			do i = 2,Nx
				div(i-1,j-1) = ( u_aux(i+1,j)-u_aux(i,j) )*deltax**(-1.0d0) + ( v_aux(i,j+1)-v_aux(i,j) )*deltay**(-1.0d0)
			enddo
		enddo
	end subroutine divergencia
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Resuelve la ecuación de Poisson con la transformada de Fourier para cada punto de la malla
	subroutine edo_fourier(deltax,deltay,Phi,Q)
		integer l
		double precision deltax, deltay, k(Nx-1)
		double precision Q(Nx-1,Ny-1), div(Nx-1,Ny-1)
		double precision Phi(Nx-1,Ny-1)
		double precision alpha, beta
		complex*16 hatPhi(Nx-1,Ny-1), hatQ(Nx-1,Ny-1)
		!Resolvemos la edo para cada l como un sistema tridiagonal con cond. periódicas.
		do l = 1, Ny-1
			!Aplicamos la Transformada de Fourier
			call trans_fourier(Q(:,l),hatQ(:,l))
		enddo
		do l = 1, Nx-1
			!Resolvemos el sistema
			k(l) = 2.0d0*(cos((2.0d0*pi)*(l-1)*(Nx-1)**(-1.0d0))-1.0d0)*deltax**(-2.0d0)
			alpha = (k(l) - 2.0d0*deltay**(-2.0d0)) !Elementos de la diag principal
			beta = deltay**(-2.0d0) !Sub diagonal y esquinas
			if (l==1) then
				alpha = (k(l) - 2.0d0*deltay**(-2.0d0) + 1.0d0) !Anclamos Phi=0.0 en el nodo (1,1) arbitrariamente
			endif
			call tridiag_complx(alpha,beta,hatQ(l,:),Ny-1,hatPhi(l,:)) !Resuelve el sistema
		enddo
		do l = 1,Ny-1
			!Transformada Inversa para Phi
			call trans_in_fourier(hatPhi(:,l),Phi(:,l))
		enddo
	end subroutine edo_fourier
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Transformada discreta de Fourier
	subroutine trans_fourier(Q,hatQ)
		integer m,l
		double precision Q(Nx-1)
		complex*16 hatQ(Nx-1), suma, aux
		complex*16, parameter :: i=(0.0d0,1.0d0)
		do m=1,Nx-1
			suma = 0.0d0
			do l = 1,Nx-1
				aux = Q(l)*exp(-i*2.0d0*pi*(m-1)*(l-1)*(Nx-1)**(-1.0d0))
				suma = suma + aux
			enddo
			hatQ(m) = suma*(Nx-1)**(-1.0d0)
		enddo
	end subroutine trans_fourier
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Transformada inversa discreta de Fourier
	subroutine trans_in_fourier(hatPhi,Phi)
		integer m,l
		double precision Phi(Nx-1)
		complex*16 hatPhi(Nx-1), suma, aux
		complex*16, parameter :: i=(0.0d0,1.0d0)
		do m=1,Nx-1
			suma = 0.0d0
			do l = 1,Nx-1
				aux = hatPhi(l)*exp(i*2.0d0*pi*(m-1)*(l-1)*(Nx-1)**(-1.0d0))
				suma = suma + aux
			enddo
			Phi(m) = suma
		enddo		
	end subroutine trans_in_fourier
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Validación del solver de Poisson
!Q		=	-4*Pi**2*( sin(2*pi*x) + sin(2*pi*y) )
!Phi	=	sin(2*pi*x) + sin(2*pi*y)
	subroutine validar_poisson(deltax,deltay)
		integer i,j
		double precision deltax, deltay, aux_x, aux_y, error
		double precision Phi(Nx-1,Ny-1), teo(Nx-1,Ny-1)
		double precision Q(Nx-1,Ny-1)
		do i = 1, Nx-1
			aux_x = (i-1.0d0)*deltax
			do j= 1, Ny-1
				aux_y = (j-1.0d0)*deltay
				Q(i,j) = -4.0d0*pi*pi*(sin(2.0d0*pi*aux_x) + sin(2.0d0*pi*aux_y))
				teo(i,j) = sin(2.0d0*pi*aux_x) + sin(2.0d0*pi*aux_y)
			enddo
		enddo

		call edo_fourier(deltax,deltay,Phi,Q)

		open(unit=10,file='plot/val_teo.txt')
		do j = 1, Ny - 1
			write(10,*) (teo(i,j),i=1,Nx-1)
		enddo
		close(10)
		open(unit=40, file="style.gnu")
		write(40,*) "set pm3d map interpolate 0,0"
		write(40,*) "set terminal pngcairo font 'Verdana,12"
		write(40,*) "set palette maxcolors 200"
		write(40,*) "set palette defined (0 'blue', 100 'white', 200 'red')"
		write(40,*) "unset key"
		write(40,*) "set xrange [0:1]"
		write(40,*) "set yrange [0:1]"
		write(40,*) "set cbrange [-2:2]"
		write(40,*) "set output 'plot/val_teo.png'"
		write(40,*) "set title 'Solución Teórica'"
		write(40,*) "splot 'plot/val_teo.txt' matrix using ($1/",(Nx-2),"):($2/",(Ny-2),"):3"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu plot/val_teo.txt')

		open(unit=20,file='plot/val_sol.txt')
		do j = 1, Ny - 1
			write(20,*) (Phi(i,j),i=1,Nx-1)
		enddo
		close(20)
		open(unit=40, file="style.gnu")
		write(40,*) "set pm3d map interpolate 0,0"
		write(40,*) "set terminal pngcairo font 'Verdana,12"
		write(40,*) "set palette maxcolors 200"
		write(40,*) "set palette defined (0 'blue', 100 'white', 200 'red')"
		write(40,*) "unset key"
		write(40,*) "set xrange [0:1]"
		write(40,*) "set yrange [0:1]"
		write(40,*) "set cbrange [-2:2]"
		write(40,*) "set output 'plot/val_sol.png'"
		write(40,*) "set title 'Solución Aproximada'"
		write(40,*) "splot 'plot/val_sol.txt' matrix using ($1/",(Nx-2),"):($2/",(Ny-2),"):3"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu plot/val_sol.txt')

		write(*,*) 'Se ha ejecutado la subrutina para validad el solver de Poisson'
		write(*,*) 
	end subroutine validar_poisson
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_presion
