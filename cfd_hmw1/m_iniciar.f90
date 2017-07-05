!Módulo para inicializar el problema y sus condiciones de contorno. En este módulo:
!Se calculan los delta_x y delta_y
!Se definen las mallas para u,v y P (todas son distintas)
!Se evalua la condicion inicial para en la malla de u. Si se le agrega la perturbacion en los parámetros del problema, esta se suma
!Contiene la subrutina auxiliares que amplia una matriz para sumarle bordes periodicos
!Como extra, contiene la subrutina de la barra de progreso
module m_iniciar
use m_constantes
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine inicializacion(deltax,deltay,u,v,P)
		integer i,j
		double precision deltax, deltay, auxi
		double precision u(Nx-1,Ny-1), v(Nx-1,Ny-1), P(nx-1,ny-1)
		double precision mallay(Ny-1)
		!Creación de Mallas + Inicialización del Problema
		deltax = (Nx-1)**(-1.0d0)
		deltay = lambda**(-1.0)*(Ny-1)**(-1.0d0)
		do i=1,Ny-1
			mallay(i)=(i-0.5d0)*deltay
			u(:,i)=0.5d0*(1.0d0+tanh(10.0d0*(1.0d0-abs(lambda*mallay(i)-0.5d0)*4.0d0)))
		enddo
		v = 0.0d0		
		call perturbar(u,v,deltax)
		P = 0.0d0
	end subroutine inicializacion
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine perturbar(u,v,deltax)
	integer i,j
	double precision u(Nx-1,Ny-1), malla, deltax, v(Nx-1,Ny-1), auxi
	do j = 1, Ny-1
		do i = 1, Nx-1
			malla = (i-1.0d0)*deltax
			if (perturbacion == 1) then
				u(i,j) = u(i,j) + 0.25d0*u(i,j)*sin(4.0d0*pi*malla)
			endif
			if (perturbacion == 2) then
				u(i,j) = u(i,j)*(1.0d0 + (0.5d0 + 0.5d0*malla)*0.25d0*sin(4.0d0*pi*malla))		!Perturbación Alternativa 2
			endif
			if (perturbacion == 3) then
				u(i,j) = u(i,j)*(1.0d0 + (0.5d0 + 0.5d0*malla)*0.25d0*sin(4.0d0*pi*malla))		!Perturbación Alternativa 3
			if (perturbacion == 3) then
				auxi = (i-0.5d0)*deltax
				v(i,j) = v(i,j) + 0.25d0*0.25d0*sin(2.0d0*pi*auxi)
			endif
			endif
		enddo
	enddo
	end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine auxiliares(u,v,u_aux,v_aux)
		integer i, j
		double precision u(Nx-1,Ny-1), v(Nx-1,Ny-1)
		double precision u_aux(Nx+1,Ny+1), v_aux(Nx+1,Ny+1)
		u_aux(2:Nx,2:Ny) = u
		v_aux(2:Nx,2:Ny) = v

		u_aux(1,2:Ny)			=	u(Nx-1,:)
		u_aux(Nx+1,2:Ny)		=	u(1,:)
		u_aux(2:Nx,1)			=	u(:,Ny-1)
		u_aux(2:Nx,Ny+1)		=	u(:,1) 

		v_aux(1,2:Ny)			=	v(Nx-1,:)
		v_aux(Nx+1,2:Ny)		=	v(1,:)
		v_aux(2:Nx,1)			=	v(:,Ny-1)
		v_aux(2:Nx,Ny+1)		=	v(:,1) 
	end subroutine auxiliares
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine progress(j)
		 implicit none
		 integer(kind=4)::j,k
		 character(len=17)::bar="???% |          |"
		 write(unit=bar(1:3),fmt="(I3)") j+1
		 do k=1, int(j/10)+1
			bar(6+k:6+k)="*"
		 end do
		 call system ('bash red.sh')
		 write(unit=6,fmt="(a1,a17)",advance="no") char(13), bar
		 if (j/=100)then
			flush(unit=6)
		 else
			write(unit=6,fmt=*)
		 endif
		 call system ('bash norm.sh')
		 return
	end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_iniciar
