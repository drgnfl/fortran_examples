!Módulo que calcula el lado derecho de la ecuación de cantidad de movimiento. Acá se hace:
!Se calcula el termino advectivo tanto para u y para v, y al final se ponderan como Adams-Bashford. Nota: La primera iteración se hace con un Euler Explícito (al realizar u_ant = u y v_ant = v en la inicializacion).
!Se calcula el término difusivo.
!Se calculan los gradientes de presión para P en el tiempo n.
!Se suman todos al final.
module m_rhs
use m_constantes
implicit none
contains      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine adveccion(u_aux,v_aux,deltax,deltay,adv_x,adv_y)
	integer i, j
	double precision adv_x(Nx-1,Ny-1),adv_y(Nx-1,Ny-1),u_aux(Nx+1,Ny+1),v_aux(Nx+1,Ny+1)
	double precision aux1,aux2, deltax, deltay
	!Recorremos cada punto de la malla
	adv_x=0.0d0
	adv_y=0.0d0
	do j=2,Ny
		do i =2,Nx
			adv_x(i-1,j-1) = -0.5d0*(u_aux(i+1,j)**2.0d0-u_aux(i-1,j)**2.0d0)*deltax**(-1.0d0)
			aux1 = 0.25d0*(u_aux(i,j)+u_aux(i,j+1))*(v_aux(i-1,j+1)+v_aux(i,j+1))
			aux2 = 0.25d0*(u_aux(i,j)+u_aux(i,j-1))*(v_aux(i-1,j)+v_aux(i,j))
			adv_x(i-1,j-1) = adv_x(i-1,j-1) - (aux1 - aux2)*deltay**(-1.0d0)

			adv_y(i-1,j-1) = -0.5d0*(v_aux(i,j+1)**2.0d0-v_aux(i,j-1)**2.0d0)*deltay**(-1.0d0)
			aux1 = 0.25d0*(u_aux(i+1,j)+u_aux(i+1,j-1))*(v_aux(i,j)+v_aux(i+1,j))
			aux2 = 0.25d0*(u_aux(i,j)+u_aux(i,j-1))*(v_aux(i,j)+v_aux(i-1,j))
			adv_y(i-1,j-1) = adv_y(i-1,j-1) - (aux1 - aux2)*deltax**(-1.0d0)
		enddo
	enddo
	end subroutine adveccion
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine laplac(u_aux,deltax,deltay,dif)
	integer i,j 
	double precision u_aux(Nx+1,Ny+1), deltax, deltay
	double precision dif(Nx-1,Ny-1)
	!Recorremos cada punto de la malla
	dif=0.0d0
	do j=2,Ny
		do i = 2,Nx
			dif(i-1,j-1) = ( u_aux(i-1,j)-2.0d0*u_aux(i,j)+u_aux(i+1,j) )*deltax**(-2.0d0)
			dif(i-1,j-1) = dif(i-1,j-1) + ( u_aux(i,j-1)-2.0d0*u_aux(i,j)+u_aux(i,j+1) )*deltay**(-2.0d0)
		enddo
	enddo
	end subroutine laplac
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine rhs(P,deltax,deltay,rhs_x,rhs_y,dt,u_aux,v_aux,u_ant_aux,v_ant_aux)
	integer i,j
	double precision deltax, deltay, rhs_x(Nx-1,Ny-1), rhs_y(Nx-1,Ny-1), dt
	double precision u_aux(nx+1,ny+1), v_aux(nx+1,ny+1), u_ant_aux(nx+1,ny+1), v_ant_aux(nx+1,ny+1), p(nx-1,ny-1)
	double precision dif(Nx-1,Ny-1)
	double precision adv_x(Nx-1,Ny-1),adv_y(Nx-1,Ny-1)
	double precision gradp_x(Nx-1,Ny-1),gradp_y(Nx-1,Ny-1)
	rhs_x = 0.0d0       
	rhs_y = 0.0d0
	!Gradiente de Presión
	do j=1,(Ny-1)
		do i=1,(Nx-1)
			!En x
			if (i==1) then !primer borde en x
				gradp_x(i,j) = (P(i,j)-P(Nx-1,j))*deltax**(-1.0d0)
			else
				gradp_x(i,j) = (P(i,j)-P(i-1,j))*deltax**(-1.0d0)
			endif
			!En y
			if (j==1) then !primer borde en y
				gradp_y(i,j) = (P(i,j)-P(i,Ny-1))*deltay**(-1.0d0)
			else
				gradp_y(i,j) = (P(i,j)-P(i,j-1))*deltay**(-1.0d0)
			endif
		enddo
	enddo
	!Comenzamos a sumar cada componente del lado derecho:
	!1. Presión
	rhs_x = rhs_x - gradp_x      
	rhs_y = rhs_y - gradp_y
	!2. Advección
	call adveccion(u_aux,v_aux,deltax,deltay,adv_x,adv_y)
	rhs_x = rhs_x + 1.5d0*adv_x     
	rhs_y = rhs_y + 1.5d0*adv_y
	call adveccion(u_ant_aux,v_ant_aux,deltax,deltay,adv_x,adv_y)
	rhs_x = rhs_x - 0.5d0*adv_x     
	rhs_y = rhs_y - 0.5d0*adv_y	
	!3. Difusión
	call laplac(u_aux,deltax,deltay,dif)
	rhs_x = rhs_x + Re**(-1.0d0)*dif
	call laplac(v_aux,deltax,deltay,dif)
	rhs_y = rhs_y + Re**(-1.0d0)*dif
	!FINALMENTE
	rhs_x = dt*rhs_x
	rhs_y = dt*rhs_y
   end subroutine rhs  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_rhs
