!Módulo para la última etapa del programa. Acá se realiza:
!Se calcula el gradiente de Phi.
!Se despeja la velocidad en el tiempo n+1 y se actualizan los datos de la velocidad en el tiempo n-1.
module m_desp_veloc
use m_constantes
use m_iniciar
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine despejar_velocidad(uu,vv,u,v,u_aux,v_aux,u_ant_aux,v_ant_aux,Phi,deltax,deltay,dt)
		integer i,j
		double precision uu(Nx-1,Ny-1),vv(Nx-1,Ny-1),u(Nx-1,Ny-1),v(Nx-1,Ny-1),u_ant_aux(Nx+1,Ny+1),v_ant_aux(Nx+1,Ny+1), Phi(Nx-1,Ny-1)
		double precision u_aux(Nx+1,Ny+1), v_aux(Nx+1,Ny+1)
		double precision deltax, deltay,dt
		double precision grad_x(Nx-1,Ny-1), grad_y(Nx-1,Ny-1)
		!Gradientes de Phi
		do j=1,(Ny-1)
			do i=1,(Nx-1)
				!En x
				if (i==1) then !primer borde en x
					grad_x(i,j) = (Phi(1,j)-Phi(Nx-1,j))*deltax**(-1.0d0)
				else
					grad_x(i,j) = (Phi(i,j)-Phi(i-1,j))*deltax**(-1.0d0)
				endif
				!En y
				if (j==1) then !primer borde en y
					grad_y(i,j) = (Phi(i,j)-Phi(i,Ny-1))*deltay**(-1.0d0)
				else
					grad_y(i,j) = (Phi(i,j)-Phi(i,j-1))*deltay**(-1.0d0)
				endif
			enddo
		enddo
		!Actualizando y calculando:
		u_ant_aux = u_aux
		u = uu - dt*grad_x

		v_ant_aux = v_aux
		v = vv - dt*grad_y
		
		call auxiliares(u,v,u_aux,v_aux)
	end subroutine despejar_velocidad
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_desp_veloc
