module m_solve_uv
use m_cons
use m_tdma
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine solve_uu(u,v,pp,uu, u_ae,u_aw,u_an,u_as,u_anb,dt)
		implicit none
		double precision u(nx+1,ny), v(nx,ny+1), pp(nx,ny)
		double precision u_ae(nx,ny), u_aw(nx,ny), u_as(nx,ny), u_an(nx,ny), u_anb(nx,ny)
		double precision uu(nx+1,ny), dt

		integer i, j
		double precision rhs(nx,ny), uuu(nx+1,ny)
		double precision d1(nx-1), dp(nx-1), d2(nx-1)
		double precision error
		uuu 	= u
		uu		= 0.0d0
		error 	= 1.0d0

		do while (error > tol_u)
			!Lado derecho de la ecuación de momentum
			do j = 1 , ny
				do i = 2, nx
					!Condición de pared en los bordes
					if (j==1) then
						rhs(i,j) = ( pp(i-1,j) - pp(i,j) )*dy + u(i,j)*dx*dy/dt - u_as(i,j)*uuu(i,j) + u_an(i,j)*uuu(i,j+1)
					elseif (j==ny) then
						rhs(i,j) = ( pp(i-1,j) - pp(i,j) )*dy + u(i,j)*dx*dy/dt + u_as(i,j)*uuu(i,j-1) - u_an(i,j)*uuu(i,j)
					!Puntos interiores
					else
						rhs(i,j) = ( pp(i-1,j) - pp(i,j) )*dy + u(i,j)*dx*dy/dt + u_as(i,j)*uuu(i,j-1) + u_an(i,j)*uuu(i,j+1)
					endif
				enddo
			enddo
							
			!Resolución del sistema pentadiagonal por iteración
			do j=1,ny
				!Condicion de borde a la entrada
				rhs(2,j) = rhs(2,j) + u(1,j)*u_aw(2,j)
				!Diagonal 1
				d1 = -u_aw(2:nx,j)
				!Diagonal Principal
				dp = dx*dy/dt + u_anb(2:nx,j)
				!Condicion de salida
				dp(nx-1) = dp(nx-1) - u_ae(nx,j)
				!Diagonal 2
				d2 = -u_ae(2:nx,j)
				!Condición de obstáculo
				if (sq == 1) then 
					if (j >= ny/2+1-my .and. j <= ny/2+my) then
						do i = mx*5+1, mx*7+1
							dp(i-1) =  dp(i-1) + 10.0d0**80.0d0
						enddo
					endif
				endif

				call tdma( d1,dp,d2,rhs(2:nx,j),nx-1,uu(2:nx,j) )
			enddo
			!Condiciones de Borde
			uu(1,:) = 1.0d0
			uu(nx+1,:) = uu(nx,:)

			error = maxval( abs( uuu-uu ) )
			uuu = uu
		enddo

	end subroutine solve_uu
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine solve_vv(u,v,pp,vv, v_ae,v_aw,v_an,v_as,v_anb,dt)
		implicit none
		double precision u(nx+1,ny), v(nx,ny+1), pp(nx,ny)
		double precision v_ae(nx,ny), v_aw(nx,ny), v_as(nx,ny), v_an(nx,ny), v_anb(nx,ny)
		double precision vv(nx,ny+1), dt

		integer i, j
		double precision rhs(nx,ny), vvv(nx,ny+1)
		double precision d1(nx), dp(nx), d2(nx)
		double precision error
		vvv 	= v
		vv		= 0.0d0
		error 	= 1.0d0

		do while (error > tol_v)
			!Lado derecho de la ecuación de momentum
			do j = 2 , ny
				do i = 1, nx
					rhs(i,j) = ( pp(i,j-1) - pp(i,j) )*dx + v(i,j)*dx*dy/dt + v_as(i,j)*vvv(i,j-1) + v_an(i,j)*vvv(i,j+1)
				enddo
			enddo
			!Resolución del sistema pentadiagonal por iteración
			do j=2,ny
				!Diagonal 1
				d1 = -v_aw(1:nx,j)
				!Diagonal Principal
				dp = dx*dy/dt + v_anb(1:nx,j)
				!Condicion de borde a la entrada
				dp(1) = dp(1) + v_aw(1,j)
				!Condicion de salida
				dp(nx) = dp(nx) - v_ae(nx,j)
				!Diagonal 2
				d2 = -v_ae(1:nx,j)

				!Condición de obstáculo
				if (sq == 1) then
					if (j >= ny/2+1-my .and. j <= ny/2+1+my) then
						do i = mx*5+1, mx*7
							dp(i) =  dp(i) + 10.0d0**80.0d0
						enddo
					endif
				endif
				call tdma( d1,dp,d2,rhs(1:nx,j),nx,vv(1:nx,j) )
			enddo
			!Condiciones de Borde
			vv(:,1) 	= 0.0d0
			vv(:,ny+1) 	= 0.0d0

			error = maxval( abs( vvv-vv ) )
			vvv = vv
		enddo

	end subroutine solve_vv
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_solve_uv
