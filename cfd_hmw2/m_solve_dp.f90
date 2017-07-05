module m_solve_dp
use m_cons
use m_tdma
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine solve_dp(u,v,uu,vv,dp, div, beta_x,beta_y,beta_p )
		implicit none
		double precision u(nx+1, ny), v(nx,ny+1)
		double precision uu(nx+1, ny), vv(nx,ny+1)
		double precision dp (nx,ny)

		integer i, j
		double precision beta_x(nx+1,ny), beta_y(nx,ny+1), beta_p(nx,ny)
		double precision rhs(nx,ny), dpp(nx,ny)
		double precision div(nx,ny)
		double precision d1(nx), diag_p(nx), d2(nx)
		double precision error_p

		dpp 		= 0.0d0
		dp			= 0.0d0
		error_p 	= 1.0d0

		do while (error_p > tol_p)
			!Lado derecho de la ecuación de la presión
			rhs = -div
			do j = 1 , ny
				do i = 1, nx
					if (j==1) then
						rhs(i,j) = rhs(i,j) + beta_y(i,j+1)*dpp(i,j+1)
					elseif (j==ny) then
						rhs(i,j) = rhs(i,j) + beta_y(i,j)*dpp(i,j-1)
					else
						rhs(i,j) = rhs(i,j) + beta_y(i,j)*dpp(i,j-1) + beta_y(i,j+1)*dpp(i,j+1)
					endif
				enddo
			enddo
			!Resolución del sistema pentadiagonal por iteración
			do j=1,ny
				!Diagonal 1
				d1 = -beta_x(1:nx,j)
				!Diagonal Principal
				diag_p = beta_p(:,j)
				!Condicion de salida se fija como referencia
				!Condición de obstáculo
				if (sq == 1) then
					if (j >= ny/2+1-my .and. j <= ny/2+my) then
						do i = mx*5+1, mx*7
							diag_p(i) =  diag_p(i) + 10.0d0**60.0d0
						enddo
					endif
				endif
				!Diagonal 2
				d2 = -beta_x(2:nx+1,j)
				call tdma( d1,diag_p,d2,rhs(:,j),nx,dp(:,j) )
			enddo

			error_p = maxval( abs( dpp-dp) )
			dpp = dp
		enddo

	end subroutine solve_dp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine solve_uv_new(uu,vv,dp,beta_x,beta_y,u_new,v_new)
		double precision uu(nx+1,ny), vv(nx,ny+1), dp(nx,ny)
		double precision beta_x(nx+1,ny), beta_y(nx,ny+1)
		double precision u_new(nx+1,ny), v_new(nx,ny+1)

		integer i,j

		!Calculo de velocidades nuevas
		do j = 1 , ny
			do i = 2, nx
				u_new(i,j) = uu(i,j) + beta_x(i,j)/dy*(dp(i-1,j)-dp(i,j))
			enddo
		enddo
		u_new(1,:) = 1.0d0
		u_new(nx+1,:) = u_new(nx,:)

		do j = 2 , ny
			do i = 1, nx
				v_new(i,j) = vv(i,j) + beta_y(i,j)/dx*(dp(i,j-1)-dp(i,j))
			enddo
		enddo
		v_new(:,1) 		= 0.0d0
		v_new(:,ny+1) 	= 0.0d0

	end subroutine solve_uv_new
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine mass_cons(u,v,div)
		implicit none
		double precision u(nx+1,ny), v(nx,ny+1), div(nx,ny)
		integer i, j
		div = 0.0d0
		do j = 1, ny
			do i = 1, nx
				div(i,j) = ( u(i+1,j) - u(i,j) )*dy + ( v(i,j+1) - v(i,j) )*dx
			enddo
		enddo
	end subroutine mass_cons
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine vort(u,v,w)
		implicit none
		integer i,j
		double precision u(Nx+1,Ny),v(Nx,Ny+1),w(Nx,Ny)
		w = 0.0d0
		do j = 1,Ny
			do i = 1,Nx
					if (i==1) then
						w(i,j) = 0.333d0*( v(i+1,j)+v(i+1,j+1) - 0.0d0 - 0.0d0 ) / dx 
					elseif (i==nx) then
						w(i,j) = 0.25d0*( v(i,j)+v(i,j+1) - v(i-1,j) - v(i-1,j+1) ) / dx
					else
						w(i,j) = 0.25d0*( v(i+1,j)+v(i+1,j+1) - v(i-1,j) - v(i-1,j+1)) / dx
					endif

					if (j==1) then
						w(i,j) = w(i,j) - 0.333d0*( u(i,j+1)+u(i+1,j+1)-0.0d0-0.0d0 ) / dy 
					elseif (j==ny) then
						w(i,j) = w(i,j) - 0.333d0*( 0.0d0+0.0d0-u(i,j-1)-u(i+1,j-1) ) / dy 
					else
						w(i,j) = w(i,j) - 0.25*( u(i,j+1)+u(i+1,j+1)-u(i,j-1)-u(i+1,j-1) ) / dy 
					endif
			enddo
		enddo
	end subroutine vort
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_solve_dp
