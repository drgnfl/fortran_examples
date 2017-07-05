module m_coef
use m_cons
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine u_coef_calc( u,v, u_ae,u_aw,u_an,u_as,u_anb,beta_x,dt )
		implicit none
		double precision u(nx+1,ny), v(nx,ny+1), dt
		double precision u_ae(nx,ny), u_aw(nx,ny), u_as(nx,ny), u_an(nx,ny), u_anb(nx,ny), beta_x(nx+1,ny)

		integer i, j
		double precision Fe, Fw, Fs, Fn
		double precision u_e, u_w, v_s, v_n
		
		beta_x = 0.0d0
		do j = 1 , ny
			do i = 2, nx
				u_e = 0.5d0*( u(i+1,j) + u(i,j) )
				u_w = 0.5d0*( u(i-1,j) + u(i,j) )
				v_s = 0.5d0*( v(i-1,j) + v(i,j) )
				v_n = 0.5d0*( v(i-1,j+1) + v(i,j+1) )

				Fe = u_e*dy
				Fw = u_w*dy
				Fn = v_n*dx
				Fs = v_s*dx

				u_ae(i,j) =  gam*dy/dx + max(-Fe,0.0d0)
				u_aw(i,j) =  gam*dy/dx + max( Fw,0.0d0)
				u_an(i,j) =  gam*dx/dy + max(-Fn,0.0d0)
				u_as(i,j) =  gam*dx/dy + max( Fs,0.0d0)
				u_anb(i,j) = u_ae(i,j) + u_aw(i,j) + Fe - Fw + u_an(i,j) + u_as(i,j) + Fn - Fs

				beta_x(i,j) = (dy**2.0d0)/( dx*dy/dt + u_anb(i,j) )
			enddo
		enddo
		beta_x(nx+1,:) = beta_x(nx,:)
		!Condición de Obstáculo
		if (sq == 1) then
			do j = ny/2+1-my, ny/2+my
				do i = mx*5+1, mx*7+1
					beta_x(i,j) =  0.0d0
				enddo
			enddo
		endif
	end subroutine u_coef_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine v_coef_calc( u,v, v_ae,v_aw,v_an,v_as,v_anb,beta_y,dt )
		implicit none
		double precision u(nx+1,ny), v(nx,ny+1), dt
		double precision v_ae(nx,ny), v_aw(nx,ny), v_as(nx,ny), v_an(nx,ny), v_anb(nx,ny),beta_y(nx,ny+1)

		integer i, j
		double precision Fe, Fw, Fs, Fn
		double precision u_e, u_w, v_s, v_n

		beta_y = 0.0d0
		do j = 2 , ny
			do i = 1, nx
				u_e = 0.5d0*( u(i+1,j) + u(i+1,j-1) )
				u_w = 0.5d0*( u(i,j) + u(i,j-1) )
				v_s = 0.5d0*( v(i,j) + v(i,j-1) )
				v_n = 0.5d0*( v(i,j) + v(i,j+1) )

				Fe = u_e*dy
				Fw = u_w*dy
				Fn = v_n*dx
				Fs = v_s*dx

				v_ae(i,j) =  gam*dy/dx + max(-Fe,0.0d0)
				v_aw(i,j) =  gam*dy/dx + max( Fw,0.0d0)
				v_an(i,j) =  gam*dx/dy + max(-Fn,0.0d0)
				v_as(i,j) =  gam*dx/dy + max( Fs,0.0d0)
				v_anb(i,j) = v_ae(i,j) + v_aw(i,j) + Fe - Fw + v_an(i,j) + v_as(i,j) + Fn - Fs

				beta_y(i,j) = (dx**2.0d0)/( dx*dy/dt + v_anb(i,j) )
			enddo
		enddo
		!Condición de Obstáculo
		if (sq == 1) then
			do j = ny/2+1-my, ny/2+my+1
				do i = mx*5+1, mx*7
					beta_y(i,j) =  0.0d0
				enddo
			enddo
		endif
	end subroutine v_coef_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine p_coef_calc( beta_x,beta_y,beta_p )
		implicit none
		double precision beta_x(nx+1,ny),beta_y(nx,ny+1),beta_p(nx,ny)
		integer i, j

		do j = 1, ny
			do i = 1, nx
				beta_p(i,j) = beta_x(i,j) + beta_x(i+1,j) + beta_y(i,j) + beta_y(i,j+1)
			enddo
		enddo
	end subroutine p_coef_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_coef
