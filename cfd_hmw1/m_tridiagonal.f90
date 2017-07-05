!Módulo para resolver un sistema lineal tridiagonal simétrica con condiciones de borde periodicas. Alpha es la diagonal principal y Beta es la diagonal secundaria.
module m_tridiagonal
use m_constantes
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine tridiag(alpha,beta,rhs,N,sol)
		integer i,N
		double precision alpha,beta,rhs(N),sol(N)
		double precision u(N),v(N),yy(N),y(N),qq(N),q(N),vv(N)
		u = 0.0d0
		u(1) = -alpha
		u(N) = beta
		vv = 0.0d0
		vv(1) = 1.0d0
		vv(N) = -beta*alpha**(-1.0d0)
	!RESOLVIENDO LA TRIDIAGONAL MODIFICADA, 1: LOWER
		do i=1,N
			if (i==1) then
				yy(i) = rhs(i)
				qq(i) = u(i)
				v(i) = 2.0d0*alpha
			else
				yy(i) = rhs(i) - beta*yy(i-1)*v(i-1)**(-1.0d0)
				qq(i) = u(i) - beta*qq(i-1)*v(i-1)**(-1.0d0)
				v(i) = alpha - beta**2.0d0*v(i-1)**(-1.0d0)
			endif
		enddo
		v(N) = (alpha + beta**2.0d0*alpha**(-1.0d0)) - beta**2.0d0*v(N-1)**(-1.0d0)
	!RESOLVIENDO LA TRIDIAGONAL MODIFICADA, 2: UPPER
		do i=N,1,-1
			if (i==N) then
				y(i) = yy(i)*v(i)**(-1.0d0)
				q(i) = qq(i)*v(i)**(-1.0d0)
			else
				y(i) = (yy(i)-beta*y(i+1))*v(i)**(-1.0d0)
				q(i) = (qq(i)-beta*q(i+1))*v(i)**(-1.0d0)
			endif
		enddo
		sol = y - (dot_product(vv,y)*(1.0d0 + dot_product(vv,q))**(-1.0d0))*q
	end subroutine tridiag
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine tridiag_complx(alpha,beta,rhs,N,sol)
		integer i,N
		double precision alpha,beta
		complex*16 rhs(N),sol(N)
		complex*16 u(N),v(N),yy(N),y(N),qq(N),q(N),vv(N)
		u = 0.0d0
		u(1) = -alpha
		u(N) = beta
		vv = 0.0d0
		vv(1) = 1.0d0
		vv(N) = -beta*alpha**(-1.0d0)
	!RESOLVIENDO LA TRIDIAGONAL MODIFICADA, 1: LOWER
		do i=1,N
			if (i==1) then
				yy(i) = rhs(i)
				qq(i) = u(i)
				v(i) = 2.0d0*alpha
			else
				yy(i) = rhs(i) - beta*yy(i-1)*v(i-1)**(-1.0d0)
				qq(i) = u(i) - beta*qq(i-1)*v(i-1)**(-1.0d0)
				v(i) = alpha - beta**2.0d0*v(i-1)**(-1.0d0)
			endif
		enddo
		v(N) = (alpha + beta**2.0d0*alpha**(-1.0d0)) - beta**2.0d0*v(N-1)**(-1.0d0)
	!RESOLVIENDO LA TRIDIAGONAL MODIFICADA, 2: UPPER
		do i=N,1,-1
			if (i==N) then
				y(i) = yy(i)*v(i)**(-1.0d0)
				q(i) = qq(i)*v(i)**(-1.0d0)
			else
				y(i) = (yy(i)-beta*y(i+1))*v(i)**(-1.0d0)
				q(i) = (qq(i)-beta*q(i+1))*v(i)**(-1.0d0)
			endif
		enddo
		sol = y - (dot_product(vv,y)*(1.0d0 + dot_product(vv,q))**(-1.0d0))*q
	end subroutine tridiag_complx
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_tridiagonal
