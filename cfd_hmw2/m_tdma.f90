module m_tdma
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine tdma(d1,dp,d2,rhs,n,sol)
		implicit none
		integer k, n
		double precision d1(n), dp(n), d2(n), rhs(n)
		double precision m, sol(n)
		!Forward
		do k = 2, n
			if (dp(k-1) == 0.0d0) then
				dp(k-1) = 1.0d-10
			endif
			m		=	d1(k) / dp(k-1)
			dp(k) 	= 	dp(k) - m*d2(k-1)
			rhs(k) 	= 	rhs(k) - m*rhs(k-1)
		enddo
		!Backward
		sol(n) = rhs(n) / dp(n)
		do k = n-1, 1, -1
			sol(k)	=	(rhs(k) - d2(k)*sol(k+1)) / dp(k)
		enddo
	end subroutine tdma
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_tdma
