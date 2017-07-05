module m_graph
use m_cons
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine graficar(u,v,p,w,cont,t_iter)

		implicit none
		double precision u(nx+1,ny), v(nx,ny+1), p(nx,ny), w(nx,ny), t_iter

		integer i, j, cont
		integer vv, pp, ww
		double precision aux_i, aux_j
		double precision aux_mag
		double precision u_aux(nx,ny), v_aux(nx,ny)
		character(len = 100) aux1, aux_nx, aux_ny, aux_re, aux_t, aux_cfl
		character(len = 100) grafname

		!Formato de auxiliares
		if (Nx<10) then
			write(aux_nx,'(I1)') Nx
		elseif (Nx<100) then
			write(aux_nx,'(I2)') Nx
		else
			write(aux_nx,'(I3)') Nx
		endif
		if (Ny<10) then
			write(aux_ny,'(I1)') Ny
		elseif (Ny<100) then
			write(aux_ny,'(I2)') Ny
		else
			write(aux_ny,'(I3)') Ny
		endif

		write(aux1,'(I6.6)') cont
		write(aux_t,'(F8.4)') t_iter
		write(aux_re, '(F5.0)') re
		write(aux_cfl, '(F4.2)') cfl

		!Escritura de tablas
		open(unit=10,file='plot/v.txt')			
		do j=1,Ny
			aux_j = dy*0.5d0 + dy*(j-1)
			do i = 1,Nx
				aux_i = dx*0.5d0 + dx*(i-1)
				u_aux(i,j) = 0.5d0*(u(i,j)+u(i+1,j))
				v_aux(i,j) = 0.5d0*(v(i,j)+v(i,j+1))
				aux_mag = sqrt(u_aux(i,j)*u_aux(i,j)+v_aux(i,j)*v_aux(i,j))
				if (aux_mag < 0.001d0) then
					aux_mag = 0.001d0
				endif
				write(10,*) aux_i, aux_j, 0.7d0*lx*u_aux(i,j)/nx, 0.7d0*lx*v_aux(i,j)/nx, aux_mag
			enddo
		end do
		close(10)

		open(unit=30,file='plot/p.txt')	
		do j=1,ny
			write(30,*) (p(i,j),i=1,nx)
		end do
		close(30)
		open(unit=30,file='plot/w.txt')											
		do j=1,ny
			write(30,*) (w(i,j),i=1,nx)
		end do
		close(30)

		grafname = 'plot/v/graf'//trim(aux1)//'.png'													
		open(unit=40, file="style.gnu")
		write(40,*) "set terminal pngcairo size 1200, 600 font 'Verdana, 12'"
		write(40,*) "set title 'Velocidad* | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
		&" | CFL = ",trim(aux_cfl)," | t* = ",trim(aux_t),"'"
		write(40,*) "set output '",trim(grafname),"'"
		write(40,*) "set xrange [0:",Lx,"]"
		write(40,*) "set yrange [0:",Ly,"]"
		write(40,*) "unset key"
		if (sq == 1) then
			write(40,*) "set object rect from 2.5,",ly/2.0d0 - 0.5d0," to 3.5,",ly/2.0d0 + 0.5d0," behind "
		endif
		write(40,*) 'set grid'
		write(40,*) 'set isosamples 2, 2'
		!write(40,*) 'set cbrange [0:1.8]'
		write(40,*) "set palette maxcolors 2000"
		write(40,*) "set palette defined ( 0 'blue',500 'cyan',1000 'green',1500 'yellow', 2000 'red')"
		write(40,*) "plot 'plot/v.txt' using 1:2:3:4:5 with vectors filled linecolor palette z"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu')

		grafname = 'plot/p/graf'//trim(aux1)//'.png'
		open(unit=40, file="style.gnu")
		write(40,*) "set terminal pngcairo size 800, 400 font 'Verdana,12"
		write(40,*) "set palette maxcolors 2000"
		write(40,*) "set autoscale fix"
		write(40,*) "set palette defined ( 0 'blue',500 'cyan',1000 'green',1500 'yellow', 2000 'red')"
		write(40,*) "unset key"
		if (sq == 1) then
			write(40,*) "set object rect from 2.5,",ly/2.0d0 - 0.5d0," to 3.5,",ly/2.0d0 + 0.5d0," front "
		endif
		write(40,*) "set title 'Presión* | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
		&" | CFL = ",trim(aux_cfl)," | t* = ",trim(aux_t),"'"
		!write(40,*) "set cbrange [-1.0:1.5]"
		write(40,*) "set output '",trim(grafname),"'"
		write(40,*) "plot 'plot/p.txt' matrix using (($1+0.5)*",lx/nx,"):($2+0.5)*",ly/ny,":3 with image"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu')

		grafname = 'plot/w/graf'//trim(aux1)//'.png'
		open(unit=40, file="style.gnu")	
		!write(40,*) "set terminal pngcairo size 800, 400 font 'Verdana,12'"
		write(40,*) "set terminal pngcairo size 400, 200 font 'Verdana,12'" !Para tumblr
		write(40,*) "set palette maxcolors 2000"
		write(40,*) "set autoscale fix"
		write(40,*) "set palette defined (0 'blue', 1000 'white', 2000 'red')"
		write(40,*) "unset key"
		if (sq == 1) then
			write(40,*) "set object rect from 2.5,",ly/2.0d0 - 0.5d0," to 3.5,",ly/2.0d0 + 0.5d0," front "
		endif
		!write(40,*) "set title 'Vorticidad* | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
		!&" | CFL = ",trim(aux_cfl)," | t* = ",trim(aux_t),"'"
		!write(40,*) "set cbrange [-15:15]"
		write(40,*) "unset colorbox"
		write(40,*) "set output '",trim(grafname),"'"
		write(40,*) "plot 'plot/w.txt' matrix using (($1+0.5)*",lx/nx,"):($2+0.5)*",ly/ny,":3 with image"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu')

		call system('rm plot/v.txt plot/p.txt plot/w.txt')

	end subroutine graficar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine graficar_cfl(dt_cfl,n,tiempo)
	
	implicit none
	integer j, n
	double precision dt_cfl(20000), tiempo(20000)
	character(len = 100) aux_nx, aux_cfl, aux_ny, aux_re, aux_sq, aux_h
		!Formato de auxiliares
		if (Nx<10) then
			write(aux_nx,'(I1)') Nx
		elseif (Nx<100) then
			write(aux_nx,'(I2)') Nx
		else
			write(aux_nx,'(I3)') Nx
		endif
		if (Ny<10) then
			write(aux_ny,'(I1)') Ny
		elseif (Ny<100) then
			write(aux_ny,'(I2)') Ny
		else
			write(aux_ny,'(I3)') Ny
		endif

		write(aux_re,'(F8.0)') Re
		write(aux_cfl,'(F4.2)') cfl
		write(aux_sq,'(I1.1)') sq
		write(aux_h, '(f3.1)') ly

		!Creación de tabla para guardar los dt utilizados
		open(unit=10,file='plot/cfl.txt')
			do j=1,n
				write(10,*) tiempo(j), dt_cfl(j)
			end do
		close(10)

		open(unit=40, file="style.gnu")
		write(40,*) "set terminal pngcairo dashed font 'Verdana,12'"
		write(40,*) "set output 'plot/cfl.png'"
		write(40,*) "set grid"
		write(40,*) "unset key"
		write(40,*) "set title 'dt vs t* | sq = ",trim(aux_sq)," | Ly = ",trim(aux_h)," | Re = "&
		&,trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny)," | CFL = ",trim(aux_cfl),"'"
		write(40,*) "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 # --- blue"
		write(40,*) "plot 'plot/cfl.txt' w lines ls 1"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu plot/cfl.txt')

	end subroutine graficar_cfl
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine armar_gif
		
		implicit none
		write(*,*)
		write(*,'(a15)') 'Armando gifs...'
		call system('convert -delay 10 -loop 0 plot/p/*.png plot/p.gif')
		call system('convert -delay 10 -loop 0 plot/v/*.png plot/v.gif')
		call system('convert -delay 10 -loop 0 plot/w/*.png plot/w.gif')

	end subroutine armar_gif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine teo_profile(u)

		implicit none
		double precision u(nx+1,ny)
		double precision aux_j
		character(len = 100) aux_re, aux_nx, aux_ny
		integer j
		open(unit=10,file='plot/teo.txt')		
		do j=1,ny
			aux_j = dy*0.5d0 + dy*(j-1)
			write(10,*) aux_j, u(nx+1,j) ,6.0d0*aux_j/ly*( 1.0d0 - aux_j/ly )
		end do
		close(10)

		write(aux_re, '(F5.0)') re

		if (Nx<10) then
			write(aux_nx,'(I1)') Nx
		elseif (Nx<100) then
			write(aux_nx,'(I2)') Nx
		else
			write(aux_nx,'(I3)') Nx
		endif
		if (Ny<10) then
			write(aux_ny,'(I1)') Ny
		elseif (Ny<100) then
			write(aux_ny,'(I2)') Ny
		else
			write(aux_ny,'(I3)') Ny
		endif

		open(unit=40, file="style.gnu")	
		write(40,*) "set terminal pngcairo size 800, 600 font 'Verdana,12'"
		write(40,*) "set autoscale fix"
		write(40,*) "set title 'Validación Flujo Poiseuille | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),"'"
		write(40,*) "set output 'plot/valid.png'"
		write(40,*) "set grid"
		write(40,*) "plot 'plot/teo.txt' using 2:1 title 'simulación' pt 7 ,\"
		write(40,*) "'plot/teo.txt' using 3:1 title 'teórica' with lines"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu plot/teo.txt')

	end subroutine teo_profile
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_graph
