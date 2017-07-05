!Módulo para la creación de gráficos y escritura de tablas
module m_graficos
use m_constantes
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine graficar(u_aux,v_aux,P,vort,cont,t_iter,deltax,deltay,vv,pp,ww)
		integer i, j, cont
		integer vv, pp, ww
		double precision u_aux(Nx+1,Ny+1), v_aux(Nx+1,Ny+1), P(Nx-1,Ny-1), vort(Nx-1,Ny-1), aux_i, aux_j, t_iter
		double precision deltax, deltay,aux_mag
		double precision u(nx-1,ny-1), v(nx-1,ny-1)
		character(len = 100) aux1, aux_nx, aux_cfl, aux_ny, aux_re, aux_t, aux_pert
		character(len = 100) grafname

		!Formato de auxiliares
		if (Nx<10) then
			write(aux_nx,'(I1)') Nx
			write(aux_ny,'(I1)') Ny
		elseif (Nx<100) then
			write(aux_nx,'(I2)') Nx
			write(aux_ny,'(I2)') Ny
		else
			write(aux_nx,'(I3)') Nx
			write(aux_ny,'(I3)') Ny
		endif
		write(aux1,'(I6.6)') cont
		write(aux_re,'(F8.0)') Re
		write(aux_t,'(F8.4)') t_iter
		write(aux_cfl,'(F4.2)') cfl
		write(aux_pert,'(I1.1)') perturbacion

		!Escritura de tablas
		open(unit=10,file='plot/v.txt')																	!Velocidad
		do j=2,Ny
			aux_j = deltay*0.5d0 + deltay*(j-2)
			do i = 2,Nx
				aux_i = deltax*0.5d0 + deltax*(i-2)
				u(i-1,j-1) = 0.5d0*(u_aux(i,j)+u_aux(i+1,j))
				v(i-1,j-1) = 0.5d0*(v_aux(i,j)+v_aux(i,j+1))
				aux_mag = sqrt(u(i-1,j-1)*u(i-1,j-1)+v(i-1,j-1)*v(i-1,j-1))
				write(10,*) aux_i, aux_j, 0.85d0*u(i-1,j-1)/aux_mag/Nx, 0.85d0*v(i-1,j-1)/aux_mag/(Ny*lambda), aux_mag
			enddo
		end do
		close(10)
		open(unit=30,file='plot/p.txt')																	!Presión
		do j=1,Ny-1
			write(30,*) (p(i,j),i=1,Nx-1)
		end do
		close(30)
		open(unit=50,file='plot/vort.txt')																!Vorticidad
		do j=1,Ny-1
			write(50,*) (vort(i,j),i=1,Nx-1)
		end do
		close(50)

		!Creación de Gráficos
		if (vv == 1) then
			grafname = 'plot/v/graf'//trim(aux1)//'.png'													!Velocidad
			open(unit=40, file="style.gnu")		!Formato del gráfico
			write(40,*) "set terminal pngcairo size 1920, 1080 font 'Verdana, 10'"
			write(40,*) "set title 'Velocidad* | p = ",trim(aux_pert)," | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
			&" | CFL = ",trim(aux_cfl)," | t* = ",trim(aux_t),"'"
			write(40,*) "set output '",trim(grafname),"'"
			write(40,*) "set xrange [0:1]"
			write(40,*) "set yrange [0:",lambda**(-1.0d0),"]"
			write(40,*) "unset key"
			write(40,*) 'set grid'
			write(40,*) 'set isosamples 2, 2'
			write(40,*) 'set cbrange [0:1.5]'
			write(40,*) "set palette maxcolors 2000"
			write(40,*) "set palette defined ( 0 'white', 250 'red', 2000 'black')"
			write(40,*) "plot 'plot/v.txt' using 1:2:3:4:5 with vectors filled linecolor palette z"
		 	close(unit=40)
			call system('gnuplot style.gnu')
			call system('rm style.gnu')
		endif

		if (pp == 1) then		
			grafname = 'plot/p/graf'//trim(aux1)//'.png'													!Presión
			open(unit=40, file="style.gnu")
			write(40,*) "set pm3d map interpolate 0,0"
			write(40,*) "set terminal pngcairo font 'Verdana,10"
			write(40,*) "set palette maxcolors 200"
			write(40,*) "set palette defined (0 'blue', 100 'white', 200 'red')"
			write(40,*) "unset key"
			write(40,*) "set xrange [0:1]"
			write(40,*) "set yrange [0:",lambda**(-1.0d0),"]"
			write(40,*) "set title 'Presión* | p = ",trim(aux_pert)," | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
			&" | CFL = ",trim(aux_cfl)," | t* = ",trim(aux_t),"'"
			write(40,*) "set cbrange [-0.3:0.3]"
			write(40,*) "set output '",trim(grafname),"'"
			write(40,*) "splot 'plot/p.txt' matrix using ($1/",(Nx-2),"):($2/",(Ny-2)*lambda,"):3"
		 	close(unit=40)
			call system('gnuplot style.gnu')
			call system('rm style.gnu')
		endif

		if (ww == 1) then
			grafname = 'plot/vort/graf'//trim(aux1)//'.png'													!Vorticidad
			open(unit=40, file="style.gnu")	
			write(40,*) "set pm3d map interpolate 0,0"
			write(40,*) "set terminal pngcairo dashed font 'Verdana,10'"
			write(40,*) "set palette maxcolors 2000"
			write(40,*) "set palette defined (0 'blue', 1000 'white', 2000 'red')"
			!write(40,*) "set palette defined (0 'yellow', 800 'yellow', 900 'blue', 1000 'white', 1100 'red', 1200 'green', 2000 'green')"
			write(40,*) "unset key"
			write(40,*) "set xrange [0:1]"
			write(40,*) "set yrange [0:",lambda**(-1.0d0),"]"
			write(40,*) "set title 'Vorticidad* | p = ",trim(aux_pert)," | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
			&" | CFL = ",trim(aux_cfl)," | t* = ",trim(aux_t),"'"
			write(40,*) "set cbrange [-22:22]"
			!write(40,*) "set cbrange [-250:250]"
			!write(40,*) "unset colorbox"
			write(40,*) "set output '",trim(grafname),"'"
			write(40,*) "splot 'plot/vort.txt' matrix using ($1/",Nx-2,"):($2/",(Ny-2)*lambda,"):3"
		 	close(unit=40)
			call system('gnuplot style.gnu')
			call system('rm style.gnu')
		endif

		call system('rm plot/v.txt plot/p.txt plot/vort.txt')
	end subroutine graficar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine graficar_cfl(dt_cfl,n,tiempo)
	integer j, n
	double precision dt_cfl(20000), tiempo(20000)
	character(len = 100) aux_nx, aux_cfl, aux_ny, aux_re, aux_pert
		!Formato de auxiliares
		if (Nx<10) then
			write(aux_nx,'(I1)') Nx
			write(aux_ny,'(I1)') Ny
		elseif (Nx<100) then
			write(aux_nx,'(I2)') Nx
			write(aux_ny,'(I2)') Ny
		else
			write(aux_nx,'(I3)') Nx
			write(aux_ny,'(I3)') Ny
		endif
		write(aux_re,'(F8.0)') Re
		write(aux_cfl,'(F4.2)') cfl
		write(aux_pert,'(I1.1)') perturbacion
		!Creación de tabla para guardar los dt utilizados
		open(unit=10,file='plot/cfl.txt')
			do j=1,n
				write(10,*) tiempo(j), dt_cfl(j)
			end do
		close(10)
		open(unit=40, file="style.gnu")
		write(40,*) "set terminal pngcairo dashed font 'Verdana,10'"
		write(40,*) "set output 'plot/cfl_",trim(aux_nx),".png'"
		write(40,*) "set grid"
		write(40,*) "unset key"
		write(40,*) "set title 'dt vs t* | p = ",trim(aux_pert)," | Re = ",trim(aux_re)," | ",trim(aux_nx)," x ",trim(aux_ny),&
		&" | CFL = ",trim(aux_cfl),"'"
		write(40,*) "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 # --- blue"
		write(40,*) "plot 'plot/cfl.txt' w lines ls 1"
	 	close(unit=40)
		call system('gnuplot style.gnu')
		call system('rm style.gnu plot/cfl.txt')
	end subroutine graficar_cfl
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine armar_gif
		write(*,*)
		write(*,*)
		write(*,'(a15)') 'Armando gifs...'
		write(*,*)
		if (gp == 1) then
			call system('convert -delay 10 -loop 0 plot/p/*.png plot/pres.gif')
		endif
		if (gv == 1) then
			call system('convert -delay 10 -loop 0 plot/v/*.png plot/vel.gif')
		endif
		if (gw == 1) then
			call system('convert -delay 10 -loop 0 plot/vort/*.png plot/vort.gif')
		endif
	end subroutine armar_gif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module m_graficos
