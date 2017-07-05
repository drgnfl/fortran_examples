!Tarea 2, CFD, UTFSM 2016-2
!Pablo Cárdenas Zamorano
program main
!Preambulo: Módulos a utilizar
use m_cons
use m_coef
use m_tdma
use m_solve_uv
use m_solve_dp
use m_graph
!Declaración de variables
implicit none
include 'glob_var.f90'
!Se inicia el programa
write(*,'(a25)') 'Ejecutando el programa...'
write(*,*)
!Inicializando system_clock
CALL system_clock(count_rate=cr)
CALL system_clock(count_max=cm)
rate = real(cr)
call system_clock(c1)
!Creación de directorios
call system('rm -r plot')
call system('mkdir plot')
call system('mkdir plot/v')
call system('mkdir plot/p')
call system('mkdir plot/w')
!Valores Iniciales
u		=	0.0d0
u(1,:)	=	1.0d0
v		=	0.0d0
p		=	0.0d0
!Inicio de Algoritmo
pp 		= 	0.0d0
dp		=	0.0d0
error_pp = 	1.0d0
t_iter = 0.0d0
n = 0
!Comienzo de las iteraciones temporales
do k = 1, kmax
	!Barra de carga
	call progress(floor(k*(kmax**(-1.0d0))*100.0d0)-1)
	!Cálculo del paso temporal
	n = n + 1																
	dt_cfl(n) = cfl*(maxval(abs(u))*dx**(-1.0d0)+maxval(abs(v))*dy**(-1.0d0))**(-1.0d0)
	tiempo(n) = t_iter											
	dt = dt_cfl(n)						
	!Calcular coeficientes para la iteración
	call u_coef_calc( u,v, u_ae,u_aw,u_an,u_as,u_anb,beta_x,dt )
	call v_coef_calc( u,v, v_ae,v_aw,v_an,v_as,v_anb,beta_y,dt )
	call p_coef_calc( beta_x, beta_y, beta_p )
	!Implementación de SIMPLE
	qq = 0			!Contador para revisar la divergencia del programa
	do while (error_pp > tol_div)
		qq = qq + 1
		if (qq == 50) then
			write(*,*)
			write(*,*)
			write(*,'(a35)') 'El programa divergió, disminuir cfl!'
			write(*,*)
			goto 10
		endif
		!1. Suponer p*
		pp = pp + alpha_p*dp
		!2. Calcular u* y v*
		call solve_uu(u,v,pp,uu, u_ae,u_aw,u_an,u_as,u_anb,dt)
		call solve_vv(u,v,pp,vv, v_ae,v_aw,v_an,v_as,v_anb,dt)
		!3. Cálculo de la corrección de presión (dp)
		call mass_cons(uu,vv,div)
		call solve_dp(u,v,uu,vv,dp, div, beta_x, beta_y, beta_p)
		!4. Cálculo del error basado en la divergencia
		error_pp = maxval( abs(div) )
		!5. Cálculo de u^{n+1} y v^{n+1}
		call solve_uv_new(uu,vv,dp,beta_x,beta_y,u_new,v_new)
	enddo
	t_iter = t_iter + dt													!Actualizamos al tiempo actual
	!Creación de gráficos (100 frames)	
	!if (mod(k,kmax/100) == 0 .or. k == 1) then
		!Cálculo de vorticidad
		call vort(u_new,v_new,w)
		call graficar(u_new,v_new,pp+alpha_p*dp,w,k,t_iter)
	!endif
	!Variables para iniciar la próxima iteración
	error_pp = 1.0d0
	u = u_new
	v = v_new
enddo
write(*,*)
!Creación de animaciones
10 if (gif==1) then
	call armar_gif
endif
call graficar_cfl(dt_cfl,n,tiempo)
!Validación con perfil teorico
if (tp == 1) then
	call teo_profile(u)
endif
!Fin del programa
call system_clock(c2)
call system ('bash red.sh')
write(*,*)
write(*,'(a36)') 'El programa finalizó correctamente!'
write(*,'(a21,f8.2,a4)')'Tiempo transcurrido: ', (c2 - c1)/rate/60.0d0, ' [min]'
write(*,*)
call system ('bash norm.sh')
call system('aplay -q end.wav')	
end program
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
