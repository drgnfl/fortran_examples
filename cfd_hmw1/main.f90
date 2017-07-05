!Tarea 1, CFD, UTFSM 2016-2
!Pablo Cárdenas Zamorano
program main
!Preambulo: Módulos a utilizar
use m_constantes
use m_iniciar
use m_rhs
use m_tridiagonal
use m_presion
use m_desp_veloc
use m_graficos

!Declaración de variables
implicit none
include 'variables.f90'
double precision start, finish
double precision tiempo(20000)
double precision aux_Phi(Nx+1,Ny+1), lap_phi(Nx-1,Ny-1)

start = 0.0d0
finish = 0.0d0
call cpu_time(start)		!Mide el tiempo

!Creación de directorios
call system('rm -r plot') 	!Borra si ya está creado
call system('mkdir plot')
call system('mkdir plot/v')
call system('mkdir plot/p')
call system('mkdir plot/vort')

!Inicialización de los valores máximos y contarodes
max_u 		= 0.0d0
max_v 		= 0.0d0
max_p 		= 0.0d0
max_vort 	= 0.0d0
dt 			= 0.0d0
t_iter		= 0.0d0 	!Iniciamos al tiempo cero
n 			= 0 		!Cuenta el número de iteraciones
cont 		= 0 		!Cuenta el número de gráficos

!Inicialización de los valores iniciales
call inicializacion(deltax,deltay,u,v,p)
call auxiliares(u,v,u_aux,v_aux)		!Se amplia u y v para que sea mas facil evaluar las condiciones periodicas

u_ant_aux = u_aux						!Se genera un euler explícito en el primer paso de tiempo
v_ant_aux = v_aux

call vorticidad(u_aux,v_aux,deltax,deltay,vort)
call graficar(u_aux,v_aux,P,vort,cont,t_iter,deltax,deltay,gp,gv,gw)

if (vp == 1) then
	call validar_poisson(deltax,deltay)
endif

!Comienzo de la iteración
write(*,'(a26)') 'Integrando en el tiempo...'
write(*,*)
do while (t_iter<t_mod)
	call progress(floor(cont*aux_n**(-1.0d0)*100.0d0)-1)												!Agrega una barra de carga
	
	n = n + 1																							!Contador de iteración
	dt_cfl(n) = cfl*(maxval(abs(u))*deltax**(-1.0d0)+maxval(abs(v))*deltay**(-1.0d0))**(-1.0d0)			!Se guarda el dt a utilizar
	tiempo(n) = t_iter																					!Se almacena el tiempo
	dt = dt_cfl(n)																						!Se asigna el dt

	call rhs(P,deltax,deltay,rhs_x,rhs_y,dt,u_aux,v_aux,u_ant_aux,v_ant_aux)							!Cálculo del RHS ecuación de momentum para v*

	!Se usa el método ADI:
	do j=1,Ny-1																							!Matriz tridiag en x
		call tridiag(1.0d0+dt/(Re*deltax**2.0d0),-dt/(2.0d0*Re*deltax**2.0d0),rhs_x(:,j),Nx-1,uu(:,j))	!1er ADI en x
		call tridiag(1.0d0+dt/(Re*deltax**2.0d0),-dt/(2.0d0*Re*deltax**2.0d0),rhs_y(:,j),Nx-1,vv(:,j))	!1er ADI en y
	enddo
	do i=1,Nx-1																							!Matriz tridiag en y	
		call tridiag(1.0d0+dt/(Re*deltay**2.0d0),-dt/(2.0d0*Re*deltay**2.0d0),uu(i,:),Ny-1,uu(i,:))		!2do ADI en x
		call tridiag(1.0d0+dt/(Re*deltay**2.0d0),-dt/(2.0d0*Re*deltay**2.0d0),vv(i,:),Ny-1,vv(i,:))		!2do ADI en y
	enddo
	
	!Se despeja v*
	uu = uu + u
	vv = vv + v

	call divergencia(uu,vv,deltax,deltay,div)
	call edo_fourier(deltax,deltay,Phi,div*dt**(-1.0d0))													!Cálculo de Phi
	P = P + Phi - 0.5d0*div*Re**(-1.0d0)																!Cálculo de la presión en el tiempo n+1

	if (ep == 1) then
		call auxiliares(u,Phi,u_aux,aux_Phi)
		call laplac(aux_phi,deltax,deltay,lap_phi)
		write(*,*)'Error Poisson:',maxval(abs(dt*lap_phi - div))
	endif

	call despejar_velocidad(uu,vv,u,v,u_aux,v_aux,u_ant_aux,v_ant_aux,Phi,deltax,deltay,dt)				!Cálculo de v en el tiempo n+1
	call vorticidad(u_aux,v_aux,deltax,deltay,vort)														!Cálculo de la vorticidad

	!Cálculo de máximos para todo el programa
	if (max_u<maxval(u)) then
		max_u = maxval(u)
	endif
	if (max_v<maxval(v)) then
		max_v = maxval(v)
	endif
	if (max_p<maxval(p)) then
		max_p = maxval(p)
	endif
	if (max_vort<maxval(vort)) then
		max_vort = maxval(vort)
	endif
	
	!Mecanismo para detener el programa si los valores se escapan (u > 100)
	if (maxval(abs(u))>100.0d0) then
		write(*,*) 
		write(*,*) 'EL CÓDIGO EXPLOTÓ'
		write(*,*) 
		exit
	endif

	t_iter = t_iter + dt												!Actualizamos al tiempo actual
	!Mecanismo para hacer gráficos
	if (t_iter*t_mod**(-1.0d0)*aux_n>cont) then
		cont = cont + 1
		call graficar(u_aux,v_aux,P,vort,cont,t_iter,deltax,deltay,gp,gv,gw)
	endif
enddo

call armar_gif															!Arma los archivos .gif
call graficar_cfl(dt_cfl,n,tiempo)										!Creación de gráfico para dt

call cpu_time(finish)	
!Escribe por pantalla los máximos obtenidos 
write(*,'(a19)') 'Resultados Finales:'
write(*,'(a40,F8.5)') 'Velocidad Horizontal Máxima:            ',max_u
write(*,'(a40,F8.5)') 'Velocidad Vertical Máxima:              ',max_v
write(*,'(a41,F8.5)') 'Presión Máxima:                         ',max_p
write(*,'(a40,F8.4)') 'Vorticidad Máxima:                      ',max_vort
write(*,*)
write(*,'(a39,F10.4,a4)') 'Tiempo Transcurrido:                   ',finish-start,' [s]'
write(*,*)
call system ('bash red.sh')
write(*,'(a36)') 'El programa finalizó correctamente!'
write(*,*)
call system ('bash norm.sh')
call system('aplay -q end.wav')											!Fin del programa, reproduce un audio
end program
