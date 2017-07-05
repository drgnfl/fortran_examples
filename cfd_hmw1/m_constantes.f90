!Este módulo contiene los parámetros del problema, además de sumarle o no la perturbación.
module m_constantes
implicit none
integer, parameter :: Nx = 64, Ny = 64, perturbacion = 1, aux_n = 100							!aux_n: numero de frames
integer, parameter :: gp = 1, gv = 1, gw = 1													!Activa los gráficos de P, V y Vort
integer, parameter :: ep = 0																	!Muestra el error de Poisson en Pantalla
integer, parameter :: vp = 0																	!Activa el validador de Poisson
double precision, parameter :: Re = 1000.0d0, cfl=0.5d0, pi = 4.0d0*atan(1.0d0), t_mod = 5.0d0, lambda = 1.0d0
end module m_constantes
