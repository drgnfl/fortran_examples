module m_cons
implicit none

!Parámetros del Problema
integer, parameter :: llx = 10, lly = 5 !Largo del dominio
integer, parameter :: mx = 8, my = 8    !Factores de aumento de malla (Nx,Ny)=(lx*mx*2,ly*my*2)
double precision, parameter :: re = 2000.0d0, cfl = 1.0d0
double precision, parameter :: alpha_p = 0.8d0	!Factor de sobrerelajación
integer, parameter :: sq = 1	!Activa o desactiva el cuadrado
integer, parameter :: kmax = 2000	!Número máximo de iteraciones
double precision, parameter :: tol_u = 1.0d-3, tol_v = 1.0d-3, tol_p = 1.0d-3, tol_div = 1.0d-3	!Tolerancias para las ecs.
integer, parameter :: tp = 0 !Activa la comparación con perfil teórico
integer, parameter :: gif = 0 !Activa la creación de animaciones

!Parámetros Calculados
double precision, parameter :: Lx = 1.0d0*llx, Ly = 1.0d0*lly
integer, parameter :: nx = mx*2*llx, ny = my*2*lly
double precision, parameter :: dx = Lx/nx, dy = Ly/ny
double precision, parameter :: gam = 1.0d0/re
end module m_cons
