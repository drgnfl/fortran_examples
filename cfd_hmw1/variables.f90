!Archivo de texto que contiene las variables globales del programa.
integer i, j, cont, n
double precision deltax, deltay, dt
double precision u(Nx-1,Ny-1), v(Nx-1,Ny-1), P(Nx-1,Ny-1)
double precision u_aux(Nx+1,Ny+1), v_aux(Nx+1,Ny+1)
double precision u_ant_aux(Nx+1,Ny+1), v_ant_aux(Nx+1,Ny+1)
double precision u_ant(Nx-1,Ny-1), v_ant(Nx-1,Ny-1)
double precision rhs_x(Nx-1,Ny-1), rhs_y(Nx-1,Ny-1)
double precision uu(Nx-1,Ny-1), vv(Nx-1,Ny-1)
double precision Phi(Nx-1,Ny-1), div(Nx-1,Ny-1), vort(Nx-1,Ny-1)
double precision max_p, max_v, max_u, max_vort
double precision dt_cfl(20000), t_iter
