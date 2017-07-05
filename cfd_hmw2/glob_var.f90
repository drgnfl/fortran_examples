!Valores del flujo
double precision u(nx+1,ny), v(nx,ny+1), p(nx,ny)
double precision uu(nx+1,ny), vv(nx,ny+1), pp(nx,ny)
double precision u_ae(nx,ny), u_aw(nx,ny), u_as(nx,ny), u_an(nx,ny), u_anb(nx,ny), beta_x(nx+1,ny)
double precision v_ae(nx,ny), v_aw(nx,ny), v_as(nx,ny), v_an(nx,ny), v_anb(nx,ny), beta_y(nx,ny+1)
double precision beta_p(nx,ny)
double precision div(nx,ny)
double precision dp(nx,ny)
double precision w(nx,ny)
double precision u_new(nx+1,ny), v_new(nx,ny+1)
double precision error_pp, max_div
integer k, n
!Variables para el clock
double precision rate
integer c1, c2, cr, cm
integer qq
!Variables asociadas a dt
double precision dt, dt_cfl(20000), tiempo(20000), t_iter

