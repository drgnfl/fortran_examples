objects = m_constantes.o m_iniciar.o m_rhs.o m_tridiagonal.o m_presion.o m_desp_veloc.o m_graficos.o main.o 
modules = m_constantes.mod m_iniciar.mod m_rhs.mod m_tridiagonal.mod m_presion.mod m_desp_veloc.mod m_graficos.mod

all:
	make run --quiet

run: cfd1
	clear
	make clean
	clear
	./cfd1
	rm -f cfd1

cfd1: main.o
	gfortran -o cfd1 $(objects)
main.o: main.f90 $(modules)
	gfortran -c main.f90

m_desp_veloc.mod: m_desp_veloc.o
	gfortran -c m_desp_veloc.f90
m_desp_veloc.o: m_desp_veloc.f90 m_constantes.mod m_iniciar.mod
	gfortran -c m_desp_veloc.f90

m_presion.mod: m_presion.o
	gfortran -c m_presion.f90
m_presion.o: m_presion.f90 m_tridiagonal.mod m_iniciar.mod
	gfortran -c m_presion.f90

m_tridiagonal.mod: m_tridiagonal.o
	gfortran -c m_tridiagonal.f90
m_tridiagonal.o: m_tridiagonal.f90 m_constantes.mod
	gfortran -c m_tridiagonal.f90

m_rhs.mod: m_rhs.o
	gfortran -c m_rhs.f90
m_rhs.o: m_rhs.f90 m_constantes.mod
	gfortran -c m_rhs.f90

m_iniciar.mod: m_iniciar.o
	gfortran -c m_iniciar.f90
m_iniciar.o: m_iniciar.f90 m_constantes.mod
	gfortran -c m_iniciar.f90

m_graficos.mod: m_graficos.o
	gfortran -c m_graficos.f90
m_graficos.o: m_graficos.f90 m_constantes.mod
	gfortran -c m_graficos.f90

m_constantes.mod: m_constantes.o
	gfortran -c m_constantes.f90
m_constantes.o: m_constantes.f90
	gfortran -c m_constantes.f90

clean:
	rm -f *.o *.mod
	echo Se ha limpiado la carpeta!
