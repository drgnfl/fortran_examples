objects = m_cons.o m_tdma.o m_solve_uv.o m_solve_dp.o m_graph.o m_coef.o main.o 
modules = m_cons.mod m_tdma.mod m_solve_uv.mod m_solve_dp.mod m_graph.mod m_coef.mod

all:
	make run --quiet

run: cfd2
	make clean
	clear
	./cfd2
	rm -f cfd2

cfd2: main.o
	gfortran -o cfd2 $(objects)
main.o: main.f90 $(modules)
	gfortran -c main.f90

m_solve_dp.mod: m_solve_dp.o
	gfortran -c m_solve_dp.f90
m_solve_dp.o: m_solve_dp.f90 m_cons.mod m_tdma.mod
	gfortran -c m_solve_dp.f90

m_solve_uv.mod: m_solve_uv.o
	gfortran -c m_solve_uv.f90
m_solve_uv.o: m_solve_uv.f90 m_cons.mod m_tdma.mod
	gfortran -c m_solve_uv.f90

m_coef.mod: m_coef.o
	gfortran -c m_coef.f90
m_coef.o: m_coef.f90 m_cons.mod
	gfortran -c m_coef.f90

m_graph.mod: m_graph.o
	gfortran -c m_graph.f90
m_graph.o: m_graph.f90 m_cons.mod
	gfortran -c m_graph.f90

m_tdma.mod: m_tdma.o
	gfortran -c m_tdma.f90
m_tdma.o: m_tdma.f90
	gfortran -c m_tdma.f90

m_cons.mod: m_cons.o
	gfortran -c m_cons.f90
m_cons.o: m_cons.f90
	gfortran -c m_cons.f90

clean:
	rm -f *.o *.mod
