msd: msd.f90
	gfortran msd.f90 -o bin/msd.x

unwrap: unwrap_PBC.f90
	gfortran unwrap_PBC.f90 -o bin/unwrap_PBC.x -O3

unwrap_cell: unwrap_PBC_cell.f90
	gfortran unwrap_PBC_cell.f90 -o bin/unwrap_PBC_cell.x -O3

wrap: wrap_PBC.f90
	gfortran wrap_PBC.f90 -o bin/wrap_PBC.x -O3

clean:
	rm bin/*.x
