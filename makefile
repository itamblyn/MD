msd: msd.f90
	ifort msd.f90 -o bin/msd.x -O3 -stand f95

unwrap: unwrap_PBC.f90
	ifort unwrap_PBC.f90 -o bin/unwrap_PBC.x -O3 -stand f95

wrap: wrap_PBC.f90
	ifort wrap_PBC.f90 -o bin/unwrap_PBC.x -O3 -stand f95

clean:
	rm bin/*.x
