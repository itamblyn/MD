!       program replicate_PBC.f
!***********************************************************
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, i, j, ix, iy, iz
        parameter (n=1000000)
        integer natom
        real*4 x, y, z
        integer nx, ny, nz
        real*4 ax, ay, az
        character*2 typat
        character*128 fin

        ! Get xyz fname & lattice constants from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fin
        write(*,*) 'Number of images: nx, ny, nz'
        read(*,*) nx, ny, nz

        ! Open input and output xyz files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='replicated.xyz')

        ! Loop over the file
        do i=1,n

          ! Read and write natom and timestep
          read(1,*,END=100) natom
          read(1,*) ax, ay, az
          write(2,*) natom*nx*ny*nz
          write(2,*) ax, ay, az

          do j=1,natom

            ! Read in coordinates from xyz file
            read(1,*,END=100) typat, x, y, z

            do ix=1,nx
              do iy=1,ny
                do iz=1,nz
                  write(2,*) typat, (ix-1)*ax + x,(iy-1)*ay + y,(iz-1)*az + z
                enddo
              enddo
            enddo
          enddo

        enddo

 100    continue

        close(1)

 90     continue

        END
