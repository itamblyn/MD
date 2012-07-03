!       program msd_com.f
!***********************************************************
!       Calculate the mean squared displacement from
!       an xyz file
!         
!       This is an edit of B.Boates's msd code. This edit
!       IGNORES the timestep in the xyz and subtracts off the COM
!       (since some MD codes fail to do this automatically)
!
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer n, m, q, i, j, k
        parameter (n=1000000,m=1024)
        integer natoms(n), timesteps(n)
        real(8) Rx0(m), Ry0(m), Rz0(m)
        real(8) Rx(m), Ry(m), Rz(m)
        real(8) dRx, dRy, dRz, sum_msd, msd(n)
        real(8) sx, sy, sz
        real(8) D(n)
        character(128) fin, dummy

        write(6,*) 'Name of unwrapped xyz file'
        read(5,*) fin

        open(1,file=fin,status='old',ERR=90)

        sx = 0
        sy = 0
        sz = 0

        read(1,*) natoms(1)
        read(1,*) 
        do j=1,natoms(1)
          read(1,*) dummy, Rx0(j), Ry0(j), Rz0(j)
          sx = sx + Rx0(j)
          sy = sy + Ry0(j)
          sz = sz + Rz0(j)
        enddo        

        sx = sx/natoms(1)
        sy = sy/natoms(1)
        sz = sz/natoms(1)

        do j=1,natoms(1)
          Rx0(j) = Rx0(j) - sx
          Ry0(j) = Ry0(j) - sy
          Rz0(j) = Rz0(j) - sz
        enddo

        open(2,file='msd.dat')
        write(2,*) '# timestep msd (center of mass has been subtracted)'
        do i=2,n

          read(1,*,END=100) natoms(i)
          read(1,*) 

          sum_msd = 0.0

          sx = 0
          sy = 0
          sz = 0

          do j=1,natoms(i)

            read(1,*) dummy, Rx(j), Ry(j), Rz(j)
            sx = sx + Rx(j)
            sy = sy + Ry(j)
            sz = sz + Rz(j)
          enddo

          sx = sx/natoms(i)
          sy = sy/natoms(i)
          sz = sz/natoms(i)

          do j=1,natoms(i)

            dRx = (Rx(j)-sx) - Rx0(j)
            dRy = (Ry(j)-sy) - Ry0(j)
            dRz = (Rz(j)-sz) - Rz0(j)
            sum_msd = sum_msd + (dRx**2 + dRy**2 + dRz**2)

          enddo
          
          msd(i) = sum_msd / natoms(i)

          ! Use timestep of 32 au in ps
!          D(i) = msd(i) / (6.0 * timesteps(i)*0.00077404298)

!          D(1) = 0.0

          write(2,*) i, msd(i) !, D(i)
!          write(2,*) i - 1, msd(i) !, D(i)

        enddo

 100    continue

        close(1)

        close(2)

 90     continue

        END
