      program molecules

        implicit none
   
        real, parameter :: pi = 3.1415927
        real, parameter :: e0 = 8.8541878e-12
        real, parameter :: coxygen = -0.8
        real, parameter :: chydrogen = 0.4
        real, parameter :: echarge = 1.6021765e-19
        real, parameter :: blength = 95.84e-12
        real, parameter :: anglehoh = 104.45
        integer, parameter :: mmolecules = 3
        integer, parameter :: apm = 3 
        real :: box
        integer :: nsimulations, i, count
        real :: tforce, avgforce
        real :: force

        call System_Clock(count)
        call random_seed(count)

        write(*,*) 'Size of the box in meters:'
        read(*,*) box
        write(*,*) 'Number of simulations:'
        read(*,*) nsimulations

        tforce = 0.0

        do i = 1, nsimulations
          force = 0.0
          call simwatermolecules(box, force)
          tforce = tforce + force
        END do

        avgforce = tforce / nsimulations

        write(*,*) 'The average force experienced by an atom:', avgforce

      contains
        
        subroutine simwatermolecules(box, force)
          real, intent(in) :: box
          real, intent(out) :: force
          real, dimension(mmolecules * apm) :: xpos, ypos, charges
          integer :: i, j, atom
          real :: dx, dy, r, fout, ox, oy, angle

          xpos = 0.0
          ypos = 0.0
          force = 0.0
          charges = 0.0

          do i = 1, mmolecules * apm
            if(mod(i-1, apm) == 0) then
              charges(i) = coxygen
            else
              charges(i) = chydrogen
            END if
          END do

          do atom = 1, mmolecules * apm, apm
            call random_number(ox)
            call random_number(oy)
            ox = ox * box
            oy = oy * box
            xpos(atom) = ox
            ypos(atom) = oy

            do j = 1, 2
              call random_number(angle)
              angle = angle * 2 * pi
              xpos(atom + j) = ox + blength * cos(angle)
              ypos(atom + j) = oy + blength * sin(angle)
            END do
          END do

        force = 0.0

        do i = 1 , mmolecules * apm
          do j = i+1, mmolecules * apm
              dx = xpos(j) - xpos(i)
              dy = ypos(j) - ypos(i)
              r = sqrt(dx**2 + dy**2)

              if(r .gt. 0.0) then
                call coulomb(charges(i) * echarge, 
     &          charges(j) * echarge, r, fout)

                force = force + fout
              END if
            END do
          END do

        force = force / (mmolecules * apm * 
     &  (mmolecules * apm - 1) / 2)

        END subroutine simwatermolecules

        subroutine coulomb(q1, q2, r, fout)
          real, intent(in) :: q1, q2, r
          real, intent(out) :: fout
          fout = (q1 * q2) / (4.0 * pi * e0 * r**2)
        END subroutine coulomb

      END program molecules
