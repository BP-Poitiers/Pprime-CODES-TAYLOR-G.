!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module stratified

   use decomp_2d
   use decomp_2d_mpi
   use variables
   use param

   implicit none

   ! Temperature difference between both walls
   real(mytype), parameter :: deltaT = 1._mytype
   real(mytype), parameter :: temp_l = deltaT/2._mytype

   private ! All functions/subroutines private by default
   public :: init_strat, boundary_conditions_strat, postprocess_strat, scalar_forcing_strat, inflow_strat, interp1d_N, outflow_strat

contains

   !############################################################################
   !!
   !! Init the stratified case
   !!
   !############################################################################
   subroutine init_strat(ux1, uy1, uz1, ep1, phi1)

      USE decomp_2d_io
      USE MPI
      use variables

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
      !real(mytype), dimension(xsize(2)) :: Nn

      integer :: n, statut, unit, iii
      real(mytype),dimension(:),allocatable :: Yi, Ni
      !real(mytype), allocatable :: Yy(:), Nn(:)
      logical :: existe


      real(mytype) :: x
      integer :: i, j, k

      !if (irestart==1) then
      !        print*, "restart passe bien par là"
      !endif


      ! This does not apply in case of restart
      if (irestart == 0) then
         ! Velocity is zero
         ux1 = zero
         uy1 = zero
         uz1 = zero

         ! Linear temperature profile
         if (numscalar >= 1) then
            do i = 1, xsize(1)
               phi1(i, :, :, 1) = 0
            end do
         end if

      end if

     ! ============= READ THE FILE TO SET UP THE BrV FREQ LATTICE =========================

      inquire(file=trim(N_file), exist=existe)
      if (.not. existe) then
              print*, "Fichier introuvable"
              stop
      elseif(existe) then
         if (mod(itime, ilist).eq.0) then
              print*, "Le fichier est trouvé : ", N_file
              print*, '------------------------'
         end if
      endif

      n = 10000
      allocate(Yi(n), Ni(n))
      unit=5555
      open(unit=unit, file=N_file, status='old', action='read', iostat=statut)

      iii=0
      do
        iii=iii+1
        read(unit, *, iostat=statut) Yi(iii), Ni(iii)
        if (statut /=0) exit
      enddo

      allocate(Yy(iii-1), Nn(iii-1))

      Yy(:) = Yi(1:iii-1)
      Nn(:) = Ni(1:iii-1)

      deallocate(Yi, Ni)

      close(unit)

      !if (size(Nn).eq.xsize(2)) then
      !        if (mod(itime, ilist).eq.0) then
      !                print*, "N(y) est de la bonne taille", size(Nn), xsize(2)
      !                print*, Nn(1), Nn(xsize(2))
      !                print*, "-------------------------"
      !        end if
      !else if (size(Nn).ne.xsize(2)) then
      !        print*, "Pas de la bonne taille", size(Nn), xsize(2)
              !stop
      !end if

   end subroutine init_strat

   !!----------------------------------------------------------------------------
   !! Boundary conditions for the stratified case
   !!----------------------------------------------------------------------------
 
   subroutine boundary_conditions_strat(ux, uy, uz, phi)

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux, uy, uz
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi

      integer :: i, j, k

      ! Velocity
      IF (nclx1 == 2) THEN
      END IF
      IF (nclxn == 2) THEN
      END IF
      IF (ncly1 == 2) THEN
      !        call outflow_strat(ux,uy,uz,phi)
      END IF
      IF (nclyn == 2) THEN
              call inflow_strat(ux,uy,uz,phi)
      END IF
      IF (nclz1 == 2) THEN
      END IF
      IF (nclzn == 2) THEN
      END IF

      ! Scalar
      if (numscalar >= 1) then
         ! x=0
         if (nclxS1 == 2) then
            phi(1, :, :, 1) = 0
         end if
         ! x=Lx
         if (nclxSn == 2) then
            phi(xsize(1), :, :, 1) = 0
         end if
         ! y=0
         if (nclyS1 == 2 .and. xstart(2) == 1) then
            do i = 1, xsize(1)
               phi(i, 1, :, 1) = 0
            end do
         end if
         ! y=Ly
         if (nclySn == 2 .and. xend(2) == ny) then
            do i = 1, xsize(1)
               phi(i, xsize(2), :, 1) = 0
            end do
         end if
      end if

   end subroutine boundary_conditions_strat

   !!--------------------------------------------------------------------------------------
   !! Inflow in y=Ly - displacement
   !!-------------------------------------------------------------------------------------

   subroutine inflow_strat (ux,uy,uz,phi)

      USE param
      USE variables
      USE MPI
      USE var, only: ux_inflow, uy_inflow, uz_inflow

      implicit none

      real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
      real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
      real(mytype) :: ut1, k_sin, raideur, vphi, ccentre1, ccentre2

      integer :: i,j,k

      k_sin = (k_wave*twopi/xlx)
      raideur = 0.05
      vphi=pulse_sin / k_sin
      ccentre2=xlx * centre2
      ccentre1=xlx * centre1
      !centreT=vphi*itime*dt

      

      do k=1, xsize(3)
          do i=1, xsize(1)
              byyn(i, k) = ampli_sin*sin(k_sin * xp(i) - pulse_sin * itime * dt) * (-1 /(1+exp(-raideur * (xp(i)-ccentre2))) + 1) * (1 /(1+exp(-0.05 * (xp(i)-ccentre1))))
              
          enddo
      enddo

   end subroutine inflow_strat


  !*******************************************************************************
  !
   subroutine scalar_forcing_strat(uy1,dphi1,phi1)
  !
  !*******************************************************************************

      USE param
      USE variables

      implicit none

      real(mytype),dimension(xsize(1),xsize(2),xsize(3), ntime) :: dphi1
      real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uy1, phi1

      integer :: i, j, k, indice
      real(mytype) :: y_glob, dy_lat, alpha
      real(mytype), dimension(xsize(2)) :: N_from_y

      dy_lat = (Yy(2) - Yy(1))

      do j=1, xsize(2)
        y_glob= (xstart(2) + j - 1)*dy
        
        indice = int(y_glob / dy_lat +1)
        alpha = (y_glob - (indice-1) * dy_lat) / dy_lat

        if (indice<1) then
                N_from_y(j) = Nn(1)
        else if (indice >= size(Yy)) then
                N_from_y(j) = Nn(size(Nn))
        else
                N_from_y(j) = Nn(indice) + alpha * (Nn(indice+1) - Nn(indice))
        end if
 
      end do

      ! Terms from decomposition
      do i=1, xsize(1)
        do j=1, xsize(2)
            do k=1, xsize(3)
                !y_glob = (xstart(2) + j - 1)*dy
                !call interp1d_N(y_glob, N_from_y)
                dphi1(i,j,k,1) = dphi1(i,j,k,1) - (1/ri(1))*N_from_y(j)*N_from_y(j)*uy1(i,j,k)
            end do
        end do
      end do

   end subroutine scalar_forcing_strat


   !!----------------------------------------------------------------------------
   !! Interpolation of the BrV lattice for a given y_global - Archaic
   !!----------------------------------------------------------------------------
   
   subroutine interp1d_N(y_pos, N_from_y)

      USE param
      USE variables
      USE MPI

      implicit none

      real(mytype), intent(in) :: y_pos
      real(mytype), intent(out) :: N_from_y
      integer :: iii

      Ny = 0.0d0

      do iii=1, size(Nn)
        if (y_pos >= Yy(iii) .and. y_pos <= Yy(iii+1)) then
                N_from_y = Nn(iii) + ( Nn(iii+1)-Nn(iii) ) * (y_pos - Yy(iii)) / (Yy(iii+1) - Yy(iii))
                return
        end if
      end do

      if (y_pos<Yy(1)) then
              N_from_y = Nn(1)
      else if (y_pos>Yy(size(Nn))) then 
              N_from_y = Nn(size(Nn))
      end if

   end subroutine interp1d_N


   !##########################################################################
   !##########################################################################

   subroutine outflow_strat (ux,uy,uz,phi)

      USE param
      USE variables
      USE MPI

      implicit none

      integer :: i,k,code
      real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
      real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
      real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cy,uymin,uymax

      udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

      uymax=-1609._mytype
      uymin=1609._mytype
      do k=1, xsize(3)
        do i=1, xsize(1)
          if (uy(i, 2, k).gt.uymax) then
                  uymax=uy(i, 2, k)
          end if
          if (uy(i, 2, k).lt.uymin) then
                  uymin=uy(i, 2, k)
          end if
        enddo
      enddo

      !call MPI_ALLREDUCE(MPI_IN_PLACE,uymax,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
      !call MPI_ALLREDUCE(MPI_IN_PLACE,uymin,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)
        
      !do k=1, xsize(3)
      !  do i=1, xsize(1)
      !    cy = uy(i, 2, k) * gdt(itr)*udy
      !    byx1(i,k)=ux(i, 1, k)-cy*(ux(i, 1, k)-ux(i, 2, k))
      !    byy1(i,k)=uy(i, 1, k)-cy*(uy(i, 1, k)-uy(i, 2, k))
      !    byz1(i,k)=uz(i, 1, k)-cy*(uz(i, 1, k)-uz(i, 2, k))
      !    if (iscalar.eq.1) then
      !      phi(i, 1, k, :) = phi(i, 1, k, :)-cy*(phi(i, 1, k, :)-phi(i, 2, k, :))
      !    endif
      !  enddo
      !enddo


   end subroutine outflow_strat


   !############################################################################
   !!
   !! Post-processing for the stratified case
   !!
   !############################################################################
   subroutine postprocess_strat(ux1, uy1, uz1, phi1)

      implicit none

      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

   end subroutine postprocess_strat

end module stratified
