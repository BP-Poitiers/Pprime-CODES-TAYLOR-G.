!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module cavity

   use decomp_2d
   use decomp_2d_mpi
   use variables
   use param

   implicit none

   ! Temperature difference between both walls
   real(mytype), parameter :: deltaT = 1._mytype
   real(mytype), parameter :: temp_l = deltaT/2._mytype

   private ! All functions/subroutines private by default
   public :: init_cavity, boundary_conditions_cavity, postprocess_cavity

contains

   !############################################################################
   !!
   !! Init the cavity case
   !!
   !############################################################################
   subroutine init_cavity(ux1, uy1, uz1, ep1, phi1)

      USE decomp_2d_io
      USE MPI

      implicit none

      real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, ep1
      real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

      real(mytype) :: x
      integer :: i

      ! This does not apply in case of restart
      if (irestart == 0) then
         ! Velocity is zero
         ux1 = zero
         uy1 = zero
         uz1 = zero

         ! Linear temperature profile
         if (numscalar >= 1) then
            do i = 1, xsize(1)
               phi1(i, :, :, 1) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
               phi1(i, :, :, 2) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if

      end if

   end subroutine init_cavity

   !############################################################################
   !!
   !! Boundary conditions for the cavity case
   !!
   !############################################################################
   subroutine boundary_conditions_cavity(ux, uy, uz, phi)

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
             ! call outflow_cav(ux,uy,uz,phi)
      END IF
      IF (nclyn == 2 .and. xend(2) == ny) THEN
             ! call inflow_cav(ux,uy,uz,phi)
              ! call outflow_cav(ux, uy, uz, phi)
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
            phi(1, :, :, 2) = 0
         end if
         ! x=Lx
         if (nclxSn == 2) then
            phi(xsize(1), :, :, 1) = 0
            phi(xsize(1), :, :, 2) = 0
         end if
         ! y=0
         if (nclyS1 == 2 .and. xstart(2) == 1) then
            do i = 1, xsize(1)
               phi(i, 1, :, 1) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
               phi(i, 1, :, 2) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if
         ! y=Ly
         if (nclySn == 2 .and. xend(2) == ny) then
            do i = 1, xsize(1)
               phi(i, xsize(2), :, 1) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
               phi(i, xsize(2), :, 2) = temp_l - deltaT*(i - 1)/real(xsize(1) - 1, kind=mytype)
            end do
         end if
      end if

      ! Clip
      if (numscalar >= 1) then
         do k = 1, xsize(3)
            do j = 1, xsize(2)
               do i = 1, xsize(1)
                  if (phi(i,j,k,1) > temp_l) phi(i,j,k,1) = temp_l
                  if (phi(i,j,k,1) < -temp_l) phi(i,j,k,1) = -temp_l
                  if (numscalar >=2) then
                         if (phi(i,j,k,2) > temp_l ) phi(i,j,k,2) = temp_l
                         if (phi(i,j,k,2) < -temp_l) phi(i,j,k,2) = -temp_l
                  endif
               enddo
            enddo
         enddo
      endif

   end subroutine boundary_conditions_cavity


   ! ##########################################################################
   ! INFLOW AND OUTFLOW REVERSED - INFLOW IN y=Ly
   ! ##########################################################################
  
   subroutine inflow_cav (ux,uy,uz,phi)

      USE param
      USE variables
      USE MPI
      USE var, only: ux_inflow, uy_inflow, uz_inflow

      implicit none

      real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
      real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
      real(mytype) :: ut1

      integer :: i,j,k

      !if (mod(itime, ilist).eq.0) then
      !      print*, "-------------------------"
      !      print*, "Prints de verifications", activate_q
      !      print*, 'Lx star ', ' - dx ', '- xp(nx)'
      !      print*, xlx, dx, xp(xsize(1))
      !      print*, '-=-=-=-=-=-=-'
      !      print*, 'pi ', 'Lx'
      !      print*, pi, xlx
      !      print*, '-=-=-=-=-=-=-'
      !      print*, 'ampli ', '- omega ', ' - period', ' - kx'
      !      print*, ampli_sin, pulse_sin, k_wave, (k_wave*twopi/xlx)
      !      print*, '-=-=-=-=-=-=-'
      !      print*, sin(xp(1)), sin(xp(100)), sin(xp(xsize(1)))
      !      print*, sin(0*dx), sin((100-1)*dx), sin((xsize(1)-1)*dx)
      !      print*, '-------------------------'
      !end if


      do k=1, xsize(3)
          do i=1, xsize(1)
              byyn(i, k) = ampli_sin*sin((k_wave*twopi/xlx)  * xp(i) - pulse_sin * itime * dt)
          enddo
      enddo

      ut1=zero
      do k=1,xsize(3)
          do j=1,xsize(1)
             ut1=ut1+byyn(j,k)   ! In en ly
          enddo
      enddo

      ut1=ut1/(real(nx*nz,mytype))


      if (mod(itime, ilist)==0) then
        print*, '------------ Input parameters ------------'
        print*, 'sinus', byyn(1,1), byyn(100,1), byyn(xsize(1), 1)
        print*, 'debit in', ut1
      end if


   end subroutine inflow_cav


   ! #######################################################################
   ! OUTFLOW IN y=0
   ! #######################################################################

   subroutine outflow_cav (ux,uy,uz,phi)

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
         !         print*, 'Changement max'
          end if
          if (uy(i, 2, k).lt.uymin) then
                  uymin=uy(i, 2, k)
               !   print*, "cahngement min"
          end if
        enddo
      enddo

      call MPI_ALLREDUCE(MPI_IN_PLACE,uymax,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uymin,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)
        
      if (byyn(200, 1)>=0) then
              cy=0.5*(uymax+uymin)*gdt(itr)*udy
              if (mod(itime, ilist)==0) then
                     ! print*, '------- Sinus parameters ------'
                      print*, "Positif, cy", cy
              endif
      else if (byyn(200,1) < 0) then
              cy=-0.5*(uymax+uymin)*gdt(itr)*udy
              if (mod(itime, ilist)==0) then
                     ! print*, '------- Sinus parameters ------'
                      print*, "Negatif, cy ", cy
              endif
      end if

      ! cy=0.5*(uymax+uymin)*gdt(itr)*udy
      do k=1, xsize(3)
        do i=1, xsize(1)
          byx1(i,k)=ux(i, 1, k)-cy*(ux(i, 1, k)-ux(i, 2, k))
          byy1(i,k)=uy(i, 1, k)-cy*(uy(i, 1, k)-uy(i, 2, k))
          byz1(i,k)=uz(i, 1, k)-cy*(uz(i, 1, k)-uz(i, 2, k))
          if (iscalar.eq.1) then
            phi(i, 1, k, :) = phi(i, 1, k, :)-cy*(phi(i, 1, k, :)-phi(i, 2, k, :))
          endif
        enddo
      enddo

      if ((nrank==0).and.(mod(itime,ilist)==0)) then
        print*, '--------------------'
       ! write(*,*) 'sinus original', byyn(200,1)
        write(*,*) 'cy', cy, 0.5*(uymax+uymin), gdt(itr)
       ! write(*,*) 'BY1 BC', real(byy1(200,1),4), real(byy1(100,1),4), real(byy1(256,1),4), real(byy1(400,1),4)
      end if


   end subroutine outflow_cav

   !############################################################################
   !!
   !! Post-processing for the cavity case
   !!
   !############################################################################
   subroutine postprocess_cavity(ux1, uy1, uz1, phi1)

      implicit none

      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

   end subroutine postprocess_cavity

end module cavity
