! <DefPhotolysis_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 

!+ Photolysis coefficients
!-------------------------------------------------------------------------------

     module DefPhotolysis_ml
!-------------------------------------------------------------------------------

!    Data needed  for photolysis calculation.  NPHODIS is the number of
!    tabulated rates from the Phodis model. NRCPHOT (<=NPHODIS) is the
!    number of photolysis rats needed by the model
!
!   10/10/01 - corrected and tidied up by jej.
!   11/10/01 - NDISS removed, minor F90 changes and docs
!              added. NLAT and CLOUDTOP added.
!-------------------------------------------------------------------------------

   use CheckStop_ml,      only: CheckStop
   use GridValues_ml    , only : gb
   use Io_ml,           only : IO_DJ, open_file, ios
   use Met_ml           , only : cc3d,cc3dmax,z_bnd
   use ModelConstants_ml,    only: KMAX_MID, KCHEMTOP, NPROC
   use Par_ml      ,    only : me,MAXLIMAX,MAXLJMAX
   use LocalVariables_ml, only : Grid  ! => izen
   implicit none
   private

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  integer, public, parameter :: &
             NRCPHOT      = 17   ! Number of photolytic reactions
   
   real, public, dimension(NRCPHOT,KCHEMTOP:KMAX_MID), &
        save :: rcphot       ! photolysis rates    -   main output

   real, public, save :: sum_rcphot     ! was jej's sum1, for debug only
   logical, public, parameter :: DEBUG_DJ = .false.

   integer, parameter, private ::  &
               HORIZON    = 90     & ! Integer solar zenith angle at sunset
             , CLOUDTOP   = 6        ! k-value above which clear-sky dj assumed
                                     ! (since..... Joffen?)

   integer, parameter, private :: &
               NPHODIS = 17       &  ! Max possible NRCPHOT
              ,NLAT    = 6           ! No. latitude outputs

    real, private, dimension(NPHODIS,KCHEMTOP:KMAX_MID,HORIZON,NLAT) :: dj

    real, private, dimension(NPHODIS,KCHEMTOP:KMAX_MID,HORIZON) :: &
                   djcl1        &
                  ,djcl3

!  Indices of photolysis rates as available from Phodis files:

    integer, public, parameter ::  &
      IDAO3    =  1 , IDBO3    =  2 , IDNO2    =  3 , &
      IDH2O2   =  4 , IDHNO3   =  5 , IDACH2O  =  6 , &
      IDBCH2O  =  7 , IDCH3CHO =  8 , IDCH3COX =  9 , &
      IDCH3COY = 10 , IDHCOHCO = 11 , IDRCOHCO = 12 , &
      IDNO3    = 13 , IDN2O5   = 14 , IDCH3O2H = 15 , &
      IDHO2NO2 = 16 , IDACETON = 17

 !/ subroutines: 

  public :: readdiss
  public :: setup_phot

 contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine readdiss(newseason)

      integer ::  newseason

      integer ::  k     &  ! help index
                 ,izn   &  ! integer zenith angle
                 ,info  &  ! used for broadcast
                 ,nr    &  ! numbering of photolytic reactions
                 ,la       ! counting every 10 deg. latitude
        real myz
        character*20 fname1, fname2, fname3





!    Open, read and broadcast clear sky rates
!---------------

        if(me == 0)then
           write(fname1,fmt='(''jclear'',i2.2,''.dat'')') newseason
           call open_file(IO_DJ,"r",fname1,needed=.true.)
           call CheckStop(ios,"DefPhotolysis: ios error in jclear ")
        endif


!       Format of input data from Phodis - careful with "17" and NPHODIS
999     FORMAT(1x, f8.3, 17(1x, 1pe8.2))


        if(me == 0)then

          do la = 1,NLAT
            do izn = 1,HORIZON
              do k = 1,KCHEMTOP
                read(IO_DJ,999) myz,(dj(nr,KCHEMTOP,izn,la),nr=1,NPHODIS)
              end do
              do k = KCHEMTOP+1,KMAX_MID
                read(IO_DJ,999) myz,(dj(nr,k,izn,la),nr=1,NPHODIS)
              end do   ! k
            end do    ! izn
          end do     ! la
          close(IO_DJ)
        endif  ! me = 0

        CALL MPI_BCAST(dj  ,8*NPHODIS*(KMAX_MID-KCHEMTOP+1)*HORIZON*NLAT,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 




!    Open, read and broadcast light cloud rates
!---------------

        if(me == 0)then
           write(fname2,fmt='(''jcl1km'',i2.2,''.dat'')') newseason
           call open_file(IO_DJ,"r",fname2,needed=.true.)
           call CheckStop(ios,"DefPhotolysis: ios error in jcl1km ")
        endif


        if(me == 0)then

          do izn = 1,HORIZON
            do k = 1,KCHEMTOP
              read(IO_DJ,999) myz,(djcl1(nr,KCHEMTOP,izn),nr=1,NPHODIS)
            end do
            do k = KCHEMTOP+1,KMAX_MID
              read(IO_DJ,999) myz,(djcl1(nr,k,izn),nr=1,NPHODIS)
            end do
          end do  ! izn

          do izn = 1,HORIZON
            do k = KCHEMTOP,KMAX_MID
              do nr=1,NPHODIS
                djcl1(nr,k,izn)=djcl1(nr,k,izn)/dj(nr,k,izn,3)-1.0
              enddo ! nr
            end do ! k
          end do  ! izn
          close(IO_DJ)
        endif   ! me = 0

        CALL MPI_BCAST(djcl1  ,8*NPHODIS*(KMAX_MID-KCHEMTOP+1)*HORIZON,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 




!    Open, read and broadcast dense cloud rates
!---------------

        if(me == 0)then
           write(fname3,fmt='(''jcl3km'',i2.2,''.dat'')') newseason
           call open_file(IO_DJ,"r",fname3,needed=.true.)
           call CheckStop(ios,"DefPhotolysis: ios error in jcl3km ")
        endif


        if(me == 0)then

          do izn = 1,HORIZON
            do k = 1,KCHEMTOP
              read(IO_DJ,999) myz,(djcl3(nr,KCHEMTOP,izn),nr=1,NPHODIS)
            end do
            do k = KCHEMTOP+1,KMAX_MID
              read(IO_DJ,999) myz,(djcl3(nr,k,izn),nr=1,NPHODIS)
            end do  ! k
          end do   ! izn
          close(IO_DJ)
          do izn = 1,HORIZON
            do k = KCHEMTOP,KMAX_MID
              do nr=1,NPHODIS
                djcl3(nr,k,izn)=djcl3(nr,k,izn)/dj(nr,k,izn,3)-1.
              enddo  ! nr
            end do  ! k
          end do   ! izn
       endif      !  me = 0

        CALL MPI_BCAST(djcl3  ,8*NPHODIS*(KMAX_MID-KCHEMTOP+1)*HORIZON,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 


!       if(me == 0) then
!        do k=1,KMAX_MID
!           write(6,*) 'jverdi i niv. k',k
!           write(6,*) (dj(3,1,k,nr),nr=1,4)
!           write(6,*) (djcl1(1,k,nr),nr=1,4)
!           write(6,*) (djcl3(1,k,nr),nr=1,4)
!         end do
!       end if

        return

        end subroutine readdiss
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        subroutine setup_phot(i,j,errcode)

!       input
        integer :: i,j
        integer :: errcode

!       local
        integer la     &  ! counting every 10 deg. latitude
               ,n      &  ! help index
               ,k      &  ! vertical index
               ,base   &  ! cloud base
               ,top    &  ! cloud top
               ,iclcat    ! cloud type

        real clear        ! clear sky fraction
        real sum1

!---- assign photolysis rates ------------------------------------------------

        errcode = 0


        if ( Grid%izen > 90 ) then     ! Photolysis rates zero when the sun is
                                  ! below the horizon
             rcphot(:,:) = 0.0 

        else !! (izen < 90)  -- sun above horizon


        !/ first find cloud base and cloud top

           iclcat = 0
           if(cc3dmax(i,j,KMAX_MID) > 1.e-4) then

              k = KMAX_MID
              do while(cc3d(i,j,k) < 1.e-4 .and. k >= CLOUDTOP)
                 k = k-1
              end do
              base = k+1

             ! if all cc3d so far are <1.e-4 we are done

              if( base < CLOUDTOP ) then 

                !  we have found a k>=CLOUDTOP with cc3d>=1.e-4, now search for top

                 k = CLOUDTOP
                 do while(cc3d(i,j,k) < 1.0e-4)
                     k = k+1
                 end do
                 top = k

                 if(top >= base) then
                   print *,'top,base'
                   errcode = 17
                   return
                 endif
                 iclcat = 1

                 if(z_bnd(i,j,top)-z_bnd(i,j,base) > 1.5e3) iclcat = 2

              end if  ! base<CLOUDTOP
            end if   ! end cc3dmax
4732        continue




            la = max(1,int(0.1*gb(i,j)-2.0001))

            if(iclcat == 0)then
              do k = KCHEMTOP,KMAX_MID
                do n=1,NRCPHOT
                  rcphot(n,k)  = dj(n,k,Grid%izen,la)
                enddo
              enddo
            else if(iclcat == 1)then
              clear = cc3dmax(i,j,KMAX_MID)
              do k = KCHEMTOP,KMAX_MID
                do n=1,NRCPHOT
                  rcphot(n,k)  = (1. +          &
                               clear*djcl1(n,k,Grid%izen)) * dj(n,k,Grid%izen,la)
                enddo  !  n
              enddo   !  k


            else
              clear = cc3dmax(i,j,KMAX_MID)
              do k = KCHEMTOP,KMAX_MID
                do n=1,NRCPHOT
                  rcphot(n,k)  = (1. +            &
                               clear*djcl3(n,k,Grid%izen))*dj(n,k,Grid%izen,la)
                 enddo
              enddo
            endif
  
            if ( DEBUG_DJ ) then
                sum_rcphot = sum_rcphot + &
                      sum ( rcphot(1:NRCPHOT, KCHEMTOP:KMAX_MID) )
            end if


          end if   !  end izen <  90 (daytime)  test

    end subroutine setup_phot
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 end module DefPhotolysis_ml
