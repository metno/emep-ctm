! <DefPhotolysis_ml.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.17>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2018 met.no
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
! <DefPhotolysis_ml.f90 - A component of the EMEP MSC-W Chemical transport Model>
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
   use DerivedFields_ml,  only: d_3d, f_3d
   use GridValues_ml   , only : glat, A_bnd, B_bnd
   use Io_ml,            only : IO_DJ, open_file, ios
   use LocalVariables_ml,only : Grid  ! => izen
   use MetFields_ml    , only : cc3d,cc3dmax,z_bnd,ps
   use Config_module,     only: TXTLEN_FILE, KMAX_MID, KCHEMTOP, NPROC, IOU_INST,&
                                   jcl1kmFile,jcl3kmFile,jclearFile,num_lev3d,lev3d
   use MPI_Groups_ml   , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_INTEGER&
                                     ,MPI_COMM_CALC, IERROR
   use Par_ml      ,     only : me,LIMAX,LJMAX
   use SmallUtils_ml,     only: key2str, find_index
   use Functions_ml,      only: StandardAtmos_km_2_kPa
   implicit none
   private

   integer, public, parameter :: &
             NRCPHOT      = 17   ! Number of photolytic reactions
   
   integer, public, parameter:: NzPHODIS=20 !number of heights defined in the input files
   real, save, public :: zPHODIS(NzPHODIS) !heights of the input files, in km assumed constants
   real, save, public :: P_PHODIS(NzPHODIS) !heights of the input files, in Pa

   real, allocatable,save,public, dimension(:,:) &
         :: rcphot       ! photolysis rates    -   main output

   real, public, save :: sum_rcphot     !  for debug only
   logical, public, parameter :: DEBUG_DJ = .false.

   integer, parameter, private ::  &
               HORIZON    = 90     & ! Integer solar zenith angle at sunset
             , CLOUDTOP   = 6        ! k-value above which clear-sky dj assumed
                                     ! (since..... Joffen?)
   integer, parameter, private :: KMAX20=20
   integer, parameter, private :: &
               NPHODIS = 17       &  ! Max possible NRCPHOT
              ,NLAT    = 6           ! No. latitude outputs

    real, allocatable,save, private, dimension(:,:,:,:) :: dj

    real, allocatable,save, private, dimension(:,:,:) :: &
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
    integer, public, parameter ::  IDRCOCHO  = IDRCOHCO ! Just tmp

    integer, save :: photo_out_ix = -1

 !/ subroutines: 

  public :: readdiss
  public :: setup_phot

 contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine readdiss(newseason)

      integer ::  newseason

      integer ::  kz     &  ! help index
                 ,izn   &  ! integer zenith angle
                 ,nr    &  ! numbering of photolytic reactions
                 ,la       ! counting every 10 deg. latitude
        real myz
        character(len=TXTLEN_FILE) fname1
        character(len=3) ::season3

        logical,save:: first_call=.true.

        if(newseason==1)season3='jan'
        if(newseason==2)season3='apr'
        if(newseason==3)season3='jul'
        if(newseason==4)season3='oct' 

        if(first_call)then
           if(.not.(allocated(rcphot)))allocate(rcphot(NRCPHOT,KCHEMTOP:KMAX_MID))
           allocate(dj(NPHODIS,NzPHODIS,HORIZON,NLAT))
           allocate(djcl1(NPHODIS,NzPHODIS,HORIZON))
           allocate(djcl3(NPHODIS,NzPHODIS,HORIZON))

           photo_out_ix = 0
           if(allocated(f_3d)) &
             photo_out_ix = find_index("D3_J(NO2)", f_3d(:)%subclass)
           if(photo_out_ix>0 .and. me==0)write(*,*)'will output J(NO2)'
           
        end if
!    Open, read and broadcast clear sky rates
!---------------

        if(me == 0)then
           fname1 = key2str(jclearFile,'SEASON',season3)
           call open_file(IO_DJ,"r",fname1,needed=.true.)
           call CheckStop(ios,"DefPhotolysis: ios error in jclear ")
        end if


!       Format of input data from Phodis - careful with "17" and NPHODIS
999     FORMAT(1x, f8.3, 17(1x, 1pe8.2)) !Format imposed by file


        if(me == 0)then

          do la = 1,NLAT
            do izn = 1,HORIZON        
              do kz = 1,NzPHODIS
                 !we assume that zPHODIS(kz) is constant, and simply overwrite
                read(IO_DJ,999) zPHODIS(kz),(dj(nr,kz,izn,la),nr=1,NPHODIS)
              end do   ! kz
            end do    ! izn
          end do     ! la
          close(IO_DJ)
        end if  ! me = 0

        CALL MPI_BCAST(zPHODIS  ,8*NzPHODIS,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
        CALL MPI_BCAST(dj  ,8*NPHODIS*NzPHODIS*HORIZON*NLAT,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
!        write(*,*)'zPHODIS, P_PHODIS ',(kz,zPHODIS(kz),P_PHODIS(kz),kz=1,NzPHODIS)

        if(first_call)then
           !convert from z in meters to pressure
           do kz = 1,NzPHODIS             
              P_PHODIS(kz) = 1000*StandardAtmos_km_2_kPa(zPHODIS(kz)) ! P_PHODIS in Pa
              !write(*,*)'PHODIS LEVELS ',kz,P_PHODIS(kz),zPHODIS(kz)
           enddo
        endif

        !   Open, read and broadcast light cloud rates
!---------------

        if(me == 0)then
           fname1 = key2str(jcl1kmFile,'SEASON',season3)
           call open_file(IO_DJ,"r",fname1,needed=.true.)
           call CheckStop(ios,"DefPhotolysis: ios error in jcl1km ")
        end if


        if(me == 0)then

          do izn = 1,HORIZON
             do kz = 1,NzPHODIS
              read(IO_DJ,999) myz,(djcl1(nr,kz,izn),nr=1,NPHODIS)
            end do
          end do  ! izn

          do izn = 1,HORIZON
            do kz = 1,NzPHODIS
              do nr=1,NPHODIS
                djcl1(nr,kz,izn)=djcl1(nr,kz,izn)/dj(nr,kz,izn,3)-1.0
              end do ! nr
            end do ! k
          end do  ! izn
          close(IO_DJ)
        end if   ! me = 0

        CALL MPI_BCAST(djcl1  ,8*NPHODIS*NzPHODIS*HORIZON,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 


!    Open, read and broadcast dense cloud rates
!---------------

        if(me == 0)then
           fname1 = key2str(jcl3kmFile,'SEASON',season3)
           call open_file(IO_DJ,"r",fname1,needed=.true.)
           call CheckStop(ios,"DefPhotolysis: ios error in jcl3km ")
        end if


        if(me == 0)then

           do izn = 1,HORIZON
              do kz = 1,NzPHODIS
                 read(IO_DJ,999) myz,(djcl3(nr,kz,izn),nr=1,NPHODIS)
              end do  ! kz
           end do   ! izn
           close(IO_DJ)
           do izn = 1,HORIZON
              do kz = 1,NzPHODIS
                 do nr=1,NPHODIS
                    djcl3(nr,kz,izn)=djcl3(nr,kz,izn)/dj(nr,kz,izn,3)-1.
                 end do  ! nr
              end do  ! k
           end do   ! izn
        end if      !  me = 0
        
        CALL MPI_BCAST(djcl3  ,8*NPHODIS*NzPHODIS*HORIZON,MPI_BYTE,0,MPI_COMM_CALC,IERROR) 
        
!        do k=1,KMAX_MID
!           write(6,*) 'jverdi i niv. k',k
!           write(6,*) (dj(3,1,k,nr),nr=1,4)
!           write(6,*) (djcl1(1,k,nr),nr=1,4)
!           write(6,*) (djcl3(1,k,nr),nr=1,4)
!         end do
!       end if
        first_call=.false.
        return

        end subroutine readdiss
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        subroutine setup_phot(i,j,errcode)

!       input
        integer :: i,j,iz
        integer :: errcode

!       local
        integer la     &  ! counting every 10 deg. latitude
               ,n      &  ! help index
               ,k      &  ! vertical index
               ,base   &  ! cloud base
               ,top    &  ! cloud top
               ,iclcat    ! cloud type

        real clear        ! clear sky fraction
        integer k2P_PHODIS(KMAX_MID)! converts model level k into PHODIS height index
        logical, save ::firstc=.true.

! make conversion between P_PHODIS heights and model level (k) for this (i,j)
! should be fast and does not need to be too accurate (smooth field) -> no interpolation
        do k=KMAX_MID,1,-1
           k2P_PHODIS(k) = NzPHODIS    
           do iz = NzPHODIS,1,-1
              if(A_bnd(k)+B_bnd(k)*ps(i,j,1)>P_PHODIS(iz))exit
              k2P_PHODIS(k) = iz       
           enddo
!           if(firstc and me==0)write(*,*)'PHODIS level conversion ',k,k2P_PHODIS(k),A_bnd(k)+B_bnd(k)*ps(i,j,1),P_PHODIS(k2P_PHODIS(k))
        enddo
        
        firstc=.false.

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
                 end if
                 iclcat = 1

                 if(z_bnd(i,j,top)-z_bnd(i,j,base) > 1.5e3) iclcat = 2

              end if  ! base<CLOUDTOP
            end if   ! end cc3dmax




            la = max(1,int(0.1*abs(glat(i,j))-2.0001))

            if(iclcat == 0)then
              do k = KCHEMTOP,KMAX_MID
                do n=1,NRCPHOT
                  rcphot(n,k)  = dj(n,k2P_PHODIS(k),Grid%izen,la)
                end do
              end do
            else if(iclcat == 1)then
              clear = cc3dmax(i,j,KMAX_MID)
              do k = KCHEMTOP,KMAX_MID
                do n=1,NRCPHOT
                  rcphot(n,k)  = (1. + clear*djcl1(n,k2P_PHODIS(k),Grid%izen)) &
                                * dj(n,k2P_PHODIS(k),Grid%izen,la)
                end do  !  n
              end do   !  k


            else
              clear = cc3dmax(i,j,KMAX_MID)
              do k = KCHEMTOP,KMAX_MID
                do n=1,NRCPHOT
                  rcphot(n,k)  = (1. +     clear*djcl3(n,k2P_PHODIS(k),Grid%izen)) &
                               * dj(n,k2P_PHODIS(k),Grid%izen,la)
                 end do
              end do
            end if
  
            if ( DEBUG_DJ ) then
                sum_rcphot = sum_rcphot + &
                      sum ( rcphot(1:NRCPHOT, KCHEMTOP:KMAX_MID) )
            end if


          end if   !  end izen <  90 (daytime)  test

          if(photo_out_ix>0) d_3d(photo_out_ix,i,j,1:num_lev3d,IOU_INST) = &
               rcphot(IDNO2,lev3d(1:num_lev3d))

    end subroutine setup_phot
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 end module DefPhotolysis_ml
