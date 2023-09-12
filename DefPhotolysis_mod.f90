! <DefPhotolysis_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version v5.0>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2023 met.no
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
! <DefPhotolysis_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
!+ Photolysis coefficients
!-------------------------------------------------------------------------------

  module DefPhotolysis_mod
    !-------------------------------------------------------------------------------
    
    !    Data needed  for photolysis calculation.  NPHODIS is the number of
    !    tabulated rates from the Phodis model. NRCPHOT (<=NPHODIS) is the
    !    number of photolysis rats needed by the model
    !
    !   10/10/01 - corrected and tidied up by jej.
    !   11/10/01 - NDISS removed, minor F90 changes and docs
    !              added. NLAT and CLOUDTOP added.
    !-------------------------------------------------------------------------------
    
       use CheckStop_mod,      only: CheckStop
       use Config_module,     only: KMAX_MID, KCHEMTOP, NPROC, IOU_INST,&
                                       MasterProc, uses, &
                                       jcl1kmFile,jcl3kmFile,jclearFile,num_lev3d,lev3d
       use DerivedFields_mod,  only: d_3d, f_3d
       use Functions_mod,      only: StandardAtmos_km_2_kPa
       use GridValues_mod   , only : glat, A_bnd, B_bnd
       use Io_mod,            only : IO_DJ, open_file, ios
       use LocalVariables_mod,only : Grid  ! => izen
       use MetFields_mod    , only : cc3d,cc3dmax,z_bnd,ps
       use MPI_Groups_mod   , only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_INTEGER&
                                         ,MPI_COMM_CALC, IERROR
       use OwnDataTypes_mod ,  only: TXTLEN_FILE
       use Par_mod      ,     only : me,LIMAX,LJMAX
       use SmallUtils_mod,     only: key2str, find_index
       use ZchemData_mod,      only: rcphot
       implicit none
       private

       integer, public, parameter :: &
       NRCPHOT          = 17  & ! Number of photolytic reactions
      ,NRCPHOTextended  = 172   ! currently max number, used for CRIv2R5

       integer, public, parameter :: & ! cloud-j index declarations
       CLJ1 = 1,CLJ2 = 2,CLJ3 = 3,CLJ4 = 4,CLJ5 = 5,CLJ6 = 6,CLJ7 = 7,CLJ8 = 8, CLJ9 = 9, & 
       CLJ10=10,CLJ11=11,CLJ12=12,CLJ13=13,CLJ14=14,CLJ15=15,CLJ16=16,CLJ17=17, CLJ18=18, &
       CLJ19=19,CLJ20=20,CLJ21=21,CLJ22=22,CLJ23=23,CLJ24=24,CLJ25=25,CLJ26=26, CLJ27=27, &
       CLJ28=28,CLJ29=29,CLJ30=30,CLJ31=31,CLJ32=32,CLJ33=33,CLJ34=34,CLJ35=35, CLJ36=36, &
       CLJ37=37,CLJ38=38,CLJ39=39,CLJ40=40,CLJ41=41,CLJ42=42,CLJ43=43,CLJ44=44, CLJ45=45, &
       CLJ46=46,CLJ47=47,CLJ48=48,CLJ49=49,CLJ50=50,CLJ51=51,CLJ52=52,CLJ53=53, CLJ54=54, &
       CLJ55=55,CLJ56=56,CLJ57=57,CLJ58=58,CLJ59=59,CLJ60=60,CLJ61=61,CLJ62=62, CLJ63=63, &
       CLJ64=64,CLJ65=65,CLJ66=66,CLJ67=67,CLJ68=68,CLJ69=69,CLJ70=70,CLJ71=71, CLJ72=72, &
       CLJ73=73,CLJ74=74,CLJ75=75,CLJ76=76,CLJ77=77,CLJ78=78,CLJ79=79,CLJ80=80, CLJ81=81, &
       CLJ82=82,CLJ83=83,CLJ84=84,CLJ85=85,CLJ86=86,CLJ87=87,CLJ88=88,CLJ89=89, CLJ90=90, &
       CLJ91=91,CLJ92=92,CLJ93=93,CLJ94=94,CLJ95=95,CLJ96=96,CLJ97=97,CLJ98=98, CLJ99=99, &
       CLJ100=100,CLJ101=101,CLJ102=102,CLJ103=103,CLJ104=104,CLJ105=105,CLJ106=106, &
       CLJ107=107,CLJ108=108,CLJ109=109,CLJ110=110,CLJ111=111,CLJ112=112,CLJ113=113, &
       CLJ114=114,CLJ115=115,CLJ116=116,CLJ117=117,CLJ118=118,CLJ119=119,CLJ120=120, &
       CLJ121=121,CLJ122=122,CLJ123=123,CLJ124=124,CLJ125=125,CLJ126=126,CLJ127=127, &
       CLJ128=128,CLJ129=129,CLJ130=130,CLJ131=131,CLJ132=132,CLJ133=133,CLJ134=134, &
       CLJ135=135,CLJ136=136,CLJ137=137,CLJ138=138,CLJ139=139,CLJ140=140,CLJ141=141, &
       CLJ142=142,CLJ143=143,CLJ144=144,CLJ145=145,CLJ146=146,CLJ147=147,CLJ148=148, &
       CLJ149=149,CLJ150=150,CLJ151=151,CLJ152=152,CLJ153=153,CLJ154=154,CLJ155=155, &
       CLJ156=156,CLJ157=157,CLJ158=158,CLJ159=159,CLJ160=160,CLJ161=161,CLJ162=162, &
       CLJ163=163,CLJ164=164,CLJ165=165,CLJ166=166,CLJ167=167,CLJ168=168,CLJ169=169, &
       CLJ170=170,CLJ171=171,CLJ172=172
    
       ! EmChem19 definitions matching CloudJ indices from CloudJ_EmChem19 input folder/file
      integer, public, parameter :: &
      IDO3_O3P  = CLJ3 , IDO3_O1D = CLJ4 , IDNO2    = CLJ9 , IDHCHO_H = CLJ5 , &
      IDHCHO_H2 = CLJ6 , IDH2O2   = CLJ7 , IDCH3O2H = CLJ8 , IDNO3    = CLJ74, &
      IDHO2NO2  = CLJ16, IDRCOCHO = CLJ64, IDCH3COY = CLJ72, IDMEK    = CLJ75, &
      IDCH3CHO  = CLJ54, IDGLYOXA = CLJ65, IDGLYOXB = CLJ66, IDGLYOXC = CLJ67, &
      IDPAN     = CLJ52, IDHNO3   = CLJ15, IDHONO   = CLJ14, IDCHOCHO = CLJ73, &
      IDNO3_NO  = CLJ11, IDNO3_NO2= CLJ12, IDACETON = CLJ68, IDN2O5   = CLJ13
      
      ! in CRIv2R5 IDO3_O1D = CLJ3, IDNO2 = CLJ11
      ! in EmChem switch CLJ9 and CLJ11, and 3 and 4

      ! copies with different names for situational use
      integer, public, parameter :: &
      IDAO3 = IDO3_O3P, IDBO3 = IDO3_O1D, IDBCH2O = IDHCHO_H2, IDACH2O = IDHCHO_H, &
      IDHCOHCO = IDCHOCHO, IDCH3COX = CLJ75 

       
       integer, public, parameter:: NzPHODIS=20 !number of heights defined in the input files
       real, save, public :: zPHODIS(NzPHODIS) !heights of the input files, in km assumed constants
       real, save, public :: P_PHODIS(NzPHODIS) !heights of the input files, in Pa
    
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
    
        ! Indices of photolysis rates as available from Phodis files:
        integer, public, parameter ::  &
          IDAO3_tb    =  1 , IDBO3_tb    =  2 , IDNO2_tb    =  3 , &
          IDH2O2_tb   =  4 , IDHNO3_tb   =  5 , IDACH2O_tb  =  6 , &
          IDBCH2O_tb  =  7 , IDCH3CHO_tb =  8 , IDCH3COX_tb =  9 , &
          IDCH3COY_tb = 10 , IDHCOHCO_tb = 11 , IDRCOHCO_tb = 12 , &
          IDNO3_tb    = 13 , IDN2O5_tb   = 14 , IDCH3O2H_tb = 15 , &
          IDHO2NO2_tb = 16 , IDACETON_tb = 17
    
    ! Corresponding MCM-photolysis rates:
    ! IDAO3    =  1 = IDO3_O3P = MCM_J2
    ! IDBO3    =  2 = IDO3_O1D = MCM_J1
    ! IDNO2    =  3 = MCM_J4
    ! IDH2O2   =  4 = MCM_J3
    ! IDHNO3   =  5 = MCM_J8
    ! IDACH2O  =  6 = MCM_J11
    ! IDBCH2O  =  7 = MCM_J12
    ! IDCH3CHO =  8 = MCM_J13
    ! IDCH3COX =  9 = MCM_J22
    ! IDCH3COY = 10 = IDBIACET = MCM_J35 ! Based on the assumption that EMEP use the notation CH3COY in the same way as the OSLO CTM3
    ! IDHCOHCO = 11 = ? MCM_J31+MCM_J32+MCM_J33 ?
    ! IDRCOHCO = 12 = IDMGLYOX = MCM_J34
    ! IDNO3    = 13 = ? MCM_J5+MCM_J6 ?
    ! IDN2O5   = 14 = NOT_IN_MCM?
    ! IDCH3O2H = 15 = MCM_J41
    ! IDHO2NO2 = 16 = NOT_IN_MCM?
    ! IDACETON = 17 = MCM_J21
    
    !       ,IDCHOCHO_2CHO  = IDHCOHCO & ! Just name change CHECK, TMP!!!
    !       ,IDCHOCHO_2CO   = IDHCOHCO & ! Just name change CHECK, TMP!!!
    !       ,IDCHOCHO_HCHO  = IDHCOHCO & ! Just name change CHECK, TMP!!!
    !    integer, public, save :: IDCHOCHO
    ! In EMEP code, we use IDBO3 for O1D production and IDAO3 for OP production
    ! MCM indices:
    
    !NEEDS FIXING. Changed from ESX to try to match above, but eg NO3 is difficult
    !  integer, public, parameter :: & 
    !   IDO3_O1D   = 2,IDO3_O3P  = 1, & !:BUG FIX RB Apr25
    !    IDNO3_NO  = IDNO3  &
    !   ,IDNO3_NO2  = IDNO3 & !HONO NEEDS FIXING!
    !   IDHCHO_H  = 6 & ! HCHO -> CO + 2 HO2
    !   ,IDHCHO_H2 = 7 !&  ! HCHO -> CO + H2
    !   ,MCM_J18   = 18, MCM_J20   = 20 &
    !   ,MCM_J22   = 22 , IDMEK     = 22 &
    !   ,IDCHOCHO_2CO  = 31 &! .. 33 QUERY
    !   ,IDCHOCHO_HCHO = 32 &! .. 33 QUERY
    !   ,IDCHOCHO_2CHO = 33 &! .. 33 QUERY
    !   ,IDCH3COCH3  = IDACETON
    !NEEDS FIXING. ESX has:
    !  integer, public, parameter :: & 
    !    IDO3_O1D   = 1,IDO3_O3P  = 2,IDH2O2    = 3,IDNO2     = 4 ,IDNO3_NO  = 5 &
    !   ,IDNO3_NO2  = 6,IDHONO    = 7,IDHNO3    = 8,IDHCHO_H  = 11,IDHCHO_H2 = 12&
    !   ,IDCH3CHO  = 13 &
    !   ,MCM_J18   = 18, MCM_J20   = 20 &
    !   ,MCM_J22   = 22 , IDMEK     = 22 &
    !   ,IDCHOCHO_2CO  = 31 &! .. 33 QUERY
    !   ,IDCHOCHO_HCHO = 32 &! .. 33 QUERY
    !   ,IDCHOCHO_2CHO = 33 &! .. 33 QUERY
    !   ,IDCH3COCH3  = 34, IDRCOCHO  = 34 &!  MGLOX, etc
    !   ,IDCH3O2H  = 41 &
    !   ,IDiC3H7ONO2  = 54 &! Used for ISON
    !   ,IDCH3COO2H  = -1 &! QUERY ???
    !   ,IDN2O5    = -1 ! .. QUERY in PHUX; not MCM?
    
        integer, save, public :: photo_out_ix_no2 = -1
        integer, save, public :: photo_out_ix_o3a = -1
        integer, save, public :: photo_out_ix_o3b = -1
        integer, save, public :: photo_out_ix_h2o2 = -1
        integer, save, public :: photo_out_ix_hno3 = -1
        integer, save, public :: photo_out_ix_ach2o = -1
        integer, save, public :: photo_out_ix_bch2o = -1
        integer, save, public :: photo_out_ix_hono = -1
        integer, save, public :: photo_out_ix_ho2no2 = -1
        integer, save, public :: photo_out_ix_no3 = -1
        integer, save, public :: photo_out_ix_ch3o2h = -1 
        integer, save, public :: photo_out_ix_MEK = -1 
        integer, save, public :: photo_out_ix_N2O5 = -1 
        integer, save, public :: photo_out_ix_GLYOX = -1 
        integer, save, public :: photo_out_ix_CH3CHO = -1 
        integer, save, public :: photo_out_ix_ACETON = -1 
        integer, save, public :: photo_out_ix_MGLYOX = -1 
        integer, save, public :: photo_out_ix_BIACET = -1 
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
              allocate(dj(NPHODIS,NzPHODIS,HORIZON,NLAT))
              allocate(djcl1(NPHODIS,NzPHODIS,HORIZON))
              allocate(djcl3(NPHODIS,NzPHODIS,HORIZON))
    
              if(allocated(f_3d).and.(.not.uses%cloudj)) then
                photo_out_ix_no2 = find_index("D3_J(NO2)", f_3d(:)%subclass)
                photo_out_ix_o3a = find_index("D3_J(O3a)", f_3d(:)%subclass)
                photo_out_ix_o3b = find_index("D3_J(O3b)", f_3d(:)%subclass)
                photo_out_ix_h2o2 = find_index("D3_J(H2O2)", f_3d(:)%subclass)
                photo_out_ix_hno3 = find_index("D3_J(HNO3)", f_3d(:)%subclass)
                photo_out_ix_ach2o = find_index("D3_J(ACH2O)", f_3d(:)%subclass)
                photo_out_ix_bch2o = find_index("D3_J(BCH2O)", f_3d(:)%subclass)
                photo_out_ix_hono = find_index("D3_J(HONO)", f_3d(:)%subclass)
                photo_out_ix_ho2no2 = find_index("D3_J(HO2NO2)", f_3d(:)%subclass)
                photo_out_ix_no3 = find_index("D3_J(NO3)", f_3d(:)%subclass)
                photo_out_ix_ch3o2h = find_index("D3_J(CH3O2H)", f_3d(:)%subclass)
                photo_out_ix_MEK = find_index("D3_J(MEK)", f_3d(:)%subclass)
                photo_out_ix_N2O5 = find_index("D3_J(N2O5)", f_3d(:)%subclass)
                photo_out_ix_GLYOX = find_index("D3_J(GLYOX)", f_3d(:)%subclass)
                photo_out_ix_CH3CHO = find_index("D3_J(CH3CHO)", f_3d(:)%subclass)
                photo_out_ix_ACETON = find_index("D3_J(ACETON)", f_3d(:)%subclass)
                photo_out_ix_MGLYOX = find_index("D3_J(MGLYOX)", f_3d(:)%subclass)
                photo_out_ix_BIACET = find_index("D3_J(BIACET)", f_3d(:)%subclass)
                if(me==0)write(*,*) 'Outputting tabulated J-values specified in config.'
              end if              
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
                  do while(cc3d(i,j,k,1) < 1.e-4 .and. k >= CLOUDTOP)
                     k = k-1
                  end do
                  base = k+1
    
                 ! if all cc3d so far are <1.e-4 we are done
                  if(USES%CLEARSKYTAB) then
                  
                    if( base < CLOUDTOP ) then 
    
                    !  we have found a k>=CLOUDTOP with cc3d>=1.e-4, now search for top
    
                    k = CLOUDTOP
                    do while(cc3d(i,j,k,1) < 1.0e-4)
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

                else
                  
                  if( base > CLOUDTOP ) then
    
                    !  we have found a k>=CLOUDTOP with cc3d>=1.e-4, now search for top
    
                    k = CLOUDTOP
                    do while(cc3d(i,j,k,1) < 1.0e-4)
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
    
                  end if ! base > cloudtop 
                end if ! uses%clearskytab
              end if ! end cc3dmax
    
    
    
    
                la = max(1,int(0.1*abs(glat(i,j))-2.0001))
    
                if(iclcat == 0)then
                  do k = KCHEMTOP,KMAX_MID
                      ! match tabulated rates with indices in EMEP
                      rcphot(IDO3_O3P,k)  = dj(IDAO3_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDO3_O1D,k)  = dj(IDBO3_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDNO2,k)     = dj(IDNO2_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDH2O2,k)    = dj(IDH2O2_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDHNO3,k)    = dj(IDHNO3_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDHCHO_H,k)  = dj(IDACH2O_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDHCHO_H2,k) = dj(IDBCH2O_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDCH3CHO,k)  = dj(IDCH3CHO_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDMEK,k)     = dj(IDCH3COX_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDCH3COY,k)  = dj(IDCH3COY_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDCHOCHO,k)  = dj(IDHCOHCO_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDRCOCHO,k)  = dj(IDRCOHCO_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDNO3,k)     = dj(IDNO3_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDN2O5,k)    = dj(IDN2O5_tb,k2P_PHODIS(k),Grid%izen,la)  
                      rcphot(IDCH3O2H,k)  = dj(IDCH3O2H_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDHO2NO2,k)  = dj(IDHO2NO2_tb,k2P_PHODIS(k),Grid%izen,la)
                      rcphot(IDACETON,k)  = dj(IDACETON_tb,k2P_PHODIS(k),Grid%izen,la)

                  end do
                else if(iclcat == 1)then
                  clear = cc3dmax(i,j,KMAX_MID)
                  do k = KCHEMTOP,KMAX_MID
                      ! match tabulated rates with indices in EMEP
                      rcphot(IDO3_O3P,k)  = (1. + clear*djcl1(IDAO3_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDAO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDO3_O1D,k)  = (1. + clear*djcl1(IDBO3_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDBO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDNO2,k)     = (1. + clear*djcl1(IDNO2_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDNO2_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDH2O2,k)    = (1. + clear*djcl1(IDH2O2_tb,k2P_PHODIS(k),Grid%izen))   &
                                          * dj(IDH2O2_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDHNO3,k)    = (1. + clear*djcl1(IDHNO3_tb,k2P_PHODIS(k),Grid%izen))   &
                                          * dj(IDHNO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDHCHO_H,k)  = (1. + clear*djcl1(IDACH2O_tb,k2P_PHODIS(k),Grid%izen))  &
                                          * dj(IDACH2O_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDHCHO_H2,k) = (1. + clear*djcl1(IDBCH2O_tb,k2P_PHODIS(k),Grid%izen))  &
                                          * dj(IDBCH2O_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCH3CHO,k)  = (1. + clear*djcl1(IDCH3CHO_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3CHO_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDMEK,k)     = (1. + clear*djcl1(IDCH3COX_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3COX_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCH3COY,k)  = (1. + clear*djcl1(IDCH3COY_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3COY_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCHOCHO,k)  = (1. + clear*djcl1(IDHCOHCO_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDHCOHCO_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDRCOCHO,k)  = (1. + clear*djcl1(IDRCOHCO_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDRCOHCO_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDNO3,k)     = (1. + clear*djcl1(IDNO3_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDNO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDN2O5,k)    = (1. + clear*djcl1(IDN2O5_tb,k2P_PHODIS(k),Grid%izen))   &
                                          * dj(IDN2O5_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCH3O2H,k)  = (1. + clear*djcl1(IDCH3O2H_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3O2H_tb,k2P_PHODIS(k),Grid%izen,la)
                                          
                      rcphot(IDHO2NO2,k)  = (1. + clear*djcl1(IDHO2NO2_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDHO2NO2_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDACETON,k)  = (1. + clear*djcl1(IDACETON_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDACETON_tb,k2P_PHODIS(k),Grid%izen,la)
                    
                  end do   !  k
    
    
                else
                  clear = cc3dmax(i,j,KMAX_MID)
                  do k = KCHEMTOP,KMAX_MID
                      ! match tabulated rates with indices in EMEP
                      rcphot(IDO3_O3P,k)  = (1. + clear*djcl3(IDAO3_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDAO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDO3_O1D,k)  = (1. + clear*djcl3(IDBO3_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDBO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDNO2,k)     = (1. + clear*djcl3(IDNO2_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDNO2_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDH2O2,k)    = (1. + clear*djcl3(IDH2O2_tb,k2P_PHODIS(k),Grid%izen))   &
                                          * dj(IDH2O2_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDHNO3,k)    = (1. + clear*djcl3(IDHNO3_tb,k2P_PHODIS(k),Grid%izen))   &
                                          * dj(IDHNO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDHCHO_H,k)  = (1. + clear*djcl3(IDACH2O_tb,k2P_PHODIS(k),Grid%izen))  &
                                          * dj(IDACH2O_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDHCHO_H2,k) = (1. + clear*djcl3(IDBCH2O_tb,k2P_PHODIS(k),Grid%izen))  &
                                          * dj(IDBCH2O_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCH3CHO,k)  = (1. + clear*djcl3(IDCH3CHO_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3CHO_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDMEK,k)     = (1. + clear*djcl3(IDCH3COX_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3COX_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCH3COY,k)  = (1. + clear*djcl3(IDCH3COY_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3COY_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCHOCHO,k)  = (1. + clear*djcl3(IDHCOHCO_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDHCOHCO_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDRCOCHO,k)  = (1. + clear*djcl3(IDRCOHCO_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDRCOHCO_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDNO3,k)     = (1. + clear*djcl3(IDNO3_tb,k2P_PHODIS(k),Grid%izen))    &
                                          * dj(IDNO3_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDN2O5,k)    = (1. + clear*djcl3(IDN2O5_tb,k2P_PHODIS(k),Grid%izen))   &
                                          * dj(IDN2O5_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDCH3O2H,k)  = (1. + clear*djcl3(IDCH3O2H_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDCH3O2H_tb,k2P_PHODIS(k),Grid%izen,la)
                                          
                      rcphot(IDHO2NO2,k)  = (1. + clear*djcl3(IDHO2NO2_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDHO2NO2_tb,k2P_PHODIS(k),Grid%izen,la)

                      rcphot(IDACETON,k)  = (1. + clear*djcl3(IDACETON_tb,k2P_PHODIS(k),Grid%izen)) &
                                          * dj(IDACETON_tb,k2P_PHODIS(k),Grid%izen,la)

                  end do
                end if
      
                if ( DEBUG_DJ ) then
                    sum_rcphot = sum_rcphot + &
                          sum ( rcphot(1:NRCPHOT, KCHEMTOP:KMAX_MID) )
                end if
    
               ! adding HONO
                rcphot(IDHONO,:)  =  0.22* rcphot(IDNO2,:)
                rcphot(IDNO3_NO,:) = 0.127 * rcphot(IDNO3,:)
                rcphot(IDNO3_NO2,:) = 0.873 * rcphot(IDNO3,:)
    
              end if   !  end izen <  90 (daytime)  test
    
        end subroutine setup_phot
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
     end module DefPhotolysis_mod