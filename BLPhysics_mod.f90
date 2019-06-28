! <BLPhysics_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2019 met.no
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
module BLPhysics_mod

 ! Collection of boundary layer met routines. Will move more code from tiphys
 !  here in future. Try to keep 1-D or elemental to allow use in offline codes 
 ! (*No* routines in use, except for testing)

 use Config_module,    only : KMAX_MID, KMAX_BND, KWINDTOP, PT, PBL
 use Landuse_mod,           only : Landcover, water_fraction
 use PhysicalConstants_mod, only : KARMAN, GRAV
 implicit none
 private

! minimum value now generally calculated as z_mid(19), but we
!   keep a fixed value for smoothing. 
! real, parameter, public :: PBL_ZiMIN=100.   ! EMEP/TI and smooth(zi)
! real, parameter, public :: PBL_ZiMAX=3000.  ! EMEP/TI

!NWPHmix now from Config_emep, NWP%HmixMethod
!NWPHMIX! Choose one Hmix method here (not needed for NWP?)
!NWPHMIX character(len=4), parameter, public :: HmixMethod = &
!NWPHMIX     "JcRb"   ! Jericevic
!NWPHMIX    !"SbRb"   ! Seibert
!NWPHMIX    !"TIZi"   ! Original from Trond Iversen tiphysics

! Choose one Kz method here. Prefered method is likely to use O'brien
! in convective, Jericevic in Stable.
 logical, parameter, public  :: NWP_Kz=.false. ! hb 23.02.2010 Kz from meteo 
 logical, parameter, public  :: USE_MIN_KZ =.false. ! "fix"
  character(len=2), parameter, public :: KzMethod = &
     "--"   ! Set U, S separately, preferred? :
  !   "JG"   ! Jericevic/Grisogono - both unstable + stable
  character(len=2), parameter, public :: UnstableKzMethod = &
     "OB"   ! O'Brien
  character(len=2), parameter, public :: StableKzMethod = &
     "JG"   ! Jericevic/Grisogono
    !"BW"   ! Brost Wynngard
    !"Sb"   ! Seibert

!Movbed
!  character(len=4), parameter, public :: FluxPROFILE = &
!     "Iter"   ! 
! !     "Ln95"   ! ! will use Launiainen1995 

 logical, parameter, public :: PIELKE = .true.
 real, public, parameter :: KZ_MINIMUM = 0.001 ! m2/s
 real, public, parameter :: KZ_MAXIMUM = 1.0e3 ! m2/s - as old kzmax
 real, public, parameter :: KZ_SBL_LIMIT = 0.1 ! m2/s - Defines stable BL height
 ! TI code had unstable if delq > 0.00001, ca. fh >10-8 so excluded neutral
 real, public, parameter :: OB_invL_LIMIT =  -1.0e-10
 real, public, parameter :: MIN_USTAR_LAND = 0.1 ! ms - Defines stable BL height

 real, parameter, private :: EPS=0.01  !prevents div by zero for WS

! - Hmix routines
public :: SeibertRiB_Hmix
public :: SeibertRiB_Hmix_3d
public :: JericevicRiB_Hmix
public :: JericevicRiB_Hmix0  ! Now allow mixing heights based upon surface T
public :: Venkatram_Hmix
public :: VogelezangHoltslag_Hmix
public :: Zilitinkevich_Hmix
public :: TI_Hmix

! Kz routines
public :: PielkeBlackadarKz
public :: BrostWyngaardKz
public :: JericevicKz
public :: O_BrienKz

! Misc
public :: fake_zmid
public :: fake_zbnd
public :: risig1

! Test all
public :: Test_BLM

! Conversions
public :: SigmaKz_2_m2s ! hb 23.02.2010 Kz from meteo
 private :: SigmaKz_2_m2s_scalar  ! function to get factor
 private :: SigmaKz_2_m2s_arrays  ! subrouitne for 3d arrays
public :: Kz_m2s_toSigmaKz
public :: Kz_m2s_toEtaKz

! Conversion of Kz in sigma coordinates to m2/s,
!  Kz(sigma)=Kz*ro**2*(GRAV/p*)**2
! We can call this conversion routine as either scalar or array
  interface SigmaKz_2_m2s
     module procedure SigmaKz_2_m2s_scalar
     module procedure SigmaKz_2_m2s_arrays
  end interface SigmaKz_2_m2s

contains

 !----------------------------------------------------------------------------

function BrostWyngaardKz(z,h,ustar,invL,Kdef) result(Kz)
  real, intent(in) :: z    ! height
  real, intent(in) :: h    ! Boundary layer depth 
  real, intent(in) :: ustar!  u*
  real, intent(in) :: invL !  1/L  
  real, intent(in) :: Kdef !  1/L  
  real :: Kz

     if ( z < h ) then
        Kz = KARMAN * ustar * z * (1-z/h)**1.5 / (1+5*z*invL)
     else
       Kz =  Kdef
     end if

end function BrostWyngaardKz

function JericevicKz(z,h,ustar,Kdef) result(Kz)
  real, intent(in) :: z    ! height
  real, intent(in) :: h    ! Boundary layer depth 
  real, intent(in) :: ustar, Kdef !  u*, default Kz
  real :: Kz
  real :: Kmax, zmax


     if ( z < h ) then
        Kmax = 0.05 * h * ustar
        zmax = 0.21 * h
        Kz = 0.39 * ustar * z * exp( -0.5*(z/zmax)**2 )
     !OS_TEST_Hmix1 

     else ! open-source had this Kz=0.0 line. Not sure why
     !OS_TEST_Hmix1    
     !Kz= 0.0
       Kz =  Kdef
     end if

end function JericevicKz


!----------------------------------------------------------------------------
! Two Rib-based mixing height methods.
! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
! Jericevic et al., ACP, 2009,  eqn (17): 
! Boundary layer height is calculated with the bulk Richardson number
! method with critical value of Ric=0.25
!----------------------------------------------------------------------------

subroutine SeibertRiB_Hmix_3d (u,v, zm, theta, pzpbl)
  real, dimension(:,:,:), intent(in) :: u,v ! winds
  real, dimension(:,:,:), intent(in) :: zm ! mid-cell height
  real, dimension(:,:,:), intent(in) :: theta !pot. temp
  real, intent(out) :: pzpbl(:,:)
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real,dimension(size(pzpbl,1),size(pzpbl,2)) :: RiB, Theta1 ! pot temp of lowest cell

! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
! Although we should use virtual pot temp, not just theta

     Theta1(:,:) = theta(:,:,KMAX_MID)
     pzpbl = -999.999
     KLOOP: do k=KMAX_MID-1, KWINDTOP, -1 

        where ( pzpbl < 0 ) 
           Rib(:,:) =      GRAV * zm(:,:,k) &
             *(theta(:,:,k)-Theta1(:,:) ) / &
             ( Theta1(:,:) * ( u(:,:,k)**2 + v(:,:,k)**2 )+EPS )
           where (Rib(:,:) >= Ric) 
              pzpbl(:,:) = zm(:,:,k)
           end where 
         end where 

      end do KLOOP

end subroutine SeibertRiB_Hmix_3d

 !----------------------------------------------------------------------------
subroutine SeibertRiB_Hmix (u,v, zm, theta, pzpbl)
  real, dimension(KMAX_MID), intent(in) :: u,v ! winds
  real, dimension(KMAX_MID), intent(in) :: zm ! mid-cell height
  real, dimension(KMAX_MID), intent(in) :: theta !pot. temp
  real, intent(out) :: pzpbl
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real :: Rib  ! bulk Richardson number
  real :: Theta1   ! pot temp of lowest cell

! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
! Although we should use virtual pot temp, not just theta

   Theta1 = theta(KMAX_MID)
   do k=KMAX_MID-1, KWINDTOP, -1
      Rib =      GRAV * zm(k) &
             *(theta(k)-Theta1 ) / &
             ( Theta1 * ( u(k)**2 + v(k)**2 )+EPS )
             !print *, k, zm(k), theta(k), sqrt(( u(k)**2 + v(k)**2 )),  RiB
      if(Rib >= Ric) then
             pzpbl = zm(k)
             exit
      end if
   end do

end subroutine SeibertRiB_Hmix

 !----------------------------------------------------------------------------

subroutine JericevicRiB_Hmix (u,v, zm, theta, zi)
  real, dimension(KMAX_MID), intent(in) :: u,v ! winds
  real, dimension(KMAX_MID), intent(in) :: zm ! mid-cell height
  real, dimension(KMAX_MID), intent(in) :: theta !pot. temp
  real, intent(out) :: zi
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real :: Rib  ! bulk Richardson number
  real :: Theta1, z1  ! pot temp  and height of lowest cell

! Jericevic et al., ACP, 2009, pp1001-,  eqn (17): 

   Theta1 = theta(KMAX_MID)
   z1     = zm(KMAX_MID)
   zi     = z1  ! start val

   do k=KMAX_MID-1, KWINDTOP, -1

       Rib =   GRAV * ( zm(k) - z1 ) &
             * (theta(k)-Theta1 ) / &
       ( 0.5*(theta(k)+Theta1) * ( u(k)**2 + v(k)**2 )+EPS )
       if(Rib >= Ric) then
              zi = zm(k)
              exit
       end if
    end do

end subroutine JericevicRiB_Hmix

 !----------------------------------------------------------------------------
subroutine JericevicRiB_Hmix0 (u,v, zm, theta, zi)
 !- as above, but allow test for surface SBL
  real, dimension(KMAX_MID), intent(in) :: u,v ! winds
  real, dimension(KMAX_MID), intent(in) :: zm ! mid-cell height
  real, dimension(KMAX_MID), intent(in) :: theta !pot. temp
  real, intent(out) :: zi
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real :: Rib  ! bulk Richardson number
  real :: Theta1, z1  ! pot temp  and height of lowest cell

! Jericevic et al., ACP, 2009, pp1001-,  eqn (17):

   Theta1 = theta(KMAX_MID)
   z1     = zm(KMAX_MID)
   zi     = z1  ! start val

   do k=KMAX_MID-1, KWINDTOP, -1

       Rib =   GRAV * ( zm(k) - z1 ) &
             * (theta(k)-Theta1 ) / &
       ( 0.5*(theta(k)+Theta1) * ( u(k)**2 + v(k)**2 )+EPS )
       if(Rib >= Ric) then
              zi = zm(k)
              exit
       end if
    end do

end subroutine JericevicRiB_Hmix0

subroutine VogelezangHoltslag_Hmix (u,v, zm, theta, q, ustar, pzpbl)
  real, dimension(KMAX_MID), intent(in) :: u,v ! winds
  real, dimension(KMAX_MID), intent(in) :: zm ! mid-cell height
  real, dimension(KMAX_MID), intent(in) :: theta !pot. temp
  real, dimension(KMAX_MID), intent(in) :: q ! spec humid
  real, intent(in) :: ustar
  real, intent(out) :: pzpbl
  real, dimension(KMAX_MID):: tv !virtual pot. temp
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real :: Rig  ! bulk Richardson number
  real :: Theta1, u1, v1, z1, bu2 !  values for lowest grid cell

! Vogelezang, D. & Holtslag, A.,  BLM, 1996, 81, 245-269, Eqn. (3)
! Although we should use virtual pot temp, not just theta


  !     Calculate virtual temperature

   tv(:) = theta(:) * (1.0+0.622*q(:))

   Theta1 = tv(KMAX_MID)
   z1     = zm(KMAX_MID)
   bu2    = 100.0 * ustar * ustar
   u1     = u(KMAX_MID)
   v1     = v(KMAX_MID)

   do k=KMAX_MID-1, KWINDTOP, -1
      Rig =      GRAV * ( zm(k) - z1 )  &
             *(tv(k)-Theta1 ) / &
             ( Theta1 * ( ( u(k) - u1) **2 + ( v(k) - v1) **2 )+ bu2 + EPS )
             !print *, k, zm(k), theta(k), sqrt(( u(k)**2 + v(k)**2 )),  RiB
      if(Rig >= Ric) then
             pzpbl = zm(k)
             exit
      end if
   end do

end subroutine VogelezangHoltslag_Hmix
 !----------------------------------------------------------------------------
function Venkatram_Hmix (ustar) result(zi)
  real, intent(in) :: ustar
  real :: zi
  ! From Venkatram, 1980, simple method!

     zi = 2.4e3 * ustar**1.5

end function Venkatram_Hmix
 !----------------------------------------------------------------------------

function Zilitinkevich_Hmix (ustar,invL,lat) result(zi)
 ! Not equator-proof yet!!
   real, intent(in) :: ustar, invL, lat
   real :: f, pi = 4.0*atan(1.0)
   real :: zi
   f = 1.46e-4 * sin(lat*pi/180.0)
   !zi = 0.2 *  ustar/f
   if( invL > f/ustar ) then
     zi = 0.4 * sqrt( ustar/(f*invL) )
   else
     zi = 0.2 *  ustar/f
   end if
end function Zilitinkevich_Hmix

 !----------------------------------------------------------------------------

subroutine PielkeBlackadarKz (u,v, zm, zb, th, Kz, Pielke_flag, debug_flag)

!  Use of Pielke/Blackadar for Kz as in orig EMEP and MET.NO NWP model 
!  (Derived from T.Iversen's tiphysics code)
!
  real, dimension(:), intent(in) :: u,v ! wind-vels
  real, dimension(:), intent(in) :: zm, zb ! heights of mid and bnd
  real, dimension(:), intent(in) :: th
  real, dimension(:), intent(out) :: Kz ! Kz  (m2/s)
  logical, intent(in) :: Pielke_flag, debug_flag

  real, dimension(KMAX_BND) :: Ris   ! Richardson number, sigma (???) coords

  integer :: k, km
  real :: xl2    ! mixing length squared
  real, save :: zmmin = 200.0 !QUERY, Pielke or EMEP?
  real :: Ric, Ric0 !critical Richardson number variables
  real :: dvdz

 !..the following variables in sigmas-levels:

   do k=2,KMAX_MID

      km=k-1

      !..wind sheer (first keep as squared)

      dvdz = (u(km)-u(k))**2 + (v(km)-v(k))**2 + EPS 

      Ris(k)=(2*GRAV/(th(km)+th(k))) * (th(km)-th(k))*(zm(km)-zm(k)) &
                  /dvdz

      !........................
      !..mixing length squared:
      !
       xl2=(KARMAN*min(zb(k),zmmin))**2

      !..............................
      !..critical richardsons number:
      !
       Ric0=0.115*((zm(km)-zm(k))*100.0)**0.175
       Ric=max(0.25,Ric0)

       dvdz = sqrt(dvdz)/(zm(km)-zm(k))

      !..................................................................
      !..exchange coefficient (Pielke,...)
      if ( Pielke_flag ) then
         !if( debug_flag ) write(6,*) "BLPielke Ris ",k, Ris(k), Ric
         if (Ris(k) > Ric ) then
            Kz(k) = KZ_MINIMUM
            if( debug_flag ) write(6,"(a,i3,9es10.2)") "BLPielke Kmin ",k, Ris(k), Ric, Kz(k)
         else
            Kz(k) = 1.1 * (Ric-Ris(k)) * xl2 * dvdz /Ric
           if( debug_flag ) write(6,"(a,i3,es10.2,f6.3,9es10.2)") "BLPielke Ks ",k, &
                Ris(k), Ric, xl2, dvdz, Kz(k)
         end if
      else

         !..exchange coefficient (Blackadar, 1979; Iversen & Nordeng, 1987):
         !
         if(Ris(k) <= 0.0) then
            Kz(k)=xl2*dvdz*sqrt(1.1-87.*Ris(k))
         elseif(Ris(k) <= 0.5*Ric) then
                Kz(k)=xl2*dvdz*(1.1-1.2*Ris(k)/Ric)
         elseif(Ris(k) <= Ric) then
                Kz(k)=xl2*dvdz*(1.-Ris(k)/Ric)
         else
                Kz(k)=KZ_MINIMUM
         end if
      end if ! Pielke or Blackadar

   end do ! k
end subroutine PielkeBlackadarKz

 !----------------------------------------------------------------------------

 !----------------------------------------------------------------------------
subroutine Test_BLM (mm,dd,hh,fH,u,v, zm, zb, pb, exnm, &
          th, q,  Kz, Kz_nwp, invL, ustar, zi )
  integer, intent(in) :: mm, dd, hh        ! date
  real, intent(in)               :: fh     ! heart flux, -ve = Unstable
  real, dimension(:), intent(in) :: u,v    ! winds
  real, dimension(:), intent(in) :: exnm   ! mid-cell exner function (CP*)
  real, dimension(:), intent(in) :: zm     ! mid-cell height
  real, dimension(:), intent(in) :: zb     ! cell boundary height
  real, dimension(:), intent(in) :: pb     ! pressure at boundaries
  real, dimension(:), intent(in) :: th     ! pot. temp
  real, dimension(:), intent(in) :: q      !  specific humid ! TEST Vogel
  real, dimension(:), intent(in) :: Kz     ! Kz  (m2/s) 
  real, dimension(:), intent(in) :: Kz_nwp ! Kz from NWP if available/used
  real, intent(in)               :: ustar  ! m/s
  real, intent(in)               :: invL   ! 1/m
  real, intent(in)               :: zi     ! pot. temp
  integer :: k

  real, dimension(size(Kz)) :: &
    Kz_OB   &! Kz  O'Brien
   ,Kz_AJ   &! Kz  Amela Jericevic
   ,Kz_BW   &! Kz  Brost-Wynaargd, unstable
   ,Kz_PBT  &! Kz  Pielke+Blackader, Pielke flag=T
   ,Kz_PBF   ! Kz    "   "   flag=F
  real :: ziSeibert, ziJericevic, ziVenki, ziTI, ziVH

    write(*,*)"HmixMETHOD "//PBL%HmixMethod
    write(*,*)"KzMETHOD "//KzMethod//"-U:"//UnstableKzMethod// &
              "-S:"//StableKzMethod

    call PielkeBlackadarKz (u,v, zm, zb, th, Kz_PBT, &
              Pielke_flag=.true., debug_flag=.false.)

    call PielkeBlackadarKz (u,v, zm, zb, th, Kz_PBF, &
              Pielke_flag=.false., debug_flag=.false.)

  !/ Hmix methods *************************************

    call SeibertRiB_Hmix( u,v, zm, th, ziSeibert)

    call JericevicRiB_Hmix (u,v, zm, th, ziJericevic)

    call VogelezangHoltslag_Hmix (u,v, zm, th, q,  ustar, ziVH)

    ziVenki = Venkatram_Hmix(ustar)

    call TI_Hmix(Kz_PBT, zm, zb, fh, th, exnm, pb, ziTI, debug_flag=.true.)

    write(*,"(a,3i3,f8.2,10(a,f5.0))") "TEST_BLM fh:", mm, dd, hh, fh, &
          " zi ", zi, " ziS: ", ziSeibert, " ziJ: ", ziJericevic, &
          " ziV: ", ziVenki, " ziTI: ", ziTI, " ziVH: ", ziVH

  !/ Kz *************************************************
    write(*,"(a,4a3,2a7,a9,9a10)") "DEBUG_Kz: ", "mm", "dd", "hh", "k", &
           "fh", "u*", "zb", "pzpbl", &
           "Kz_m2s", "Kz_nwp", "Kz_PBT", "Kz_PBF", "Kz_OB", "KBW", "KAJ"

    Kz_OB(:) = Kz_PBT(:) ! sim to orig emep
    Kz_BW(:) = Kz_PBT(:) ! sim to orig emep

    if( fh < 0 ) then
       if ( invL > 0 ) print *, "TEST BLM SIGN ERRORR!!", fh, invL
       call O_BrienKz( zi, zb, ustar, invL, Kz_OB(:),.false.)
    end if

    do k = 2, size(th)

       if ( zb(k) < zi  ) then
          Kz_AJ(k) = JericevicKz( zb(k), ziJericevic, ustar, -8.888  )
       else
          Kz_AJ(k) = -9.999
       end if


       if( fh > 0 ) then  ! Query choices above zi?
           if( zb(k) < ziSeibert ) then ! Query
               Kz_BW(k) = BrostWyngaardKz( zb(k), ziSeibert, ustar, invL, -8.888  )
           else
               Kz_BW(k) = KZ_MINIMUM
           end if
       end if

       write(*,"(a,4i3,f7.1,f7.3,2f8.0,9es10.2)") "DEBUG_Kz: ", &
        mm, dd, hh, k, fh, ustar, zb(k), zi, Kz(k), Kz_nwp(k), Kz_PBT(k), &
           Kz_PBF(k), Kz_OB(k), Kz_BW(k), Kz_AJ(k)

    end do

    ! Write out as arrays for plotting:
    write(*,"(a,3i3,25f7.1)") "KZMAT Kz_nwp ", mm, dd, hh, zi, &
            (Kz_nwp(k),k=5,20)
    write(*,"(a,3i3,25f7.1)") "KZMAT Kz_PBT ", mm, dd, hh, zi, &
            (Kz_PBT(k),k=5,20)
    write(*,"(a,3i3,25f7.1)") "KZMAT Kz_AJ  ", mm, dd, hh, ziJericevic, &
            (Kz_AJ(k), k=5,20)

end subroutine Test_BLM

subroutine TI_Hmix (Kz, zm, zb, fh, th, exnm, pb, zi, debug_flag)
  real, dimension(:), intent(in) :: Kz     ! Kz (m2/s units) 
  real, dimension(:), intent(in) :: zm, zb ! heights of mid and bnd
  real, intent(in) :: fh  ! surface flux of sensible heat, W/m2
  real, dimension(:), intent(in) :: th
  real, dimension(:), intent(in) :: exnm, pb ! exner_mid and p_bnd
  real, intent(out) :: zi  ! Final height of ABL
  logical, intent(in) :: debug_flag

  real, dimension(KMAX_BND) :: &
      xksm   ! spacially smoothed Kz in z direction, m2/s.

  integer :: k, km, km1, km2, kp, nh1, nh2

  real, parameter :: DTZ = 3600.0 !time interval for integration of surface 
                                  !heat fluxes in ABL-height calculations, s
  real :: zis    ! height of the stable ABL, m

  ! For unstable BL
  real :: delq   ! available heat flux for developing the unstable ABL, J/m2
                 ! (heat-input per m2 from the ground during unstable BL)
  real :: thsrf, xdth, dthdzm, dthc, xdthdz
  real, dimension(KMAX_MID) :: thc, dthdz ! Assures th increases with height
  real :: dpidth !heat increasement in accordance with temp. increasement, J/m2
  real :: pidth  !heat used to adjust air temperature, J/m2
  integer :: trc !help variable telling whether or not unstable ABL exists
  integer :: kabl
  real :: ziu    ! height of the unstable ABL, m

  ! Smooth Kz verticall over 3 levels

  k=2
  km=1
  kp=3
  xksm(k)=( (zm(km)-zm(k))*Kz(k) + (zm(k)-zm(kp))*Kz(kp) )&
              / ( zm(km) - zm(kp) )

  k=KMAX_MID
  km2=k-2
  km1=k-1
  xksm(k)=( (zm(km2)-zm(km1))*Kz(km1) + (zm(km1)-zm(k))*Kz(k) )&
               / ( zm(km2) - zm(k) )

  do k = 3,KMAX_MID-1
      km1=k-1
      km2=k-2
      kp=k+1
      xksm(k)=(  (zm(km2)-zm(km1))*Kz(km1) + (zm(km1)-zm(k))*Kz(k)&
            + (zm(k)-zm(kp))*Kz(kp) ) / ( zm(km2) - zm(kp) )
  end do ! k

 !............................................................
 !..The height of the stable BL is the lowest level for which:
 !..xksm .le. 1 m2/s (this limit may be changed):
 
  zis = PBL%ZiMIN
  nh1 = KMAX_MID
  nh2 = 1

  do k=KMAX_MID,2,-1

     if(xksm(k) >= KZ_SBL_LIMIT .and. nh2 == 1) then
        nh1=k   ! Still unstable
     else
        nh2=0   ! Now stable
     end if
  end do

  k=nh1
  if(zb(nh1) >=  PBL%ZiMIN) then

      if( abs(xksm(k)-xksm(k-1)) > eps) then

           zis=((xksm(k)-KZ_SBL_LIMIT )*zb(k-1) &
               + (KZ_SBL_LIMIT -xksm(k-1))*zb(k))&
                     /(xksm(k)-xksm(k-1))
      else
          zis= PBL%ZiMIN
      end if

   end if


   zi = zis

  !..height of unstable boundary layer:
  !
  !..assuring that th is increasing with height.
  !..adjusted th-sounding is assigned to thc-array.
  !..This adjusted th is not meant to be used in
  !..other parts of the model program
  !
  dthdzm = 1.e-4
  thc(KMAX_MID)=th(KMAX_MID)
  do k=KMAX_MID-1,1,-1

      dthc = (th(k)-th(k+1))&
             / (zm(k)-zm(k+1))

      dthdz(k)=max(dthc,dthdzm)

      thc(k)=thc(k+1)+dthdz(k)*(zm(k)-zm(k+1))
        ! if ( debug_flag ) write(6,"(a,i3,3es10.3)") "DEBUG THC ", 
        !  k, th(k), dthc, thc(k)

  end do

  !..estimated as the height to which an hour's input
  !..of heat from the ground is vertically distributed,
  !..assuming dry adiabatic adjustment.

  delq=-min(fh,0.0)*DTZ
     thsrf=0.0
     ziu=0.0
     if ( debug_flag ) write(6,"(a,3es10.3)") "DEBUG fH (-ve => U) ", fh

     !.................................

     trc=0                        !..   =0 for stable BL (delq=0):
     if(delq >  0.00001) trc=1    !..trc=1 for unstable BL (delq>0):

     !------------------------------------------------------------
     ! calculating the height of unstable ABL

    kabl = KMAX_MID
    do while( trc == 1)
        kabl = kabl-1
        pidth=0.

        do k=KMAX_MID,kabl,-1
           xdth = thc(kabl)-thc(k)
           dpidth = exnm(k)*xdth*(pb(k+1)-pb(k))/GRAV
           pidth = pidth + dpidth
          ! if ( debug_flag ) write(6,"(a,2i3,6es11.3,i4)") "DEBUG PID ",&
          ! kabl, k, xdth, exnm(k), pb(k),  dpidth, pidth 
        end do

       if(pidth >= delq.and.trc == 1  ) then

          !  at level kabl or below level kabl and above level kabl+1

          thsrf = thc(kabl)- (thc(kabl)-thc(KMAX_MID))    &
                * (pidth-delq)/pidth

          xdthdz = (thc(kabl)-thc(kabl+1)) / (zm(kabl)-zm(kabl+1))

          ziu = zm(kabl+1) + (thsrf-thc(kabl+1))/xdthdz

          trc=0  
          !  if ( debug_flag ) write(6,"(a,i3,2es10.3,i4)") "DEBUG PICTH ",&
          !      kabl, delq, pidth, trc 

        end if


        if ( debug_flag ) write(6,"(a,i3,es10.3,i5)") "DEBUG mid ", &
             kabl, delq, trc 
        if(kabl <= 4 .and. trc == 1  ) then

           write(6,*)'PBL ziu calculations failed!', fh
           !if ( debug_flag ) then
              do k=KMAX_MID,kabl,-1
                 write(6,"(a,i3,3es10.3)") "DEBUG ZIU ", k, th(k), thc(k)
              end do
           !end if
           

           ziu=PBL%ZiMAX

           trc=0 
        end if

     end do ! while

  zi = max( ziu, zis)
  zi = min( PBL%ZiMAX, zi)


end subroutine TI_Hmix

 !----------------------------------------------------------------------------

    !************************************************************************!
      subroutine O_BrienKz( zi, zs_bnd, ustar, invL , Kz, debug_flag )       !
    !**********************************************************************!
    !..exchange coefficients for convective boundary layer:
    !..o'brien's profile formula:

    real,intent(in) :: zi
    real,intent(in), dimension(:) :: zs_bnd
    real,intent(in) :: ustar, invL
    real, intent(inout), dimension(:) :: Kz
    logical, intent(in) :: debug_flag

    real, parameter :: FHSL = 0.04 ! surface layer as fraction zi

    real :: ux3 ,hsl  &
         ,zimhs ,zimz ,zmhs
    real  ::  Kzhs  ,dKzdz ,Kzzi ,hs
    integer :: k

    !c..exchange parameter and its vertical derivative at z = hs

      Kzhs=0.    ! Kz at hs
      dKzdz=0.   ! d Kz(hs) /dz - derivateive at hs (or KB)
      Kzzi=0.    ! Kz at top of BL, set to zero or KZ_MINIMUM?

      !...................................................................
      !..air density at ground level is always calculated diagnostically:

      ux3 = ustar*ustar*ustar

     !..........................
     !..unstable surface-layer.:

     hs=FHSL*zi          !  height of surface layer

     !c..hsl=hs/l where l is the monin-obhukov length
     !hsl = KARMAN*GRAV*hs*fh*KAPPA /(ps*ux3)

     hsl = hs*invL

     ! Garratt \Phi function - take from MicroMet

     Kzhs = ustar*KARMAN*hs*sqrt(1.0-16.0*hsl)  ! /Pr=1.00
     dKzdz = Kzhs*(1.-0.5*16.0*hsl/(1.0-16.0*hsl))/hs

     Kz(KMAX_MID)=Kzhs  ! QUERY - should be at z_bnd(20)?
     if ( debug_flag ) write(*,"(a,f7.2,3es12.3)") "OBRIEN Kz20 ", &
        hs, invL, Kzhs, dKzdz
     if ( 16.0*hsl >= 1.0 ) then
        write(*,"(a,f7.2,4es12.3)") "OBRIEN NEG ", hs, hsl,  Kzhs, dKzdz
     end if

   !..exchange parameter at z = ziu

   !do  k=1,KMAX_MID !had -ve values at k=1, and not needed
   do  k=2,KMAX_MID

         if( zs_bnd(k) >= zi) then

            Kzzi=Kz(k)  ! values above zi  stored

            if ( debug_flag ) write(*,"(a,i3,es12.3)") "OBRIEN Kzzi ", k,  Kzzi

         else 
            !.....................................................
            !..the obrien-profile for z<ziu                      .
            !.....................................................
            !
            if(zs_bnd(k) <= hs) then
               Kz(k)=zs_bnd(k)*Kzhs/hs
               if ( debug_flag ) &
                   write(*,"(a,i3,es12.3)") "OBRIEN Kzhs ", k, Kz(k)
            else !! if( zs_bnd(k) < zi) then
               zimhs = zi-hs
               zimz  =zi-zs_bnd(k)
               zmhs  =zs_bnd(k)-hs
               Kz(k) = Kzzi+(zimz/zimhs)*(zimz/zimhs)  &
                    *(Kzhs-Kzzi+zmhs*(dKzdz     &
                    + 2.*(Kzhs-Kzzi)/zimhs))
               if ( debug_flag ) &
                   write(*,"(a,i3,es12.3)") "OBRIEN Kz(k) ", k,  Kz(k)
            end if

         end if

   end do

  end subroutine O_BrienKz
 !----------------------------------------------------------------------------
  function risig1(ws,th1,th2,z) ! Amela Jericevic, not used in new version
   !calculates the bulk Richardson number
      implicit none
      real, intent(in) :: ws   ! wind-speed
      real, intent(in) :: th1, th2
      real, intent(in) :: z
      real :: risig1
      real :: dvdz

      dvdz = ws*ws+0.001
      risig1=(2.0*GRAV/(th1+th2))*(th2-th1)*z/dvdz
   end function risig1
 !----------------------------------------------------------------------------

subroutine fake_zmid(z)
   
 ! Data for testing. Note that z here is really geopotetntial
 ! heught, but near the ground it doesn't matter so much
 ! from 1 time-period, see DEBUG_Z in Met_mod
   real, dimension(20), intent(out) :: z
   z( 1) =   15244.7131
   z( 2) =   13531.2981
   z( 3) =   12177.9740
   z( 4) =   11004.0286
   z( 5) =    9788.6193
   z( 6) =    8430.8546
   z( 7) =    7013.5860
   z( 8) =    5700.0226
   z( 9) =    4598.2893
   z(10) =    3692.3391
   z(11) =    2934.3594
   z(12) =    2296.2730
   z(13) =    1761.6936
   z(14) =    1317.4565
   z(15) =     951.4678
   z(16) =     655.8842
   z(17) =     424.7443
   z(18) =     254.0498
   z(19) =     137.4284
   z(20) =      45.6146
end subroutine fake_zmid

subroutine fake_zbnd(z)
   
 ! Data for testing. Note that z here is really geopotetntial
 ! heught, but near the ground it doesn't matter so much
 ! from 1 time-period, see DEBUG_Z in Met_mod
   real, dimension(21), intent(out) :: z
  z( 1) =    16276.236
  z( 2) =    14319.868
  z( 3) =    12815.771
  z( 4) =    11598.220
  z( 5) =    10461.131
  z( 6) =     9168.470
  z( 7) =     7747.153
  z( 8) =     6324.928
  z( 9) =     5103.225
  z(10) =     4110.656
  z(11) =     3286.077
  z(12) =     2591.897
  z(13) =     2007.792
  z(14) =     1520.443
  z(15) =     1116.948
  z(16) =      787.332
  z(17) =      525.163
  z(18) =      324.823
  z(19) =      183.602
  z(20) =       91.302
  z(21) =       91.302
end subroutine fake_zbnd
subroutine SigmaKz_2_m2s_scalar (roa,ps,Kz_fac)   ! hb
  real, intent(in) :: roa
  real, intent(in) :: ps
  real :: fac
  real :: Kz_fac

     fac= (ps - PT)/(GRAV*roa)
     Kz_fac= fac*fac

end subroutine SigmaKz_2_m2s_scalar

subroutine SigmaKz_2_m2s_arrays (SigmaKz,roa,ps,Kz)
  real, intent(in), dimension(:,:,:) :: SigmaKz, roa
  real, intent(in), dimension(:,:)   :: ps
  real, intent(out), dimension(:,:,:) :: Kz 
  real :: fac
  integer :: i,j,k

 ! Kz has dim 1:KMAX_MID, whereas SigmaKz has 1:KMAX_BND
  do k = 1, size(Kz,3)
    do j = 1, size(Kz,2)
      do i = 1, size(Kz,1)
         fac= (ps(i,j) - PT)/(GRAV*roa(i,j,k))
         Kz(i,j,k) = fac*fac*SigmaKz(i,j,k)
      end do
    end do
  end do

end subroutine SigmaKz_2_m2s_arrays

subroutine Kz_m2s_toSigmaKz (Kz,roa,ps,SigmaKz)
  real, intent(in), dimension(:,:,:) :: Kz, roa
  real, intent(in), dimension(:,:)   :: ps
  real, intent(out), dimension(:,:,:) :: SigmaKz 
  real :: fac
  integer :: i,j,k

 ! Kz has dim 1:KMAX_MID, whereas SigmaKz has 1:KMAX_BND
  do k = 1, size(Kz,3)
    do j = 1, size(Kz,2)
      do i = 1, size(Kz,1)
         fac= (GRAV*roa(i,j,k))/(ps(i,j) - PT)
         SigmaKz(i,j,k) = fac*fac*Kz(i,j,k)
      end do
    end do
  end do
  k=size(Kz,3)+1  ! k=21 for Sigma's
  SigmaKz(:,:,k) = 0.0

end subroutine Kz_m2s_toSigmaKz

subroutine Kz_m2s_toEtaKz (Kz,roa,ps,EtaKz,Eta_mid,A_mid,B_mid)
  real, intent(in), dimension(:,:,:) :: Kz, roa
  real, intent(in), dimension(:,:)   :: ps
  real, intent(out), dimension(:,:,:) :: EtaKz 
  real, intent(in), dimension(:)   :: Eta_mid,A_mid,B_mid
  real :: fac
  integer :: i,j,k

! Kz has dim 1:KMAX_MID, whereas EtaKz has 1:KMAX_BND
! Kz defined at middle of level, EtaKz defined at level boundaries
! EtaKz = Kz*(roa*g* d(Eta)/d(P) )**2
  do k = 2, size(Kz,3)
    do j = 1, size(Kz,2)
      do i = 1, size(Kz,1)
         fac= (GRAV*(roa(i,j,k)+roa(i,j,k-1))*0.5)*&
         (Eta_mid(k)-Eta_mid(k-1))/(A_mid(k)+B_mid(k)*ps(i,j)-A_mid(k-1)-B_mid(k-1)*ps(i,j))
         EtaKz(i,j,k) = fac*fac*Kz(i,j,k)
      end do
    end do
  end do
  k=size(Kz,3)+1  ! surface
  EtaKz(:,:,k) = 0.0
  EtaKz(:,:,1) = 0.0!top

end subroutine Kz_m2s_toEtaKz
end module BLPhysics_mod
