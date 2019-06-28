! <Radiation_mod.f90 - A component of the EMEP MSC-W Chemical transport Model, version rv4.33>
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
module Radiation_mod

  !>
  !! Collection of routines to calculate radiation terms, also for
  !! canopies. IMPORTANT - Most routines expect SolarSetup to 
  !! have been called first.
  !
  !  Module usable by stand-alone deposition code.

  use PhysicalConstants_mod  , only: PI, DEG2RAD, RAD2DEG, DAY_ZEN, DAY_COSZEN
  use TimeDate_mod  , only: julian_date, day_of_year
  implicit none
  private

  !/ Subroutines:
  public :: SolarSetup    !> => decl, sindecl, eqt_h, etc., + daytime, solarnoon
  public :: ZenithAngle   !> => CosZen=cos(Zen), Zen=zenith angle (degrees) 
  public :: ZenithAngleS  !>    (simpler version)
  public :: ClearSkyRadn  !> => irradiance (W/m2), clear-sky
  public :: WeissNormanPAR   ! OCT/NOV  2017 after RIVM found bug
  public :: CloudAtten    !> => Cloud-Attenuation factor
  public :: fCloudAtten    !> => Cloud-Attenuation factor
  public :: CanopyPAR     !> => sun & shade PAR  values, and LAIsunfrac
  public :: ScaleRad      !>  Scales modelled radiation where observed values
                          !!  available.

  !/ Functions:
  public :: daytime       !> true if zen < 89.9 deg
  public :: daylength     !> Lenght of day, hours
  public :: solarnoon     !> time of solarnoon, hours
  public :: sunrise       !> time of sunrise, hours
  public :: sunset        !> time of sunset, hours


  real, public, parameter :: &
      PARfrac = 0.45,   &           !> approximation to fraction (0.45 to 0.5)
                        !! of total radiation in PAR waveband (400-700nm)
      Wm2_uE  = 4.57,   &           !> converts from W/m^2 to umol/m^2/s
      Wm2_2uEPAR= PARfrac * Wm2_uE  !> converts from W/m^2 to umol/m^2/s PAR


  ! Some variables which are dependent only on day of year and GMT time
  ! - hence they do not vary from grid to grid and can safely be stored here

  real, private, save :: rdecl,sinrdecl,cosrdecl
  real, private, save :: eqtime, eqt_h, tan_decl

logical, private, parameter :: DEBUG = .false.


  !================== Ashrae/Iqbal stuff  -- see ClearSkyRad subroutine

    type, public :: ashrae_tab
        integer :: nday
        real    :: a
        real    :: b
        real    :: c
    end type ashrae_tab

    type(ashrae_tab), save, public ::  Ashrae    ! Current values

    type(ashrae_tab), parameter, dimension(14), private ::  ASHRAE_REV = (/ &
    !               nday    a      b      c
         ashrae_tab(  1, 1203.0, 0.141, 0.103 ) &
        ,ashrae_tab( 21, 1202.0, 0.141, 0.103 ) &
        ,ashrae_tab( 52, 1187.0, 0.142, 0.104 ) &
        ,ashrae_tab( 81, 1164.0, 0.149, 0.109 ) &
        ,ashrae_tab(112, 1130.0, 0.164, 0.120 ) &
        ,ashrae_tab(142, 1106.0, 0.177, 0.130 ) &
        ,ashrae_tab(173, 1092.0, 0.185, 0.137 ) &
        ,ashrae_tab(203, 1093.0, 0.186, 0.138 ) &
        ,ashrae_tab(234, 1107.0, 0.182, 0.134 ) &
        ,ashrae_tab(265, 1136.0, 0.165, 0.121 ) &
        ,ashrae_tab(295, 1136.0, 0.152, 0.111 ) &
        ,ashrae_tab(326, 1190.0, 0.144, 0.106 ) &
        ,ashrae_tab(356, 1204.0, 0.141, 0.103 ) &
        ,ashrae_tab(366, 1203.0, 0.141, 0.103 ) &
      /)
  !=============================
  

contains

 !<===========================================================================
  subroutine SolarSetup(year,month,day,hour)

    ! Sets up decelention and related terms, as well as Ashrae coefficients
    ! Should be called before other routines.

    integer, intent(in) :: year,month,day
    real, intent(in) :: hour

    real :: d,ml,rml,w,wr,ec,epsi,yt,pepsi,cww
    real :: sw,ssw, eyt, feqt1, feqt2, feqt3, feqt4, feqt5, &
            feqt6, feqt7, feqt,ra,reqt  
    real :: dayinc
    integer :: i

    logical, parameter :: MY_DEBUG = .false.


!* count days from dec.31,1973 

    d = julian_date(year,month,day) - julian_date(1973,12,31) + 1 
    d = d + hour/24.0


!* calc geom mean longitude

    ml = 279.2801988 + 0.9856473354*d + 2.267e-13*d*d
    rml = ml*DEG2RAD

!* calc equation of time in sec
!*  w = mean long of perigee
!*  e = eccentricity
!*  epsi = mean obliquity of ecliptic

    w = 282.4932328 + 4.70684e-5*d + 3.39e-13*d*d
    wr = w*DEG2RAD
    ec = 1.6720041e-2 - 1.1444e-9*d - 9.4e-17*d*d
    epsi = 23.44266511 - 3.5626e-7*d - 1.23e-15*d*d
    pepsi = epsi*DEG2RAD
    yt = tan(pepsi/2.0)
    yt = yt*yt
    cww = cos(wr)
    sw = sin(wr)
    ssw = sin(2.0*wr)
    eyt = 2.0*ec*yt
    feqt1 = sin(rml)*(-eyt*cww - 2.0*ec*cww)
    feqt2 = cos(rml)*(2.0*ec*sw - eyt*sw)
    feqt3 = sin(2.0*rml)*(yt - (5.0*ec*ec/4.0)*(cww*cww-sw*sw))
    feqt4 = cos(2.0*rml)*(5.0*ec*ec*ssw/4.0) 
    feqt5 = sin(3.0*rml)*(eyt*cww)
    feqt6 = cos(3.0*rml)*(-eyt*sw)
    feqt7 = -sin(4.0*rml)*(0.5*yt*yt)
    feqt = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7

    eqtime = feqt*13751.0

!* equation of time in hrs:

    eqt_h = eqtime/3600.0

!* convert eq of time from sec to deg

    reqt = eqtime/240.0

!* calc right ascension in rads

    ra = ml - reqt

!* calc declination in rads, deg

    tan_decl = 0.43360*sin( ra * DEG2RAD )

    rdecl    = atan(tan_decl)
    sinrdecl = sin(rdecl)
    cosrdecl = cos(rdecl)

  !-----------------------------------------------------------------

     ! Find coefficients for Iqbal/Ashrae algorith, used in ClearSkyRadn routine
     ! first, perform the table look up

      d = day_of_year(year,month,day)
      do i = 1, 14
        if (d <=   ASHRAE_REV(i)%nday ) exit
      end  do

      if ( DEBUG .and. i<1.or.i>14) write(unit=6,fmt=*) "solbio: index err!"

      if ( ASHRAE_REV(i)%nday == 1) then
        Ashrae%a = ASHRAE_REV(1)%a
        Ashrae%b = ASHRAE_REV(1)%b
        Ashrae%c = ASHRAE_REV(1)%c
      else
        dayinc = real( d-ASHRAE_REV(i-1)%nday ) / &
                 real( ASHRAE_REV(i)%nday-ASHRAE_REV(i-1)%nday )
        Ashrae%a = ASHRAE_REV(i-1)%a + &
             ( ASHRAE_REV(i)%a - ASHRAE_REV(i-1)%a )*dayinc
        Ashrae%b = ASHRAE_REV(i-1)%b + &
             ( ASHRAE_REV(i)%b - ASHRAE_REV(i-1)%b )*dayinc
        Ashrae%c = ASHRAE_REV(i-1)%c + &
             ( ASHRAE_REV(i)%c - ASHRAE_REV(i-1)%c )*dayinc
      end if

 end subroutine SolarSetup


 !<===========================================================================
  elemental subroutine ZenithAngle(thour, latitude, longitude, Z, CosZ )
  !   IMPORTANT - Call SolarSetup before use to get decl terms and eqt_H

    real, intent(in) :: thour
    real, intent(in) :: latitude
    real, intent(in) :: longitude
    real, intent(out) :: Z                  ! Zenith Angle (degrees)
    real, intent(out) :: CosZ               ! cos(Z)

    real :: rlt,lzgmt,zpt,lbgmt

        rlt = latitude *DEG2RAD
        lbgmt = 12.0 - eqt_h - longitude *24.0/360.0
        lzgmt = 15.0*(thour - lbgmt)
        zpt = lzgmt*DEG2RAD
        CosZ = sin(rlt)*sinrdecl + cos(rlt)*cosrdecl*cos(zpt)
    
        CosZ = min( 1.0, CosZ)
        CosZ = max(-1.0, CosZ)
        Z = acos(CosZ)*RAD2DEG

 end subroutine ZenithAngle
  ! ======================================================================
  subroutine ZenithAngleS(lon, lat, daynr, nydays, hr, Z, CosZ )
  ! ======================================================================
  ! (Non-elemental version of ZenithAngle, used for debugging sometimes)
  !  routine determines (approximate) cos(zen), where "zen" denotes the zenith 
  !  angle, (in accordance with Iversen and Nordeng (1987, p.28))
  !  dnmi  29-9-95  Hugo Jakobsen, modified Dave, 2002-2004
  !
  ! arguments:              
    real,    intent(in) ::  lon      ! longitude (degrees), east is positive
    real,    intent(in) ::  lat      ! latitude (degrees), north is positive
    integer, intent(in) ::  daynr    ! day nr. (1..366)
    integer, intent(in) ::  nydays   ! number of days per year (365 or 366)
    real,    intent(in) ::  hr       ! hour  (0-24, gmt)  ! ds - was integer
    real,    intent(out) ::  Z       ! zenith angle (degrees) 
    real,    intent(out) ::  CosZ
    
     !/ Local....
     real    :: lonr, latr, arg, decl, tangle

     lonr=lon*DEG2RAD             ! convert to radians
     latr=lat*DEG2RAD             ! convert to radians

     arg = ((daynr - 80.0)/nydays) * 2.0 * PI

     decl = 23.5 * sin(arg) * DEG2RAD
     
     tangle = lonr + (hr/12.0-1.0)*PI      !no eqtime correction
     CosZ   =(sin(decl)*sin(latr)+cos(decl)*cos(latr)*cos(tangle))
     Z      = acos(CosZ) * RAD2DEG 
  end subroutine ZenithAngleS

 !=============================================================================
 ! After Weiss & Norman 1985 (WN below), but use notation 'dbh' for direct beam
 ! on horizontal surface to avoid confusion with irradiance (normal to beams)
 ! PAR is here in W/m2

  elemental subroutine WeissNormanPAR(p,CosZ,fcloud,Svdbh,Svdif)
    real, intent(in) :: p                  ! Pressure, Pa
    real, intent(in) :: CosZ               ! Cos(Zenith)
    real, intent(in) :: fcloud   ! cloud factor (RATIO in WN), clear sky=1.
  ! Outputs clear-sky values of PAR:
    real, intent(out) ::  Svdbh   ! Direct beam on horizontal surface,  W/m2
    real, intent(out) ::  Svdif   ! Diffuse, W/m2

   ! R terms represent potential (clear-sky) radiation:
    real ::  Rvdbh                ! Direct PAR beam on horizontal surface, W/m2
    real ::  Rvdif                ! Diffuse PAR, W/m2
    real ::  Rvh, Svh             ! Total visible on horizontal surface, W/m2
    real ::  fvdbh                ! fraction vis, direct beam
    real ::  m                    ! optical air mass, 1/cos(theta)
    real, parameter  :: PRES0 = 101300.0  ! std sea-level pressure (Pa)
    real, parameter  :: A=0.9, B = 0.7    ! WN coefficients
  
    if ( CosZ > DAY_COSZEN ) then

      m=1/cosZ
      Rvdbh = 600 * exp( -0.185*(p/PRES0) *m ) * cosZ   !WN1
      Rvdif = 0.4*( 600 - Rvdbh ) * cosZ                !WN2
      Rvh   = Rvdbh + Rvdif                             !WN9

      if ( fcloud <  A ) then
        fvdbh = Rvdbh/Rvh * ( 1- ((A-fcloud)/B)**0.6667 ) ! !WN11
      else
        fvdbh = Rvdbh/Rvh    !WN11
      end if

      Svh   = fcloud * Rvh      ! C7

      Svdbh = fvdbh *  Svh      ! C12
      Svdif = Svh - Svdbh       ! C12
      
    else
      Svdbh  = 0.0
      Svdif  = 0.0
    end if

  end subroutine WeissNormanPAR
 !=============================================================================

  elemental subroutine ClearSkyRadn(p,CosZ,Idirect,Idiffuse)

  !     Computes the radiation terms needed for stomatal and BVOC calculations.
  !     Methodology for this calculation taken from Iqbal, M., 1983,
  !     An introduction to solar radiation, Academic Press, New York, 
  !     pp. 202-210.
  !
  !  history:
  !
  !     From T. Pierce's SolBio code:
  !     Development of this routine was prompted by the need for a
  !     horizontal rather than an actinic flux calculation (which had
  !     been performed by Soleng). Furthermore, Soleng computed total
  !     radiation only out to the near-ir spectrum. This program
  !     is designed only for approximate radiation estimates to be used
  !     for stomatal calculations.
  !
  !     8/90    initial development of SolBio by T. Pierce
  !     95       modified by Hugo Jakobsen, 29/9-95
  !     04-05    modified by Dave Simpson for EMEP and DO3SE models

    real, intent(in) :: p                  ! Pressure, Pa
    real, intent(in) :: CosZ               ! Cos(Zenith)

  ! Calculates clear-sky values of:

    real, intent(out) :: Idirect  ! total direct solar radiation (W/m2)
    real, intent(out) :: Idiffuse ! diffuse solar radiation (W/m)

  ! Local:
    real ::  Idrctn               ! direct normal solar radiation (W/m2)
    real, parameter  :: PRES0 = 101300.0  ! std sea-level pressure (Pa)

    real, parameter  :: cn = 1.0
  !        cn    - clearness number (defined as the ratio of normal
  !                incident radiation, under local mean water vapour,
  !                divided by the normal incident radiation, for
  !                water vapour in a basic atmosphere)  
  !              - currently, this value set equal to 1.0, but eventually
  !                may vary as a function of latitude and month pending further
  !                literature review.

    if ( CosZ > DAY_COSZEN ) then
      Idrctn     = cn * Ashrae%a * exp(- Ashrae%b * (p/PRES0)/CosZ)
      Idiffuse   = Ashrae%c * Idrctn
      Idirect    = Idrctn *  CosZ
    else
      Idirect    = 0.0
      Idiffuse   = 0.0
    end if

   !X solar  = Idirect + Idiffuse ! total solar radiation, diff.+direct (W/m2)

  end subroutine ClearSkyRadn
!===========================================================================
  elemental subroutine CloudAtten(cl,a,b)!,c)
    ! Routine applies a cloud-attenuation factor to arguments, which could
    ! be say, Idrctt,Idfuse,solar, or just solar: the last 2 arguments are
    ! optional

    ! Agument
      real, intent(in)  :: cl               ! cloud fraction   (0-1)
      real, intent(inout)           :: a     
      real, intent(inout), optional :: b!,c

      real :: f           ! cloud attenuation factor

      f = 1.0 - 0.75*cl**3.4     !(source: Kasten & Czeplak (1980)) 

      a = a * f

      if( present(b) ) b = b * f
!      if( present(c) ) c = c * f

  end subroutine CloudAtten
!===========================================================================
  elemental function fCloudAtten(cl) result(f)
    real, intent(in) :: cl
    integer :: icl
    real :: f 
    real, parameter, dimension(0:100) :: fcloud = [ (1-0.75*(0.01*icl)**3.4, icl=0,100) ]
      ! Varies :
      !  10  0.999701440    30  0.987489522    50  0.928950787    
      !  60  0.867938757    70  0.776953936    80  0.648789823    
      !  90  0.475813627   100  0.250000000    
    f = fcloud( int(100*cl) )
  end function fCloudAtten
!===========================================================================
    !WN17 subroutine CanopyPAR(LAI,sinB,Idrctt,Idfuse,&
    !WN17 PARdbh, PARdif is in W/m2
    subroutine CanopyPAR(LAI,sinB,PARdbh,PARdif,&
                            PARsun,PARshade,LAIsunfrac)
!===========================================================================
!
!    Calculates sun and shade components of incoming PAR (W/m2 on horizontal
!    surface, using methodology as described in Emberson et  al. (2000,
!    EMEP Note 6/00), eqns.  6-11, based upon sun/shade method of Norman (1979,1982)
!    (N82 eqns. 14-17).

!     input arguments:

    real, intent(in)  :: LAI       ! leaf area index (m^2/m^2), one-sided
    real, intent(in)  :: sinB      ! B = solar elevation angle; sinB = CosZen
    real, intent(in)  :: PARdbh, PARdif     ! Direct, diffuse Radn, W/m2
    real, intent(out) :: PARsun, PARshade   ! vegetation components of PAR, W/m2
    real, intent(out) :: LAIsunfrac


!     internal variables:

    real :: LAIsun    ! sunlit LAI

    real, parameter :: cosA    = 0.5   ! A = mean leaf inclination (60 deg.), 
     ! where it is assumed that leaf inclination has a spherical distribution



    LAIsun = (1.0 - exp(-0.5*LAI/sinB) ) * sinB/cosA
    LAIsunfrac = LAIsun/LAI

  ! PAR flux densities (W/m2) evaluated using method of Norman (1982, p.79): 
  ! "conceptually, 0.07 represents a scattering coefficient"  

    PARshade = PARdif * exp(-0.5*LAI**0.7) +  &
               0.07 * PARdbh  * (1.1-0.1*LAI)*exp(-sinB)   

    PARsun = PARdbh *cosA/sinB + PARshade

  end subroutine CanopyPAR

!--------------------------------------------------------------------

  subroutine ScaleRad(ObsRad, Idrctt,Idfuse)
      real, intent(in) :: ObsRad
      real, intent(inout) :: Idrctt, Idfuse

      real :: Scale
      logical :: MY_DEBUG = .false.
              
   !-----  
   ! Observation frequently don't have PAR, but instead have global radiation.
   ! We scale the modelled global radiation (solar) by the observed, and
   ! distribute Idrctt and Idrctn according to this.

   ! This routine is from  a stand-alone version of the drydep, so not used
   ! in emepctm. The subroutine is retained in order to keep the standalone 
   ! and emepctm code the same

          Scale = -999 ! gfortran compilation gives a warning that scale is 
                       ! not set unless  ObsRad > 0.0 
      
          if ( ObsRad > 0.0 ) Scale =  ObsRad/(Idrctt+Idfuse)

          if (MY_DEBUG) then
              if (Scale <0.1 .or. Scale>10.0) then
                 print  "(a35,2f10.3)","Obs and Mod Radiation large diff", &
                     ObsRad,Idrctt+Idfuse
              end if
           end if
           Idrctt=Scale*Idrctt
           Idfuse=Scale*Idfuse
   end subroutine ScaleRad

  !=============================================================================
  ! FUNCTIONS
  !=============================================================================

  ! Define function for daytime, to keep definitions consistent throughout code.
  ! NB - older code had a check for zen>1.0e-15 --- Why?!

  elemental function daytime(zen) result (day)
       real, intent(in) :: zen    ! Zenith angle (degrees)
       logical :: day
 
       if( zen < DAY_ZEN ) then
            day = .true.
       else
            day = .false.
       end if
  end function daytime

  !-----------------------------------------------------------------
  ! Calculate length of day, following Jones, Appendix 7:
  !   IMPORTANT - Call SolarSetup before use to get tan_decl

  elemental function daylength(lat) result (len)
       real, intent(in) :: lat    !  Latitude, deg.
       real :: len, arg
       real, parameter :: TIMEFAC = RAD2DEG * 2.0 /15.0

        arg =  -tan_decl * tan( DEG2RAD*lat )

        if( arg <= -1.0 ) then 
           len = 24.0  !! Polar summer
        else if( arg >= 1.0 ) then 
           len = 0.0   !! Polar night
        else
           len = acos( -tan_decl * tan( DEG2RAD*lat ) ) &
                              * TIMEFAC  ! eqn A7.2, Jones
        end if
  end function daylength
  !-----------------------------------------------------------------
  ! Calculate solar noon, following Jones, Appendix 7:
  !   IMPORTANT - Call SolarSetup before use to get eqt_t

  elemental function solarnoon(long) result (noon)
       real, intent(in) :: long    !  Longitude, deg.
       real :: noon

        noon = 12.0 - eqt_h - long*24.0/360.0

  end function solarnoon
  !-----------------------------------------------------------------
  ! sunrise/sunset
  elemental function sunrise(lat,lon) result (hour)
    real, intent(in) :: lat,lon    !  latitude,Longitude, deg.
    real :: hour
    hour = solarnoon(lon)-daylength(lat)*0.5
   !hour = min(max(hour,0.0),24.0)
  end function sunrise
  elemental function sunset(lat,lon) result (hour)
    real, intent(in) :: lat,lon    !  latitude,Longitude, deg.
    real :: hour
    hour = solarnoon(lon)+daylength(lat)*0.5
   !hour = min(max(hour,0.0),24.0)
  end function sunset
!===============================================================
end module Radiation_mod
!===============================================================
!TSTEMX program testr
!TSTEMX use Radiation_mod
!TSTEMX use PhysicalConstants_mod, only:  DEG2RAD
!TSTEMX use TimeDate_mod  , only: day_of_year
!TSTEMX implicit none
!TSTEMX integer :: iZen,icl, beta, difFrac
!TSTEMX real :: CosZ,Idirect,Idiffuse,cl,fcloud, PARdbh, PARdif
!TSTEMX real :: PARdbh_cs, PARdif_cs    !clear-sky values
!TSTEMX call SolarSetup(2012,6,21,0.0)  ! Day 173
!TSTEMX print *, "DAY ", day_of_year(2012,6,21)
!
!TSTEMX do iZen = 90, 0, -10
!TSTEMX   print "(/,a3,a6,a3,5a10)", 'run', 'cl', 'Z', 'PARdbh_cs', 'PARdif_cs', 'PARdbh', 'PARdif', 'diff/tot'
!TSTEMX do icl = 0, 10
!TSTEMX    cl = 0.1 * icl
!TSTEMX    beta = 90.0 - iZen   ! Use solar elev angle, better for plots
!TSTEMX    CosZ = cos(iZen*DEG2RAD)
!TSTEMX    call ClearSkyRadn(1.0e5,CosZ,Idirect,Idiffuse)
!TSTEMX    PARdbh_cs = Idirect*PARfrac; PARdif_cs = Idiffuse *PARfrac! save for printout
!TSTEMX    !ORIG METHOD:
!TSTEMX    call CloudAtten(cl,Idirect,Idiffuse) !BUG SCALED BOTH DIR & DIF IDENTICALLY
!TSTEMX    difFrac = 100   ! Only for printout, but limit as night approaches
!TSTEMX    if( Idiffuse > 0.0) difFrac = nint( 100*Idiffuse/(Idirect+ Idiffuse))
!TSTEMX    fcloud = fCloudAtten(cl)
!TSTEMX    print "(a3,2f6.1,i3,4f10.3,i5)", 'old', cl, 100*fcloud, int(beta), PARdbh_cs, PARdif_cs, Idirect*PARfrac, Idiffuse*PARfrac, difFrac
!TSTEMX    !NEW  METHOD! Note that weissNorman uses PAR directly
!TSTEMX    call WeissNormanPAR(1.0e5,CosZ,fcloud,PARdbh,PARdif)
!TSTEMX    difFrac = 100   ! Only for printout, but limit as night approaches
!TSTEMX    if( PARdif >0) difFrac = nint(100*PARdif/(PARdbh+ PARdif))
!TSTEMX    print "(a3,2f6.1,i3,4f10.3,i5)", 'new', cl, 100*fcloud, int(beta), PARdbh_cs, PARdif_cs, PARdbh, PARdif, difFrac
!TSTEMX   end do
!TSTEMX end do
!TSTEMX print *, DEG2RAD
!TSTEMX end program testr
