! <Aqueous_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Aqueous_ml

!-----------------------------------------------------------------------
! Aqueous scavenging routines.
!
! A "minimal" version of the aqueous reactions.
! The scavenging of soluble compounds is based upon the work of Berge 
! and Jakobsen (1998) and Eliassen and Saltbones (1983).
! Simple scavenging coefficients are used. A distinction is made
! between in-cloud and sub-cloud scavenging.
!
! Usage:
!
! Uses My_WetDep_ml, where the scavenging coefficients, divided by
! relevant factors, are specified.
!
! Setup_Clouds(i,j)  called from Runchem_ml
! WetDeposition(i,j) called from Runchem_ml if prec. clouds are present
!
! Refs:
! Berge, E., 1993, Coupling of wet scavenging of sulphur to clouds in a
!    numerical weather prediction model, Tellus, 45B, 1-22
! Berge, E. and Jakobsen, H.A., 1998, {A regional scale multi-layer model for
!    the calculation of long-term transport and deposition of air pollution in
!    Europe, Tellus, 50,205-223
! Eliassen, A. and Saltbones, J., 1983, Modelling of long-range transport of
!    sulphur over Europe: a two year model run and some experiments, Atmos.
!    Environ., 17, 1457-1473
! Seland, O. and T. Iversen (1999) A scheme for black carbon and
!    sulphate aerosols tested in a hemispheric scale, Eulerian dispersion
!    model. Atm.  Env. Vol. 33, pp.2853-2879.
!-----------------------------------------------------------------------

  use My_WetDep_ml,  only : WetDep, NWETDEP, WetDep_Budget, WDEP_PREC
  use Derived_ml,    only : IOU_INST & ! Index: instantaneous values
                            ,d_2d      ! Contains Wet deposition fields
  use GridValues_ml, only : gridwidth_m,xm2,xmd,carea
  use ModelConstants_ml, only: &
      CHEMTMIN, CHEMTMAX       &       ! -> range of temperature 
     ,KMAX_MID                 &       ! -> ground, k=20
     ,KUPPER                   &       ! -> top of cloud-chemistry, k=6
     ,KCHEMTOP                 &       ! -> top of chemistry, now k=2
     ,dt => dt_advec           &       ! -> model timestep
     ,PT, ATWAIR                       ! -> pressure at top, atw. air
  use Met_ml,            only : pr, roa, z_bnd, cc3d, ps, lwc
  use Par_ml,            only : me     ! for DEBUG
  use Setup_1dfields_ml, only : xn_2d, amk

  implicit none
  private

! Subroutines:

  public  :: init_aqueous
  public  :: Setup_Clouds       ! characterises clouds and calls
                                ! WetDeposition if rain
  public  :: WetDeposition      ! simplified setup_wetdep
  private :: tabulate_aqueous
  private :: get_frac
  private :: setup_aqurates

! Outputs:

  logical, public, save, dimension(KUPPER:KMAX_MID) :: &
           incloud              ! True for in-cloud k values 


! Variables used in module:

  logical, private, save :: DEBUG_AQ = .false.

  real, private, save, dimension(KUPPER:KMAX_MID) :: &
        pr_acc                  ! Accumulated precipitation

  integer, private, save  :: kcloudtop   ! k-level of highest-cloud
  integer, private, save  :: ksubcloud   ! k-level just below cloud

  real, private, parameter :: & ! Define limits for "cloud"
      PR_LIMIT = 1.0e-7  &      ! for accumulated precipitation
     ,CW_LIMIT = 1.0e-10 &      ! for cloud water, kg(H2O)/kg(air)
     ,B_LIMIT  = 1.0e-3         ! for cloud cover (fraction)

  real, private, save :: &      ! Set in init below
      INV_Hplus          &      ! = 1.0/Hplus       (1/H+)
     ,INV_Hplus0p4              ! = INV_Hplus**0.4  (1/H+)**0.4

! The Henry's law coefficients, K, given in units of M or M atm-1,
! are calculated as effective. A factor K1fac = 1+K1/H+ is defined
! here also.

  integer, public, parameter :: &
  NHENRY  = 3, &  ! No. of species with Henry's law applied
  NK1     = 1, &  ! No. of species needing effective Henry's calc.
  IH_SO2  = 1, &
  IH_H2O2 = 2, &
  IH_O3   = 3

! Aqueous fractions:

  real, public, dimension ( NHENRY, KUPPER:KMAX_MID ), save  :: frac_aq
  real, private, dimension (NHENRY, CHEMTMIN:CHEMTMAX), save :: H
  real, private, dimension (NK1,CHEMTMIN:CHEMTMAX), save     :: K1fac

! Aqueous reaction rates for usage in gas-phase chemistry:

  integer, private, parameter :: NAQUEOUS = 4   ! No. aqueous rates
  integer, private, parameter :: NAQRC    = 3   ! No. constant rates

  real, public, dimension(NAQUEOUS,KCHEMTOP:KMAX_MID), save :: aqrck

  real, private, dimension(NAQRC), save :: aqrc ! constant rates for
                                                ! so2 oxidn.
  real, private, dimension(2), save :: vw       ! constant rates for
  logical, public,save :: prclouds_present      ! true if precipitating
                                                ! clouds

  integer, public, parameter :: &
      ICLOHSO2  = 1   &  ! for [oh] + [so2]
     ,ICLRC1    = 2   &  ! for [h2o2] + [so2]
     ,ICLRC2    = 3   &  ! for [o3] + [so2]
     ,ICLRC3    = 4      ! for [o3] + [o2] (Fe catalytic)


! Incloud scavenging: (only dependant on precipitation (not cloud water)
!-----------------------------------------------------------------------
! The parameterization of the scavenging of soluble chemical 
! components is scaled to the precipitation in each layer. For 
! incloud scavenging it is based on the parameterization described 
! in Berge (1998). The incloud scavenging of a soluble component 
! X is given by the expression:
!
!       X * f * pr_acc * W_sca
! Q =  ------------------------
!         Z_sca * rho_water
!
! where pr_acc is the accumulated precipitation in the layer, 
! Z_sca is the scavenging depth (scale=1000m) and rho_water is the
! density of water.
! f is an efficiency parameter. The module My_WetDep_ml should
! return the user-defined array WetDep%W_sca, with W_Sca representing
! the value of f.W_Sca/Z_sca/rho_water


! Sub cloud scavenging:
!-----------------------------------------------------------------------
! The sub-cloud scavenging distinguishes between particulate and
! gas-phase components. The scavenging of gases is calculated as
! Q = vw * P, where P is the accumulated precipitation (pr_acc, m)
! from all the layers above. (From setup_1d)
! For particles the scavenging is believed to be much less effective,
! as they follow the air-current around the droplets (Berge, 1993).
! Scavenging for particles is calculated as 
!    Q = A.e.P/v
! where A is 5.2 m3 kg-1 s-1, v is the fall-speed of the droplets,
! 5 m s-1, and e is the scavenging efficiency, 0.1.


contains

!-----------------------------------------------------------------------
subroutine Setup_Clouds(i,j)

!-----------------------------------------------------------------------
! DESCRIPTION
! Define incloud and precipitating clouds.
! The layer must contain at least 1.e-7 kgwater/kg air to 
! be considered a cloud. 
!
! Also calculates 
!  pr_acc - the accumulated precipitation for each k
!  b      - fractional cloud cover for each k
!-----------------------------------------------------------------------

  integer, intent(in) ::  i,j

  real, dimension(KUPPER:KMAX_MID) :: &
       b           &  !  Cloud-area (fraction)
      ,cloudwater     !  Cloud-water (volume mixing ratio) 
                      !  cloudwater = 1.e-6 same as 1.g m^-3

  integer :: k

! Add up the precipitation in the column:
  pr_acc(KUPPER) = sum ( pr(i,j,1:KUPPER) ) ! prec. from above 
  do k= KUPPER+1, KMAX_MID
     pr_acc(k) = pr_acc(k-1) + pr(i,j,k)
     pr_acc(k) = max( pr_acc(k), 0.0 )
  end do

  prclouds_present = .false.  
  if ( pr_acc(KMAX_MID) > PR_LIMIT )  prclouds_present = .true. 
                             ! --> precipitation at the surface

! initialise with .false. and 0:
  incloud(:)  = .false.
  cloudwater(:) = 0.


! Loop starting at surface finding the cloud base:

  ksubcloud = KMAX_MID+1       ! k-coordinate of sub-cloud limit

  do k = KMAX_MID, KUPPER, -1  
     if ( lwc(i,j,k) >  CW_LIMIT ) exit
     ksubcloud = k
  end do
  if ( ksubcloud == 0 ) return ! No cloud water found below level 6
                               ! Cloud above level 6 are likely thin
                               ! cirrus clouds, and if included may
                               ! need special treatment...
                               ! ==> assume no cloud 

! Define incloud part of the column requiring that both cloud water
! and cloud fractions are above limit values 

  kcloudtop = -1               ! k-level of cloud top
  do k = KUPPER, ksubcloud-1

     b(k) = cc3d(i,j,k)

! Units: kg(w)/kg(air) * kg(air(m^3) / density of water 10^3 kg/m^3
! ==> cloudwater (volume mixing ratio of water to air in cloud 
! (when devided by cloud fraction b )
! cloudwater(k) = 1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) / b(k)   

     if ( lwc(i,j,k) >  CW_LIMIT ) then
        cloudwater(k) = lwc(i,j,k) / b(k) ! value of cloudwater in the
                                          ! cloud fraction of the grid
        incloud(k) = .true.
        if ( kcloudtop < 0 ) kcloudtop = k
     end if

  end do

  if ( prclouds_present .and. kcloudtop == -1 ) then
     if ( DEBUG_AQ ) write(6,"(a20,i3,2i4,3es12.4)") &
        "ERROR prclouds sum_cw", &
         me, i,j, maxval(lwc(i,j,KUPPER:KMAX_MID),1) , &
         maxval(pr(i,j,:)), pr_acc(KMAX_MID)
     kcloudtop = KUPPER ! for safety
  end if

! sets up the aqueous phase reaction rates (SO2 oxidation) and the 
! fractional solubility 

  call setup_aqurates(b ,cloudwater,incloud)

  if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
     write(6,*) "DEBUG_AQ me presetn", me, prclouds_present
     write(6,*) "(a15,i3,2i4,es14.4)", "DEBUG_AQ me ", me, &
                kcloudtop, ksubcloud, pr_acc(KMAX_MID)
  end if

end subroutine Setup_Clouds


!-----------------------------------------------------------------------
subroutine init_aqueous() 

!-----------------------------------------------------------------------
! DESCRIPTION
! Calls initial tabulations, sets frac_aq to zero above cloud level, and
! sets constant rates.
! MTRLIM represents mass transport limitations between the clouds
! and the remainder of the grid-box volume. (so2 will be rapidly
! depleted within the clouds, and must be replenished from the
! surrounding cloudfree volume.
!-----------------------------------------------------------------------

  use PhysicalConstants_ml, only: AVOG ! Avogadro's No.

  real, parameter :: &
     Hplus = 5.0e-5                    ! H+ (Hydrogen ion concentration)
  real, parameter :: MASSTRLIM = 1.0   ! Mass transport limitation

  INV_Hplus    = 1.0/Hplus             ! 1/H+
  INV_Hplus0p4 = INV_Hplus**0.4        ! (1/H+)**0.4

! tabulations
!========================
  call tabulate_aqueous()
!========================

! Constant rates: The rates given in Berge (1993) are in mol-1 l.
! These need to be multiplied by 1.0e3/AVOG/Vf,so we perform the
! 1.0e3/AVOG scaling here.

! so2aq + h2o2   ---> so4, ref: Möller 1980
  aqrc(1) = 8.3e5 * 1.0e3/AVOG * MASSTRLIM

! (so2aq + hso3-) + H+ + o3 ---> so4, ref: Martin & Damschen 1981
  aqrc(2) = 1.8e4 * 1.0e3/AVOG * MASSTRLIM

! (so2aq + hso3-) + o2 ( + Fe ) --> so4, see documentation below
  aqrc(3) = 3.3e-10  * MASSTRLIM  

! Regarding aqrc(3):
! catalytic oxidation with Fe. The assumption is that 2% of SIV
! is oxidised per hour inside the droplets, corresponding to a
! conversion rate of 5.6^-6 (units s^-1  -- Therfore no conversion
! from mol l^-1)

! 5.6e-6 * 0.5e-6 (liquid water fraction) /8.5e-3 (fso2 at 10deg C)

! multiply with the assumed liquid water fraction from Seland and
! Iversen (0.5e-6) and with an assumed fso2 since the reaction is
! scaled by the calculated value for these parameters later.

end subroutine init_aqueous


!-----------------------------------------------------------------------
subroutine tabulate_aqueous()

!-----------------------------------------------------------------------
! DESCRIPTION
! Tabulates Henry's law coefficients over the temperature range
! defined in Tabulations_ml.
! For SO2, the effective Henry's law is given by
!   Heff = H * ( 1 + K1/H+ )  
! where k2 is omitted as it is significant only at high pH.
! We tabulate also the factor 1+K1/H+ as K1fac.
!-----------------------------------------------------------------------

  real, dimension(CHEMTMIN:CHEMTMAX) :: t, tfac ! Temperature, K, factor
  integer :: i

  t(:)           = (/ ( real(i), i=CHEMTMIN, CHEMTMAX ) /)
  tfac(:)        = 1.0/t(:) -  1.0/298.0

  H (IH_SO2 ,:)  = 1.23    * exp(3020.0*tfac(:) ) 
  H (IH_H2O2,:)  = 7.1e4   * exp(6800.0*tfac(:) )
  H (IH_O3  ,:)  = 1.13e-2 * exp(2300.0*tfac(:) )

! Need  effective Henry's coefficient for SO2:
  K1fac(IH_SO2  ,:)  =  &
     ( 1.0 + 1.23e-2 * exp(2010.0*tfac(:) ) * INV_Hplus)

  H (IH_SO2 ,:)  = H(IH_SO2,:) * K1fac(IH_SO2,:)

end subroutine tabulate_aqueous


!-----------------------------------------------------------------------
subroutine setup_aqurates(b ,cloudwater,incloud)

!-----------------------------------------------------------------------
! DESCRIPTION
! sets the rate-coefficients for thr aqueous-phase reactions
!-----------------------------------------------------------------------

  use Setup_1dfields_ml, only : &
     itemp          ! temperature (K)

  real, dimension(KUPPER:KMAX_MID) :: &
     b            & ! cloud-aread (fraction)
    ,cloudwater     ! cloud-water
  logical, dimension(KUPPER:KMAX_MID) :: &
           incloud  ! True for in-cloud k values 

! Outputs  -> aqurates

! local
  real, dimension(KUPPER:KMAX_MID) :: &
         fso2grid & ! f_aq * b = {f_aq}
        ,fso2aq   & ! only so2.h2o part (not hso4-)
        ,caqh2o2  & ! rate of oxidation of so2 with H2O2
        ,caqo3    & ! rate of oxidation of so2 with H2O2
        ,caqsx      ! rate of oxidation of so2 with o2 ( Fe )

  integer k

  call get_frac(cloudwater,incloud) ! => frac_aq
    
! initialize:
  aqrck(:,:)=0.


! Gas phase ox. of SO2 is "default"
! in cloudy air, only the part remaining in gas phase (not
! dissolved) is oxidized 

  aqrck(ICLOHSO2,:) = 1.0 

  do k = KUPPER,KMAX_MID
     if ( incloud(k) ) then ! Vf > 1.0e-10) ! lwc > CW_limit
        fso2grid(k) = b(k) * frac_aq(IH_SO2,k)
        fso2aq  (k) = fso2grid(k) / K1fac(IH_SO2,itemp(k))
        caqh2o2 (k) = aqrc(1) * frac_aq(IH_H2O2,k) / cloudwater(k)
        caqo3   (k) = aqrc(2) * frac_aq(IH_O3,k) / cloudwater(k)
        caqsx   (k) = aqrc(3) / cloudwater(k)  

! oh + so2 gas-phase
        aqrck(ICLOHSO2,k) = ( 1.0-fso2grid(k) ) ! now correction factor!
        aqrck(ICLRC1,k)   = caqh2o2(k) * fso2aq(k)
        aqrck(ICLRC2,k)   = caqo3(k) * INV_Hplus0p4 * fso2grid(k)
        aqrck(ICLRC3,k)   = caqsx(k) *  fso2grid(k) 

     end if
  enddo

end subroutine setup_aqurates


!-----------------------------------------------------------------------
subroutine get_frac(cloudwater,incloud)

!-----------------------------------------------------------------------
! DESCRIPTION
! Calculating pH dependant solubility fractions: Calculates the fraction
! of each soluble gas in the aqueous phase, frac_aq
!-----------------------------------------------------------------------

! intent in from used modules  : cloudwater and logical incloud
! intent out to rest of module : frac_aq

  use Setup_1dfields_ml, only : &
      temp       & ! temperature (K)
     ,itemp        ! temperature (K)
  use PhysicalConstants_ml, only: RGAS_ATML ! Gas-constant

! local
  real, dimension (KUPPER:KMAX_MID) :: &
      cloudwater   ! Volume fraction - see notes above.
  logical, dimension(KUPPER:KMAX_MID) :: &
      incloud      ! True for in-cloud k values
  real, dimension (KUPPER:KMAX_MID) :: VfRT ! Vf * Rgas * Temp

  integer :: ih, k ! index over species with Henry's law, vertical level k

! Make sure frac_aq is zero outside clouds:
  frac_aq(:,:) = 0.

  do k = KUPPER, KMAX_MID
     if ( incloud(k) ) then

        VfRT(k) = cloudwater(k) * RGAS_ATML * temp(k)

! Get aqueous fractions:
        do ih = 1, NHENRY
           frac_aq(ih,k) = 1.0 / ( 1.0+1.0/( H(ih,itemp(k))*VfRT(k) ) )
        end do

     end if
  end do

end subroutine get_frac


!-----------------------------------------------------------------------
subroutine WetDeposition(i,j)

!-----------------------------------------------------------------------
! DESCRIPTION
! Calculates wet deposition and changes in xn concentrations
! WetDeposition called from RunChem if precipitation reach the surface
!-----------------------------------------------------------------------

! input
  integer, intent(in) ::  i,j

! local
  integer :: itot        ! index in xn_2d arrays
  integer :: spec        ! species index from WetDep array
  integer :: n,k

  real    :: invgridarea ! xm2/(h*h)
  real    :: f_rho       ! Factors in rho calculation
  real    :: rho(KUPPER:KMAX_MID)
  real, dimension(KUPPER:KMAX_MID) :: vw ! Svavenging rates (tmp. array)

  real                    :: loss    ! conc. loss due to scavenging
  real,dimension(NWETDEP) :: sumloss ! sum conc. loss due to scavenging

! Loop starting from above:
  f_rho  = xmd(i,j)*(ps(i,j,1) - PT)/ATWAIR
  do k=kcloudtop, KMAX_MID           ! No need to go above cloudtop
     rho(k)  = f_rho * carea(k) / amk(k)
  end do
  sumloss(:) = 0.0

  invgridarea = xm2(i,j)/( gridwidth_m*gridwidth_m )

! calculate concentration after wet deposition and sum up the vertical
! column of the depositions for the fully soluble species.

  if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
     Write(6,*) "(a15,i3,2i4,es14.4)", "DEBUG_WDEP2", me, &
                   kcloudtop, ksubcloud, pr_acc(KMAX_MID)
  end if ! DEBUG

  do spec = 1, NWETDEP

! Put both in- and sub-cloud scavenging ratios in the array vw:

     vw(kcloudtop:ksubcloud-1) = WetDep(spec)%W_sca ! Scav. for incloud
     vw(ksubcloud:KMAX_MID  )  = WetDep(spec)%W_sub ! Scav. for subcloud

     if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
        Write(6,*) "(a15,i3,2es14.4)", "DEBUG_WDEP Sca", &
                   spec, WetDep(spec)%W_sca, WetDep(spec)%W_sub
     end if ! DEBUG

     itot = WetDep(spec)%itot

     do k = kcloudtop, KMAX_MID
        loss = xn_2d(itot,k) * ( 1.0 - exp( -vw(k)*pr_acc(k)*dt ) )
        xn_2d(itot,k) = xn_2d(itot,k) - loss
        sumloss(spec) = sumloss(spec) + loss * rho(k)

        if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
           Write(6,*)"(a30,4i4)", "DEBUG_WDEP me, k, itot, spec", &
                     me, k, itot, spec
           Write(6,*)"(a30,i4, 2es12.4, f8.2)", &
                     "DEBUG_WDEP me, vw, pr, dt ", me,  &
                     vw(k), pr_acc(k), dt
           Write(6,*)"(a30,3es12.4)", "DEBUG_WDEP loss, xn, sumloss", &
                     loss, xn_2d(itot,k), sumloss(spec)
        end if ! DEBUG
     end do ! k loop
  end do ! spec loop

  d_2d(WDEP_PREC,i,j,IOU_INST) = sum ( pr(i,j,:) ) * dt
                                              ! Same for all models

! add other losses into twetdep and wdep arrays:
  call WetDep_Budget(i,j,sumloss,invgridarea) ! Model-specific

end subroutine WetDeposition

!-----------------------------------------------------------------------
end module Aqueous_ml
