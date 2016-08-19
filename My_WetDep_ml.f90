! <My_WetDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module My_WetDep_ml
 use MassBudget_ml,     only : totwdep
 use ModelConstants_ml, only : atwS, atwN, atwPM
 use Derived_ml,  only : f_2d, d_2d, IOU_INST
 use SmallUtils_ml,  only : find_index


  use GenSpec_tot_ml          ! SO2, SO4, etc.
  use GenSpec_adv_ml          ! IXADV_SO2, IXADV_SO4, etc.
  implicit none
  private

  public :: Init_WetDep       ! Call from Unimod
  public :: WetDep_Budget     ! Call from Aqueous_ml

  type, public :: WScav
     integer  :: itot     !ds may05 - was adv - confusing
     real  :: W_sca       ! Scavenging ratio/z_Sca/rho = W_sca/1.0e6
     real  :: W_sub       ! same for subcloud
  end type WScav
  

  integer, public, parameter :: NWETDEP =  14  !SeaS 11  ! Number of solublity classes
  type(WScav), public, dimension(NWETDEP), save  :: WetDep
  
  integer, public, save  :: WDEP_PREC   ! Used in Aqueous_ml
  integer, private, save :: WDEP_SOX, WDEP_OXN, WDEP_RDN

contains

  subroutine Init_WetDep()

  !/ INCLOUDFAC is A/v where A is 5.2 m3 kg-1 s-1, !  and v is the fallspeed (5 m/s). 
    real, parameter ::  FALLSPEED = 5.0                            ! m/s 
    real, parameter ::  SUBCLFAC = 5.2 / FALLSPEED

  !/ e is the scavenging efficiency (0.1 for fine particles, 0.4 for course)

    real, parameter ::  EFF25 = 0.1*SUBCLFAC  & 
                      , EFFCO = 0.4*SUBCLFAC  ! collection efficiency b/clouds - coarse

   !/.. setup the scavenging ratios for in-cloud and sub-cloud. For
   !    gases, sub-cloud = 0.5 * incloud. For particles, sub-cloud=
   !    efficiency * INCLOUDFAC. See also notes in Aqueous_ml.

   !/..                        W_Sca  W_sub
  
    WetDep(1)   = WScav(SO2,    0.3,  0.15)   ! Berge+Jakobsen, issh
    WetDep(2)   = WScav(SO4,    1.0,  EFF25)  ! Berge+Jakobsen, jej
    WetDep(3)   = WScav(aNH4,   1.0,  EFF25)
    WetDep(4)   = WScav(NH3,    1.4,  0.5 )  ! subcloud = 1/3 of cloud for gases
    WetDep(5)   = WScav(aNO3,   1.0,  EFF25)  
    WetDep(6)   = WScav(HNO3,   1.4,  0.5)   ! 
    WetDep(7)   = WScav(H2O2,   1.4,  0.5)   ! 
    WetDep(8)   = WScav(HCHO,   0.1,  0.03)  ! 
    WetDep(9)   = WScav(pNO3,   1.0,  EFFCO) !!
    WetDep(10)  = WScav(PM25,   1.0,  EFF25)
    WetDep(11)  = WScav(PMCO,   1.0,  EFFCO)
    WetDep(12)  = WScav(SSFI,   1.0,  EFF25)   !SeaS
    WetDep(13)  = WScav(SSCO,   1.0,  EFFCO)   !SeaS
    WetDep(14)  = WScav(Pb210,  1.0,  EFF25)   !

   !####################### ds NEW define indices here ##########

     WDEP_PREC= find_index("WDEP_PREC",f_2d(:)%name)
     WDEP_SOX = find_index("WDEP_SOX",f_2d(:)%name)
     WDEP_OXN = find_index("WDEP_OXN",f_2d(:)%name)
     WDEP_RDN = find_index("WDEP_RDN",f_2d(:)%name)
   !####################### ds END define indices here ##########

  end subroutine Init_WetDep

  subroutine WetDep_Budget(i,j,sumloss,invgridarea)
     integer,            intent(in) ::  i,j
     real, dimension(:), intent(in) :: sumloss
     real, intent(in)  :: invgridarea
     real :: wdeps, wdepox, wdepred, wdeppm25, wdeppmco


      !wdeps = sumloss(SO2) + sumloss(SO4)
       wdeps = sumloss(1) + sumloss(2)

      !wdepred = sumloss(NH3)  + sumloss(NH4) & 
       wdepred = sumloss(4)  + sumloss(3) !
  
      !wdepox  = sumloss(HNO3) + sumloss(aNO3) + pNO3
       wdepox  = sumloss(6) + sumloss(5) + sumloss(9)
      wdeppm25= sumloss(7) 
      wdeppmco= sumloss(8) 

       totwdep(IXADV_SO4)  = totwdep(IXADV_SO4) + wdeps
       totwdep(IXADV_HNO3) = totwdep(IXADV_HNO3) + wdepox
       totwdep(IXADV_NH3)  = totwdep(IXADV_NH3)  + wdepred
       totwdep(IXADV_PM25)  = totwdep(IXADV_PM25)  + wdeppm25
       totwdep(IXADV_PMco)  = totwdep(IXADV_PMco)  + wdeppmco


       d_2d(WDEP_SOX,i,j,IOU_INST) = wdeps * atwS * invgridarea 
       d_2d(WDEP_OXN,i,j,IOU_INST) = wdepox * atwN * invgridarea 
       d_2d(WDEP_RDN,i,j,IOU_INST) = wdepred * atwN * invgridarea 


  end subroutine WetDep_Budget
end module My_WetDep_ml
