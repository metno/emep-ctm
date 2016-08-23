! <Output_hourly.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
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
!***********************************************************************
subroutine hourly_out() !!  spec,ofmt,ix1,ix2,iy1,iy2,unitfac)
!***********************************************************************
!**    DESCRIPTION:
!       Calculates and
!       Outputs hourly concentration (or met) values for a sub-set of the grid.
!
!**    REVISION HISTORY:
!      Extended to produce new file, Hourly.mmyy, every month, 10/5/01 ds
!      stop_test used instead of stop_all, su, 05/01
!      Extended for variable format, met, xn_adv or xn_shl, ds, and to use
!       Asc2D type 19/4/01
!      Corrected for IRUNBEG, etc., su, 4/01
!      New, ds, 5/3/99
!
!*************************************************************************
!
  use My_Outputs_ml,    only: NHOURLY_OUT,    & ! No. outputs
                              NLEVELS_HOURLY, & ! No. output levels
                              hr_out,         & ! Required outputs
                              LEVELS_HOURLY ! Output selected model levels
            !NML SELECT_LEVELS_HOURLY, LEVELS_HOURLY ! Output selected model levels

  use CheckStop_ml,     only: CheckStop
  use Chemfields_ml,    only: xn_adv,xn_shl,cfac,PM25_water,PM25_water_rh50,AOD
  use ChemGroups_ml,    only: chemgroups
  use Derived_ml,       only: num_deriv2d        ! D2D houtly output type
  use DerivedFields_ml, only: f_2d,d_2d          ! D2D houtly output type
  use OwnDataTypes_ml,  only: Asc2D, Deriv
  use ChemSpecs_shl_ml ,only: NSPEC_SHL          ! Maps indices
  use ChemChemicals_ml ,only: species            ! Gives names
  use GridValues_ml,    only: i_fdom, j_fdom,&   ! Gives emep coordinates
                              debug_proc, debug_li,debug_lj
  use Io_ml,            only: IO_HOURLY
  use ModelConstants_ml,only: KMAX_MID, MasterProc, &
                              IOU_INST, IOU_HOUR, IOU_YEAR, IOU_HOUR_PREVIOUS, &
                              DEBUG => DEBUG_OUT_HOUR,runlabel1,HOURLYFILE_ending,&
                              FORECAST
  use ModelConstants_ml,only: SELECT_LEVELS_HOURLY !NML
  use ModelConstants_ml,only: MFAC! TESTSHL
  use MetFields_ml,     only: t2_nwp,th, roa, surface_precip, ws_10m ,rh2m,&
                              pzpbl, ustar_nwp, Kz_m2s, &
                              Idirect, Idiffuse, z_bnd, z_mid,ps
  use NetCDF_ml,        only: Out_netCDF, CloseNetCDF, Init_new_netCDF, fileName_hour, &
                              Int1, Int2, Int4, Real4, Real8  !Output data type to choose
  use OwnDataTypes_ml,  only: TXTLEN_DERIV,TXTLEN_SHORT
  use Par_ml,           only: MAXLIMAX, MAXLJMAX, GIMAX,GJMAX,        &
                              me, IRUNBEG, JRUNBEG, limax, ljmax
! FUTURE  use Pollen_ml,        only: heatsum, pollen_left, AreaPOLL
  use TimeDate_ml,      only: current_date
  use TimeDate_ExtraUtil_ml,only : date2string
  use Units_ml,         only: Group_Units

  implicit none

!*.. Components of  hr_out
!* character(len=TXTLEN_DERIV):: name   ! netCDF variable name
!* character(len=TXTLEN_SHORT):: type   ! "ADVppbv" or "ADVugm3" or "SHLmcm3"
!* character(len=9) :: ofmt     ! Output format (e.g. es12.4)
!* integer          :: spec     ! Species number in xn_adv or xn_shl array or other arrays
!* integer          :: ix1,ix2  ! bottom-left & upper-right x
!* integer          :: iy1,iy2  ! bottom-left & upper-right y
!* integer          :: nk       ! number of vertical levels
!* character(len=TXTLEN_SHORT) :: unit   ! netCDF unit attribute
!* real             :: unitconv !  conv. factor
!* real             :: max      ! Max allowed value for output

! local variables
  logical, save     :: my_first_call = .true. ! Set false after file opened
  integer msnr                        ! Message number for rsend
  real hourly(MAXLIMAX,MAXLJMAX)      ! Local hourly value  (e.g. ppb)
  real ghourly(GIMAX,GJMAX)           ! Global hourly value (e.g. ppb)
  real :: arrmax                      ! Maximum value from array
  real :: unit_conv                   ! Unit conversion (ppb ug etc.)
  real :: g                           ! tmp - saves value of ghourly(i,j)
  integer, dimension(2) :: maxpos     ! Location of max value
  integer i,j,ih,ispec,itot,iadv      ! indices
  integer :: k,ik,iik                 ! Index for vertical level
  integer ist,ien,jst,jen             ! start and end coords
  character(len=TXTLEN_DERIV) :: name ! For output file, species names
  character(len=4)  :: suffix         ! For date "mmyy"
  integer, save :: prev_month = -99   ! Initialise with non-possible month
  type(Deriv) :: def1           ! for NetCDF
  real :: scale                 ! for NetCDF
  integer ::CDFtype,nk,klevel   ! for NetCDF
  character(len=TXTLEN_SHORT)    :: hr_out_type=""      ! hr_out%type
  integer                        :: hr_out_nk=0         ! hr_out%nk
  integer, pointer, dimension(:) :: gspec=>null()       ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null()  ! & unit conv. factors

  character(len=len(fileName_hour)) :: filename
  logical, save :: debug_flag      ! = ( MasterProc .and. DEBUG )
  logical       :: surf_corrected  ! to get 3m values

  logical       :: file_exist=.false.

  if(NHOURLY_OUT<= 0) then
    if(my_first_call.and.MasterProc.and.DEBUG) &
       write(*,*) "DEBUG Hourly_out: nothing to output!"
    my_first_call = .false.
    return
  endif

  if(my_first_call) then
    debug_flag=(debug_proc.and.DEBUG)

  !/ Ensure that domain limits specified in My_Outputs lie within
  !  model domain. In emep coordinates we have:
    do ih = 1, NHOURLY_OUT
      hr_out(ih)%ix1 = max(IRUNBEG,hr_out(ih)%ix1)
      hr_out(ih)%iy1 = max(JRUNBEG,hr_out(ih)%iy1)
      hr_out(ih)%ix2 = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
      hr_out(ih)%iy2 = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)
!     hr_out(ih)%nk  = min(KMAX_MID,hr_out(ih)%nk)
      if(debug_flag) write(*,*) "DEBUG Hourly nk ", ih, hr_out(ih)%nk
    enddo ! ih
  endif  ! first_call

  filename=trim(runlabel1)//date2string(HOURLYFILE_ending,current_date)
  inquire(file=filename,exist=file_exist)
  if(my_first_call.or..not.file_exist)then
    if(debug_flag) write(*,*) "DEBUG ",HOURLYFILE_ending,"-Hourlyfile ", trim(filename)
    call Init_new_netCDF(trim(filename),IOU_HOUR)

  !! Create variables first, without writing them (for performance purposes)   
    do ih=1,NHOURLY_OUT
      def1%name=hr_out(ih)%name
      def1%unit=hr_out(ih)%unit
      def1%class=hr_out(ih)%type
      ist = hr_out(ih)%ix1
      jst = hr_out(ih)%iy1
      ien = hr_out(ih)%ix2
      jen = hr_out(ih)%iy2
      nk  = hr_out(ih)%nk
      CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
      scale=1.
      if(any(hr_out(ih)%type==(/"ADVppbv     ","ADVugXX     ","ADVugXXgroup",&
                                "COLUMN      ","COLUMNgroup ","D2D         "/)))nk=1
      select case(nk)
      case(1)       ! write as 2D
        call Out_netCDF(IOU_HOUR,def1,2,1,hourly,scale,CDFtype,ist,jst,ien,jen,&
          create_var_only=.true.)
      case(2:)      ! write as 3D
        call Out_netCDF(IOU_HOUR,def1,3,1,hourly,scale,CDFtype,ist,jst,ien,jen,1,&
          create_var_only=.true.,chunksizes=(/ien-ist+1,jen-jst+1,1,1/))
      endselect
    enddo
  endif

  my_first_call = .false.
!......... Uses concentration/met arrays from Chem_ml or Met_ml ..................
!
!        real xn_adv(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)
!        real cfac(NSPEC_ADV,MAXLIMAX,MAXLJMAX)
! or...
!        real xn_shl(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)
! or...
!        real temp2m(MAXLIMAX,MAXLJMAX)
!
!..........................................................................

  hourly(:,:) = 0.0 ! initialize

  HLOOP: do ih = 1, NHOURLY_OUT

    hr_out_type=trim(hr_out(ih)%type)
    hr_out_nk=hr_out(ih)%nk
    if(any(hr_out_type==(/"ADVppbv     ","ADVugXX     ","ADVugXXgroup",&
                          "COLUMN      ","COLUMNgroup ","D2D         "/)))hr_out_nk=1

    KVLOOP: do k = 1,hr_out_nk

      msnr  = 3475 + ih
      ispec = hr_out(ih)%spec
      name  = hr_out(ih)%name
      if(debug_flag) &
        write(*,'(a,2i4,1X,a,/a,1X,2a,i3)')"DEBUG DERIV HOURLY", ih, ispec, &
          trim(name),"INTO HOUR TYPE:", &
          trim(hr_out(ih)%type) // " "//trim(hr_out(ih)%name), " nk:", hr_out_nk

      if(any(hr_out_type==(/"COLUMN     " ,"COLUMNgroup"/)))then
        ik=KMAX_MID-hr_out(ih)%nk+1  ! top of the column
        if(ik>=KMAX_MID)ik=1         ! 1-level column does not make sense
      else
        ik=KMAX_MID-k+1              ! all levels from model bottom are outputed,
        if(debug_flag) write(*,*)"SELECT LEVELS? ", ik, SELECT_LEVELS_HOURLY
        if(SELECT_LEVELS_HOURLY)then ! or the output levels are taken
          ik=LEVELS_HOURLY(k)        ! from LEVELS_HOURLY array (default)
          hr_out_type=hr_out(ih)%type
          if(debug_flag) write(*,*)"DEBUG SELECT LEVELS", ik, hr_out_type
          surf_corrected = (ik==0)  ! Will implement cfac
          if(debug_flag.and.surf_corrected) &
            write(*,*)"DEBUG HOURLY Surf_correction", ik, k
!TESTHH QUERY: see below
          if(ik==0)then
            ik=KMAX_MID              ! surface/lowermost level
            if(debug_flag) write(*,*)"DEBUG LOWEST LEVELS", ik, hr_out_type
            if(any(hr_out_type==(/"BCVppbv     ","BCVugXX     ",&
                                  "BCVugXXgroup"/)))&
!                                 "BCVugXXgroup","Out3D       "/)))&
              hr_out_type(1:3)="ADV" ! ensure surface output
            if(any(hr_out_type==(/"PMwater"/)))&
              hr_out_type=trim(hr_out_type)//"SRF"
          else
!TESTHH QUERY:
            ik=KMAX_MID-ik+1         ! model level to be outputed
            if(any(hr_out_type==(/"ADVppbv     ","ADVugXX     ","ADVugXXgroup"/)))&
              ik=KMAX_MID            ! all ADV* types represent surface output
          endif
        endif
      endif

   !----------------------------------------------------------------
   ! Multi-layer output.
   !  Specify NLEVELS_HOURLY here, and in hr_out defs use either:
   !    ADVppbv to get surface concentrations (onyl relevant for
   !            layer k=20 of course - gives meaningless number
   !            for higher levels).
   !  Or,
   !    BCVppbv to get grid-centre concentrations (relevant for all layers).
   !
   !  For ug output (ug/m3, ugC/m3, ugS/m3, ugN/m3) use
   !    ADVugXX to get surface concentrations (only lowermost model level).
   !  Or,
   !    BCVppbv to get grid-centre concentrations (model levels).
   !  For ug/m3     output, set hr_out%unitconv=to_ug_ADV(ixadv).
   !  For ugX/m3    output, set hr_out%unitconv=to_ug_X(ixadv).
   !
   !  For Column integrated output use COLUMN
   !  For ug/m2     output, set hr_out%unitconv=to_ug_ADV(ixadv).
   !  For ugX/m2    output, set hr_out%unitconv=to_ug_X(ixadv).
   !  For molec/cm2 output, set hr_out%unitconv=to_molec_cm2.
   !----------------------------------------------------------------

      if(debug_flag) write(*,"(5a,i4)") "DEBUG Hourly MULTI ",&
              trim(hr_out(ih)%name), " case ", trim(hr_out_type), " k: ",  ik
      OPTIONS: select case(hr_out_type)

      case("ADVppbv")
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                      * cfac(ispec,i,j) &    ! 50m->3m conversion
                      * unit_conv            ! Units conv.
        endforall

      case("Out3D")
        itot = ispec
        iadv  = ispec - NSPEC_SHL
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        if(itot < NSPEC_SHL .and.  debug_proc) write(*,"(a,a,4es12.3)") &
          "OUT3D MAR22 ", trim(name), unit_conv, &
           MFAC,roa(2,2,ik,1)*MFAC, xn_shl(ispec,2,2,ik)

        !MAR22 Added SHL option.
        if( itot <= NSPEC_SHL ) then !TESTSHL
         ! CRUDE units fix. Sort out later.
         !   Inverse of No. air mols/cm3 = 1/M
         ! where M =  roa (kgair m-3) * MFAC  when ! scale in ug,  else 1
         !inv_air_density3D(i,j,k) = 1.0/( roa(i,j,k,1) * MFAC )

          forall(i=1:limax,j=1:ljmax)
             hourly(i,j) = xn_shl(ispec,i,j,ik)
          end forall
          if(index(hr_out(ih)%unit,"ppt")>0) then
            forall(i=1:limax,j=1:ljmax)
                hourly(i,j) = hourly(i,j)/( roa(i,j,ik,1) * MFAC ) * 1.0e9
            end forall
          else
            call CheckStop("SHL Out3D option not coded yet")
          end if
        
        else  ! ADV:, original code
  
          if(index(hr_out(ih)%unit,"ppb")>0) then
            forall(i=1:limax,j=1:ljmax)
              hourly(i,j) = xn_adv(iadv ,i,j,ik) &
                          * unit_conv            ! Units conv.
            end forall
          else if(index(hr_out(ih)%unit,"ug")>0) then
            forall(i=1:limax,j=1:ljmax)
              hourly(i,j) = xn_adv(iadv ,i,j,ik) &
                        * roa(i,j,ik,1)        & ! density.
                        * unit_conv              ! Units conv.
            end forall
          else
            call CheckStop("ERROR: Output_hourly  unit problem"//trim(name) )
          endif
        end if ! MAR22 ADV/SHL split 

        if(surf_corrected.and.ik==KMAX_MID.and.itot>NSPEC_SHL) then
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = hourly(i,j)*cfac(iadv,i,j) ! 50m->3m conversion
          endforall
        endif

        if(debug_flag) then
          i=debug_li; j=debug_lj
          write(*,'(A,2I4,1X,L2,2f10.4)')"Out3D K-level"//trim(name), ik,  &
            itot, surf_corrected, hourly(i,j), cfac(ispec-NSPEC_SHL,i,j)
        endif

      case("BCVppbv")
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall(i=1:limax,j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,ik) &
                      * unit_conv            ! Units conv.
        endforall
        if(DEBUG) &
          write(*,'(A,I0,1X,L2)')"K-level", ik, trim(name), itot, surf_corrected

      case("ADVugXX")  !ug/m3, ugX/m3 output at the surface
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                      * cfac(ispec,i,j) &     ! 50m->3m conversion
                      * unit_conv       &     ! Units conv.
                      * roa(i,j,KMAX_MID,1)   ! density.
        endforall

      case("BCVugXX")  ! ug/m3, ugX/m3 output at model mid-levels
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_adv(ispec,i,j,ik) &
                !BCV  * cfac(ispec,i,j) &     ! 50m->3m conversion
                      * unit_conv       &     ! Units conv.
                      * roa(i,j,ik,1)         ! density.
        endforall
        if(DEBUG) &
          write(*,'(a,i5,a10,i5)')"K-level", ik, trim(name), itot

      case("ADVugXXgroup")  ! GROUP output in ug/m3, ugX/m3 at the surface
        call Group_Units(hr_out(ih),gspec,gunit_conv,debug_flag,name)
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = dot_product(xn_adv(gspec,i,j,KMAX_MID),&
                                    cfac(gspec,i,j) & ! 50m->3m conversion
                                   *gunit_conv(:))  & ! Units conv.
                      * roa(i,j,KMAX_MID,1)           ! density.
        endforall
        if(DEBUG) &
          write(*,'(2a10,99i5)')"Surface", trim(name), gspec+NSPEC_SHL
        deallocate(gspec,gunit_conv)

      case("BCVugXXgroup")  ! GROUP output in ug/m3, ugX/m3 at model mid-levels
        call Group_Units(hr_out(ih),gspec,gunit_conv,debug_flag,name)
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = dot_product(xn_adv(gspec,i,j,ik), &
                                    gunit_conv(:))  & ! Units conv.
                      * roa(i,j,ik,1)                 ! density.
        endforall
        if(DEBUG) &
          write(*,'(a10,i7,a10,i7)')"K-level", ik, trim(name), gspec+NSPEC_SHL
        deallocate(gspec,gunit_conv)

      case("PMwater")    ! PM water content in ug/m3 at model mid-levels
        if(hr_out(ih)%unit/="ug/m3")hr_out(ih)%unit="ug"
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = PM25_water(i,j,ik)

      case("PMwaterSRF")  ! PM water content in ug/m3 at surface level
        if(hr_out(ih)%unit/="ug/m3")hr_out(ih)%unit="ug"
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = PM25_water_rh50(i,j)

      case("Z","Z_MID")
        name = "Z_MID"
        unit_conv =  hr_out(ih)%unitconv
        if(surf_corrected)then
          forall(i=1:limax,j=1:ljmax) hourly(i,j) = 0.0
        else
          forall(i=1:limax,j=1:ljmax) hourly(i,j) = z_mid(i,j,ik)*unit_conv
        endif

      case("AOD")
        name = "AOD 550nm"
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = AOD(i,j)

      case("COLUMN")    ! Column output in ug/m2, ugX/m2, molec/cm2
        itot = NSPEC_SHL + ispec
        name = species(itot)%name
        unit_conv =  hr_out(ih)%unitconv
        if(ih>0) hourly(:,:) = 0.0
        do iik=ik,KMAX_MID
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = hourly(i,j)                     &
                      + xn_adv(ispec,i,j,iik)             &
                      *                  unit_conv        & ! Units conv.
                      * roa(i,j,iik,1)                    & ! density.
                      * (z_bnd(i,j,iik)-z_bnd(i,j,iik+1))   ! level thickness
          endforall
        enddo

      case("COLUMNgroup")! GROUP Column output in ug/m2, ugX/m2, molec/cm2
        call Group_Units(hr_out(ih),gspec,gunit_conv,debug_flag,name)
        if(ih>1) hourly(:,:) = 0.0
        do iik=ik,KMAX_MID
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = hourly(i,j)                     &
                      + dot_product(xn_adv(gspec,i,j,iik),&
                                    gunit_conv(:))        & ! Units conv.
                      * roa(i,j,iik,1)                    & ! density.
                      * (z_bnd(i,j,iik)-z_bnd(i,j,iik+1))   ! level thickness
          endforall
        enddo
        if(DEBUG) &
          write(*,'(a10,i7,a10,i7)')"K-level", ik, trim(name), gspec+NSPEC_SHL
        deallocate(gspec,gunit_conv)

      case("SHLmcm3")        ! No cfac for short-lived species
        itot = ispec
        name = species(itot)%name
        forall ( i=1:limax, j=1:ljmax)
          hourly(i,j) = xn_shl(ispec,i,j,KMAX_MID) &
                      * hr_out(ih)%unitconv  ! Units conv.
        endforall

      case("T2_C")        ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = t2_nwp(i,j,1) - 273.15

      case("rh2m")        ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = rh2m(i,j,1)*100

      case("ws_10m")      ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = ws_10m(i,j,1)
   
      case("theta")       ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = th(i,j,KMAX_MID,1)

      case("PRECIP")      ! No cfac for surf.variable; Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = surface_precip(i,j)

      case("Idirect")     ! Direct radiation (W/m2); Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = Idirect(i,j)

      case("Idiffus")     ! Diffuse radiation (W/m2); Skip Units conv.
        forall(i=1:limax,j=1:ljmax) hourly(i,j) = Idiffuse(i,j)

      case ( "D2D" )
       ! Here ispec is the index in the f_2d arrays
        call CheckStop(ispec<1.or.ispec>num_deriv2d,&
          "ERROR-DEF! Hourly_out: "//trim(hr_out(ih)%name)//", wrong D2D id!")
        if(hr_out(ih)%unit=="") hr_out(ih)%unit = f_2d(ispec)%unit
        unit_conv = hr_out(ih)%unitconv*f_2d(ispec)%scale
        if(f_2d(ispec)%avg)then           ! non accumulated variables
          if( debug_flag ) write(*,*) " D2Davg ",&
          trim(hr_out(ih)%name), ih, ispec, trim(f_2d(ispec)%name), f_2d(ispec)%avg
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = d_2d(ispec,i,j,IOU_INST) * unit_conv
          endforall
        else                              ! hourly accumulated variables
          if(debug_flag) then
            i=debug_li
            j=debug_lj
            write(*,"(2a,2i4,a,3g12.3)") "OUTHOUR D2Dpre ",&
              trim(hr_out(ih)%name), ih, ispec,trim(f_2d(ispec)%name),&
              d_2d(ispec,i,j,IOU_YEAR), d_2d(ispec,i,j,IOU_HOUR_PREVIOUS),&
              unit_conv
          endif
  
          forall(i=1:limax,j=1:ljmax)
            hourly(i,j) = (d_2d(ispec,i,j,IOU_YEAR)&
                          -d_2d(ispec,i,j,IOU_HOUR_PREVIOUS)) * unit_conv
            d_2d(ispec,i,j,IOU_HOUR_PREVIOUS)=d_2d(ispec,i,j,IOU_YEAR)
          endforall
        endif
        if(debug_flag) &
          write(*,'(a,2i3,2es12.3)')"HHH DEBUG D2D", ispec, ih, &
            hr_out(ih)%unitconv, hourly(debug_li,debug_lj)

      case DEFAULT
        call CheckStop( "ERROR-DEF! Hourly_out: '"//trim(hr_out(ih)%type)//&
                        "' hourly type not found!")

      endselect OPTIONS

      if(debug_flag) &
        write(*,"(a,3i4,2g12.3)")"DEBUG-HOURLY-OUT:"//trim(hr_out(ih)%name),&
           me,ih,ispec, hourly(debug_li,debug_lj), hr_out(ih)%unitconv

      !/ Get maximum value of hourly array
      hourly(limax+1:MAXLIMAX,:) = 0.0
      hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0
      arrmax = maxval(hourly)

      if((hr_out(ih)%max>0.0).and.(arrmax>hr_out(ih)%max)) then
        write(*,*) "Hourly value too big!: ", ih, trim(hr_out(ih)%type), arrmax
        write(*,*) "Species : ", trim(name)," : ",  " ispec ", ispec
        write(*,*) "max allowed is : ",  hr_out(ih)%max
        write(*,*) "unitconv was   : ", hr_out(ih)%unitconv
        write(*,*) " me, limax, ljmax, MAXLIMAX,MAXLJMAX : ",  me, &
                          limax, ljmax ,MAXLIMAX,MAXLJMAX
        maxpos = maxloc(hourly)
        write(*,*) "Location is i=", maxpos(1), " j=", maxpos(2)
        write(*,*) "EMEP coords ix=", i_fdom(maxpos(1)), " iy=", j_fdom(maxpos(2))
        write(*,*) "hourly is ", hourly(maxpos(1),maxpos(2))
        if(hr_out(ih)%type(1:3)=="ADV") then
          write(*,*) "xn_ADV is ", xn_adv(ispec,maxpos(1),maxpos(2),KMAX_MID)
          write(*,*) "cfac   is ",   cfac(ispec,maxpos(1),maxpos(2))
        endif
        call CheckStop("Error, Output_hourly/hourly_out: too big!")
      endif

!NetCDF hourly output
      def1%name=hr_out(ih)%name
      def1%unit=hr_out(ih)%unit
      def1%class=hr_out(ih)%type
      ist = max(IRUNBEG,hr_out(ih)%ix1)
      jst = max(JRUNBEG,hr_out(ih)%iy1)
      ien = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
      jen = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)
      nk  = min(KMAX_MID,hr_out_nk)
      CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
      scale=1.

      select case(nk)
      case(1)       ! write as 2D
        call Out_netCDF(IOU_HOUR,def1,2,1,hourly,scale,CDFtype,ist,jst,ien,jen)
      case(2:)      ! write as 3D
        klevel=ik
        if(nk<KMAX_MID)  klevel=KMAX_MID-ik+1 !count from ground and up
        if(SELECT_LEVELS_HOURLY) klevel=k     !order is defined in LEVELS_HOURLY
        call Out_netCDF(IOU_HOUR,def1,3,1,hourly,scale,CDFtype,ist,jst,ien,jen,klevel)
      !case default   ! no output
      endselect
    enddo KVLOOP
    
  enddo HLOOP
  !write out surface pressure at each time step to define vertical coordinates
  if(NHOURLY_OUT>0)then
      def1%name='PS'
      def1%unit='hPa'
      def1%class='Surface pressure'
      CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
      scale=1.
      call Out_netCDF(IOU_HOUR,def1,2,1,ps(:,:,1)*0.01,scale,CDFtype,ist,jst,ien,jen)     
  endif

!Not closing seems to give a segmentation fault when opening the daily file
!Probably just a bug in the netcdf4/hdf5 library.
  call CloseNetCDF

! Write text file wich contents
  if(.not.(FORECAST.and.MasterProc))return
  i=index(filename,'.nc')-1;if(i<1)i=len_trim(filename)
  open(IO_HOURLY,file=filename(1:i)//'.msg',position='append')
  write(IO_HOURLY,*)date2string('FFF: YYYY-MM-DD hh',current_date)
  close(IO_HOURLY)
endsubroutine hourly_out
