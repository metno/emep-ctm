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
  use My_Outputs_ml,    only: NHOURLY_OUT, &     ! No. outputs
                              NLEVELS_HOURLY, &  ! No. output levels
                              FREQ_HOURLY, &     ! No. hours between outputs
                              Asc2D, hr_out, &   ! Required outputs
                              Hourly_ASCII, &    ! ASCII output or not
                              to_ug_ADV, to_ug_C, to_ug_S, to_ug_N, &
                              SELECT_LEVELS_HOURLY, LEVELS_HOURLY !Output selected model levels
  use CheckStop_ml,     only: CheckStop
  use Chemfields_ml,    only: xn_adv,xn_shl, cfac, PM25_water_rh50
  use ChemGroups_ml,    only: chemgroups
  use Derived_ml,       only: num_deriv2d        ! D2D houtly output type
  use DerivedFields_ml, only: f_2d,d_2d          ! D2D houtly output type
  use OwnDataTypes_ml,  only: Deriv
  use ChemSpecs_shl_ml ,only: NSPEC_SHL          ! Maps indices
  use ChemChemicals_ml ,only: species            ! Gives names
  use GridValues_ml,    only: i_fdom, j_fdom     ! Gives emep coordinates
  use Io_ml,            only: IO_HOURLY
  use ModelConstants_ml,only: KMAX_MID, DEBUG_i, DEBUG_j, MasterProc, &
                              IOU_INST, IOU_HOUR, IOU_YEAR, IOU_HOUR_PREVIOUS
  use MetFields_ml,     only: t2_nwp,th, roa, surface_precip,         &
                              Idirect, Idiffuse, z_bnd
  use NetCDF_ml,        only: Out_netCDF,                             &
                              Int1, Int2, Int4, Real4, Real8  !Output data type to choose
  use OwnDataTypes_ml,  only: TXTLEN_DERIV,TXTLEN_SHORT
  use Par_ml,           only: MAXLIMAX, MAXLJMAX, GIMAX,GJMAX,        &
                              me, IRUNBEG, JRUNBEG, limax, ljmax
  use TimeDate_ml,      only: current_date

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
  logical, save     :: debug_flag = .false.
  integer, save     :: i_debug, j_debug       ! Coords matching i,j
  integer msnr                        ! Message number for rsend
  real hourly(MAXLIMAX,MAXLJMAX)      ! Local hourly value  (e.g. ppb)
  real ghourly(GIMAX,GJMAX)           ! Global hourly value (e.g. ppb)
  real :: arrmax                      ! Maximum value from array
  real :: unit_conv                   ! Unit conversion (ppb ug etc.)
  real :: g                           ! tmp - saves value of ghourly(i,j)
  integer, dimension(2) :: maxpos     ! Location of max value
  integer i,j,ih,ispec,itot           ! indices
  integer :: k,ik,iik                 ! Index for vertical level
  integer ist,ien,jst,jen             ! start and end coords
  character(len=TXTLEN_DERIV) :: name ! For output file, species names
  character(len=4)  :: suffix         ! For date "mmyy"
  integer, save :: prev_month = -99   ! Initialise with non-possible month
  logical, parameter :: DEBUG = .false.
  type(Deriv) :: def1           ! for NetCDF
  real :: scale                 ! for NetCDF
  integer ::CDFtype,nk,klevel   ! for NetCDF
  character(len=TXTLEN_SHORT)        :: hr_out_type ! hr_out%type
  integer                            :: hr_out_nk   ! hr_out%nk
  integer, allocatable, dimension(:) :: gspec       ! group array of indexes
  real,    allocatable, dimension(:) :: gunit_conv  ! group array of unit conv. factors

  if ( NHOURLY_OUT <= 0 ) then
    if ( MasterProc .and. DEBUG ) print *,"DEBUG Hourly_out: nothing to output!"
    return
  endif

  if ( my_first_call ) then

  !/ Ensure that domain limits specified in My_Outputs lie within
  !  model domain. In emep coordinates we have:

    do ih = 1, NHOURLY_OUT
      hr_out(ih)%ix1 = max(IRUNBEG,hr_out(ih)%ix1)
      hr_out(ih)%iy1 = max(JRUNBEG,hr_out(ih)%iy1)
      hr_out(ih)%ix2 = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
      hr_out(ih)%iy2 = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)
      hr_out(ih)%nk  = min(KMAX_MID,hr_out(ih)%nk)
    enddo ! ih

    if ( DEBUG ) then
      do j = 1, ljmax
        do i = 1, limax
          if ( i_fdom(i)==DEBUG_i .and. j_fdom(j)==DEBUG_j) then
            debug_flag = .true.
            i_debug = i
            j_debug = j
            !print *, "DEBUG FOUNDIJ me ", me, " IJ ", i, j
          endif
        enddo
      enddo
    endif ! DEBUG
    my_first_call = .false.
  endif  ! first_call

  if( MasterProc .and. Hourly_ASCII .and. current_date%month/=prev_month ) then
    if(prev_month>0) close(IO_HOURLY)      ! Close last-months file
    prev_month = current_date%month

    !/.. Open new file for write-out
    write(suffix,"(2i2.2)")current_date%month,modulo(current_date%year,100)
    name = "Hourly." // suffix
    open(file=name,unit=IO_HOURLY,action="write")

    ! Write summary of outputs to top of Hourly file
    write(IO_HOURLY,"(1(I0,1X,A))")          &
      NHOURLY_OUT,    "Outputs",             &
      FREQ_HOURLY,    "Hours betwen outputs",&
      NLEVELS_HOURLY, "Max Level(s)"
    write(IO_HOURLY,"(1(a21,a21,a10,i4,5i4,a21,es12.5,es10.3))")hr_out
  endif

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

    hr_out_type=hr_out(ih)%type
    hr_out_nk=hr_out(ih)%nk
    if(any(hr_out_type==(/"ADVppbv     ","ADVugXX     ","ADVugXXgroup",&
                          "COLUMN      ","COLUMNgroup "/)))hr_out_nk=1

    KVLOOP: do k = 1,hr_out_nk

      msnr  = 3475 + ih
      ispec = hr_out(ih)%spec
      name  = hr_out(ih)%name
      if ( DEBUG .and. debug_flag ) &
        print "(A,2(1X,I0),1X,A,/A,1X,A)",&
          "DEBUG OH", me, ispec, trim(name),&
          "INTO HOUR TYPE", trim(hr_out(ih)%type)

      if(any(hr_out_type==(/"COLUMN     " ,"COLUMNgroup"/)))then
        ik=KMAX_MID-hr_out(ih)%nk+1  ! top of the column
        if(ik>=KMAX_MID)ik=1         ! 1-level column does not make sense
      else
        ik=KMAX_MID-k+1              ! all levels from model bottom are outputed,
        if(SELECT_LEVELS_HOURLY)then ! or the output levels are taken
          ik=LEVELS_HOURLY(k)        ! from LEVELS_HOURLY array (default)
          hr_out_type=hr_out(ih)%type
          if(ik==0)then
            ik=KMAX_MID              ! surface/lowermost level
            if(any(hr_out_type==(/"BCVppbv     ","BCVugXX     ","BCVugXXgroup"/)))&
              hr_out_type(1:3)="ADV" ! ensure surface output
          else
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

      OPTIONS: select case ( trim(hr_out_type) )
        case ( "ADVppbv" )
          itot = NSPEC_SHL + ispec
          name = species(itot)%name
          unit_conv =  hr_out(ih)%unitconv
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                        * cfac(ispec,i,j) &    ! 50m->3m conversion
                        * unit_conv            ! Units conv.
          end forall

        case ( "BCVppbv" )
          itot = NSPEC_SHL + ispec
          name = species(itot)%name
          unit_conv =  hr_out(ih)%unitconv
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = xn_adv(ispec,i,j,ik) & !BCV:KMAX_MID) &
                   !BCV * cfac(ispec,i,j) &    ! 50m->3m conversion
                        * unit_conv            ! Units conv.
          end forall
          if ( DEBUG .and. debug_flag ) &
            print "(2(A,'=',I0,1X))", "K-level", ik, trim(name), itot

        case ( "ADVugXX" )  !ug/m3, ugX/m3 output at the surface
          itot = NSPEC_SHL + ispec
          name = species(itot)%name
          unit_conv =  hr_out(ih)%unitconv
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                        * cfac(ispec,i,j) &     ! 50m->3m conversion
                        * unit_conv       &     ! Units conv.
                        * roa(i,j,KMAX_MID,1)   ! density.
          end forall

        case ( "BCVugXX" )  ! ug/m3, ugX/m3 output at model mid-levels
          itot = NSPEC_SHL + ispec
          name = species(itot)%name
          unit_conv =  hr_out(ih)%unitconv
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = xn_adv(ispec,i,j,ik) &
                  !BCV  * cfac(ispec,i,j) &     ! 50m->3m conversion
                        * unit_conv       &     ! Units conv.
                        * roa(i,j,ik,1)         ! density.
          end forall
          if ( DEBUG .and. debug_flag ) &
            print "(2(A,'=',I0,1X))", "K-level", ik, trim(name), itot

        case ( "ADVugXXgroup" )  ! GROUP output in ug/m3, ugX/m3 at the surface
          call group_setup()
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = dot_product(xn_adv(gspec,i,j,KMAX_MID),&
                                      cfac(gspec,i,j) & ! 50m->3m conversion
                                     *gunit_conv(:))  & ! Units conv.
                        * roa(i,j,KMAX_MID,1)           ! density.
          end forall
          if ( DEBUG .and. debug_flag ) &
            print "(A,1X,A,'=',30(I0,:,'+'))", "Surface", trim(name), gspec+NSPEC_SHL
          deallocate(gspec,gunit_conv)

        case ( "BCVugXXgroup" )  ! GROUP output in ug/m3, ugX/m3 at model mid-levels
          call group_setup()
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = dot_product(xn_adv(gspec,i,j,ik), &
                                      gunit_conv(:))  & ! Units conv.
                        * roa(i,j,ik,1)                 ! density.
          end forall
          if ( DEBUG .and. debug_flag ) &
            print "(A,'=',I0,1X,A,'=',30(I0,:,'+'))", "K-level", ik, trim(name), gspec+NSPEC_SHL
          deallocate(gspec,gunit_conv)

        case ( "PMwater" )  ! PM water content in ug/m3
          if(trim(hr_out(ih)%unit)/="ug/m3")hr_out(ih)%unit="ug"
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = PM25_water_rh50(i,j)
          end forall

        case ( "COLUMN" )    ! Column output in ug/m2, ugX/m2, molec/cm2
          itot = NSPEC_SHL + ispec
          name = species(itot)%name
          unit_conv =  hr_out(ih)%unitconv
          if(ih>0) hourly(:,:) = 0.0
          do iik=ik,KMAX_MID
            forall ( i=1:limax, j=1:ljmax)
              hourly(i,j) = hourly(i,j)                     &
                        + xn_adv(ispec,i,j,iik)             &
                        * roa(i,j,iik,1)                    & ! density.
                        * (z_bnd(i,j,iik)-z_bnd(i,j,iik+1))   ! level thickness
            end forall
          enddo

        case ( "COLUMNgroup" )! GROUP Column output in ug/m2, ugX/m2, molec/cm2
          call group_setup()
          if(ih>1) hourly(:,:) = 0.0
          do iik=ik,KMAX_MID
            forall ( i=1:limax, j=1:ljmax)
              hourly(i,j) = hourly(i,j)                     &
                        + dot_product(xn_adv(gspec,i,j,iik),&
                                      gunit_conv(:))        & ! Units conv.
                        * roa(i,j,iik,1)                    & ! density.
                        * (z_bnd(i,j,iik)-z_bnd(i,j,iik+1))   ! level thickness
            end forall
          enddo
          if ( DEBUG .and. debug_flag ) &
            print "(A,'=',I0,1X,A,'=',30(I0,:,'+'))", "K-level", ik, trim(name), gspec+NSPEC_SHL
          deallocate(gspec,gunit_conv)

        case ( "SHLmcm3" )        ! No cfac for short-lived species
          itot = ispec
          name = species(itot)%name
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = xn_shl(ispec,i,j,KMAX_MID) &
                        * hr_out(ih)%unitconv  ! Units conv.
          end forall

        case ( "T2_C   " )        ! No cfac for surf.variable
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = t2_nwp(i,j,1) - 273.15     ! Skip Units conv.
          end forall

        case ( "theta  " )        ! No cfac for met.variable
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = th(i,j,KMAX_MID,1)  ! Skip Units conv.
          end forall

        case ( "PRECIP " )        ! No cfac for surf.variables
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = surface_precip(i,j)     ! Skip Units conv.
          end forall

        case ( "Idirect" )        !  Direct radiation (W/m2)
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = Idirect(i,j)    ! Skip Units conv.
          end forall

        case ( "Idiffus" )        !  Diffuse radiation (W/m2)
          forall ( i=1:limax, j=1:ljmax)
            hourly(i,j) = Idiffuse(i,j)    ! Skip Units conv.
          end forall

        case ( "D2D" )
          call CheckStop(ispec<1.or.ispec>num_deriv2d,"ERROR-DEF! Hourly_out: "&
                         //trim(hr_out(ih)%name)//", wrong D2D id!")
          if(hr_out(ih)%unit=="") hr_out(ih)%unit = f_2d(ispec)%unit
          unit_conv = hr_out(ih)%unitconv*f_2d(ispec)%scale
          if(f_2d(ispec)%avg)then           ! non accumulated variables
            forall ( i=1:limax, j=1:ljmax)
              hourly(i,j) = d_2d(ispec,i,j,IOU_INST) * unit_conv
            end forall
          else                              ! hourly accumulated variables
            forall ( i=1:limax, j=1:ljmax)
              hourly(i,j) = (d_2d(ispec,i,j,IOU_YEAR)&
                            -d_2d(ispec,i,j,IOU_HOUR_PREVIOUS)) * unit_conv
              d_2d(ispec,i,j,IOU_HOUR_PREVIOUS)=d_2d(ispec,i,j,IOU_YEAR)
            end forall
          endif
          if( DEBUG  .and. debug_flag) &
            print "(a,2i3,2es12.3)","HHH DEBUG D2D", ispec, ih, &
              hr_out(ih)%unitconv, hourly(i_debug,j_debug)

        case DEFAULT
          call CheckStop( "ERROR-DEF! Hourly_out: "//trim(hr_out(ih)%type)//&
                          " hourly type not found!")

      end select OPTIONS

      if(DEBUG .and. debug_flag ) then
        i = i_debug
        j = j_debug
        print *,"DEBUG-HOURLY-TH ",me,ih,ispec,hourly(i,j),&
                hr_out(ih)%unitconv
      endif

      !/ Get maximum value of hourly array
      hourly(limax+1:MAXLIMAX,:) = 0.0
      hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0
      arrmax = maxval(hourly)

      if ((hr_out(ih)%max>0.0).and.(arrmax>hr_out(ih)%max)) then
        write(6,*) "Hourly value too big!: ", ih, trim(hr_out(ih)%type), arrmax
        write(6,*) "Species : ", trim(name)," : ",  " ispec ", ispec
        write(6,*) "max allowed is : ",  hr_out(ih)%max
        write(6,*) "unitconv was   : ", hr_out(ih)%unitconv
        write(6,*) " me, limax, ljmax, MAXLIMAX,MAXLJMAX : ",  me, &
                          limax, ljmax ,MAXLIMAX,MAXLJMAX
        maxpos = maxloc(hourly)
        write(6,*) "Location is i=", maxpos(1), " j=", maxpos(2)
        write(6,*) "EMEP coords ix=", i_fdom(maxpos(1)), " iy=", j_fdom(maxpos(2))
        write(6,*) "hourly is ", hourly(maxpos(1),maxpos(2))
        if ( hr_out(ih)%type(1:3) == "ADV" ) then
          write(6,*) "xn_ADV is ", xn_adv(ispec,maxpos(1),maxpos(2),KMAX_MID)
          write(6,*) "cfac   is ",   cfac(ispec,maxpos(1),maxpos(2))
        end if
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
      nk = min(KMAX_MID,hr_out_nk)
      CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
      scale=1.

      if (nk == 1) then       !write as 2D
        call Out_netCDF(IOU_HOUR,def1,2,1,hourly,scale,CDFtype,ist,jst,ien,jen)
      elseif( nk > 1 ) then   !write as 3D
        klevel=ik
        if(nk<KMAX_MID)  klevel=KMAX_MID-ik+1 !count from ground and up
        if(SELECT_LEVELS_HOURLY) klevel=k     !order is defined in LEVELS_HOURLY
        call Out_netCDF(IOU_HOUR,def1,3,1,hourly,scale,CDFtype,ist,jst,ien,jen,klevel)
        !else nk<1 : no output
      endif

      if(Hourly_ASCII)then

      !/ Send to ghourly
        call local2global(hourly,ghourly,msnr)

        if ( MasterProc ) then
          !....   write out for a sub-section of the grid:

          !/** We need to correct for small run-domains and the asked-for
          !    output domain. We can only print out the intersection of
          !    these two rectangles.

          write(IO_HOURLY,"('Spec ',i3,' Var ',i2,' = ',2a12,'Lev: ',i2,' Date:',i5,3i3)")  &
            ispec, ih, name, hr_out(ih)%name, ik,                        &
            current_date%year,current_date%month,current_date%day,current_date%hour

          if ( DEBUG .and. debug_flag ) print *, "TTTHOUR ISTS", me, ist, ien, jst, jen

          !/ In model coordinates we have:
          ist = max(1,hr_out(ih)%ix1-IRUNBEG+1)
          jst = max(1,hr_out(ih)%iy1-JRUNBEG+1)
          ien = min(GIMAX,hr_out(ih)%ix2-IRUNBEG+1)
          jen = min(GJMAX,hr_out(ih)%iy2-JRUNBEG+1)

          do i = ist,ien
            do j = jst,jen
              g = ghourly(i,j)
              if ( g /= 0.0 ) then
                write(IO_HOURLY,trim(hr_out(ih)%ofmt) ) g
              else ! Save disc-space used by thousands of  0.00000
                write(IO_HOURLY,"(i1)" ) 0
              endif
            enddo ! j
          enddo   ! i

        endif  ! me loop

      endif !Hourly_ASCII

    enddo KVLOOP
  enddo HLOOP

  contains

  subroutine group_setup()
    if(allocated(gspec)) deallocate(gspec)
    select case ( hr_out(ih)%spec )
      case (1:size(chemgroups))
        name = trim(chemgroups(hr_out(ih)%spec)%name)//"_"//trim(hr_out(ih)%unit)
        allocate(gspec(size(chemgroups(hr_out(ih)%spec)%ptr)))
        gspec=chemgroups(hr_out(ih)%spec)%ptr-NSPEC_SHL
      case DEFAULT
        call CheckStop( "ERROR-DEF! Hourly_out: "//&
                        " hourly variable "//trim(hr_out(ih)%name)//&
                        " type "//trim(hr_out(ih)%type)//", wrong group id!")
    end select
    if(allocated(gunit_conv)) deallocate(gunit_conv)
    allocate(gunit_conv(size(gspec)))
    select case (trim(hr_out(ih)%unit(1:3)))
      case("ppb" )     ; gunit_conv(:)=unit_conv  ! this option makes no sence
      case("ug","ug/" ); gunit_conv(:)=to_ug_ADV(gspec)
      case("ugC")      ; gunit_conv(:)=to_ug_C(gspec)
      case("ugS")      ; gunit_conv(:)=to_ug_S(gspec)
      case("ugN")      ; gunit_conv(:)=to_ug_N(gspec)
      case default     ; gunit_conv(:)=unit_conv
        call CheckStop(index(hr_out(ih)%unit,"molec")==0,&
                        "ERROR-DEF! Hourly_out: "//trim(hr_out(ih)%type)//&
                        "hourly, wrong group unit='"//trim(hr_out(ih)%unit)//"'!")
    end select
    if ( DEBUG .and. debug_flag ) print "(A,'=',30(A,':',I0,:,'+'))",&
      trim(name),(trim(species(gspec(ispec)+NSPEC_SHL)%name),gspec(ispec),ispec=1,size(gspec))
  end subroutine group_setup

end subroutine hourly_out
