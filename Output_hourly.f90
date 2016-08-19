! <Output_hourly.f90 - A component of the EMEP MSC-W Unified Eulerian
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
   use My_Outputs_ml,    only : NHOURLY_OUT, &      ! No. outputs
                                 NLEVELS_HOURLY, &  ! ds rv1_8_2 
                                 FREQ_HOURLY, &     ! No. hours between outputs
                                 Asc2D, hr_out, &   ! Required outputs
                                 Hourly_ASCII       ! ASCII output or not

   use CheckStop_ml,     only : CheckStop
   use Chemfields_ml ,   only : xn_adv,xn_shl, cfac
   use Derived_ml,    only : d_2d, IOU_INST,IOU_HOUR,Deriv
   use GenSpec_shl_ml ,  only : NSPEC_SHL        ! Maps indices
   use GenChemicals_ml , only : species          ! Gives names
   use GridValues_ml,    only : i_fdom, j_fdom   ! Gives emep coordinates
   use Io_ml,            only : IO_HOURLY
   use ModelConstants_ml,only : NPROC,KMAX_MID,DEBUG_i,DEBUG_j,identi,runlabel1
   use Met_ml,           only : t2_nwp,th, roa, surface_precip, &
                                   Idirect, Idiffuse
   use NetCDF_ml,        only : Out_netCDF,Init_new_netCDF &
                                ,Int1,Int2,Int4,Real4,Real8 !Output data type to choose
   use Par_ml ,          only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX    &
                                ,li0,li1,lj0,lj1 & 
                                ,me,IRUNBEG,JRUNBEG,limax,ljmax
   use TimeDate_ml      ,only : current_date

   implicit none

   !*.. Components of  hr_out
   !*  character(len=3) :: type   ! "ADVp" or "ADVu" or "SHL" or "T2 "
   !*  integer          :: spec   ! Species number in xn_adv or xn_shl array
   !* character(len=12) :: ofmt   ! Output format (e.g. es12.4)
   !*  integer          :: ix1    ! bottom-left x
   !*  integer          :: iy1    ! bottom-left y
   !*  integer          :: ix2    ! upper-right x
   !*  integer          :: iy2    ! upper-right y
   !*  real             :: unitconv   !  conv. factor
   !*  real             :: max    ! max allowed value

   ! local variables
   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE),INFO
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
   integer :: ik                       ! Index for vertical level
   integer ist,ien,jst,jen             ! start and end coords
   character(len=50) :: errmsg = "ok"  ! For  consistecny check
   character(len=20) :: name           ! For output file, species names
   character(len=120) :: netCDFName    ! For netCDF output filename
   character(len=4)  :: suffix         ! For date "mmyy"
   integer, save :: prev_month = -99   ! Initialise with non-possible month
   logical, parameter :: DEBUG = .false.
   integer :: NLEVELS_HOURLYih
   type(Deriv) :: def1 !for NetCDF
   real :: scale !for NetCDF
   integer ::CDFtype,nk,klevel!for NetCDF

    if ( my_first_call ) then

      !/ Ensure that domain limits specified in My_Outputs lie within
      !  model domain. In emep coordinates we have:

        do ih = 1, NHOURLY_OUT

           hr_out(ih)%ix1 = max(IRUNBEG,hr_out(ih)%ix1)
           hr_out(ih)%iy1 = max(JRUNBEG,hr_out(ih)%iy1)
           hr_out(ih)%ix2 = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
           hr_out(ih)%iy2 = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)
           hr_out(ih)%nk = min(KMAX_MID,hr_out(ih)%nk)

        end do ! ih
        if ( DEBUG ) then
           do j = 1, ljmax
              do i = 1, limax
                  if ( i_fdom(i)==DEBUG_i .and. j_fdom(j)==DEBUG_j) then
                       debug_flag = .true.
                       i_debug = i
                       j_debug = j
                       !print *, "DEBUG FOUNDIJ me ", me, " IJ ", i, j
                  end if
              end do
            end do
        end if ! DEBUG
        my_first_call = .false.
    end if  ! first_call

   !     hourly(:,:) = 0.0      ! Initialise (ljmax+1:MAXLJMAX, limax+1:LIMAX
   !                            !  would have done,  but this is simpler)
   ! else                        
   ! Mask the edges of the hourly array, so that we can use maxval later
   ! This makes the code a bit neater below, but costs some CPU time here,
   ! and in evaluating maxval over the whole MAXLIMAX*MAXLJMAX dimension.

        !u7.5vg FIX hourly(limax+1:MAXLIMAX,:) = 0.0
        !u7.5vg FIX hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0

     hourly(:,:) = 0.0

    !end if  ! first_call

   if(me == 0 .and. current_date%month /= prev_month ) then

        if ( prev_month > 0 .and. Hourly_ASCII) close(IO_HOURLY)      ! Close last-months file

       !/.. Open new file for write-out

        write(suffix,fmt="(2i2.2)") current_date%month, &
                           modulo ( current_date%year, 100 )
        if(Hourly_ASCII)then
           name = "Hourly" // "." // suffix
           open(file=name,unit=IO_HOURLY,action="write")
        endif

        prev_month = current_date%month

!        netCDFName =trim(runlabel1)//"_hour" // "."// suffix // ".nc"
!        call Init_new_netCDF(netCDFName,IOU_HOUR)

        if(Hourly_ASCII)then
       !ds rv1.6.2: Write summary of outputs to top of Hourly file
       !  - remember - with corrected domain limits here
        write(IO_HOURLY,*) NHOURLY_OUT, " Outputs"
        write(IO_HOURLY,*) FREQ_HOURLY, " Hours betwen outputs"
        write(IO_HOURLY,*) NLEVELS_HOURLY, "Max Level(s)"    !ds rv1_8_2

        do ih = 1, NHOURLY_OUT
           write(IO_HOURLY,fmt="(a12,a8,a10,i4,5i4,a13,es12.5,es10.3)") hr_out(ih)
        end do
        endif !Hourly_ASCII
   end if


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


   HLOOP: do ih = 1, NHOURLY_OUT
      NLEVELS_HOURLYih=hr_out(ih)%nk
   KVLOOP: do ik = KMAX_MID, KMAX_MID-NLEVELS_HOURLYih+1, -1

      msnr  = 3475 + ih
      ispec = hr_out(ih)%spec 
      name  = hr_out(ih)%name   !ds rv1.6.1 
      if ( DEBUG .and. debug_flag ) print *, "DEBUG OH ", me, ispec, name,  &
         hr_out(ih)%type

       if(DEBUG .and. debug_flag ) print *, "INTO HOUR TYPE ",hr_out(ih)%type

   !----------------------------------------------------------------
   ! Multi-layer output. 
   !  Specify NLEVELS_HOURLY here, and in hr_out defs use either:
   !
   !      ADVppbv to get surface concentrations (onyl relevant for
   !              layer k=20 of course - gives meaningless  number f
   !               or higher levels.
   ! Or,
   !      BCVppbv to get grid-centre concentrations (relevant for
   !      all layers.
   !----------------------------------------------------------------


       OPTIONS: select case ( hr_out(ih)%type ) 
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
            if ( DEBUG .and. debug_flag ) print *, "K-level", ik, name, itot

         case ( "ADVugm3" )
            itot = NSPEC_SHL + ispec 
            name = species(itot)%name
            unit_conv =  hr_out(ih)%unitconv * species(itot)%molwt
            forall ( i=1:limax, j=1:ljmax)
                  hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                                 * cfac(ispec,i,j) &     ! 50m->3m conversion
                                 * unit_conv       &     ! Units conv.
                                 * roa(i,j,KMAX_MID,1)   ! density.
            end forall

          case ( "SHLmcm3" )        ! No cfac for short-lived species
            itot = ispec 
            name = species(itot)%name
            forall ( i=1:limax, j=1:ljmax)
                     hourly(i,j) = xn_shl(ispec,i,j,KMAX_MID) &
                                    * hr_out(ih)%unitconv  ! Units conv.
            end forall

          case ( "T2_C   " )        ! No cfac for short-lived species
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = t2_nwp(i,j,1) - 273.15     ! Skip Units conv.
            end forall

          case ( "theta  " )        ! No cfac for short-lived species
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = th(i,j,KMAX_MID,1)  ! Skip Units conv.
            end forall

          case ( "PRECIP " )        ! No cfac for short-lived species
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

          case ( "D2D" )        ! No cfac for short-lived species

            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = d_2d(ispec,i,j,IOU_INST) * hr_out(ih)%unitconv
            end forall

    	    if( DEBUG  .and. debug_flag) &
               write(6,"(a12,2i3,2es12.3)") "HHH DEBUG", ispec, ih, &
                 hr_out(ih)%unitconv, hourly(i_debug,j_debug)

          case DEFAULT 
             errmsg = "ERROR-DEF! Hourly_out: " // hr_out(ih)%type
             call CheckStop( errmsg  // "hourly type not found!")

       end select OPTIONS 

	if(DEBUG .and. debug_flag ) then
             i = i_debug
             j = j_debug
             print *,"DEBUG-HOURLY-TH ",me,ih,ispec,hourly(i,j),&
                      hr_out(ih)%unitconv
        end if


     !ds rv1.6.2 ---- why needed?
        hourly(limax+1:MAXLIMAX,:) = 0.0
        hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0

      !/ Get maximum value of hourly array

       arrmax = maxval(hourly)
       if ( arrmax  >   hr_out(ih)%max ) then
            write(6,*) "Hourly value too big!: ", ih, hr_out(ih)%type, arrmax
            write(6,*) "Species : ", name," : ",  " ispec ", ispec
            write(6,*) "max allowed is : ",  hr_out(ih)%max
            write(6,*) "unitconv was   : ", hr_out(ih)%unitconv
            write(6,*) " me, limax, ljmax, MAXLIMAX,MAXLJMAX : ",  me, &
                             limax, ljmax ,MAXLIMAX,MAXLJMAX
            maxpos = maxloc(hourly)
            write(6,*) "Location is i=", maxpos(1), " j=", maxpos(2)
            write(6,*) "EMEP coords ix=", i_fdom(maxpos(1)), " iy=", j_fdom(maxpos(2))
            write(6,*) "hourly is ", hourly(maxpos(1),maxpos(2))
            if ( hr_out(ih)%type == "ADV" ) then
              write(6,*) "xn_ADV is ", xn_adv(ispec,maxpos(1),maxpos(2),KMAX_MID)
              write(6,*) "cfac   is ",   cfac(ispec,maxpos(1),maxpos(2))
            end if

              call CheckStop("Error, Output_hourly/hourly_out: too big!")

       endif

!NetCDF hourly output
       def1%name=hr_out(ih)%name
       def1%unit=hr_out(ih)%unit
       def1%class=hr_out(ih)%type 
       ist = max(1,hr_out(ih)%ix1-IRUNBEG+1)
       jst = max(1,hr_out(ih)%iy1-JRUNBEG+1)
       ien = min(GIMAX,hr_out(ih)%ix2-IRUNBEG+1)
       jen = min(GJMAX,hr_out(ih)%iy2-JRUNBEG+1)
       nk = min(KMAX_MID,hr_out(ih)%nk)
       CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
       scale=1.

       if (nk == 1) then !write as 2D
          call Out_netCDF(IOU_HOUR,def1,2 &
            ,1,hourly(:,:),scale,CDFtype,ist,jst,ien,jen)

       else if( nk > 1 ) then   !write as 3D
          !CHANGED 23 Mar 2007  klevel=KMAX_MID-ik+1
          klevel=ik
          call Out_netCDF(IOU_HOUR,def1,3 &
            ,1,hourly(:,:),scale,CDFtype,ist,jst,ien,jen,klevel)
          !else nk<1 : no output
       endif

        if(Hourly_ASCII)then

      !/ Send to ghourly

       call local2global(hourly,ghourly,msnr)

       if (me ==  0) then

            !....   write out for a sub-section of the grid:

            !/** We need to correct for small run-domains and the asked-for
            !    output domain. We can only print out the intersection of
            !    these two rectangles.

            !ds!/ In emep coordinates we have:

            !dsist = max(IRUNBEG,hr_out(ih)%ix1)
            !dsjst = max(JRUNBEG,hr_out(ih)%iy1)
            !dsien = min(GIMAX+IRUNBEG-1,hr_out(ih)%ix2)
            !dsjen = min(GJMAX+JRUNBEG-1,hr_out(ih)%iy2)

            !ds rv1_8_2 Extra info:
            write(IO_HOURLY,"('Spec ',i3,' Var ',i2,' = ',2a12,'Lev: ',i2,' Date:',i5,3i3)")  &
                  ispec,  ih, name, hr_out(ih)%name                          &
                 ,ik  &  ! ds rv1_8_2 
                 ,current_date%year,current_date%month,current_date%day       &
                 ,current_date%hour !ds                                           &
                 !ds ,ist, ien, jst, jen,         &
                 !ds  unit_conv

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
                    write(IO_HOURLY, fmt=hr_out(ih)%ofmt ) g
                else ! Save disc-space used by thousands of  0.00000
                    write(IO_HOURLY, fmt="(i1)" ) 0
                end if

              end do ! j
            end do   ! i

       end if  ! me loop

       endif !Hourly_ASCII

      end do KVLOOP
      end do HLOOP

  end subroutine hourly_out
