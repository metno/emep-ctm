! <OutputChem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                     module OutputChem_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  use CheckStop_ml,      only: CheckStop 
  use Derived_ml,        only: IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY,  &
                               f_2d, d_2d, LENOUT2D                   &
                              ,f_3d, d_3d, nav_3d, nav_2d, LENOUT3D   & 
                              , num_deriv2d, num_deriv3d &
                              ,ResetDerived, Deriv
  use My_Outputs_ml,     only: NBDATES, wanted_dates_inst,            & 
                               Ascii3D_WANTED
  use Io_ml,             only: IO_WRTCHEM
  use ModelConstants_ml, only: nprint, END_OF_EMEPDAY, KMAX_MID
  use NetCDF_ml,         only: CloseNetCDF, Out_netCDF
  use Par_ml,            only: MAXLIMAX,MAXLJMAX,GIMAX,GJMAX,me,      &
                               IRUNBEG,JRUNBEG
  use TimeDate_ml      , only: current_date, max_day  ! days in month

  implicit none

 !/* subroutines:

  public :: Wrtchem
  public :: Output_fields   ! (iotyp)
  public :: Output_f2d      ! (iotyp, dim, nav, def, dat)
  public :: Output_f3d      ! (iotyp, dim, nav, def, dat)

  logical, private, parameter :: MY_DEBUG = .false.

  contains 

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine Wrtchem(numt)

     !---------------------------------------------------------------------
     ! DESCRIPTION:
     !   Writes out data fields as NetCDF 
     !
     !   END_OF_EMEPDAY = 6am, i.g. EMEP sampling period from 6am to 6am
     !   Daily outputs for "EMEP" days (which end between 0 and 6am) are
     !   dated by the date of sampling start, i.e. date of the previous day.
     !   
     !   Thus, the first output should occur just as Jan 2nd starts (e.g.
     !   at 6am on 2nd Jan); logical Jan_1st helps dealing with this and 
     !   it also marks end of a year run.
     !   (For runs starting in other months, one partial write-out will
     !   occur at 6am of the 1st day, but this should be over-written
     !   as soon as a full day of data is available).
     !----------------------------------------------------------------------


     integer, intent(in) ::  numt

     real, dimension(MAXLIMAX, MAXLJMAX) :: local_2d  !local 2D array
     real, dimension(GIMAX, GJMAX)       :: glob_2d   !array for whole domain
     integer :: i,j,n,k,icmp,msnr1
     integer :: nyear,nmonth,nday,nhour,nmonpr
     integer :: mm_out, dd_out  
     logical :: Jan_1st, End_of_Run
     real    :: scale
     character*30 outfilename
   !------------------------------

     nyear  = current_date%year
     nmonth = current_date%month
     nday   = current_date%day
     nhour  = current_date%hour

     dd_out = nday
     mm_out = nmonth
     Jan_1st    = ( nmonth == 1 .and. nday == 1 )
     End_of_Run = ( mod(numt,nprint) == 0       )

     if(me==0 .and. MY_DEBUG) write(6,"(a12,i5,5i4)") "DAILY DD_OUT ",          &
          numt, nmonth, mm_out, nday, dd_out, nhour


     if ( END_OF_EMEPDAY  <= 7 ) then 

        !. END_OF_EMEPDAY = 6am - end of EMEP daily sampling period
        !. Daily outputs are dated with the start of sampling period 

        dd_out = nday - 1     ! only used for daily outputs

        if(me==0 .and. MY_DEBUG) write(6,"(a12,i5,5i4)") "DAILY SET ",         &
                  numt, nmonth, mm_out, nday, dd_out, nhour 

        if(dd_out == 0) then
           mm_out = nmonth - 1

           if(nmonth == 1)  mm_out = 12

             dd_out = max_day(mm_out, nyear)  !  Last day of month

              if(me==0 .and. MY_DEBUG) write(6,"(a12,i5,4i4)") "DAILY FIX ",     &
                             numt, nmonth, mm_out, nday, dd_out 
        end if
     end if      ! for END_OF_EMEPDAY <= 7


     !== Instantaneous results output ====
     !   Possible actual array output for specified days and hours 
     !   is defined in wanted_dates_bi array in My_Outputs

     do n = 1, NBDATES
        if ( wanted_dates_inst(n)%month == nmonth .and. &
             wanted_dates_inst(n)%day   == nday   .and. &
             wanted_dates_inst(n)%hour  == nhour ) then

           call Output_fields(IOU_INST)

        end if 
     end do


     !== Daily output ====

     if (nhour ==  END_OF_EMEPDAY ) then

        if ( numt > 1 .and. .not. Jan_1st ) then   ! Doesn't write out 1 Jan.

           call Output_fields(IOU_DAY)

        end if

        call ResetDerived(IOU_DAY)    ! For daily averaging, reset also 1 Jan.

     end if

     !== Output at the end of the run

     if ( End_of_Run ) then

        call Output_fields(IOU_DAY)   ! Daily outputs
        call Output_fields(IOU_YEAR)  ! Yearly outputs

     end if


     !/ NEW MONTH

     if (nday == 1 .and. nhour == 0) then
        nmonpr = nmonth-1

        if (nmonpr.eq.0) nmonpr=12

        !== Monthly output ====

        call Output_fields(IOU_MON)  

        !== ASCII output of 3D fields (if wanted)

        if(Ascii3D_WANTED) then

           if (num_deriv3d > 0) then
              msnr1 = 2000

              do n = 1, num_deriv3d

                 if( me == 0 ) then

                    write(outfilename,fmt='(a,a5,i2.2)')   &
                         trim( f_3d(n)%name ), ".out.",  nmonpr
                    open (IO_WRTCHEM,file=outfilename)
                    write(IO_WRTCHEM,fmt="(4i4)") IRUNBEG, GIMAX+IRUNBEG-1, &
                          JRUNBEG, GJMAX+JRUNBEG-1 ! domain
                 end if

                 if (nav_3d(n,IOU_MON) == 0 ) then
                    write(IO_WRTCHEM,*) "ERROR in 3D ASCII output: nav=0"
                 else 

                    do k = 1, KMAX_MID

                      local_2d(:,:) = d_3d(n,:,:,k,IOU_MON)/nav_3d(n,IOU_MON) 
                      call local2global(local_2d,glob_2d,msnr1)

                       if (me  ==  0) then
                          do j=1,GJMAX
                             do i=1,GIMAX
                                write(IO_WRTCHEM,"(es10.3)") glob_2d(i,j)
                             end do
                          end do
                       end if  ! me loop

                    end do ! k
                 end if ! nav == 0

                 if( me == 0 ) close(IO_WRTCHEM)

              end do  ! 3D-variables loop num_deriv3d
           end if

        endif       ! Ascii3D_WANTED

        call ResetDerived(IOU_MON)

     endif              ! End of NEW MONTH

  end subroutine Wrtchem

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine Output_fields(iotyp)

    integer, intent(in) :: iotyp

  !     
  !*** 2D fields, e.g. surface SO2, SO4, NO2, NO3 etc.; AOT, fluxes 
  !--------------------

     if(num_deriv2d > 0) call Output_f2d(iotyp,num_deriv2d,nav_2d,f_2d,d_2d)
  


  !*** 3D concentration fields, e.g. O3
  !--------------------

     if(num_deriv3d > 0) call Output_f3d(iotyp,num_deriv3d,nav_3d,f_3d,d_3d)

     call CloseNetCDF

  end subroutine Output_fields

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine Output_f2d (iotyp, dim, nav, def, dat)

  !=========================================
  ! Sends fields to NetCDF output routines
  !=========================================

    integer,                         intent(in) :: iotyp
    integer,                         intent(in) :: dim ! No. fields
    integer, dimension(dim,LENOUT2D),intent(in) :: nav ! No. items averaged
    type(Deriv), dimension(dim),     intent(in) :: def ! Definition of fields
    real, dimension(dim,MAXLIMAX,MAXLJMAX,LENOUT2D), intent(in) :: dat

    logical :: wanted     ! Set true for required year, month, day or inst.
    integer :: icmp       ! component index
    real    :: scale      ! Scaling factor
  !--------------------------------------

     do icmp = 1, dim
        
        wanted = .false.
        if( iotyp == IOU_YEAR) wanted = def(icmp)%year
        if( iotyp == IOU_MON ) wanted = def(icmp)%month
        if( iotyp == IOU_DAY ) wanted = def(icmp)%day
        if( iotyp == IOU_INST) wanted = def(icmp)%inst
        
        if ( wanted ) then 
           
          scale  = def(icmp)%scale
           if (iotyp /= IOU_INST )    &
             scale = scale / max(1,nav(icmp,iotyp))

             call Out_netCDF(iotyp,def(icmp),2,1,dat(icmp,:,:,iotyp),scale)
           
        endif     ! wanted
     enddo        ! component loop
     
   end subroutine Output_f2d

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine  Output_f3d (iotyp, dim, nav, def, dat)

  !=========================================
  ! Sends fields to NetCDF output routines
  !=========================================
     
   use Derived_ml,        only: IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY, &
                                Deriv,LENOUT3D
   use ModelConstants_ml, only: KMAX_MID
   use NetCDF_ml,         only: Out_netCDF
   use Par_ml,            only: MAXLIMAX, MAXLJMAX 


    implicit none

    integer,                         intent(in) :: iotyp
    integer,                         intent(in) :: dim ! No. fields
    integer, dimension(dim,LENOUT3D),intent(in) :: nav ! No. items averaged
    type(Deriv), dimension(dim),     intent(in) :: def ! definition of fields
    real, dimension(dim,MAXLIMAX,MAXLJMAX,KMAX_MID,LENOUT3D), intent(in):: dat

    logical :: wanted     ! Set true for required year, month, day or inst.
    integer :: icmp       ! component index
    real    :: scale      ! Scaling factor
 !------------------------------------------------

     do icmp = 1, dim

        wanted = .false.
        if( iotyp == IOU_YEAR) wanted = def(icmp)%year
        if( iotyp == IOU_MON ) wanted = def(icmp)%month
        if( iotyp == IOU_DAY ) wanted = def(icmp)%day
        if( iotyp == IOU_INST) wanted = def(icmp)%inst
        
        if ( wanted ) then 
           
          scale = def(icmp)%scale
           if (iotyp /= IOU_INST)   &
               scale = scale /max(1,nav(icmp,iotyp))

                call Out_netCDF(iotyp, def(icmp), 3, KMAX_MID,   &
                                dat(icmp,:,:,:,iotyp), scale)
           
           endif     ! wanted
        enddo        ! component loop
     
   end subroutine Output_f3d

 end module OutputChem_ml
