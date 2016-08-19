! <PhyChem_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module PhyChem_ml
!
!     physical and chemical routine calls within one advection step
!     driven from here
!   
!     Output of hourly data
!
!-----------------------------------------------------------------------------
   use My_Outputs_ml , only : NHOURLY_OUT, FREQ_SITE, FREQ_SONDE, FREQ_HOURLY
   use My_Timing_ml,   only : Code_timer, Add_2timing, tim_before, tim_after  

   use Advection_ml,   only: advecdiff,advecdiff_poles,adv_int
   use Chemfields_ml,  only : xn_adv,cfac,xn_shl
   use Derived_ml,     only : IOU_INST, DerivedProds, Derived, &
                               num_deriv2d,d_2d, f_2d
   use DryDep_ml,      only : drydep,init_drydep
   use Emissions_ml,   only : EmisSet  
   use GridValues_ml,  only : debug_proc, debug_li,debug_lj,& !ds jun2005
                             gl, gb, projection, Poles
   use Met_ml ,        only : roa,z_bnd,z_mid,metint, ps, cc3dmax, &
                               zen,coszen,Idirect,Idiffuse
   use ModelConstants_ml, only : KMAX_MID, nmax, nstep &
                        ,dt_advec  &    ! time-step for phyche/advection
                        ,END_OF_EMEPDAY ! (usually 6am)
   use Nest_ml,        only : readxn, wrtxn
   use Par_ml,         only : me, MAXLIMAX, MAXLJMAX
   use TimeDate_ml,       only : date,daynumber,day_of_year, add_secs, &
                                 current_date, timestamp,  &
                                 make_timestamp, make_current_date
   use Trajectory_ml,  only : trajectory_out     ! 'Aircraft'-type  outputs
   use Radiation_ml,   only : SolarSetup,       &! sets up radn params
                             ZenithAngle,      &! gets zenith angle
                             ClearSkyRadn,     &! Idirect, Idiffuse
                             CloudAtten         ! 
   use Runchem_ml,     only : runchem   ! Calls setup subs and runs chemistry
   use Sites_ml,       only: siteswrt_surf, siteswrt_sondes    ! outputs
   use Timefactors_ml, only : NewDayFactors  
!-----------------------------------------------------------------------------
implicit none
private

public :: phyche


contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 subroutine phyche(numt)
   integer, intent(in) ::  numt

   integer ::  i,j,k,n
   logical, parameter :: DEBUG = .false.
   logical, save :: End_of_Day = .false.

   integer :: ndays
   real :: thour


   type(timestamp) :: ts_now !date in timestamp format

    !------------------------------------------------------------------
    !     start of inner time loop which calls the physical and
    !     chemical routines.


    DO_OUTER: do nstep = 1,nmax

       !     Hours since midnight at any time-step
       !	using current_date we have already nstep taken into account

       thour = real(current_date%hour) + current_date%seconds/3600.0 & 
                   + 0.5*dt_advec/3600.0 

       if ( DEBUG .and. debug_proc ) then
           write(6,*) "PhyChe debug ", me, thour,  &
                      current_date%hour, current_date%seconds

           if ( current_date%hour == 12 ) then

               ndays = day_of_year(current_date%year,current_date%month, &
                                    current_date%day)
               write(6,*) 'thour,ndays,nstep,dt', thour,ndays,nstep,dt_advec
           endif

        endif

        if (me == 0) write(6,"(a15,i6,f8.3)") 'timestep nr.',nstep,thour

        call wrtxn(current_date,.false.) !Write xn_adv for future nesting
        call readxn(current_date) !Read xn_adv from earlier runs
 
!        ==================
         call Code_timer(tim_before)

        call EmisSet(current_date)
        call Add_2timing(15,tim_after,tim_before,"phyche:EmisSet")

!       For safety we initialise add instant. values here to zero.
!       Usually not needed, but sometimes
!       ==================
         d_2d(:,:,:,IOU_INST) = 0.0
!       ==================


       !===================================

        call SolarSetup(current_date%year,current_date%month, &
                           current_date%day,thour)

        call ZenithAngle(thour, gb, gl, zen, coszen )

        if( DEBUG .and. debug_proc  ) then
          write(*,*) "PhyChem ZenRad ", current_date%day, current_date%hour, &
               thour, gl(debug_li,debug_lj),gb(debug_li,debug_lj), &
                     zen(debug_li,debug_lj),coszen(debug_li,debug_lj)
        end if

        call ClearSkyRadn(ps(:,:,1),coszen,Idirect,Idiffuse)

        call CloudAtten(cc3dmax(:,:,KMAX_MID),Idirect,Idiffuse)

       !===================================
        call Add_2timing(16,tim_after,tim_before,"phyche:ZenAng")


        !================
        if( (Poles(1)==1.or.Poles(2)==1).and. &
                              trim(projection)==trim('lon lat'))then
           call advecdiff_poles
        else
           call advecdiff
        endif

        call Add_2timing(17,tim_after,tim_before,"phyche:advecdiff")
        !================

        call Code_timer(tim_before)


       !/ See if we are calculating any before-after chemistry productions:

          !=============================
          if ( nstep == nmax ) call DerivedProds("Before",dt_advec)
          !=============================

          call Add_2timing(26,tim_after,tim_before,"phyche:MACHO-prod")

         !===================================
           call init_drydep()
         !===================================

         !=========================================================

          call runchem(numt)   !  calls setup subs and runs chemistry

          call Add_2timing(28,tim_after,tim_before,"Runchem")

          !=========================================================
                           

         !/ See if we are calculating any before-after chemistry productions:

          !=============================
          if ( nstep == nmax ) call DerivedProds("After",dt_advec)
          !=============================

          call Code_timer(tim_before)
          !=============================
          call Add_2timing(34,tim_after,tim_before,"phyche:drydep")



          !=============================
          ! this output needs the 'old' current_date_hour

           call trajectory_out
          !=============================

!	the following partly relates to end of time step - hourly output
!	partly not depends on current_date
!	=> add dt_advec to current_date already here


          !====================================
          ts_now = make_timestamp(current_date)

          call add_secs(ts_now,dt_advec)

          current_date = make_current_date(ts_now)

          !====================================



          End_of_Day = (current_date%seconds == 0 .and. &
                        current_date%hour    == END_OF_EMEPDAY)

          if( End_of_Day .and. me == 0 ) then
              write(*,"(a20,2i4,i6)") "END_OF_EMEPDAY, Hour,seconds=", &
                END_OF_EMEPDAY, current_date%hour,current_date%seconds
          endif

          call Derived(dt_advec,End_of_Day)


         ! Hourly Outputs:
          if ( current_date%seconds == 0 ) then

              if ( modulo(current_date%hour, FREQ_SITE) == 0 )  &
                                 call siteswrt_surf(xn_adv,cfac,xn_shl)

              if ( modulo(current_date%hour, FREQ_SONDE) == 0 ) &
                  call siteswrt_sondes(xn_adv,xn_shl)

              if ( NHOURLY_OUT > 0 .and.  &
                     modulo(current_date%hour, FREQ_HOURLY) == 0 ) &
                         call hourly_out()

          end if

          call Add_2timing(35,tim_after,tim_before,"phyche:outs")


          call metint


          call adv_int


          call Add_2timing(36,tim_after,tim_before,"phyche:ints")

      enddo DO_OUTER

   end subroutine phyche
!-----------------------------------------------------------------------------
end module PhyChem_ml
