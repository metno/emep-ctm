! <Volcanos_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Volcanos_ml
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!
!                    module Volcanos_ml
!
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !-----------------------------------------------------------------------!
  ! Processes SOx emission heights from volcanoes 
  ! In volcanos_ml, only the height of the volcanos are read from volcanos.dat,
  ! the emissions themselves comes from gridSOx, so this is to ensure that if
  ! suddenly e.g. Iceland starts to report volcano emissions and these are in
  ! gridSOx, the program will discover this.
  ! Note that nvolc is set according to the emission input from the gridSOx 
  ! files.  It counts the number of grids with the "country code" for volcanoes
  ! that are in the gridSOx file.
  !
  ! Note - EmisGet and other routines use "restricted" coords, which introduces
  ! some complications here. Hopefully we can tidy up one day.
  !-----------------------------------------------------------------------!

 use CheckStop_ml,          only : CheckStop
 use My_Emis_ml,            only : QRCVOL,molwt
 use EmisDef_ml,            only : NSECTORS,ISNAP_NAT
 use GridValues_ml,         only : sigma_bnd, i_fdom, j_fdom,i_local, j_local
 use Io_ml,                 only : ios, open_file, check_file, IO_VOLC
 use ModelConstants_ml,     only : KMAX_BND,KMAX_MID,PT, NPROC
 use Met_ml,                only : ps, roa
 use Par_ml,                only : IRUNBEG, JRUNBEG, me, li0,lj0,li1,lj1  &
                                  ,gi0, gi1, gj0, gj1 !TEST
 use PhysicalConstants_ml,  only : GRAV, AVOG

 implicit none
 private


 !/* subroutines:

  public :: VolcGet
  public :: Set_Volc     
  public :: Scale_Volc   


  integer, public, parameter  :: NMAX_VOLC = 3  ! Max number of volcanoes
  integer, public, save       :: nvolc = 0    & ! No. grids with volcano 
                                                ! emissions in gridSOx
                                ,volc_no = 0 
  integer, public, save, dimension(NMAX_VOLC):: &
                        height_volc,          &  ! Height of volcanoes
                        i_volc, j_volc           ! Volcano's EMEP coordinates
  real, private, save, dimension(NMAX_VOLC)::   &
                        rcemis_volc0 ! Emissions part varying every hour
  real, public, save, dimension(NMAX_VOLC) ::   &
                        rcemis_volc, & ! Emissions part varying every time-step 
                        emis_volc = 0.0 ! Volcanoes' emissions

  logical, private, parameter :: DEBUG_VULC = .false.

contains 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine VolcGet(height_volc)
!***********************************************************************
 !-----------------------------------------------------------------------!
 ! Reads volcanoes' coorinates (i,j) and height level(height) 
 ! Returns to EmisSet with height of volcanos (height_volc)
 ! Input file: Volcanoes.dat
 !-----------------------------------------------------------------------!

  integer, intent(out), dimension(NMAX_VOLC) :: height_volc
  integer            :: nvolc_read,height,i,j          ! Local variables
  character (len=13) :: fname
  logical            :: fexist 
     
     fname = "Volcanoes.dat"  

     ios=0    ! Start with  assumed ok status

     call open_file(IO_VOLC,"r",fname,needed=.true.,skip=1)

     call CheckStop(ios,"VolcGet: problems with Volcanoes.dat ")


     height_volc(:)=0.0
     nvolc_read=0

     READVOLC: do
           read(IO_VOLC,*,iostat=ios) i,j,height

           if (DEBUG_VULC) write(*,*)'found i,j,heigh',i,j,height
           if ( ios /= 0 ) exit READVOLC

  !/** Read (i,j) are given for the full EMEP polar-stereographic domain
  !    Convert them to actual run domain
           i = i -IRUNBEG+1    
           j = j -JRUNBEG+1    

  !/** Set the volcano number to be the same as in emission data (gridSOx)

           do volc_no=1,nvolc
               if ((i_volc(volc_no)==i) .and. (j_volc(volc_no)==j)) then
                  height_volc(volc_no)=height
                  nvolc_read=nvolc_read+1
                  if (DEBUG_VULC) write(*,*)'Found volcano with height k=',height
               endif
           enddo
     enddo READVOLC

     write(6,*) nvolc_read,' volcanos on volcanos.dat &
             & match volcanos on emislist.sox'
     write(6,*) nvolc,' volcanos found in emislist.sox'

     call CheckStop(nvolc_read < nvolc, "Volc missing in Volcanos.dat")
  
     close(IO_VOLC)

    end subroutine VolcGet
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine Set_Volc

 !-----------------------------------------------------------------------!
 ! Starts converting emission units from kg/m2/s to.... (hourly)
 !-----------------------------------------------------------------------!

    !**Local variables
    integer            :: k,i,j, i1,i2,j1,j2
    real               :: unit_conv1  

    rcemis_volc0(:) = 0.0
    unit_conv1      = 0.0

    !/** Set volcano
    do volc_no=1,nvolc
       k=height_volc(volc_no)
       i=i_volc(volc_no) +IRUNBEG-1   !NEW
       j=j_volc(volc_no) +JRUNBEG-1   !NEW

       if ( DEBUG_VULC ) &
       write(6,'(a20/4i6/6i6/4i6)')'Volcan: check1 ',  &
       i,j, i_volc(volc_no),j_volc(volc_no),           &
       i_local(i),j_local(j), li0, li1, lj0, lj1,      &
       gi0,gi1,gj0,gj1

       if ( (i_local(i) >= li0) .and. (i_local(i) <= li1 )  .and.  &
            (j_local(j) >= lj0) .and. (j_local(j) <= lj1) ) then 

          unit_conv1 = GRAV* 0.001*AVOG/ &
                      (sigma_bnd(KMAX_BND-k+1) - sigma_bnd(KMAX_BND-k))

          rcemis_volc0(volc_no) = emis_volc(volc_no)  &
                                * unit_conv1 / molwt(QRCVOL)

         if ( DEBUG_VULC ) &
           write(*,*)'rc_emis_volc is ',rcemis_volc(volc_no)

       endif
    enddo ! volc_no
   end subroutine Set_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   subroutine Scale_Volc

 !-----------------------------------------------------------------------!
 ! Finishing converting volcano emissions to molecules/cm3/s 
 ! (every advection timestep)
 !-----------------------------------------------------------------------!

   integer i,j,k,i_l,j_l, i1,i2,j1,j2
   real unit_conv2 


  do volc_no=1,nvolc

     k=height_volc(volc_no)
     i=i_volc(volc_no) +IRUNBEG-1   !NEW
     j=j_volc(volc_no) +JRUNBEG-1   !NEW
    ! i=i_volc(volc_no)
    ! j=j_volc(volc_no)

     if ( DEBUG_VULC ) &
     write(6,'(a20/4i6/6i6/4i6)')'Volcan: check2 ', &
     i,j, i_volc(volc_no),j_volc(volc_no),           &
     i_local(i),j_local(j), li0, li1, lj0, lj1,      &
     gi0,gi1,gj0,gj1

       if ( (i_local(i) >= li0) .and. (i_local(i) <= li1 )  .and.  &
            (j_local(j) >= lj0) .and. (j_local(j) <= lj1) ) then 

        i_l = i_local(i) !local i
        j_l = j_local(j) !local j


        if ( DEBUG_VULC ) &
        write(6,'(a30,4i8)')'Volcan: check 3: ',   &
        i_l, j_l, i_volc(volc_no)-gi0+1, j_volc(volc_no)-gj0+1

        unit_conv2 = roa(i_l,j_l,KMAX_BND-k,1) / (ps(i_l,j_l,1)-PT)

        rcemis_volc(volc_no) = rcemis_volc0(volc_no) * unit_conv2

        if ( DEBUG_VULC ) &
           write(*,*)'rc_emis_volc is ',rcemis_volc(volc_no)

     endif
  enddo ! volc_no

  end subroutine Scale_Volc
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module Volcanos_ml
