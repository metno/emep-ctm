! <AirEmis_mod.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
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
module AirEmis_mod

  ! NOx emissions from lightning
  ! Emissions converted from Kg/month/gridcell to flux 
  ! in molecules cm-3 s-1 on T21 (5.65x5.65deg) resolution 
  ! 
  ! Variable listing is given below
  
    
   use Io_mod                , only : IO_AIRN, IO_LIGHT, ios, open_file
   use GridValues_mod        , only : glon,glat, GRIDWIDTH_M
   use MetFields_mod         , only : z_bnd  
   use Config_module    ,      only : KCHEMTOP, KMAX_MID, KMAX_BND, NPROC, &
                                     USES, lightningFile
   use MPI_Groups_mod,         only : MPI_BYTE, MPI_DOUBLE_PRECISION, MPI_REAL8, MPI_INTEGER, MPI_LOGICAL, &
                             MPI_MIN, MPI_MAX, MPI_SUM, &
                             MPI_COMM_CALC, MPI_COMM_WORLD, MPISTATUS, IERROR, ME_MPI, NPROC_MPI
   use OwnDataTypes_mod,        only: TXTLEN_FILE
   use Par_mod               , only : LIMAX, LJMAX, limax,ljmax, me
   use PhysicalConstants_mod , only : AVOG
   use SmallUtils_mod,         only : key2str
   use TimeDate_mod,           only : current_date

   implicit none
   private

   real, public, dimension(:,:,:), save,allocatable :: &
                  airn                 & ! aircraft  NOx emissions
                 ,airlig                 ! lightning NOx emissions

   public :: lightning

   private :: air_inter !interpolate the data into required grid
 
   integer,private ,parameter :: ILEV=18
   logical, parameter :: MY_DEBUG = .false.

 contains
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


    subroutine lightning()

      integer, parameter ::  ILON = 64, IGL = 16, GGL = 32
      real,    parameter ::  DLON =  8.4375 - 2.8125, RLON0 = -2.8125

      integer i, j, k, nlon, ngl, nlev ,level
      real    zrmin, zfak, secmonth, sumnox



! Definition of the grid   ~ t21 :
! Data read in from N --> S  and from longitude 0  ( not from +-180 )  

      real, dimension(GGL)    ::     ygrida       ! grid mid. pt. N-S
      real, dimension(IGL)    ::     ygrdum   &   ! grid mid. pt. N-S
                                    ,area         ! grid area N-S
      real, dimension(ILON+1) ::     rlon         !  

      integer, dimension(ILON,GGL)         :: intnox ! global emission 
                                                     ! kg/month/cell  
      real,    dimension(ILON,GGL,-1:ILEV) :: flux   ! emission converted
                                                     ! to flux
                                                     ! molecules/cm3/s 
            

      character(len=TXTLEN_FILE) :: fname

      data ygrdum / 85.76058712, 80.26877907, 74.74454037, &
                69.21297617, 63.67863556, 58.14295405, &
            52.60652603, 47.06964206, 41.53246125, &
            35.99507841, 30.45755396, 24.91992863, &
            19.38223135, 13.84448373,  8.30670286, & 
             2.76890301/

      data area  /  268516310010.37,  64778953547.94, 101133318591.76, &
             136519495343.02, 170627187541.37, 203140714921.94, &
            233757180440.66, 262190945894.60, 288176619318.18, &
             311471618334.64, 331858463070.53, 349146817322.73, &
            363175270173.37, 373812845114.06, 380960223957.41, &
            384550674664.23/

! ---- Defines the grid ----------------------------------------------



      secmonth = 1.
      flux(:,:,:) = 0.

      
      if(.not.allocated(airlig))then
         allocate(airlig(KCHEMTOP:KMAX_MID,LIMAX,LJMAX))
         airlig=0.0
      end if

! --- Read Emission data received from DLR 

     if(me == 0)then
     sumnox = 0.
!     write(fname,fmt='(''lightning'',i2.2,''.dat'')')     &
!        current_date%month
     fname = key2str(lightningFile, 'MM', current_date%month)

! - open and read 1 line of header
  
         call open_file(IO_LIGHT,"r",fname,needed=.true.,skip=1)
             if (ios /= 0 )   WRITE(*,*) 'MPI_ABORT: ', "ioserror: lightning" 
                if (ios /= 0 ) call  MPI_ABORT(MPI_COMM_CALC,9,IERROR) 
             end if

      if(me == 0)then
     read(IO_LIGHT,'(3i4,2e22.13)') nlon,ngl,nlev,zrmin,zfak
     if (MY_DEBUG) write(6,*) nlon,ngl,nlev,zrmin,zfak
         do k = 1,nlev
         read(IO_LIGHT,'(i2)') level
         read(IO_LIGHT,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
             do i = 1,ngl
                do j = 1,nlon
               flux(j,i,k)=(float(intnox(j,i))*zfak)+zrmin
                     sumnox = sumnox + flux(j,i,k)
                end do
             end do
         end do

         close(IO_LIGHT)

         write(6,*) 'Sum of NOx emissions from lightning: ',sumnox

      end if

      call air_inter(ILON   ,IGL    ,GGL  ,1     ,      & 
                     flux   ,airlig              ,      &
             ygrdum ,ygrida ,DLON ,RLON0 ,      &
             rlon   ,area   ,secmonth)


  end subroutine lightning

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine air_inter(ILON   ,IGL    ,GGL    ,iktop ,  &
                       flux   ,airem  ,                 &
                       ygrdum ,ygrida ,DLON   ,RLON0 ,  &
                       rlon   ,area   ,secmonth  )


      integer, intent(in) :: ILON,IGL,GGL,iktop
      real, intent(in)    :: area(IGL), ygrdum(IGL),DLON,RLON0,secmonth


      real, intent(inout) ::   flux(ILON,GGL,-1:ILEV)

      real, dimension(KCHEMTOP:KMAX_MID,LIMAX,LJMAX), intent(out) :: airem
      real, intent(out)   :: ygrida(GGL)
      real, intent(out)   :: rlon(ILON+1)

      !    local
      integer lon,lat,i,j,ig,jg,kg,k, i_sh
      integer la_tst1, la_tst2, lo_tst1, lo_tst2   !  test area for sums
      real    height,     &  !  height of the emission levels
              atwno2,     &  !  atomic weight of NO2
              vol            !  volume of model grid boxes
      real    frac, above, below, glij
      real sum,sumnox,volcm,sum2,sum2_out
      integer, dimension(LIMAX,LJMAX) :: ixn  & !  mapping of emission 
                                              ,jxn    !  grid to model grid
      integer, dimension(KMAX_MID)          :: ilevel
      real    fraca(KMAX_MID), fracb(KMAX_MID)


      height = 1.e5
      atwno2 = 46.

      !  print out values on a sub-domain for comparison with direct model input

      la_tst1 = 7
      la_tst2 = 13
      lo_tst1 = 1
      lo_tst2 = 5

      if(me == 0)then
         sum = 0.
         sumnox = 0.
 
         do k = iktop,ILEV
            do lat = la_tst1,la_tst2
               do lon = lo_tst1,lo_tst2
                  sum = sum + flux(lon,lat,k)
               end do
            end do

            do lat=1,GGL
                if(lat<=IGL)then
                   volcm = area(lat)*1.e4*height
                else
      ! -- area not defined for Southern Hemisphere
                   volcm = area(GGL-lat+1)*1.e4*height
                end if
        
                do lon=1,ILON
                   sumnox = sumnox + flux(lon,lat,k)
                   flux(lon,lat,k)=flux(lon,lat,k)*1.e3*AVOG &
                                     /volcm/secmonth/atwno2
                end do !lon
            end do    !lat 
         end do       !k
            
         if(MY_DEBUG)write(6,*) 'SUMNOX, ANCAT:',sumnox
      end if        !me=0


      CALL MPI_BCAST(flux(1,1,iktop), 8*GGL*ILON*(ILEV+1-iktop), MPI_BYTE, 0,&
          MPI_COMM_CALC, IERROR)

      ! -- N/S
      ygrida(1) = 90.
      ygrida(GGL) = -90.

      do i=2,IGL
         i_sh = GGL + 1 - i
         ygrida(i) = (ygrdum(i-1)+ygrdum(i))*0.5
         ygrida(i_sh) = - ygrida(i)
      end do

      ! -  E/W
      rlon(1) = RLON0

      do i=2,ILON
         rlon(i) = rlon(i-1) + DLON
      end do
      rlon(ILON+1) = rlon(1)+360.

      ! -- Assign gridpoints to the EMEP grid


      jg = GGL-1
      do j = 1,ljmax
         do i = 1,limax
            if(abs(glat(i,j)-90.0)<0.001) then
               ixn(i,j) = 1
               jxn(i,j) = 1
            else

               do while(glat(i,j)<ygrida(jg+1))
                  jg = jg+1
               end do

               do while(glat(i,j)>=ygrida(jg))
                  jg = jg-1
               end do
 
               jxn(i,j) = jg
               glij = glon(i,j)
               if(glij<=rlon(1)) glij = glij+360.
               ig = int((glij-rlon(1))/DLON)+1
               if(ig>ILON) ig=ig-ILON 
               ixn(i,j) = ig
            end if
         end do   !i
      end do      !j


      do j = 1,ljmax
         do i = 1,limax

         kg = 1
        above = 1.e3
        below = 0.

            do k = KMAX_MID,KCHEMTOP,-1

               do while(z_bnd(i,j,k+1)>below+1.e3) 
                  below = below+1.e3
               end do

               do while (z_bnd(i,j,k)>above) 
                   kg = kg+1
                   above = above+1.e3
               end do

               ilevel(k) = kg
               fraca(k) = 1.
               if(above-below>1.1e3)   &
                  fraca(k) = (z_bnd(i,j,k)-(above-1.e3))        &
                               /(z_bnd(i,j,k) - z_bnd(i,j,k+1))
               fracb(k) = 0.
               if(above-below>2.1e3)   &
                   fracb(k) = (below+1.e3 - z_bnd(i,j,k+1))      &
                              /(z_bnd(i,j,k) - z_bnd(i,j,k+1))
            end do  ! k

            lon = ixn(i,j)
            lat = jxn(i,j)

            do k = KCHEMTOP,KMAX_MID
                 frac = 1. - fraca(k) - fracb(k)
                 kg = ilevel(k)
                 if(kg<=ILEV)then
                    airem(k,i,j) = flux(lon,lat,kg)*fraca(k)             &
                         + flux(lon,lat,kg-1)*frac              &
                         + flux(lon,lat,kg-2)*fracb(k)
                 else
                    !zero emissions
                    airem(k,i,j) = 0.0
                 end if
            end do

           ! surface emissions

             if(iktop == 0)        &
               airem(KMAX_MID,i,j) = airem(KMAX_MID,i,j)+flux(lon,lat,0)
         end do
      end do


      !! Print out on a limited part of the domain both raw data ( flux ) and 
      !! the re-gridded emissions ( airem ).  Expect only an approximate match! 
      !! Summation of aircraft emission and lightning are on different grids 
      !! and hence for different domains. 

      sum2 = 0.
      do j = 1,ljmax
         do i = 1,limax

            if(glon(i,j)>rlon(lo_tst1) .and. glon(i,j)<rlon(lo_tst2+1) .and.     &
               glat(i,j)<ygrida(la_tst1) .and. glat(i,j)>ygrida(la_tst2+1)) then
                do k=KCHEMTOP,KMAX_MID
                   vol = GRIDWIDTH_M*GRIDWIDTH_M                             &
                        *(z_bnd(i,j,k)-z_bnd(i,j,k+1))*1.e6
                   sum2 = sum2 + airem(k,i,j)*atwno2/AVOG*             &
                          vol*secmonth*1.e-3
                end do
            end if
         end do
      end do

        CALL MPI_ALLREDUCE(sum2, sum2_out,1, &
        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_CALC, IERROR) 
      if(me == 0.and.MY_DEBUG) write(6,*) 'ancat on limited area:',sum,sum2_out

  end subroutine air_inter
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module AirEmis_mod

