! <AirEmis_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module AirEmis_ml
  ! NOx Emissions from Aircraft and Lightning.
  ! ANCAT emissions converted from Kg/month/gridcell to flux 
  ! in molecules cm-3 s-1 on ANCAT grid (2.8 X 2.8 degree). 
  ! Emissions on (finer) model grid then assigned from the ANCAT grid where
  ! model grid falls within. Note that aircraft emissions is given on a t42 
  ! and lightning on t21.
  !
  ! Ref: Gardner et.al., 1997, The ANCAT/EC global inventory of NOx emission
  !      from aircraft, Atm. Environ., 31(12), 1751-1766.
  !
  ! 
  ! Variable listing is given below
  !
  !  
   use Par_ml               , only : MAXLIMAX, MAXLJMAX, limax,ljmax, me
   use ModelConstants_ml    , only : KCHEMTOP, KMAX_MID, KMAX_BND, NPROC
   use Io_ml                , only : IO_AIRN, IO_LIGHT, ios, open_file
   use GridValues_ml        , only : gl,gb, GRIDWIDTH_M
   use PhysicalConstants_ml , only : AVOG
   use Met_ml               , only : z_bnd  
   use TimeDate_ml,           only : current_date

   implicit none
   private

   real, public, dimension(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), save :: &
                  airn                 & ! aircraft NOx emissions
                 ,airlig                 ! lightning NOx emissions

   public :: aircraft_nox !reads in the raw data
   public :: lightning

   private :: air_inter !interpolate the data into required grid
 
   include 'mpif.h'
   
   integer STATUS(MPI_STATUS_SIZE),INFO
   real MPIbuff
   integer,private ,parameter :: ILEV=18
   logical, parameter :: MY_DEBUG = .false.

 contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
  subroutine aircraft_nox(newseason)


!input
      integer,intent(in):: newseason

!local
      integer, parameter ::  ILON = 128, IGL = 32, GGL = 64
      real,    parameter ::  DLON = 4.21875-1.40625,RLON0 = -1.40625

      integer i, j, k, nlon, ngl, nlev, level,nlevb
      real zrmin, zfak, secmonth

! Definition of the ancat grid   ~ t42 :
! Data read in from N --> S  and from longitude 0  ( not from +-180 )  

      real, dimension(GGL) ::     ygrida       ! grid mid. pt. N-S
      real, dimension(IGL) ::     ygrdum   &   ! grid mid. pt. N-S
                                 ,area         ! grid area N-S
      real, dimension(ILON+1) ::   rlon        !  

      integer, dimension(ILON,GGL)         :: intnox ! global emission kg/month/cell  
      real,    dimension(ILON,GGL,-1:ILEV) :: flux   ! emission converted to flux
                                                     ! molecules/cm3/s 

      character*20 fname

      data ygrdum / 87.86379884, 85.09652699, 82.31291295, 79.52560657, &
                    76.73689968, 73.94751515, 71.15775201, 68.36775611, & 
                    65.57760700, 62.78735180, 59.99702011, 57.20663153, &
	            54.41619953, 51.62573360, 48.83524097, 46.04472663, &
	            43.25419467, 40.46364818, 37.67308960, 34.88252099, &
	            32.09194388, 29.30135962, 26.51076933, 23.72017390, &
	            20.92957425, 18.13897099, 15.34836476, 12.55775612, &
                     9.76714556,  6.97653355,  4.18592053,  1.39530/	 
	

      data area    / 3.4123374E+09,  8.2558925E+09,  1.2956985E+10, 	 &
		     1.7624316E+10,  2.2249359E+10,  2.6821497E+10,      &
		     3.1329966E+10,  3.5764105E+10,  4.0113402E+10,	 &
		     4.4367553E+10,  4.8516461E+10,  5.2550296E+10, 	 &
		     5.6459481E+10,  6.0234756E+10,  6.3867154E+10,	 &
		     6.7348070E+10,  7.0669238E+10,  7.3822790E+10, 	 &
		     7.6801245E+10,  7.9597535E+10,  8.2205024E+10, 	 &
		     8.4617535E+10,  8.6829343E+10,  8.8835195E+10, 	 &
		     9.0630349E+10,  9.2210528E+10,  9.3572006E+10, 	 &
		     9.4711529E+10,  9.5626412E+10,  9.6314483E+10,	 &
		     9.6774103E+10,  9.7004184E+10/


! ---- Defines the ANCAT grid ----------------------------------------------
!   WARNING!!!  This is not the correct t42 grid, but the aircraft 
!               emissions are defined in this grid
!
!-----Variable listing---------------
!
!MAXLIMAX      ==> Maximum no. of local points in longitude
!MAXLJMAX      ==> Maximum no. of local points in latitude
!limax         ==> Actual number of local points in longitude
!ljmax         ==> Actual number of local points in latitude
!NPROC         ==> Total no. of processors for parallel computation
!me            ==> Address of processer, host=0    (numbering starts at 0 
!                   in south-west corner of ground level
!KCHEMTOP      ==> Topmost level where chemistry is performed (k=2, chemistry is not done for k=1)
!KMAX_MID      ==> Number of levels in vertical (=20)
!KMAX_BND      ==> Number of levels in vertical + 1 (=KMAX_MID+1)
!current_date  ==> derived type containing date info
!IO_AIRN       ==> Input variable for Aircraft_nox emission from input emission file
!IO_LIGHT      ==> Input variable for lightning emission from input emission file
!ios           ==> I/O error status number
!open_file     ==> Checks that file exists and opens if required
!gl            ==> Geographical longitude of EMEP grid center
!gb            ==> Geographical latitude of EMEP grid  center 
!GRIDWIDTH_M   ==> Width of grid at 60N, in meters
!AVOG          ==> Avogadro's No.
!z_bnd         ==> Height of full layers
!GGL           ==> No. of latitudes (=64 on T42 grid and =32 on T21 grid)
!IGL           ==> No. of latitudes in NH (=32 for T42 and =16 for T21. Data read from N to S)
!ILON          ==> No. of longitudes (=128 for T42 and =64 for T21)
!DLON          ==> Delta longitude
!RLON0         ==> First longitude point
!


	secmonth = 3600.*24.*31.
	flux(:,:,:) = 0.



! --- Open and read ancat data (originally from DLR, EU project POLINAT)
! --- Commercial aircraft emissions every season
! --- Military aircraft emission read in as annual data

      if(me == 0)then

          write(fname,fmt='(''ancat'',i2.2,''.dat'')') newseason
          call open_file(IO_AIRN,"r",fname,needed=.true.,skip=1)
          if (ios /= 0)   WRITE(*,*) 'MPI_ABORT: ', "ioserror: ancat" 
          if (ios /= 0) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 

      end if ! me == 0


      if(me == 0)then

         read(IO_AIRN,'(3i4,2e22.13)') nlon,ngl,nlev,zrmin,zfak
 	 if (MY_DEBUG)write(6,*) nlon,ngl,nlev,zrmin,zfak
              do k = 0,NLEV-1
                 read(IO_AIRN,'(i2)') level
	         read(IO_AIRN,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
                 do i = 1,ngl
                    do j = 1,nlon
                       flux(j,i,k)=(float(intnox(j,i))*zfak)+zrmin
                    end do
                 end do
              end do

         close(IO_AIRN)

         call open_file(IO_AIRN,"r","ancatmil.dat",needed=.true.,skip=1)
         if (ios /= 0)   WRITE(*,*) 'MPI_ABORT: ', "ioserror: ancatmil" 
         if (ios /= 0) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 

      end if ! me == 0


      if(me == 0)then
 	 read(IO_AIRN,'(3i4,2e22.13)') nlon,ngl,nlevb,zrmin,zfak
	 if (MY_DEBUG) write(6,*) nlon,ngl,nlevb,zrmin,zfak
              do k = 1,nlevb
 	         read(IO_AIRN,'(i2)') level
	         read(IO_AIRN,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
                 do i = 1,ngl
                    do j = 1,nlon
	               flux(j,i,k)=flux(j,i,k)+(float(intnox(j,i))*zfak)+zrmin
                    end do
                 end do
              end do

         close(IO_AIRN)

      endif ! me == 0

      call air_inter(ILON   ,IGL    ,GGL  ,0     ,  &
                     flux   ,airn                ,  &
		     ygrdum ,ygrida ,DLON ,RLON0 ,  &
                     rlon   ,area   ,secmonth )
  
  end subroutine aircraft_nox

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine lightning()

      integer, parameter ::  ILON = 64, IGL = 16, GGL = 32
      real,    parameter ::  DLON =  8.4375 - 2.8125, RLON0 = -2.8125

      integer i, j, k, nlon, ngl, nlev ,level
      real    zrmin, zfak, secmonth, sumnox



! Definition of the ancat grid   ~ t21 :
! Data read in from N --> S  and from longitude 0  ( not from +-180 )  
! NB!!  note the difference between lightning and aircraft emission grid


      real, dimension(GGL)    ::     ygrida       ! grid mid. pt. N-S
      real, dimension(IGL)    ::     ygrdum   &   ! grid mid. pt. N-S
                                    ,area         ! grid area N-S
      real, dimension(ILON+1) ::     rlon         !  

      integer, dimension(ILON,GGL)         :: intnox ! global emission kg/month/cell  
      real,    dimension(ILON,GGL,-1:ILEV) :: flux   ! emission converted to flux
                                                     ! molecules/cm3/s 
            

      character*20 fname

      data ygrdum / 85.76058712, 80.26877907, 74.74454037,             &
	            69.21297617, 63.67863556, 58.14295405,	       &
		    52.60652603, 47.06964206, 41.53246125,	       &
		    35.99507841, 30.45755396, 24.91992863,	       &
		    19.38223135, 13.84448373,  8.30670286,	       & 
		     2.76890301/

      data area  /  268516310010.37,  64778953547.94, 101133318591.76, &
 		    136519495343.02, 170627187541.37, 203140714921.94, &
		    233757180440.66, 262190945894.60, 288176619318.18, &
	 	    311471618334.64, 331858463070.53, 349146817322.73, &
		    363175270173.37, 373812845114.06, 380960223957.41, &
		    384550674664.23/

! ---- Defines the ANCAT grid ----------------------------------------------



      secmonth = 1.
      flux(:,:,:) = 0.

! --- Read Emission data received from DLR 

      if(me == 0)then
	 sumnox = 0.
	 write(fname,fmt='(''lightn'',i2.2,''.dat'')') 	&
		current_date%month

! - open and read 1 line of header
  
         call open_file(IO_LIGHT,"r",fname,needed=.true.,skip=1)
             if (ios /= 0 )   WRITE(*,*) 'MPI_ABORT: ', "ioserror: lightning" 
                if (ios /= 0 ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
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

         write(6,*) 'lightning sum nox in AirEmis',sumnox

      endif

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


      integer, parameter  :: KMAX_BND_AIR = 21 
      integer, intent(in) :: ILON,IGL,GGL,iktop
      real, intent(in)    :: area(IGL), ygrdum(IGL),DLON,RLON0,secmonth


      real, intent(inout) ::   flux(ILON,GGL,-1:ILEV)

      real, dimension(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), intent(out) :: airem
      real, intent(out)   :: ygrida(GGL)
      real, intent(out)   :: rlon(ILON+1)

      !	local
      integer info
      integer lon,lat,i,j,ig,jg,kg,k, i_sh
      integer la_tst1, la_tst2, lo_tst1, lo_tst2   !  test area for sums
      real    height,     &  !  height of the emission levels
              atwno2,     &  !  atomic weight of NO2
              vol            !  volume of model grid boxes
      real    frac, above, below, glij
      real sum,sumnox,volcm,sum2
      integer, dimension(MAXLIMAX,MAXLJMAX) :: ixn  & !  mapping of emission 
                                              ,jxn    !  grid to model grid
      integer, dimension(KMAX_MID)          :: ilevel
      real    fraca(KMAX_MID), fracb(KMAX_MID)


      height = 1.e5
      atwno2 = 46.

      !  print out values on a sub-domain for comparison with direct model input
      !  NB!  due to different resolution the subdomain will be different for 
      !       aircraft and lightning emissions

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
                endif
        
                do lon=1,ILON
                   sumnox = sumnox + flux(lon,lat,k)
                   flux(lon,lat,k)=flux(lon,lat,k)*1.e3*AVOG &
                                     /volcm/secmonth/atwno2
                end do !lon
            end do    !lat 
         end do       !k
            
         write(6,*) 'SUMNOX, ANCAT:',sumnox
      endif		!me=0


      CALL MPI_BCAST(flux(1,1,iktop), 8*GGL*ILON*(ILEV+1-iktop), MPI_BYTE, 0,&
          MPI_COMM_WORLD, INFO)

      ! -- N/S
      ygrida(1) = 90.
      ygrida(GGL) = -90.

      do i=2,IGL
         i_sh = GGL + 1 - i
         ygrida(i) = (ygrdum(i-1)+ygrdum(i))*0.5
         ygrida(i_sh) = - ygrida(i)
      enddo

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
            if(abs(gb(i,j)-90.0)<0.001) then
               ixn(i,j) = 1
               jxn(i,j) = 1
            else

               do while(gb(i,j)<ygrida(jg+1))
                  jg = jg+1
               enddo

               do while(gb(i,j)>=ygrida(jg))
                  jg = jg-1
               enddo
 
               jxn(i,j) = jg
               glij = gl(i,j)
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
               if(above-below>1.1e3) 			      &
                  fraca(k) = (z_bnd(i,j,k)-(above-1.e3))        &
                               /(z_bnd(i,j,k) - z_bnd(i,j,k+1))
               fracb(k) = 0.
               if(above-below>2.1e3)			      &
                   fracb(k) = (below+1.e3 - z_bnd(i,j,k+1))      &
                              /(z_bnd(i,j,k) - z_bnd(i,j,k+1))
            end do  ! k

            lon = ixn(i,j)
            lat = jxn(i,j)

            do k = KCHEMTOP,KMAX_MID
                 frac = 1. - fraca(k) - fracb(k)
                 kg = ilevel(k)
                 airem(k,i,j) = flux(lon,lat,kg)*fraca(k)             &
                              + flux(lon,lat,kg-1)*frac		      &
                              + flux(lon,lat,kg-2)*fracb(k)
            end do

           ! surface emissions

             if(iktop == 0)		&
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

            if(gl(i,j)>rlon(lo_tst1) .and. gl(i,j)<rlon(lo_tst2+1) .and.     &
               gb(i,j)<ygrida(la_tst1) .and. gb(i,j)>ygrida(la_tst2+1)) then
                do k=KCHEMTOP,KMAX_MID
                   vol = GRIDWIDTH_M*GRIDWIDTH_M                             &
                        *(z_bnd(i,j,k)-z_bnd(i,j,k+1))*1.e6
                   sum2 = sum2 + airem(k,i,j)*atwno2/AVOG*		     &
                          vol*secmonth*1.e-3
                end do
            end if
         end do
      end do

        MPIbuff=sum2
        CALL MPI_ALLREDUCE(MPIbuff,sum2, 1, &
        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      if(me == 0) write(6,*) 'ancat on limited area:',sum,sum2

  end subroutine air_inter
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module AirEmis_ml

