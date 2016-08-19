! <MassBudget_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 module   MassBudget_ml

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! DESCRIPTION 
! Routine to cross check the mass balance of the model
! 29/10/02 - output formatting and descriptions improved, ds.
! 1/10/01 - code for derived fields removed. MY_MASS_PRINT ADDED, ds
! Oct, 2001 - ** new ** mass budget method by jej
! Nov.2001, ds, "use" statements moved to top, much code moved to REMOVED 
! section at end
!_____________________________________________________________________________

!! use DryDep_ml, only : NDRYDEP_ADV, DDepMap , DryDep_Budget

 use ChemChemicals_ml, only : species       ! species identifier
 use ChemSpecs_tot_ml,  only : NSPEC_TOT     ! No. species (long-lived)
 use ChemSpecs_adv_ml,  only : NSPEC_ADV     ! No. species (long-lived)
 use ChemSpecs_shl_ml,  only : NSPEC_SHL     ! No. species (shorshort-lived)
 use Chemfields_ml ,  only : xn_adv        ! advective flag
 use GridValues_ml ,  only : carea,xmd     ! cell area, 1/xm2 where xm2 is 
                                           ! the area factor in the middle 
                                           ! of the cell 
 use Io_ml         ,  only : IO_RES        ! =25
 use MetFields_ml  ,  only : ps            ! surface pressure  
 use ModelConstants_ml,                 &
                      only : KMAX_MID   &  ! Number of levels in vertical
                            ,MasterProc &  ! Master processor
                            ,NPROC      &  ! No. processors
                            ,PT         &  ! Pressure at top
                            ,ATWAIR     &  ! Mol. weight of air(Jones,1992)
                            ,TXTLEN_NAME&
                            ,EXTENDEDMASSBUDGET
 use Par_ml,          only : MAXLIMAX   & 
                            ,MAXLJMAX   &  
                            ,li0,li1    &
                            ,lj0,lj1    &
                            ,limax,ljmax&
                            ,gi0, gj0   &
                            ,GIMAX,GJMAX
 use Setup_1dfields_ml, only : amk ! Air concentrations 

!Variable listing
!MAXLIMAX    ==> Maximum number of local points in longitude
!MAXLJMAX    ==> Maximum number of local points in latitude
!li0         ==> First local index in longitude when 
!                outer boundary is excluded
!li1         ==> Last local index in longitude when 
!                outer boundary is excluded
!lj0         ==> First local index in latitude when 
!                outer boundary is excluded
!lj1         ==> Last local index in latitude when 
!                outer boundary is excluded
!NPROC       ==> Total no. of processors for parallel computation
!limax       ==> Actual number of local points in longitude
!ljmax       ==> Actual number of local points in latitude
!me          ==> Address of processer, host=0 (numbering starts at 0
!                in south-west corner of ground level
!gi0         ==> Global address of longitude start point
!gj0         ==> Global address of latitute start point
!GIMAX = 132 ==> Number of global points in longitude
!GJMAX = 111 ==> Number of global points in latitude

 
implicit none
private
   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE),INFO
   real    MPIbuff(NSPEC_ADV*KMAX_MID)

! Some work arrays used in Aqueous_ml and (in future) DryDry
! Use tot index for convenience
  real, public, save, dimension(NSPEC_TOT) ::   &
      wdeploss, & 
      ddeploss

! The following parameters are used to check the global mass budget:
! Initialise here also.

  real, public, save, dimension(NSPEC_ADV) ::   &
      sumint   = 0.0   & !  initial mass
     ,fluxin   = 0.0   & !  mass in  across lateral boundaries
     ,fluxout  = 0.0   & !  mass out across lateral boundaries
     ,totddep  = 0.0   & !  total dry dep
     ,totwdep  = 0.0   & !  total wet dep
     ,totem    = 0.0   & !  total emissions
     ,totox    = 0.0   & !  total oxidation
     ,totldep  = 0.0     !  local deposition (not in use - Lagrangian)

  real, public, save, dimension(NSPEC_ADV) ::  &
      amax = -2.0   &  ! maximum concentration in field -2
     ,amin =  2.0      ! minimum concentration in field  2

  public :: Init_massbudget
  public :: massbudget
!  public :: DryDep_Budget

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine Init_massbudget()
    ! Initialise mass-budget - calculate mass of concentrations fields
    ! within 3-D grid, after boundary conditions
    !
    !----------------------------------------------------------------------

    integer i, j, k, n, info    ! lon,lat,lev indexes
                                ! n - No. of species
                                ! info - printing info
    real rwork

    do k=2,KMAX_MID   
      do j=lj0,lj1
        do i=li0,li1
            rwork = carea(k)* xmd(i,j)*(ps(i,j,1) - PT)
            sumint(:) = sumint(:) + xn_adv(:,i,j,k)*rwork  ! sumint in kg      
        enddo
      enddo
    enddo

      MPIbuff(1:NSPEC_ADV)= sumint (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, sumint , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 

    if(MasterProc.and.EXTENDEDMASSBUDGET)then
         do n = 1,NSPEC_ADV
       if(sumint(n) >  0. ) then
             write(IO_RES,"(a15,i4,4x,e10.3)") "Initial mass",n,sumint(n) 
             write(6,"(a15,i4,4x,e10.3)") "Initial mass",n,sumint(n) 
           end if
         enddo
    end if

 end subroutine Init_massbudget


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine massbudget()

  ! Converted from old masbud.f by ds, March 2001
  ! sums over all sulphur and nitrogen, so is model independant.


  integer ::  i, j, k, n, nn, info               ! lon,lat,lev indexes
                                                 ! n - No. of species
                                                 ! nn - Total no. of short
                                                 ! lived and advected species
                                                 ! info - printing info
  integer :: ifam                                ! family index
  real, dimension(NSPEC_ADV,KMAX_MID) ::  sumk   ! total mass in each layer
  integer, parameter :: NFAMILIES = 3            ! No. of families         
  character(len=8), dimension(NFAMILIES), save :: family_name = &
           (/ "Sulphur ", "Nitrogen", "Carbon  " /)

  real, dimension(NFAMILIES) ::family_init  & ! initial total mass of 
                                              ! species family
                              ,family_mass &  ! total family mass at the 
                                              ! end of the model run
                              ,family_inflow &! total family mass flowing in  
                              ,family_outflow&! total family mass flowing out 
                              ,family_ddep&   ! total family mass dry dep. 
                              ,family_wdep&   ! total family mass wet dep. 
                              ,family_em  &   ! total family mass emitted 
                              ,family_input & ! total family mass input
                              ,family_fracmass  ! mass fraction (should be 1.0)

  real, dimension(NSPEC_ADV) :: &
        xmax, xmin, & ! min and max value for the individual species
        sum_mass,   & ! total mass of species
        frac_mass,  & ! mass budget fraction (should=1) for groups of species
        gfluxin, gfluxout,  & ! flux in  and out
        gtotem,    & ! total emission
        gtotddep, gtotwdep, & ! total dry and wet deposition
        gtotldep, & ! local dry deposition
        gtotox      ! oxidation of SO2 

  real :: totdiv,helsum, natoms

    sum_mass(:)     = 0.
    frac_mass(:)    = 0.
    xmax(:)         = -2.
    xmin (:)        = 2.
    gfluxin(:)      = fluxin(:)
    gfluxout(:)     = fluxout(:)
    gtotem(:)       = totem(:)
    gtotddep(:)     = totddep(:)
    gtotwdep(:)     = totwdep(:)
    gtotldep(:)     = totldep(:)
    gtotox(:)       = totox(:)


    sumk(:,:) = 0.

    do k = 1,KMAX_MID
      do j = lj0,lj1
        do i = li0,li1

            helsum  = carea(k)*xmd(i,j) * (ps(i,j,1) - PT)

            xmax(:) = amax1(xmax(:),xn_adv(:,i,j,k))
            xmin(:) = amin1(xmin(:),xn_adv(:,i,j,k))

            sumk(:,k) = sumk(:,k) + xn_adv(:,i,j,k)*helsum

        enddo
      enddo
    enddo



      MPIbuff(1:NSPEC_ADV)= xmax(1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, xmax, NSPEC_ADV,&
      MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= xmin   (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, xmin   , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gfluxin (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gfluxin , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gfluxout (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gfluxout , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gtotem (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gtotem , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gtotddep (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gtotddep , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gtotwdep (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gtotwdep , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gtotldep (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gtotldep , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      MPIbuff(1:NSPEC_ADV)= gtotox (1:NSPEC_ADV) 
      CALL MPI_ALLREDUCE(MPIbuff, gtotox , NSPEC_ADV, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
      j=0
      do k=1,KMAX_MID
         do i=1,NSPEC_ADV
            j=j+1
            MPIbuff(j)= sumk (i,k)
         enddo
      enddo
      CALL MPI_ALLREDUCE(MPIbuff, sumk , NSPEC_ADV*KMAX_MID, &
      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 

!     make some temporary variables used to hold the sum over all 
!     domains. Remember that sumint already holds the sum over all
!     domains, see inass
!
    amax(:) = max( amax(:), xmax(:) )
    amin(:) = min( amin(:), xmin(:) )
    do k = 2,KMAX_MID
      sum_mass(:) = sum_mass(:)+sumk(:,k)
    enddo



    do n = 1,NSPEC_ADV

      totdiv = sumint(n) + gtotem(n) + gfluxin(n)
      frac_mass(n) = sum_mass(n)  + (gtotddep(n)+gtotwdep(n))*ATWAIR &
                   + gfluxout(n) 

      if(totdiv >  0.0 ) frac_mass(n) = frac_mass(n)/totdiv


    end do


   if ( MasterProc ) then   ! printout from node 0

    do n = 1,NSPEC_ADV
      if (gtotem(n) > 0.0 .and. EXTENDEDMASSBUDGET) write(6,*)   &
                           'tot. emission of species ',n,gtotem(n)
    end do

    family_init(:)  = 0.
    family_mass(:) = 0.
    family_inflow(:)  = 0.
    family_outflow(:)  = 0.
    family_input(:)  = 0.
    family_fracmass(:) = 0.
    family_ddep(:) = 0.
    family_wdep(:) = 0.
    family_em(:)   = 0.
    natoms = 0.0

    write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'      

     do ifam = 1, 3

       write(6,"(a8,i3,a12)") 'family ', ifam,  family_name(ifam) 
       do n = 1, NSPEC_ADV

         nn = NSPEC_SHL + n


           if ( ifam == 1 ) natoms = real(species(nn)%sulphurs)

           if ( ifam == 2 ) natoms = real(species(nn)%nitrogens)

           if ( ifam == 3 ) natoms = real(species(nn)%carbons)

           if (natoms > 0) then
             family_init(ifam) = family_init(ifam) + sumint(n)*natoms
             family_mass(ifam) = family_mass(ifam) + sum_mass(n)*natoms
             family_inflow(ifam) = family_inflow(ifam) + gfluxin(n)*natoms
             family_outflow(ifam) = family_outflow(ifam) + gfluxout(n)*natoms
             family_ddep(ifam) = family_ddep(ifam) + gtotddep(n)*natoms
             family_wdep(ifam) = family_wdep(ifam) + gtotwdep(n)*natoms
             family_em(ifam) = family_em(ifam) + gtotem(n)*natoms
           end if
       end do  ! NSPEC_ADV

       family_input(ifam) = family_init(ifam) &
                         + family_inflow(ifam) &
                         + family_em(ifam)   

       if (family_input(ifam) > 0.0 ) &
              family_fracmass(ifam) = (family_mass(ifam) &
                                 +  family_outflow(ifam)  &
                                 +  family_ddep(ifam)*ATWAIR  & 
                                 +  family_wdep(ifam)*ATWAIR) & 
                                 / family_input(ifam)


      write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'      
      write(6,*)

      write(6,"(a9,5a12)") "family", "sumint", "summas", &
                               "fluxout","fluxin", "fracmass"
      write(6,"(a9,5es12.4)") family_name(ifam), &
             family_init(ifam), family_mass(ifam),family_outflow(ifam), &
             family_inflow(ifam), family_fracmass(ifam)

      write(6,*)
      write(6,"(a9,3a14)") "ifam", "totddep","totwdep","totem"
      write(6,"(i9,3es14.3)") ifam, family_ddep(ifam)*ATWAIR  &
                        , family_wdep(ifam)*ATWAIR  &
                        , family_em(ifam) 
      write(6,*)
      write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'

    end do  ! ifam = 1,3
    

   end if

  if ( MasterProc .and. EXTENDEDMASSBUDGET) then     ! printout from node 0

     !/.. now use species array which is set in My_MassBudget_ml
     do n = 1,NSPEC_ADV
        write(6,*)
        write(IO_RES,*)
        do k = 1,KMAX_MID
           write(6,950)      n,species(n+NSPEC_SHL)%name, k,sumk(n,k)
           write(IO_RES,950) n,species(n+NSPEC_SHL)%name, k,sumk(n,k)
        end do
     enddo
950  format(' Spec ',i3,2x,a12,5x,'k= ',i2,5x,es12.5)
     
     do n = 1,NSPEC_ADV

        write(6,*)
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'      
        write(6,*)

        write(6,"(a3,6a12)") " n ", "Spec", "sumint", "summas", &
                               "fluxout","fluxin", "fracmass"
        write(6,"(i3,1x,a11,5es12.4)") n,species(n+NSPEC_SHL)%name, &
                  sumint(n),sum_mass(n), gfluxout(n),gfluxin(n), frac_mass(n)

        write(6,*)
        write(6,"(a3,6a12)")  " n ", "species", &
                      "totox", "totddep", "totwdep", "totem", "totldep"
        write(6,"(i3,1x,a11,5es12.4)") n, species(n+NSPEC_SHL)%name, &
                 gtotox(n), gtotddep(n), gtotwdep(n), gtotem(n), gtotldep(n)
        write(6,*)
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
                 
    enddo
!              

   end if  ! MasterProc

 end subroutine massbudget



!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!  subroutine DryDep_Budget(i,j,Loss,convfac)
!     !use ChemSpecs_adv_ml
!
!      real, dimension(NSPEC_ADV), intent(in) :: Loss
!      real, dimension(NSPEC_ADV)             :: DryLoss
! 
!     real, intent(in)  :: convfac
!     integer           :: n,nadv,i,j  ! index in IXADV_  arrays
!
!     DryLoss(:)=Loss(:)* convfac /amk(KMAX_MID)   !molec/cm3->mix ratio 
!
!      do n = 1, NDRYDEP_ADV 
!         nadv    = DDepMap(n)%ind
!         totddep( nadv ) = totddep (nadv) + DryLoss(nadv) 
!
!      enddo
!  end subroutine DryDep_Budget
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 end module MassBudget_ml
!--------------------------------------------------------------------------
