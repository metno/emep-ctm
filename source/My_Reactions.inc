! <My_Reactions.inc - A component of the EMEP MSC-W Unified Eulerian
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


!    do iter = 1, toiter(k)
!
! The chemistry is iterated several times, more close to the ground
! than aloft.
! For some reason, it proved faster to include files as given below
! than to use loops. 
if(k>=KCHEMTOP)then

   include 'My_FastReactions.inc'

endif
if(k>=6)then

   include 'My_FastReactions.inc'

endif

if(k>=KEMISTOP)then

   include 'My_FastReactions.inc'

endif
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   end do !! End iterations
 ! Just before SO4
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!-> SO4 

      P =   &
      2e-12*aqrck(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
       + aqrck(ICLRC1,K)*XNEW(H2O2) * xnew(SO2 )  &
       + aqrck(ICLRC2,K)*XNEW(O3) * xnew(SO2 )  &
       + aqrck(ICLRC3,K) * xnew(SO2 )  &
      + 0.05*rcemis(QRCSO2,K) 
      ! L = 0.0
 

       xnew(SO4)=  amax1(0.0,  xold(SO4) + dt2 * P)


!-> pNO3

      P =   &
        rcmisc(10,k) * xnew(HNO3)
      L = 0.0


       xnew(pNO3)=  amax1(0.0,  xold(pNO3) + dt2 * P)


!-> NH3 

      P =   &
     rcemis(QRCNH3,K) 
      ! L = 0.0
 

       xnew(NH3)=  amax1(0.0,  xold(NH3) + dt2 * P)

!-> PM25

      P =   &
     rcemis(QRCPM25,K) 
      ! L = 0.0
       xnew(PM25)=  amax1(0.0,  xold(PM25) + dt2 * P)


!-> PMCO

      P =   &
     rcemis(QRCPMCO,K) 
      ! L = 0.0
       xnew(PMCO)=  amax1(0.0,  xold(PMCO) + dt2 * P)

!-> Rn222

      P =  rc_Rn222(k)
      L = 1.0/(5.51614*24.0*3600.0)
       xnew(Rn222)=  amax1(0.0, ( xold(Rn222) + dt2 * P))  &
                     /(1.0 + dt2*L ) 

!-> Pb210
      P = 1.0/(5.51614*24.0*3600.0) * xnew(Rn222)
      ! L = 0.0
       xnew(Pb210)=  max(0.0,  xold(Pb210) + dt2 * P)

!SeaS
   if (SEASALT) then

!-> Sea salt fine

      P =   &
     rcss(QSSFI,K) 
      ! L = 0.0
       xnew(SSfi)=  amax1(0.0,  xold(SSfi) + dt2 * P)

!-> Sea salt coarse

      P =   &
     rcss(QSSCO,K) 
      ! L = 0.0
       xnew(SSco)=  amax1(0.0,  xold(SSco) + dt2 * P)
   endif

!-> AMSU 
      ! P = 0.0
 
      ! L = 0.0
 

!       xnew(AMSU)=  amax1(0.0, xold(AMSU))


!-> AMNI 
      ! P = 0.0
 
      ! L = 0.0
 

!       xnew(AMNI)=  amax1(0.0, xold(AMNI))

