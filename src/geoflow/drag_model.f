C*******************************************************************
C* Copyright (C) 2003 University at Buffalo
C*
C* This software can be redistributed free of charge.  See COPYING
C* file in the top distribution directory for more details.
C*
C* This software is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C*
C* Author: 
C* Description: 
C*
C*******************************************************************
C* $Id:$
C*

C***********************************************************************
      subroutine calc_drag_force (uvec, vsolid, vfluid, den_solid,
     &             den_fluid, vterminal, drag, tiny)
C***********************************************************************

      implicit none
      integer i, j
      double precision uvec(6), vsolid(2), vfluid(2)
      double precision den_fluid, den_solid
      double precision drag(4), vterminal, tiny

!     function calls
      double precision sgn

!     local variables and constants
      double precision exponant, temp, volf, denfrac
      double precision delv(2)
      parameter (exponant = 3.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     model in pitman-le paper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      volf = uvec(2)/uvec(1)
!      if (volf.gt.0.85) volf=0.85
      temp = uvec(1)*volf*(1.-volf)**(-exponant)
      denfrac = den_fluid/den_solid

       
      do 10 i=1,4
        drag(i)=0.
10    continue
 
      ! fluid vel - solid vel
      if ( uvec(1).gt.tiny ) then
        do 20 j=1,2
          delv(j) = (vfluid(j)-vsolid(j))/vterminal
20      continue

!       cap (u-v)/Vt to 0.1 
!        do 30 j=1,2
!          if (dabs(delv(j)).gt.0.1) then
!            delv(j)=0.1*sgn(delv(i), tiny)
!          endif
!30      continue

!       compute individual drag-forces
        drag(1) = (1.-denfrac)*temp*(1.-volf)*delv(1)
        drag(2) = (1.-denfrac)*temp*(1.-volf)*delv(2)
        drag(3) = (1.-denfrac)*temp*delv(1)/denfrac
        drag(4) = (1.-denfrac)*temp*delv(2)/denfrac
      endif

      end subroutine calc_drag_force

C***************************************************************************
C*
C*
C*
C*
C***************************************************************************
      subroutine calc_drag_force2 (uvec, vsolid, vfluid, den_solid,
     &             den_fluid, avgdia, viscosity, drag, tiny)
C***************************************************************************

      implicit none
      integer i, j

!     declare function arguments
      double precision uvec(6), vsolid(2), vfluid(2)
      double precision den_fluid, den_solid, avgdia, viscosity
      double precision drag(4), tiny

!     local variables
      double precision coeff, re, vf, vs
      double precision beta, kterm
      double precision relvel(2), relvmag, denfrac

      do 50 i=1, 4
50      drag(i)=0.

      if ( uvec(1) .gt. tiny ) then

!       volume fractions
        vs = uvec(2)/uvec(1)
        vf = 1.0 - vs

!       relative velocities
        relvmag=0.
        do 60 i=1,2
          relvel(i) = vfluid(i) - vsolid(i)
60        relvmag = relvmag + relvel(i)**2
        relvmag = dsqrt(relvmag)

!       reynold's number
        re = den_fluid*vf*avgdia*relvmag/viscosity

!       drag coeff
        coeff = (0.63 + 4.8/dsqrt(re))**2

!       compute beta
        beta = 3.7 - 0.65*dexp(-0.5*(1.5 - dlog(re))**2)

!       calculate K
        kterm = 0.75*coeff*den_fluid*relvmag*vs*(vf**(2-beta))/avgdia

!       finally the drag
        drag(1) = (1.-denfrac)*kterm*relvel(1)
        drag(2) = (1.-denfrac)*kterm*relvel(2)
        drag(3) = (1.-denfrac)*kterm*relvel(1)/denfrac
        drag(4) = (1.-denfrac)*kterm*relvel(2)/denfrac
      endif

      end subroutine calc_drag_force2   
