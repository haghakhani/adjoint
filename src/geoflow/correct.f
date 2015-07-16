C*******************************************************************
C     * Copyright (C) 2003 University at Buffalo
C     *
C     * This software can be redistributed free of charge.  See COPYING
C     * file in the top distribution directory for more details.
C     *
C     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C     *
C     * Author: 
C     * Description: 
C     *
C*******************************************************************
C     * $Id: correct.f 143 2007-06-25 17:58:08Z dkumar $ 
C     *

C***********************************************************************
      subroutine correct(uvec, uprev, fluxxp, fluxyp, fluxxm, fluxym,
     1     tiny, dtdx, dtdy, dt, dUdx, dUdy, xslope, yslope, curv,
     2     intfrictang, bedfrictang, g, dgdx ,kactxy, frict_tiny,
     3     forceint,forcebed, dragfoce ,DO_EROSION, eroded, v_solid,
     4     terminal_vel, eps, IF_STOPPED,
     5     fluxsrc)
C***********************************************************************

      implicit none
      double precision forceint, forcebed, eroded, speed
      double precision forceintx, forceinty
      double precision forcebedx, forcebedy
      double precision forcebedmax, forcebedequil, forcegrav
      double precision unitvx, unitvy, v_solid(2)
      double precision alphaxx, alphayy, alphaxy, alphaxz, alphayz
      double precision tanbed, terminal_vel, dragfoce(2)

      double precision fluxxp(3),fluxyp(3),tiny, uprev(3), ustore(3)
      double precision fluxxm(3), fluxym(3)
      double precision uvec(3), dUdx(3), dUdy(3)
      double precision h_inv, hphi_inv, curv(2), frict_tiny
      double precision intfrictang, bedfrictang, kactxy, dgdx(2)
      double precision dtdx, dtdy, dt, g(3), sgn_dudy, sgn_dvdx, tmp
      double precision dnorm, fluxsrc(3)
      double precision xslope,yslope,slope
      double precision t1, t2
      double precision erosion_rate,threshold,es,totalShear
      double precision eps, drag(4)

!     function calls
      double precision sgn

      integer i
      integer DO_EROSION, IF_STOPPED
      parameter(threshold=1.0D-02,erosion_rate=0.1)

c     initialize to zero
      forceintx=0.0
      forcebedx=0.0
      forceinty=0.0
      forcebedy=0.0
      unitvx=0.0
      unitvy=0.0
      eroded=0.0

      Ustore(1)=Uprev(1)
     1     -dtdx*(fluxxp(1)-fluxxm(1))
     2     -dtdy*(fluxyp(1)-fluxym(1))
     3     +dt*fluxsrc(1)
      Ustore(1)=dmax1(Ustore(1),0.0d0)

      Ustore(2)=Uprev(2)
     1     -dtdx*(fluxxp(2)-fluxxm(2))
     2     -dtdy*(fluxyp(2)-fluxym(2))
     3     +dt*fluxsrc(2)

      Ustore(3)=Uprev(3)
     1     -dtdx*(fluxxp(3)-fluxxm(3))
     2     -dtdy*(fluxyp(3)-fluxym(3))
     3     +dt*fluxsrc(3)


      ustore(1) = dmax1(ustore(1),0.)

      if(uvec(1).gt.tiny) then
c     Source terms ...
c     here speed is speed squared
         speed=v_solid(1)**2+v_solid(2)**2
         if(speed.gt.0.0) then
c     here speed is speed
            speed=dsqrt(speed)
            unitvx=v_solid(1)/speed
            unitvy=v_solid(2)/speed
         else
            unitvx=0.0
            unitvy=0.0
         endif
         tanbed=dtan(bedfrictang)
         h_inv = 1.d0/uvec(1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     x-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     the gravity force in the x direction
         forcegrav=g(1)*uvec(1)

c     the internal friction force
         tmp = h_inv*(dUdy(2)-v_solid(1)*dUdy(1))
         sgn_dudy = sgn(tmp, frict_tiny)
         forceintx=sgn_dudy*uvec(1)*kactxy*(g(3)*dUdy(1)
     $        +dgdx(2)*uvec(1))*dsin(intfrictang)

c     the bed friction force for fast moving flow 
         forcebedx=unitvx*
     $        dmax1(g(3)*Uvec(1)+v_solid(1)*Uvec(2)*curv(1),0.0d0)
     $        *tanbed

         if (abs(Ustore(2) + dt*forcegrav)
     $    .gt.abs(dt*(forcebedx+forceintx))) then
            Ustore(2) = Ustore(2) + dt*(forcegrav -forcebedx -forceintx)
         else
            Ustore(2)=0.d0
         endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     y-direction source terms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     the gravity force in the y direction
         forcegrav=g(2)*Uvec(1)

c     the internal friction force
         tmp = h_inv*(dudx(3)-v_solid(2)*dUdx(1))
         sgn_dvdx = sgn(tmp, frict_tiny)
         forceinty=sgn_dvdx*Uvec(1)*kactxy*(g(3)*
     $        dUdx(1)+dgdx(1)*Uvec(1))*dsin(intfrictang)

c     the bed friction force for fast moving flow 
         forcebedy=unitvy
     $        *dmax1(g(3)*Uvec(1)+v_solid(2)*Uvec(3)*curv(2),0.0d0)
     $        *tanbed

         if (abs(Ustore(3) + dt*forcegrav)
     $    .gt.abs(dt*(forcebedy+forceinty))) then
            Ustore(3) = Ustore(3) + dt*(forcegrav -forcebedy -forceinty)
         else
            Ustore(3)=0.d0
         endif

      endif

c     update the state variables
      Uvec(1) = Ustore(1)
      Uvec(2) = Ustore(2)
      Uvec(3) = Ustore(3)

      return
      end
