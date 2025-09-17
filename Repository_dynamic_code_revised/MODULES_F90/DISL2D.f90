      Module DISL2D

      REAL(8), PRIVATE :: mu1,mu2,nu1,nu2
      REAL(8), PRIVATE :: c1,c2,d2,gamm1,gamm2,kap1,kap2,pi
      REAL(8), PRIVATE :: Tp(4,4),Ts(4,4),Up(4,4),Us(4,4)

      Contains

      Subroutine DISL2D_DISPL(Ut,Ub,x0,z0,delta,Hd,x,z,Ux,Uz)
      IMPLICIT NONE
! input: Ut tensile component of the burger vector (opening)
!        Ub shear component of the burger vector (slip)
!        x0,z0 coordinates of the dislocation center (x horizontal, z vertical downward)
!        delta angle between the dislocation and z-axis (for vertical disl delta=0) positive counterclockwise
!        Hd lenght of the dislocation element
!        x,z coordinates of the point in which the displacement will be calculated
! output: Ux,Uz x and z components of the displacement induced by the dislocation
      REAL(8) :: Ut,Ub,x0,z0,delta,Hd,x,z,Ux,Uz

      Ux = uxUt(x0,z0,delta,Hd,x,z)*Ut + uxUb(x0,z0,delta,Hd,x,z)*Ub
      Uz = uzUt(x0,z0,delta,Hd,x,z)*Ut + uzUb(x0,z0,delta,Hd,x,z)*Ub
      
      End Subroutine DISL2D_DISPL



      Subroutine DISL2D_STRSS(Ut,Ub,x0,z0,delta,Hd,x,z,Sxx,Sxz,Szz)
      IMPLICIT NONE
! input: Ut tensile component of the burger vector (opening)
!        Ub shear component of the burger vector (slip)
!        x0,z0 coordinates of the dislocation center (x horizontal, z vertical downward)
!        delta angle between the dislocation and z-axis (for vertical disl delta=0) positive counterclockwise
!        Hd lenght of the dislocation element
!        x,z coordinates of the point in which the stress will be calculated
! output: Sxx,Sxz, Szz xx, xz and zz components of the stress field induced by the dislocation
      REAL(8) :: Ut,Ub,x0,z0,delta,Hd,x,z,Sxx,Sxz,Szz

      Sxx = sxxUt(x0,z0,delta,Hd,x,z)*Ut*1.D-3 + sxxUb(x0,z0,delta,Hd,x,z)*Ub*1.D-3
      Sxz = sxzUt(x0,z0,delta,Hd,x,z)*Ut*1.D-3 + sxzUb(x0,z0,delta,Hd,x,z)*Ub*1.D-3
      Szz = szzUt(x0,z0,delta,Hd,x,z)*Ut*1.D-3 + szzUb(x0,z0,delta,Hd,x,z)*Ub*1.D-3
      
      End Subroutine DISL2D_STRSS



      Subroutine SET_EPAR(mu_1,mu_2,nu_1,nu_2)
      IMPLICIT NONE

      REAL(8) :: mu_1,mu_2,nu_1,nu_2
      REAL(8) :: a1,a2,delt,gamm,cp,cm,d,e

      mu1=mu_1
      mu2=mu_2
      nu1=nu_1
      nu2=nu_2
      pi = 4.D0*atan(1.D0)

      If (mu2.eq.0.D0) Then ! homogeneous half space
      
        nu2 = 0.D0

        gamm1 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu1)
        gamm2 = 1.D0/(4.D0*pi)

        kap1  = 1.D0/(2.D0*pi*mu1)
        kap2  = 1.D0/(6.D0*pi*mu1)

        c1    = 0.D0
        c2    = 1.D0/(2.D0*pi) * mu1/(1.D0-nu1)
        d2    = 0.D0
      
      Else
      
        If ((nu1.ne.nu2).or.(mu1.ne.mu2)) Then ! bounded medium

          a1    = (3.D0-4.D0*nu1)/(4.D0*mu1*mu1)
          a2    = (3.D0-4.D0*nu2)/(4.D0*mu2*mu2)

          delt  = 1.D0/(2.D0*pi) * (mu2 /(1.D0-nu2) - mu1 /(1.D0-nu1))

          gamm1 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu1)
          gamm2 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu2)
          gamm  = gamm2 - gamm1

          kap1  = 1.D0/(2.D0*pi) * 1.D0/(mu1 + (3.D0 - 4.D0*nu1)*mu2)
          kap2  = 1.D0/(2.D0*pi) * 1.D0/(mu2 + (3.D0 - 4.D0*nu2)*mu1)

          cp    = (1.D0+(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2)
          cm    = (1.D0-(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2)

          d     = (1.D0-2.D0*nu2)/(2.D0*mu2)-(1.D0-2.D0*nu1)/(2.D0*mu1)
          e     = (1.D0-     nu2)/      mu2 +(1.D0-     nu1)/      mu1

          c1    = 1.D0/(e*e-d*d) *( delt*(a1+cp) + gamm*d)
          c2    = 1.D0/(e*e-d*d) *(-delt*(a2+cp) + gamm*d)
          d2    = 1.D0/(e*e-d*d) *( delt*cm      - gamm*e)

        Else ! homogeneous medium

          a1    = (3.D0-4.D0*nu1)/(4.D0*mu1*mu1)
          a2    = a1

          delt  = 0.D0

          gamm1 = 1.D0/(4.D0*pi) * 1.D0/(1.D0-nu1)
          gamm2 = gamm1
          gamm  = 0.D0

          kap1  = 1.D0/(2.D0*pi) * 1.D0/(mu1 + (3.D0 - 4.D0*nu1)*mu2)
          kap2  = kap1

          cp    = (1.D0+(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2)
          cm    = (1.D0-(3.D0-4.D0*nu1)*(3.D0-4.D0*nu2))/(8.D0*mu1*mu2)

          d     = 0.D0
          e     = 2.D0*(1.D0-nu1)/mu1

          c1    = 0.D0
          c2    = 0.D0
          d2    = 0.D0

        EndIf

      EndIf

      Tp(1,1)=   mu1*kap1 - mu2*kap2
      Tp(1,2)=   gamm1 - 2.D0*mu2*kap1
      Tp(1,3)=  -gamm1 + 2.D0*mu1*kap1
      Tp(1,4)=   2.D0*(gamm1 - 2.D0*mu2*kap1)

      Tp(2,1)=   mu1*kap1 - mu2*kap2
      Tp(2,2)=   gamm2 - 2.D0*mu1*kap2
      Tp(2,3)=  -gamm2 + 2.D0*mu1*kap1
      Tp(2,4)=   0.D0

      Tp(3,1)=   gamm1 - mu1*kap1 - mu2*kap2
      Tp(3,2)=   gamm1 - 2.D0*mu2*kap1
      Tp(3,3)=   gamm1 - 2.D0*mu1*kap1
      Tp(3,4)=   2.D0*(gamm1 - 2.D0*mu2*kap1)

      Tp(4,1)=   gamm2 - mu1*kap1 - mu2*kap2
      Tp(4,2)=  -gamm2 + 2.D0*mu1*kap2
      Tp(4,3)=   gamm2 - 2.D0*mu1*kap1
      Tp(4,4)=   0.D0

      Ts(1,1)=   mu2*kap2 - mu1*kap1
      Ts(1,2)=   gamm2 - 2.D0*mu1*kap2
      Ts(1,3)=  -gamm2 + 2.D0*mu2*kap2
      Ts(1,4)=   2.D0*(gamm2 - 2.D0*mu1*kap2)

      Ts(2,1)=   mu2*kap2 - mu1*kap1
      Ts(2,2)=   gamm1 - 2.D0*mu2*kap1
      Ts(2,3)=  -gamm1 + 2.D0*mu2*kap2
      Ts(2,4)=   0.D0

      Ts(3,1)=   gamm2 - mu2*kap2 - mu1*kap1
      Ts(3,2)=   gamm2 - 2.D0*mu1*kap2
      Ts(3,3)=   gamm2 - 2.D0*mu2*kap2
      Ts(3,4)=   2.D0*(gamm2 - 2.D0*mu1*kap2)

      Ts(4,1)=   gamm1 - mu2*kap2 - mu1*kap1
      Ts(4,2)=  -gamm1 + 2.D0*mu2*kap1
      Ts(4,3)=   gamm1 - 2.D0*mu2*kap2
      Ts(4,4)=   0.D0

      Up(1,1)= -(1.D0-     nu1)/(     mu1) * c2 +&
     &          (1.D0-2.D0*nu1)/(2.D0*mu1) * d2
      Up(1,2)=   1.D0/(2.D0*mu1) * (c2-d2)
      Up(1,3)=  (3.D0-4.D0*nu1)/(2.D0*mu1) * (c2-d2)
      Up(1,4)=  -1.D0/(     mu1) * (c2-d2)

      Up(2,1)=  (1.D0-     nu2)/(     mu2) * c1 +&
     &          (1.D0-2.D0*nu2)/(2.D0*mu2) * d2
      Up(2,2)=   1.D0/(2.D0*mu2) * (c1+d2)
      Up(2,3)=  -1.D0/(2.D0*mu2) * (c1-d2)
      Up(2,4)=   0.D0

      Up(3,1)=  (1.D0-     nu1)/(     mu1) * d2 -&
     &          (1.D0-2.D0*nu1)/(2.D0*mu1) * c2
      Up(3,2)=  -1.D0/(2.D0*mu1) * (c2-d2)
      Up(3,3)=  (3.D0-4.D0*nu1)/(2.D0*mu1) * (c2-d2)
      Up(3,4)=   1.D0/(     mu1) * (c2-d2)

      Up(4,1)= -(1.D0-     nu2)/(     mu2) * d2 -&
     &          (1.D0-2.D0*nu2)/(2.D0*mu2) * c1
      Up(4,2)=   1.D0/(2.D0*mu2) * (c1+d2)
      Up(4,3)=  -1.D0/(2.D0*mu2) * (c1-d2)
      Up(4,4)=   0.D0

      Us(1,1)= -(1.D0-     nu2)/(     mu2) * c1 -&
     &          (1.D0-2.D0*nu2)/(2.D0*mu2) * d2
      Us(1,2)=   1.D0/(2.D0*mu2) * (c1+d2)
      Us(1,3)=  (3.D0-4.D0*nu2)/(2.D0*mu2) * (c1+d2)
      Us(1,4)=  -1.D0/(     mu2) * (c1+d2)

      Us(2,1)=  (1.D0-     nu1)/(     mu1) * c2 -&
     &          (1.D0-2.D0*nu1)/(2.D0*mu1) * d2
      Us(2,2)=   1.D0/(2.D0*mu1) * (c2-d2)
      Us(2,3)=  -1.D0/(2.D0*mu1) * (c2+d2)
      Us(2,4)=   0.D0

      Us(3,1)= -(1.D0-     nu2)/(     mu2) * d2 -&
     &          (1.D0-2.D0*nu2)/(2.D0*mu2) * c1
      Us(3,2)=  -1.D0/(2.D0*mu2) * (c1+d2)
      Us(3,3)=  (3.D0-4.D0*nu2)/(2.D0*mu2) * (c1+d2)
      Us(3,4)=   1.D0/(     mu2) * (c1+d2)

      Us(4,1)=  (1.D0-     nu1)/(     mu1) * d2 -&
     &          (1.D0-2.D0*nu1)/(2.D0*mu1) * c2
      Us(4,2)=   1.D0/(2.D0*mu1) * (c2-d2)
      Us(4,3)=  -1.D0/(2.D0*mu1) * (c2+d2)
      Us(4,4)=   0.D0

      End Subroutine SET_EPAR

!     x comp of the displacement field due to the tensile component of the burger vector
      Function uxUt(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uxUt
      REAL(8) :: x0,z0,delta,Hd,x,z

      uxUt = uxUX(x0,z0,delta,Hd,x,z)*cos(delta) -&
     &       uxUZ(x0,z0,delta,Hd,x,z)*sin(delta)

      End Function uxUt

!     x comp of the displacement field due to the shear component of the burger vector
      Function uxUb(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uxUb
      REAL(8) :: x0,z0,delta,Hd,x,z

      uxUb = uxUX(x0,z0,delta,Hd,x,z)*sin(delta) +&
     &       uxUZ(x0,z0,delta,Hd,x,z)*cos(delta)

      End Function uxUb

!     z comp of the displacement field due to the tensile component of the burger vector
      Function uzUt(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uzUt
      REAL(8) :: x0,z0,delta,Hd,x,z

      uzUt = uzUX(x0,z0,delta,Hd,x,z)*cos(delta) -&
     &       uzUZ(x0,z0,delta,Hd,x,z)*sin(delta)

      End Function uzUt

!     z comp of the displacement field due to the shear component of the burger vector
      Function uzUb(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uzUb
      REAL(8) :: x0,z0,delta,Hd,x,z

      uzUb = uzUX(x0,z0,delta,Hd,x,z)*sin(delta) +&
     &       uzUZ(x0,z0,delta,Hd,x,z)*cos(delta)

      End Function uzUb

!     xx comp of the stress field due to the tensile component of the burger vector
      Function sxxUt(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxxUt
      REAL(8) :: x0,z0,delta,Hd,x,z

      sxxUt = sxxUX(x0,z0,delta,Hd,x,z)*cos(delta) -&
     &        sxxUZ(x0,z0,delta,Hd,x,z)*sin(delta)

      End Function sxxUt

!     xz comp of the stress field due to the tensile component of the burger vector
      Function sxzUt(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxzUt
      REAL(8) :: x0,z0,delta,Hd,x,z

      sxzUt = sxzUX(x0,z0,delta,Hd,x,z)*cos(delta) -&
     &        sxzUZ(x0,z0,delta,Hd,x,z)*sin(delta)

      End Function sxzUt

!     zz comp of the stress field due to the tensile component of the burger vector
      Function szzUt(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: szzUt
      REAL(8) :: x0,z0,delta,Hd,x,z

      szzUt = szzUX(x0,z0,delta,Hd,x,z)*cos(delta) -&
     &        szzUZ(x0,z0,delta,Hd,x,z)*sin(delta)

      End Function szzUt

!     xx comp of the stress field due to the shear component of the burger vector
      Function sxxUb(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxxUb
      REAL(8) :: x0,z0,delta,Hd,x,z

      sxxUb = sxxUX(x0,z0,delta,Hd,x,z)*sin(delta) +&
     &        sxxUZ(x0,z0,delta,Hd,x,z)*cos(delta)

      End Function sxxUb

!     xz comp of the stress field due to the shear component of the burger vector
      Function sxzUb(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxzUb
      REAL(8) :: x0,z0,delta,Hd,x,z

      sxzUb = sxzUX(x0,z0,delta,Hd,x,z)*sin(delta) +&
     &        sxzUZ(x0,z0,delta,Hd,x,z)*cos(delta)

      End Function sxzUb

!     zz comp of the stress field due to the shear component of the burger vector
      Function szzUb(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: szzUb
      REAL(8) :: x0,z0,delta,Hd,x,z

      szzUb = szzUX(x0,z0,delta,Hd,x,z)*sin(delta) +&
     &        szzUZ(x0,z0,delta,Hd,x,z)*cos(delta)

      End Function szzUb

!     x component of displacement due to x component of the burger vector
      Function uxUX(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uxUX
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      REAL(8) :: x1,z1,x2,z2,dGx
      INTEGER :: j

      uxUX = 0.D0

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_dGx (dGx,x0,z0,delta,Hd,x,z)

      uxUX = uxUX + dGx

      If (z1.ge.0.D0) Then

        Call fill_dY  (dY,x0,z0,delta,Hd,x,z)

        If (z.ge.0.D0) Then    ! formula for CASE1 within half-space 1
          Do j=1,4
          uxUX = uxUX + (-Up(1,j)*dY(1,j))
          EndDo
        Else                   ! formula for CASE1 within half-space 2
          Do j=1,4
          uxUX = uxUX + (-Up(2,j)*dY(2,j))
          EndDo
        EndIf
      Else
        If (z2.le.0.D0) Then

          Call fill_dY (dY,x0,z0,delta,Hd,x,z)

          If (z.le.0.D0) Then  ! formula for CASE2 within half-space 2
            Do j=1,4
            uxUX = uxUX + (-Us(1,j)*dY(1,j))
            EndDo
          Else                 ! formula for CASE2 within half-space 1
            Do j=1,4
            uxUX = uxUX + (-Us(2,j)*dY(2,j))
            EndDo
          EndIf
        Else                   ! formulas for mixed CASE (2-1)

          Call fill_Y(Y1,x1,z1,Hd,x,z)
          Call fill_Y(Y2,x2,z2,Hd,x,z)

          If (z.ge.0.D0) Then  ! formula for mixed CASE within half-space 1
            Do j=1,4
            uxUX = uxUX + (-Us(2,j)*Y1(2,j) + Up(1,j)*Y2(1,j))
            EndDo
          Else                 ! formula for mixed CASE within half-space 2
            Do j=1,4
            uxUX = uxUX + (-Us(1,j)*Y1(1,j) + Up(2,j)*Y2(2,j))
            EndDo
          EndIf
        EndIf
      EndIf

      End Function uxUX

!     x component of displacement due to z component of the burger vector
      Function uxUZ(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uxUZ
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      REAL(8) :: x1,z1,x2,z2,dHx
      INTEGER :: j

      uxUZ = 0.D0

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_dHx (dHx,x0,z0,delta,Hd,x,z)

      uxUZ = uxUZ - dHx

      If (z1.ge.0.D0) then

        Call fill_dY (dY,x0,z0,delta,Hd,x,z)

        If (z.ge.0.D0) then    ! formula for CASE1 within half-space 1
          Do j=1,4
          uxUZ = uxUZ - (-Tp(3,j)*dY(3,j))
          EndDo
        Else                   ! formula for CASE1 within half-space 2
          Do j=1,4
          uxUZ = uxUZ - (-Tp(4,j)*dY(4,j)) !+
          EndDo
        EndIf
      Else
        If (z2.le.0.D0) then

          Call fill_dY (dY,x0,z0,delta,Hd,x,z)

          If (z.le.0.D0) then  ! formula for CASE2 within half-space 2
            Do j=1,4
            uxUZ = uxUZ - (-Ts(3,j)*dY(3,j))
            EndDo
          Else                 ! formula for CASE2 within half-space 1
            Do j=1,4
            uxUZ = uxUZ - (-Ts(4,j)*dY(4,j)) !+
            EndDo
          EndIf
        Else                   ! formulas for mixed CASE (2-1)

          Call fill_Y(Y1,x1,z1,Hd,x,z)
          Call fill_Y(Y2,x2,z2,Hd,x,z)

          If (z.ge.0.D0) then  ! formula for mixed CASE within half-space 1
            Do j=1,4
            uxUZ = uxUZ - (-Ts(4,j)*Y1(4,j) + Tp(3,j)*Y2(3,j))
            EndDo
          Else                 ! formula for mixed CASE within half-space 2
            Do j=1,4
            uxUZ = uxUZ - (-Ts(3,j)*Y1(3,j) + Tp(4,j)*Y2(4,j))
            EndDo
          EndIf
        EndIf
      EndIf

      End Function uxUZ

!     z component of displacement due to x component of the burger vector
      Function uzUX(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uzUX
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      REAL(8) :: x1,z1,x2,z2,dGz
      INTEGER :: j

      uzUX = 0.D0

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_dGz (dGz,x0,z0,delta,Hd,x,z)

      uzUX = uzUX + dGz

      If (z1.ge.0.D0) then

        Call fill_dY (dY,x0,z0,delta,Hd,x,z)

        If (z.ge.0.D0) then    ! formula for CASE1 within half-space 1
          Do j=1,4
          uzUX = uzUX + (-Up(3,j)*dY(3,j))
          EndDo
        Else                   ! formula for CASE1 within half-space 2
          Do j=1,4
          uzUX = uzUX + (-Up(4,j)*dY(4,j))
          EndDo
        EndIf
      Else
        If (z2.le.0.D0) then

          Call fill_dY (dY,x0,z0,delta,Hd,x,z)

          If (z.le.0.D0) then  ! formula for CASE2 within half-space 2
            Do j=1,4
            uzUX = uzUX + (-Us(3,j)*dY(3,j))
            EndDo
          Else                 ! formula for CASE2 within half-space 1
            Do j=1,4
            uzUX = uzUX + (-Us(4,j)*dY(4,j))
            EndDo
          EndIf
        Else                   ! formulas for mixed CASE (2-1)

          Call fill_Y(Y1,x1,z1,Hd,x,z)
          Call fill_Y(Y2,x2,z2,Hd,x,z)

          If (z.ge.0.D0) then  ! formula for mixed CASE within half-space 1
            Do j=1,4
            uzUX = uzUX + (-Us(4,j)*Y1(4,j) + Up(3,j)*Y2(3,j))
            EndDo
          Else                 ! formula for mixed CASE within half-space 2
            Do j=1,4
            uzUX = uzUX + (-Us(3,j)*Y1(3,j) + Up(4,j)*Y2(4,j))
            EndDo
          EndIf
        EndIf
      EndIf

      End Function uzUX

!     z component of displacement due to z component of the burger vector
      Function uzUZ(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: uzUZ
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      REAL(8) :: x1,z1,x2,z2,dHz
      INTEGER :: j

      uzUZ = 0.D0

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_dHz (dHz,x0,z0,delta,Hd,x,z)

      uzUZ = uzUZ - dHz

      If (z1.ge.0.D0) then

        Call fill_dY (dY,x0,z0,delta,Hd,x,z)

        If (z.ge.0.D0) then    ! formula for CASE1 within half-space 1
          Do j=1,4
          uzUZ = uzUZ - (-Tp(1,j)*dY(1,j))
          EndDo
        Else                   ! formula for CASE1 within half-space 2
          Do j=1,4
          uzUZ = uzUZ - (-Tp(2,j)*dY(2,j))
          EndDo
        EndIf
      Else
        If (z2.le.0.D0) then

          Call fill_dY (dY,x0,z0,delta,Hd,x,z)

          If (z.le.0.D0) then  ! formula for CASE2 within half-space 2
            Do j=1,4
            uzUZ = uzUZ - (-Ts(1,j)*dY(1,j))
            EndDo
          Else                 ! formula for CASE2 within half-space 1
            Do j=1,4
            uzUZ = uzUZ - (-Ts(2,j)*dY(2,j))
            EndDo
          EndIf
        Else                   ! formulas for mixed CASE (2-1)

          Call fill_Y(Y1,x1,z1,Hd,x,z)
          Call fill_Y(Y2,x2,z2,Hd,x,z)

          If (z.ge.0.D0) then  ! formula for mixed CASE within half-space 1
            Do j=1,4
            uzUZ = uzUZ - (-Ts(2,j)*Y1(2,j) + Tp(1,j)*Y2(1,j))
            EndDo
          Else                 ! formula for mixed CASE within half-space 2
            Do j=1,4
            uzUZ = uzUZ - (-Ts(1,j)*Y1(1,j) + Tp(2,j)*Y2(2,j))
            EndDo
          EndIf
        EndIf
      EndIf

      End Function uzUZ

!     xx component of the stress field due to x component of the burger vector
      Function sxxUX(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxxUX
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      REAL(8) :: x1,z1,x2,z2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

      If (z1.ge.0.D0) Then

         Call fill_dI (dI,x0,z0,delta,Hd,x,z)

         If (z.lt.0.D0) Then    ! formula for CASE1 within half-space 2
         sxxUX =&
     &      -3.D0/(4.D0*(1.D0-nu2)) *                                   &
     &       2.D0*mu2/(3.D0*pi) *                                       &
     &   (  (z-z1)*(3.D0*(x-x1)**2+(z-z1)**2)/                          &
     &     ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &      (z-z2)*(3.D0*(x-x2)**2+(z-z2)**2)/                          &
     &     ((x-x2)**2+(z-z2)**2)**2         ) -                         &
     &   (2.D0*c1+d2)*dI(2,1) - (c1+d2)*dI(2,2) + (c1-d2)*dI(2,3)
         Else                   ! formula for CASE1 within half-space 1
         sxxUX =&
     &      -3.D0/(4.D0*(1.D0-nu1)) *                                   &
     &       2.D0*mu1/(3.D0*pi) *                                       &
     &   (  (z-z1)*(3.D0*(x-x1)**2+(z-z1)**2)/                          &
     &     ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &      (z-z2)*(3.D0*(x-x2)**2+(z-z2)**2)/                          &
     &     ((x-x2)**2+(z-z2)**2)**2         ) +                         &
     &   (2.D0*c2-d2)*dI(1,1) - (c2-d2)*dI(1,2) - 3.D0*(c2-d2)*dI(1,3) +&
     &    2.D0*(c2-d2)*dI(1,4)
         EndIf
         
      Else

         If (z2.le.0.D0) Then

            Call fill_dI (dI,x0,z0,delta,Hd,x,z)

            If (z.gt.0.D0) Then    ! formula for CASE2 within half-space 1
         sxxUX =&
     &      -3.D0/(4.D0*(1.D0-nu1)) *                                   &
     &       2.D0*mu1/(3.D0*pi) *                                       &
     &   (  (z-z1)*(3.D0*(x-x1)**2+(z-z1)**2)/                          &
     &     ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &      (z-z2)*(3.D0*(x-x2)**2+(z-z2)**2)/                          &
     &     ((x-x2)**2+(z-z2)**2)**2         ) -                         &
     &   (2.D0*c2-d2)*dI(2,1) - (c2-d2)*dI(2,2) + (c2+d2)*dI(2,3)
            Else                   ! formula for CASE2 within half-space 2
         sxxUX =&
     &      -3.D0/(4.D0*(1.D0-nu2)) *                                   &
     &       2.D0*mu2/(3.D0*pi) *                                       &
     &    ( (z-z1)*(3.D0*(x-x1)**2+(z-z1)**2)/                          &
     &     ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &      (z-z2)*(3.D0*(x-x2)**2+(z-z2)**2)/                          &
     &     ((x-x2)**2+(z-z2)**2)**2         ) +                         &
     &   (2.D0*c1+d2)*dI(1,1) - (c1+d2)*dI(1,2) - 3.D0*(c1+d2)*dI(1,3) +&
     &    2.D0*(c1+d2)*dI(1,4)
            EndIf
         Else                      ! formulas for mixed CASE (2-1)

         Call fill_I (Ih1,x1,z1,x,z)
         Call fill_I (Ih2,x2,z2,x,z)

              If (z.lt.0.D0) Then  ! formula for mixed CASE within half-space 2
         sxxUX =&
     &      -3.D0/(4.D0*(1.D0-nu2)) *                                   &
     &       2.D0*mu2/(3.D0*pi) *                                       &
     &    ( (z-z1)*(3.D0*(x-x1)**2+(z-z1)**2)/                          &
     &     ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &      (z-z2)*(3.D0*(x-x2)**2+(z-z2)**2)/                          &
     &     ((x-x2)**2+(z-z2)**2)**2         ) +                         &
     &   (2.D0*c1+d2)*(Ih1(1,1)+Ih2(2,1)) -                             &
     &   (c1+d2)*(Ih1(1,2)+3.D0*Ih1(1,3)-2.D0*Ih1(1,4)-Ih2(2,2)) -      &
     &   (c1-d2)*Ih2(2,3)
              Else                 ! formula for mixed CASE within half-space 1
         sxxUX =&
     &      -3.D0/(4.D0*(1.D0-nu1)) *                                   &
     &       2.D0*mu1/(3.D0*pi) *                                       &
     &    ( (z-z1)*(3.D0*(x-x1)**2+(z-z1)**2)/                          &
     &     ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &      (z-z2)*(3.D0*(x-x2)**2+(z-z2)**2)/                          &
     &     ((x-x2)**2+(z-z2)**2)**2         ) -                         &
     &   (2.D0*c2-d2)*(Ih2(1,1)+Ih1(2,1)) +                             &
     &   (c2-d2)*(Ih2(1,2)+3.D0*Ih2(1,3)-2.D0*Ih2(1,4)-Ih1(2,2)) +      &
     &   (c2+d2)*Ih1(2,3)
              EndIf
              
         EndIf

      EndIf

      End Function sxxUX

!     xz component of the stress field due to x component of the burger vector
      Function sxzUX (x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxzUX
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      REAL(8) :: x1,z1,x2,z2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

      If (z1.ge.0.D0) Then

           Call fill_dI (dI,x0,z0,delta,Hd,x,z)

           If (z.lt.0.D0) Then        ! formula for CASE1 within half-space 2
           sxzUX =-3.D0/(4.D0*(1.D0-nu2)) *                               &
     &             2.D0*mu2/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2   ) -                           &
     &          c1*dI(4,1) - (c1+d2)*dI(4,2) + (c1-d2)*dI(4,3)
           Else                       ! formula for CASE1 within half-space 1
           sxzUX =-3.D0/(4.D0*(1.D0-nu1)) *                               &
     &             2.D0*mu1/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2   ) -                           &
     &          c2*dI(3,1) + (c2-d2)*dI(3,2) +                            &
     &         (c2-d2)*dI(3,3) - 2.D0*(c2-d2)*dI(3,4)
           EndIf

      Else

           If (z2.le.0.D0) Then

              Call fill_dI (dI,x0,z0,delta,Hd,x,z)

              If (z.gt.0.D0) Then     ! formula for CASE2 within half-space 1
           sxzUX =-3.D0/(4.D0*(1.D0-nu1)) *                               &
     &             2.D0*mu1/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2   ) -                           &
     &          c2*dI(4,1) - (c2-d2)*dI(4,2) + (c2+d2)*dI(4,3)
              Else                    ! formula for CASE2 within half-space 2
           sxzUX =-3.D0/(4.D0*(1.D0-nu2)) *                               &
     &             2.D0*mu2/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2   ) -                           &
     &          c1*dI(3,1) + (c1+d2)*dI(3,2) +                            &
     &         (c1+d2)*dI(3,3) - 2.D0*(c1+d2)*dI(3,4)
              EndIf

           Else                       ! formulas for mixed CASE

              Call fill_I (Ih1,x1,z1,x,z)
              Call fill_I (Ih2,x2,z2,x,z)

              If (z.lt.0.D0) Then     ! formula for mixed CASE within half-space 2
           sxzUX =-3.D0/(4.D0*(1.D0-nu2)) *                               &
     &             2.D0*mu2/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2   ) -                           &
     &          c1*(Ih1(3,1)-Ih2(4,1)) +                                  &
     &         (c1+d2)*(Ih1(3,2)+Ih1(3,3)-2.D0*Ih1(3,4)+Ih2(4,2)) -       &
     &         (c1-d2)*Ih2(4,3)
              Else                    ! formula for mixed CASE within half-space 1
           sxzUX =-3.D0/(4.D0*(1.D0-nu1)) *                               &
     &             2.D0*mu1/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2   ) +                           &
     &          c2*(Ih2(3,1)-Ih1(4,1)) -                                  &
     &         (c2-d2)*(Ih2(3,2)+Ih2(3,3)-2.D0*Ih2(3,4)+Ih1(4,2)) +       &
     &         (c2+d2)*Ih1(4,3)
              EndIf

           EndIf

      EndIf

      End Function sxzUX

!     zz component of the stress field due to x component of the burger vector
      Function szzUX(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: szzUX
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      REAL(8) :: x1,z1,x2,z2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

      If (z1.ge.0.D0) Then

           Call fill_dI (dI,x0,z0,delta,Hd,x,z)

           If (z.lt.0.D0) Then        ! formula for CASE1 within half-space 2
            szzUX =-3.D0/(4.D0*(1.D0-nu2)) *                           &
     &              2.D0*mu2/(3.D0*pi) *                               &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                        &
     &            ((x-x1)**2+(z-z1)**2)**2 -                           &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                        &
     &            ((x-x2)**2+(z-z2)**2)**2   ) -                       &
     &           d2*dI(2,1) + (c1+d2)*dI(2,2) - (c1-d2)*dI(2,3)
           Else                       ! formula for CASE1 within half-space 1
            szzUX =-3.D0/(4.D0*(1.D0-nu1)) *                           &
     &              2.D0*mu1/(3.D0*pi) *                               &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                        &
     &            ((x-x1)**2+(z-z1)**2)**2 -                           &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                        &
     &            ((x-x2)**2+(z-z2)**2)**2   ) -                       &
     &           d2*dI(1,1) + (c2-d2)*dI(1,2) -                        &
     &          (c2-d2)*dI(1,3) - 2*(c2-d2)*dI(1,4)
           EndIf

      Else
           If (z2.le.0.D0) Then

              Call fill_dI (dI,x0,z0,delta,Hd,x,z)

              If (z.gt.0.D0) Then     ! formula for CASE2 within half-space 1
            szzUX =-3.D0/(4.D0*(1.D0-nu1)) *                           &
     &              2.D0*mu1/(3.D0*pi) *                               &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                        &
     &            ((x-x1)**2+(z-z1)**2)**2 -                           &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                        &
     &            ((x-x2)**2+(z-z2)**2)**2    ) +                      &
     &           d2*dI(2,1) + (c2-d2)*dI(2,2) - (c2+d2)*dI(2,3)
              Else                    ! formula for CASE2 within half-space 2
            szzUX =-3.D0/(4.D0*(1.D0-nu2)) *                           &
     &              2.D0*mu2/(3.D0*pi) *                               &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                        &
     &            ((x-x1)**2+(z-z1)**2)**2 -                           &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                        &
     &            ((x-x2)**2+(z-z2)**2)**2    ) +                      &
     &           d2*dI(1,1) + (c1+d2)*dI(1,2) -                        &
     &          (c1+d2)*dI(1,3) - 2*(c1+d2)*dI(1,4)
              EndIf

           Else                       ! formulas for mixed CASE

              Call fill_I (Ih1,x1,z1,x,z)
              Call fill_I (Ih2,x2,z2,x,z)

              If (z.lt.0.D0) Then     ! formula for mixed CASE within half-space 2
            szzUX =-3.D0/(4.D0*(1.D0-nu2)) *                           &
     &              2.D0*mu2/(3.D0*pi) *                               &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                        &
     &            ((x-x1)**2+(z-z1)**2)**2 -                           &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                        &
     &            ((x-x2)**2+(z-z2)**2)**2    ) +                      &
     &           d2*(Ih1(1,1)+Ih2(2,1)) +                              &
     &          (c1+d2)*(Ih1(1,2)-Ih1(1,3)-2.D0*Ih1(1,4)-Ih2(2,2)) +   &
     &          (c1-d2)*Ih2(2,3)
              Else                    ! formula for mixed CASE within half-space 1
            szzUX =-3.D0/(4.D0*(1.D0-nu1)) *                           &
     &              2.D0*mu1/(3.D0*pi) *                               &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                        &
     &            ((x-x1)**2+(z-z1)**2)**2 -                           &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                        &
     &            ((x-x2)**2+(z-z2)**2)**2    ) +                      &
     &           d2*(Ih2(1,1)+Ih1(2,1)) -                              &
     &          (c2-d2)*(Ih2(1,2)-Ih2(1,3)-2.D0*Ih2(1,4)-Ih1(2,2)) -   &
     &          (c2+d2)*Ih1(2,3)
              EndIf

           EndIf

      EndIf

      End Function szzUX

!     xx component of the stress field due to z component of the burger vector
      Function sxxUZ(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxxUZ
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      REAL(8) :: x1,z1,x2,z2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

      If (z1.ge.0.D0) Then

           Call fill_dI (dI,x0,z0,delta,Hd,x,z)

           If (z.lt.0.D0) Then        ! formula for CASE1 within half-space 2
           sxxUZ =-3.D0/(4.D0*(1.D0-nu2)) *                               &
     &             2.D0*mu2/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2    ) -                          &
     &         (c1+2.D0*d2)*dI(4,1) - (c1+d2)*dI(4,2) + (c1-d2)*dI(4,3)
           Else                       ! formula for CASE1 within half-space 1
           sxxUZ =-3.D0/(4.D0*(1.D0-nu1)) *                               &
     &             2.D0*mu1/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2    ) -                          &
     &         (c2-2.D0*d2)*dI(3,1) + (c2-d2)*dI(3,2) -                   &
     &          3.D0*(c2-d2)*dI(3,3) + 2.D0*(c2-d2)*dI(3,4)
           EndIf

      Else

           If (z2.le.0.D0) Then

              Call fill_dI (dI,x0,z0,delta,Hd,x,z)

              If (z.gt.0.D0) Then     ! formula for CASE2 within half-space 1
           sxxUZ =-3.D0/(4.D0*(1.D0-nu1)) *                               &
     &             2.D0*mu1/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2    ) -                          &
     &         (c2-2.D0*d2)*dI(4,1) - (c2-d2)*dI(4,2) + (c2+d2)*dI(4,3)
              Else                    ! formula for CASE2 within half-space 2
           sxxUZ =-3.D0/(4.D0*(1.D0-nu2)) *                               &
     &             2.D0*mu2/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2    ) -                          &
     &         (c1+2.D0*d2)*dI(3,1) + (c1+d2)*dI(3,2) -                   &
     &         3.D0*(c1+d2)*dI(3,3) + 2.D0*(c1+d2)*dI(3,4)
              EndIf

           Else                       ! formulas for mixed CASE

              Call fill_I (Ih1,x1,z1,x,z)
              Call fill_I (Ih2,x2,z2,x,z)

              If (z.lt.0.D0) Then     ! formula for mixed CASE within half-space 2
           sxxUZ =-3.D0/(4.D0*(1.D0-nu2)) *                               &
     &             2.D0*mu2/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2    ) -                          &
     &         (c1+2.D0*d2)*(Ih1(3,1)-Ih2(4,1)) +                         &
     &         (c1+d2)*(Ih1(3,2)-3.D0*Ih1(3,3)+2.D0*Ih1(3,4)+Ih2(4,2)) -  &
     &         (c1-d2)*Ih2(4,3)
              Else                    ! formula for mixed CASE within half-space 1
           sxxUZ =-3.D0/(4.D0*(1.D0-nu1)) *                               &
     &             2.D0*mu1/(3.D0*pi) *                                   &
     &         ( (x-x1)*((z-z1)**2-(x-x1)**2)/                            &
     &           ((x-x1)**2+(z-z1)**2)**2 -                               &
     &           (x-x2)*((z-z2)**2-(x-x2)**2)/                            &
     &           ((x-x2)**2+(z-z2)**2)**2    ) -                          &
     &         (c2-2.D0*d2)*(Ih1(4,1)-Ih2(3,1)) -                         &
     &         (c2-d2)*(Ih1(4,2)+Ih2(3,2)-3.D0*Ih2(3,3)+2.D0*Ih2(3,4)) +  &
     &         (c2+d2)*Ih1(4,3)
              Endif

           EndIf

      EndIf

      End Function sxxUZ

!     xz component of the stress field due to z component of the burger vector
      Function sxzUZ(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: sxzUZ
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      REAL(8) :: x1,z1,x2,z2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

      If (z1.ge.0.D0) Then

           Call fill_dI (dI,x0,z0,delta,Hd,x,z)

           If (z.lt.0.D0) Then        ! formula for CASE1 within half-space 2
            sxzUZ =-3.D0/(4.D0*(1.D0-nu2)) *                              &
     &              2.D0*mu2/(3.D0*pi) *                                  &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                           &
     &            ((x-x1)**2+(z-z1)**2)**2 -                              &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                           &
     &            ((x-x2)**2+(z-z2)**2)**2    ) +                         &
     &          d2*dI(2,1) + (c1+d2)*dI(2,2) - (c1-d2)*dI(2,3)
           Else                       ! formula for CASE1 within half-space 1
            sxzUZ =-3.D0/(4.D0*(1.D0-nu1)) *                              &
     &              2.D0*mu1/(3.D0*pi) *                                  &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                           &
     &            ((x-x1)**2+(z-z1)**2)**2 -                              &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                           &
     &            ((x-x2)**2+(z-z2)**2)**2    ) +                         &
     &          d2*dI(1,1) + (c2-d2)*dI(1,2) -                            &
     &          (c2-d2)*dI(1,3) + 2*(c2-d2)*dI(1,4)
           EndIf

      Else
           If (z2.le.0.D0) Then

              Call fill_dI (dI,x0,z0,delta,Hd,x,z)

              If (z.gt.0.D0) Then     ! formula for CASE2 within half-space 1
            sxzUZ =-3.D0/(4.D0*(1.D0-nu1)) *                              &
     &              2.D0*mu1/(3.D0*pi) *                                  &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                           &
     &            ((x-x1)**2+(z-z1)**2)**2 -                              &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                           &
     &            ((x-x2)**2+(z-z2)**2)**2    ) -                         &
     &          d2*dI(2,1) + (c2-d2)*dI(2,2) - (c2+d2)*dI(2,3)
              Else                    ! formula for CASE2 within half-space 2
            sxzUZ =-3.D0/(4.D0*(1.D0-nu2)) *                              &
     &              2.D0*mu2/(3.D0*pi) *                                  &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                           &
     &            ((x-x1)**2+(z-z1)**2)**2 -                              &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                           &
     &            ((x-x2)**2+(z-z2)**2)**2    ) -                         &
     &          d2*dI(1,1) + (c1+d2)*dI(1,2) -                            &
     &          (c1+d2)*dI(1,3) + 2*(c1+d2)*dI(1,4)
              EndIf

           Else                       ! formulas for mixed CASE

              Call fill_I (Ih1,x1,z1,x,z)
              Call fill_I (Ih2,x2,z2,x,z)

              If (z.lt.0.D0) Then     ! formula for mixed CASE within half-space 2
            sxzUZ =-3.D0/(4.D0*(1.D0-nu2)) *                              &
     &              2.D0*mu2/(3.D0*pi) *                                  &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                           &
     &            ((x-x1)**2+(z-z1)**2)**2 -                              &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                           &
     &            ((x-x2)**2+(z-z2)**2)**2   ) -                          &
     &          d2*(Ih1(1,1)+Ih2(2,1)) +                                  &
     &          (c1+d2)*(Ih1(1,2)-Ih1(1,3)+2.D0*Ih1(1,4)-Ih2(2,2)) +      &
     &          (c1-d2)*Ih2(2,3)
              Else                    ! formula for mixed CASE within half-space 1
            sxzUZ =-3.D0/(4.D0*(1.D0-nu1)) *                              &
     &              2.D0*mu1/(3.D0*pi) *                                  &
     &          ( (z-z1)*((z-z1)**2-(x-x1)**2)/                           &
     &            ((x-x1)**2+(z-z1)**2)**2 -                              &
     &            (z-z2)*((z-z2)**2-(x-x2)**2)/                           &
     &            ((x-x2)**2+(z-z2)**2)**2    ) -                         &
     &          d2*(Ih1(2,1)+Ih2(1,1)) +                                  &
     &          (c2-d2)*(Ih1(2,2)-Ih2(1,2)+Ih2(1,3)-2.D0*Ih2(1,4)) -      &
     &          (c2+d2)*Ih1(2,3)
              EndIf

           Endif

      EndIf

      End Function sxzUZ

!     zz component of the stress field due to z component of the burger vector
      Function szzUZ(x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: szzUZ
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),Ih1(4,4),Ih2(4,4)
      REAL(8) :: x1,z1,x2,z2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

      If (z1.ge.0.D0) Then

         Call fill_dI (dI,x0,z0,delta,Hd,x,z)

         If (z.lt.0.D0) Then        ! formula for CASE1 within half-space 2
         szzUZ = 3.D0/(4.D0*(1.D0-nu2)) *                                 &
     &         (2.D0*mu2)/(3.D0*pi)*                                      &
     &      ( (x-x1)*((x-x1)**2+3.D0*(z-z1)**2)/                          &
     &        ((x-x1)**2+(z-z1)**2)**2 -                                  &
     &        (x-x2)*((x-x2)**2+3.D0*(z-z2)**2)/                          &
     &        ((x-x2)**2+(z-z2)**2)**2         ) -                        &
     &      c1*dI(4,1) + (c1+d2)*dI(4,2) - (c1-d2)*dI(4,3)                
         Else                       ! formula for CASE1 within half-space 1
         szzUZ = 3.D0/(4.D0*(1.D0-nu1)) *                                 &
     &         (2.D0*mu1)/(3.D0*pi)*                                      &
     &      ( (x-x1)*((x-x1)**2+3.D0*(z-z1)**2)/                          &
     &        ((x-x1)**2+(z-z1)**2)**2 -                                  &
     &        (x-x2)*((x-x2)**2+3.D0*(z-z2)**2)/                          &
     &        ((x-x2)**2+(z-z2)**2)**2         ) -                        &
     &      c2*dI(3,1) - (c2-d2)*dI(3,2) - (c2-d2)*dI(3,3) -              &
     &      2.D0*(c2-d2)*dI(3,4)
         EndIf
      Else

         If (z2.le.0.D0) Then

            Call fill_dI (dI,x0,z0,delta,Hd,x,z)

            If (z.gt.0.D0) Then     ! formula for CASE2 within half-space 1
         szzUZ = 3.D0/(4.D0*(1.D0-nu1)) *                                 &
     &         (2.D0*mu1)/(3.D0*pi)*                                      &
     &      ( (x-x1)*((x-x1)**2+3.D0*(z-z1)**2)/                          &
     &        ((x-x1)**2+(z-z1)**2)**2 -                                  &
     &        (x-x2)*((x-x2)**2+3.D0*(z-z2)**2)/                          &
     &        ((x-x2)**2+(z-z2)**2)**2         ) -                        &
     &      c2*dI(4,1) + (c2-d2)*dI(4,2) - (c2+d2)*dI(4,3)
            Else                    ! formula for CASE2 within half-space 2
         szzUZ = 3.D0/(4.D0*(1.D0-nu2)) *                                 &
     &        (2.D0*mu2)/(3.D0*pi)*                                       &
     &     ( (x-x1)*((x-x1)**2+3.D0*(z-z1)**2)/                           &
     &       ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &       (x-x2)*((x-x2)**2+3.D0*(z-z2)**2)/                           &
     &       ((x-x2)**2+(z-z2)**2)**2          ) -                        &
     &      c1*dI(3,1) - (c1+d2)*dI(3,2) - (c1+d2)*dI(3,3) -              &
     &      2.D0*(c1+d2)*dI(3,4)
            EndIf
         Else                       ! formulas for mixed CASE

         Call fill_I (Ih1,x1,z1,x,z)
         Call fill_I (Ih2,x2,z2,x,z)

            If (z.lt.0.D0) Then     ! formula for mixed CASE within half-space 2
         szzUZ = 3.D0/(4.D0*(1.D0-nu2)) *                                 &
     &        (2.D0*mu2)/(3.D0*pi)*                                       &
     &     ( (x-x1)*((x-x1)**2+3.D0*(z-z1)**2)/                           &
     &       ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &       (x-x2)*((x-x2)**2+3.D0*(z-z2)**2)/                           &
     &       ((x-x2)**2+(z-z2)**2)**2         ) -                         &
     &      c1*(Ih1(3,1)-Ih2(4,1)) -                                      &
     &     (c1+d2)*(Ih1(3,2)+Ih1(3,3)+2.D0*Ih1(3,4)+Ih2(4,2)) +           &
     &     (c1-d2)*Ih2(4,3)
            Else                    ! formula for mixed CASE within half-space 1
         szzUZ = 3.D0/(4.D0*(1.D0-nu1)) *                                 &
     &        (2.D0*mu1)/(3.D0*pi)*                                       &
     &     ( (x-x1)*((x-x1)**2+3.D0*(z-z1)**2)/                           &
     &       ((x-x1)**2+(z-z1)**2)**2 -                                   &
     &       (x-x2)*((x-x2)**2+3.D0*(z-z2)**2)/                           &
     &       ((x-x2)**2+(z-z2)**2)**2        ) -                          &
     &      c2*(Ih1(4,1)-Ih2(3,1)) +                                      &
     &     (c2-d2)*(Ih1(4,2)+Ih2(3,2)+Ih2(3,3)+2.D0*Ih2(3,4)) -           &
     &     (c2+d2)*Ih1(4,3)
            EndIf
            
         EndIf

      EndIf

      End Function szzUZ


      Subroutine fill_dI(dI,x0,z0,delta,Hd,x,z)
      IMPLICIT NONE
      
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dI(4,4),I1(4,4),I2(4,4)
      REAL(8) :: x1,z1,x2,z2
      INTEGER :: n,m

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

           Call fill_I (I1,x1,z1,x,z)
           Call fill_I (I2,x2,z2,x,z)

                Do n=1,4
                   Do m=1,4
                      dI(n,m) = 0.D0
                      dI(n,m) = I1(n,m) - I2(n,m)
                   EndDo
                EndDo

      End Subroutine fill_dI


      Subroutine fill_dY(dY,x0,z0,delta,Hd,x,z)
      IMPLICIT NONE
      
      REAL(8) :: x0,z0,delta,Hd,x,z

      REAL(8) :: dY(4,4),Y1(4,4),Y2(4,4)
      REAL(8) :: x1,z1,x2,z2
      INTEGER :: n,m

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

       If (abs(z1).lt.1.D-12) z1 = 0.D0
       If (abs(z2).lt.1.D-12) z2 = 0.D0

           Call fill_Y (Y1,x1,z1,Hd,x,z)
           Call fill_Y (Y2,x2,z2,Hd,x,z)

                Do n=1,4
                   Do m=1,4
                      dY(n,m) = 0.D0
                      dY(n,m) = Y1(n,m) - Y2(n,m)
                   EndDo
                EndDo

      End Subroutine fill_dY


      Subroutine fill_I(I,xd,zd,x,z)
      IMPLICIT NONE

      REAL(8) :: I(4,4)
      REAL(8) :: xd,zd,x,z

                 I(1,1)=(zd+z)/((x-xd)**2 + (z+zd)**2)

                 I(1,2)=z*((z+zd)**2 - (x-xd)**2)/         &
     &                    ((x-xd)**2 + (z+zd)**2)**2

                 I(1,3)=zd*((z+zd)**2 - (x-xd)**2)/        &
     &                     ((x-xd)**2 + (z+zd)**2)**2

                 I(1,4)=2.D0*z*zd*(zd+z)*                  &
     &                  ((z+zd)**2 - 3.D0*(x-xd)**2)/      &
     &                  ((x-xd)**2 + (z+zd)**2)**3

                 I(2,1)=(zd-z)/((x-xd)**2 + (z-zd)**2)

                 I(2,2)=z*((z-zd)**2 - (x-xd)**2)/         &
     &                    ((x-xd)**2 + (z-zd)**2)**2

                 I(2,3)=zd*((z-zd)**2 - (x-xd)**2)/        &
     &                     ((x-xd)**2 + (z-zd)**2)**2

                 I(2,4)=2.D0*z*zd*(zd-z)*                  &
     &                  ((z-zd)**2 - 3.D0*(x-xd)**2)/      &
     &                  ((x-xd)**2 + (z-zd)**2)**3

                 I(3,1)=(x-xd)/((x-xd)**2 + (z+zd)**2)

                 I(3,2)=2.D0*z*(x-xd)*(zd+z)/              &
     &                      ((x-xd)**2 + (z+zd)**2)**2

                 I(3,3)=2.D0*zd*(x-xd)*(zd+z)/             &
     &                       ((x-xd)**2 + (z+zd)**2)**2

                 I(3,4)=2.D0*z*zd*(x-xd)*                  &
     &                  (3.D0*(zd+z)**2 - (x-xd)**2)/      &
     &                  ((x-xd)**2 + (z+zd)**2)**3

                 I(4,1)=(x-xd)/((x-xd)**2 + (z-zd)**2)

                 I(4,2)=2.D0*z*(x-xd)*(zd-z)/              &
     &                      ((x-xd)**2 + (z-zd)**2)**2

                 I(4,3)=2.D0*zd*(x-xd)*(zd-z)/             &
     &                       ((x-xd)**2 + (z-zd)**2)**2

                 I(4,4)=2.D0*z*zd*(x-xd)*                  &
     &                  (3.D0*(zd-z)**2 - (x-xd)**2)/      &
     &                  ((x-xd)**2 + (z-zd)**2)**3

      End Subroutine fill_I


      Subroutine fill_Y(Y,xd,zd,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: Y(4,4)
      REAL(8) :: xd,zd,x,z,Hd

      REAL(8) :: arctanP,arctanM,segnoP,segnoM

      If (abs(x-xd).lt.1.D-12) Then ! IMPORTANT: this will happen every time fill_Y
                                    ! is called for computing the displacement at the
                                    ! CENTER of a VERTICAL elementary disl. (that is
                                    ! needed for computing the asimmetry in the opening
                                    ! and shear displacement)
        If (zd+z.gt.0.D0) Then
          arctanP =  pi*0.5D0
        Else
          arctanP = -pi*0.5D0
        EndIf

        If (zd-z.gt.0.D0) Then
          arctanM =  pi*0.5D0
        Else
          arctanM = -pi*0.5D0
        EndIf

        If (abs(zd).gt.1.D-12) Then

          If (zd.gt.1.D-12) Then
            segnoP = +1.D0
            segnoM = +1.D0
          Else
            segnoP = -1.D0
            segnoM = -1.D0
          EndIf

        Else
          segnoP = +sign(1.D0,z)
          segnoM = -sign(1.D0,z)
        EndIf

      Else ! (when abs(x-xd) > 1.D-12)

        arctanP = atan((zd+z)/(x-xd))
        arctanM = atan((zd-z)/(x-xd))

        If (abs(zd).lt.1.D-12) Then
          segnoP = +sign(1.D0,(x-xd)*z)!+(x-xd)/abs(x-xd)*sign(1.D0,z)
          segnoM = -sign(1.D0,(x-xd)*z)!-(x-xd)/abs(x-xd)*sign(1.D0,z)
        Else
          segnoP =  sign(1.D0,(x-xd)*zd)!(x-xd)*zd/(abs((x-xd)*zd))
          segnoM =  sign(1.D0,(x-xd)*zd)!(x-xd)*zd/(abs((x-xd)*zd))
        EndIf

      EndIf

!       NUMERICAL LIMIT INSTEAD OF ANALYTICAL (SEE ABOVE)
!       zdp = zd
!       zdm = zd
! 
!       If (abs(x-xd).lt.1.D-14) Then !(disl. verticale)
!         x = x  + 1.D-12
!       EndIf
!       If (abs(z-zd).lt.1.D-14) Then !(disl. orizzontale)
!         If (abs(zd).lt.1.D-14) Then
!           zdp = zd + 1.D-12
!           zdm = zd + 1.D-12
!         EndIf
!       Else
!         If (abs(zd).lt.1.D-14) Then
!           zdp = zd + 1.D-12*z/abs(z)
!           zdm = zd - 1.D-12*z/abs(z)
!         EndIf
!       EndIf
! 
!                  Y(1,1)= -(-pi*0.5D0*(x-xd)*zdp/abs((x-xd)*zdp)
!      &                   + atan((z+zd)/(x-xd)))

                 Y(1,1)= pi*0.5D0*segnoP - arctanP

                 Y(1,2)=  z*(x-xd)/                        &
     &                      ((x-xd)**2 + (z+zd)**2)

                 Y(1,3)= zd*(x-xd)/                        &
     &                      ((x-xd)**2 + (z+zd)**2)

                 Y(1,4)=2.D0*z*zd*(z+zd)*(x-xd)/           &
     &                     (((x-xd)**2 + (z+zd)**2)**2)

!                  Y(2,1)= -(-pi*0.5D0*(x-xd)*zdm/abs((x-xd)*zdm)
!      &                   + atan((zd-z)/(x-xd)))

                 Y(2,1)= pi*0.5D0*segnoM - arctanM

                 Y(2,2)= z*(x-xd)/                         &
     &                   ((x-xd)**2 + (z-zd)**2)

                 Y(2,3)= zd*(x-xd)/                        &
     &                   ((x-xd)**2 + (z-zd)**2)           

                 Y(2,4)= 2.D0*z*zd*(zd-z)*(x-xd)/          &
     &                   (((x-xd)**2 + (z-zd)**2)**2)

                 Y(3,1)= 0.5D0*log( ((x-xd)**2 + (z+zd)**2)/    &
     &                              (Hd**2) )

                 Y(3,2)= -z*(z+zd)/                        &
     &                   ((x-xd)**2 + (z+zd)**2)

                 Y(3,3)= -zd*(z+zd)/                       &
     &                   ((x-xd)**2 + (z+zd)**2)

                 Y(3,4)= z*zd*((x-xd)**2 - (z+zd)**2)/     &
     &                          (((x-xd)**2 + (z+zd)**2)**2)

                 Y(4,1)= 0.5D0*log( ((x-xd)**2 + (z-zd)**2)/    &
     &                              (Hd**2) )

                 Y(4,2)= -z*(zd-z)/                        &
     &                   ((x-xd)**2 + (z-zd)**2)

                 Y(4,3)= -zd*(zd-z)/                       &
     &                   ((x-xd)**2 + (z-zd)**2)

                 Y(4,4)= z*zd*((x-xd)**2 - (z-zd)**2)/     &
     &                          (((x-xd)**2 + (z-zd)**2)**2)

      End Subroutine fill_Y


      Subroutine fill_dGx(dGx,x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: dGx,x0,z0,x,z,Hd,delta

      REAL(8) :: z1,z2,x1,x2
      REAL(8) :: Gx1,Gx2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_Gx (Gx1,x1,z1,delta,Hd,x,z)
      Call fill_Gx (Gx2,x2,z2,delta,Hd,x,z)

      dGx = Gx1 - Gx2

      End  Subroutine fill_dGx


      Subroutine fill_dGz(dGz,x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: dGz,x0,z0,x,z,Hd,delta

      REAL(8) :: z1,z2,x1,x2
      REAL(8) :: Gz1,Gz2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_Gz (Gz1,x1,z1,delta,Hd,x,z)
      Call fill_Gz (Gz2,x2,z2,delta,Hd,x,z)

      dGz = Gz1 - Gz2

      End Subroutine fill_dGz


      Subroutine fill_dHx(dHx,x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: dHx,x0,z0,x,z,Hd,delta

      REAL(8) :: z1,z2,x1,x2
      REAL(8) :: Hx1,Hx2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_Hx (Hx1,x1,z1,delta,Hd,x,z)
      Call fill_Hx (Hx2,x2,z2,delta,Hd,x,z)

      dHx = Hx1 - Hx2

      End Subroutine fill_dHx


      Subroutine fill_dHz(dHz,x0,z0,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: dHz,x0,z0,x,z,Hd,delta

      REAL(8) :: z1,z2,x1,x2
      REAL(8) :: Hz1,Hz2

      z1 = z0 - (Hd*0.5D0)*cos(delta)
      z2 = z0 + (Hd*0.5D0)*cos(delta)

      x1 = x0 - (Hd*0.5D0)*sin(delta)
      x2 = x0 + (Hd*0.5D0)*sin(delta)

      If (abs(z1).lt.1.D-12) z1 = 0.D0
      If (abs(z2).lt.1.D-12) z2 = 0.D0

      Call fill_Hz (Hz1,x1,z1,delta,Hd,x,z)
      Call fill_Hz (Hz2,x2,z2,delta,Hd,x,z)

      dHz = Hz1 - Hz2

      End Subroutine fill_dHz


      Subroutine fill_Gx(Gx,xd,zd,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: Gx,xd,zd,x,z,Hd,delta

      REAL(8) :: nu,Arctan

      If (z.ge.0.D0) Then
      nu = nu1
      Else
      nu = nu2
      EndIf

      If ((x-xd)*cos(delta)-(z-zd)*sin(delta).gt.1.D-12) Then
      Arctan =  pi*0.5D0 +                                           &
     &          atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/         &
     &                ((x-xd)*cos(delta)-(z-zd)*sin(delta)) )
      Else
      If ((x-xd)*cos(delta)-(z-zd)*sin(delta).lt.-1.D-12) Then
      Arctan = -pi*0.5D0 +                                           &
     &          atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/         &
     &                ((x-xd)*cos(delta)-(z-zd)*sin(delta)) )
      Else
      If ((x-xd)*sin(delta)+(z-zd)*cos(delta).gt.1.D-12) Then
      Arctan = pi
      Else
      If ((x-xd)*sin(delta)+(z-zd)*cos(delta).lt.-1.D-12) Then
      Arctan = 0.D0
      Else
      Arctan = 3.D0*pi/4.D0
      EndIf
      EndIf
      EndIf
      EndIf

      Gx = 1.D0/(2.D0*pi)*                                           &
     &     (Arctan + 1.D0/(2.D0*(1.D0-nu))*                          &
     &                (((x-xd)*(z-zd))/((x-xd)**2+(z-zd)**2)) )

      End Subroutine fill_Gx


      Subroutine fill_Gz(Gz,xd,zd,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: Gz,xd,zd,x,z,Hd,delta

      REAL(8) :: nu

      If (z.ge.0.D0) Then
      nu = nu1
      Else
      nu = nu2
      EndIf

      Gz = -1.D0/(4.D0*pi*(1.D0-nu))*                                &
     &      ( (1.D0-2.D0*nu)*0.5D0*log( ((x-xd)**2+(z-zd)**2)/       &
     &                                   (Hd**2) ) -                 &
     &        (z-zd)**2/((x-xd)**2+(z-zd)**2) )

      End Subroutine fill_Gz


      Subroutine fill_Hx(Hx,xd,zd,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: Hx,xd,zd,x,z,Hd,delta

      REAL(8) :: nu

      If (z.ge.0.D0) Then
      nu = nu1
      Else
      nu = nu2
      EndIf
                                                                     
      Hx = -1.D0/(4.D0*pi*(1.D0-nu))*                                &
     &      ( (1.D0-2.D0*nu)*0.5D0*log( ((x-xd)**2+(z-zd)**2)/       &
     &                                   (Hd**2) ) +                 &
     &        (z-zd)**2/((x-xd)**2+(z-zd)**2) )

      End Subroutine fill_Hx


      Subroutine fill_Hz(Hz,xd,zd,delta,Hd,x,z)
      IMPLICIT NONE

      REAL(8) :: Hz,xd,zd,x,z,Hd,delta

      REAL(8) :: nu,Arctan

      If (z.ge.0.D0) Then
      nu = nu1
      Else
      nu = nu2
      EndIf

      If ((x-xd)*cos(delta)-(z-zd)*sin(delta).gt.1.D-12) Then
      Arctan =  pi*0.5D0 +                                            &
     &          atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/          &
     &                ((x-xd)*cos(delta)-(z-zd)*sin(delta)) )
      Else
      If ((x-xd)*cos(delta)-(z-zd)*sin(delta).lt.-1.D-12) Then
      Arctan = -pi*0.5D0 +                                            &
     &          atan( ((x-xd)*sin(delta)+(z-zd)*cos(delta))/          &
     &                ((x-xd)*cos(delta)-(z-zd)*sin(delta)) )
      Else
      If ((x-xd)*sin(delta)+(z-zd)*cos(delta).gt.1.D-12) Then
      Arctan = pi
      Else
      If ((x-xd)*sin(delta)+(z-zd)*cos(delta).lt.-1.D-12) Then
      Arctan = 0.D0
      Else
      Arctan = 3.D0*pi/4.D0
      EndIf
      EndIf
      EndIf
      EndIf

      Hz = -1.D0/(2.D0*pi)*                                           &
     &     (Arctan - 1.D0/(2.D0*(1.D0-nu))*                           &
     &                 (((x-xd)*(z-zd))/((x-xd)**2+(z-zd)**2)) )
     
      End Subroutine fill_Hz

      End Module DISL2D
