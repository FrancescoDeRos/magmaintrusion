      Module CRACK2D
      Use DISL2D, only : DISL2D_DISPL, DISL2D_STRSS
      
      contains

      Subroutine CRACK2D_DISPL(Ut,Ub,x0z0,delta,Hd,n_dis,x,z,Ux,Uz)

      IMPLICIT NONE
! input: Ut tensile component of the burger vector of each element of the crack (opening)
!        Ub shear component of the burger vector of each element of the crack (slip)
!        x0,z0 coordinates of each dislocation center (x horizontal, z vertical downward)
!        delta angle between the dislocation and z-axis of each element of the crack (for vertical disl delta=0) positive counterclockwise
!        Hd lenght of each dislocation element
!        x,z coordinates of the point in which the displacement will be calculated
! output: Ux,Uz x and z components of the displacement induced by the boundary element crack
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:),x,z,Ux,Uz
      INTEGER :: n_dis
      REAL(8) :: x0,z0,Ux_k,Uz_k
      INTEGER :: k

      Ux=0.D0
      Uz=0.D0
      
      Do k=1,n_dis

        x0 = x0z0(k,1)
        z0 = x0z0(k,2)

        Call DISL2D_DISPL(Ut(k),Ub(k),x0,z0,delta(k),Hd(k),x,z,Ux_k,Uz_k)
        
        Ux=Ux+Ux_k
        Uz=Uz+Uz_k
        
      EndDo

      End Subroutine CRACK2D_DISPL



      Subroutine CRACK2D_STRSS(Ut,Ub,x0z0,delta,Hd,n_dis,x,z,Sxx,Sxz,Szz)

      IMPLICIT NONE
! input: Ut tensile component of the burger vector of each element of the crack (opening)
!        Ub shear component of the burger vector of each element of the crack (slip)
!        x0,z0 coordinates of each dislocation center (x horizontal, z vertical downward)
!        delta angle between the dislocation and z-axis of each element of the crack (for vertical disl delta=0) positive counterclockwise
!        Hd lenght of each dislocation element
!        x,z coordinates of the point in which the stress will be calculated
! output: Sxx,Sxz, Szz x and z components of the stress induced by the boundary element crack
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:),x,z,Sxx,Sxz,Szz
      INTEGER :: n_dis
      REAL(8) :: x0,z0,Sxx_k,Sxz_k,Szz_k
      INTEGER :: k

      Sxx=0.D0
      Sxz=0.D0
      Szz=0.D0
      
      Do k=1,n_dis

        x0 = x0z0(k,1)
        z0 = x0z0(k,2)

        Call DISL2D_STRSS(Ut(k),Ub(k),x0,z0,delta(k),Hd(k),x,z,Sxx_k,Sxz_k,Szz_k)
        
        Sxx=Sxx+Sxx_k
        Sxz=Sxz+Sxz_k
        Szz=Szz+Szz_k
        
      EndDo

      End Subroutine CRACK2D_STRSS


      Subroutine CALC_Dx0Dz0(Ut,Ub,x0z0,delta,Hd,n_dis,Dx0Dz0)

      IMPLICIT NONE
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:),Dx0Dz0(:,:)
      INTEGER :: n_dis

      REAL(8) :: x0,z0,Ux,Uz
      INTEGER :: i

      Ux=0.D0
      Uz=0.D0
      Do i=1,n_dis

          x0 = x0z0(i,1)
          z0 = x0z0(i,2)

          Call CRACK2D_DISPL(Ut,Ub,x0z0,delta,Hd,n_dis,x0,z0,Ux,Uz)
           
          Dx0Dz0(i,1) = Ux-( Ut(i)*cos(delta(i))+Ub(i)*sin(delta(i)) )*0.5D0
          Dx0Dz0(i,2) = Uz-(-Ut(i)*sin(delta(i))+Ub(i)*cos(delta(i)) )*0.5D0

      EndDo

      End Subroutine CALC_Dx0Dz0


      End Module CRACK2D
      
