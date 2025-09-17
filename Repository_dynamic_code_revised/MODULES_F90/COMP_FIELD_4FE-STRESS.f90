      Module COMP_FIELD
      Use BELEMENT,       Only : COMP_LITH_P, stress_rot
      Use EXTERNAL_FIELD, Only : EX_STRESS
      Use CRACK2D
      
      REAL(8), PRIVATE :: ax,az,bx,bz,delta_out,scl
      INTEGER, PRIVATE :: n_pnt
      LOGICAL, PUBLIC  :: movie
      LOGICAL, PRIVATE :: comp_dyke_stress,dikeRF_output
      
      REAL(8), PRIVATE :: pi=4.D0*datan(1.D0)

      contains

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_OUTPUT(ax_,az_,bx_,bz_,n_pnt_,movie_,comp_dyke_stress_,dikeRF_output_,delta_out_,scl_)

      IMPLICIT NONE
      REAL(8) :: ax_,az_,bx_,bz_,delta_out_,scl_
      INTEGER :: n_pnt_
      LOGICAL :: movie_,comp_dyke_stress_,dikeRF_output_
      
      ax = ax_
      az = az_
      bx = bx_
      bz = bz_
      
      n_pnt = n_pnt_
      movie = movie_
      
      comp_dyke_stress = comp_dyke_stress_
      dikeRF_output    = dikeRF_output_
      delta_out        = delta_out_
      scl              = scl_
      
      If (movie.eqv..true.) Then
        write (*,*) 'In your output folder, create a subfolder named "output_frames_r[run number]"'
        write (*,*) 'where run_num is 2 digit integer'
        write(*,'(" [PRESS ENTER TO CONTINUE] ")')
        read(*,*)
      EndIf

      End Subroutine SET_OUTPUT
            
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine DISPL_FIELD(Ut,Ub,x0z0,delta,Hd,n_dis,n_outf1,n_outf2)

      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:)
      INTEGER :: n_dis
      INTEGER :: n_outf1,n_outf2

      REAL(8) :: x,z,xstep,zstep,Ux,Uz
      INTEGER :: n_output,i,j
      CHARACTER(50) :: outf01
      CHARACTER(2)  :: n_outf01
      CHARACTER(3)  :: n_outf02
      
      If (comp_dyke_stress.eqv..true.) Then

      Write(n_outf01,20) n_outf1
      Write(n_outf02,21) n_outf2
      n_output=30+n_outf1
      outf01 = './output/xzUxUz_r'//n_outf01//'_i'//n_outf02//'.dat'
      Open(n_output,file = outf01,access='sequential',form='formatted',status='new')
   
      xstep=(bx-ax)/dble(n_pnt)
      zstep=(bz-az)/dble(n_pnt)

      x=ax-xstep
      Do i=1,n_pnt+1
        x =  x + xstep
        z = az - zstep
        Do j=1,n_pnt+1
           z = z + zstep

           Call CRACK2D_DISPL(Ut,Ub,x0z0,delta,Hd,n_dis,x,z,Ux,Uz)

           Write (n_output,22) x,z,Ux,Uz

         EndDo
      EndDo

20    Format (i2.2)
21    Format (i3.3)
22    Format (2x,5e16.6)

      Close(n_output)
      
      EndIf

      End Subroutine DISPL_FIELD

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine STRSS_FIELD(Ut,Ub,x0z0,delta,Hd,n_dis,n_outf1,n_outf2)

      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:)
      INTEGER :: n_dis
      INTEGER :: n_outf1,n_outf2

      REAL(8) :: x,z,xstep,zstep,Sxx,Sxz,Szz,Stress(2,2),StressR(2,2)
      REAL(8) :: delta_
      INTEGER :: n_output,i,j
      CHARACTER(50) :: outf01
      CHARACTER(2)  :: n_outf01
      CHARACTER(3)  :: n_outf02
      
      If (comp_dyke_stress.eqv..true.) Then

      Write(n_outf01,20) n_outf1
      Write(n_outf02,21) n_outf2
      n_output=30+n_outf1
      outf01 = './output/xzSxxSxzSzz_r'//n_outf01//'_i'//n_outf02//'.dat'
      Open(n_output,file = outf01,access='sequential',form='formatted',status='new')

      If (dikeRF_output.eqv..true.) Then
        delta_ = atan( (x0z0(1,1)-x0z0(n_dis,1)) / (x0z0(1,2)-x0z0(n_dis,2)) )
      Else
        delta_ = delta_out
      EndIf
      
      xstep=(bx-ax)/dble(n_pnt)
      zstep=(bz-az)/dble(n_pnt)

      x=ax-xstep
      Do i=1,n_pnt+1
        x =  x + xstep
        z = az - zstep
        Do j=1,n_pnt+1
           z = z + zstep

           Call CRACK2D_STRSS(Ut,Ub,x0z0,delta,Hd,n_dis,x,z,Sxx,Sxz,Szz)
           
           Stress(1,1)=Sxx
           Stress(1,2)=Sxz
           Stress(2,1)=Stress(1,2)
           Stress(2,2)=Szz

           Call stress_rot(Stress,delta_,StressR)

           Write (n_output,22) x,z,StressR(1,1),StressR(1,2),StressR(2,2)

         EndDo
      EndDo

20    Format (i2.2)
21    Format (i3.3)
22    Format (2x,5e16.6)

      Close(n_output)
      
      EndIf

      End Subroutine STRSS_FIELD

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine CRACK_OVERP(Ut,Ub,x0z0,delta,Hd,n_dis,OverP_av)

      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:),OverP_av
      INTEGER :: n_dis

      REAL(8) :: Sxx,Sxz,Szz,Stress(2,2),StressR(2,2)
      INTEGER :: i

      OverP_av = 0.D0

      Do i=1,n_dis

           Call CRACK2D_STRSS(Ut,Ub,x0z0,delta,Hd,n_dis,x0z0(i,1),x0z0(i,2),Sxx,Sxz,Szz)
           
           Stress(1,1)=Sxx
           Stress(1,2)=Sxz
           Stress(2,1)=Stress(1,2)
           Stress(2,2)=Szz

           Call stress_rot(Stress,delta(i),StressR)

           OverP_av = OverP_av + StressR(1,1)

      EndDo

      OverP_av = OverP_av / n_dis

      End Subroutine CRACK_OVERP

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd,n_dis,n_outf1,n_outf2)

      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd(:)
      INTEGER :: n_dis
      INTEGER :: n_outf1,n_outf2

      REAL(8) :: x_fr,z_fr,Dx0Dz0(n_dis,2)
      INTEGER :: n_output,l
      CHARACTER(50) :: outf01
      CHARACTER(2)  :: n_outf01
      CHARACTER(5)  :: n_outf02

      Write(n_outf01,20) n_outf1
      Write(n_outf02,21) n_outf2
      n_output=30+n_outf1
      If (movie.eqv..true.) Then
        outf01 = 'output/output_frames_r'//n_outf01//'/crack_shape_'//n_outf02//'.dat'
      Else
        outf01 = 'output/crack_shape_r'//n_outf01//'_i'//n_outf02//'.dat'
      EndIf
      Open(n_output,file = outf01,access='sequential',form='formatted',status='new')
      
      Call CALC_Dx0Dz0(Ut,Ub,x0z0,delta,Hd,n_dis,Dx0Dz0)
!~         Dx0Dz0=0.D0

        write (n_output,22) x0z0(1,1)+Dx0Dz0(1,1)*1.d-3*scl+   &
     &                      Hd(1)*sin(delta(1))*0.5d0,         &
     &                      x0z0(1,2)+Dx0Dz0(1,2)*1.d-3*scl+   &
     &                      Hd(1)*cos(delta(1))*0.5d0

        Do l=1,n_dis

           x_fr = x0z0(l,1) +                      &
     &           (-Ut(l)*cos(delta(l))*0.5d0       &
     &            -Ub(l)*sin(delta(l))*0.5d0       &
     &            +Dx0Dz0(l,1)       )*1.d-3*scl
           z_fr = x0z0(l,2) +                      &
     &           (+Ut(l)*sin(delta(l))*0.5d0       &
     &            -Ub(l)*cos(delta(l))*0.5d0       &
     &            +Dx0Dz0(l,2)         )*1.d-3*scl
        write (n_output,22) x_fr, z_fr

        EndDo

        write (n_output,22) x0z0(n_dis,1)+Dx0Dz0(n_dis,1)*1.d-3*scl- &
     &                      Hd(n_dis)*sin(delta(n_dis))*0.5d0,       &
     &                      x0z0(n_dis,2)+Dx0Dz0(n_dis,2)*1.d-3*scl- &
     &                      Hd(n_dis)*cos(delta(n_dis))*0.5d0

        Do l=0,n_dis-1
           x_fr = x0z0(n_dis-l,1) +                              &
     &           (+Ut(n_dis-l)*cos(delta(n_dis-l))*0.5d0         &
     &            +Ub(n_dis-l)*sin(delta(n_dis-l))*0.5d0         &
     &            +Dx0Dz0(n_dis-l,1)               )*1.d-3*scl
           z_fr = x0z0(n_dis-l,2) +                              &
     &           (-Ut(n_dis-l)*sin(delta(n_dis-l))*0.5d0         &
     &            +Ub(n_dis-l)*cos(delta(n_dis-l))*0.5d0         &
     &            +Dx0Dz0(n_dis-l,2)               )*1.d-3*scl
        write (n_output,22) x_fr, z_fr
        EndDo

        write (n_output,22) x0z0(1,1)+Dx0Dz0(1,1)*1.d-3*scl+     &
     &                      Hd(1)*sin(delta(1))*0.5d0,           &
     &                      x0z0(1,2)+Dx0Dz0(1,2)*1.d-3*scl+     &
     &                      Hd(1)*cos(delta(1))*0.5d0
         
20    Format (i2.2)
21    Format (i5.5)
22    Format (2x,6e16.6)

      Close(n_output)

      End Subroutine CRACK_SHAPE

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMPUTE_LITH_P_PROFILE(z_min,z_max,grid_step)

      IMPLICIT NONE

      REAL(8) :: z_min,z_max,grid_step
      
      REAL(8) :: z,Plit
      INTEGER :: i,n_z
      CHARACTER(50) :: outf01

      outf01 = './output/zPlit.dat'
      Open(50,file = outf01,access='sequential',form='formatted',status='new')

      n_z = int((z_max-z_min)/grid_step) - 1
      z = z_min
      Do i=1,n_z

         z = z + grid_step
         Call COMP_LITH_P(Plit,z)

         Write (50,22) z,Plit

      EndDo

22    Format (2x,5e16.6)

      Close(50)

      End Subroutine COMPUTE_LITH_P_PROFILE

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
      
      Subroutine COMPUTE_EX_STRESS_DIRECTIONS(x_min,x_max,z_min,z_max,grid_step)
!~!   compute the max and min stress directions on (x,z) plane
      IMPLICIT NONE
      
      REAL(8) :: x_min,x_max,z_min,z_max,grid_step

      REAL(8) :: S_ex(2,2)
      REAL(8) :: x_i,z_j
      REAL(8) :: EIGENVEC(2,2),EIGENVEC_max(2),EIGENVEC_min(2)
      REAL(8) :: EIGENVAL(2),  EIGENVAL_max,   EIGENVAL_min
      REAL(8) :: alpha_max_compr,alpha_min_compr
      INTEGER :: n_x,n_z,i,j

      REAL(8) :: WORK(68)
      INTEGER :: LWORK=68,INFO

      CHARACTER*50 :: outf01,outf02

      outf01 = './output/inplane_max_compr_stress.dat'
      outf02 = './output/inplane_min_compr_stress.dat'
      Open(51,file = outf01,access='sequential',form='formatted',status='new')
      Open(52,file = outf02,access='sequential',form='formatted',status='new')

      n_x = int((x_max-x_min)/grid_step) - 1
      n_z = int((z_max-z_min)/grid_step) - 1

      x_i = x_min
      Do i=1,n_x
         x_i = x_i + grid_step
         z_j = z_min
         Do j=1,n_z
            z_j = z_j + grid_step

            S_ex  = 0.D0
            Call EX_STRESS(S_ex,x_i,z_j)

            EIGENVEC = S_ex
            Call DSYEV('V','U',2,EIGENVEC,2,EIGENVAL,WORK,LWORK,INFO)

            If (EIGENVAL(1).lt.EIGENVAL(2)) Then
              EIGENVEC_min(:) = EIGENVEC(:,1)
              EIGENVEC_max(:) = EIGENVEC(:,2)
              EIGENVAL_min    = EIGENVAL(1)
              EIGENVAL_max    = EIGENVAL(2)
            Else
              EIGENVEC_min(:) = EIGENVEC(:,2)
              EIGENVEC_max(:) = EIGENVEC(:,1)
              EIGENVAL_min    = EIGENVAL(2)
              EIGENVAL_max    = EIGENVAL(1)
            EndIf
            
            If (EIGENVEC_min(2).ne.0.D0) Then
              alpha_max_compr = (180.D0/pi)*datan(EIGENVEC_min(1)/EIGENVEC_min(2)) !-EIGENVEC_min(1)
            Else
              If (EIGENVEC_min(1).gt.0.D0) Then
                alpha_max_compr =  90.D0
              Else
                alpha_max_compr = -90.D0
              EndIf
            EndIf
            
            alpha_min_compr = alpha_max_compr - dsign(90.D0,alpha_max_compr)

            Write (51,22) x_i, z_j, alpha_max_compr, EIGENVAL_min, EIGENVAL_max-EIGENVAL_min
            Write (52,22) x_i, z_j, alpha_min_compr, EIGENVAL_max, EIGENVAL_max-EIGENVAL_min
            
         EndDo
      EndDo

      Close(51)
      Close(52)

22    Format (2x,5e16.6)
      
      End Subroutine COMPUTE_EX_STRESS_DIRECTIONS

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMPUTE_EX_STRESS_DIRECTIONS_2(x_min,x_max,z_min,z_max,grid_step)
!~!   compute the max and min stress directions on (x,z) plane
      IMPLICIT NONE

      REAL(8) :: x_min,x_max,z_min,z_max,grid_step

      REAL(8) :: S_ex(2,2),S_exR(2,2)
      REAL(8) :: x_i,z_j,delta_k,delta_max,T_k,T_max
      INTEGER :: i,j,n_x,n_z,k
      CHARACTER*50 :: outf01,outf02

      outf01 = './output/xz_sig1.dat'
      outf02 = './output/xz_sig3.dat'
      Open(51,file = outf01,access='sequential',form='formatted',status='new')
      Open(52,file = outf02,access='sequential',form='formatted',status='new')

      n_x = int((x_max-x_min)/grid_step) - 1
      n_z = int((z_max-z_min)/grid_step) - 1

      x_i = x_min
      Do i=1,n_x
         x_i = x_i + grid_step
         z_j = z_min
         Do j=1,n_z
            z_j = z_j + grid_step

            S_ex  = 0.D0
            S_exR = 0.D0
            Call EX_STRESS(S_ex,x_i,z_j)

            Do k=-180,179
               delta_k = dble(k)*pi/360.D0
               If (z_j.gt.0.D0) Call stress_rot(S_ex,delta_k,S_exR)
               T_k = S_exR(1,1)
               If((T_k.gt.T_max+1.D-12).or.(k.eq.-180)) Then
                 delta_max = delta_k
                 T_max = T_k
               EndIf
            EndDo

            If (z_j.gt.0.D0) Call stress_rot(S_ex,delta_max,S_exR)
            Write (51,22) x_i,z_j,(delta_max         )*180.D0/pi, S_exR(2,2)
            Write (52,22) x_i,z_j,(delta_max-pi*0.5d0)*180.D0/pi, S_exR(1,1)

         EndDo
      EndDo
      
      Close(51)
      Close(52)

22    Format (2x,4e16.6)

      End Subroutine COMPUTE_EX_STRESS_DIRECTIONS_2

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

!~       Subroutine stress_rot(stress,delta,stressR)

!~       IMPLICIT NONE
!~       REAL(8) :: stress(:,:),stressR(:,:),delta
!~       REAL(8) :: rot(2,2),rotT(2,2),stress2(2,2)

!~       rot (1,1) =  cos(delta)
!~       rot (1,2) =  sin(delta)
!~       rot (2,1) = -sin(delta)
!~       rot (2,2) =  cos(delta)

!~       rotT(1,1) =  cos(delta)
!~       rotT(1,2) = -sin(delta)
!~       rotT(2,1) =  sin(delta)
!~       rotT(2,2) =  cos(delta)

!~       stress2=matmul(rotT,stress)
!~       stressR=matmul(stress2,rot)

!~       End Subroutine stress_rot


      End Module COMP_FIELD
