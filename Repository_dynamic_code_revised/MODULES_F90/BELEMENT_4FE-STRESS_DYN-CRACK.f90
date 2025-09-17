      Module BELEMENT
      Use EXTERNAL_FIELD
      Use CRACK2D,       Only : CALC_Dx0Dz0
      Use DISL2D,        Only : SET_EPAR, DISL2D_STRSS

!~    PARAMETERS CONSTRAINING THE MAGMA PROPERTIES
      REAL(8), PRIVATE :: Pidr,Kappa,A0,v,eta
!~    PARAMETERS CONSTRAINING THE LITHOSTATIC PRESSURE PROFILE WITHIN THE HOST ROCK
      REAL(8), PRIVATE :: Plit1,Plit2,z_dns,z_ref,Plit1_dns,Plit2_dns,Kr1,Kr2,zl
      LOGICAL, PRIVATE :: Lt_stress_on
!~    PARAMETERS CONSTRAINING THE TECTONIC STRESS WITHIN THE HOST ROCK
      REAL(8), PRIVATE :: SigXX,mx_SXX,mz_SXX,SigZZ,mx_SZZ,mz_SZZ
      LOGICAL, PRIVATE :: Tk_stress_on
      INTEGER, PRIVATE :: n_DYN=0,n_STA=0
      
      REAL(8), PRIVATE :: pi=4.D0*datan(1.D0),g=9.81D0

      REAL(8), ALLOCATABLE, PRIVATE :: F_old(:,:),F_buf(:,:)
!      REAL(8), PRIVATE :: W_closed_el ! wrong

      contains

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_IDR_PAR(rhoM,Kappa_,A0_,v_,eta_)
!~ set magma parameters
      IMPLICIT NONE

      REAL(8) :: rhoM,Kappa_,A0_,v_,eta_

      Pidr  = rhoM*g
      Kappa = Kappa_
      A0    = A0_
      v     = v_
      eta   = eta_

      End Subroutine SET_IDR_PAR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_LITH_PAR(mu1,mu2,nu1,nu2,rhoR1,rhoR2,z_dns_,z_ref_,Lt_stress_on_)
!~ set density parameters for the host rock and call the sub for setting elastic parameters in DISL2D module
      IMPLICIT NONE
      
      
      REAL(8) :: mu1,mu2,nu1,nu2,rhoR1,rhoR2,z_dns_,z_ref_,moho_depth_
      LOGICAL :: Lt_stress_on_

      Call SET_EPAR(mu1,mu2,nu1,nu2)

      Lt_stress_on = Lt_stress_on_
      
      Plit1 = rhoR1*g
      Plit2 = rhoR2*g
      
      z_dns = z_dns_
      z_ref = z_ref_

      Kr1 = 2.D0*mu1*(1.D0+nu1)/(3.D0*(1.D0-2.D0*nu1))

      If (z_dns.eq.0.D0) Then
        If (nu2.eq.0.5D0) Then
          write(*,*) 'WARNING: nu2=0.5 RESULTS IN Kr2=INF'
          write(*,*) '***** nu2 will be set to 0.49 *****'
          write(*,*) '*****[press enter to continue]*****'
          read(*,*)
          Kr2 = 2.D0*mu2*(1.D0+0.49D0)/(3.D0*(1.D0-2.D0*0.49D0))
        Else
          Kr2 = 2.D0*mu2*(1.D0+nu2)/(3.D0*(1.D0-2.D0*nu2))
        EndIf
      Else
        Kr2=Kr1
      EndIf      

      Plit1_dns = Kr1*( exp(Plit1/Kr1*(z_dns-z_ref)) - 1.D0 )
      If (Kr2.gt.0.D0) Then
        Plit2_dns = Kr2*( exp(Plit2/Kr2*(z_dns-z_ref)) - 1.D0 )
      Else
        Plit2_dns = 0.D0
      EndIf
      
      End Subroutine SET_LITH_PAR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_TK_PAR(SigXX_,mx_SXX_,mz_SXX_,SigZZ_,mx_SZZ_,mz_SZZ_,Tk_stress_on_)
!~ set parameters defining the tectonic stress acting on each element of the b. el. crack
      IMPLICIT NONE
      
      REAL(8) :: SigXX_,mx_SXX_,mz_SXX_,SigZZ_,mx_SZZ_,mz_SZZ_
      LOGICAL :: Tk_stress_on_

      Tk_stress_on = Tk_stress_on_

      SigXX  = SigXX_
      mx_SXX = mx_SXX_
      mz_SXX = mz_SXX_
      SigZZ  = SigZZ_
      mx_SZZ = mx_SZZ_
      mz_SZZ = mz_SZZ_

      End Subroutine SET_TK_PAR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_IN_CRACK_GEOM(x0z0,delta,Hd_,n_dis,x_top,z_top,delta_optimal)
!~ set the array x0z0 given delta, Hd_, n_dis, x_top, z_top
!~ when iter=0 delta, Hd, n_dis, x_top, z_top are set to the initial values
      IMPLICIT NONE
      
      REAL(8) :: x0z0(:,:),delta(:),Hd_(:),x_top,z_top
      INTEGER :: n_dis
      LOGICAL :: delta_optimal

      REAL(8) :: S_tk(2,2),S_ex(2,2),S_tkR(2,2),S_exR(2,2)
      REAL(8) :: delta_k,delta_max,T_k,T_max
      INTEGER :: i,k
      
!~ --------------------------------------------------------------------------
!~    COMPUTE THE OPTIMAL ORIENTATION FOR THE INITIAL DIRECTION OF THE DIKE
!~ --------------------------------------------------------------------------
      If (delta_optimal.eqv..true.) Then

        S_tk = 0.D0
        S_ex = 0.D0
        S_tkR = 0.D0
        S_exR = 0.D0
        Call EX_STRESS(S_ex,x_top,z_top)
        If (Tk_stress_on.eqv..true.) Call TK_STRESS(S_tk,x_top,z_top)

        Do k=-180,179
          delta_k = dble(k)*pi/360.D0
          Call stress_rot (S_tk,delta_k,S_tkR)
          Call stress_rot (S_ex,delta_k,S_exR)
          T_k = S_exR(1,1) + S_tkR(1,1)
          If((T_k.gt.T_max+1.D-12).or.(k.eq.-180)) Then
            delta_max = delta_k
            T_max = T_k
          EndIf
        EndDo

        delta = delta_max

      EndIf
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
      
      x0z0(n_dis,1) = x_top + Hd_(n_dis)*0.5D0*sin(delta(n_dis))
      x0z0(n_dis,2) = z_top + Hd_(n_dis)*0.5D0*cos(delta(n_dis))

      Do i=n_dis-1,1,-1

         x0z0(i,1) = x0z0(i+1,1) + (Hd_(i+1)*0.5D0)*sin(delta(i+1))     &
     &                           + (Hd_(i  )*0.5D0)*sin(delta(i  ))

         x0z0(i,2) = x0z0(i+1,2) + (Hd_(i+1)*0.5D0)*cos(delta(i+1))     &
     &                           + (Hd_(i  )*0.5D0)*cos(delta(i  ))

      EndDo

      End Subroutine SET_IN_CRACK_GEOM
      
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SOLVE_CRACK(Ut,Ub,x0z0_in,delta_in,Hd_in,A_iter,n_dis_in,esc)
!~ compute burger vector components on each element of the b. el. crack
!~ and exclude the elements with negative opening given the crack geometry
!~ and calling the sub for calculateing Ut and Ub
!~ Ut(:)        OUTPUT
!~ Ub(:)        OUTPUT
!~ x0z0_in(:,:)  INPUT/OUTPUT
!~ delta_in(:)   INPUT/OUTPUT
!~ Hd_in(:)      INPUT/OUTPUT
!~ A_in          INPUT/OUTPUT
!~ n_dis_in      INPUT/OUTPUT
!~ esc          OUTPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0_in(:,:),delta_in(:),Hd_in(:),A_iter
      INTEGER :: n_dis_in
      LOGICAL :: esc
      
      REAL(8), ALLOCATABLE :: Ub_closed(:),sig(:),tau(:),x0z0(:,:),delta(:),Hd_(:)
      INTEGER :: n_dis,n,l
      
      n_dis = n_dis_in
      
      ALLOCATE (Ub_closed(n_dis),sig(n_dis),tau(n_dis),x0z0(n_dis,2),delta(n_dis),Hd_(n_dis))

      delta  = delta_in
      Hd_    = Hd_in
      x0z0   = x0z0_in
      
10    Call COMP_UtUb(Ut,Ub,x0z0,delta,Hd_,A_iter,n_dis)

      If (Ut(1).le.0.D0) Then
      
        n_dis = n_dis - 1

        Do n=1,n_dis
          delta(n)  = delta(n+1)
          Hd_(n)    = Hd_(n+1)
          x0z0(n,:) = x0z0(n+1,:)
        EndDo

          delta_in(n_dis+1)  = 0.D0
          Hd_in(n_dis+1)     = 0.D0
          x0z0_in(n_dis+1,:) = 0.D0
          Ut(n_dis+1)        = 0.D0
          Ub(n_dis+1)        = 0.D0

        GoTo 10
      EndIf

      Do n=2,n_dis-1
        If (Ut(n).le.0.D0) Then

          write(22,*) ' ** WARNING ** DISLOCATION number', n, '/', n_dis, 'CLOSED'
          write(*,*)  ' ** WARNING ** DISLOCATION number', n, '/', n_dis, 'CLOSED'
          esc = .false.!.true.

          n_dis = n_dis - 1

          Do l=1,n_dis
            delta(l)  = delta(l+1)
            Hd_(l)    = Hd_(l+1)
            x0z0(l,:) = x0z0(l+1,:)
          EndDo

          delta_in(n_dis+1)  = 0.D0
          Hd_in(n_dis+1)     = 0.D0
          x0z0_in(n_dis+1,:) = 0.D0
          Ut(n_dis+1)        = 0.D0
          Ub(n_dis+1)        = 0.D0

          GoTo 10
        EndIf
      EndDo

      If (Ut(n_dis).le.0.D0) Then
        write(*,*) ' ** WARNING ** UPPER DISLOCATION IS CLOSED '
        esc = .true.
      EndIf

      Do n=1,n_dis
        delta_in(n)  = delta(n)
        Hd_in(n)     = Hd_(n)
        x0z0_in(n,:) = x0z0(n,:)
      EndDo
      
      n_dis_in = n_dis
      
      DEALLOCATE (Ub_closed,sig,tau,x0z0,delta,Hd_)

      End subroutine SOLVE_CRACK

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SOLVE_DYN_CRACK(Ut_old,Hd_old,Ut,Ub,DPvisc,x0z0_in,delta_in,Hd_in,A_iter,n_dis_in,esc)
!~ compute burger vector components on each element of the moving crack
!~ and exclude the elements with negative opening given the crack geometry
!~ and calling the sub to compute Ut and Ub
!~ Ut_old(:)     INPUT
!~ Hd_old(:)     INPUT
!~ Ut(:)        OUTPUT
!~ Ub(:)        OUTPUT
!~ DPvisc(:)    OUTPUT
!~ x0z0_in(:,:)  INPUT/OUTPUT
!~ delta_in(:)   INPUT/OUTPUT
!~ Hd_in(:)      INPUT/OUTPUT
!~ A_in          INPUT/OUTPUT
!~ n_dis_in      INPUT/OUTPUT
!~ esc          OUTPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut_old(:),Hd_old(:),Ut(:),Ub(:),DPvisc(:),x0z0_in(:,:),delta_in(:),Hd_in(:),A_iter
      INTEGER :: n_dis_in
      LOGICAL :: esc
      
      REAL(8), ALLOCATABLE :: x0z0(:,:),delta(:),Hd_(:)
      INTEGER :: n_dis,Dn,n,l
      
      CHARACTER(50) :: outf
      CHARACTER(3)  :: n_outf
      CHARACTER(3)  :: DYN_outf

      n_DYN = n_DYN+1
      Write(DYN_outf,89) n_DYN

      n_dis = n_dis_in
      
      
      ALLOCATE (x0z0(n_dis,2),delta(n_dis),Hd_(n_dis))

      delta  = delta_in
      Hd_    = Hd_in
      x0z0   = x0z0_in
      
      l=0
10    Call COMP_DYN_UtUb(Ut_old,Hd_old,Ut,Ub,x0z0,delta,Hd_,A_iter,n_dis)
      l=l+1
       
      Call COMP_DYN_DPvisc(Ut_old,Hd_old,Ut,Hd_,n_dis,DPvisc)
!~       Write(n_outf,89) l
!~       outf = './output/Ut_DYN_'//DYN_outf//'_'//n_outf//'.dat'
!~       Open(99,file = outf,access='sequential',form='formatted',status='new')
!~       Do n=1,n_dis
!~         write(99,*) n, Ut(n), DPvisc(n+l-1)
!~       EndDo
!~       Close(99)

      If (Ut(1).le.0.D0) Then ! or set a threshold for minimum opening
      
!~         Write (*,*) '**********************' 
!~         Write (*,*) 'Ut(1) =',Ut(1) 
!~         Write (*,*) 'PRESS ENTER TO CONTINUE'
!~         Write (*,*) '**********************' 
!~         Read  (*,*)
              
        n_dis = n_dis - 1

        Do n=1,n_dis
          delta(n)  = delta(n+1)
          Hd_(n)    = Hd_(n+1)
          x0z0(n,:) = x0z0(n+1,:)
        EndDo

          delta_in(n_dis+1)  = 0.D0
          Hd_in(n_dis+1)     = 0.D0
          x0z0_in(n_dis+1,:) = 0.D0
          Ut(n_dis+1)        = 0.D0
          Ub(n_dis+1)        = 0.D0

        GoTo 10
      EndIf

      If (minval(Ut(:)).lt.0.D0) Then
        Write (*,*) '***********************************'
        Write (*,*) 'NEAGTIVE OPENING IN SOLVE_DYN_CRACK'          
        Write (*,*) 'Ut_min =',minval(Ut) 
!~         do n=1,n_dis
!~             Write(114,*)n,Ut(n)
!~         enddo
        Write (*,*) 'THE RUN WILL STOP HERE'
        Write (*,*) '***********************************' 
!~         STOP
		esc=.true.
      EndIf

!~       Do n=2,n_dis-1
!~         If (Ut(n).le.0.D0) Then

!~           write(22,*) ' ** WARNING ** DISLOCATION number', n, '/', n_dis, 'CLOSED'
!~           write(*,*)  ' ** WARNING ** DISLOCATION number', n, '/', n_dis, 'CLOSED'
!~           esc = .false.!.true.

!~           n_dis = n_dis - 1

!~           Do l=1,n_dis
!~             delta(l)  = delta(l+1)
!~             Hd_(l)    = Hd_(l+1)
!~             x0z0(l,:) = x0z0(l+1,:)
!~           EndDo

!~           delta_in(n_dis+1)  = 0.D0
!~           Hd_in(n_dis+1)     = 0.D0
!~           x0z0_in(n_dis+1,:) = 0.D0
!~           Ut(n_dis+1)        = 0.D0
!~           Ub(n_dis+1)        = 0.D0

!~           GoTo 10
!~         EndIf
!~       EndDo

!~       If (Ut(n_dis).le.0.D0) Then
!~         write(*,*) ' ** WARNING ** UPPER DISLOCATION IS CLOSED '
!~         esc = .true.
!~       EndIf

      If (size(Ut_old).gt.(n_dis-1)) Then
        Dn = size(Ut_old) - (n_dis-1)
        Do n=1,n_dis
          DPvisc(n) = DPvisc(n+Dn)
        EndDo
        Do n=n_dis+1,size(DPvisc)
          DPvisc(n) = 0.D0
        EndDo
      EndIf

      Do n=1,n_dis
        delta_in(n)  = delta(n)
        Hd_in(n)     = Hd_(n)
        x0z0_in(n,:) = x0z0(n,:)
      EndDo
      
      n_dis_in = n_dis
      
      DEALLOCATE (x0z0,delta,Hd_)

89    Format (i3.3)

      End subroutine SOLVE_DYN_CRACK

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SOLVE_STA_CRACK(Ut,Ub,x0z0_in,z_load_in,delta_in,Hd_in,DPvisc,A_iter,n_dis_in,esc)
!~ compute burger vector components on each element of the stationary moving crack
!~ and exclude the elements with negative opening given the crack geometry
!~ and calling the sub to compute Ut and Ub
!~ Ut(:)         OUTPUT
!~ Ub(:)         OUTPUT
!~ DPvisc(:)     OUTPUT
!~ x0z0_in(:,:)  INPUT/OUTPUT
!~ delta_in(:)   INPUT/OUTPUT
!~ Hd_in(:)      INPUT/OUTPUT
!~ A_in          INPUT/OUTPUT
!~ n_dis_in      INPUT/OUTPUT
!~ esc           OUTPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),DPvisc(:),x0z0_in(:,:),delta_in(:),Hd_in(:),A_iter
      REAL(8) :: z_load_in
      INTEGER :: n_dis_in
      LOGICAL :: esc
      
      REAL(8), ALLOCATABLE :: x0z0(:,:),delta(:),Hd_(:),Ut_iter(:)
      INTEGER :: n_dis
      REAL(8), ALLOCATABLE :: x0z0_STA(:,:),delta_STA(:),Hd_STA(:)
      INTEGER :: n_dis_STA
      INTEGER :: n,l
      LOGICAL :: Ut1,Ut2
      CHARACTER(50) :: outf
      CHARACTER(3)  :: n_outf
      CHARACTER(3)  :: STA_outf

      n_STA = n_STA+1
      Write(STA_outf,89) n_STA

      n_dis = n_dis_in
      
      ALLOCATE (x0z0(n_dis,2),delta(n_dis),Hd_(n_dis),Ut_iter(n_dis))
      ALLOCATE (x0z0_STA(n_dis,2),delta_STA(n_dis),Hd_STA(n_dis))

!***** COMPUTE THE STATIC CRACK SHAPE ******
!*** CALLING COMP_STA_UtUb with DPvisc=0 ***
      delta  = delta_in
      Hd_    = Hd_in
      x0z0   = x0z0_in
      DPvisc = 0.D0
      zl    = z_load_in
      
10    Call COMP_STA_UtUb(Ut,Ub,DPvisc,x0z0,delta,Hd_,A_iter,n_dis)

      If (Ut(1).le.0.D0) Then
      
        n_dis = n_dis - 1

        Do n=1,n_dis
          delta(n)  = delta(n+1)
          Hd_(n)    = Hd_(n+1)
          x0z0(n,:) = x0z0(n+1,:)
        EndDo

          delta(n_dis+1)  = 0.D0
          Hd_(n_dis+1)    = 0.D0
          x0z0(n_dis+1,:) = 0.D0
          Ut(n_dis+1)     = 0.D0
          Ub(n_dis+1)     = 0.D0
        GoTo 10
      EndIf
      
      outf = './output/Ut_STATIC_'//STA_outf//'.dat'
      Open(99,file = outf,access='sequential',form='formatted',status='new')
      Do n=1,n_dis
        write(99,*) n, Ut(n), DPvisc(n)
      EndDo
      Close(99)
      
!***** END COMPUTING THE STATIC CRACK SHAPE ******

!***** COMPUTE THE STATIONARY CRACK SHAPE ******
!*** CALLING COMP_STA_UtUb and update DPvisc ***

      l = 0
11    Ut_iter(:) = Ut(:)
      l = l+1
      
      n_dis_STA = n_dis
      delta_STA = delta
      Hd_STA    = Hd_
      x0z0_STA  = x0z0

      Call COMP_STA_DPvisc(Ut,Hd_STA,n_dis_STA,DPvisc)
      
      Call DEALLOC_Fold_Fbuf()
12    Call COMP_STA_UtUb(Ut,Ub,DPvisc,x0z0_STA,delta_STA,Hd_STA,A_iter,n_dis_STA)

      If (Ut(1).le.0.D0) Then
      
        n_dis_STA = n_dis_STA - 1

        Do n=1,n_dis_STA
          delta_STA(n)  = delta_STA(n+1)
          Hd_STA(n)     = Hd_STA(n+1)
          x0z0_STA(n,:) = x0z0_STA(n+1,:)
        EndDo

          delta_STA(n_dis_STA+1)  = 0.D0
          Hd_STA(n_dis_STA+1)     = 0.D0
          x0z0_STA(n_dis_STA+1,:) = 0.D0
          Ut(n_dis_STA+1)         = 0.D0
          Ub(n_dis_STA+1)         = 0.D0

        GoTo 12
      EndIf

      If (n_dis-n_dis_STA.gt.0) Then
        Do n = n_dis,1+n_dis-n_dis_STA,-1
          Ut(n) = Ut(n-n_dis+n_dis_STA) 
          Ub(n) = Ub(n-n_dis+n_dis_STA) 
        EndDo
        Do n = n_dis-n_dis_STA,1,-1
          Ut(n) = 0.D0
          Ub(n) = 0.D0 
        EndDo

      EndIf

      Write(n_outf,89) l
      outf = './output/Ut_STA_'//STA_outf//'_'//n_outf//'.dat'
      Open(99,file = outf,access='sequential',form='formatted',status='new')
      Do n=1,n_dis
        write(99,*) n, Ut(n), DPvisc(n)
      EndDo
      Close(99)
      If (minval(Ut(:)).lt.0.D0) Then
        Write (*,*) '***********************************'
        Write (*,*) 'NEAGTIVE OPENING IN SOLVE_STA_CRACK'          
        Write (*,*) 'Ut_min =',minval(Ut) 
        Write (*,*) 'THE RUN WILL STOP HERE'
        Write (*,*) '***********************************' 
        STOP
      EndIf
      
      If (l.lt.30) GoTo 11 !use Ut_iter for the convergence condition here
      
!***** END COMPUTING THE STATIONARY CRACK SHAPE ******

13    If (Ut(1).le.0.D0) Then
      
        n_dis = n_dis - 1

        Do n=1,n_dis
          delta(n)  = delta(n+1)
          Hd_(n)    = Hd_(n+1)
          x0z0(n,:) = x0z0(n+1,:)
          Ut(n)     = Ut(n+1)
          Ub(n)     = Ub(n+1)
          DPvisc(n) = DPvisc(n+1)
        EndDo

        delta(n_dis_in+1)  = 0.D0
        Hd_(n_dis_in+1)    = 0.D0
        x0z0(n_dis_in+1,:) = 0.D0
        Ut(n_dis_in+1)     = 0.D0
        Ub(n_dis_in+1)     = 0.D0
        DPvisc(n_dis_in+1) = 0.D0     

        GoTo 13
      EndIf
      
      n_dis_in = n_dis
      delta_in = delta
      Hd_in    = Hd_
      x0z0_in  = x0z0
            
      DEALLOCATE (x0z0,delta,Hd_,x0z0_STA,delta_STA,Hd_STA,Ut_iter)

89    Format (i3.3)

      End subroutine SOLVE_STA_CRACK

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_STA_DPvisc(Ut,Hd_,n_dis,DPvisc)
!compute the viscous pressure change within a stationary moving crack 
!Ut(:)     INPUT
!Hd_(:)    INPUT
!n_dis     INPUT
!DPvisc(:) OUTPUT

      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Hd_(:),DPvisc(:)
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: INSUM(:)
      INTEGER :: j,i
      
      ALLOCATE (INSUM(n_dis))
      
      Ut(:) = Ut(:)*1.D-3
      
      INSUM(1) = -Ut(1)
      Do j=2, n_dis
        INSUM(j) = INSUM(j-1) + (Ut(j-1) - Ut(j))
      EndDo

      If (Ut(1).gt.0.D0) Then
        DPvisc(1) = ( (12.D0*v*eta*Hd_(1)) / (Ut(1)**3) ) * INSUM(1)
      Else
        DPvisc(1) = 0.D0
      EndIf
      Do i=2, n_dis
        If (Ut(i).gt.0.D0) Then
          DPvisc(i) = DPvisc(i-1) + ( (12.D0*v*eta*Hd_(i)) / (Ut(i)**3) ) * INSUM(i)
        Else
          DPvisc(i) = DPvisc(i-1)
        EndIf
      EndDo

      Ut(:) = Ut(:)*1.D3

      End Subroutine COMP_STA_DPvisc

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

!~       Subroutine SOLVE_DYN_CRACK_OLD(Ut_old,DPvisc_old,Ut,Ub,x0z0_in,delta_in,Hd_in,DPvisc,A_iter,n_dis_in,esc)
!~ !!~ compute burger vector components on each element of the b. el. crack
!~ !!~ and exclude the elements with negative opening given the crack geometry
!~ !!~ and calling the sub for calculateing Ut and Ub
!~ !!~ Ut_old(:)     INPUT
!~ !!~ DPvisc_old(:) INPUT
!~ !!~ Ut(:)        OUTPUT
!~ !!~ Ub(:)        OUTPUT
!~ !!~ x0z0_in(:,:)  INPUT/OUTPUT
!~ !!~ delta_in(:)   INPUT/OUTPUT
!~ !!~ Hd_in(:)      INPUT/OUTPUT
!~ !!~ DPvisc(:)     OUTPUT
!~ !!~ A_in          INPUT/OUTPUT
!~ !!~ n_dis_in      INPUT/OUTPUT
!~ !!~ esc          OUTPUT
!~       IMPLICIT NONE
      
!~       REAL(8) :: Ut_old(:),DPvisc_old(:),Ut(:),Ub(:),x0z0_in(:,:),delta_in(:),Hd_in(:),DPvisc(:),A_iter
!~       INTEGER :: n_dis_in
!~       LOGICAL :: esc
      
!~       REAL(8), ALLOCATABLE :: x0z0(:,:),delta(:),Hd_(:),Ut_iter(:)
!~       INTEGER :: n_dis,n,l
!~       LOGICAL :: Ut1,Ut2
!~       CHARACTER(50) :: outf
!~       CHARACTER(3)  :: n_outf
!~       CHARACTER(3)  :: DYN_outf

!~       n_DYN = n_DYN+1
!~       Write(DYN_outf,89) n_DYN

!~       n_dis = n_dis_in
      
!~       ALLOCATE (x0z0(n_dis,2),delta(n_dis),Hd_(n_dis),Ut_iter(n_dis))

!~       delta  = delta_in
!~       Hd_    = Hd_in
!~       x0z0   = x0z0_in

!~ !***** COMPUTE THE DYNAMIC CRACK SHAPE ITERATIVELY ******
!~ !***** (FIRST ITERATION WITH DPvisc = DPvisc_old)  ******
!~       l = 0
!~ 11    Ut_iter(:) = Ut(:)
!~       l = l+1

!~       If (l.eq.1) Then
!~         DPvisc(1:n_dis-1) = DPvisc_old(1:n_dis-1)
!~         DPvisc(n_dis)     = DPvisc(n_dis-1)
!~       Else
!~         Call COMP_DYN_DPvisc_OLD(Ut_old,Ut,Hd_,n_dis_in,DPvisc)
!~       EndIf

!~       Call DEALLOC_Fold_Fbuf()
!~ 12    Call COMP_STA_UtUb(Ut,Ub,DPvisc,x0z0,delta,Hd_,A_iter,n_dis)

!~       If (Ut(1).le.0.D0) Then
      
!~         n_dis = n_dis - 1

!~         Do n=1,n_dis
!~           delta(n)  = delta(n+1)
!~           Hd_(n)    = Hd_(n+1)
!~           x0z0(n,:) = x0z0(n+1,:)
!~           DPvisc(n) = DPvisc(n+1)
!~         EndDo

!~           delta(n_dis+1)  = 0.D0
!~           Hd_(n_dis+1)    = 0.D0
!~           x0z0(n_dis+1,:) = 0.D0
!~           Ut(n_dis+1)     = 0.D0
!~           Ub(n_dis+1)     = 0.D0
!~           DPvisc(n_dis+1) = 0.D0

!~         GoTo 12
!~       EndIf

!~       If (n_dis_in-n_dis.gt.0) Then
!~         Do n = n_dis_in,1+n_dis_in-n_dis,-1
!~           Ut(n) = Ut(n-n_dis_in+n_dis) 
!~           Ub(n) = Ub(n-n_dis_in+n_dis) 
!~         EndDo
!~         Do n = n_dis_in-n_dis,1,-1
!~           Ut(n) = 0.D0
!~           Ub(n) = 0.D0 
!~         EndDo

!~         delta  = delta_in
!~         Hd_    = Hd_in
!~         x0z0   = x0z0_in
!~         n_dis  = n_dis_in
!~       EndIf

!~       Write(n_outf,89) l
!~       outf = './output/Ut_DYN_'//DYN_outf//'_'//n_outf//'.dat'
!~       Open(99,file = outf,access='sequential',form='formatted',status='new')
!~       Do n=1,n_dis
!~         write(99,*) n, Ut(n), DPvisc(n)
!~       EndDo
!~       Close(99)
!~       Do n = 1, n_dis-1
!~         Ut1=.false.
!~         If (Ut(n).gt.0.D0)   Ut1=.true.
!~         Ut2=.false.
!~         If (Ut(n+1).gt.0.D0) Ut2=.true.
!~         If ( (Ut1.eqv..true.).and.(Ut2.eqv..false.) ) Then
!~           write(*,*) 'negative opening at element n/n_dis = ', n+1,'/',n_dis
!~           write(*,*) 'Ut(',n,') = ', Ut(n)
!~           write(*,*) 'Ut(',n+1,') = ', Ut(n+1)
!~           read(*,*)
!~         EndIf
!~       EndDo
      
!~       If (l.lt.31) GoTo 11 !use Ut_iter for the convergence condition here
      
!~ !***** END COMPUTING THE DYNAMIC CRACK SHAPE ******
      
!~ 13    If (Ut(1).le.0.D0) Then
      
!~         n_dis_in = n_dis_in - 1

!~         Do n=1,n_dis_in
!~           delta_in(n)  = delta_in(n+1)
!~           Hd_in(n)     = Hd_in(n+1)
!~           x0z0_in(n,:) = x0z0_in(n+1,:)
!~           Ut(n)        = Ut(n+1)
!~           Ub(n)        = Ub(n+1)
!~         EndDo

!~         delta_in(n_dis_in+1)  = 0.D0
!~         Hd_in(n_dis_in+1)     = 0.D0
!~         x0z0_in(n_dis_in+1,:) = 0.D0
!~         Ut(n_dis_in+1)        = 0.D0
!~         Ub(n_dis_in+1)        = 0.D0

!~         GoTo 13
!~       EndIf
      
!~       DEALLOCATE (x0z0,delta,Hd_,Ut_iter)

!~ 89    Format (i3.3)

!~       End subroutine SOLVE_DYN_CRACK_OLD

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_DYN_DPvisc(Ut_old,Hd_old,Ut,Hd_,n_dis,DPvisc)
!compute the viscous pressure change within a moving crack 
!Ut_old(:) INPUT
!Hd_old(:) INPUT
!Ut(:)     INPUT
!Hd_(:)    INPUT
!n_dis     INPUT
!DPvisc(:) OUTPUT

      IMPLICIT NONE
      
      REAL(8) :: Ut_old(:),Hd_old(:),Ut(:),Hd_(:),DPvisc(:)
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: FLUX_disl(:),DPvisc_disl(:),Ut_new(:)
      REAL(8) :: MULT,A,A_old,SUMFLUX,SUMDPvisc
      INTEGER :: n_dis_old,Dn
      INTEGER :: j,i
      
      n_dis_old = size(Ut_old)

      ALLOCATE (FLUX_disl(n_dis_old+1),DPvisc_disl(n_dis_old+1),Ut_new(n_dis_old+1))
      
      Ut_new(:) = 0.D0
      
      MULT = 12.D0*v*eta*1.D6 ! keep the opening in [mm] or [m] and multiply by 1.D6 this coefficient to adjust the dimention, numerically more efficient.

      If (n_dis_old.eq.(n_dis-1)) Then
        Ut_new(1:n_dis) = Ut(1:n_dis)
      Else
        If (n_dis_old.gt.(n_dis-1)) Then
          Dn = n_dis_old - (n_dis-1)
          Do i=1+Dn,n_dis_old+1
            Ut_new(i) = Ut(i-Dn) 
          EndDo
        Else
          Write (*,*) '**********************' 
          Write (*,*) 'n_dis_old is less than n_dis-1' 
          Write (*,*) 'message from SUB COMP_MOV_DPvisc' 
          Write (*,*) 'THE RUN WILL STOP HERE'
          Write (*,*) '**********************' 
          STOP
        EndIf
      EndIf

!***********************************************************************************
!************************** FOR COMPRESSIBLE FLUIDS ********************************
!***********************************************************************************
      A     = 0.D0
      A_old = 0.D0
      Do i=1,n_dis_old
        A     = A     + Ut_new(i) 
        A_old = A_old + Ut_old(i)
      EndDo
      A = A + Ut_new(n_dis_old+1)

      FLUX_disl(1) = 0.5D0*(Ut_new(1) - Ut_old(1)*A/A_old)
      Do i=2, n_dis_old
        SUMFLUX = 0.D0
        Do j=1,i-1
          SUMFLUX = SUMFLUX + (Ut_new(j) - Ut_old(j)*A/A_old)
        EndDo
        FLUX_disl(i) = SUMFLUX + 0.5D0*(Ut_new(i) - Ut_old(i)*A/A_old)
      EndDo
      SUMFLUX = 0.D0
      Do j=1,n_dis_old
        SUMFLUX = SUMFLUX + (Ut_new(j) - Ut_old(j)*A/A_old)
      EndDo
      FLUX_disl(n_dis_old+1) = SUMFLUX + 0.5D0*Ut_new(n_dis_old+1)

      Do i=1,n_dis_old
        DPvisc_disl(i) = (MULT*Hd_old(i)/Ut_old(i)**3)*FLUX_disl(i)
      EndDo
      DPvisc_disl(n_dis_old+1) = (MULT*Hd_old(n_dis_old)/Ut_old(n_dis_old)**3)*FLUX_disl(n_dis_old+1) ! Here I made the approx of using Ut_old(n_dis_old) when i=n_dis_old+1 in order to be consistent with the approx made in fill_DYN_F (COEFD1)
      
      DPvisc(1) = 0.5D0*DPvisc_disl(1)
      Do i=2,n_dis_old+1
        SUMDPvisc = 0.D0
        Do j=1,i-1
          SUMDPvisc = SUMDPvisc + DPvisc_disl(j)
        EndDo
        DPvisc(i) = SUMDPvisc + 0.5D0*DPvisc_disl(i)
      EndDo

      End Subroutine COMP_DYN_DPvisc

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------


      Subroutine COMP_UtUb(Ut,Ub,x0z0,delta,Hd_,A_iter,n_dis)
!~ !compute burger vector components on each element of the b. el. crack
!~ !given the crack geometry and calling the sub for calculateing the
!~ !stress constrains on each crack element (invert a linear system)
!~ !Ut(:)     OUTPUT
!~ !Ub(:)     OUTPUT
!~ !x0z0(:,:)  INPUT
!~ !delta(:)   INPUT
!~ !Hd_(:)     INPUT
!~ !A_iter     INPUT/OUTPUT
!~ !n_dis      INPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),A_iter
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: F(:,:),sig(:),DeZ(:),tau(:),B(:),X(:)
      REAL(8) :: A_old,z_top
      INTEGER, ALLOCATABLE :: IPIV(:)
      INTEGER :: INFO,counter,i

!~       INTEGER :: DT,COUNTI,COUNTF,COUNT_RATE
      
!~       CALL SYSTEM_CLOCK(COUNTI, COUNT_RATE)

      ALLOCATE (F(n_dis*2,n_dis*2),sig(n_dis),DeZ(n_dis),tau(n_dis))
      ALLOCATE (IPIV(n_dis*2),B(n_dis*2),X(n_dis*2))
      
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)
      Call fill_F (F,x0z0,delta,Hd_,n_dis)
!~       Call fill_F_old (F,x0z0,delta,Hd_,n_dis)

      Call DGETRF (n_dis*2,n_dis*2,F,n_dis*2,IPIV,INFO)
      If (info.ne.0) then
      write(*,*)'info DGETRF = ', info
      STOP
      EndIf

!---------------------------------------------!
!------------COMPUTE BURGER VECT--------------!
!---------------------------------------------!
      z_top = (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))
      DeZ(1:n_dis) = x0z0(1:n_dis,2) - z_top
      
      B(n_dis+1:n_dis*2) = -tau(1:n_dis) ! MPa
      
      counter = 0
10    B(1:n_dis)   = -(sig(1:n_dis) + Kappa + Pidr*DeZ(1:n_dis)*(A0/A_iter)) ! MPa
      X(1:n_dis*2) = B(1:n_dis*2)*1.D3 ! to have opening in meters
      
      Call DGETRS ('N',n_dis*2,1,F,n_dis*2,IPIV,X,n_dis*2,INFO)
      If (info.ne.0) then
      write(*,*)'info DGETRI = ', info
      STOP
      EndIf      
      
      A_old = A_iter
      A_iter= 0.D0
      Do i=1,n_dis
      A_iter = A_iter + Hd_(i)*X(i)*1.D-3 ![km^2 <- X(i) is in meters]
      EndDo

      counter = counter+1

      If ( (abs((A_old-A_iter)/A_old).gt.1.D-10).and.(counter.le.10) ) GoTo 10

      If (counter.gt.10) write (*,*) '** WARNING ** counter>10 in sub COMP_UtUb'

      Ut(1:n_dis) = X(1:n_dis)
      Ub(1:n_dis) = X(n_dis+1:2*n_dis)

      DEALLOCATE (F,sig,DeZ,tau,IPIV,B,X)

!~       CALL SYSTEM_CLOCK(COUNTF)

!~       DT=COUNTF-COUNTI
      
!~       Write(*,*) 'inversion time =', DT,' count rate =', COUNT_RATE

      End Subroutine COMP_UtUb
!

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_DYN_UtUb(Ut_old,Hd_old,Ut,Ub,x0z0,delta,Hd_,A_iter,n_dis)
!~ !compute burger vector components on each element of the moving crack
!~ !given the crack geometry and calling the sub to compute the
!~ !stress constrains on each crack element (invert a linear system)
!~ !Ut_old(:)  INPUT
!~ !Hd_old(:)  INPUT
!~ !Ut(:)     OUTPUT
!~ !Ub(:)     OUTPUT
!~ !x0z0(:,:)  INPUT
!~ !delta(:)   INPUT
!~ !Hd_(:)     INPUT
!~ !A_iter     INPUT/OUTPUT
!~ !n_dis      INPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut_old(:),Hd_old(:),Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),A_iter
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: F(:,:),sig(:),DeZ(:),tau(:),B(:),X(:)
      REAL(8) :: A_iter_old,z_top
      INTEGER, ALLOCATABLE :: IPIV(:)
      INTEGER :: INFO,counter,i,j

      REAL(8), ALLOCATABLE :: COEF1(:,:),COEF2(:)
      REAL(8) :: MULT,COEF2_j,A_old

!~       INTEGER :: DT,COUNTI,COUNTF,COUNT_RATE
      
!~       CALL SYSTEM_CLOCK(COUNTI, COUNT_RATE)

      ALLOCATE (F(n_dis*2,n_dis*2),sig(n_dis),DeZ(n_dis),tau(n_dis))
      ALLOCATE (IPIV(n_dis*2),B(n_dis*2),X(n_dis*2))
      
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)
      Call fill_DYN_F (F,x0z0,delta,Hd_,Ut_old,Hd_old,n_dis)

      Call DGETRF (n_dis*2,n_dis*2,F,n_dis*2,IPIV,INFO)
      If (info.ne.0) then
      write(*,*)'info DGETRF = ', info
      STOP
      EndIf

!---------------------------------------------!
!------------COMPUTE BURGER VECT--------------!
!---------------------------------------------!
      z_top = (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))
      DeZ(1:n_dis) = x0z0(1:n_dis,2) - z_top
      
      B(n_dis+1:n_dis*2) = -tau(1:n_dis) ! MPa
      
      counter = 0
10    B(1:n_dis)   = -( sig(1:n_dis) + Kappa + Pidr*DeZ(1:n_dis)*(A0/A_iter) ) ! MPa
      X(1:n_dis*2) = B(1:n_dis*2)*1.D3 ! to have opening in meters
      
      Call DGETRS ('N',n_dis*2,1,F,n_dis*2,IPIV,X,n_dis*2,INFO)
      If (info.ne.0) then
      write(*,*)'info DGETRI = ', info
      STOP
      EndIf      
      
      A_iter_old = A_iter
      A_iter= 0.D0
      Do i=1,n_dis
      A_iter = A_iter + Hd_(i)*X(i)*1.D-3 ![km^2 <- X(i) is in meters]
      EndDo

      counter = counter+1

      If ( (abs((A_iter_old-A_iter)/A_iter_old).gt.1.D-6).and.(counter.le.100) ) GoTo 10

      If (counter.gt.100) Then
        write (*,*) '** WARNING ** counter>100 in sub COMP_UtUb'
!~         read (*,*)
      EndIf

      Ut(1:n_dis) = X(1:n_dis)
      Ub(1:n_dis) = X(n_dis+1:2*n_dis)
      
      DEALLOCATE (F,sig,DeZ,tau,IPIV,B,X)

!~       CALL SYSTEM_CLOCK(COUNTF)

!~       DT=COUNTF-COUNTI
      
!~       Write(*,*) 'inversion time =', DT,' count rate =', COUNT_RATE

      End Subroutine COMP_DYN_UtUb

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_STA_UtUb(Ut,Ub,DPvisc,x0z0,delta,Hd_,A_iter,n_dis)
!~ !compute burger vector components on each element of the stationary moving crack
!~ !given the crack geometry and calling the sub to compute the
!~ !stress constrains on each crack element (invert a linear system)
!~ !Ut(:)     OUTPUT
!~ !Ub(:)     OUTPUT
!~ !DPvisc(:)  INPUT
!~ !x0z0(:,:)  INPUT
!~ !delta(:)   INPUT
!~ !Hd_(:)     INPUT
!~ !A_iter     INPUT/OUTPUT
!~ !n_dis      INPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),DPvisc(:),x0z0(:,:),delta(:),Hd_(:),A_iter
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: F(:,:),sig(:),DeZ(:),tau(:),B(:),X(:)
      REAL(8) :: A_old,z_top
      INTEGER, ALLOCATABLE :: IPIV(:)
      INTEGER :: INFO,counter,i

!~       INTEGER :: DT,COUNTI,COUNTF,COUNT_RATE
      
!~       CALL SYSTEM_CLOCK(COUNTI, COUNT_RATE)

      ALLOCATE (F(n_dis*2,n_dis*2),sig(n_dis),DeZ(n_dis),tau(n_dis))
      ALLOCATE (IPIV(n_dis*2),B(n_dis*2),X(n_dis*2))
      
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)
      Call fill_F (F,x0z0,delta,Hd_,n_dis)
!~       Call fill_F_old (F,x0z0,delta,Hd_,n_dis)
      Call DGETRF (n_dis*2,n_dis*2,F,n_dis*2,IPIV,INFO)
      If (info.ne.0) then
      write(*,*)'info DGETRF = ', info
      STOP
      EndIf

!---------------------------------------------!
!------------COMPUTE BURGER VECT--------------!
!---------------------------------------------!
      z_top = (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))
      DeZ(1:n_dis) = x0z0(1:n_dis,2) - z_top
      
      B(n_dis+1:n_dis*2) = -tau(1:n_dis) ! MPa
      
      counter = 0
10    B(1:n_dis)   = -(sig(1:n_dis) + Kappa + Pidr*DeZ(1:n_dis)*(A0/A_iter) + DPvisc(1:n_dis)) ! MPa
      X(1:n_dis*2) = B(1:n_dis*2)*1.D3 ! to have opening in meters
      
      Call DGETRS ('N',n_dis*2,1,F,n_dis*2,IPIV,X,n_dis*2,INFO)
      If (info.ne.0) then
      write(*,*)'info DGETRI = ', info
      STOP
      EndIf      
      
      A_old = A_iter
      A_iter= 0.D0
      Do i=1,n_dis
      A_iter = A_iter + Hd_(i)*X(i)*1.D-3 ![km^2 <- X(i) is in meters]
      EndDo

      counter = counter+1

      If ( (abs((A_old-A_iter)/A_old).gt.1.D-10).and.(counter.le.100) ) GoTo 10

      If (counter.gt.100) write (*,*) '** WARNING ** counter>100 in sub COMP_DYN_UtUb'

      Ut(1:n_dis) = X(1:n_dis)
      Ub(1:n_dis) = X(n_dis+1:2*n_dis)

      DEALLOCATE (F,sig,DeZ,tau,IPIV,B,X)

!~       CALL SYSTEM_CLOCK(COUNTF)

!~       DT=COUNTF-COUNTI
      
!~       Write(*,*) 'inversion time =', DT,' count rate =', COUNT_RATE

      End Subroutine COMP_STA_UtUb
      
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_DYN_OverP (Ut,Ub,x0z0,DPvisc,delta,Hd_,A,n_dis,OverP,OverP_av)
!~ compute deformation energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ OverP    OUTPUT
!~ OverP_av OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),OverP(:),DPvisc(:)
      REAL(8) :: A,OverP_av
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: sig(:),tau(:)
      REAL(8) :: DePk_iter,Pidr_iter,DeZn
      INTEGER :: n

      DePk_iter =     - Kappa*(A-A0)/A0
      Pidr_iter = Pidr - Pidr*(A-A0)/A
      
      OverP_av = 0.D0

      ALLOCATE (sig(n_dis),tau(n_dis))
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)

      Do n=1,n_dis

        DeZn = x0z0(n,2) - (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))

        OverP(n) = sig(n) + DePk_iter + Pidr_iter*DeZn + DPvisc(n)
        
        OverP_av = OverP_av + OverP(n)

      EndDo

      OverP_av = OverP_av/n_dis
      
      DEALLOCATE (sig,tau)

      End Subroutine COMP_DYN_OverP

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_OverP (Ut,Ub,x0z0,delta,Hd_,A,n_dis,OverP,OverP_av)
!~ compute deformation energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ OverP    OUTPUT
!~ OverP_av OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),OverP(:)
      REAL(8) :: A,OverP_av
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: sig(:),tau(:)
      REAL(8) :: DePk_iter,Pidr_iter,DeZn
      INTEGER :: n

      DePk_iter =     - Kappa*(A-A0)/A0
      Pidr_iter = Pidr - Pidr*(A-A0)/A
      
      OverP_av = 0.D0

      ALLOCATE (sig(n_dis),tau(n_dis))
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)

      Do n=1,n_dis

        DeZn = x0z0(n,2) - (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))

        OverP(n) = sig(n) + DePk_iter + Pidr_iter*DeZn
        
        OverP_av = OverP_av + OverP(n)

      EndDo

      OverP_av = OverP_av/n_dis
      
      DEALLOCATE (sig,tau)

      End Subroutine COMP_OverP

!~ DELTA ---------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_DELTA_iter(x_top,z_top,delta,Hd_,n_dis,z_start,   &
     &                           OverP_av,L_norm,                       &
     &                           dip_av,dip_ndis,dip_max,               &
     &                           Sig_tip,Sig_xz_av,Sig_xz_n,Sig_xz_max, &
     &                           DELTA1,DELTA2,DELTA3)
!~ compute DELTA parameter (according to 2 definitions)
!~ x_top      INPUT
!~ z_top      INPUT
!~ delta      INPUT
!~ Hd_        INPUT
!~ n_dis      INPUT
!~ z_start    INPUT
!~ OverP_av   INPUT
!~ L_norm     OUTPUT
!~ dip_av     OUTPUT
!~ dip_max    OUTPUT
!~ Sig_tip    OUTPUT
!~ Sig_xz_av  OUTPUT
!~ Sig_xz_max OUTPUT
!~ DELTA1     OUTPUT
!~ DELTA2     OUTPUT
      IMPLICIT NONE

      REAL(8) :: delta(:),Hd_(:),Sig_tip(2,2)
      REAL(8) :: x_top,z_top,z_start,OverP_av,L_norm,dip_av,dip_max,dip_ndis, &
     &           Sig_xz_av,Sig_xz_max,Sig_xz_n,DELTA1,DELTA2,DELTA3
      INTEGER :: n_dis

      REAL(8) :: EIGENVEC(2,2),EIGENVAL(2),EIGENVEC_1(2),WORK(68)
      REAL(8) :: dip_EIGENVEC_1,Sig_tip_av(2,2),Sig_tip_max(2,2),Sig_tip_n(2,2)
      INTEGER :: LWORK=68,INFO

      L_norm=n_dis*Hd_(1)/z_start
      Sig_tip=0.D0
      Call EX_STRESS(Sig_tip,x_top,z_top)
      
      dip_av   = sum(delta)/n_dis
      dip_ndis = delta(n_dis)! sum(delta)/n_dis
      
      EIGENVEC = Sig_tip
      Call DSYEV('V','U',2,EIGENVEC,2,EIGENVAL,WORK,LWORK,INFO)
      EIGENVEC_1(:) = EIGENVEC(:,1)
      If (EIGENVEC_1(2).ne.0.D0) Then
          dip_EIGENVEC_1 = datan(EIGENVEC_1(1)/EIGENVEC_1(2))
      Else
        If (EIGENVEC_1(1).gt.0.D0) Then
          dip_EIGENVEC_1 =  pi*0.5D0
        Else
          dip_EIGENVEC_1 = -pi*0.5D0
        EndIf
      EndIf
      dip_max = dip_EIGENVEC_1 - dsign(pi,dip_EIGENVEC_1)*0.25D0

      Call stress_rot(Sig_tip,dip_av,Sig_tip_av)
      Sig_xz_av=abs(Sig_tip_av(1,2))

      Call stress_rot(Sig_tip,dip_ndis,Sig_tip_n)
      Sig_xz_n=abs(Sig_tip_n(1,2))

      Call stress_rot(Sig_tip,dip_max,Sig_tip_max)
      Sig_xz_max=abs(Sig_tip_max(1,2))

      dip_av  = (180.D0/pi)*dip_av
      dip_ndis = (180.D0/pi)*dip_ndis
      dip_max = (180.D0/pi)*dip_max
      
      DELTA1 = Sig_xz_av /(OverP_av*L_norm)
      DELTA2 = Sig_xz_n  /(OverP_av*L_norm)
      DELTA3 = Sig_xz_max/(OverP_av*L_norm)

      End Subroutine COMP_DELTA_iter

!~ END DELTA -----------------------------------------------------------
!~ ---------------------------------------------------------------------

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_U (Ut,Ub,x0z0,delta,Hd_,A,n_dis,U)
!~ compute gravitational energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ U     OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
      REAL(8) :: A,U
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: Dx0Dz0(:,:)
      REAL(8) :: Z_CM
      INTEGER :: i

      ALLOCATE (Dx0Dz0(n_dis,2))
      Call CALC_Dx0Dz0(Ut,Ub,x0z0,delta,Hd_,n_dis,Dx0Dz0)

      Z_CM = 0.D0
      Do i = 1,n_dis
        Z_CM = Z_CM + Hd_(i)*Ut(i)*1.D-3*(x0z0(i,2))!+Dx0Dz0(i,2)*1.D-3)
      EndDo
      Z_CM = Z_CM/A

      U = -Z_CM*(Pidr*A0*1.D3)   ! MPa*km*m

!~       Write (22,*) ' x0z0(',n_dis,'|', delta(n_dis)*180.D0/(4.D0*atan(1.D0)),')', &
!~      &              ' = (', x0z0(n_dis,1),';', x0z0(n_dis,2),')'
      Write (22,*) 'U   =',U
!~       Write (22,*) 'Dz(n_dis = ',n_dis,') =',Dx0Dz0(n_dis,2)
!~       Write (*,*)  'U   =',U

      DEALLOCATE (Dx0Dz0)

      End Subroutine COMP_U

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_W (Ut,Ub,x0z0,delta,Hd_,A,n_dis,W)
!~ compute deformation energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ W     OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
      REAL(8) :: A,W
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: sig(:),tau(:),OverP(:)
      REAL(8) :: W_r,W_f,DePk_iter,Pidr_iter,DeZn
      INTEGER :: n

      W_r = 0.D0
      W_f = 0.D0

      DePk_iter =     - Kappa*(A-A0)/A0
      Pidr_iter = Pidr - Pidr*(A-A0)/A

      ALLOCATE (sig(n_dis),tau(n_dis),OverP(n_dis))
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)

      Do n=1,n_dis

        DeZn = x0z0(n,2) - (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))

        OverP(n) = sig(n) + DePk_iter + Pidr_iter*DeZn

        W_r = W_r - 0.5D0*Hd_(n)*( (-OverP(n) + 2.D0*sig(n))*Ut(n) +   &
     &                                               tau(n) *Ub(n)  )   ! MPa*km*m

      EndDo

      W_f = 0.5D0*Kappa*A*1.D3*((A-A0)/A0)**2    ! MPa*km*m

      W = W_r + W_f ! + W_closed_el

      Write (22,*) 'W_r =',W_r
      Write (22,*) 'W_f =',W_f
!~       Write (*,*) 'W_r =',W_r
!~       Write (*,*) 'W_f =',W_f

      DEALLOCATE (sig,tau,OverP)

      End Subroutine COMP_W

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_DYN_W (Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,W)
!~ compute deformation energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ DPvisc INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ W     OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),DPvisc(:)
      REAL(8) :: A,W
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: sig(:),tau(:),OverP(:)
      REAL(8) :: W_r,W_f,DePk_iter,Pidr_iter,DeZn
      INTEGER :: n

      W_r = 0.D0
      W_f = 0.D0

      DePk_iter =     - Kappa*(A-A0)/A0
      Pidr_iter = Pidr - Pidr*(A-A0)/A

      ALLOCATE (sig(n_dis),tau(n_dis),OverP(n_dis))
      Call STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)

      Do n=1,n_dis

        DeZn = x0z0(n,2) - (x0z0(n_dis,2)+Hd_(n_dis)*0.5D0*cos(delta(n_dis)))

        OverP(n) = sig(n) + DePk_iter + Pidr_iter*DeZn + DPvisc(n)

        W_r = W_r - 0.5D0*Hd_(n)*( (-OverP(n) + 2.D0*sig(n))*Ut(n) +   &
     &                                               tau(n) *Ub(n)  )   ! MPa*km*m

      EndDo

      W_f = 0.5D0*Kappa*A*1.D3*((A-A0)/A0)**2    ! MPa*km*m

      W = W_r + W_f ! + W_closed_el

      Write (22,*) 'W_r =',W_r
      Write (22,*) 'W_f =',W_f
!~       Write (*,*) 'W_r =',W_r
!~       Write (*,*) 'W_f =',W_f

      DEALLOCATE (sig,tau,OverP)

      End Subroutine COMP_DYN_W

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------


      Subroutine COMP_D(Ut_old,Hd_old,Ut,Hd_,n_dis,D)
!compute the viscous dissipation energy 
!Ut_old(:) INPUT
!Hd_old(:) INPUT
!Ut(:)     INPUT
!Hd_(:)    INPUT
!n_dis     INPUT
!D(:)      OUTPUT

      IMPLICIT NONE
      
      REAL(8) :: Ut_old(:),Hd_old(:),Ut(:),Hd_(:),D
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: Ut_new(:),Ut_avr(:),DPvisc(:),dDPvisc_dz(:)
      REAL(8) :: MULT,Hd
      INTEGER :: n_dis_old,Dn
      INTEGER :: j,i
      
      n_dis_old = size(Ut_old)

      ALLOCATE (DPvisc(n_dis_old+1),dDPvisc_dz(n_dis_old+1),Ut_new(n_dis_old+1),Ut_avr(n_dis_old+1))

      Hd = Hd_(n_dis)                   ! Will not work if Hd_(i) is not constant
      MULT = -Hd/(12.D0*eta*1.D3*v*1.D3)! keep the opening in [mm] or [m] and multiply eta and v by 1.D3, numerically more efficient.

      Ut_new(:) = 0.D0
      D         = 0.D0
      
      If (n_dis_old.eq.(n_dis-1)) Then
        Ut_new(1:n_dis) = Ut(1:n_dis)
      Else
        If (n_dis_old.gt.(n_dis-1)) Then
          Dn = n_dis_old - (n_dis-1)
          Do i=1+Dn,n_dis_old+1
            Ut_new(i) = Ut(i-Dn) 
          EndDo
        Else
          Write (*,*) '**********************' 
          Write (*,*) 'n_dis_old is less than n_dis-1' 
          Write (*,*) 'message from SUB COMP_D' 
          Write (*,*) 'THE RUN WILL STOP HERE'
          Write (*,*) '**********************' 
          STOP
        EndIf
      EndIf
      
      Call COMP_DYN_DPvisc(Ut_old,Hd_old,Ut,Hd_,n_dis,DPvisc)
      
      dDPvisc_dz(1) = DPvisc(1)/(Hd*0.5D0)
      Do i=2,n_dis_old+1
        dDPvisc_dz(i) = (DPvisc(i) - DPvisc(i-1))/Hd
      EndDo

      Do i=1,n_dis_old
        Ut_avr(i) = (Ut_new(i) + Ut_old(i))*0.5D0
      EndDo
      Ut_avr(n_dis_old+1) = Ut_new(n_dis_old+1)*0.5D0

      Do i=1,n_dis_old+1
        D = D + Hd*Ut_avr(i)**3*dDPvisc_dz(i)**2
      EndDo
      
      D = D*MULT*1.D-3 ! D is multiplied by 1.D-3 so that its units are [KPa*m^2] or [MPa*km^2]

	  DEALLOCATE (DPvisc,dDPvisc_dz,Ut_new,Ut_avr)

      End Subroutine COMP_D

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------
!~       Subroutine COMP_Ecf(Ut_old,Hd_old,Ut,Hd_,n_dis,Ecf)
!~ !compute the viscous dissipation energy 
!~ !Ut_old(:) INPUT
!~ !Hd_old(:) INPUT
!~ !Ut(:)     INPUT
!~ !Hd_(:)    INPUT
!~ !n_dis     INPUT
!~ !Ecf(:)      OUTPUT

!~       IMPLICIT NONE
      
!~       REAL(8) :: Ut_old(:),Hd_old(:),Ut(:),Hd_(:),Ecf
!~       INTEGER :: n_dis

!~       REAL(8), ALLOCATABLE :: Ut_new(:),Ut_avr(:),DPvisc(:),dDPvisc_dz(:)
!~       REAL(8) :: MULT,Hd
!~       INTEGER :: n_dis_old,Dn
!~       INTEGER :: j,i
      
!~       n_dis_old = size(Ut_old)

!~       ALLOCATE (DPvisc(n_dis_old+1),dDPvisc_dz(n_dis_old+1),Ut_new(n_dis_old+1),Ut_avr(n_dis_old+1))

!~       Hd = Hd_(n_dis)                   ! Will not work if Hd_(i) is not constant
!~       MULT = -Hd/(12.D0*eta*1.D3*v*1.D3)! keep the opening in [mm] or [m] and multiply eta and v by 1.D3, numerically more efficient.

!~       Ut_new(:) = 0.D0
!~       D         = 0.D0
      
!~       If (n_dis_old.eq.(n_dis-1)) Then
!~         Ut_new(1:n_dis) = Ut(1:n_dis)
!~       Else
!~         If (n_dis_old.gt.(n_dis-1)) Then
!~           Dn = n_dis_old - (n_dis-1)
!~           Do i=1+Dn,n_dis_old+1
!~             Ut_new(i) = Ut(i-Dn) 
!~           EndDo
!~         Else
!~           Write (*,*) '**********************' 
!~           Write (*,*) 'n_dis_old is less than n_dis-1' 
!~           Write (*,*) 'message from SUB COMP_D' 
!~           Write (*,*) 'THE RUN WILL STOP HERE'
!~           Write (*,*) '**********************' 
!~           STOP
!~         EndIf
!~       EndIf
      
!~       Call COMP_DYN_DPvisc(Ut_old,Hd_old,Ut,Hd_,n_dis,DPvisc)
      
!~       dDPvisc_dz(1) = DPvisc(1)/(Hd*0.5D0)
!~       Do i=2,n_dis_old+1
!~         dDPvisc_dz(i) = (DPvisc(i) - DPvisc(i-1))/Hd
!~       EndDo

!~       Do i=1,n_dis_old
!~         Ut_avr(i) = (Ut_new(i) + Ut_old(i))*0.5D0
!~       EndDo
!~       Ut_avr(n_dis_old+1) = Ut_new(n_dis_old+1)*0.5D0

!~       Do i=1,n_dis_old+1
!~         D = D + Hd*Ut_avr(i)**3*dDPvisc_dz(i)**2
!~       EndDo
      
!~       D = D*MULT*1.D-3 ! D is multiplied by 1.D-3 so that its units are [KPa*m^2] or [MPa*km^2]

!~ 	  DEALLOCATE (DPvisc,dDPvisc_dz,Ut_new,Ut_avr)

!~       End Subroutine COMP_Ecf

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_Fold_eq_Fbuf
 
      IMPLICIT NONE

      If (ALLOCATED(F_old).eqv..true.) DEALLOCATE (F_old)
      If (ALLOCATED(F_buf).eqv..true.) Then
        ALLOCATE (F_old(size(F_buf(:,1)),size(F_buf(:,1))))
        F_old(:,:) = F_buf(:,:)
      EndIf

      End Subroutine SET_Fold_eq_Fbuf


      Subroutine SET_Fbuf_eq_Fold
 
      IMPLICIT NONE
      
      If (ALLOCATED(F_buf).eqv..true.) DEALLOCATE (F_buf)
      If (ALLOCATED(F_old).eqv..true.) Then
        ALLOCATE (F_buf(size(F_old(:,1)),size(F_old(:,1))))
        F_buf(:,:) = F_old(:,:)
      EndIf

      End Subroutine SET_Fbuf_eq_Fold


      Subroutine DEALLOC_Fold_Fbuf()

      If (ALLOCATED(F_old).eqv..true.) DEALLOCATE (F_old)
      If (ALLOCATED(F_buf).eqv..true.) DEALLOCATE (F_buf)

      End Subroutine DEALLOC_Fold_Fbuf

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine fill_F (F,x0z0,delta,Hd_,n_dis)
!!~ Fill the matrix F with the loading coefficient of the elementary dislocation distribution
!!~ F     OUTPUT
!!~ x0z0  INPUT
!!~ delta INPUT
!!~ Hd_   INPUT
!!~ n_dis INPUT
      IMPLICIT NONE

      REAL(8) :: F(:,:),x0z0(:,:),delta(:),Hd_(:)
      INTEGER :: n_dis

      REAL(8) :: sig_Ut(2,2),sig_UtR(2,2),sig_Ub(2,2),sig_UbR(2,2)
      INTEGER :: i,j
      LOGICAL :: F_buf_on

10    F_buf_on = ALLOCATED(F_buf)
      
!~       write(*,*) F_old_on,size(F_old(:,1)),size(F(:,1))

      If (F_buf_on.eqv..false.) Then

!~       write(*,*) '******* THIS SHOULD HAPPEN ONLY ONCE'
!~       read(*,*)

      Do i=1,n_dis
        Do j=1,n_dis

              Call DISL2D_STRSS(1.D3,0.D0,x0z0(j,1),x0z0(j,2),delta(j), & ! coefficients in MPa/km if Ut (and Ub) = 1.D3, the result will be in km and X in compute_UtUb must be multiplied by 1.D3
     &                             Hd_(j),x0z0(i,1),x0z0(i,2),          & ! otherwise if Ut (and Ub) = 1.D0 the coefficients are in MP/m, the results will be in m BUT (-Kappa/A0)*Hd_(j) must be multiplied by 1.D-3 and in compute_UtUb X must not be multiplied by 1.D3
     &                          sig_Ut(1,1),sig_Ut(1,2),sig_Ut(2,2))
     
              sig_Ut(2,1) = sig_Ut(1,2)

              Call stress_rot(sig_Ut,delta(i),sig_UtR)

              F(i,j)       = sig_UtR(1,1) + (-Kappa/A0)*Hd_(j) ! MPa/km
              F(i+n_dis,j) = sig_UtR(1,2)! + sforz_tR(1,1)*attr(i-n_dis)

        EndDo
      EndDo

      Do i=1,n_dis
        Do j=n_dis+1,n_dis*2

              Call DISL2D_STRSS(0.D0,1.D3,x0z0(j-n_dis,1),x0z0(j-n_dis,2),delta(j-n_dis), &
     &                             Hd_(j-n_dis),x0z0(i,1),x0z0(i,2),          &
     &                          sig_Ub(1,1),sig_Ub(1,2),sig_Ub(2,2))
     
              sig_Ub(2,1) = sig_Ub(1,2)

              Call stress_rot(sig_Ub,delta(i),sig_UbR)

              F(i,j)       = sig_UbR(1,1)
              F(i+n_dis,j) = sig_UbR(1,2)! + sforz_bR(1,1)*attr(i-n_dis)

        EndDo
      EndDo      

      ALLOCATE(F_buf(n_dis*2,n_dis*2))
      F_buf(:,:) = F(:,:)

      Else ! (F_buf_on=.true.)

      If (size(F_buf(:,1)).eq.n_dis*2) Then ! no disl added or closed
      
      write(*,*) '******* THIS SHOULD NEVER HAPPEN 1 (sub fill_F in BELEMENT MODULE)'
      read(*,*)
      F(:,:) = F_buf(:,:)

      Else ! F_buf is loaded
      
      If (size(F_buf(:,1)).eq.n_dis*2+2) Then ! (A) one disl has been closed from the bottom
      
      F(1      :n_dis  , 1      :n_dis  ) = F_buf(2      :n_dis+1  , 2      :n_dis+1  )
      F(1      :n_dis  , n_dis+1:n_dis*2) = F_buf(2      :n_dis+1  , n_dis+3:n_dis*2+2)
      F(n_dis+1:n_dis*2, 1      :n_dis  ) = F_buf(n_dis+3:n_dis*2+2, 2      :n_dis+1  )
      F(n_dis+1:n_dis*2, n_dis+1:n_dis*2) = F_buf(n_dis+3:n_dis*2+2, n_dis+3:n_dis*2+2)
      
      DEALLOCATE (F_buf)
      ALLOCATE (F_buf(n_dis*2,n_dis*2))

      F_buf(:,:) = F(:,:)

      Else ! (A)
      
      If (size(F_buf(:,1)).eq.n_dis*2-2) Then ! (B) one disl has been added at the top
      
      F(1      :n_dis-1  , 1      :n_dis-1  ) = F_buf(1    :n_dis-1  , 1    :n_dis-1  )
      F(1      :n_dis-1  , n_dis+1:n_dis*2-1) = F_buf(1    :n_dis-1  , n_dis:n_dis*2-2)
      F(n_dis+1:n_dis*2-1, 1      :n_dis-1  ) = F_buf(n_dis:n_dis*2-2, 1    :n_dis-1  )
      F(n_dis+1:n_dis*2-1, n_dis+1:n_dis*2-1) = F_buf(n_dis:n_dis*2-2, n_dis:n_dis*2-2)
      
      Do j = 1,n_dis
      
              Call DISL2D_STRSS(1.D3,0.D0,x0z0(j,1),x0z0(j,2),delta(j), &
     &                             Hd_(j),x0z0(n_dis,1),x0z0(n_dis,2),          &
     &                          sig_Ut(1,1),sig_Ut(1,2),sig_Ut(2,2))
     
              sig_Ut(2,1) = sig_Ut(1,2)

              Call stress_rot(sig_Ut,delta(n_dis),sig_UtR)

              F(n_dis,j)   = sig_UtR(1,1) + (-Kappa/A0)*Hd_(j)
              F(n_dis*2,j) = sig_UtR(1,2)! + sforz_tR(1,1)*attr(i-n_dis)
      
      EndDo

      Do j = n_dis+1,n_dis*2
      
              Call DISL2D_STRSS(0.D0,1.D3,x0z0(j-n_dis,1),x0z0(j-n_dis,2),delta(j-n_dis), &
     &                             Hd_(j-n_dis),x0z0(n_dis,1),x0z0(n_dis,2),          &
     &                          sig_Ub(1,1),sig_Ub(1,2),sig_Ub(2,2))
     
              sig_Ub(2,1) = sig_Ub(1,2)

              Call stress_rot(sig_Ub,delta(n_dis),sig_UbR)

              F(n_dis,j)   = sig_UbR(1,1)
              F(n_dis*2,j) = sig_UbR(1,2)! + sforz_bR(1,1)*attr(i-n_dis)
      
      EndDo

      Do i=1,n_dis

              Call DISL2D_STRSS(1.D3,0.D0,x0z0(n_dis,1),x0z0(n_dis,2),delta(n_dis), &
     &                             Hd_(n_dis),x0z0(i,1),x0z0(i,2),          &
     &                          sig_Ut(1,1),sig_Ut(1,2),sig_Ut(2,2))
     
              sig_Ut(2,1) = sig_Ut(1,2)

              Call stress_rot(sig_Ut,delta(i),sig_UtR)

              F(i,n_dis)       = sig_UtR(1,1) + (-Kappa/A0)*Hd_(n_dis)
              F(i+n_dis,n_dis) = sig_UtR(1,2)! + sforz_tR(1,1)*attr(i-n_dis)

      EndDo

      Do i=1,n_dis

              Call DISL2D_STRSS(0.D0,1.D3,x0z0(n_dis,1),x0z0(n_dis,2),delta(n_dis), &
     &                             Hd_(n_dis),x0z0(i,1),x0z0(i,2),          &
     &                          sig_Ub(1,1),sig_Ub(1,2),sig_Ub(2,2))
     
              sig_Ub(2,1) = sig_Ub(1,2)

              Call stress_rot(sig_Ub,delta(i),sig_UbR)

              F(i,n_dis*2)       = sig_UbR(1,1)
              F(i+n_dis,n_dis*2) = sig_UbR(1,2)! + sforz_bR(1,1)*attr(i-n_dis)

      EndDo      

      DEALLOCATE (F_buf)
      ALLOCATE (F_buf(n_dis*2,n_dis*2))
      
      F_buf(:,:) = F(:,:)
      
      Else ! (B) more than one disl closed or added, should not happen if disl at the tail close one by one
      
      write(*,*) '******* THIS SHOULD NEVER HAPPEN 2 (sub fill_F in BELEMENT MODULE)'
      read(*,*)

      DEALLOCATE (F_buf)
      GoTo 10

      EndIf ! (B)
      EndIf ! (A)
      EndIf
      EndIf
      
      End subroutine fill_F

!---------------------------------------------------------------------
!---------------------------------------------------------------------

      Subroutine fill_DYN_F (F,x0z0,delta,Hd_,Ut_old,Hd_old,n_dis)
!!~ Fill the matrix F with the loading coefficient of the elementary dislocations
!!~ it also compute the coefficients due to the viscous flow through the dislocations
!!~ F     OUTPUT
!!~ x0z0   INPUT
!!~ delta  INPUT
!!~ Hd_    INPUT
!!~ Ut_old INPUT
!!~ Hd_old INPUT
!!~ n_dis  INPUT
      IMPLICIT NONE

      REAL(8) :: F(:,:),x0z0(:,:),delta(:),Hd_(:),Ut_old(:),Hd_old(:)
      INTEGER :: n_dis

      REAL(8) :: sig_Ut(2,2),sig_UtR(2,2),sig_Ub(2,2),sig_UbR(2,2)
      INTEGER :: i,j
      LOGICAL :: F_buf_on

      REAL(8), ALLOCATABLE :: COEFA1(:,:),COEFA2(:)
      REAL(8), ALLOCATABLE :: COEFB1(:,:),COEFB2(:)
      REAL(8), ALLOCATABLE :: COEFC1(:,:),COEFC2(:)
      REAL(8), ALLOCATABLE :: COEFD1(:,:),COEFD2(:)
      REAL(8) :: MULT,A_old,SUMA2_j,SUMB2_j,SUMC2_j
      INTEGER :: n_dis_old,Dn

10    F_buf_on = ALLOCATED(F_buf)
      
!~ !      write(*,*) F_old_on,size(F_old(:,1)),size(F(:,1))

      If (F_buf_on.eqv..false.) Then

!~ !      write(*,*) '******* THIS SHOULD HAPPEN ONLY ONCE'
!~ !      read(*,*)

      Do i=1,n_dis
        Do j=1,n_dis

              Call DISL2D_STRSS(1.D3,0.D0,x0z0(j,1),x0z0(j,2),delta(j), & ! coefficients in MPa/km if Ut (and Ub) = 1.D3, the result will be in km and X in compute_UtUb must be multiplied by 1.D3
     &                             Hd_(j),x0z0(i,1),x0z0(i,2),          & ! otherwise if Ut (and Ub) = 1.D0 the coefficients are in MP/m, the results will be in m BUT (-Kappa/A0)*Hd_(j) must be multiplied by 1.D-3 and in compute_UtUb X must not be multiplied by 1.D3
     &                          sig_Ut(1,1),sig_Ut(1,2),sig_Ut(2,2))
     
              sig_Ut(2,1) = sig_Ut(1,2)

              Call stress_rot(sig_Ut,delta(i),sig_UtR)

              F(i,j)       = sig_UtR(1,1) + (-Kappa/A0)*Hd_(j) ! MPa/km
              F(i+n_dis,j) = sig_UtR(1,2)! + sforz_tR(1,1)*attr(i-n_dis)

        EndDo
      EndDo

      Do i=1,n_dis
        Do j=n_dis+1,n_dis*2

              Call DISL2D_STRSS(0.D0,1.D3,x0z0(j-n_dis,1),x0z0(j-n_dis,2),delta(j-n_dis), &
     &                             Hd_(j-n_dis),x0z0(i,1),x0z0(i,2),          &
     &                          sig_Ub(1,1),sig_Ub(1,2),sig_Ub(2,2))
     
              sig_Ub(2,1) = sig_Ub(1,2)

              Call stress_rot(sig_Ub,delta(i),sig_UbR)

              F(i,j)       = sig_UbR(1,1)
              F(i+n_dis,j) = sig_UbR(1,2)! + sforz_bR(1,1)*attr(i-n_dis)

        EndDo
      EndDo
      
      ALLOCATE(F_buf(n_dis*2,n_dis*2))
      F_buf(:,:) = F(:,:)

      Else ! (F_buf_on=.true.)

      If (size(F_buf(:,1)).eq.n_dis*2) Then ! no disl added or closed
      
      write(*,*) '******* THIS SHOULD NEVER HAPPEN 1 (sub fill_F in BELEMENT MODULE)'
      read(*,*)
      
      DEALLOCATE(F_buf)
      GoTo 10
!~ !      F(:,:) = F_buf(:,:)

      Else ! F_buf is loaded
      
      If (size(F_buf(:,1)).eq.n_dis*2+2) Then ! (A) one disl has been closed from the bottom
      
      F(1      :n_dis  , 1      :n_dis  ) = F_buf(2      :n_dis+1  , 2      :n_dis+1  )
      F(1      :n_dis  , n_dis+1:n_dis*2) = F_buf(2      :n_dis+1  , n_dis+3:n_dis*2+2)
      F(n_dis+1:n_dis*2, 1      :n_dis  ) = F_buf(n_dis+3:n_dis*2+2, 2      :n_dis+1  )
      F(n_dis+1:n_dis*2, n_dis+1:n_dis*2) = F_buf(n_dis+3:n_dis*2+2, n_dis+3:n_dis*2+2)
      
      DEALLOCATE (F_buf)
      ALLOCATE (F_buf(n_dis*2,n_dis*2))

      F_buf(:,:) = F(:,:)

      Else ! (A)
      
      If (size(F_buf(:,1)).eq.n_dis*2-2) Then ! (B) one disl has been added at the top
      
      F(1      :n_dis-1  , 1      :n_dis-1  ) = F_buf(1    :n_dis-1  , 1    :n_dis-1  )
      F(1      :n_dis-1  , n_dis+1:n_dis*2-1) = F_buf(1    :n_dis-1  , n_dis:n_dis*2-2)
      F(n_dis+1:n_dis*2-1, 1      :n_dis-1  ) = F_buf(n_dis:n_dis*2-2, 1    :n_dis-1  )
      F(n_dis+1:n_dis*2-1, n_dis+1:n_dis*2-1) = F_buf(n_dis:n_dis*2-2, n_dis:n_dis*2-2)
      
      Do j = 1,n_dis
      
              Call DISL2D_STRSS(1.D3,0.D0,x0z0(j,1),x0z0(j,2),delta(j), &
     &                             Hd_(j),x0z0(n_dis,1),x0z0(n_dis,2),          &
     &                          sig_Ut(1,1),sig_Ut(1,2),sig_Ut(2,2))
     
              sig_Ut(2,1) = sig_Ut(1,2)

              Call stress_rot(sig_Ut,delta(n_dis),sig_UtR)

              F(n_dis,j)   = sig_UtR(1,1) + (-Kappa/A0)*Hd_(j)
              F(n_dis*2,j) = sig_UtR(1,2)! + sforz_tR(1,1)*attr(i-n_dis)
      
      EndDo

      Do j = n_dis+1,n_dis*2
      
              Call DISL2D_STRSS(0.D0,1.D3,x0z0(j-n_dis,1),x0z0(j-n_dis,2),delta(j-n_dis), &
     &                             Hd_(j-n_dis),x0z0(n_dis,1),x0z0(n_dis,2),          &
     &                          sig_Ub(1,1),sig_Ub(1,2),sig_Ub(2,2))
     
              sig_Ub(2,1) = sig_Ub(1,2)

              Call stress_rot(sig_Ub,delta(n_dis),sig_UbR)

              F(n_dis,j)   = sig_UbR(1,1)
              F(n_dis*2,j) = sig_UbR(1,2)! + sforz_bR(1,1)*attr(i-n_dis)
      
      EndDo

      Do i=1,n_dis

              Call DISL2D_STRSS(1.D3,0.D0,x0z0(n_dis,1),x0z0(n_dis,2),delta(n_dis), &
     &                             Hd_(n_dis),x0z0(i,1),x0z0(i,2),          &
     &                          sig_Ut(1,1),sig_Ut(1,2),sig_Ut(2,2))
     
              sig_Ut(2,1) = sig_Ut(1,2)

              Call stress_rot(sig_Ut,delta(i),sig_UtR)

              F(i,n_dis)       = sig_UtR(1,1) + (-Kappa/A0)*Hd_(n_dis)
              F(i+n_dis,n_dis) = sig_UtR(1,2)! + sforz_tR(1,1)*attr(i-n_dis)

      EndDo

      Do i=1,n_dis

              Call DISL2D_STRSS(0.D0,1.D3,x0z0(n_dis,1),x0z0(n_dis,2),delta(n_dis), &
     &                             Hd_(n_dis),x0z0(i,1),x0z0(i,2),          &
     &                          sig_Ub(1,1),sig_Ub(1,2),sig_Ub(2,2))
     
              sig_Ub(2,1) = sig_Ub(1,2)

              Call stress_rot(sig_Ub,delta(i),sig_UbR)

              F(i,n_dis*2)       = sig_UbR(1,1)
              F(i+n_dis,n_dis*2) = sig_UbR(1,2)! + sforz_bR(1,1)*attr(i-n_dis)

      EndDo      

      DEALLOCATE (F_buf)
      ALLOCATE (F_buf(n_dis*2,n_dis*2))
      
      F_buf(:,:) = F(:,:)
      
      Else ! (B) more than one disl closed or added, should not happen if disl at the tail close one by one
      
      write(*,*) '******* THIS SHOULD NEVER HAPPEN 2 (sub fill_F in BELEMENT MODULE)'
      read(*,*)

      DEALLOCATE (F_buf)
      GoTo 10

      EndIf ! (B)
      EndIf ! (A)
      EndIf
      EndIf
      
!***************** MOVING CRACK *******************

      If (eta.gt.0.D0) Then

      n_dis_old = size(Ut_old)
      
      ALLOCATE (COEFA1(n_dis_old+1,n_dis_old+1),COEFA2(n_dis_old+1))
      ALLOCATE (COEFB1(n_dis_old+1,n_dis_old+1),COEFB2(n_dis_old+1))
      ALLOCATE (COEFC1(n_dis_old+1,n_dis_old+1),COEFC2(n_dis_old+1))
      ALLOCATE (COEFD1(n_dis_old+1,n_dis_old+1),COEFD2(n_dis_old+1))

!~ !      Ut_old(:) = Ut_old(:)*1.D-3
      MULT = 12.D0*v*eta*1.D9 ! keep the opening in [mm] or [m] and multiply by 1.D9 the coefficient to adjust the dimention, numerically more efficient.

      A_old = 0.D0
      Do i = 1,n_dis_old
        A_old = A_old + Ut_old(i)*Hd_old(i)
      EndDo 

!~ !********* COMPUTE NEW COEFFICIENTS *********

! **** A1 ****
      COEFA1(:,:) = 0.D0
      COEFA1(1,1) = MULT*Hd_old(1)/Ut_old(1)**3
      Do i=2,n_dis_old
        Do j=1,i
          COEFA1(i,j) = COEFA1(i-1,j) + MULT*Hd_old(i)/Ut_old(i)**3
        EndDo
      EndDo
      COEFA1(n_dis_old+1,:) = COEFA1(n_dis_old,:)                       ! Exact, no approx (other than using Ut_old(i) instead of Ut_averege(i)

! **** A2 ****
      COEFA2(1) = -MULT*Hd_old(1)/Ut_old(1)**2
      Do i=2,n_dis_old
        SUMA2_j = 0.D0
        Do j=1,i
          SUMA2_j = SUMA2_j + Ut_old(j)/Ut_old(i)**3
        EndDo
        COEFA2(i) = COEFA2(i-1) - MULT*Hd_old(i)*SUMA2_j
      EndDo
      COEFA2(n_dis_old+1) = COEFA2(n_dis_old)                           ! Exact, no approx (other than using Ut_old(i) instead of Ut_averege(i)

! **** B1 ****
      COEFB1(:,:) = 0.D0
      Do i=1,n_dis_old
        Do j=1,i
          COEFB1(i,j) = -0.5D0*MULT*Hd_old(i)/Ut_old(i)**3
        EndDo
      EndDo
      COEFB1(n_dis_old+1,:) = 0.D0                                      ! Exact, no approx (other than using Ut_old(i) instead of Ut_averege(i))

! **** B2 ****
      Do i=1,n_dis_old
        SUMB2_j = 0.D0
        Do j=1,i
          SUMB2_j = SUMB2_j + Hd_old(j)*Ut_old(j)
        EndDo
        COEFB2(i) = 0.5D0*MULT*SUMB2_j/Ut_old(i)**3
      EndDo
      COEFB2(n_dis_old+1) = 0.D0                                        ! Exact, no approx (other than using Ut_old(i) instead of Ut_averege(i))

! **** C1 ****
      COEFC1(:,:) = 0.D0
      Do i=2,n_dis_old+1
        Do j=1,i-1
          COEFC1(i,j) = -0.5D0*MULT*Hd_old(j)/Ut_old(j)**3
        EndDo
      EndDo

! **** C2 ****
      COEFC2(1) = 0.D0
      Do i=2,n_dis_old+1
        SUMC2_j = 0.D0
        Do j=1,i-1
          SUMC2_j = SUMC2_j + Hd_old(j)/Ut_old(j)**2
        EndDo
        COEFC2(i) = 0.5D0*MULT*SUMC2_j                                   
      EndDo

! **** D1 ****
      COEFD1(:,:) = 0.D0
      Do i=1,n_dis_old
        COEFD1(i,i) = -0.25D0*MULT*Hd_old(i)/Ut_old(i)**3
      EndDo
      COEFD1(n_dis_old+1,n_dis_old+1) = COEFD1(n_dis_old,n_dis_old)     ! Here I made the approx of using Ut_old(n_dis_old) when i=n_dis_old+1

! **** D2 ****
      Do i=1,n_dis_old
        COEFD2(i) = 0.25D0*MULT*Hd_old(i)/Ut_old(i)**2
      EndDo
      COEFD2(n_dis_old+1) = 0.D0                                        ! Exact, no approx (other than using Ut_old(i) instead of Ut_averege(i))

!~ !**************** UPDATE F ******************

      Dn = n_dis_old - (n_dis-1) ! if Dn = 0 -> the crack is growing and no tail elements are closing
                                 ! if Dn > 0 -> the tail element(s) of the crack are closing
      If (Dn.lt.0) Then          ! if Dn < 0 -> something is wrong
        Write (*,*) '**********************' 
        Write (*,*) 'n_dis_old is less than n_dis-1' 
        Write (*,*) 'message from SUB fill_MOV_F' 
        Write (*,*) 'THE RUN WILL STOP HERE'
        Write (*,*) '**********************' 
        STOP
      EndIf

      Do i=1,n_dis
        Do j=1,i
          F(i,j) = F(i,j) + COEFA1(i+Dn,j+Dn) + COEFB1(i+Dn,j+Dn) + COEFC1(i+Dn,j+Dn) + COEFD1(i+Dn,j+Dn)
        EndDo
      EndDo

      Do i=1,n_dis
        Do j=1,n_dis
          F(i,j) = F(i,j) + (Hd_(i)/A_old)*(COEFA2(i+Dn) + COEFB2(i+Dn) + COEFC2(i+Dn) + COEFD2(i+Dn))
        EndDo
      EndDo

!~ !      Ut_old(:) = Ut_old(:)*1.D3

      DEALLOCATE (COEFA1,COEFA2,COEFB1,COEFB2,COEFC1,COEFC2,COEFD1,COEFD2)
      
      EndIf

!*************** END MOVING CRACK *****************
          
      End subroutine fill_DYN_F

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine STRESS_CONSTR(sig,tau,x0z0,Hd_,delta,n_dis)
!~ set shear stress and confining pressure acting on each element of the crack.
!~ will call subs that compute lithostatic, tectonic and topographic stress
!~ sig   OUTPUT
!~ tau   OUTPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ n_dis INPUT
      IMPLICIT NONE
      
      REAL(8) :: sig(:),tau(:),x0z0(:,:),Hd_(:),delta(:)
      INTEGER :: n_dis

      REAL(8), ALLOCATABLE :: P_lit(:)     
      INTEGER :: i
!~       REAL(8) :: S_tk(2,2),S_tkR(2,2),S_ex(2,2),S_exR(2,2)
      REAL(8) :: sig_N,sig_T,z_tip                                       ! new

!~       S_tkR = 0.D0
!~       S_exR = 0.D0

      ALLOCATE (P_lit(n_dis))

      Call LITH_P(P_lit,x0z0,n_dis)
      
      Do i=1,n_dis
      
!~       Call EX_STRESS(S_ex,x0z0(i,1),x0z0(i,2))
!~       Call stress_rot(S_ex,delta(i),S_exR)

!~       If (Tk_stress_on.eqv..true.) Then
!~         Call TK_STRESS(S_tk,x0z0(i,1),x0z0(i,2))
!~         Call stress_rot(S_tk,delta(i),S_tkR)
!~       EndIf

!~       sig(i) = S_tkR(1,1) + S_exR(1,1) - P_lit(i)
!~       tau(i) = S_tkR(1,2) + S_exR(1,2)
	  z_tip = x0z0(n_dis,2) - Hd_(n_dis)*cos(delta(n_dis))*0.5d0
      Call COMP_sigN_sigT(sig_N,sig_T,x0z0(i,1),x0z0(i,2),z_tip,delta(i))     ! new

      sig(i) = sig_N - P_lit(i)                                         ! new
      tau(i) = sig_T                                                    ! new
      
      EndDo

      DEALLOCATE (P_lit)
      
      End Subroutine STRESS_CONSTR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_sigN_sigT_old(sig_N,sig_T,x,z,delta)
!~ compute normal and shear stress in x,z with respect to a plane dipping as delta, due to the external stress and tectonic stress
!~ will not account for the lithostatic pressure
!~ sig_N OUTPUT
!~ sig_T OUTPUT
!~ x,z   INPUT
!~ delta INPUT
      IMPLICIT NONE
      
      REAL(8) :: sig_N,sig_T,x,z,delta
      INTEGER :: n_dis

      REAL(8) :: S_tk(2,2),S_tkR(2,2),S_ex(2,2),S_exR(2,2)

      S_tk  = 0.D0
      S_ex  = 0.D0
      S_tkR = 0.D0
      S_exR = 0.D0

      Call EX_STRESS(S_ex,x,z)
      Call stress_rot(S_ex,delta,S_exR)

      If (Tk_stress_on.eqv..true.) Then
        Call TK_STRESS(S_tk,x,z)
        Call stress_rot(S_tk,delta,S_tkR)
      EndIf

      sig_N = S_tkR(1,1) + S_exR(1,1)
      sig_T = S_tkR(1,2) + S_exR(1,2)
      
      End Subroutine COMP_sigN_sigT_old
      !~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_sigN_sigT(sig_N,sig_T,x,z,ztop,delta)
!~ compute normal and shear stress in x,z with respect to a plane dipping as delta, due to the external stress and tectonic stress
!~ will not account for the lithostatic pressure
!~ sig_N OUTPUT
!~ sig_T OUTPUT
!~ x,z   INPUT
!~ delta INPUT
      IMPLICIT NONE
      
      REAL(8) :: sig_N,sig_T,x,z,ztop,delta
      INTEGER :: n_dis

      REAL(8) :: S_tk(2,2),S_tkR(2,2),S_ex(2,2),S_exR(2,2)

      S_tk  = 0.D0
      S_ex  = 0.D0
      S_tkR = 0.D0
      S_exR = 0.D0
!~ 	  write(*,*)ztop, zl
	  if (ztop.LE.zl) then
          Call EX_STRESS(S_ex,x,z)
          Call stress_rot(S_ex,delta,S_exR)
      endif
      
      If (Tk_stress_on.eqv..true.) Then
        Call TK_STRESS(S_tk,x,z)
        Call stress_rot(S_tk,delta,S_tkR)
      EndIf
	  
      sig_N = S_tkR(1,1) + S_exR(1,1)
      sig_T = S_tkR(1,2) + S_exR(1,2)

      
      End Subroutine COMP_sigN_sigT

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine LITH_P(P_lit,x0z0,n_dis)
!~ set lithostatic pressure acting on each element of the b. el. crack
      IMPLICIT NONE
      
      REAL(8) :: P_lit(:),x0z0(:,:)
      INTEGER :: n_dis

!~       REAL(8) :: z
      REAL(8) :: P_lit_const                                            ! new
      INTEGER :: k

      If (Lt_stress_on.eqv..true.) Then

        Do k=1,n_dis

          Call COMP_LITH_P(P_lit(k),x0z0(k,2))                          ! new
!~           z  = x0z0(k,2)

!~           If (z.le.z_dns) Then

!~             If (Kr2.eq.0.D0) Then
!~               write (*,*) 'DIKES DONT FLY :('
!~               stop
!~             EndIf

!~             P_lit(k) = Kr2*( exp(Plit2/Kr2*(z-z_ref)) - 1.D0 )

!~           Else ! (z.gt.dnsty_tr)

!~             P_lit(k) = Kr1*( exp(Plit1/Kr1*(z-z_ref)) - 1.D0 ) +     &
!~      &                 (Plit2_dns - Plit1_dns) +                     &
!~      &                 (Plit2_dns - Plit1_dns)*(z - z_dns)*Plit1/Kr1

!~           EndIf

        EndDo

      Else
      
        Call COMP_LITH_P(P_lit_const,0.D0)                              ! new
        
        P_lit = P_lit_const                                             ! new
!~         z  = z_ref

!~         P_lit = Kr1*( exp(Plit1/Kr1*z) - 1.D0 )

      EndIf

      End Subroutine LITH_P

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_LITH_P(Plit,z)
!~ compute lithostatic pressure at depth z
      IMPLICIT NONE
      
      REAL(8) :: Plit,z

      If (Lt_stress_on.eqv..true.) Then

          If (z.le.z_dns) Then

            If (Kr2.eq.0.D0) Then
              write (*,*) '********************************************************'
              write (*,*) '|In subroutine COMP_LITH_P (BELEMENT_4FE-STRESS.f90)   |'
              write (*,*) '|the depth (z) for lithostatic pressure calculation    |'
              write (*,*) '|is shallower than z_dns (density-transition depth)    |'
              write (*,*) '|but the density of the shallower layer is set to ZERO.|'
              write (*,*) '*************  The simulation STOPPED  *****************'
              stop
            EndIf

            Plit = Kr2*( exp(Plit2/Kr2*(z-z_ref)) - 1.D0 )

          Else ! (z.gt.dnsty_tr)

            Plit = Kr1*( exp(Plit1/Kr1*(z-z_ref)) - 1.D0 ) +     &
     &                 (Plit2_dns - Plit1_dns) +                     &
     &                 (Plit2_dns - Plit1_dns)*(z - z_dns)*Plit1/Kr1

          EndIf

      Else
      
        Plit = Kr1*( exp(Plit1/Kr1*z_ref) - 1.D0 )

      EndIf

      End Subroutine COMP_LITH_P

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine TK_STRESS(S_tk,x,z)
!~ compute tectonic stress tensor acting in (x,z)
      IMPLICIT NONE
      
      REAL(8) :: S_tk(2,2),x,z

!~         If (z.lt.0.D0) Then
        S_tk(1,1) = SigXX + mx_SXX*x + mz_SXX*z
        S_tk(1,2) = 0.D0
        S_tk(2,1) = 0.D0
        S_tk(2,2) = SigZZ + mx_SZZ*x + mz_SZZ*z
!~         Else
!~         S_tk(1,1) = 0.D0
!~         S_tk(1,2) = 0.D0
!~         S_tk(2,1) = 0.D0
!~         S_tk(2,2) = 0.D0
!~         EndIf

      End Subroutine TK_STRESS

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine stress_rot(stress,delta,stressR)

      IMPLICIT NONE
      REAL(8) :: stress(:,:),stressR(:,:),delta
      REAL(8) :: rot(2,2),rotT(2,2),stress2(2,2)

      rot (1,1) =  cos(delta)
      rot (1,2) =  sin(delta)
      rot (2,1) = -sin(delta)
      rot (2,2) =  cos(delta)

      rotT(1,1) =  cos(delta)
      rotT(1,2) = -sin(delta)
      rotT(2,1) =  sin(delta)
      rotT(2,2) =  cos(delta)

      stress2=matmul(rotT,stress)
      stressR=matmul(stress2,rot)

      End Subroutine stress_rot
      
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine gradient(DeE,Diss,dJ)
!~ Estimation of the gradient of the cost function J=(DeE-Diss)^2
!~ DeE 		input
!~ Diss		input
!~ dJ 		output

      IMPLICIT NONE
      REAL(8) :: DeE,Diss,J,dJ,Ed,vel,velocity
      
      J=(DeE-Diss)**2 		! cost function J = energy release - viscous dissipation
      Ed=Diss/v
      dJ=-2*Ed*(DeE-Diss)	! gradient of J
                  
      End Subroutine gradient
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine update_v(DeE,Diss,vel)
!~ Update the velocity according to the gradient of the cost function J
!~ v = v - alpha dJ/dv
!~ DeE 		input
!~ Diss		input
!~ DJ 		input
!~ vel		output

      IMPLICIT NONE
      REAL(8) :: DeE,Diss,vel,Ed

	  Ed=Diss/v 			! viscous dissipation without velocity
      v=DeE/Ed
      vel=v

      End Subroutine update_v
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine update_v_cst(vel)
!~ Send vel to PROPAGATION module
!~ vel		output

      IMPLICIT NONE
      REAL(8) :: vel
      
      vel=v

      End Subroutine update_v_cst

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      End Module BELEMENT
