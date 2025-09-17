      Module PROPAGATION
      Use COMP_FIELD
      Use BELEMENT

      INTEGER, PRIVATE :: n_iter
      INTEGER, PRIVATE :: n_dir_in
      REAL(8), PRIVATE :: alfa
      REAL(8), PRIVATE :: x_stop_min,x_stop_max,z_stop,z_load,Dz_snap
      REAL(8), PRIVATE :: E_threshold_1,E_threshold_2,E_threshold_int
      LOGICAL, PRIVATE :: int_dir,check_int,vel_var,z_ana
      
      REAL(8), ALLOCATABLE, PRIVATE :: x0z0_old(:,:),delta_old(:),Hd_old(:)
! SEV-ANALOG-TRAJECTORY
      REAL(8), ALLOCATABLE, PRIVATE :: z0delta_ana(:,:)
      INTEGER, PRIVATE :: n_delta
! END SEV-ANALOG-TRAJECTORY
      REAL(8), ALLOCATABLE, PRIVATE :: Ut_old(:),Ub_old(:)
      REAL(8), PRIVATE :: E_old
      INTEGER, PRIVATE :: n_dis_old

      REAL(8), PRIVATE :: pi=4.D0*datan(1.D0)
      
      contains

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_PROP_PAR(n_iter_,n_dir_,alfa_,x_stop_min_,x_stop_max_,z_stop_, z_load_,  &
     &                        Dz_snap_,E_threshold_1_,E_threshold_2_,E_threshold_int_,vel_var_,z_ana_)

      IMPLICIT NONE
      
      INTEGER :: n_iter_
      INTEGER :: n_dir_
      REAL(8) ::  alfa_
      REAL(8) :: x_stop_min_,x_stop_max_,z_stop_,z_load_,Dz_snap_
      REAL(8) :: E_threshold_1_,E_threshold_2_,E_threshold_int_
      LOGICAL :: vel_var_,z_ana_

      n_iter  = n_iter_

      n_dir_in = n_dir_
      alfa  = alfa_

      x_stop_min  = x_stop_min_
      x_stop_max  = x_stop_max_
      z_stop      = z_stop_
! SEV - STRESS FIXED DEPTH
	  z_load      = z_load_
! END - SEV - STRESS FIXED DEPTH
      Dz_snap     = Dz_snap_

      E_threshold_1   = E_threshold_1_
      E_threshold_2   = E_threshold_2_
      E_threshold_int = E_threshold_int_

      int_dir = .false.
      If ((E_threshold_1.eq.E_threshold_2).and.(E_threshold_1.eq.E_threshold_int)) Then
        check_int = .false.
      Else
        check_int = .true.
      EndIf
      
      vel_var = vel_var_
      z_ana = z_ana_

      End Subroutine SET_PROP_PAR

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine START_PROP(x0z0_in,delta_in,Hd_in,A0,n_dis,run,xf_top,zf_top)
!~ x0z0_in  INPUT
!~ delta_in INPUT
!~ Hd_in    INPUT
!~ A0       INPUT
!~ n_dis    INPUT/OUTPUT
!~ run      INPUT
!~ xf_top  OUTPUT
!~ zf_top  OUTPUT
      IMPLICIT NONE

      REAL(8) :: x0z0_in(:,:),delta_in(:),Hd_in(:)
      REAL(8) :: A0,xf_top,zf_top
      INTEGER :: n_dis,run

      REAL(8), ALLOCATABLE :: Ut(:),    Ub(:),    DPvisc(:), x0z0(:,:), delta(:), Hd_(:)
      REAL(8) :: A,A_old,E,D,DeE,Diss,OverP_av,OverP_av2,Ub_max,Ub_min
      REAL(8) :: x_top,z_top,x_start,z_start,z_snap
!~ DELTA -------------------------------------------------------------
      REAL(8) :: L_norm,dip_av,dip_ndis,dip_max,Sig_tip(2,2),           &
     &           Sig_xz_av,Sig_xz_n,Sig_xz_max,DELTA1,DELTA2,DELTA3
!~ -------------------------------------------------------------------
!~ Plit sig_T sig_N ------------------------------------------------
      REAL(8) :: Plit,sig_N,sig_T
!~ -------------------------------------------------------------------
!~ Cost function -----------------------------------------------------
      REAL(8) :: dJ,vel
!~ -------------------------------------------------------------------
      INTEGER :: n_dir,iter,i,iv
      LOGICAL :: esc,esc_old
      CHARACTER(50) :: outf01,outf02,outf03,outf04,outf05,outf06,outf07,outf08
      CHARACTER(2)  :: n_outf
      CHARACTER(3)  :: iter_outf

      REAL(8), ALLOCATABLE :: P_lit(:),OverP(:)
      
      ! REAL(8) :: vel,vit

      Write(n_outf,20) run
      outf01 = './output/XtZtDeE_'//n_outf//'.dat'
      outf02 = './output/RUNINFO_'//n_outf//'.dat'
!~ DELTA -------------------------------------------------------------
      outf03 = './output/XtZtDELTAINF_'//n_outf//'.txt'
      outf04 = './output/XtZtDELTA_'//n_outf//'.dat'
!~ -------------------------------------------------------------------
      outf05 = './output/XtZtPlitSigNSigT_'//n_outf//'.dat'
      outf06 = './output/XtZtvel_'//n_outf//'.dat'
      

      Open(21,file = outf01,access='sequential',form='formatted',status='new')
      Open(22,file = outf02,access='sequential',form='formatted',status='new')
!~ DELTA -------------------------------------------------------------
      Open(23,file = outf03,access='sequential',form='formatted',status='new')
      Open(24,file = outf04,access='sequential',form='formatted',status='new')
!~ -------------------------------------------------------------------
      Open(25,file = outf05,access='sequential',form='formatted',status='new')
      Open(26,file = outf06,access='sequential',form='formatted',status='new')

      E = 0.D0
      A = A0

      ALLOCATE (Ut(n_dis),Ub(n_dis),DPvisc(n_dis),x0z0(n_dis,2),delta(n_dis),Hd_(n_dis))

      x0z0  = x0z0_in
      delta = delta_in
      Hd_   = Hd_in

      Write (22,*) 'n_dis input  =', n_dis
      Write (*,*)  'n_dis input  =', n_dis
      Write (22,*) 'A0  =', A
      Write (*,*)  'A0  =', A

!~       Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,esc)
      Call SOLVE_STA_CRACK(Ut,Ub,x0z0,z_load,delta,Hd_,DPvisc,A,n_dis,esc)
      
!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE (Ut_old(n_dis),Hd_old(n_dis))
      Ut_old(:) = 0.D0
      Hd_old(:) = 0.D0
!!!!!!!!!!!!!!!!!!!!!!!!  
!~       Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
      Call COMP_DYN_ENERGY(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D)
!!!!!!!!!!!!!!!!!!!!!!!!
      DEALLOCATE(Ut_old, Hd_old)
!!!!!!!!!!!!!!!!!!!!!!!!
      Write (22,*) 'Etot =',E, 'D=', D
            
      x_top = x0z0(n_dis,1) - Hd_(n_dis)*sin(delta(n_dis))*0.5d0
      z_top = x0z0(n_dis,2) - Hd_(n_dis)*cos(delta(n_dis))*0.5d0
      
      x_start = x_top
      z_start = z_top

      z_snap = z_top - Dz_snap

      Do i=1,n_dis
        Write (21,21) x0z0(i,1),x0z0(i,2),0.D0
      EndDo
      Write (22,*) 'n_dis output =', n_dis
      Write (*,*)  'n_dis output =', n_dis
      Write (22,*) 'A  =', A
      Write (*,*)  'A  =', A
      Write (22,*) 'delta_in  =', delta(n_dis)/pi*180.D0
      Write (*,*)  'delta_in  =', delta(n_dis)/pi*180.D0
      Write (22,*) '(x_start;z_start)=(',x_start,';',z_start,')'
      Write (*,*)  '(x_start;z_start)=(',x_start,';',z_start,')'
      
      Call DISPL_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,0)
      Call STRSS_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,0)
      Call CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd_,n_dis,run,0)

!~ CALCOLO OVERPRESSURE ----------------------------------------------
      ALLOCATE (OverP(n_dis))
      Call COMP_DYN_OverP(Ut,Ub,x0z0,DPvisc,delta,Hd_,A,n_dis,OverP,OverP_av)
      Write (22,*) 'OverP_av  =', OverP_av
      Write (*,*)  'OverP_av  =', OverP_av
      Write (22,*) 'OverP_max =', OverP(n_dis)
      Write (*,*)  'OverP_max =', OverP(n_dis)
      DEALLOCATE (OverP)
!~ FINE CALCOLO OVERPRESSURE -----------------------------------------

!~ DELTA -------------------------------------------------------------
        Call COMP_DELTA_iter(x_top,z_top,delta,Hd_,n_dis,z_start,       &
     &                           OverP_av,L_norm,                       &
     &                           dip_av,dip_ndis,dip_max,               &
     &                           Sig_tip,Sig_xz_av,Sig_xz_n,Sig_xz_max, &
     &                           DELTA1,DELTA2,DELTA3)
        Write (24,22) x_top,z_top,OverP_av,L_norm,                      &
       &              Sig_xz_av,Sig_xz_n,Sig_xz_max,                    &
       &              DELTA1,DELTA2,DELTA3
        Write (23,*) '--------------------------------'
        Write (23,*) '---------- ITER = 0 -------'
        Write (23,*) '--------------------------------'
        Write (23,*) 'x_top, z_top    = ',x_top,z_top
        Write (23,*) 'S_xx, S_zz, Sxz = ',Sig_tip(1,1),Sig_tip(2,2),Sig_tip(1,2)
        Write (23,*) 'dip_av          = ',dip_av
        Write (23,*) 'dip_ndis        = ',dip_ndis
        Write (23,*) 'dip_max         = ',dip_max
        Write (23,*) 'S_xz_av         = ',Sig_xz_av
        Write (23,*) 'S_xz_n          = ',Sig_xz_n
        Write (23,*) 'S_xz_max        = ',Sig_xz_max
        Write (23,*) 'OverP_av        = ',OverP_av
        Write (23,*) 'L_norm          = ',L_norm
        Write (23,*) 'DELTA1, DELTA2, DELTA3  = ',DELTA1,DELTA2,DELTA3
        Write (23,*) ''
!~ END DELTA ---------------------------------------------------------

!~ COMPUTE P_lit sig_N sig_T -----------------------------------------
        Call COMP_LITH_P (Plit,x0z0(n_dis,2))
        Call COMP_sigN_sigT (sig_N,sig_T,x0z0(n_dis,1),x0z0(n_dis,2),x0z0(n_dis,2),delta(n_dis)) 
        Write (25,21) x0z0(n_dis,1),x0z0(n_dis,2),Plit,sig_N,sig_T
!~ END COMPUTE P_lit sig_N sig_T -------------------------------------
		Write (26,21) x0z0(n_dis,1),x0z0(n_dis,2)

      Do iter = 1,n_iter
      
		iv=1

        Write (22,*) '**************************************************'
        Write (*,*)  '**************************************************'
        Write (22,*) 'ITER =',iter,'/',n_iter
        Write (*,*)  'ITER =',iter,'/',n_iter

!~         If (iter.le.1) Then
!~           n_dir = 61
!~         Else
!~           If (iter.le.10) Then
!~             n_dir = 21
!~           Else
            n_dir = n_dir_in
!~           EndIf
!~         EndIf

        ALLOCATE (Ut_old(n_dis),Ub_old(n_dis),x0z0_old(n_dis,2),delta_old(n_dis),Hd_old(n_dis))
      
        Ut_old(1:n_dis)     = Ut(1:n_dis)
        Ub_old(1:n_dis)     = Ub(1:n_dis)
        x0z0_old(1:n_dis,:) = x0z0(1:n_dis,:)
        delta_old(1:n_dis)  = delta(1:n_dis)
        Hd_old(1:n_dis)     = Hd_(1:n_dis)
        
        A_old     = A
        E_old     = E
        n_dis_old = n_dis
        

11      DEALLOCATE (Ut,Ub,x0z0,delta,Hd_,DPvisc)
        ALLOCATE (Ut(n_dis+1),Ub(n_dis+1),x0z0(n_dis+1,2),delta(n_dis+1),Hd_(n_dis+1),DPvisc(n_dis+1))
 
        x0z0(1:n_dis,:) = x0z0_old(1:n_dis,:)
        delta(1:n_dis)  = delta_old(1:n_dis)
        Hd_(1:n_dis)    = Hd_old(1:n_dis)

!~         If (check_int.eqv..false.) Then
		If (z_ana.eqv..false.) then
			Call TEST_DIR(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D,n_dir,alfa,esc)
		Else
			Call TEST_DIR_ANALOG(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D,n_dir,alfa,esc)
		EndIf

!~         Call TEST_DIR(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D,n_dir,alfa,esc)


!~         Else
!~           If (sign(1.D0,(z_top-Hd_(n_dis)*sin(delta(n_dis))))/sign(1.D0,z_top).lt.0.D0) Then
!~             Call TEST_DIR_INT(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,esc)
!~           Else
!~             Call TEST_DIR(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,esc)
!~           EndIf
!~         EndIf

        x_top = x0z0(n_dis,1) - Hd_(n_dis)*sin(delta(n_dis))*0.5d0
        z_top = x0z0(n_dis,2) - Hd_(n_dis)*cos(delta(n_dis))*0.5d0

!~         If ((delta(n_dis).eq.pi*0.D0).and.(int_dir.eqv..true.)) Then
!~           DeE = (E_old - E)/Hd_(n_dis) - E_threshold_int ! MPa*m ( E = U+W [MPa*km^2*10^3] )
!~         Else
!~           If (z_top.gt.0.D0) Then
!~             DeE = (E_old - E)/Hd_(n_dis) - E_threshold_1 ! MPa*m ( E = U+W [MPa*km^2*10^3] )
!~           Else
!~             DeE = (E_old - E)/Hd_(n_dis) - E_threshold_2 ! MPa*m ( E = U+W [MPa*km^2*10^3] )
!~           EndIf
!~         EndIf

			DeE = (E_old - E)/Hd_(n_dis) - E_threshold_1 ! MPa*m ( E = U+W [MPa*km^2*10^3] )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Diss = -D*1e3/Hd_(n_dis)
		
		If (vel_var.eqv..true.) Then
		
			If ((DeE.le.0.D0).and.(sqrt((z_top-z_load)*(z_top-z_load)).lt.Hd_(n_dis))) Then
				Call update_v_cst(vel)
!~ 				Write(113,*)iter,x_top,z_top,E_old,E,DeE,Diss,vel, iv, 'loop 1'
			Else
				Call gradient(DeE,Diss,dJ)
				If (sqrt(dJ*dJ) .GT. 1.0D-3) Then
					n_dis=n_dis_old
					Call update_v(DeE,Diss,vel)
					Call DEALLOC_Fold_Fbuf()
!~ 					Write(113,*)iter,x_top,z_top,E_old,E,DeE,Diss,vel,iv, 'loop 2'
					iv=iv+1
					GoTo 11
				Else 
					Call update_v_cst(vel)
!~ 					Write(113,*)iter,x_top,z_top,E_old,E,DeE,Diss,vel, iv, 'loop 3'
				EndIf  
			EndIf
		Else 
			Call update_v_cst(vel)
!~ 			Write(113,*)iter,x_top,z_top,E_old,E,DeE,Diss,vel, iv, 'loop 3'
		EndIf
               
        Write(iter_outf,30) iter
        outf07 = './output/Ut_'//iter_outf//'_'//n_outf//'.dat'
        Open(27,file = outf07,access='sequential',form='formatted',status='new')
        Do i=1,n_dis
          Write(27,*) i, Ut(i)
        EndDo
        Close(27)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        Write (22,*) '--------------------------------------------------'
        Write (*,*)  '--------------------------------------------------'
        Write (22,*) 'n_dis =',n_dis
        Write (22,*) 'DeE  (iter=',iter, ') =',DeE
        Write (22,*) 'delta(     ',n_dis,') =',delta(n_dis)*180.D0/pi

        ALLOCATE (P_lit(n_dis))
        Call LITH_P(P_lit,x0z0,n_dis)
        Write (22,*) 'Plit (     ',n_dis,') =',P_lit(n_dis)
        Write (22,*) 'Plit (     ',n_dis-1,') =',P_lit(n_dis-1)
        Write (22,*) '[...]'
        Write (22,*) 'Plit (     ',2,') =',P_lit(2)
        Write (22,*) 'Plit (     ',1,') =',P_lit(1)
        DEALLOCATE (P_lit)

        Write (22,*) 'x0z0(1,1) = ',x0z0(1,1)
        Write (22,*) 'x0z0(1,2) = ',x0z0(1,2)

        Write (22,*) 'A_iter =',A
        Write (22,*) 'x_top  =',x_top
        Write (22,*) 'z_top  =',z_top
        Write (*,*)  'n_dis=',n_dis
        Write (*,*)  'delta(',n_dis,') =',delta(n_dis)*180.D0/pi
        Write (*,*)  'A(iter=',iter,') =',A
        Write (*,*)  'x_top =',x_top,'  z_top =',z_top

        ALLOCATE (OverP(n_dis))
        Call COMP_DYN_OverP(Ut,Ub,x0z0,DPvisc,delta,Hd_,A,n_dis,OverP,OverP_av)
!---SEV overpressure profile
        Write(iter_outf,31) iter
        outf08 = './output/Overpressure_'//iter_outf//'_'//n_outf//'.dat'
        Open(28,file = outf08,access='sequential',form='formatted',status='new')
        Do i=1,n_dis
          write(28,*) i, OverP(i), DPvisc(i)
        EndDo
        Close(28) 
!---END SEV overpressure profile
        Write (22,*) 'OverP_av =', OverP_av
        Write (*,*)  'OverP_av =', OverP_av
        DEALLOCATE (OverP)

!~ DELTA -------------------------------------------------------------
        Call COMP_DELTA_iter(x_top,z_top,delta,Hd_,n_dis,z_start,       &
     &                           OverP_av,L_norm,                       &
     &                           dip_av,dip_ndis,dip_max,               &
     &                           Sig_tip,Sig_xz_av,Sig_xz_n,Sig_xz_max, &
     &                           DELTA1,DELTA2,DELTA3)
!~ -------------------------------------------------------------------

        Write (22,*) 'Ut_max =', maxval(Ut), 'Ut_last_el =', Ut(1)
        Write (*,*) 'Ut_max =', maxval(Ut), 'Ut_last_el =', Ut(1)
        Ub_max = maxval(Ub)
        Ub_min = minval(Ub)
        If (abs(Ub_max).lt.abs(Ub_min)) Ub_max = Ub_min
        Write (22,*) 'Ub_max =', Ub_max, 'Ub_last_el =', Ub(1)
        Write (*,*) 'Ub_max =', Ub_max, 'Ub_last_el =', Ub(1)
		Write (22,*) 'vel =', vel
		Write (*,*)  'vel =', vel

        If ( (esc.eqv..true.).or.(DeE.le.0.D0)) Then
			If (sqrt((z_top-z_load)*(z_top-z_load)).lt.Hd_(n_dis)) Then
				GoTo 12
			else
				If (esc.eqv..true.) write(*,*) '***** TEST DISLOCATIONS CLOSED *****'
				If (DeE.le.0.D0)    write(*,*) '** ENERGY RELEASE UNDER THRESHOLD **'
            
				Call DISPL_FIELD(Ut_old,Ub_old,x0z0_old,delta_old,Hd_old,n_dis_old,run,iter)
				Call STRSS_FIELD(Ut_old,Ub_old,x0z0_old,delta_old,Hd_old,n_dis_old,run,iter)
				Call CRACK_SHAPE(Ut_old,Ub_old,x0z0_old,delta_old,Hd_old,n_dis_old,run,iter)

				DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
				xf_top = x_top
				zf_top = z_top
				Call DEALLOC_Fold_Fbuf()
				EXIT
			endif
        EndIf
		
12      Write (21,21) x_top,z_top,DeE,D/Hd_(n_dis) !, vit !!!!!!!!!!!!!!!!!!
!~ DELTA -------------------------------------------------------------
        Write (24,22) x_top,z_top,OverP_av,L_norm,                      &
       &              Sig_xz_av,Sig_xz_n,Sig_xz_max,                    &
       &              DELTA1,DELTA2,DELTA3
        
        Write (23,*) '--------------------------------'
        Write (23,*) '------- ITER = ',iter,' --------'
        Write (23,*) '--------------------------------'
        Write (23,*) 'x_top, z_top    = ',x_top,z_top
        Write (23,*) 'S_xx, S_zz, Sxz = ',Sig_tip(1,1),Sig_tip(2,2),Sig_tip(1,2)
        Write (23,*) 'dip_av          = ',dip_av
        Write (23,*) 'dip_ndis        = ',dip_ndis
        Write (23,*) 'dip_max         = ',dip_max
        Write (23,*) 'S_xz_av         = ',Sig_xz_av
        Write (23,*) 'S_xz_n          = ',Sig_xz_n
        Write (23,*) 'S_xz_max        = ',Sig_xz_max
        Write (23,*) 'OverP_av        = ',OverP_av
        Write (23,*) 'L_norm          = ',L_norm
        Write (23,*) 'DELTA1, DELTA2, DELTA3  = ',DELTA1,DELTA2,DELTA3
        Write (23,*) ''
!~ END DELTA ---------------------------------------------------------
!~ COMPUTE P_lit sig_N sig_T -----------------------------------------
        Call COMP_LITH_P (Plit,x0z0(n_dis,2))
        Call COMP_sigN_sigT (sig_N,sig_T,x0z0(n_dis,1),x0z0(n_dis,2),x0z0(n_dis,2),delta(n_dis)) 
        Write (25,21) x0z0(n_dis,1),x0z0(n_dis,2),Plit,sig_N,sig_T
!~ END COMPUTE P_lit sig_N sig_T -------------------------------------
		Write (26,*) x0z0(n_dis,1),x0z0(n_dis,2),vel,maxval(Ut)

        If (movie.eqv..true.) Call CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)
      
        If ( ( z_top.le.z_stop     + Hd_(n_dis)*cos(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( x_top.le.x_stop_min + Hd_(n_dis)*sin(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( x_top.ge.x_stop_max - Hd_(n_dis)*sin(delta(n_dis)) ) &
     &              .or.                                            &
     &       ( z_top.le.z_snap     + Hd_(n_dis)*cos(delta(n_dis)) ) &
     &              .or.                                            &
     &        ( iter.eq.n_iter )                                    ) Then

          Call DISPL_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)
          Call STRSS_FIELD(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)
          If (movie.eqv..false.) Call CRACK_SHAPE(Ut,Ub,x0z0,delta,Hd_,n_dis,run,iter)

          If (z_top.le.z_snap+Hd_(n_dis)*cos(delta(n_dis))) Then
            z_snap=z_snap-Dz_snap
            DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
            CYCLE
          EndIf

          DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)
          xf_top = x_top
          zf_top = z_top
          Call DEALLOC_Fold_Fbuf()
          EXIT
        EndIf

        DEALLOCATE (Ut_old,Ub_old,x0z0_old,delta_old,Hd_old)

      EndDo

      DEALLOCATE (Ut,Ub,x0z0,delta,Hd_,DPvisc)
      xf_top = x_top
      zf_top = z_top
      Call DEALLOC_Fold_Fbuf()

20    Format (i2.2)
21    Format (2x,5f20.10)
22    Format (2x,10d16.6)
30    Format (i3.3)
31    Format (i3.3,2f16.6)

      Close(21)
      Close(22)
      Close(23)
      Close(24)
      Close(25)

      End Subroutine START_PROP

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine TEST_DIR(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D,n_dir,alfa,esc)
!~ !test different directions of propagation and give as output the crack
!~ !opening in the direction that minimise the total energie of the system
!~ !Ut     INPUT/OUTPUT
!~ !Ub     INPUT/OUTPUT
!~ !x0z0   INPUT/OUTPUT
!~ !delta  INPUT/OUTPUT
!~ !Hd_    INPUT/OUTPUT
!~ !DPvisc OUTPUT
!~ !A      INPUT/OUTPUT
!~ !n_dis  INPUT/OUTPUT
!~ !E      OUTPUT
!~ !n_dir  INPUT
!~ !alfa   INPUT
!~ !esc    OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),DPvisc(:)
      REAL(8) :: A,E,alfa,D
      INTEGER :: n_dis,n_dir
      LOGICAL :: esc

      REAL(8), ALLOCATABLE :: Ut_dir(:),Ub_dir(:),x0z0_dir(:,:),delta_dir(:),Hd_dir(:),DPvisc_dir(:)
      REAL(8) :: A_dir,E_min,angle
      REAL(8), ALLOCATABLE :: E_(:),E_tot(:),D_(:)
      INTEGER :: n_dis_in,dir,i
      INTEGER :: n_dis_min,delta_n_dis
      INTEGER, ALLOCATABLE :: n_dis_(:)
      LOGICAL, ALLOCATABLE :: esc_dir(:)
      LOGICAL :: redo_test_dir
      
      !REAL(8) :: vel

      redo_test_dir = .false.
      
      Call SET_Fold_eq_Fbuf

10    If (n_dir.eq.1) Then
      
      dir=1
      n_dis = n_dis+1

      delta(n_dis)  = delta(n_dis-1) - ((n_dir+1)*0.5D0-dir)*alfa
      Hd_(n_dis)    = Hd_(n_dis-1)
      x0z0(n_dis,1) = x0z0(n_dis-1,1) -                                  &
     &                          0.5D0*Hd_(n_dis-1)*sin(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *sin(delta(n_dis))
      x0z0(n_dis,2) = x0z0(n_dis-1,2) -                                  &
     &                          0.5D0*Hd_(n_dis-1)*cos(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *cos(delta(n_dis))
        
      esc = .false.

      Call SOLVE_DYN_CRACK(Ut_old,Hd_old,Ut,Ub,DPvisc,x0z0,delta,Hd_,A,n_dis,esc)
!~       Call SOLVE_STA_CRACK(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,esc)
!~       Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,esc)
      Call COMP_DYN_ENERGY(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D)
!~       Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
      
      Else

      Write (22,*) '--------------------------------------------------'
      Write (22,*) '       ** TEST DIRECTIONS OF PROPAGATION **       '
      Write (22,*) '--------------------------------------------------'
      Write (*,*)  '--------------------------------------------------'

      ALLOCATE (E_(n_dir),D_(n_dir),E_tot(n_dir),n_dis_(n_dir),esc_dir(n_dir))

      ALLOCATE (x0z0_dir(n_dis+1,2),delta_dir(n_dis+1),Hd_dir(n_dis+1), &
     &          DPvisc_dir(n_dis+1),Ut_dir(n_dis+1),Ub_dir(n_dis+1))


!~       Call omp_set_num_threads(5) 
!~ !$omp parallel do private(A_dir,n_dis_,x0z0_dir,delta_dir,Hd_dir,esc_dir,Ut_dir,Ub_dir,E_,E_tot,angle)
      Do dir=1,n_dir

        A_dir     = A
        n_dis_(dir) = n_dis+1

        x0z0_dir(1:n_dis,:) = x0z0(1:n_dis,:)
        delta_dir(1:n_dis)  = delta(1:n_dis)
        Hd_dir(1:n_dis)     = Hd_(1:n_dis)

        delta_dir(n_dis+1)  = delta_dir(n_dis) - ((n_dir+1)*0.5D0-dir)*alfa
        Hd_dir(n_dis+1)     = Hd_dir(n_dis)
        x0z0_dir(n_dis+1,1) = x0z0_dir(n_dis,1) -                             &
     &                          0.5D0*Hd_dir(n_dis)  *sin(delta_dir(n_dis)) - &
     &                          0.5D0*Hd_dir(n_dis+1)*sin(delta_dir(n_dis+1))
        x0z0_dir(n_dis+1,2) = x0z0_dir(n_dis,2) -                             &
     &                          0.5D0*Hd_dir(n_dis)  *cos(delta_dir(n_dis)) - &
     &                          0.5D0*Hd_dir(n_dis+1)*cos(delta_dir(n_dis+1))
        
        esc_dir(dir) = .false.

        Call SOLVE_DYN_CRACK(Ut_old,Hd_old,Ut_dir,Ub_dir,DPvisc_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_(dir),esc_dir(dir))
!~         Call SOLVE_DYN_CRACK_OLD(Ut_old,DPvisc_old,Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,DPvisc_dir,A_dir,n_dis_(dir),esc_dir(dir))
!~         Call SOLVE_STA_CRACK(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,DPvisc_dir,A_dir,n_dis_(dir),esc_dir(dir))
!~         Call SOLVE_CRACK(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_(dir),esc_dir(dir))
        If (esc_dir(dir).eqv..false.) Then
!~           Call COMP_ENERGY(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_(dir),E_(dir))
          Call COMP_DYN_ENERGY(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,DPvisc_dir,A_dir,n_dis_(dir),E_(dir),D_(dir))
        Else
          E_(dir) = E_old
        EndIf
        Call SET_Fbuf_eq_Fold

!~         If (int_dir.eqv..false.) Then
!~           If (x0z0_dir(n_dis_(dir),2)-0.5D0*Hd_dir(n_dis_(dir))*cos(delta_dir(n_dis_(dir))).gt.0.D0) Then
        E_tot(dir) = E_(dir) + E_threshold_1*Hd_dir(n_dis_(dir))
!~           Else
!~             E_tot(dir) = E_(dir) + E_threshold_2*Hd_dir(n_dis_(dir))
!~           EndIf
!~         Else
!~           E_tot(dir) = E_(dir) + E_threshold_int*Hd_dir(n_dis_(dir))
!~         EndIf

        angle=delta_dir(n_dis_(dir))*180.D0/pi
        Write(22,*)'Etot(',dir,'/',n_dir,')=',E_(dir),' delta=',angle,' n_dis=',n_dis_(dir),' esc_dir=',esc_dir(dir)
        Write(*,*) 'test direction (',dir,'/',n_dir,') =',delta_dir(n_dis_(dir))*180.D0/pi

      EndDo
!~ !$omp end parallel do

      DEALLOCATE (x0z0_dir,delta_dir,Hd_dir,Ut_dir,Ub_dir,DPvisc_dir)

      E_min = E_tot(1)
      Do i=1,n_dir
        If (E_tot(i).le.E_min) Then
          E     = E_(i)
          E_min = E_tot(i)
          D		= D_(i)
          dir   = i
        EndIf
      EndDo

!~    IF THE NUMBER OF DISL. EL. CHANGES FOR DIFFERENT TEST DIRECTION THE ENERGY RELEASE CAN BE AFFECTED BY THAT,
!~    THEREFORE I RE-TEST THE DIRECTIONS USING THE MIN NUMBER OF DISL. EL. FOUND TESTING THE DIRECTIONS OF PROPAGATION
      If (minval(n_dis_).ne.maxval(n_dis_)) Then
      If ( ( (dir.ne.1).and.(dir.ne.n_dir).and.(n_dis_(dir-1).ne.n_dis_(dir+1)) ) &
     &                                    .or.                                    &
     &                        ( (dir.eq.1).and.(n_dis_(dir).ne.n_dis_(dir+1)) )   &
     &                                    .or.                                    &
     &                    ( (dir.eq.n_dir).and.(n_dis_(dir).ne.n_dis_(dir-1)) ) ) Then
     
        If (redo_test_dir.eqv..true.) Then
          Write (*,*) 'ERROR IN SUBROUTINE TEST_DIR - MODULE PROPAGATION'
          Stop
        EndIf
      
        redo_test_dir = .true.
      
        n_dis_min   = minval(n_dis_)
        delta_n_dis = (n_dis_old+1) - n_dis_min
      
        n_dis = n_dis_min-1
        x0z0(1:n_dis,:) = x0z0_old(delta_n_dis+1:n_dis_old,:)
        delta(1:n_dis)  = delta_old(delta_n_dis+1:n_dis_old)
        Hd_(1:n_dis)    = Hd_old(delta_n_dis+1:n_dis_old)
        
        Call DEALLOC_Fold_Fbuf()
        DEALLOCATE (E_,E_tot,n_dis_,esc_dir)
        Call SOLVE_DYN_CRACK(Ut_old,Hd_old,Ut,Ub,DPvisc,x0z0,delta,Hd_,A,n_dis,esc)
!~         Call SOLVE_DYN_CRACK_OLD(Ut_old,DPvisc_old,Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,esc)
!~         Call SOLVE_STA_CRACK(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,esc)
!~         Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,esc)
        
        GoTo 10
        
      EndIf
      EndIf
      
      If (redo_test_dir.eqv..true.) Then
        n_dis           = n_dis_old
        x0z0(1:n_dis,:) = x0z0_old(1:n_dis,:)
        delta(1:n_dis)  = delta_old(1:n_dis)
        Hd_(1:n_dis)    = Hd_old(1:n_dis)
        Call DEALLOC_Fold_Fbuf()
      EndIf      
      
      n_dis = n_dis+1

      delta(n_dis)  = delta(n_dis-1) - ((n_dir+1)*0.5D0-dir)*alfa
      Hd_(n_dis)    = Hd_(n_dis-1)
      x0z0(n_dis,1) = x0z0(n_dis-1,1) -                               &
     &                          0.5D0*Hd_(n_dis-1)*sin(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *sin(delta(n_dis))
      x0z0(n_dis,2) = x0z0(n_dis-1,2) -                               &
     &                          0.5D0*Hd_(n_dis-1)*cos(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *cos(delta(n_dis))
        
      esc = .false.

      Call SOLVE_DYN_CRACK(Ut_old,Hd_old,Ut,Ub,DPvisc,x0z0,delta,Hd_,A,n_dis,esc)
!~       Call SOLVE_DYN_CRACK_OLD(Ut_old,DPvisc_old,Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,esc)
!~       Call SOLVE_STA_CRACK(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,esc)
!~       Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,esc)
!~       write(115,*)redo_test_dir, 'true',esc, 'false'
      If ((redo_test_dir.eqv..true.).and.(esc.eqv..false.)) Then
!~         Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
        Call COMP_DYN_ENERGY(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D)
      EndIf      
!~ 	  write(115,*)E,D
      DEALLOCATE (E_,E_tot,n_dis_,esc_dir)
      
      EndIf

      End Subroutine TEST_DIR
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine TEST_DIR_ANALOG(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D,n_dir,alfa,esc)
!~ !Impose direction of propagation of analog experiment and give as 
!~ !output the crack opening in the direction of the analog experiment 
!~ !Ut     INPUT/OUTPUT
!~ !Ub     INPUT/OUTPUT
!~ !x0z0   INPUT/OUTPUT
!~ !delta  INPUT/OUTPUT
!~ !Hd_    INPUT/OUTPUT
!~ !A      INPUT/OUTPUT
!~ !n_dis  INPUT/OUTPUT
!~ !E      OUTPUT
!~ !n_dir  INPUT
!~ !alfa   INPUT
!~ !esc    OUTPUT
      IMPLICIT NONE

      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),DPvisc(:)
      REAL(8) :: A,E,alfa, z,D
      INTEGER :: n_dis,n_dir
      LOGICAL :: esc

      INTEGER :: n_dis_in,i
      LOGICAL :: redo_test_dir

      redo_test_dir = .false.
   
      n_dis = n_dis+1
      Hd_(n_dis) = Hd_(n_dis-1)
      z = x0z0(n_dis-1,2) - 0.5*Hd_(n_dis-1)*cos(delta(n_dis-1))
     
      do i=1,n_delta
          if (((z-z0delta_ana(i,1)) .LE. 1e-5) .AND. (z-z0delta_ana(i,1)) .GE. -1e-5 ) then
              delta(n_dis)=z0delta_ana(i,2)
          endif
      enddo
      
      x0z0(n_dis,1) = x0z0(n_dis-1,1) -                                  &
     &                          0.5D0*Hd_(n_dis-1)*sin(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *sin(delta(n_dis))
      x0z0(n_dis,2) = x0z0(n_dis-1,2) -                                  &
     &                          0.5D0*Hd_(n_dis-1)*cos(delta(n_dis-1)) - &
     &                          0.5D0*Hd_(n_dis)  *cos(delta(n_dis))
         
      esc = .false.    
      
      Call SOLVE_DYN_CRACK(Ut_old,Hd_old,Ut,Ub,DPvisc,x0z0,delta,Hd_,A,n_dis,esc)
!~       Call SOLVE_STA_CRACK(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,esc)
!~       Call SOLVE_CRACK(Ut,Ub,x0z0,delta,Hd_,A,n_dis,esc)
      Call COMP_DYN_ENERGY(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D)
!~       Call COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)


      End Subroutine TEST_DIR_ANALOG 
      
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------   

      Subroutine SET_ANALOG_PATH(nb_ana)
      
      CHARACTER*50 :: inpf00
      
      inpf00 = './input/z0delta_analog.txt'
      Open(40,file = inpf00,access='sequential',form='formatted')
      
      ALLOCATE(z0delta_ana(nb_ana,2))
      
      Read(40,*) n_delta
      Do i=1,nb_ana
          Read(40,*) z0delta_ana(i,1), z0delta_ana(i,2)
      enddo
      Close(40)
      End Subroutine SET_ANALOG_PATH
!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

!~       Subroutine TEST_DIR_INT(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E,n_dir,alfa,esc)
!~ !!~ test different directions of propagation and give as output the crack
!~ !!~ opening in the direction that minimise the total energie of the system
!~ !!~ Ut     INPUT/OUTPUT
!~ !!~ Ub     INPUT/OUTPUT
!~ !!~ x0z0   INPUT/OUTPUT
!~ !!~ delta  INPUT/OUTPUT
!~ !!~ Hd_    INPUT/OUTPUT
!~ !!~ A      INPUT/OUTPUT
!~ !!~ n_dis  INPUT/OUTPUT
!~ !!~ E      OUTPUT
!~ !!~ n_dir  INPUT
!~ !!~ alfa   INPUT
!~ !!~ esc    OUTPUT
!~       IMPLICIT NONE
!~ 
!~       REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
!~       REAL(8) :: A_old,E_old,A,E,alfa
!~       INTEGER :: n_dis,n_dir
!~       LOGICAL :: esc
!~ 
!~       REAL(8), ALLOCATABLE :: Ut_dir(:),Ub_dir(:),x0z0_dir(:,:),delta_dir(:),Hd_dir(:)
!~       REAL(8), ALLOCATABLE :: x0z0_in(:,:),delta_in(:),Hd_in(:)
!~       REAL(8) :: A_dir,A_in,E_dir
!~       REAL(8) :: c_1,c_2,Ddelta,delta_1,E_,E_dir_,h_n_dir
!~       INTEGER :: n_dis_dir,n_dis_in,dir
!~       LOGICAL :: esc_dir
!~ 
!~       esc=.true.
!~ 
!~       ALLOCATE (x0z0_dir(n_dis+1,2),delta_dir(n_dis+1),Hd_dir(n_dis+1),Ut_dir(n_dis+1),Ub_dir(n_dis+1))
!~       ALLOCATE (x0z0_in(n_dis,2),   delta_in(n_dis),   Hd_in(n_dis))
!~ 
!~       A_in = A
!~       n_dis_in = n_dis
!~       x0z0_in(:,:) = x0z0(1:n_dis,:)
!~       delta_in(:)  = delta(1:n_dis)
!~       Hd_in(:)     = Hd_(1:n_dis)
!~ 
!~       h_n_dir = (n_dir-1)*0.5D0
!~ 
!~       If (tan(delta(n_dis)).gt.0.D0) Then
!~         Ddelta =-alfa*2.D0
!~         c_1 = delta(n_dis)
!~         c_2 = delta(n_dis)-h_n_dir*alfa*2.D0
!~       Else
!~         Ddelta = alfa*2.D0
!~         c_1 = delta(n_dis)+h_n_dir*alfa*2.D0
!~         c_2 = delta(n_dis)
!~       EndIf
!~       
!~       If (sin(delta(n_dis)).gt.0.D0) Then
!~         delta_1 = 90.D0*pi/180.D0
!~       Else
!~         delta_1 =-90.D0*pi/180.D0
!~       EndIf
!~ 
!~       Write (22,*) '--------------------------------------------------'
!~       Write (22,*) '** TEST INTERFACE DIRECTION **'
!~       dir = 0
!~       Do While (c_1.gt.c_2)
!~ 
!~         dir = dir+1
!~ 
!~         A_dir     = A_in
!~         n_dis_dir = n_dis_in+1
!~ 
!~         x0z0_dir(1:n_dis_in,:) = x0z0_in(:,:)
!~         delta_dir(1:n_dis_in)  = delta_in(:)
!~         Hd_dir(1:n_dis_in)     = Hd_in(:)
!~ 
!~         delta_dir(n_dis_dir)  = delta_1 + Ddelta*(dir-1)
!~         Hd_dir(n_dis_dir)     = Hd_dir(n_dis_dir-1)
!~         x0z0_dir(n_dis_dir,1) = x0z0_dir(n_dis_dir-1,1) -                               &
!~      &                          0.5D0*Hd_dir(n_dis_dir-1)*sin(delta_dir(n_dis_dir-1)) - &
!~      &                          0.5D0*Hd_dir(n_dis_dir)  *sin(delta_dir(n_dis_dir))
!~         x0z0_dir(n_dis_dir,2) = x0z0_dir(n_dis_dir-1,2) -                               &
!~      &                          0.5D0*Hd_dir(n_dis_dir-1)*cos(delta_dir(n_dis_dir-1)) - &
!~      &                          0.5D0*Hd_dir(n_dis_dir)  *cos(delta_dir(n_dis_dir))
!~         
!~         esc_dir = .false.
!~ 
!~         Write (22,*) '--------------------------------------------------'
!~         Write (*,*)  '--------------------------------------------------'
!~ 
!~         Call SOLVE_CRACK(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_dir,esc_dir)
!~         Call COMP_ENERGY(Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir,A_dir,n_dis_dir,E_dir)
!~ 
!~         If (dir.eq.1) Then
!~           E_dir_ = E_dir + E_threshold_int*Hd_dir(n_dis_dir)
!~         Else
!~           If (x0z0_dir(n_dis_dir,2)-0.5D0*Hd_dir(n_dis_dir)*cos(delta_dir(n_dis_dir)).gt.0.D0) Then
!~             E_dir_ = E_dir + E_threshold_1*Hd_dir(n_dis_dir)
!~           Else
!~             E_dir_ = E_dir + E_threshold_2*Hd_dir(n_dis_dir)
!~           EndIf
!~         EndIf
!~ 
!~         Write (22,*) 'Etot (',dir,'/',n_dir,')=',E_dir_,' @ delta=',delta_dir(n_dis_dir)*180.D0/pi
!~         Write (*,*)  'Etot (',dir,'/',n_dir,')=',E_dir_,' @ delta=',delta_dir(n_dis_dir)*180.D0/pi
!~ 
!~         If(((dir.eq.1).or.(E_dir_.lt.E_)).and.(esc_dir.eqv..false.)) Then
!~           Ut    = Ut_dir
!~           Ub    = Ub_dir
!~           x0z0  = x0z0_dir
!~           delta = delta_dir
!~           Hd_   = Hd_dir
!~           A     = A_dir
!~           n_dis = n_dis_dir
!~           E     = E_dir
!~           E_    = E_dir_
!~           esc   = esc_dir
!~         EndIf
!~ 
!~         If (tan(delta(n_dis)).gt.0.D0) Then
!~           c_1 = delta_dir(n_dis_dir)
!~         Else
!~           c_2 = delta_dir(n_dis_dir)
!~         EndIf
!~  
!~       EndDo
!~ 
!~       DEALLOCATE (Ut_dir,Ub_dir,x0z0_dir,delta_dir,Hd_dir)
!~       
!~       If (delta(n_dis).eq.delta_1) int_dir=.true.
!~ 
!~       End Subroutine TEST_DIR_INT

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_ENERGY(Ut,Ub,x0z0,delta,Hd_,A,n_dis,E)
!~ compute total energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ DPvisc INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ E     OUTPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:)
      REAL(8) :: A,E
      INTEGER :: n_dis

      REAL(8) :: U,W

      Call COMP_U(Ut,Ub,x0z0,delta,Hd_,A,n_dis,U)
      Call COMP_W(Ut,Ub,x0z0,delta,Hd_,A,n_dis,W)
      
      E = U + W ! [MPa*km^2*10^3] it will be devided by the disl lenght in km so that the energy release will be in MPa*m

      End Subroutine COMP_ENERGY

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine COMP_DYN_ENERGY(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,E,D)
!~ compute total energy of the system
!~ Ut    INPUT
!~ Ub    INPUT
!~ x0z0  INPUT
!~ delta INPUT
!~ Hd_   INPUT
!~ DPvisc INPUT
!~ A     INPUT
!~ n_dis INPUT
!~ E     OUTPUT
      IMPLICIT NONE
      
      REAL(8) :: Ut(:),Ub(:),x0z0(:,:),delta(:),Hd_(:),DPvisc(:)
      REAL(8) :: A,E
      INTEGER :: n_dis

      REAL(8) :: U,W,D

      Call COMP_U(Ut,Ub,x0z0,delta,Hd_,A,n_dis,U)
      Call COMP_DYN_W(Ut,Ub,x0z0,delta,Hd_,DPvisc,A,n_dis,W)
!~       Call COMP_D_OLD(Ut,Ut_old,x0z0,x0z0_old,Hd_,D)
      Call COMP_D(Ut_old,Hd_old,Ut,Hd_,n_dis,D)
      
      E = U + W ! [MPa*km^2*10^3] it will be devided by the disl lenght in km so that the energy release will be in MPa*m

      End Subroutine COMP_DYN_ENERGY

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      End Module PROPAGATION
