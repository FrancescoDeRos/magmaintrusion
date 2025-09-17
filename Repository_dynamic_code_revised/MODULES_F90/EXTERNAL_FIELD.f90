      Module EXTERNAL_FIELD

      REAL(8), ALLOCATABLE, PRIVATE :: STRESS_MATRIX(:,:,:)
      REAL(8), PRIVATE :: x_min,x_max,z_min,z_max,grd
      contains

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine SET_EX_STRESS_MATRIX(x_min_,x_max_,z_min_,z_max_,grd_,input_stress)

      IMPLICIT NONE

      REAL(8):: x_min_,x_max_,z_min_,z_max_,grd_
      CHARACTER*50 :: input_stress
      
      INTEGER :: n_pnt,n_pnt_x,n_pnt_z,i,j
      REAL(8), ALLOCATABLE :: X(:),Z(:),Sxx(:),Szz(:),Sxz(:)
      INTEGER :: IO

      x_min = x_min_
      x_max = x_max_
      z_min = z_min_
      z_max = z_max_
      grd   = grd_
      Open(49,file=input_stress, access='sequential',form='formatted')

      n_pnt = 0
      Do
        Read(49,*,IOSTAT=IO)
        If (IO.ne.0) Exit
        n_pnt = n_pnt + 1
      EndDo
      Rewind(49)

      ALLOCATE (X(n_pnt),Z(n_pnt),Sxx(n_pnt),Szz(n_pnt),Sxz(n_pnt))

      Do i=1,n_pnt
        Read(49,*) X(i), Z(i), Sxx(i), Szz(i), Sxz(i)
        If (IO.gt.0) Then
          Write(*,*) 'I/O error when reading ', input_stress, 'at line', i
          Stop
        EndIf
      EndDo
      Rewind(49)

      Close(49)

      n_pnt_x = int((x_max-x_min)/grd) + 1
      n_pnt_z = int((z_max-z_min)/grd) + 1
      
      If ((n_pnt_x*n_pnt_z).ne.n_pnt) Then
        Write (*,*) '(n_pnt_x*n_pnt_z) is not equal to n_pnt'
        Stop
      EndIf

      ALLOCATE (STRESS_MATRIX(n_pnt_x,n_pnt_z,5))
      
      Do j=1,n_pnt_z
        Do i=1,n_pnt_x

          STRESS_MATRIX(i,j,1) =   X(n_pnt_x*(j-1)+i)
          STRESS_MATRIX(i,j,2) =   Z(n_pnt_x*(j-1)+i)
          STRESS_MATRIX(i,j,3) = Sxx(n_pnt_x*(j-1)+i)
          STRESS_MATRIX(i,j,4) = Szz(n_pnt_x*(j-1)+i)
          STRESS_MATRIX(i,j,5) = Sxz(n_pnt_x*(j-1)+i)
          
        EndDo
      EndDo

      DEALLOCATE (X,Z,Sxx,Szz,Sxz)

      End Subroutine SET_EX_STRESS_MATRIX

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine EX_STRESS(S_ex,x,z)

      IMPLICIT NONE
      
      REAL(8) :: S_ex(2,2),x,z
      INTEGER :: i,j

      If (ALLOCATED(STRESS_MATRIX).eqv..false.) Then
        write(*,*) 'STRESS_MATRIX has not been set'
        stop
      EndIf

      i = nint((x-x_min)/grd) + 1
      j = nint((z_max-z)/grd) + 1

      S_ex = 0.D0
      S_ex(1,1) = STRESS_MATRIX(i,j,3)      
      S_ex(2,2) = STRESS_MATRIX(i,j,4)      
      S_ex(1,2) = STRESS_MATRIX(i,j,5)      
      S_ex(2,1) = S_ex(1,2)      

!~       write(*,*) 'x=',x,'   x_matrix=',STRESS_MATRIX(i,j,1)
!~       write(*,*) 'z=',z,'   z_matrix=',STRESS_MATRIX(i,j,2)
!~       write(*,*) 'i=',i,'   j=',j
!~       write(*,*) ''


      End Subroutine EX_STRESS

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      Subroutine DEALLOC_EX_STRESS_MATRIX()

      If (ALLOCATED(STRESS_MATRIX).eqv..true.) DEALLOCATE (STRESS_MATRIX)

      End Subroutine DEALLOC_EX_STRESS_MATRIX

!~ ---------------------------------------------------------------------
!~ ---------------------------------------------------------------------

      End Module EXTERNAL_FIELD
