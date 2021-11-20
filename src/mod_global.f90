!======================================================================
! Written by Mark A. Herndon
! Lehigh University, Department of Mechanical Engineering and Mechanics
!======================================================================
MODULE mod_global
    IMPLICIT NONE
    ! CONSTANTS
    REAL(KIND=8), PARAMETER :: pi = 4.0*ATAN(1.0)

    ! USER-SPECIFIED PARAMETERS
    INTEGER      :: nt   !< # of time steps
    INTEGER      :: nv   !< # of vortices in real plane
    INTEGER      :: nk   !< # of wavenumbers
    INTEGER      :: nvt  !< Total # of vortices in ground-image system
    REAL(KIND=8) :: kspan !< Wavenumber span (kb)
    LOGICAL      :: GE   !< In Ground Effect Logical

    REAL(KIND=8) :: dt  !< Time step
    REAL(KIND=8) :: b_0 !< Initial vortex separation
    REAL(KIND=8) :: a   !< Vortex core radius 
    REAL(KIND=8) :: ka  !< Vortex core radius time wavenumber k
    REAL(KIND=8) :: omega !< Self induced rotation frequency 

    ! USER-SPECIFIED INITIAL CONDITIONS
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Y_0    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Z_0    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: eta_0  
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: zeta_0 
    !$OMP THREADPRIVATE(eta_0, zeta_0,Y_0,Z_0)

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: omega_array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: wave_nums
    ! VORTEX POSITION AND PERTURBATION AMPLITUDE ARRAYS
    ! DIMENSION(nvt,nt) --> Ex. Y(vortex index, time index) == position
    ! vortex (vortex index) at t = (time index)
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Y    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: Z    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: eta  
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: zeta 

    !$OMP THREADPRIVATE(Y,Z,eta,zeta)
    ! VORTEX CIRCULATION STRENGTH AND ORIENTATION
    ! DIMENSION(nvt) --> Ex. GAM(vortex index 1) ...
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: GAM   
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: a_vec   

    ! DERIVED PARAMETERS
    REAL(KIND=8) :: b     
    REAL(KIND=8) :: h     
    INTEGER      :: m     !< Dimension of VORT array
    ! GLOBAL VARIABLES
    INTEGER :: n    !< Time integration indexing integer
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Vc_0
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Vc_new
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Vp_0
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Vp_new
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tau
    !$OMP THREADPRIVATE(Vc_0, Vc_new, Vp_0, Vp_new)    
    ! PROPAGATOR MATRIX
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: s
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: PHI
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)   :: C_MAT    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)   :: V    
    !$OMP THREADPRIVATE(s, PHI, C_MAT)

    ABSTRACT INTERFACE
        FUNCTION mutual_induction(beta)
            REAL(KIND=8), INTENT(IN) :: beta
            REAl(KIND=8)             :: mutual_induction
        END FUNCTION
    END INTERFACE
CONTAINS

!======================================================================
!======================================================================
SUBROUTINE ALLOCATE_VARIABLES
    IMPLICIT NONE

    ! ALLOCATE INITIAL POSITION ARRAYS
    ALLOCATE(Y_0(nv))
    ALLOCATE(Z_0(nv))
    ALLOCATE(eta_0(nv))
    ALLOCATE(zeta_0(nv))

    ! ALLOCATE VORTEX POSITION AND PERTURBATION AMPLITUDE HISTORY ARRAYS
    ALLOCATE(Y(nv,nt))
    ALLOCATE(Z(nv,nt))
    ALLOCATE(eta(nv,nt))
    ALLOCATE(zeta(nv,nt))

    ! ALLOCATE VORTEX CIRCULATION ARRAY
    ALLOCATE(GAM(nvt))
    ALLOCATE(a_vec(nvt))

    ! ALLOCATE VORT ARRAYS FOR RK5 INTEGRATOR
    m = nv*2
    ALLOCATE(Vc_0(m))
    ALLOCATE(Vc_new(m))
    ALLOCATE(Vp_0(m))
    ALLOCATE(Vp_new(m))
    ALLOCATE(tau(nt))
    
    ALLOCATE(PHI(m,m,nt,nk))
    ALLOCATE(C_MAT(m,m,nk))
    
    ! omega / wavenumber arrays
    ALLOCATE(omega_array(nk))
    ALLOCATE(wave_nums(nk))
    ALLOCATE(s(nk,nt))
    ALLOCATE(V(m,nt,nk))

END SUBROUTINE ALLOCATE_VARIABLES
!======================================================================
!======================================================================
END MODULE mod_global
