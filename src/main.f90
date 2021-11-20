!=======================================================================
!=======================================================================
PROGRAM MAIN
    USE mod_file_io,   ONLY : read_input_data, WRITE_TRAJECTORIES,      &
                              WRITE_SOLUTION_FILE, input_fname
    USE mod_global,    ONLY : wave_nums, ka, omega, omega_array, nk, a, &
                              mutual_induction, GE, s
    USE omp_lib
    IMPLICIT NONE
    INTEGER :: ii, jj, kk
    PROCEDURE(mutual_induction) :: PSI, PHI
    REAL(KIND=8) :: tstart, tend, tfinal, ts, te, tf
    CALL GET_COMMAND_ARGUMENT(1,input_fname)
    CALL read_input_data
    CALL GENERATE_WAVE_NUMS
    
    ! Generate omega array
    omega_array(1) = 0.d0
    DO ii = 2, nk
        ka = wave_nums(ii)
        CALL CALC_OMEGA
        omega_array(ii) = omega
    END DO
    
    OPEN(1,FILE='DATA/OMEGA.x',FORM='UNFORMATTED',ACCESS='STREAM',ACTION='WRITE',STATUS='REPLACE')
    WRITE(1) nk
    WRITE(1) wave_nums
    WRITE(1) omega_array
    CLOSE(1)

    CALL CALC_CROW
    ! Set ground effect constraints if GE == TRUE  
    IF ( GE == .TRUE. ) THEN
        CALL SET_GROUND_EFFECT
    END IF

    ! Calculate trajectories of vortices ( or set constant position )
    CALL CALC_TRAJECTORIES
    CALL WRITE_TRAJECTORIES(nk)
    tstart = OMP_get_wtime()
    DO ii = 1, nk
        ka = wave_nums(ii)
        omega = omega_array(ii)
        CALL CALC_PERTURBATIONS(ii)
        CALL CALC_SINGULAR_VALUES(ii)
        tend = OMP_get_wtime()
        tf   = tend - tstart
        IF (MOD(ii,20) .EQ. 0 ) THEN
            WRITE(*,'(A,I4,X,A,X,I4,3X,A,3X,F12.6,3X,A)') 'COMPLETED ITERATION ', ii, '/',nk,'. . .', tf, 'SECONDS ELAPSED'
        END IF
    END DO    
    CALL WRITE_SOLUTION_FILE(nk)


END PROGRAM MAIN
!=======================================================================
!=======================================================================
SUBROUTINE CALC_SINGULAR_VALUES(ii)
    USE mod_global, ONLY : m, phi, s, nt, V
    USE LAPACK95, ONLY : GESVD
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    INTEGER :: i, j, ni, nj, nn
    INTEGER(KIND=8) :: INFO
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A_tmp, U, VT
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: s_array
    CHARACTER(LEN=1) :: JOB = 'U'
    ALLOCATE(A_tmp(m,m))
    ALLOCATE(s_array(m))
    ALLOCATE(VT(m,m))
    
    DO nn = 1, nt
        A_tmp(:,:) = phi(:,:,nn,ii)
!        WRITE(*,'("Matrix A"/(<m>F8.4))'), ((A_tmp(i,j), i = 1, m), j = 1, m)
        CALL GESVD(A=A_tmp,S=s_array,VT=VT,JOB=JOB,INFO=INFO)
    !    WRITE(*,'("Matrix V"/(<m>F8.4))'), ((VT(i,j), i = 1, m), j = 1, m)
!        WRITE(*,'("Singular Values"/(1F8.4))'), (s_array(i), i = 1, m)
        s(ii,nn) = s_array(1)
        DO j = 1, m
            V(j,nn,ii) = VT(j,1)
        !    WRITE(*,*) VT(j,1)
        END DO
    END DO
    
    
    DEALLOCATE(A_tmp); DEALLOCATE(s_array); DEALLOCATE(VT)
END SUBROUTINE CALC_SINGULAR_VALUES
!=======================================================================
!=======================================================================
SUBROUTINE CALC_PERTURBATIONS(ii)
    USE mod_numerical_routines, ONLY : DERIVATIVE, RK5
    USE mod_global,             ONLY : nv, m, dt, eta, zeta, eta_0, zeta_0, &
                                       Vp_0, Vp_new, tau, n, nt, PHI
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: ii
    PROCEDURE(DERIVATIVE) :: VORTEX_DERIV
    INTEGER, SAVE               :: ii_2, jj, i, kk, mm
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), SAVE   :: I_mat
    ! Create identity matrix for initial conditions on PHI
    ALLOCATE(I_mat(m,m))
    ii_2 = ii;
    I_mat = 0.d0
    !WRITE(*,'("Matrix I"/(<m>F8.4))'), ((I_mat(mm,jj), mm = 1, m), jj = 1, m)
    !WRITE(*,*) '--------------------'
    DO jj = 1, m
        I_mat(jj,jj) = 1.d0
    END DO   
    !WRITE(*,'("Matrix I"/(<m>F8.4))'), ((I_mat(mm,jj), mm = 1, m), jj = 1, m)
    DO jj = 1, m
        Vp_0(:)   = I_mat(:,jj)
        eta_0(:)  = I_mat(1:m/2,jj)
        zeta_0(:) = I_mat(m/2 +1:m,jj)
        eta(:,1)  = I_mat(1:m/2,jj)
        zeta(:,1) = I_mat(m/2 +1:m,jj)
        DO n = 1, nt
            CALL Rk5(Vp_0,Vp_new,dt,m,VORTEX_DERIV)
            DO i = 1, nv
                Vp_0(i)         = Vp_new(i)
                Vp_0(i+nv)      = Vp_new(i+nv)
                eta(i,n+1)      = Vp_new(i)
                zeta(i,n+1)     = Vp_new(i+nv)
            END DO
        END DO
        DO mm = 1, m/2
            PHI(mm,jj,:,ii_2)   = eta(mm,:)
            PHI(mm+2,jj,:,ii_2) = zeta(mm,:)
        END DO
    END DO

    DEALLOCATE(I_mat)
END SUBROUTINE CALC_PERTURBATIONS
!=======================================================================
!=======================================================================
SUBROUTINE CALC_TRAJECTORIES
    USE mod_numerical_routines, ONLY : RK5, DERIVATIVE
    USE mod_global,             ONLY : nt, dt, Y_0, Z_0, Y, Z, Vc_0, &
                                       Vc_new, m, nv, tau
    IMPLICIT NONE
    PROCEDURE(DERIVATIVE) :: VORTEX_T_DERIV
    INTEGER :: i, n
    DO i = 1, nv
        Vc_0(i)    = Y_0(i)
        Vc_0(i+nv) = Z_0(i)
        Y(i,1)     = Y_0(i)
        Z(i,1)     = Z_0(i)
        !Y(i,:)     = Y_0(i)
        !Z(i,:)     = Z_0(i)
    END DO
    
    tau(1) = 0.0

    DO n = 1, nt
       CALL Rk5(Vc_0,Vc_new,dt,m,VORTEX_T_DERIV)
       DO i = 1, nv
           Vc_0(i)    = Vc_new(i)
           Vc_0(i+nv) = Vc_new(i+nv)
           Y(i,n+1)   = Vc_new(i)
           Z(i,n+1)   = Vc_new(i+nv)
       END DO
       tau(n+1) = dt*REAL(n,KIND=8)
   END DO

END SUBROUTINE CALC_TRAJECTORIES
!=======================================================================
!=======================================================================
FUNCTION VORTEX_DERIV(x_0,m,h,ch)
    USE mod_global, ONLY : GE, pi, nv, nvt, GAM, mutual_induction, ka, omega, a, Y, Z, n, b_0
    IMPLICIT NONE
    INTEGER,                    INTENT(IN)    :: m
    REAL(KIND=8),               INTENT(IN)    :: h, ch
    REAL(KIND=8), DIMENSION(m), INTENT(IN)    :: x_0
    REAL(KIND=8), DIMENSION(m)                :: VORTEX_DERIV
    ! FUNCTION SPECIFIC VARIABLES AND PARAMETERS
    PROCEDURE(mutual_induction)               :: PSI, PHI    !< Special functions
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_temp      !< Temporary y array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_temp      !< Temporary z array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: eta_temp    !< Temporary eta array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: zeta_temp   !< Temporary zeta array

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_deriv     !< y derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_deriv     !< z derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: eta_deriv   !< eta derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: zeta_deriv  !< zeta derivative array
    INTEGER                                   :: i, j        !< Loop index integers
    REAL(KIND=8)                              :: y_mn, z_mn  !< Relative y and z coordinates
    REAL(KIND=8)                              :: r_mn        !< Relative radius
    REAL(KIND=8)                              :: sum_y       !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_z       !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_eta     !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_zeta    !< Sum holder for y eq
    REAL(KIND=8)                              :: V1_mn       !< First term for eta equation
    REAL(KIND=8)                              :: V2_mn       !< Second term for eta equation
    REAL(KIND=8)                              :: V3_mn       !< Third term for eta equation
    REAL(KIND=8)                              :: V4_mn       !< Fourth term for eta equation
    REAL(KIND=8)                              :: W1_mn       !< First term for zeta equation
    REAL(KIND=8)                              :: W2_mn       !< Second term for zeta equation
    REAL(KIND=8)                              :: W3_mn       !< Third term for zeta equation
    REAL(KIND=8)                              :: W4_mn       !< Fourth term for zeta equation
    REAL(KIND=8)                              :: k
    
    k = ka/a;

    ALLOCATE(eta_temp(nvt))
    ALLOCATE(zeta_temp(nvt))
    ALLOCATE(y_temp(nvt))
    ALLOCATE(z_temp(nvt))
    ALLOCATE(eta_deriv(nv))
    ALLOCATE(zeta_deriv(nv))

    ! Initialize sum holder values to 0.0
    sum_eta = 0.0; sum_zeta = 0.0;

    ! Map x_0 values into temporary arrays
    DO i = 1, nv
        y_temp(i)    = Y(i,n) + ch*(Y(i,n+1)-Y(i,n))/h 
        z_temp(i)    = Z(i,n) + ch*(Z(i,n+1)-Z(i,n))/h 
        eta_temp(i)  = x_0(i)
        zeta_temp(i) = x_0(i+nv)
    END DO

    IF ( GE == .TRUE. ) THEN
     DO i = nv+1, nvt
         y_temp(i)     =  y_temp(i-nv)
         z_temp(i)     = -z_temp(i-nv)
         eta_temp(i)   =  x_0(i-nv)
         zeta_temp(i)  = -x_0(i)
     END DO
    END IF

    DO i = 1, nv
        DO j = 1, nvt
            IF ( i .EQ. j ) THEN
                CYCLE
            ELSEIF ( i .NE. j ) THEN
                z_mn  = z_temp(j) - z_temp(i)
                y_mn  = y_temp(j) - y_temp(i)
                r_mn  = SQRT(y_mn**2 + z_mn**2)

                V1_mn =   GAM(j)*2.d0*y_mn*z_mn/(r_mn**4)
                V2_mn = -(GAM(j)*2.d0*y_mn*z_mn/(r_mn**4))*PHI(k*r_mn) 
                V3_mn = -(GAM(j)/(r_mn**2))*(1.d0 - (2.d0*z_mn**2/r_mn**2))
                V4_mn =  (GAM(j)/(r_mn**2))*(PSI(k*r_mn) - ((2.d0*z_mn**2)/r_mn**2)*PHI(k*r_mn))

                W1_mn =  -GAM(j)*2.d0*y_mn*z_mn/(r_mn**4)
                W2_mn =  (GAM(j)*2.d0*y_mn*z_mn/(r_mn**4))*PHI(k*r_mn) 
                W3_mn =  (GAM(j)/(r_mn**2))*(1.d0 - (2.d0*y_mn**2/r_mn**2))
                W4_mn = -(GAM(j)/(r_mn**2))*(PSI(k*r_mn) - ((2.d0*y_mn**2)/r_mn**2)*PHI(k*r_mn))

                sum_eta  = sum_eta  + V1_mn*eta_temp(i)  + V2_mn*eta_temp(j)  + V3_mn*zeta_temp(i) + V4_mn*zeta_temp(j)
                sum_zeta = sum_zeta + W1_mn*zeta_temp(i) + W2_mn*zeta_temp(j) + W3_mn*eta_temp(i)  + W4_mn*eta_temp(j)
            END IF
        END DO
        eta_deriv(i)  = sum_eta  + (GAM(i)/(a**2))*omega*zeta_temp(i) 
        zeta_deriv(i) = sum_zeta - (GAM(i)/(a**2))*omega*eta_temp(i) 
        sum_eta       = 0.d0
        sum_zeta      = 0.d0
    END DO

    DO i = 1, nv
        VORTEX_DERIV(i) = eta_deriv(i)
        VORTEX_DERIV(i+nv) = zeta_deriv(i)
    END DO
    

    DEALLOCATE(y_temp)
    DEALLOCATE(z_temp)
    DEALLOCATE(eta_temp)
    DEALLOCATE(zeta_temp)
    DEALLOCATE(eta_deriv)
    DEALLOCATE(zeta_deriv)
END FUNCTION VORTEX_DERIV
!=======================================================================
!=======================================================================
FUNCTION VORTEX_T_DERIV(x_0,m,h,ch)
    USE mod_global, ONLY : GE, pi, nv, nvt, GAM, ka, omega, a
    IMPLICIT NONE
    INTEGER,                    INTENT(IN)    :: m
    REAL(KIND=8),               INTENT(IN)    :: h, ch
    REAL(KIND=8), DIMENSION(m), INTENT(IN)    :: x_0
    REAL(KIND=8), DIMENSION(m)                :: VORTEX_T_DERIV
    ! FUNCTION SPECIFIC VARIABLES AND PARAMETERS
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_temp   !< Temporary y array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_temp   !< Temporary z array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: y_deriv  !< y derivative array
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:)   :: z_deriv  !< z derivative array
    INTEGER                                   :: i, j     !< Loop index integers
    REAL(KIND=8)                              :: y_mn     !< Relative y coordinates
    REAL(KIND=8)                              :: z_mn     !< Relative z coordinates
    REAL(KIND=8)                              :: r_mn     !< Relative radius
    REAL(KIND=8)                              :: sum_y    !< Sum holder for y eq
    REAL(KIND=8)                              :: sum_z    !< Sum holder for y eq

    ALLOCATE(y_temp(nvt))
    ALLOCATE(z_temp(nvt))

    ALLOCATE(y_deriv(nv))
    ALLOCATE(z_deriv(nv))

    ! Initialize sum holder values to 0.0
    sum_y = 0.0; sum_z = 0.0;

    ! Map x_0 values into temporary arrays
    DO i = 1, nv
        y_temp(i)    = x_0(i)
        z_temp(i)    = x_0(i+nv)
    END DO
    IF ( GE == .TRUE. ) THEN
     DO i = nv+1, nvt
         y_temp(i)    =  x_0(i-nv)
         z_temp(i)    = -x_0(i)
     END DO
    END IF

    DO i = 1, nv
        DO j = 1, nvt
            IF ( i .EQ. j ) THEN
                CYCLE
            ELSEIF ( i .NE. j ) THEN
                z_mn  = z_temp(j) - z_temp(i)
                y_mn  = y_temp(j) - y_temp(i)
                r_mn  = SQRT(y_mn**2 + z_mn**2)
                
                sum_y = sum_y + GAM(j)*z_mn/r_mn**2
                sum_z = sum_z - GAM(j)*y_mn/r_mn**2
            END IF
        END DO
        y_deriv(i)    = sum_y
        z_deriv(i)    = sum_z
        sum_y         = 0.d0
        sum_z         = 0.d0
    END DO
    DO i = 1, nv
        VORTEX_T_DERIV(i)    = y_deriv(i)
        VORTEX_T_DERIV(i+nv) = z_deriv(i)
    END DO
    
    DEALLOCATE(y_temp)
    DEALLOCATE(z_temp)
    DEALLOCATE(y_deriv)
    DEALLOCATE(z_deriv)
END FUNCTION VORTEX_T_DERIV
!=======================================================================
!=======================================================================
SUBROUTINE SET_GROUND_EFFECT
    USE mod_global, ONLY : nv, nvt, GAM
    IMPLICIT NONE
    INTEGER :: j

    DO j = nv+1, nvt
        GAM(j) = -GAM(j-nv)
    END DO

END SUBROUTINE SET_GROUND_EFFECT
!=======================================================================
!=======================================================================
FUNCTION PSI(BETA)
    USE special_function_interface, ONLY : BESSELK0, BESSELK1
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: BETA
    REAL(KIND=8)             :: PSI
    
    IF ( BETA == 0.d0 ) THEN
        PSI = 1.d0
    ELSE IF ( BETA .GE. 600 ) THEN
        PSI = 0.d0
    ELSE
        PSI = (BETA**2)*BESSELK0(ABS(BETA)) + & 
              ABS(BETA)*BESSELK1(ABS(BETA))
    END IF
END FUNCTION PSI
!=======================================================================
!=======================================================================
FUNCTION PHI(BETA)
    USE special_function_interface, ONLY : BESSELKN
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: beta
    REAL(KIND=8)             :: PHI

    IF ( BETA == 0.d0 ) THEN
        PHI = 1.d0
    ELSE IF ( BETA .GE. 600.d0 ) THEN
        PHI = 0.d0
    ELSE 
        PHI = (1.d0/2.d0)*(BETA**2)*BESSELKN(2,ABS(BETA))
    END IF
END FUNCTION PHI
!=======================================================================
!=======================================================================
SUBROUTINE CALC_OMEGA
    USE mod_numerical_routines, ONLY : BISECTION_METHOD, root_function
    USE mod_global, ONLY : omega, ka
    IMPLICIT NONE
    PROCEDURE(root_function) :: BESSEL_ROOT
    PROCEDURE(root_function) :: DISPERSION
    REAL(KIND=8) :: x !< return value for bisection root finding method
    REAL(KIND=8) :: a !< Left intial end point for bisection method
    REAL(KIND=8) :: b !< Right intial end point for bisection method
    REAL(KIND=8) :: tol = 1E-14 !< Tolerance for root approximation
    REAL(KIND=8) :: eps = 1E-10; !< Adjust bessl_root_val by eps
    REAL(KIND=8) :: bessel_root_val !< clearly indicate root of bessel
    
    a = 0.01d0; b = 5.d0; ! First root of bessel J in this interval
    CALL BISECTION_METHOD(BESSEL_ROOT, x, a, b, tol)
    bessel_root_val = x;
    b = bessel_root_val - eps;
    CALL BISECTION_METHOD(DISPERSION, x, a, b, tol)
    omega = ((2*ka/SQRT(ka**2 + x**2)) - 1.d0);
    
END SUBROUTINE CALC_OMEGA
!=======================================================================
!=======================================================================
FUNCTION DISPERSION(BETA)
    USE special_function_interface, ONLY : BESSELJ0, BESSELJ1, BESSELJN, &
                                           BESSELK0, BESSELK1, BESSELKN
    USE mod_global, ONLY : ka
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: BETA
    REAL(KIND=8)             :: DISPERSION
    REAL(KIND=8)             :: J0, J1, J2, J1_p 
    REAL(KIND=8)             :: K0, K1, K2, K1_p 

    J0 = BESSELJ0(BETA); J1 = BESSELJ1(BETA); J2 = BESSELJN(2,BETA);
    K0 = BESSELK0(ka);   K1 = BESSELK1(ka);   K2 = BESSELKN(2,ka);

    J1_p =  (J0 - J2)/2.d0;
    K1_p = -(K0 + K2)/2.d0;

    DISPERSION = (1.d0/beta)*(J1_p/J1) + K1_p/(ka*K1) + &
                 SQRT(beta**2 + (ka)**2)/(ka*beta**2) 
END FUNCTION DISPERSION
!=======================================================================
!=======================================================================

FUNCTION BESSEL_ROOT(x)
    USE special_function_interface, ONLY : BESSELJ1
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: x
    REAL(KIND=8)             :: BESSEL_ROOT
    REAL(KIND=8)             :: y
    
    y = BESSELJ1(x)
    BESSEL_ROOT = y

END FUNCTION BESSEL_ROOT
!=======================================================================
!=======================================================================
SUBROUTINE GENERATE_WAVE_NUMS
    USE mod_global, ONLY : wave_nums, nk, kspan
    IMPLICIT NONE
    INTEGER      :: k
    REAL(KIND=8) :: dk !< Wave number step size
    
    dk = kspan/REAL(nk,KIND=8)
    DO k = 1, nk
        wave_nums(k) = REAL(k,KIND=8)*dk
    END DO
    DO k = 1, nk
        wave_nums(k) = wave_nums(k) + wave_nums(5);
    END DO

!    wave_nums(1) = 0.1432d0
    !wave_nums(1) = 0.1095d0
END SUBROUTINE GENERATE_WAVE_NUMS
!=======================================================================
!=======================================================================
SUBROUTINE CALC_CROW
    USE mod_global, ONLY : wave_nums, nk, a, ka, omega, omega_array, &
                           mutual_induction
    
    IMPLICIT NONE
    PROCEDURE(mutual_induction)               :: PSI, PHI    !< Special functions
    INTEGER :: ii, jj
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: alpha
    REAL(KIND=8) :: beta, CHI
    ALLOCATE(alpha(nk))

    DO jj = 1, nk
        beta = wave_nums(jj)/a
        omega = omega_array(jj)/(a**2)
        CHI = -(PSI(beta) - 2.d0*PHI(beta))
        alpha(jj) = ((1.d0 - PSI(beta) + omega)*(1.d0 + CHI - omega))
    END DO
    
    OPEN(1,FILE='DATA/alpha.x',FORM='UNFORMATTED',ACCESS='STREAM',ACTION='WRITE',STATUS='REPLACE')
    WRITE(1) nk
    WRITE(1) wave_nums
    WRITE(1) alpha
    WRITE(1) omega_array
    CLOSE(1)
END SUBROUTINE CALC_CROW
