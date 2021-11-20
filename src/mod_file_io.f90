!=======================================================================
! Written by Mark A. Herndon
! Lehigh University, Department of Mechanical Engineering and Mechanics
!=======================================================================
MODULE mod_file_io
    IMPLICIT NONE
    INTEGER :: SOL_WRITE_FREQ, s
    CHARACTER(:), ALLOCATABLE :: input_file_contents
    CHARACTER(LEN=100) :: input_fname

CONTAINS

SUBROUTINE read_input_data
    USE mod_global, ONLY : nt, dt, nv, nvt, nk, kspan, GE,         &
                           Y_0, Z_0, eta_0, zeta_0, GAM, a, a_vec, &
                           ALLOCATE_VARIABLES
    IMPLICIT NONE
    INTEGER :: i, i2

    NAMELIST /CODE_DATA/ nt, dt, nv, nvt, nk, kspan, GE
    NAMELIST /VORTEX_DATA/ Y_0, Z_0, eta_0, zeta_0, GAM, a, a_vec

    !INQUIRE ( FILE = 'input_parameters.dat', SIZE=s )
    INQUIRE ( FILE = TRIM(input_fname), SIZE=s )
    ALLOCATE ( CHARACTER(LEN=s) :: input_file_contents )

    !OPEN(UNIT=3,FILE='input_parameters.dat',ACCESS='STREAM',ACTION='READ',STATUS='OLD')
    OPEN(UNIT=3,FILE=TRIM(input_fname),ACCESS='STREAM',ACTION='READ',STATUS='OLD')
    READ(3) input_file_contents
    CLOSE(3)

    DO i = 1, s
        IF (input_file_contents(i:i) == '!' ) THEN
           i2 = i + INDEX(input_file_contents(i:),ACHAR(10))-2
            input_file_contents(i:i2) = ''
        END IF
    END DO

    READ ( input_file_contents, NML=CODE_DATA )

    CALL ALLOCATE_VARIABLES

    READ ( input_file_contents, NML=VORTEX_DATA )

END SUBROUTINE read_input_data

!=======================================================================
!=======================================================================
SUBROUTINE WRITE_SOLUTION_FILE(ii)
    USE mod_global, ONLY : nv, nvt, nt, Y, Z, tau, eta, zeta, GE, PHI, nk, s, wave_nums, V
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    INTEGER :: zz
    CHARACTER(LEN=60) :: fname_1, fname_2
    zz = INT(Z(1,1));
    IF ( GE == .FALSE. ) THEN
        WRITE(fname_2,'("DATA/perturbations-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) wave_nums
        WRITE(1) V
        CLOSE(1)

        WRITE(fname_2,'("DATA/perturbations_2-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        !WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) wave_nums
        WRITE(1) V
        CLOSE(1)
    ELSE
        WRITE(fname_2,'("DATA/perturbations-GE-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) wave_nums
        WRITE(1) V
        CLOSE(1)

        WRITE(fname_2,'("DATA/perturbations_2-GE-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_2,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt, nk
        !WRITE(1) PHI
        WRITE(1) tau
        WRITE(1) s
        WRITE(1) wave_nums
        WRITE(1) V
        CLOSE(1)
    END IF
END SUBROUTINE WRITE_SOLUTION_FILE
!=======================================================================
!=======================================================================
SUBROUTINE WRITE_TRAJECTORIES(ii)
    USE mod_global, ONLY : nv, nt, Y, Z, tau, GE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ii
    INTEGER :: zz
    CHARACTER(LEN=60) :: fname_1

    zz = INT(Z(1,1));
    IF ( GE == .FALSE. ) THEN
        WRITE(fname_1,'("DATA/vortex_trajectories-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_1,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt
        WRITE(1) Y
        WRITE(1) Z
        WRITE(1) tau
        CLOSE(1)
    ELSE
        WRITE(fname_1,'("DATA/vortex_trajectories-GE-",I4.4,"-",I3.3,".x")'), ii, zz
        OPEN(1,FILE=fname_1,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='REPLACE',ACTION='WRITE')
        WRITE(1) nv, nt
        WRITE(1) Y
        WRITE(1) Z
        WRITE(1) tau
        CLOSE(1)
    END IF
    
END SUBROUTINE WRITE_TRAJECTORIES
END MODULE mod_file_io
