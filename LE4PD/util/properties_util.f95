! Original FORTRAN codes by Jeremy Copperman & Marina Guenza
! Author: Pablo Romano
! Github: @pgromano
! Last Modified: 06/29/2016

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE a MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_a_MATRIX(a, M, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nres
    DOUBLE PRECISION,INTENT(in) :: M(1:nres,1:nres)
    DOUBLE PRECISION,INTENT(inout) :: a(1:nres-1,1:nres)

    INTEGER(kind=8) :: i,j,k,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    k = 0
    n = 1
    DO i=1,nres
        IF(i.eq.n)THEN
            n=n+nres
        ELSE
            k=k+1
            DO j=1,nres
                a(k,j) = M(i,j)
            END DO
        END IF
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE ALPHA CARBON BOND VECTORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_BOND_VECTORS(bond, length, C, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: C(1:nconf, 1:nres, 1:3)
    DOUBLE PRECISION,INTENT(inout) :: bond(1:nconf, 1:nres-1, 1:3)
    DOUBLE PRECISION,INTENT(inout) :: length(1:nconf, 1:nres-1)

    INTEGER(kind=8) :: i,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO n=1,nconf
        DO i=1,nres-1
            bond(n,i,1) = C(n,i+1,1)-C(n,i,1)
            bond(n,i,2) = C(n,i+1,2)-C(n,i,2)
            bond(n,i,3) = C(n,i+1,3)-C(n,i,3)
            length(n,i) = (bond(n,i,1)**2+bond(n,i,2)**2+bond(n,i,3)**2)**0.5
        END DO
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCUALTE PROTEIN ALPHA-CARBON CENTER OF MASS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_CA_COM(COM,C,nconf,nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER,INTENT(in) :: nconf,nres
    DOUBLE PRECISION,INTENT(in) :: C(1:nconf,1:nres,1:3)
    DOUBLE PRECISION,INTENT(inout) :: COM(1:nconf,1:3)

    INTEGER :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    COM = 0
    DO i=1,nconf
        DO j=1,nres
            COM(i,1) = COM(i,1) + C(i,j,1)
            COM(i,2) = COM(i,2) + C(i,j,2)
            COM(i,3) = COM(i,3) + C(i,j,3)
        END DO
        COM(i,1) = COM(i,1)/nres
        COM(i,2) = COM(i,2)/nres
        COM(i,3) = COM(i,3)/nres
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE GAUSSIAN NETWORK MODEL CONTACTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_CONTACTS(contacts, eig, eigv, R, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: R(1:nconf, 1:nres, 1:nres)
    DOUBLE PRECISION,INTENT(inout) :: contacts(1:nconf, 1:nres, 1:nres)
    DOUBLE PRECISION,INTENT(inout) :: eig(1:nconf, 1:nres)
    DOUBLE PRECISION,INTENT(inout) :: eigv(1:nconf, 1:nres, 1:nres)

    INTEGER(kind=8) :: i,j,n,io
    DOUBLE PRECISION :: csum(1:nres), rc, work(1:nres*nres)
    DOUBLE PRECISION :: eig_in(1:nres), eigv_in(1:nres, 1:nres)
    rc = 0.7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contacts = 0.0
    DO n=1,nconf
        DO i=1,nres
            DO j=1,nres
                IF(R(n,i,j).ne.0.0)THEN
                    IF((i.ne.j).AND.(1.0/R(n,i,j).le.rc))THEN
                        contacts(n,i,j) = -1.0
                    END IF
                END IF
            END DO
        END DO
    END DO

    DO n=1,nconf
        csum = 0.0
        DO i=1,nres
            DO j=1,nres
                IF(i.ne.j)THEN
                    csum(i) = csum(i)-contacts(n,i,j)
                END IF
            END DO
        END DO

        DO i=1,nres
            contacts(n,i,i) = csum(i)
        END DO
    END DO

    !WRITE(*,*)"io = 0:  (success)"
    !WRITE(*,*)"io < 0:  (illegal value)"
    !WRITE(*,*)"io > 0:  (no convergence)"
    DO n=1,nconf
        eig_in = 0
        eigv_in = contacts(n,1:nres,1:nres)
        CALL DSYEV("V","U",nres,eigv_in,nres,eig_in,work,nres*nres,io)
        !WRITE(*,*)"ssyev: ",io
        DO i=1,nres
            eig(n,i) = eig_in(i)
            DO j=1,nres
                eigv(n,i,j) = eigv_in(i,j)
            END DO
        END DO
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE H MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_H_MATRIX(H, R, fp, fratio, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nres
    DOUBLE PRECISION,INTENT(in) :: R(1:nres,1:nres), fp(1:nres), fratio
    DOUBLE PRECISION,INTENT(inout) :: H(1:nres,1:nres)

    INTEGER(kind=8) :: i,j,n
    DOUBLE PRECISION :: avg_fppr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    avg_fppr = SUM(fp)/FLOAT(nres)  ! Average protein-contributed friction per residue
    DO i=1,nres
        H(i,i) = avg_fppr/fp(i)
        IF(i.ne.1)THEN
            DO j=1,i-1
                H(i,j) = fratio*R(i,j)
                H(j,i) = H(i,j)
            END DO
        END IF
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE ATOMIC MODE LENGTH ARRAY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_MODE_LENGTH_ARRAY(mlen_out, mlen_in, apr, natoms, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: natoms, nres
    INTEGER,INTENT(in) :: apr(1:nres) ! atoms per residue
    DOUBLE PRECISION,INTENT(in) :: mlen_in(1:nres-1,1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: mlen_out(1:nres-1,1:natoms)

    INTEGER(kind=8) :: i,j,k,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO k=1,nres-1                       ! MODE LOOP
        n=0
        DO i=1,nres                     ! RESIDUE LOOP
            DO j=1,apr(i)               ! ATOM LOOP
                n=n+1
                IF(i.eq.1)THEN
                    mlen_out(k,n) = 0.d0
                ELSE
                    mlen_out(k,n) = mlen_in(k,i-1)
                END IF
            END DO
        END DO
    END DO
    END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE MODE TRAJECTORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_MODE_TRAJ(mode, bonds, QI, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: bonds(0:nconf-1,1:nres-1,1:3), QI(1:nres-1,1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: mode(0:nconf-1,1:nres-1,1:2)

    INTEGER(kind=8) :: a,n,t
    DOUBLE PRECISION :: pi, mode_vec(0:nconf-1,1:nres-1,1:3),mode_mag(0:nconf-1,nres-1)
    pi = 4.D0*DATAN(1.D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    mode = 0.d0
    mode_vec = 0.d0
    mode_mag = 0.d0
    DO t=0,nconf-1
        ! Calculate instantaneous mode vector
        DO a=1,nres-1
            DO n=1,nres-1
                mode_vec(t,a,1) = QI(a,n)*bonds(t,n,1)+mode_vec(t,a,1) + mode_vec(t,a,1)
                mode_vec(t,a,2) = QI(a,n)*bonds(t,n,2)+mode_vec(t,a,2) + mode_vec(t,a,2)
                mode_vec(t,a,3) = QI(a,n)*bonds(t,n,3)+mode_vec(t,a,3) + mode_vec(t,a,3)
            END DO
            mode_mag(t,a) = (mode_vec(t,a,1)**2 +&
                            &mode_vec(t,a,2)**2 +&
                            &mode_vec(t,a,3)**2)**0.5
        END DO

        ! Calculate Phi & Theta
        DO a=1,nres-1
            ! Calculate Phi
            mode(t,a,1) = DATAN(mode_vec(t,a,2)/mode_vec(t,a,1))
            IF(mode_vec(t,a,1).lt.0.d0)THEN
                mode(t,a,1) = mode(t,a,1)+pi
            END IF
            IF(mode(t,a,1).lt.0.d0)THEN
                mode(t,a,1) = mode(t,a,1)+2.d0*pi
            END IF

            ! Calculate Theta
            mode(t,a,2) = DACOS(mode_vec(t,a,3)/mode_mag(t,a))
        END DO
    END DO
    END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE MEAN SQUARED FLUCTUATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_MSF(MSF, contacts, eig, temp, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: contacts(1:nconf, 1:nres, 1:nres)
    DOUBLE PRECISION,INTENT(in) :: temp, eig(1:nconf, 1:nres)
    DOUBLE PRECISION,INTENT(inout) :: MSF(1:nconf, 1:nres, 1:nres)

    INTEGER(kind=8) :: i,j,k,n
    DOUBLE PRECISION :: Rb, c

    Rb = 0.0019872
    c = (3.0*Rb*temp)/6.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO n=1,nconf
        DO i=1,nres
            DO j=1,nres
                DO k=2,nres
                    MSF(n,i,j) = MSF(n,i,j)+&
                    &(contacts(n,i,k)*contacts(n,j,k))/eig(n,k)
                END DO
            END DO
        END DO
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROJECT BOND DYNAMICS ONTO NH BOND VECTORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_NH_MATRIX(mat,C_bond,NH_bond,C_length,NH_length,nconf,nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: C_bond(1:nconf, 1:nres-1, 1:3)
    DOUBLE PRECISION,INTENT(in) :: C_length(1:nconf, 1:nres-1)
    DOUBLE PRECISION,INTENT(in) :: NH_bond(1:nconf, 1:nres-1, 1:3)
    DOUBLE PRECISION,INTENT(in) :: NH_length(1:nconf, 1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: mat(1:nconf, 1:nres-1, 1:nres-1)

    INTEGER(kind=8) :: i,j,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO n=1,nconf
        DO i=1,nres-1
            DO j=1,nres-1
                mat(n,i,j) = (NH_bond(n,i,1)*C_bond(n,j,1) +&
                            & NH_bond(n,i,2)*C_bond(n,j,2) +&
                            & NH_bond(n,i,3)*C_bond(n,j,3)) /&
                            &(NH_length(n,i)*C_length(n,j))
            END DO
        END DO
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE NH BOND VECTORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_NH_VECTORS(bond, length, N, H, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: N(1:nconf, 1:nres, 1:3)
    DOUBLE PRECISION,INTENT(in) :: H(1:nconf, 1:nres, 1:3)
    DOUBLE PRECISION,INTENT(inout) :: bond(1:nconf, 1:nres-1, 1:3)
    DOUBLE PRECISION,INTENT(inout) :: length(1:nconf, 1:nres-1)

    INTEGER(kind=8) :: i,nn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO nn=1,nconf
        DO i=1,nres-1
            bond(nn,i,1) = H(nn,i+1,1)-N(nn,i+1,1)
            bond(nn,i,2) = H(nn,i+1,2)-N(nn,i+1,2)
            bond(nn,i,3) = H(nn,i+1,3)-N(nn,i+1,3)
            length(nn,i) = (bond(nn,i,1)**2+bond(nn,i,2)**2+bond(nn,i,3)**2)**0.5
        END DO
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE T1, T2, & NOE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_NMR_OBSERVABLES(T1,T2,NOE,P2,time_in,NHfactor,timescale,nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nres, timescale
    DOUBLE PRECISION,INTENT(in) :: NHfactor
    DOUBLE PRECISION,INTENT(in) :: P2(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION,INTENT(in) :: time_in(0:1000*timescale)
    DOUBLE PRECISION,INTENT(inout) :: T1(1:nres-1), T2(1:nres-1), NOE(1:nres-1)

    INTEGER(kind=8) :: n,t
    DOUBLE PRECISION :: c2, di, d2, rnhin3, w0, w1, w2, w3, w4
    DOUBLE PRECISION :: uo, hp, wh, wn, gh, gn, dn, pi
    DOUBLE PRECISION :: time(0:1000*timescale), dt(1:1000*timescale)
    DOUBLE PRECISION :: fj0(1:nres-1), fint0(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fints0(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fj1(1:nres-1), fint1(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fints1(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fj2(1:nres-1), fint2(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fints2(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fj3(1:nres-1), fint3(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fints3(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fj4(1:nres-1), fint4(0:1000*timescale,1:nres-1)
    DOUBLE PRECISION :: fints4(0:1000*timescale,1:nres-1)

! ORIGINAL COMMENT (???):
!   All the parameters below (save for pi) should be changeable. How could we do that?
! UPDATE (pgromano):
!   Super easy to set this up with F2PY... but what are these values?
!   I will eventually come back and rewrite everything with more intuitive names
!   or at the very least add comments to describe what is being done.

    uo = 1.256637D-06
    hp = 6.62608D-34
    wh = 599.98D+06
    wn = 60.8D+06
    gh = 26.7519D+07
    gn = -2.7126D+07
    dn = 160D-06
    pi = 3.141592654d0
    rnhin3 = 9.90688385D+29

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Convert time into picoseconds
    rnhin3=rnhin3*NHfactor
    DO t=0,1000*timescale
        time(t) = time_in(t)*1D-12
        IF(t.gt.0)THEN
            dt(t) = time(t)-time(t-1)
        END IF
    END DO

    w0 = 0.d0
    w1 = wn
    w2 = wh
    w3 = wh+wn
    w4 = wh-wn

    DO n=1,nres-1
        fint0(0,n) = 1.d0
        fints0(0,n) = 0.d0
        fj0(n) = 0.d0
        DO t=1,1000*timescale
            fint0(t,n) = P2(t,n)*DCOS(w0*time(t)*2*pi)
            fints0(t,n) = fints0(t-1,n) + &
                        & 0.5*dt(t)*(fint0(t,n)+fint0(t-1,n))
        END DO
        fj0(n) = 0.4*fints0(1000*timescale,n)
    END DO

    DO n=1,nres-1
        fint1(0,n) = 1.d0
        fints1(0,n) = 0.d0
        fj1(n) = 0.d0
        DO t=1,1000*timescale
            fint1(t,n) = P2(t,n)*DCOS(w1*time(t)*2*pi)
            fints1(t,n) = fints1(t-1,n) + &
                        & 0.5*dt(t)*(fint1(t,n)+fint1(t-1,n))
        END DO
        fj1(n) = 0.4*fints1(1000*timescale,n)
    END DO

    DO n=1,nres-1
        fint2(0,n) = 1.d0
        fints2(0,n) = 0.d0
        fj2(n) = 0.d0
        DO t=1,1000*timescale
            fint2(t,n) = P2(t,n)*DCOS(w2*time(t)*2*pi)
            fints2(t,n) = fints2(t-1,n) + &
                        & 0.5*dt(t)*(fint2(t,n)+fint2(t-1,n))
        END DO
        fj2(n) = 0.4*fints2(1000*timescale,n)
    END DO

    DO n=1,nres-1
        fint3(0,n) = 1.d0
        fints3(0,n) = 0.d0
        fj3(n) = 0.d0
        DO t=1,1000*timescale
            fint3(t,n) = P2(t,n)*DCOS(w3*time(t)*2*pi)
            fints3(t,n) = fints3(t-1,n) + &
                        & 0.5*dt(t)*(fint3(t,n)+fint3(t-1,n))
        END DO
        fj3(n) = 0.4*fints3(1000*timescale,n)
    END DO

    DO n=1,nres-1
        fint4(0,n) = 1.d0
        fints4(0,n) = 0.d0
        fj4(n) = 0.d0
        DO t=1,1000*timescale
            fint4(t,n) = P2(t,n)*DCOS(w4*time(t)*2*pi)
            fints4(t,n) = fints4(t-1,n) + &
                        & 0.5*dt(t)*(fint4(t,n)+fint4(t-1,n))
        END DO
        fj4(n) = 0.4*fints4(1000*timescale,n)
    END DO

    di = uo*uo*hp*hp*gh*gh*gn*gn                ! (uo*hp*gh*gn)**2
    d2 = di*rnhin3*rnhin3/(64*pi*pi*pi*pi)      ! (di*(rnhin3**2))/(8*(pi**2))**2
    c2 = 4*pi*pi*wn*wn*dn*dn/3                  ! ((2*pi*wn*dn)**2)/3

    DO n=1,nres-1
        ! Calculate T1 per residue
        T1(n) = (c2*fj1(n)+0.25*d2*(3*fj1(n)+6*fj3(n)+fj4(n)))**(-1)

        ! Calculate T2 per residue
        T2(n) = ((c2*(4*fj0(n)+3*fj1(n))/6) + &
                & d2*(4*fj0(n)+3*fj1(n)+6*fj2(n)+6*fj3(n)+fj4(n))/8)**(-1)

        ! Calculate NOE per residue
        NOE(n) = 1.d0+0.25*d2*(gh/gn)*(6*fj3(n)-fj4(n))*T1(n)
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE Q MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_NU_EIGENVALUES(w_mu, QI, U, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nres
    DOUBLE PRECISION,INTENT(in) :: QI(1:nres-1,1:nres-1), U(1:nres-1,1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: w_mu(1:nres-1)

    INTEGER(kind=8) :: i,j,k,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate mu eigenvalues and mode lengths
    DO k=1,nres-1
        DO i=1,nres-1
            DO j=1,nres-1
                w_mu(k) = w_mu(k)+QI(k,i)*U(i,j)*QI(k,j)
            END DO
        END DO
    END DO

! Invert mu eigenvalues
    DO n=1,nres-1
        w_mu(n) = 1.0/w_mu(n)
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE M1 BOND AUTOCORRELATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_P2(P2, time, tau, tau_m1, barriers, modelength, lam, Q, QI,&
    &   NH, blsq, mu, sigma, temp, nres, timescale)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nres, temp, timescale
    DOUBLE PRECISION,INTENT(in) :: mu(1:nres-1), sigma, blsq
    DOUBLE PRECISION,INTENT(in) :: NH(1:nres-1, 1:nres-1)
    DOUBLE PRECISION,INTENT(in) :: Q(1:nres-1,1:nres-1), QI(1:nres-1,1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: tau(1:nres-1), tau_m1(1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: barriers(1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: modelength(1:nres-1,1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: lam(1:nres-1)
    DOUBLE PRECISION,INTENT(inout) :: time(1:1000*timescale+1)
    DOUBLE PRECISION,INTENT(inout) :: P2(1:1000*timescale+1,1:nres-1)

    INTEGER(kind=8) :: i, j, k, n, t, io, index, t_final, ti, tf, dt
    DOUBLE PRECISION :: AMPsum(1000*timescale+1,1:nres), modeAmp(1:nres-1,1:nres-1)
    DOUBLE PRECISION :: w1, w2, w3, X2, ES, t_red(1:1000*timescale+1), fricorr(1:nres-1)
    DOUBLE PRECISION :: norm(1:nres-1), modenorm(1:nres-1)
    DOUBLE PRECISION :: tau_run_avg(1:1000*timescale+2,1:nres-1)
    DOUBLE PRECISION :: tau_m1_run_avg(1:1000*timescale+2,1:nres-1)
    DOUBLE PRECISION :: Rb, eps,  NHAmp, pi

    pi = 3.141592654
    Rb = 0.00198
    eps = 6.42
    NHAmp = 0.02
    t_final = 1000*timescale+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculate the Mode fluctuations
    DO i=1,nres-1
        DO j=1,nres-1
            modelength(j,i) = (((blsq*(Q(j,i)**2)/mu(i))**0.5)*10) !(Angstroms)
        END DO
    END DO

! Calculate barriers in modes from "GLASSY MYSTERY FUNCTION"
    barriers = 0.d0
    fricorr = 0.d0
    DO i=4,nres-1
        barriers(i) = ((blsq/mu(i))**0.5)*eps
        fricorr(i) = DEXP(barriers(i)/(Rb*temp))
    END DO

! Calculate the Mode Amplitudes (modeAmp) in M1 relaxation
! First 3 modes use NH orientations expressed from CA basis
    modeAmp = 0.d0
    DO i=1,nres-1
        DO j=1,3 ! MODE LOOP
            DO k=1,nres-1 ! Dummy loop over CA vector basis set
                modeAmp(i,j) = modeAmp(i,j)+QI(j,k)*NH(i,k)
            END DO
            modeAmp(i,j) = (modeAmp(i,j)**2)*mu(j)
        END DO
    END DO

! Internal modes
    DO i=1,nres-1
        DO j=4,nres-1
            modeAmp(i,j) = (Q(i,j)**2)/mu(j)
        END DO
    END DO

! Normalize the sum of the first 3 modes to the alpha carbon bond representation
! and remove the fluctuations from NH basis.
    norm = 0.d0
    modenorm = 0.d0
    DO i=1,nres-1
        DO j=1,3
            norm(i) = norm(i)+(Q(i,j)**2)/mu(j)
            modenorm(i) = modenorm(i)+modeAmp(i,j)
        END DO
        norm(i) = norm(i)-NHAmp
        DO j=1,3
            modeAmp(i,j) = modeAmp(i,j)*(norm(i)/modenorm(i))
        END DO
    END DO

! Add rigid body eigen combinations
    w1 = 0.5*(lam(1)+lam(2))
    w2 = 0.5*(lam(1)+lam(3))
    w3 = 0.5*(lam(2)+lam(3))
    lam(1) = w1
    lam(2) = w2
    lam(3) = w3

! Here since we are done with the dummy norm array, we will use the first entry
! as a different normalization term.
!   norm(1) : normalize mode amplitude of CA basis
    DO i=1,nres-1
        norm(1) = 0.d0
        DO j=1,nres-1
            norm(1) = norm(1)+modeAmp(i,j)
        END DO
        norm(1)=norm(1)+NHAmp

        DO j=1,nres-1
            modeAmp(i,j) = modeAmp(i,j)/norm(1)
        END DO
    END DO

! Build reduced time spacings to inteligently build statistics by sampling a lot
! at short times and sample less at long times. Here we set up to sample given
! timescale.

    index = 1
    t_red(1) = 0.d0
    DO t=1,timescale
        ti = INT(10.d0**(t-1))
        tf = INT(10.d0**(t+2))
        dt = INT(10.d0**(t-1))

        IF(t.ne.timescale)THEN
            DO n=ti,tf,dt
                index = index+1
                t_red(index+1) = t_red(index)+FLOAT(dt)*0.001d0
            END DO
        ELSE IF(t.eq.timescale)THEN
            DO n=ti,tf-1,dt
                index = index+1
                t_red(index+1) = t_red(index)+FLOAT(dt)*0.001d0
            END DO
        END IF
    END DO

! Calculate the M1 bond auto-correlation. This can be found in equation 2.7 from
! Jeremy Coppermans's thesis.
    AMPsum = 0.d0
    tau_run_avg = 0.d0
    tau_m1_run_avg = 0.d0
    DO n=1,nres-1                       ! Residue Loop
        AMPsum(1,n) = 0.d0
        DO t=1,t_final                  ! Time Loop
            DO k=1,nres-1               ! Mode Loop
                IF(k.le.3)THEN
                    ES = DEXP(-lam(k)*t_red(t))
                ELSE
                    ES = DEXP(-(lam(k)*t_red(t))/fricorr(k))
                END IF
                AMPsum(t,n) = AMPsum(t,n)+modeAmp(n,k)*ES
            END DO
            AMPsum(t,n) = AMPsum(t,n)+NHAmp*DEXP(-t_red(t)/(sigma*0.2))

! Change of variables such that the second order Legendre polynomial of the time
! dependent bond orientation function (P2) can be expressed as a function of M1.
! This can be found in equation 2.8 from Jeremy Coppermans's thesis.
!       X2 = x**2 = (1-M1**2)/(M1**2)
!       P2 = 1-3(x**2-(2/pi) * x**3 * (1-(2/pi) * atan(x)))

            P2(t,n) = 0.d0
            IF(AMPsum(t,n).gt.0.255d0)THEN
                IF(AMPsum(t,n).gt.0.99999)THEN
                    X2 = 0.d0
                ELSE
                    X2 = (1.d0-(AMPsum(t,n)**2))/(AMPsum(t,n)**2)
                END IF

                P2(t,n) = 1.d0-(3.d0*X2)+(3.d0*(DSQRT(X2)**3) * &
                        & (1.d0/(2.d0/pi)-DATAN(DSQRT(X2))))

                IF(P2(t,n).le.10E-9)THEN
                    P2(t,n) = 1E-9
                END IF
            ELSE
                P2(t,n) = ((3.d0/5.d0)*(AMPsum(t,n)**2)) + &
                        & ((6.d0/35.d0)*(AMPsum(t,n)**2)**2) + &
                        & ((8.d0/105.d0)*(AMPsum(t,n)**2)**3) + &
                        & ((16.d0/385.d0)*(AMPsum(t,n)**2)**4)
            END IF

! Rescale reduced times to real time by sigma factor
            time(t) = t_red(t)/sigma

            IF(t.gt.1)THEN
! Calculate bond autocorrelation relaxation time
                tau_run_avg(t+1,n) = 0.5 * (t_red(t)-t_red(t-1)) * &
                                    &(P2(t,n)+P2(t-1,n)) + tau_run_avg(t,n)

! Calculate M1 bond autocorrelation relaxation time
                tau_m1_run_avg(t+1,n) = 0.5 * (t_red(t)-t_red(t-1)) * &
                                    & (AMPsum(t,n)+AMPsum(t-1,n)) + &
                                    & tau_m1_run_avg(t,n)
            END IF

            IF(P2(t,n).le.1E-8)THEN
                P2(t,n) = 1E-8
            END IF
        END DO
        tau(n) = tau_run_avg(t_final,n)/sigma
        tau_m1(n) = tau_m1_run_avg(t_final,n)/sigma
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE R MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_R_MATRIX(R, C, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: C(1:nconf, 1:nres, 1:3)
    DOUBLE PRECISION,INTENT(inout) :: R(1:nconf, 1:nres, 1:nres)

    INTEGER(kind=8) :: i,j,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO n=1,nconf
        DO i=1,nres
            DO j=1,nres
                IF(i.ne.j)THEN
                    R(n,i,j) = ((C(n,i,1)-C(n,j,1))**2+ &
                               &(C(n,i,2)-C(n,j,2))**2+ &
                               &(C(n,i,3)-C(n,j,3))**2)**0.5
                    R(n,i,j) = 1.0/R(n,i,j)
                ELSE IF(i.eq.j)THEN
                    R(n,i,j) = 0.0
                END IF
            END DO
        END DO
    END DO
    END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE U MATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CALCULATE_U_MATRIX(U, bond, MSF, temp, nconf, nres)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    INTEGER(kind=8),INTENT(in) :: nconf, nres
    DOUBLE PRECISION,INTENT(in) :: temp, bond(1:nconf, 1:nres-1, 1:3)
    DOUBLE PRECISION,INTENT(in) :: MSF(1:nconf, 1:nres, 1:nres)
    DOUBLE PRECISION,INTENT(inout) :: U(1:nconf, 1:nres-1, 1:nres-1)

    INTEGER(kind=8) :: i,j,n
    DOUBLE PRECISION :: Rb, c

    Rb = 0.001987204
    c = (3.0*Rb*temp)/6.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO n=1,nconf
        DO i=1,nres-1
            DO j=1,nres-1
                U(n,i,j) = bond(n,i,1)*bond(n,j,1) + &
                        & bond(n,i,2)*bond(n,j,2) + &
                        & bond(n,i,3)*bond(n,j,3) + &
                        & c*(MSF(n,i,j) + MSF(n,i+1,j+1) - &
                        & MSF(n,i+1,j) - MSF(n,i,j+1))
            END DO
        END DO
    END DO
    END
