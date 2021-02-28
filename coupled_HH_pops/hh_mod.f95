!! This Program is the module used in the simulation of following paper
!! 
!!Transmission delays and frequency detuning can regulate information flow between brain regions
!! Aref Pariz, Ingo Fischer , Alireza Valizadeh, Claudio Mirasso
!! PLOS Comp. Bio. 2021

!! Compiled on 21 Feb. 2021 by gfortran 10.2.0 on Ubuntu 20.10 as following lines

!! gfortran -m64 -O3 hh_mod.f95 main_hh.f95 -o a1.out
!! echo $deeout $dI | ./a1.out

!! Written by Aref Pariz. Email: pariz.aref@gmail.com

MODULE hh_mod
    IMPLICIT NONE
    !-----------------
    INTEGER, PARAMETER 		:: DP = KIND(4.0D0)
    INTEGER 			:: neurons = 100, numnet = 2
    INTEGER                     :: r = 2000 ! Poissonian rate
    !INTEGER, DIMENSION(:)       :: selected_neurons
    REAL(KIND = DP) :: PI = 4*ATAN(1.d0)
    REAL(KIND = DP) :: omega = 3, amp1 = 0.1, amp2=0.1, selected_neurons_percent = 0.5
    REAL(KIND = DP) :: deeout, deiout !deeout=1, deiout=1
    REAL(KIND = DP), PARAMETER 	:: dt = 0.01
	REAL(KIND = DP)         	:: xg
    REAL(KIND = DP), PARAMETER 	:: gL_bar = 0.3, gK_bar = 36, gNa_bar = 120
    REAL(KIND = DP), PARAMETER 	:: E_L = -54.387, E_K = -77, E_Na = 50
    REAL(KIND = DP), PARAMETER 	:: Pee1 = 0.1, Pei1 = 0.1, Pie1 = 0.1, Pii1 = 0.1, Peeout = 0.05, PeiNet=0.
    REAL(KIND = DP), PARAMETER 	:: Pee2 = 0.1, Pei2 = 0.1, Pie2 = 0.1, Pii2 = 0.1
    REAL(KIND = DP), PARAMETER 	:: dee = 1, dei = 1, die = 1, dii = 1
    REAL(KIND = DP), PARAMETER 	:: sec = 6    !REAL(KIND = DP) :: gee1, gee2, gei1, gei2, gie1, gie2, gii1, gii2, 
    REAL(KIND = DP), PARAMETER  :: geeout = 0.0025, geiout= 0
    REAL(KIND = DP), PARAMETER 	:: gee1 = 0.0025, gei1 = 0.005, gie1 = 0.01, gii1 = 0.01
    REAL(KIND = DP), PARAMETER 	:: gee2 = 0.0025, gei2 = 0.005, gie2 = 0.01, gii2 = 0.01
    REAL(KIND = DP), PARAMETER 	:: taud = 3.0_DP, taur = 0.5_DP
    
!    REAL(KIND = DP)            	:: gii
    
CONTAINS
    SUBROUTINE initialize_seed()
        INTEGER :: i
        INTEGER, DIMENSION(8) :: dtt
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(SIZE = i)
        ALLOCATE(seed(i))
        CALL DATE_AND_TIME(VALUES = dtt) !dtt = 1
        CALL RANDOM_SEED(GET = seed)
        seed(i) = dtt(8)
        seed(1) = dtt(8) * dtt(7) * dtt(6)
        CALL RANDOM_SEED(PUT = seed)
    END SUBROUTINE initialize_seed
    FUNCTION funv(vv, nn, mm, hh, synI, extI)
        IMPLICIT NONE
        REAL(KIND = DP), PARAMETER :: gL_bar = 0.3, gK_bar = 36, gNa_bar = 120
        REAL(KIND = DP), PARAMETER :: E_L = -54.387, E_K = -77, E_Na = 50
        REAL(KIND = DP), INTENT(IN) :: vv, nn, mm, hh, synI, extI
        REAL(KIND = DP) :: funv
        funv = extI + synI  - gK_bar * (nn * nn * nn * nn)*(vv - E_K) &
        -gNa_bar * (mm * mm * mm) * hh * (vv - E_Na) - gL_bar * (vv - E_L)
    END FUNCTION funv

    FUNCTION funn(vv, nn)
        IMPLICIT NONE
        REAL(KIND = DP), INTENT(IN) :: vv, nn
        REAL(KIND = DP) :: alpha_n, beta_n, funn
        alpha_n = 0.01_DP * (vv + 55.0_DP)/(1.0_DP - exp(-0.1_DP * (vv + 55.0_DP)))
        beta_n = 0.125_DP * exp(-0.0125_DP * (vv + 65.0_DP))
        funn = alpha_n * (1.0_DP - nn) - beta_n * nn
    END FUNCTION funn

    FUNCTION funm(vv, mm)
        IMPLICIT NONE
        REAL(KIND = DP), INTENT(IN) :: vv, mm
        REAL(KIND = DP) :: alpha_m, beta_m, funm
        alpha_m = 0.1_DP * ((vv + 40.0_DP)/(1.0_DP - exp(-0.1_DP * (vv + 40.0_DP))))
        beta_m = 4.0_DP * exp(-0.0556_DP * (vv + 65.0_DP))
        funm = alpha_m * (1.0_DP - mm) - beta_m * mm
    END FUNCTION funm

    FUNCTION funh(vv, hh)
        IMPLICIT NONE
        REAL(KIND = DP), INTENT(IN) :: vv, hh
        REAL(KIND = DP) :: alpha_h, beta_h, funh
        alpha_h = 0.07_DP * exp(-0.05_DP * (vv + 65.0_DP))
        beta_h = 1.0_DP/(1.0_DP + exp(-0.1_DP * (vv + 35_DP)))
        funh = alpha_h * (1.0_DP - hh) - beta_h * hh
    END FUNCTION funh

    FUNCTION intrnd(a, b, c) result(intrnd_r)
        integer, parameter :: pesrandnum = 1000
        integer :: a, b, c
        real, dimension(:) :: x(pesrandnum)
        integer, dimension(:) :: y(pesrandnum)
        integer, dimension(:) :: intrnd_r(c)
        integer, dimension(:,:) :: pz(pesrandnum)
        integer :: ii, jj
        call random_number(x)
        y = a + int((b - a) * x)
        pz(1:pesrandnum) = 0
        DO ii = 1, pesrandnum
            DO jj = ii, pesrandnum
                IF (y(ii) .EQ. y(jj) .AND. ii .ne. jj) THEN
                    pz(ii) = 1
                END IF
            END DO
        END DO
        jj = 0
        DO ii = 1, pesrandnum
            IF (pz(ii) .EQ. 0 .AND. jj .lt. c) THEN
                jj = jj + 1
                intrnd_r(jj) = y(ii)
            ENDIF
        END DO
    END function intrnd

    FUNCTION creatA(N, Ne, Ni, numnet, Pee1, Pei1, Pie1, Pii1, Pee2, Pei2, Pie2, Pii2, PeiNet, Peeout, con)
        IMPLICIT NONE
        INTEGER :: nPee, nPei, nPie, nPii, nPeiNet
        INTEGER, DIMENSION(:,:) :: con
        !        REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: xpee, xpei, xpie, xpii
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ypee, ypei, ypie, ypii
        INTEGER :: Ne, Ni, N, numnet
        REAL(KIND = DP) :: Pee1, Pee2, Pei1, Pei2, Pie1, Pie2, Pii1, Pii2, Peeout, PeiNet
        INTEGER :: ii, jj, c1, c2, st1, en1, st2, en2, bb
        INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:) :: creatA
        ALLOCATE(creatA(N * numnet, N * numnet))
        creatA(1:numnet * N, 1:numnet * N) = 0

        
!        ALLOCATE(xpee(nPee), xpei(nPei), xpie(nPie), xpii(nPii))
!        ALLOCATE(ypee(nPee), ypei(nPei), ypie(nPie), ypii(nPii))
        DO c1 = 1, numnet ! Post, Rwos
            DO c2 = 1, numnet ! Pre, Columns
                IF (con(c1, c2) .EQ. 1) THEN
                    IF (c1 .EQ. c2) THEN
                        IF  (c1 .EQ. 1) THEN
                            nPee = INT(Ne * Pee1)
                            nPei = INT(Ne * Pei1)
                            nPie = INT(Ni * Pie1)
                            nPii = INT(Ni * Pii1)
                            bb = 1
                        ELSEIF (c1 .EQ. 2) THEN
                            nPee = INT(Ne * Pee2)
                            nPei = INT(Ne * Pei2)
                            nPie = INT(Ni * Pie2)
                            nPii = INT(Ni * Pii2)
                            bb = 1
                        END IF
                    ELSE
                        nPee = INT(Ne * Peeout)
                        nPeiNet = INT(Ne * PeiNet)
                        bb = 2
                    END IF
                    !! E -> E
                    st1 = (c1 - 1) * Ne + 1
                    en1 = c1 * Ne
                    st2 = (c2 - 1) * Ne + 1
                    en2 = c2 * Ne
                    DO ii = st1, en1 ! Rows
                        ypee = intrnd(st2, en2, nPee)
                        DO jj = 1, size(ypee) ! Columns
                            creatA(ii, ypee(jj)) = bb
                        END DO
                    END DO
                    !! E -> I
                    IF (c1 .EQ. c2) THEN
                        st1 = numnet * Ne + (c1 - 1) * Ni + 1
                        en1 = numnet * Ne + c1 * Ni
                        st2 = (c2 - 1) * Ne + 1
                        en2 = c2 * Ne
                        DO ii = st1, en1 ! Rows
                            ypei = intrnd(st2, en2, nPei)
                            DO jj = 1, size(ypei) ! Columns
                                creatA(ii, ypei(jj)) = bb
                            END DO
                        END DO

                        !! I -> E
                        st1 = (c1 - 1) * Ne + 1
                        en1 = c1 * Ne
                        st2 = numnet * Ne + (c2 - 1) * Ni + 1
                        en2 = numnet * Ne + c2 * Ni
                        DO ii = st1, en1 ! Rows
                            ypie = intrnd(st2, en2, nPie)
                            DO jj = 1, size(ypie)
                                creatA(ii, ypie(jj)) = 1
                            END DO
                        END DO
                        !! I -> I
                        st1 = numnet * Ne + (c1 - 1) * Ni + 1
                        en1 = numnet * Ne + c1 * Ni
                        st2 = numnet * Ne + (c2 - 1) * Ni + 1
                        en2 = numnet * Ne + c2 * Ni
                        DO ii = st1, en1
                            ypii = intrnd(st2, en2, nPii)
                            DO jj = 1, size(ypii)
                                creatA(ii, ypii(jj)) = 1
                            END DO
                        END DO
                    ELSEIF (c1 .NE. c2) THEN
                        ! E1&2 -> I2&1
                        st1 = (c1 - 1) * Ne + 1
                        en1 = c1 * Ne
                        st2 = numnet * Ne + (c2 - 1) * Ni + 1
                        en2 = st2 + Ni-1
                        DO ii = st2, en2 ! Rows
                            ypei = intrnd(st1, en1, nPeiNet)
                           
                            DO jj = 1, size(ypei) ! Columns
                                creatA(ii, ypei(jj)) = bb
                                !print*, "C1=",c1," c2=",c2, "nPeiNet=",nPeiNet
                            END DO
                        END DO
                    END IF
                END IF
            END DO
        END DO
        DO ii=1, numnet * neurons
            creatA(ii,ii) = 0
        END DO
    END FUNCTION creatA

    FUNCTION create_delay(dee, dei, die, dii, deeout, deiout, neurons, Ne, Ni, numnet, A, dt)
        implicit none
        INTEGER :: ii, jj, Ne, Ni, numnet, neurons
        REAL(KIND = DP) :: dee, dei, die, dii, deeout, deiout, dt
        INTEGER(KIND = 1), ALLOCATABLE, DIMENSION(:,:) :: A
        INTEGER(KIND = DP), DIMENSION(:,:) :: create_delay(numnet * neurons, numnet * neurons)
        ! E <-> E
        create_delay(1:Ne,1:Ne)=dee/dt;
        create_delay(Ne+1:numnet*Ne,Ne+1:numnet*Ne)=dee/dt;
        create_delay(1:Ne,Ne+1:numnet*Ne)=deeout/dt;
        create_delay(Ne+1:numnet*Ne,1:Ne)=deeout/dt;
        
        ! E -> I
        create_delay(numnet*Ne+1:numnet*Ne+Ni,1:Ne)=dei/dt;
        create_delay(numnet*Ne+Ni+1:numnet*neurons,Ne+1:numnet*Ne)=dei/dt;
        create_delay(numnet*Ne+Ni+1:numnet*neurons,1:Ne)=deiout/dt;
        create_delay(numnet*Ne+1:numnet*Ne+Ni,Ne+1:numnet*Ne)=deiout/dt;
        
        ! I -> E
        create_delay(1:Ne,numnet*Ne+1:numnet*Ne+Ni)=die/dt;
        create_delay(Ne+1:numnet*Ne,numnet*Ne+Ni+1:numnet*neurons)=die/dt;
        
        ! I -> I
        create_delay(numnet*Ne+1:numnet*Ne+Ni,numnet*Ne+1:numnet*Ne+Ni)=dii/dt;
        create_delay(numnet*Ne+Ni+1:numnet*neurons,numnet*Ne+Ni+1:numnet*neurons)=dii/dt;
        
!        create_delay = 0
!        DO ii = 1, numnet * Ne ! Rows, Pre
!            DO jj = 1, numnet * Ne ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    create_delay(ii, jj) = int(dee/dt)
!                END IF
!            END DO
!        END DO
!        DO ii = numnet * Ne + 1, numnet * neurons ! Rows, Pre
!            DO jj = 1, numnet * Ne ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    create_delay(ii, jj) = int(dei/dt)
!                END IF
!            END DO
!        END DO
!        DO ii = 1, numnet * Ne ! Rows, Pre
!            DO jj = numnet * Ne + 1, numnet * neurons ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    create_delay(ii, jj) = int(die/dt)
!                END IF
!            END DO
!        END DO
!        DO ii = numnet * Ne + 1, numnet * neurons ! Rows, Pre
!            DO jj = numnet * Ne + 1, numnet * neurons ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    create_delay(ii, jj) = int(dii/dt)
!                END IF
!            END DO
!        END DO
    END function create_delay

    FUNCTION synapses(gee1, gei1, gie1, gii1, gee2, gei2, gie2, gii2, geeout, geiout, neurons, Ne, Ni, numnet, A)
        IMPLICIT NONE
        INTEGER :: ii, jj, Ne, Ni, numnet, neurons
        REAL(KIND = DP) :: gee1, gei1, gie1, gii1, gee2, gei2, gie2, gii2, geeout, geiout
        INTEGER(KIND = 1), DIMENSION(:,:) :: A
        REAL(KIND = DP), ALLOCATABLE, DIMENSION (:,:) :: synapses
        ALLOCATE(synapses(numnet * neurons, numnet * neurons))
	synapses = 0
        ! E <-> E
        synapses(1:Ne,1:Ne)=gee1;
        synapses(Ne+1:numnet*Ne,Ne+1:numnet*Ne)=gee2;
        synapses(1:Ne,Ne+1:numnet*Ne)=geeout;
        synapses(Ne+1:numnet*Ne,1:Ne)=geeout;
        
        ! E -> I
        synapses(numnet*Ne+1:numnet*Ne+Ni,1:Ne)=gei1;
        synapses(numnet*Ne+Ni+1:numnet*neurons,Ne+1:numnet*Ne)=gei2;
        synapses(numnet*Ne+Ni+1:numnet*neurons,1:Ne)=geiout;
        synapses(numnet*Ne+1:numnet*Ne+Ni,Ne+1:numnet*Ne)=geiout;
        
        ! I -> E
        synapses(1:Ne,numnet*Ne+1:numnet*Ne+Ni)=gie1;
        synapses(Ne+1:numnet*Ne,numnet*Ne+Ni+1:numnet*neurons)=gie2;
        
        ! I -> I
        synapses(numnet*Ne+1:numnet*Ne+Ni,numnet*Ne+1:numnet*Ne+Ni)=gii1;
        synapses(numnet*Ne+Ni+1:numnet*neurons,numnet*Ne+Ni+1:numnet*neurons)=gii2;
!        DO ii = 1, numnet * Ne ! Rows, Pre
!            DO jj = 1, numnet * Ne ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    synapses(ii, jj) = gee
!                END IF
!            END DO
!        END DO
!        DO ii = numnet * Ne + 1, numnet * neurons ! Rows, Pre
!            DO jj = 1, numnet * Ne ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    synapses(ii, jj) = gei
!                END IF
!            END DO
!        END DO
!        DO ii = 1, numnet * Ne ! Rows, Pre
!            DO jj = numnet * Ne + 1, numnet * neurons ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    synapses(ii, jj) = gie
!                END IF
!            END DO
!        END DO
!        DO ii = numnet * Ne + 1, numnet * neurons ! Rows, Pre
!            DO jj = numnet * Ne + 1, numnet * neurons ! Columns, Post
!                IF (A(ii, jj) .EQ. 1) THEN
!                    synapses(ii, jj) = gii
!                END IF
!            END DO
!        END DO
    END function synapses


    FUNCTION Poissonsp(numnet, neurons, time, r, taud, dt)
        IMPLICIT NONE
        INTEGER :: ii, jj, numnet, neurons, time, r
        INTEGER(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: Xpoissonsp
        REAl(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: Poissonsp
        REAL(KIND = DP) :: dt, x, taud
        ALLOCATE(Poissonsp(numnet * neurons, time), Xpoissonsp(numnet * neurons, time))
        Poissonsp = 0
        DO ii = 1, time
            DO jj = 1, numnet * neurons
                CALL random_number(x)
                IF (x .lt. r * dt / 1000) THEN
                    Xpoissonsp(jj, ii) = 1
                ELSE
                    Xpoissonsp(jj, ii) = 0
                END IF
            END DO
        END DO
        DO jj = 1, numnet * neurons
            DO ii = 1, time
                IF (Xpoissonsp(jj, ii) .EQ. 0) THEN
                    Poissonsp(jj, ii + 1) = Poissonsp(jj, ii) - dt * Poissonsp(jj, ii) / taud
                ELSE
                    Poissonsp(jj, ii + 1) = Poissonsp(jj, ii) + 1
                END IF
            END DO
        END DO
    END FUNCTION Poissonsp
END MODULE hh_mod
