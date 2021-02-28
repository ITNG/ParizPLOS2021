!! This Program is the Main program used in the simulation of following paper
!! 
!!Transmission delays and frequency detuning can regulate information flow between brain regions
!! Aref Pariz, Ingo Fischer , Alireza Valizadeh, Claudio Mirasso
!! PLOS Comp. Bio. 2021

!! Compiled on 21 Feb. 2021 by gfortran 10.2.0 on Ubuntu 20.10 as following lines

!! gfortran -m64 -O3 hh_mod.f95 main_hh.f95 -o a1.out
!! echo $deeout $dI | ./a1.out

!! Written by Aref Pariz. Email: pariz.aref@gmail.com


PROGRAM hh_double

    USE hh_mod

    IMPLICIT NONE
    INTEGER 	:: Eneuron, Ineuron
    INTEGER 	:: ii, jj, ll, kk, ierr
    INTEGER	:: time, f1, f2, ensemble=10, ens
    INTEGER(KIND = DP)	:: delay0
    INTEGER, DIMENSION(:), ALLOCATABLE              :: selected_neurons
    INTEGER,            DIMENSION(:,:), ALLOCATABLE :: con
    INTEGER(KIND = DP), DIMENSION(:),   ALLOCATABLE :: rhosum
    INTEGER(KIND = DP), DIMENSION(:),   ALLOCATABLE :: st1, st2, st3
    !INTEGER(KIND = 1), 	DIMENSION(:),	ALLOCATABLE :: fired
    INTEGER(KIND = 1),  DIMENSION(:,:), ALLOCATABLE :: A
    INTEGER(KIND = 1),  DIMENSION(:,:), ALLOCATABLE :: rho
    INTEGER(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: delay
    INTEGER(KIND = 1) :: xdeeout
    INTEGER(KIND = 1) :: dI
    INTEGER(KIND = DP), DIMENSION(:),   ALLOCATABLE :: dIneuron
    
    REAL :: start, finish
    REAL, PARAMETER  :: NeRate = 0.8
    REAL(KIND = DP) :: xv, xn, xm, xh
    REAL(KIND = DP) :: xsyn, xIext, xnoise
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: vth
    REAL(KIND = DP), PARAMETER  	:: Esynout = 0
    REAL(KIND = DP), PARAMETER 		:: geeIext = 0.004
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: C, s0, S, Isyn, S1, S2, S3
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: m, n, h, v
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE :: taurise,taudecay 
    REAL(KIND = DP), DIMENSION(:), ALLOCATABLE	:: Esyn
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: synapse
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: SS
    !REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Vsave
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: vrnd
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: xIsyn
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: IextSave
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: Isignal
    REAL(KIND = DP), ALLOCATABLE, DIMENSION(:) :: I0
    !REAL(KIND = DP), ALLOCATABLE, DIMENSION(:,:) :: Inoise
    CHARACTER(len=1024)         					:: filename1,filename2
    
    REAL(KIND = DP), DIMENSION(:,:), ALLOCATABLE :: Iext  
    time = INT(1000 * sec/dt)
	
    ALLOCATE(synapse(numnet * neurons, numnet * neurons), &
    !        fired(numnet * neurons), &
            taudecay(numnet*neurons), &
            taurise(numnet*neurons), &
            C(numnet*neurons), &
            I0(numnet*neurons), &
            dIneuron(numnet*neurons), &
            vth(numnet*neurons), &
            A(numnet * neurons, numnet * neurons), &
            delay(numnet * neurons, numnet * neurons), &
            SS(neurons * numnet, time), &
            Isyn(neurons * numnet), &
            Isignal(time), &
            !Inoise(numnet * neurons, time), &
            xIsyn(neurons * numnet, neurons * numnet), &
            m(neurons * numnet), &
            n(neurons * numnet), &
            h(neurons * numnet), &
            v(neurons * numnet), &
            s0(neurons * numnet), &
            Esyn(numnet * neurons), &
            S(neurons * numnet), &
            S1(neurons * numnet), &
            S2(neurons * numnet), &
            S3(neurons * numnet), &
            st1(neurons * numnet), &
            st2(neurons * numnet), &
            st3(neurons * numnet), &
            Iext(numnet * neurons, time), &
            IextSave(numnet * neurons, time), &
            !Vsave(numnet * neurons, time), &
            rho(numnet * neurons, time), &
            rhosum(time), &
            con(numnet,numnet), &
            vrnd(numnet*neurons), &
            selected_neurons(INT(selected_neurons_percent * Eneuron)), &
            STAT = ierr)
    IF (ierr /= 0) THEN
            PRINT*, "Error :: in ALLOCATION, check the DIMENSION of matrices"
            STOP
    END IF
    PRINT*, "Allocation succeed"
    READ(*,*) xdeeout, dI
    deeout = 0.0 + 0.25 * xdeeout
    deiout = 0.0
    CALL cpu_time(start)
    Eneuron = INT(NeRate * neurons)
    Ineuron = INT(neurons - Eneuron)

    !!-----------------Periodic Signal-------------START----------------

        !selected_neurons = intrnd(1, Eneuron, INT(selected_neurons_percent * Eneuron))
        Isignal(:)=0
        OPEN(unit=1,file='Isignalx.txt')
        READ(1,*) (Isignal(ii), ii=1,time)
        !DO ll=1, Eneuron
        !    Isignal(ll,:)=Isignal(1,:)
        !    PRINT*, ll
        !END DO
        CLOSE(1)
        !PRINT*, Isignal(50,:)
                
        !DO ii = 0, time-1
        !    DO jj=1, size(selected_neurons) 
        !        Isignal(selected_neurons(jj),ii) = amp1 * sin( omega * 2 * PI * ii * dt / 1000 ) &
        !        + amp2 * sin( SQRT(omega) * 2 * PI * ii * dt / 1000 )
        !    END DO
        !END DO
        !PRINT*, "Isignal Created"
        !DO ii = 1, numnet * neurons
        !WRITE(900, *) (Isignal(ii, jj), jj = 1, time)
        !END DO
        !CLOSE(900)
        !PRINT*, "Isignal has written"
        !STOP
    !!-----------------Periodic Signal-------------END------------------
        Esyn(1:numnet * Eneuron) = 0.0
        Esyn(numnet * Eneuron + 1:numnet * neurons) = -80.0

        taudecay(1:numnet*Eneuron)=3;                  ! AMPA  DECAY TIME
        taudecay(numnet*Eneuron+1:numnet*neurons)=3;   ! GABAa DECAY TIME
        taurise(1:numnet*Eneuron)=0.5;                   ! AMPA  RISE  TIME
        taurise(numnet*Eneuron+1:numnet*neurons)=0.5;    ! GABAa RISE  TIME
        con(1,1) = 1            ! Net1 -> Net1
        con(1,2) = 1            ! Net2 -> Net1
        con(2,1) = 1            ! Net1 -> Net2
        con(2,2) = 1            ! Net2 -> Net2
!!---------------------OPENING FILES SECTION------------------------

!	    OPEN(UNIT = 100, FILE  ="rho.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!	    OPEN(UNIT = 200, FILE = "A.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!	    OPEN(UNIT = 300, FILE = "Isyn.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!	    OPEN(UNIT = 400, FILE = "V.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!	    OPEN(UNIT = 500, FILE = "synaptic_weight.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!	    OPEN(UNIT = 600, FILE = "delay.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!	    OPEN(UNIT = 700, FILE = "externalc.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!       OPEN(UNIT = 800, FILE = "SS.txt", ACTION = "WRITE", STATUS = "UNKNOWN")
!       OPEN(UNIT = 900, FILE = "Isgnal.txt", ACTION = "WRITE", STATUS = "UNKNOWN")

!!---------------------End of OPENING FILES SECTION-----------------

!!-%%%%%%%%%%%%%%%%%%-MAIN PART OF PROGRAM STARS-%%%%%%%%%%%%%%%%%%%%%%%-
                !!-----------------GENERATING CONCCETION MATRIX START---------------
        A = 0
        delay = 0
        synapse = 0
        I0(1:Eneuron) = 8.5
        I0(1+Eneuron:numnet * Eneuron) = 9.5
        I0(numnet * Eneuron + 1:numnet * Eneuron + Ineuron) = 8.5
        I0(numnet * Eneuron + Ineuron + 1:numnet * neurons) = 9.5
        dIneuron = 0
        dIneuron(1:Eneuron) = 1
        dIneuron(numnet * Eneuron + 1 : numnet * Eneuron + Ineuron) = 2
        A = creatA(neurons, Eneuron, Ineuron, numnet, Pee1, Pei1, Pie1, Pii1, Pee2, Pei2, Pie2, Pii2, PeiNet, peeout, con)
        PRINT*, "Connectivity Matrix"
        !DO ii = 1, numnet * neurons
        !    WRITE(200, *) (A(ii, jj), jj = 1, numnet * neurons)
        !END DO
        !CLOSE(200)
        !PRINT*, "A has written"
        !STOP
        
!!-----------------GENERATING CONNECTION MATRIX END-----------------

!!------------------GENERATING DELAY SECTION----START---------------

        delay = create_delay(dee, dei, die, dii, deeout, deiout, neurons, Eneuron, Ineuron, numnet, A, dt)
        Do ii=1,numnet*neurons
                Do jj=1,numnet*neurons
                        IF (A(ii,jj) .EQ. 0) THEN
                                delay(ii,jj) = 0
                        END IF
                END DO
        END DO
        print*, "Delay Matrix"
        !DO ii = 1, numnet * neurons
        !    WRITE(600, *) (delay(ii, jj), jj = 1, numnet * neurons)
        !END DO
        !CLOSE(600)
        !PRINT*, "Delay has written"
!        STOP

!!------------------GENERATING DELAY SECTION----END-----------------

!!------------------SYNAPTIC WEIGHT SECTION----START----------------

        synapse = synapses(gee1, gei1, gie1, gii1, gee2, gei2, gie2, gii2, geeout, geiout, neurons, Eneuron, Ineuron, numnet, A)
        Do ii=1,numnet*neurons
                Do jj=1,numnet*neurons
                        IF (A(ii,jj) .EQ. 0) THEN
                                synapse(ii,jj)=0
                        END IF
                END DO
        END DO
        PRINT*, "Synaptic Weight Matrix"
        !DO ii = 1, numnet * neurons
        !    WRITE(500, *) (synapse(ii, jj), jj = 1, numnet * neurons)
        !END DO
        !CLOSE(500)
        !PRINT*, "synaptic weight has written"
        
!!------------------SYNAPTIC WEIGHT SECTION----END------------------


!PRINT*, 'it is ok till here 1'
!!------------------GENERATION INOISE ---------START----------------
		!DO ii = 1, time -1
        !    DO ll = 1, numnet * neurons
        !        call random_number(xnoise)
        !        Inoise(ll, ii) = 0.2*(-1+0.5*xnoise)/sqrt(dt)
        !    END DO
        !END DO
!!------------------GENERATION INOISE ---------END------------------        
DO ens = 1, ensemble
        PRINT*, "ens=", INT(ens), ", deeout=", INT(deeout), '(ms)', ", dI", dI
        rho = 0
        
        SS = 0
        S1 = 0
        S2 = 0
        !S3 = 0
        !delay0 = INT(MAX(dee, dei, die, dii, deeout, deiout)/dt)
!        PRINT*, delay0
        st1 = -100000
        st2 = -100000
        !st3 = -100000
        
        s0 = [ - 59.8977, 0.0536, 0.5925, 0.3192 ];
        call random_number(vrnd)
        v 	= s0(1) + 5 * vrnd
        m 	= s0(2)
        h 	= s0(3)
        n 	= s0(4)
        Isyn = 0.0
        DO ii=1, numnet * neurons
            call random_number(xnoise)
            vth(ii)=20.0+5.0*(-1+0.5*xnoise)
        END DO
        !!------------------GENERATING IEXT SECTION----START----------------

        Iext = Poissonsp(numnet, neurons, time, r, taud, dt)
        PRINT*, "External Current is created"
        
!        WRITE(filename1,'(A4,I3,A7)') "Iext",INT(gii*1000),"E-3.txt" 
!        OPEN(UNIT=mn+200, FILE=TRIM(filename1), ACTION = "WRITE", STATUS = "UNKNOWN")
        !DO ii = 1, numnet * neurons
        !    WRITE(700, *) (Iext(ii, jj), jj = 1, time)
        !END DO
        !CLOSE(700)
        !PRINT*, "Iext has written"

!!------------------GENERATING IEXT SECTION----END------------------



!! ----------LOOP OVER TIME AND NEURONS--------START----------------
        PRINT*, "Main part is started"
        DO ii = 1, time - 1
            DO ll = 1, numnet * neurons
                xv = v(ll)
                xn = n(ll)
                xm = m(ll)
                xh = h(ll)
                xsyn = Isyn(ll)
				xIext = -geeIext*(xv-Esynout)*Iext(ll,ii)  + I0(ll)
				
                IF (dIneuron(ll) .EQ. 1) THEN
                    xIext = xIext  + (dI * 0.05) +Isignal(ii)
                ELSEIF (dIneuron(ll) .EQ. 2) THEN
                    xIext = xIext  + (dI * 0.05)
                END IF
                IextSave(ll,ii) = xIext
                v(ll) = xv + dt * funv(xv, xn, xm, xh, xsyn, xIext)
                m(ll) = xm + dt * funm(xv, xm)
                h(ll) = xh + dt * funh(xv, xh)
                n(ll) = xn + dt * funn(xv, xn)
                IF ((xv .lt. vth(ll)) .AND. (v(ll) .gt. vth(ll))) THEN
                    st1(ll) = st2(ll)
                    !st2(ll) = st3(ll)
                    st2(ll) = ii
                    !fired(ll) = 1
                    rho(ll, ii + 1) = 1
                ELSE
                    rho(ll, ii + 1) = 0
                END IF
                S1(ll) = (EXP(-dt*(ii-st1(ll))/taur)-EXP(-dt*(ii-st1(ll))/taudecay(ll))) &
                         /((taurise(ll)/taudecay(ll))**(taurise(ll)/(taudecay(ll)-taurise(ll))) &
                         -(taurise(ll)/taudecay(ll))**(taudecay(ll)/(taudecay(ll)-taurise(ll))))
                        
                S2(ll) = (EXP(-dt*(ii-st2(ll))/taur)-EXP(-dt*(ii-st2(ll))/taudecay(ll))) &
                         /((taurise(ll)/taudecay(ll))**(taurise(ll)/(taudecay(ll)-taurise(ll))) &
                         -(taurise(ll)/taudecay(ll))**(taudecay(ll)/(taudecay(ll)-taurise(ll))))
                        
                !S3(ll) = (EXP(-dt*(ii-st3(ll))/taurise(ll))-EXP(-dt*(ii-st3(ll))/taudecay(ll))) &
                !         /((taurise(ll)/taudecay(ll))**(taurise(ll)/(taudecay(ll)-taurise(ll))) &
                !         -(taurise(ll)/taudecay(ll))**(taudecay(ll)/(taudecay(ll)-taurise(ll))))
                        
                S(ll) =  S1(ll) + S2(ll)! + S3(ll)
                SS(ll, ii + 1) = S(ll)
            END DO 
            DO jj = 1, numnet * neurons ! Post
                DO kk = 1, numnet * neurons ! Pre
                    IF (A(jj, kk) .NE. 0) THEN
                        IF ((ii .GT. delay(jj,kk))) THEN
                            xIsyn(jj,kk) = SS(kk, ii - delay(jj, kk)) * synapse(jj,kk) * (v(jj) - Esyn(kk))
                        ELSE
                            xIsyn(jj,kk) = 0
                        END IF
                    END IF
                END DO
                Isyn(jj) = 0
                DO kk=1, numnet * neurons
                    Isyn(jj) = Isyn(jj) + xIsyn(jj,kk)
                END DO
            END DO           
        END DO

        !!-----------LOOP OVER TIME AND NEURONS----END----------------------

        PRINT*, "ready to write"

        WRITE(filename1,'(A11,I3,A4,I3,A3,I3,A4)') "rho_deeout_",INT(xdeeout),"_ens",INT(ens),"_dI",INT(dI),".txt" 
        OPEN(UNIT=ens, FILE=TRIM(filename1), ACTION = "WRITE", STATUS = "UNKNOWN")

        !!---------- WRITE TO FILE, SECTION----------------------------------


        DO ii = 1,time
                WRITE(ens,*) (rho(jj,ii) , jj = 1, numnet * neurons)
        END DO
        CLOSE(ens)
        PRINT*, "rho has written"
        
        CALL cpu_time(finish)
        PRINT*, 'Time = ', finish-start, ' seconds.'

!        DO ii = 1,time
!            WRITE(700, *) (IextSave(jj, ii), jj = 1, numnet * neurons )
!        END DO
!        CLOSE(700)
!        PRINT*,"Iext has written to disk"
!
!        DO ii = 1, time
!            WRITE(800, *) (SS(jj, ii), jj = 1, numnet * neurons)
!        END DO
!        CLOSE(800)
!        PRINT*,"SS has written to disk"DO ii = 1, time

!
    END DO

END PROGRAM hh_double
