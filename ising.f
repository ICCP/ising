      PROGRAM ISING
      ! 2D Monte Carlo simulation of Ising Model with Metropolis algorithm
      ! Program written by Tridip Das (dastridi)
      ! This program calculates magnetization with temperature variation 
      ! Susceptibility, Energy and Specific heat with temperature
      ! Also finds the critical temperature at which magnetization is
      ! lost. Scaled temperature is used, Boltzmann const, kB assumed as 1
      !
        IMPLICIT NONE
      ! Variable Declaration
        INTEGER :: I, J, M, N ! counters
        INTEGER, ALLOCATABLE :: A(:,:),AOLD(:,:),ANEW(:,:) ! 2D array containing spins
        INTEGER :: NROWS, NCOLS ! Defining # of rows and columns
      !
        REAL :: TEMPERATURE, BETA, MAGNETIZATION, ENERGY
        REAL :: MAG2_AVG, ENERGY2_AVG, SUSP, Cv
        REAL :: DELTA_E, MAG_AVG, RAND_NUM, ENERGY_AVG
        INTEGER :: TEMP_ITER, ITER_AVG, MC_ITER, TEMP_VAR
        INTEGER ::  NO_OF_MC_STEPS, NO_OF_ITER, TEMP_STEPS
        REAL :: T_INITIAL, T_FINAL, TEMP_INTERVAL, START, FINISH
      !
      !
        PRINT *, " MONTE CARLO SIMULATION OF 2D ISING MODEL "
        PRINT *, "   SIMULATION WITH METROPOLIS ALGORITHM  "
        PRINT *, "   APPLYING PERIODIC BOUNDARY CONDITION  "
        PRINT *, "--------Code developed by Tridip Das------"
        CALL CPU_TIME(START)
        PRINT *, " Start time ", START
        PRINT *, "------------------------------------------"
      ! INPUT PARAMETERS FROM inp_file
      !  OPEN(UNIT=11,FILE='inpfile',STATUS='OLD',ACTION='READ')
      ! Uncomment below lines to get input from user
      !  PRINT *, " PROVIDE THE NUMBER OF ROWS FOR ISING MODEL"
      !  READ(*,*) NROWS
      !  PRINT *, " PROVIDE THE NUMBER OF COLS FOR ISING MODEL"
      !  READ(*,*) NCOLS
      !  PRINT *, " NUMBER OF SIMULATION STEPS TO RUN FOR MC             &
!     &             EQULIBRIATION "
      !  READ(*,*) NO_OF_MC_STEPS
      !  PRINT *, " NUMBER OF ITERATIONS TO PERFORM AT EACH TEMPERATURE  &
!     &             STEP TO CALCULATE AVERAGE OVER MC PREDICTION"
      !  READ(*,*) NO_OF_ITER
      !  PRINT *, " PROVIDE INITIAL, FINAL TEMPERATURE AND INTERVAL"
      !  PRINT *, " SEPARATED BY ENTER "
      !  READ(*,*) T_INITIAL
      !  READ(*,*) T_FINAL
      !  READ(*,*) TEMP_INTERVAL
      !  READ(11,*);READ(11,*) NO_OF_STEPS
      !  READ(11,*);READ(11,*) NO_OF_ITER
      !  READ(11,*);READ(11,*) T_INITIAL
      !  READ(11,*);READ(11,*) T_FINAL
      !  READ(11,*);READ(11,*) TEMP_INTERVAL
      !  CLOSE(11)
      !
        NROWS = 10 ; NCOLS = 10; 
        NO_OF_MC_STEPS = 1000000; NO_OF_ITER = 1
        T_INITIAL = 1.0; T_FINAL = 4.0; TEMP_INTERVAL = 0.1 
        TEMP_STEPS = INT((T_FINAL - T_INITIAL)/TEMP_INTERVAL)
      !
        ALLOCATE(A(NROWS+2,NCOLS+2))
        ALLOCATE(AOLD(NROWS+2,NCOLS+2))
        ALLOCATE(ANEW(NROWS+2,NCOLS+2))
      !
      ! Initialize Spin array
        A(:,:) = 1
      ! Write in a OUTPUT file
      !
        OPEN(UNIT=31,FILE='Spin_states.out',STATUS='REPLACE',           &
     &       ACTION='WRITE')
      !
        OPEN(UNIT=33,FILE='Magnetization.out',STATUS='REPLACE',         &
     &       ACTION='WRITE')
      !  WRITE(32,*) T_INITIAL
      !  WRITE(32,*) TEMP_INTERVAL
      !  WRITE(32,*) T_FINAL
      !
        WRITE(33,*) "Temperature ", " Iteration ", " Magnetization"
      !
        OPEN(UNIT=32,FILE='temperature_sweep.dat',STATUS='REPLACE',     &
     &       ACTION='WRITE')
        WRITE(32,*) " Temperature ", " Magnetization ",                 &
     &              " Susceptibility ", " Energy ", " Cv "
      !
        MAGNETIZATION = 0.0; ENERGY = 0.0;
      !
        TEMP_SCAN: DO TEMP_ITER = 0, TEMP_STEPS
      !   
          TEMPERATURE = T_INITIAL + TEMP_ITER * TEMP_INTERVAL
          BETA = 1.0/TEMPERATURE
          MAG_AVG = 0.0; MAG2_AVG = 0.0
          ENERGY_AVG = 0.0; ENERGY2_AVG = 0.0  
      !
          TAKE_AVG: DO ITER_AVG = 1, NO_OF_ITER
               AOLD = A
      ! 
            MONTE_CARLO: DO MC_ITER = 1, NO_OF_MC_STEPS
                ANEW = AOLD
      !         Ramdomly change the spin
                CALL RANDOM_NUMBER(RAND_NUM)
                M = NINT((NROWS-1)*RAND_NUM+2)      ! Choose a random row
      !
                CALL RANDOM_NUMBER(RAND_NUM)      
                N = NINT((NCOLS-1)*RAND_NUM+2)      ! Choose a random col
                TEMP_VAR = -1*ANEW(M,N)
                ANEW(M,N) = TEMP_VAR
      !  Calculating change in energy due to new configuration
      !
                DELTA_E = -2*ANEW(M,N)*(ANEW(M-1,N)+ANEW(M+1,N)         &
     &                               + ANEW(M,N-1)+ANEW(M,N+1))
      ! 
                CALL RANDOM_NUMBER(RAND_NUM)
      !
                 IF  (EXP(-BETA*DELTA_E)> RAND_NUM) THEN
                    AOLD = ANEW
                      IF (M == 2) AOLD(NROWS+2,N) = TEMP_VAR
                      IF (M == NROWS+1) AOLD(1,N) = TEMP_VAR
                      IF (N == 2) AOLD(M,NCOLS+2) = TEMP_VAR
                      IF (N == NCOLS+1) AOLD(M,1) = TEMP_VAR
                ELSE
                    CONTINUE    
                ENDIF
      !
                MAGNETIZATION = ABS(SUM(AOLD(2:NROWS+1,2:NCOLS+1))/     &
     &                                 (NROWS*NCOLS*1.0))
      !
                IF ((MOD(NINT(TEMPERATURE*10),2) == 0) .AND.            &
     &               MOD(MC_ITER,1000) == 0) THEN
                        WRITE(33,*) TEMPERATURE, MC_ITER, MAGNETIZATION
                ENDIF
      !
            ENDDO MONTE_CARLO
      !         
           MAG_AVG = MAG_AVG + MAGNETIZATION
           MAG2_AVG = MAG2_AVG + MAGNETIZATION**2
      !
           DO I = 2, NROWS + 1
             DO J = 2, NCOLS + 1
               ENERGY = ENERGY -ANEW(I,J)*(ANEW(I-1,J)+ANEW(I+1,J)      &
     &                                    +ANEW(I,J-1)+ANEW(I,J+1))
             ENDDO
           ENDDO
           ENERGY = ENERGY/(NROWS*NCOLS*2)
           ENERGY_AVG = ENERGY_AVG + ENERGY
           ENERGY2_AVG = ENERGY2_AVG + ENERGY**2
      !
          ENDDO TAKE_AVG
      !
           MAG_AVG = MAG_AVG/NO_OF_ITER
           MAG2_AVG = MAG2_AVG/NO_OF_ITER
           ENERGY_AVG = ENERGY_AVG/NO_OF_ITER
           ENERGY2_AVG = ENERGY2_AVG/NO_OF_ITER
      !
          SUSP = BETA*(MAG2_AVG - MAG_AVG**2)
          Cv   = (BETA/TEMPERATURE)*(ENERGY2_AVG - ENERGY_AVG**2)
      !
          WRITE(32,*) TEMPERATURE, MAG_AVG, SUSP, ENERGY, Cv
      !      
      !    WRITE(31,*) TEMPERATURE
      !    WRITE(31,"(10I2)") AOLD(2:NROWS+1,2:NCOLS+1) + 1
      !    WRITE(31,*) "..................................."
!       
        ENDDO TEMP_SCAN
            WRITE(31,*) AOLD(2:NROWS+1,2:NCOLS+1)
        CLOSE(31)
        CLOSE(32)
        CLOSE(33)
        print "(12I2)",  AOLD
        PRINT *, "........................................"
        PRINT *, ENERGY/(NROWS*NCOLS*2)
        CALL CPU_TIME(FINISH)
        PRINT *, " End time ", FINISH
      END PROGRAM ISING
        
