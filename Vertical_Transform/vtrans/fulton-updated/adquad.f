c=======================================================================
c adquad.f
c
c Double precision Fortran 77 version of "adquad.f" as adapted from
c various versions of "adquad.f" that were obtained from various
c versions of the NCAR Spectral Transform Shallow Water Model (STSWM).
c=======================================================================
c
C PACKAGE ADQUAD         DESCRIPTION OF USER ENTRIES ADQUAD AND QUAD
C                        FOLLOWS THIS PACKAGE DESCRIPTION.
C
C LATEST REVISION        DECEMBER 2002
C
C PURPOSE                TO CALCULATE AUTOMATICALLY THE INTEGRAL OF F(X)
C                        OVER THE FINITE INTERVAL (A,B) WITH RELATIVE
C                        ERROR LESS THAN EPSIL.  THE PACKAGE CONTAINS
C                        TWO ROUTINES, ADQUAD AND QUAD.  ADQUAD IS
C                        RECOMMENDED AS THE PRIMARY USER ENTRY, THOUGH
C                        USERS DESIRING MORE CONTROL MAY WISH TO CALL
C                        SUBROUTINE QUAD THEMSELVES.
C
C SPECIAL CONDITIONS     IF THE RANGE INCLUDES DISCONTINUITIES OR
C                        SINGULARITIES, SPECIAL TREATMENT MAY BE
C                        NECESSARY.  DISCONTINUITIES MAY BE TREATED BY
C                        SUBDIVIDING THE INTERVAL AT EACH DISCONTINUITY.
C                        THE ROUTINE WILL COPE WITH A SINGULARITY OF THE
C                        TYPE X**P, WHERE 0 .GT. P .GT. -.9 AT EITHER OR
C                        BOTH ENDS OF THE RANGE, BUT THE CALCULATION IS
C                        TIME CONSUMING, INVOLVING MANY FUNCTION
C                        EVALUATIONS.
C
C I/O                    ERROR MESSAGES ARE WRITTEN DIRECTLY.
C
C PRECISION              SINGLE
C
C LANGUAGE               FORTRAN
C
C HISTORY                THIS PACKAGE IS A FORTRAN 66 STANDARDIZED
C                        VERSION OF ACM ALGORITHM 468.  IT WAS
C                        CODED AT NCAR IN THE EARLY 1970'S
C                        BY MEMBERS OF THE SCIENTIFIC COMPUTING
C                        DIVISION IN BOULDER, COLORADO.
C                        MODIFIED TO SAME INTERFACE AS NAG ROUTINE
C
C REFERENCES             SEE PUBLICATION
C                          ALGORITHM FOR AUTOMATIC NUMERICAL
C                          INTEGRATION OVER A FINITE INTERVAL
C                        BY T.N.L. PATTERSON, 1973
C                        COMMUNICATIONS OF THE ACM, VOL. 16,
C                        PP. 694-699.
C***********************************************************************
C
C FUNCTION ADQUAD(A,B,F,EPSIL,NPTS,RELERR,IER)
C
C PURPOSE                TO CALCULATE AUTOMATICALLY THE INTEGRAL OF F(X)
C                        OVER THE FINITE INTERVAL (A,B) WITH RELATIVE
C                        ERROR LESS THAN EPSIL.
C
C USAGE                  R = ADQUAD(A,B,F,EPSIL,NPTS,RELERR,IER)
C                        WHERE R IS THE VALUE OF THE INTEGRAL.
C
C ARGUMENTS
C
C ON INPUT               A
C                          LOWER LIMIT OF INTEGRATION.
C
C                        B
C                          UPPER LIMIT OF INTEGRATION.
C
C                        F
C                          F(X) IS USER WRITTEN FUNCTION TO CALCULATE
C                          THE INTEGRAND.  F MUST BE DECLARED EXTERNAL
C                          IN THE CALLING PROGRAM.
C
C                        EPSIL
C                          REQUIRED RELATIVE ERROR OF THE INTEGRAL.  IT
C                          IS RECOMMENDED THAT EPSIL BE LESS THAN .001.
C
C ON OUTPUT              NPTS
C                          NUMBER OF INTEGRAND EVALUATIONS PERFORMED.
C
C                        RELERR
C                          CRUDE ESTIMATE OF RELATIVE ERROR OBTAINED BY
C                          SUMMING THE ABSOLUTE VALUES OF THE ERRORS
C                          PRODUCED BY QUAD ON EACH SUBINTERVAL AND
C                          DIVIDING BY THE CALCULATED VALUE OF THE
C                          INTEGRAL.
C
C                        IER
C                          ERROR PARAMETER.
C                          = 0  REQUIRED ACCURACY ACHIEVED.
C                          = 1  RELAXED CONVERGENCE FOR AT LEAST ONE
C                               SUBINTERVAL.  ACCURACY SHOULD BE CHECKED
C                               BY EXAMINING RELERR.  IF A SUBINTERVAL
C                               DOES NOT CONVERGE WITH RELATIVE ERROR
C                               EPSIL, A RELAXED CONVERGENCE CRITERION
C                               IS APPLIED.  THIS IS
C                                 ABS(RESULT(K)-RESULT(K-1)) .LE.
C                                 ESTIM*EPSIL
C                               WHERE ESTIM IS THE ESTIMATE OF THE
C                               INTEGRAL OBTAINED BY APPLYING QUAD TO
C                               THE WHOLE INTERVAL.  THIS ALLOWS FOR THE
C                               SITUATION WHERE NEARLY ALL THE ERROR MAY
C                               BE CONCENTRATED IN ONE SUBINTERVAL AS
C                               FOR A SINGULARITY.
C                          = 2  INTERVAL STACK OVERFLOW.  CHECK
C                               ACCURACY.  THE INTERVAL STACK ALLOWS FOR
C                               50 INTERVALS NEEDING FURTHER
C                               SUBDIVISION.  IF THE STACK IS FULL WHEN
C                               ANOTHER INTERVAL NEEDS TO BE ADDED, THAT
C                               INTERVAL IS ACCEPTED AS IT STANDS, EVEN
C                               THOUGH CONVERGENCE HAS FAILED.  THIS
C                               ENTAILS LOSS OF ACCURACY.
C
C METHOD                 ADQUAD CALLS SUBROUTINE QUAD TO INTEGRATE
C                        OVER A FINITE INTERVAL EMPLOYING SUCCESSIVE
C                        ADAPTIVE SUBDIVISION IF CONVERGENCE TO THE
C                        REQUIRED ACCURACY IS NOT ACHIEVED.
C***********************************************************************
C
C     SUBROUTINE QUAD(A,B,F,RESULT,EPSIL,NPTS,K,IER)
C
C DIMENSION OF           RESULT(8)
C ARGUMENTS
C
C PURPOSE                TO CALCULATE AUTOMATICALLY THE INTEGRAL OF F(X)
C                        OVER THE FINITE INTERVAL (A,B) WITH RELATIVE
C                        ERROR LESS THAN EPSIL.
C
C USAGE                  CALL QUAD(A,B,F,RESULT,EPSIL,NPTS,K,IER)
C
C ARGUMENTS
C
C ON INPUT               A
C                          LOWER LIMIT OF INTEGRATION.
C
C                        B
C                          UPPER LIMIT OF INTEGRATION.
C
C                        F
C                          F(X) IS USER WRITTEN FUNCTION TO CALCULATE
C                          THE INTEGRAND.  F MUST BE DECLARED EXTERNAL
C                          IN THE CALLING PROGRAM.
C
C                        EPSIL
C                          REQUIRED RELATIVE ERROR OF THE INTEGRAL.
C                          THIS SHOULD BE LESS THAN .001.
C
C ON OUTPUT              RESULT
C                          ARRAY IN WHICH ARE RETURNED THE RESULTS OF
C                          APPLYING 1,3,7,15,31,63,127, AND
C                          255-POINT FORMULAE SUCCESSIVELY.  THE NUMBER
C                          OF FORMULAE ACTUALLY USED WILL DEPEND UPON
C                          EPSIL.  THE DIMENSION OF RESULT IN THE
C                          CALLING PROGRAM SHOULD BE AT LEAST 8.
C
C                        NPTS
C                          NUMBER OF INTEGRAND EVALUATIONS PERFORMED.
C
C                        K
C                          RESULT(K) CONTAINS THE VALUE OF THE INTEGRAL
C                          TO THE REQUIRED ACCURACY.  K IS DETERMINED
C                          FROM THE CONVERGENCE CRITERION:
C                            ABS(RESULT(K)-RESULT(K-1)) .LE.
C                            EPSIL*ABS(RESULT(K))
C
C                        IER
C                          ERROR PARAMETER
C                          = 0  REQUIRED ACCURACY ACHIEVED.
C                          = 1  REQUIRED ACCURACY NOT ACHIEVED AFTER
C                               WORKING THROUGH ALL EIGHT FORMULAE.
C
C TIMING                 IN ALL BUT THE MOST TRIVIAL EXAMPLES, THE
C                        TIMING IS DETERMINED BY THE TIME, T, NEEDED FOR
C                        AN INTEGRAND EVALUATION BY THE USER PROVIDED
C                        FUNCTION F(X).  THE TOTAL NUMBER, NPTS, OF
C                        INTEGRAND EVALUATIONS DEPENDS CRITICALLY UPON
C                        THE BEHAVIOR OF F(X) WITHIN THE INTERVAL (A,B)
C                        AND MAY RUN INTO THOUSANDS FOR PATHOLOGICAL
C                        CASES.  THE TOTAL TIME IS APPROXIMATELY NPTS*T.
C
C METHOD                 QUAD IMPLEMENTS THE BASIC ALGORITHM WHICH
C                        INTEGRATES OVER THE WHOLE INTERVAL USING A
C                        SEQUENCE OF INTERLEAVING 1, 3, 7, 15, 31, 63,
C                        127 AND 255-POINT EXTENDED GAUSS-TYPE
C                        QUADRATURE FORMULAE.  SINCE EACH SUCCESSIVE
C                        FORMULA EMPLOYS ALL POINTS USED BY ITS
C                        PREDECESSOR, NO INTEGRAND VALUES ARE WASTED
C                        WHEN THE ORDER OF THE INTEGRATION FORMULA IS
C                        INCREASED.  THE PROCESS IS DEEMED TO CONVERGE
C                        WHEN THE RELATIVE DIFFERENCE BETWEEN VALUES
C                        OF THE INTEGRAL OBTAINED FROM TWO SUCCESSIVE
C                        FORMULAE IS LESS THAN EPSIL.
C
C***********************************************************************
C
      SUBROUTINE QUAD(A,B,F,RESULT,EPSIL,NPTS,K,IER)
C
      IMPLICIT NONE
      DOUBLE PRECISION  A, B, EPSIL
      DOUBLE PRECISION  F
      DOUBLE PRECISION  RESULT(8)
      INTEGER           NPTS, K, IER
C
      EXTERNAL AQDAT1
C
      DOUBLE PRECISION  FUNCT(127), P(381)
      DOUBLE PRECISION  SUM, DIFF, FZERO, ACUM, X
      INTEGER           I, IOLD, INEW, J
C
      COMMON /ADQD1/ P
C
      IER = 0.0D+00
C
C CHECK FOR TRIVIAL CASE
C
      IF (A .EQ. B) GO TO 107
C
C SCALE FACTORS
C
      SUM  = (B+A)*5.0D-01
      DIFF = (B-A)*5.0D-01
C
C 1-POINT FORMULA
C
      FZERO = F(SUM)
      RESULT(1) = 2.0D+00*FZERO*DIFF
      I = 0
      IOLD = 0
      INEW = 1
      K = 2
      ACUM = 0.0D+00
      GO TO 103
  101 CONTINUE
      IF (K .EQ. 8) GO TO 105
      K = K+1
      ACUM = 0.0D+00
C
C CONTRIBUTION FROM FUNCTION VALUES ALREADY COMPUTED
C
      DO 102 J=1,IOLD
         I = I+1
         ACUM = ACUM+P(I)*FUNCT(J)
  102 CONTINUE
C
C CONTRIBUTION FROM NEW FUNCTION VALUES
C
  103 CONTINUE
      IOLD = IOLD+INEW
      DO 104 J=INEW,IOLD
         I = I+1
         X = P(I)*DIFF
         FUNCT(J) = F(SUM+X)+F(SUM-X)
         I = I+1
         ACUM = ACUM+P(I)*FUNCT(J)
  104 CONTINUE
      INEW = IOLD+1
      I = I+1
      RESULT(K) = (ACUM+P(I)*FZERO)*DIFF
C
C CHECK FOR CONVERGENCE
C
      IF (ABS(RESULT(K)-RESULT(K-1))-ABS(EPSIL*RESULT(K))) 106,106,101
C
C CONVERGENCE NOT ACHIEVED
C
  105 CONTINUE
      IER = 1
      IF (EPSIL .LT. 0.0D+00) GO TO 106
         WRITE(*,108)
  108    FORMAT(/,' QUAD: REQUIRED ACCURACY NOT ACHIEVED')
C
C NORMAL TERMINATION
C
  106 CONTINUE
      NPTS = INEW+IOLD
      RETURN
C
C TRIVIAL CASE
C
  107 CONTINUE
      K = 2
      RESULT(1) = 0.0D+00
      RESULT(2) = 0.0D+00
      NPTS = 0
      RETURN
      END
C
      FUNCTION ADQUAD(A,B,F,EPSIL,NPTS,RELERR,IER)
C
      IMPLICIT NONE
      DOUBLE PRECISION  ADQUAD
      DOUBLE PRECISION  A, B, EPSIL, RELERR
      DOUBLE PRECISION  F
      EXTERNAL          F
      INTEGER           NPTS, IER
C
      DOUBLE PRECISION  RESULT(8), STACK(10000)
      DOUBLE PRECISION  COMP, ESTIM, SUB1, SUB2, SUB3
      INTEGER           ISMAX
      INTEGER           K, IC, IS, ISTACK, NF, IIER
C
      DATA ISMAX/10000/
C
      CALL QUAD(A,B,F,RESULT,-EPSIL,NPTS,K,IER)
      ADQUAD = RESULT(K)
      RELERR = 0.0D+00
      IF (ADQUAD .NE. 0.0D+00)
     &   RELERR = ABS((RESULT(K)-RESULT(K-1))/ADQUAD)
C
C CHECK IF SUBDIVISION IS NEEDED
C
      IF (IER .EQ. 0) RETURN
C
C SUBDIVIDE
C
      ESTIM = ABS(ADQUAD*EPSIL)
      RELERR = 0.0D+00
      ADQUAD = 0.0D+00
      IS = 1
      IC = 0
      ISTACK = 0
      SUB1 = A
      SUB3 = B
  101 CONTINUE
      SUB2 = (SUB1+SUB3)*5.0D-01
      CALL QUAD(SUB1,SUB2,F,RESULT,-EPSIL,NF,K,IER)
      NPTS = NPTS+NF
      COMP = ABS(RESULT(K)-RESULT(K-1))
C
C CHECK FOR INTERVAL CONVERGENCE
C
      IF (IER .EQ. 0) GO TO 103
C
C CHECK FOR RELAXED CONVERGENCE
C
      IF (COMP .LE. ESTIM) GO TO 110
C
C CHECK FOR FULL STACK
C
      IF (IS .GE. ISMAX) GO TO 102
C
C SUBINTERVAL(SUB1,SUB2) DID NOT CONVERGE. STACK FOR FUTURE PROCESSING.
C
      STACK(IS) = SUB1
      IS = IS+1
      STACK(IS) = SUB2
      IS = IS+1
      GO TO 104
C
C STACK FULL, CONVERGENCE FAILURE IGNORED. ISTACK SET TO 1 TO INDICATE
C THIS.
C
  102 CONTINUE
      ISTACK = 1
C
C FIRST HALF INTERVAL CONVERGED
C
  103 CONTINUE
      ADQUAD = ADQUAD+RESULT(K)
      RELERR = RELERR+COMP
C
C PROCESS SECOND HALF INTERVAL
C
  104 CONTINUE
      CALL QUAD(SUB2,SUB3,F,RESULT,-EPSIL,NF,K,IER)
      NPTS = NPTS+NF
      COMP = ABS(RESULT(K)-RESULT(K-1))
C
C CHECK FOR CONVERGENCE
C
      IF (IER .EQ. 0) GO TO 105
C
C CHECK FOR RELAXED CONVERGENCE.
C
      IF (COMP .LE. ESTIM) GO TO 111
C
C NO CONVERGENCE. SUBDIVIDE INTERVAL(SUB2,SUB3) AND REPEAT.
C
      SUB1 = SUB2
      GO TO 101
C
C SECOND HALF INTERVAL CONVERGED
C
  105 CONTINUE
      ADQUAD = ADQUAD+RESULT(K)
      RELERR = RELERR+COMP
C
C CHECK FOR ALL DELINQUENT INTERVALS NOW PROCESSED.
C
      IF (IS .EQ. 1) GO TO 106
C
C SUBDIVIDE DELINQUENT INTERVAL LAST STACKET).
C
      IS = IS-1
      SUB3 = STACK(IS)
      IS = IS-1
      SUB1 = STACK(IS)
      GO TO 101
C
C FINAL RESULT WITH SUBDIVISION
C
  106 CONTINUE
      IER = IC
      IF(ISTACK.NE.0)IER=2
      IIER = IER+1
      GO TO (109,107,108),IIER
  107 CONTINUE
         WRITE(*,112) 
  112    FORMAT(/,' ADQUAD: RELAXED CONVERGENCE')
      GO TO 109
  108 CONTINUE
         WRITE(*,113)
  113    FORMAT(/,' ADQUAD: INTERVAL STACK OVERFLOW, CHECK ACCURACY')
  109 CONTINUE
      IF (ADQUAD .NE. 0.0D+00) RELERR = RELERR/ABS(ADQUAD)
      RETURN
C
C RELAXED CONVERGENCE
C
  110 CONTINUE
      IC = 1
      GO TO 103
  111 CONTINUE
      IC = 1
      GO TO 105
      END
C
      BLOCKDATA AQDAT1
C
      IMPLICIT NONE
      DOUBLE PRECISION  P(381)
      COMMON /ADQD1/ P
C
      DATA P(1),P(2),P(3)/
     & 7.74596669241483D-01, 5.55555555555557D-01, 8.88888888888889D-01/
      DATA P(4),P(5),P(6)/
     & 2.68488089868333D-01, 9.60491268708019D-01, 1.04656226026467D-01/
      DATA P(7),P(8),P(9)/
     & 4.34243749346802D-01, 4.01397414775962D-01, 4.50916538658474D-01/
      DATA P(10),P(11),P(12)/
     & 1.34415255243784D-01, 5.16032829970798D-02, 2.00628529376989D-01/
      DATA P(13),P(14),P(15)/
     & 9.93831963212756D-01, 1.70017196299402D-02, 8.88459232872258D-01/
      DATA P(16),P(17),P(18)/
     & 9.29271953151245D-02, 6.21102946737228D-01, 1.71511909136392D-01/
      DATA P(19),P(20),P(21)/
     & 2.23386686428967D-01, 2.19156858401588D-01, 2.25510499798206D-01/
      DATA P(22),P(23),P(24)/
     & 6.72077542959908D-02, 2.58075980961766D-02, 1.00314278611795D-01/
      DATA P(25),P(26),P(27)/
     & 8.43456573932111D-03, 4.64628932617579D-02, 8.57559200499902D-02/
      DATA P(28),P(29),P(30)/
     & 1.09578421055925D-01, 9.99098124967666D-01, 2.54478079156187D-03/
      DATA P(31),P(32),P(33)/
     & 9.81531149553739D-01, 1.64460498543878D-02, 9.29654857429739D-01/
      DATA P(34),P(35),P(36)/
     & 3.59571033071293D-02, 8.36725938168868D-01, 5.69795094941234D-02/
      DATA P(37),P(38),P(39)/
     & 7.02496206491528D-01, 7.68796204990037D-02, 5.31319743644374D-01/
      DATA P(40),P(41),P(42)/
     & 9.36271099812647D-02, 3.31135393257977D-01, 1.05669893580235D-01/
      DATA P(43),P(44),P(45)/
     & 1.12488943133187D-01, 1.11956873020953D-01, 1.12755256720769D-01/
      DATA P(46),P(47),P(48)/
     & 3.36038771482077D-02, 1.29038001003512D-02, 5.01571393058995D-02/
      DATA P(49),P(50),P(51)/
     & 4.21763044155885D-03, 2.32314466399103D-02, 4.28779600250078D-02/
      DATA P(52),P(53),P(54)/
     & 5.47892105279628D-02, 1.26515655623007D-03, 8.22300795723591D-03/
      DATA P(55),P(56),P(57)/
     & 1.79785515681282D-02, 2.84897547458336D-02, 3.84398102494556D-02/
      DATA P(58),P(59),P(60)/
     & 4.68135549906281D-02, 5.28349467901166D-02, 5.59784365104763D-02/
      DATA P(61),P(62),P(63)/
     & 9.99872888120358D-01, 3.63221481845531D-04, 9.97206259372224D-01/
      DATA P(64),P(65),P(66)/
     & 2.57904979468569D-03, 9.88684757547428D-01, 6.11550682211726D-03/
      DATA P(67),P(68),P(69)/
     & 9.72182874748583D-01, 1.04982469096213D-02, 9.46342858373402D-01/
      DATA P(70),P(71),P(72)/
     & 1.54067504665595D-02, 9.10371156957005D-01, 2.05942339159128D-02/
      DATA P(73),P(74),P(75)/
     & 8.63907938193691D-01, 2.58696793272147D-02, 8.06940531950218D-01/
      DATA P(76),P(77),P(78)/
     & 3.10735511116880D-02, 7.39756044352696D-01, 3.60644327807826D-02/
      DATA P(79),P(80),P(81)/
     & 6.62909660024781D-01, 4.07155101169443D-02, 5.77195710052045D-01/
      DATA P(82),P(83),P(84)/
     & 4.49145316536321D-02, 4.83618026945841D-01, 4.85643304066732D-02/
      DATA P(85),P(86),P(87)/
     & 3.83359324198731D-01, 5.15832539520484D-02, 2.77749822021825D-01/
      DATA P(88),P(89),P(90)/
     & 5.39054993352661D-02, 1.68235251552208D-01, 5.54814043565595D-02/
      DATA P(91),P(92),P(93)/
     & 5.63443130465928D-02, 5.62776998312542D-02, 5.63776283603847D-02/
      DATA P(94),P(95),P(96)/
     & 1.68019385741038D-02, 6.45190005017574D-03, 2.50785696529497D-02/
      DATA P(97),P(98),P(99)/
     & 2.10881524572663D-03, 1.16157233199551D-02, 2.14389800125039D-02/
      DATA P(100),P(101),P(102)/
     & 2.73946052639814D-02, 6.32607319362634D-04, 4.11150397865470D-03/
      DATA P(103),P(104),P(105)/
     & 8.98927578406411D-03, 1.42448773729168D-02, 1.92199051247278D-02/
      DATA P(106),P(107),P(108)/
     & 2.34067774953141D-02, 2.64174733950583D-02, 2.79892182552381D-02/
      DATA P(109),P(110),P(111)/
     & 1.80739564445388D-04, 1.28952408261042D-03, 3.05775341017553D-03/
      DATA P(112),P(113),P(114)/
     & 5.24912345480885D-03, 7.70337523327974D-03, 1.02971169579564D-02/
      DATA P(115),P(116),P(117)/
     & 1.29348396636074D-02, 1.55367755558440D-02, 1.80322163903913D-02/
      DATA P(118),P(119),P(120)/
     & 2.03577550584721D-02, 2.24572658268161D-02, 2.42821652033366D-02/
      DATA P(121),P(122),P(123)/
     & 2.57916269760242D-02, 2.69527496676331D-02, 2.77407021782797D-02/
      DATA P(124),P(125),P(126)/
     & 2.81388499156271D-02, 9.99982430354891D-01, 5.05360952078625D-05/
      DATA P(127),P(128),P(129)/
     & 9.99598799671912D-01, 3.77746646326985D-04, 9.98316635318407D-01/
      DATA P(130),P(131),P(132)/
     & 9.38369848542380D-04, 9.95724104698407D-01, 1.68114286542147D-03/
      DATA P(133),P(134),P(135)/
     & 9.91495721178104D-01, 2.56876494379402D-03, 9.85371499598521D-01/
      DATA P(136),P(137),P(138)/
     & 3.57289278351730D-03, 9.77141514639705D-01, 4.67105037211432D-03/
      DATA P(139),P(140),P(141)/
     & 9.66637851558417D-01, 5.84344987583563D-03, 9.53730006425761D-01/
      DATA P(142),P(143),P(144)/
     & 7.07248999543356D-03, 9.38320397779592D-01, 8.34283875396818D-03/
      DATA P(145),P(146),P(147)/
     & 9.20340025470011D-01, 9.64117772970252D-03, 8.99744899776941D-01/
      DATA P(148),P(149),P(150)/
     & 1.09557333878379D-02, 8.76513414484705D-01, 1.22758305600827D-02/
      DATA P(151),P(152),P(153)/
     & 8.50644494768350D-01, 1.35915710097655D-02, 8.22156254364980D-01/
      DATA P(154),P(155),P(156)/
     & 1.48936416648152D-02, 7.91084933799848D-01, 1.61732187295777D-02/
      DATA P(157),P(158),P(159)/
     & 7.57483966380512D-01, 1.74219301594641D-02, 7.21423085370098D-01/
      DATA P(160),P(161),P(162)/
     & 1.86318482561388D-02, 6.82987431091078D-01, 1.97954950480975D-02/
      DATA P(163),P(164),P(165)/
     & 6.42276642509760D-01, 2.09058514458120D-02, 5.99403930242243D-01/
      DATA P(166),P(167),P(168)/
     & 2.19563663053178D-02, 5.54495132631931D-01, 2.29409642293877D-02/
      DATA P(169),P(170),P(171)/
     & 5.07687757533716D-01, 2.38540521060385D-02, 4.59130011989833D-01/
      DATA P(172),P(173),P(174)/
     & 2.46905247444876D-02, 4.08979821229888D-01, 2.54457699654648D-02/
      DATA P(175),P(176),P(177)/
     & 3.57403837831532D-01, 2.61156733767061D-02, 3.04576441556714D-01/
      DATA P(178),P(179),P(180)/
     & 2.66966229274503D-02, 2.50678730303482D-01, 2.71855132296248D-02/
      DATA P(181),P(182),P(183)/
     & 1.95897502711100D-01, 2.75797495664819D-02, 1.40424233152560D-01/
      DATA P(184),P(185),P(186)/
     & 2.78772514766137D-02, 8.44540400837110D-02, 2.80764557938172D-02/
      DATA P(187),P(188),P(189)/
     & 2.81846489497457D-02, 2.81763190330167D-02, 2.81888141801924D-02/
      DATA P(190),P(191),P(192)/
     & 8.40096928705192D-03, 3.22595002508787D-03, 1.25392848264749D-02/
      DATA P(193),P(194),P(195)/
     & 1.05440762286332D-03, 5.80786165997757D-03, 1.07194900062519D-02/
      DATA P(196),P(197),P(198)/
     & 1.36973026319907D-02, 3.16303660822264D-04, 2.05575198932735D-03/
      DATA P(199),P(200),P(201)/
     & 4.49463789203206D-03, 7.12243868645840D-03, 9.60995256236391D-03/
      DATA P(202),P(203),P(204)/
     & 1.17033887476570D-02, 1.32087366975291D-02, 1.39946091276191D-02/
      DATA P(205),P(206),P(207)/
     & 9.03727346587510D-05, 6.44762041305726D-04, 1.52887670508776D-03/
      DATA P(208),P(209),P(210)/
     & 2.62456172740443D-03, 3.85168761663987D-03, 5.14855847897819D-03/
      DATA P(211),P(212),P(213)/
     & 6.46741983180368D-03, 7.76838777792199D-03, 9.01610819519566D-03/
      DATA P(214),P(215),P(216)/
     & 1.01788775292361D-02, 1.12286329134080D-02, 1.21410826016683D-02/
      DATA P(217),P(218),P(219)/
     & 1.28958134880121D-02, 1.34763748338165D-02, 1.38703510891399D-02/
      DATA P(220),P(221),P(222)/
     & 1.40694249578135D-02, 2.51578703842806D-05, 1.88873264506505D-04/
      DATA P(223),P(224),P(225)/
     & 4.69184924247851D-04, 8.40571432710723D-04, 1.28438247189701D-03/
      DATA P(226),P(227),P(228)/
     & 1.78644639175865D-03, 2.33552518605716D-03, 2.92172493791781D-03/
      DATA P(229),P(230),P(231)/
     & 3.53624499771678D-03, 4.17141937698409D-03, 4.82058886485126D-03/
      DATA P(232),P(233),P(234)/
     & 5.47786669391895D-03, 6.13791528004137D-03, 6.79578550488277D-03/
      DATA P(235),P(236),P(237)/
     & 7.44682083240758D-03, 8.08660936478883D-03, 8.71096507973207D-03/
      DATA P(238),P(239),P(240)/
     & 9.31592412806942D-03, 9.89774752404876D-03, 1.04529257229060D-02/
      DATA P(241),P(242),P(243)/
     & 1.09781831526589D-02, 1.14704821146939D-02, 1.19270260530193D-02/
      DATA P(244),P(245),P(246)/
     & 1.23452623722438D-02, 1.27228849827324D-02, 1.30578366883530D-02/
      DATA P(247),P(248),P(249)/
     & 1.33483114637252D-02, 1.35927566148124D-02, 1.37898747832410D-02/
      DATA P(250),P(251),P(252)/
     & 1.39386257383068D-02, 1.40382278969086D-02, 1.40881595165083D-02/
      DATA P(253),P(254),P(255)/
     & 9.99997596379750D-01, 6.93793643241083D-06, 9.99943996207055D-01/
      DATA P(256),P(257),P(258)/
     & 5.32752936697805D-05, 9.99760490924434D-01, 1.35754910949228D-04/
      DATA P(259),P(260),P(261)/
     & 9.99380338025023D-01, 2.49212400482998D-04, 9.98745614468096D-01/
      DATA P(262),P(263),P(264)/
     & 3.89745284473282D-04, 9.97805354495956D-01, 5.54295314930373D-04/
      DATA P(265),P(266),P(267)/
     & 9.96514145914890D-01, 7.40282804244503D-04, 9.94831502800622D-01/
      DATA P(268),P(269),P(270)/
     & 9.45361516858527D-04, 9.92721344282788D-01, 1.16748411742996D-03/
      DATA P(271),P(272),P(273)/
     & 9.90151370400771D-01, 1.40490799565515D-03, 9.87092527954033D-01/
      DATA P(274),P(275),P(276)/
     & 1.65611272815445D-03, 9.83518657578632D-01, 1.91971297101387D-03/
      DATA P(277),P(278),P(279)/
     & 9.79406281670862D-01, 2.19440692536384D-03, 9.74734459752401D-01/
      DATA P(280),P(281),P(282)/
     & 2.47895822665757D-03, 9.69484659502459D-01, 2.77219576459345D-03/
      DATA P(283),P(284),P(285)/
     & 9.63640621569812D-01, 3.07301843470258D-03, 9.57188216109859D-01/
      DATA P(286),P(287),P(288)/
     & 3.38039799108691D-03, 9.50115297521293D-01, 3.69337791702565D-03/
      DATA P(289),P(290),P(291)/
     & 9.42411565191083D-01, 4.01106872407503D-03, 9.34068436157727D-01/
      DATA P(292),P(293),P(294)/
     & 4.33264096809299D-03, 9.25078932907077D-01, 4.65731729975685D-03/
      DATA P(295),P(296),P(297)/
     & 9.15437587155765D-01, 4.98436456476553D-03, 9.05140358813263D-01/
      DATA P(298),P(299),P(300)/
     & 5.31308660518706D-03, 8.94184568335557D-01, 5.64281810138445D-03/
      DATA P(301),P(302),P(303)/
     & 8.82568840247341D-01, 5.97291956550816D-03, 8.70293055548114D-01/
      DATA P(304),P(305),P(306)/
     & 6.30277344908575D-03, 8.57358310886234D-01, 6.63178124290190D-03/
      DATA P(307),P(308),P(309)/
     & 8.43766882672707D-01, 6.95936140939044D-03, 8.29522194637402D-01/
      DATA P(310),P(311),P(312)/
     & 7.28494798055382D-03, 8.14628787655138D-01, 7.60798966571904D-03/
      DATA P(313),P(314),P(315)/
     & 7.99092290960843D-01, 7.92794933429486D-03, 7.82919394118284D-01/
      DATA P(316),P(317),P(318)/
     & 8.24430376303287D-03, 7.66117819303759D-01, 8.55654356130769D-03/
      DATA P(319),P(320),P(321)/
     & 7.48696293616938D-01, 8.86417320948252D-03, 7.30664521242183D-01/
      DATA P(322),P(323),P(324)/
     & 9.16671116356077D-03, 7.12033155362253D-01, 9.46368999383007D-03/
      DATA P(325),P(326),P(327)/
     & 6.92813769779114D-01, 9.75465653631741D-03, 6.73018830230419D-01/
      DATA P(328),P(329),P(330)/
     & 1.00391720440569D-02, 6.52661665410019D-01, 1.03168123309476D-02/
      DATA P(331),P(332),P(333)/
     & 6.31756437711193D-01, 1.05871679048852D-02, 6.10318113715188D-01/
      DATA P(334),P(335),P(336)/
     & 1.08498440893373D-02, 5.88362434447664D-01, 1.11044611340069D-02/
      DATA P(337),P(338),P(339)/
     & 5.65905885423653D-01, 1.13506543159806D-02, 5.42965666498311D-01/
      DATA P(340),P(341),P(342)/
     & 1.15880740330440D-02, 5.19559661537457D-01, 1.18163858908302D-02/
      DATA P(343),P(344),P(345)/
     & 4.95706407918762D-01, 1.20352707852796D-02, 4.71425065871658D-01/
      DATA P(346),P(347),P(348)/
     & 1.22444249816120D-02, 4.46735387662029D-01, 1.24435601907140D-02/
      DATA P(349),P(350),P(351)/
     & 4.21657686626164D-01, 1.26324036435421D-02, 3.96212806057616D-01/
      DATA P(352),P(353),P(354)/
     & 1.28106981638774D-02, 3.70422087950079D-01, 1.29782022395374D-02/
      DATA P(355),P(356),P(357)/
     & 3.44307341599437D-01, 1.31346900919602D-02, 3.17890812068477D-01/
      DATA P(358),P(359),P(360)/
     & 1.32799517439305D-02, 2.91195148518247D-01, 1.34137930851101D-02/
      DATA P(361),P(362),P(363)/
     & 2.64243372410927D-01, 1.35360359349562D-02, 2.37058845589829D-01/
      DATA P(364),P(365),P(366)/
     & 1.36465181025713D-02, 2.09665238243181D-01, 1.37450934430019D-02/
      DATA P(367),P(368),P(369)/
     & 1.82086496759252D-01, 1.38316319095064D-02, 1.54346811481378D-01/
      DATA P(370),P(371),P(372)/
     & 1.39060196013255D-02, 1.26470584372302D-01, 1.39681588065169D-02/
      DATA P(373),P(374),P(375)/
     & 9.84823965981194D-02, 1.40179680394566D-02, 7.04069760428552D-02/
      DATA P(376),P(377),P(378)/
     & 1.40553820726499D-02, 4.22691647653637D-02, 1.40803519625536D-02/
      DATA P(379),P(380),P(381)/
     & 1.40938864107825D-02, 1.40928450691604D-02, 1.40944070900962D-02/
C
C
C REVISION HISTORY---
C
C FEBRUARY 1977    CHANGED DATA STATEMENTS FOR P TO MEET REQUIREMENTS OF
C                  NEW COMPILER DATA STATEMENT PROCESSOR.
C JUNE 1977        CHANGED ALL DATA STATEMENTS TO ENHANCE PORTABILITY.
C
C JANUARY 1978     DELETED REFERENCES TO THE  *COSY  CARDS
C                  AND MOVED THE REVISION HISTORIES TO APPEAR BEFORE
C                  THE FINAL END CARD
C JANUARY 1979     ADDED EXTERNAL DECLARATION OF BLOCK DATA SUBROUTINE
C                  TO MAKE SURE IT GETS LOADED FROM A BINARY LIBRARY.
C-----------------------------------------------------------------------
      END
