C --- STAGECOACH, version 2.2b
C --- This program computes various age-based life history parameters
C --- for stage-structured populations.  See Cochran and Ellner 1992,
C --- Ecological Monographs 62:345-364.
      DOUBLE PRECISION B(70,70),BJ(70),C(70,140),D(70,70),DOMEIG,
     $E(70,70),EIG(70),EIGL(70),F(70,70),G(70),P(70,70),S(70,70),T,TLX,
     $VB,WV,X(70),Y,YT(70),YV(70),Z(70),ST(70,70)
      DIMENSION IT(70,71),LF(10),LX(10),N(70)
C
C --- Dimension of some arrays may be reduced, if necessary:  minimum
C     size for B,D,E,F,P,and S is (NS,NS) or (20,20), whichever is 
C     larger--all must have the same size; for BJ,EIG,EIGL,G,N,X,YT,YV,
C     and Z is (NS); for C is (NS,2*NS); for IT is (NS,NF); and for LF 
C     and LX is (10).
C
      CHARACTER*32 ALABEL
      OPEN(8,FILE="STGCCH.dat")
      REWIND(8)
      OPEN(4,FILE="STGCCH.rst")
      REWIND(4)
      READ(8,86) ALABEL
86    FORMAT(A32)
      READ(8,*) NS,NF,NT,INTYP,NBE,MAX10,TLX,ISCRN
      
C
C --- NS=# of stage classes
C --- NF=# of fate classes (NS+1,the additional fate being death)
C --- NT=# of replicates with the same NS,NF and format
C --- INTYP determines the type of transition data
C     INTYP = 0, frequency data (FORMAT I3,19(I4)) 
C           = 1, probabilities (FORMAT 6(D11.5,1X))
C --- NBE determines if fecundities are given as "newborn equivalents"
C     NBE = 0, newborn equivalents (gamma(i))
C         = 1, raw numbers (beta(i))
C --- MAX10 = (maximum # of timesteps for stable age distribution, lx,
C              fx, and Vx/V1)/10
C --- TLX = stopping value for lx
C --- ISCRN determines if markers and warnings are printed to the sceen
C     ISCRN = 0, no output on screen
C           = 1, markers and warnings printed to screen
C --- ALABEL= 32 character label
C
      WRITE(4,112)
      IF(ISCRN.EQ.0) GO TO 95
      WRITE(*,112)
      WRITE(*,*) 'YOU HAVE CHOSEN TO MONITOR LABELS'
      WRITE(*,100) NS,NF,NT,INTYP,NBE,MAX10,TLX,ISCRN,ALABEL
   95 WRITE(4,100) NS,NF,NT,INTYP,NBE,MAX10,TLX,ISCRN,ALABEL
      NS2=NS*2
      CALL INITM(P,70,70,NS,NS)
      DO 1 K=1,NT
      IF(INTYP.EQ.1) GO TO 10
      IF(ISCRN.NE.0) WRITE(*,*) 'Survival frequencies'
      WRITE(4,*) 'Survival frequencies'
      DO 2 I=1,NS
      READ(8,*) (IT(I,J), J=1,NF)
C
C --- IT(I,J)=# of individuals moving from stage I to stage J
C
      WRITE(4,101) (IT(I,J), J=1,NF)
    2 CONTINUE
C --- Compute the transition probabilities and store in array P
      CALL PROB1(70,71,NF,IT,N,P)
      GO TO 5
   10 CONTINUE
      IF(ISCRN.NE.0) WRITE(*,*) 'Survival matrix' 
      WRITE(4,*) 'Survival matrix'
      DO 13 I=1,NS
      READ(8,*) (P(I,J), J=1,NS)
C
C --- P(I,J)=probability of survival from J to I
C
      WRITE(4,106) (P(I,J), J=1,NS)
   13 CONTINUE
C --- Input fecundities
    5 CALL INITM(B,70,70,NS,NS)
      CALL INITM(F,70,70,NS,NS)
C
C --- NBIRTH=# of transitions with nonzero fecundity
C --- NNBC=# of newborn classes
C --- N(I)=newborn class I
C --- B(I,J)=births into I from J
C
      READ(8,*) NBIRTH,NNBC
      WRITE(4,127)
      IF(ISCRN.NE.0) WRITE(*,127)
      IF(NBIRTH.EQ.0) GO TO 6
      READ(8,*) (N(I), I=1,NNBC)
      DO 21 L=1,NBIRTH
      READ(8,*) I,J,B(I,J)
      WRITE(4,117) I,J,B(I,J)
   21 CONTINUE
    6 CONTINUE
C --- Input fission
      READ(8,*) NBIRTH
C
C --- NBIRTH=# of nonzero transitions due to fission
C --- F(I,J)=expected number of stage I individuals produced by fission
C            of an average individual in stage J
C
      IF(NBIRTH.EQ.0) GO TO 35
      DO 22 L=1,NBIRTH
      READ(8,*) I,J,F(I,J)
   22 CONTINUE
      WRITE(4,128)
      IF(ISCRN.NE.0) WRITE(*,128)
C --- Determine if the C matrix is convergent.  (If not, the sum of
C     C**t does not converge.)
      DO 34 I=1,NS
      WRITE(4,106) (F(I,J), J=1,NS)
      DO 34 J=1,NS
      S(I,J)=P(I,J)+F(I,J)
      E(I,J)=S(I,J)
C --- Store P+F in F
      F(I,J)=S(I,J)
   34 CONTINUE
      CALL EIGVEC(70,NS,S,X,Z,D,0)
C --- Find the dominant eigenvalue of C
      CALL DOMINANT(70,X,Z,YV,NS,INDW) 
C --- INDC indicates if C is convergent (0=convergent; 1=not convergent)
      IF(X(INDW).GE.1.) INDC=1
      IF(INDC.EQ.1) WRITE(4,130)
      IF(INDC.EQ.1.AND.ISCRN.NE.0) WRITE(*,130)   
C --- Store the complete A matrix in S and E
      GO TO 36
   35 CALL STORE(70,70,70,70,NS,NS,F,P)
   36 DO 3 I=1,NS
      DO 3 J=1,NS
      S(I,J)=F(I,J)+B(I,J)
      E(I,J)=S(I,J)
    3 CONTINUE
      WRITE(4,129)
      IF(ISCRN.NE.0) WRITE(*,129)
      DO 44 I=1,NS
      WRITE(4,106) (S(I,J), J=1,NS)
   44 CONTINUE
C --- Find the right eigenvalues and eigenvectors of A using EISPACK.
C --- In this call, X(I)=real part of the Ith eigenvalue,
C --- Z(I)=imaginary part of the Ith eigenvalue,
C --- D(J,I)=the Jth component of the Ith eigenvector (if Z(I)=0)
C --- See the documentation of subroutine EIGVEC if Z(I) is nonzero.
      CALL EIGVEC(70,NS,S,X,Z,D,0)
      CALL DOMINANT(70,X,Z,YV,NS,INDW)
      WRITE(4,*) 'Right Eigenvalues:'
C --- Find the dominant eigenvalue and its associated eigenvector 
      WRITE(4,104) (X(I),Z(I), I=1,NS)
      WRITE(4,*) 'Right eigenvectors in rows:'
      DO 42 I=1,NS
      IF(NS.GT.6) WRITE(4,115) I,(D(J,I), J=1,NS)
      IF(NS.LE.6) WRITE(4,119) I,(D(J,I), J=1,NS)
   42 CONTINUE
C --- Normalize the eigenvector and store in EIG(J)
      CALL NORMLZ(70,NS,D,INDW,EIG)       
      DOMEIG=X(INDW)
      WRITE(4,*) 'Dominant right eigenvector = stable stage dist.'
      WRITE(4,102) DOMEIG,(EIG(J), J=1,NS)
      IF(ISCRN.NE.0) WRITE(*,102) DOMEIG,(EIG(J), J=1,NS)

C --- Find the left eigenvectors.  First replace the values in array S.
      CALL STORE(70,70,70,70,NS,NS,S,E)
      do i=1,70
	do j=1,70
		st(i,j)=s(j,i)
	end do
      end do

      CALL EIGVEC(70,NS,ST,X,Z,D,0)
C --- Find the dominant left eigenvector (reproductive value)
      CALL DOMINANT(70,X,Z,YV,NS,INDV)
      WRITE(4,*) 'Left Eigenvalues:'
      WRITE(4,104) (X(I),Z(I), I=1,NS)
      WRITE(4,*) 'Left eigenvectors in rows:'
      WV=0.0D0
      CALL NORMLZ(70,NS,D,INDV,EIGL) 
      INDG=0
      DO 11 I=1,NS
      IF(INDG.GE.1) GO TO 11
      IF(EIGL(I).NE.0) INDG=I
   11 CONTINUE
      DO 83 I=1,NS
C --- G(I) = fecundity of stage I in newborn equivalents, weighted by first
C            class with nonzero reproductive value
C --- WV = <W,V>
      WV=WV+EIG(I)*EIGL(I)
      IF(NS.GT.6) WRITE(4,115) I,(D(J,I), J=1,NS)
      IF(NS.LE.6) WRITE(4,119) I,(D(J,I), J=1,NS)
      G(I)=B(1,I)*EIGL(1)/EIGL(INDG)
      DO 83 J=2,NS
      G(I)=G(I)+B(J,I)*EIGL(J)/EIGL(INDG)
   83 CONTINUE
C --- IF NBE=1, use G(I) = unweighted fecundity
      IF(NBE.EQ.1) CALL COLSUM(B,70,70,NS,NS,G,70)
      WRITE(4,*) 'Dominant left eigenvector = reproductive value'
      WRITE(4,102) DOMEIG,(EIGL(J), J=1,NS)
C --- Calculate eigenvalue sensitivity
      WRITE(4,*) 'Sensitivity matrix, dL/daij'
      IF(ISCRN.EQ.1) WRITE(*,*) 'Sensitivity matrix, dL/daij'
      DO 93 I=1,NS
      DO 94 JJ=1,NS
      S(I,JJ)=EIGL(I)*EIG(JJ)/WV
   94 CONTINUE
      IF(NS.GT.6) WRITE(4,115) I,(S(I,J), J=1,NS)
      IF(NS.LE.6) WRITE(4,119) I,(S(I,J), J=1,NS)
   93 CONTINUE
C --- Calculate eigenvalue elasticity
      WRITE(4,*) 'Elasticity matrix'
      IF(ISCRN.EQ.1) WRITE(*,*) 'Elasticity matrix'
      DO 91 I=1,NS
      DO 92 J=1,NS
      S(I,J)=S(I,J)*E(I,J)/DOMEIG
   92 CONTINUE
      IF(NS.GT.6) WRITE(4,115) I,(S(I,J), J=1,NS)
      IF(NS.LE.6) WRITE(4,119) I,(S(I,J), J=1,NS)
   91 CONTINUE
      CALL INITV(Z,70,NS)
      T=0.0D0
C --- Calculuate the contribution of each newborn class
      DO 30 NN=1,NNBC
      Z(N(NN))=B(N(NN),1)*EIG(1)
      DO 31 J=2,NS
      Z(N(NN))=Z(N(NN))+B(N(NN),J)*EIG(J)
   31 CONTINUE
      T=T+Z(N(NN))
   30 CONTINUE 
C --- Store for future use in BJ(J)
      DO 32 J=1,NS
      BJ(J)=Z(J)/T
   32 CONTINUE
C --- B is now available to use for storage
C --- Compute age using inverse of I-C/DOMEIG.
C --- The subroutine GAUSS takes C/DOMEIG and forms the matrix I-C/DOMEIG,
C     then inverts it and returns the inverse in the second half of C
      DO 4 I=1,NS
      DO 4 J=1,NS
      C(I,J)=F(I,J)/DOMEIG
    4 CONTINUE
      CALL GAUSS(C,NS,70,140,.1E-15)
C --- Square the inverse and store it in the first half of C (S,D TEMP.)
      DO 24 I=1,NS
      DO 24 J=1,NS
      JJ=J+NS
      S(I,J)=C(I,JJ)
   24 CONTINUE
      CALL SMMULT(S,S,D,E,70,NS)
      CALL STORE(70,140,70,70,NS,NS,C,D)
C --- Cube the inverse, store it in D
      CALL SMMULT(S,D,D,E,70,NS)
C --- Premultiply INVERSE CUBED by I+(P+F)/DOMEIG, store in D
      DO 20 I=1,NS
      DO 20 J=1,NS
      S(I,J)=F(I,J)/DOMEIG
   20 CONTINUE
      DO 33 I=1,NS
      S(I,I)=S(I,I)+1.0D0
   33 CONTINUE
      CALL SMMULT(S,D,D,E,70,NS)
C --- Store YI denominators in YT(I), numerators IN Z(I),
C --- Variance numerator in X(I)
      DO 28 I=1,NS
      YT(I)=C(I,NF)*BJ(1)
      Z(I)=C(I,1)*BJ(1)
      X(I)=D(I,1)*BJ(1)
      DO 27 J=2,NS
      JJ=J+NS
      YT(I)=YT(I)+C(I,JJ)*BJ(J)
      Z(I)=Z(I)+C(I,J)*BJ(J)
      X(I)=X(I)+D(I,J)*BJ(J)
   27 CONTINUE
      Z(I)=Z(I)/YT(I)
      X(I)=X(I)/YT(I)-Z(I)**2
      IF(X(I).LT.0.0D0.AND.X(I).GT.-.5D-7) X(I)=0.0D0
      IF(X(I).GE.0.0D0) X(I)=X(I)**.5
C --- Shift S back to (P+F)/DOMEIG
      S(I,I)=S(I,I)-1.0D0
   28 CONTINUE 
      WRITE(4,111) (I,Z(I),X(I), I=1,NS)
C --- store (P+F)/DOMEIG in D
      CALL STORE(70,70,70,70,NS,NS,D,S)
C --- Compute ABAR for population (STORE IN T)
      T=Z(1)*EIG(1)*G(1)
      Y=EIG(1)*G(1)
      DO 40 I=2,NS
      T=T+Z(I)*EIG(I)*G(I)
      Y=Y+EIG(I)*G(I)
   40 CONTINUE
      T=T/Y
      IF(ISCRN.EQ.1) WRITE(*,*) 'Population ABAR'
      WRITE(4,125) T
C --- Compute the stable age distribution
C --- Denominators are already in YT(I), (P+F)/DOMEIG in D
      CALL INITM(C,70,140,NS,NS2)
      T=1.0D0
      WRITE(4,137)
      IF(ISCRN.NE.0) WRITE(*,137)
      DO 59 I=1,NS
      C(2,I)=D(I,1)*BJ(1)
      DO 59 J=2,NS
      C(2,I)=C(2,I)+D(I,J)*BJ(J)
   59 CONTINUE
      C(1,1)=BJ(1)*EIG(1)/YT(1)
      C(2,1)=C(2,1)*EIG(1)/YT(1)
      DO 65 I=2,NS
      C(1,1)=C(1,1)+BJ(I)*EIG(I)/YT(I)
      C(2,1)=C(2,1)+C(2,I)*EIG(I)/YT(I)
   65 CONTINUE
      DO 68 L=3,10
      CALL SMMULT(D,S,S,E,70,NS) 
      DO 66 I=1,NS
      C(L,I)=S(I,1)*BJ(1)
      DO 66 J=2,NS
      C(L,I)=C(L,I)+S(I,J)*BJ(J)
   66 CONTINUE
      C(L,1)=C(L,1)*EIG(1)/YT(1)
      DO 68 I=2,NS
      C(L,1)=C(L,1)+C(L,I)*EIG(I)/YT(I)
   68 CONTINUE
      WRITE(4,136) (L,C(L,1), L=1,10)
      DO 79 LL=1,10
      T=T-C(LL,1)
   79 CONTINUE
      LMAX=MAX10-1
      DO 67 L=1,LMAX
      IF(T.LT..1D-7) GO TO 67
      DO 70 LL=1,10
      LX(LL)=L*10+LL
      CALL SMMULT(D,S,S,E,70,NS) 
      DO 69 I=1,NS
      C(LL,I)=S(I,1)*BJ(1)
      DO 69 J=2,NS
      C(LL,I)=C(LL,I)+S(I,J)*BJ(J)
   69 CONTINUE
      C(LL,1)=C(LL,1)*EIG(1)/YT(1)
      DO 70 I=2,NS
      C(LL,1)=C(LL,1)+C(LL,I)*EIG(I)/YT(I)
   70 CONTINUE
      WRITE(4,136) (LX(LL),C(LL,1), LL=1,10)
      DO 78 LL=1,10
      T=T-C(LL,1)
   78 CONTINUE
   67 CONTINUE
C --- If fission occurs, SKIP MEAN REMAINING LIFESPAN AND MTFP
      IF(NBIRTH.NE.0) GO TO 25
C --- Compute mean remaining lifespan using the inverse of I-P
C --- Store matrix in C.
      CALL STORE(70,140,70,70,NS,NS,C,P)
C --- Invert I-P and store it in the second half of C
      CALL GAUSS(C,NS,70,140,.1E-15)
      WRITE(4,116)
      DO 9 J=1,NS
      JJ=J+NS
      YT(J)=C(1,JJ)
      DO 8 I=2,NS
      YT(J)=YT(J)+C(I,JJ)
    8 CONTINUE
    9 CONTINUE
C --- Compute the variance.
C --- Square the inverse and store it in the first half of C
      DO 43 I=1,NS
      DO 43 J=1,NS
      JJ=J+NS
      D(I,J)=C(I,JJ)
   43 CONTINUE
      CALL SMMULT(D,D,D,E,70,NS)
      CALL STORE(70,140,70,70,NS,NS,C,D)
C --- Premultiply by I+P, STORE IN D
      CALL STORE(70,70,70,70,NS,NS,S,P)
      DO 41 I=1,NS
      S(I,I)=S(I,I)+1.0D0
   41 CONTINUE
      CALL SMMULT(S,D,D,E,70,NS)
C --- Sum over I
      DO 50 J=1,NS
      YV(J)=D(1,J)
      DO 49 I=2,NS
      YV(J)=YV(J)+D(I,J)
   49 CONTINUE
      YV(J)=YV(J)-YT(J)**2
      IF(YV(J).GT.0.0D0) YV(J)=YV(J)**.5
      WRITE(4,114) J,YT(J),YV(J)
   50 CONTINUE
C --- Compute mean time to first passage from J to II
      DO 47 II=1,NS
      WRITE(4,105) II
C --- Compute the matrix DI 
      CALL STORE(70,70,70,70,NS,NS,D,P)
      DO 48 I=1,NS
      DO 48 J=1,NS
      IF(J.EQ.II) D(I,J)=0.0D0
   48 CONTINUE
C --- Compute I-DI, store in C and invert into last half of C.
      CALL STORE(70,140,70,70,NS,NS,C,D)
      CALL GAUSS(C,NS,70,140,.1E-15)
C --- Square the inverse and store in B, then first half of C.
C --- Cube the inverse and store in S.
      DO 56 I=1,NS
      DO 56 J=1,NS
      JJ=J+NS
      S(I,J)=C(I,JJ)
   56 CONTINUE
      CALL SMMULT(S,S,B,E,70,NS)
      CALL SMMULT(S,B,S,E,70,NS)
      CALL STORE(70,140,70,70,NS,NS,C,B)
      CALL STORE(70,70,70,70,NS,NS,B,D)
      DO 57 I=1,NS
C --- Compute I+DI, store in B
      B(I,I)=B(I,I)+1.0D0
   57 CONTINUE
C --- Premultiply INV(I-DI) CUBED TIMES I+DI, STORE IN S
      CALL SMMULT(B,S,S,E,70,NS)
C --- Store mean in X, variance numerator in Z 
      DO 62 J=1,NS
      JJ=J+NS
C --- If stage II is unreachable from stage J, X(J) = -1
      X(J)=-1.0D0 
      IF(C(II,JJ).EQ.0.0D0) GO TO 62
      X(J)=X(J)+C(II,J)/C(II,JJ)
      Z(J)=(S(II,J)-(C(II,J)**2/C(II,JJ)))/C(II,JJ)
      IF(X(J).LE.0.0D0) Z(J)=0.0D0
      IF(Z(J).LT.0.0D0) Z(J)=0.0D0
      IF(Z(J).GT.0.0D0) Z(J)=Z(J)**.5
   62 CONTINUE
      WRITE(4,103) (X(J), J=1,NS)
      WRITE(4,107) (Z(J), J=1,NS)
C --- Compute lifespan (X(J)) and variance (Z(J))
      DO 76 J=1,NS
      IF(X(J).EQ.-1.0D0) GO TO 76
      X(J)=1.0D0+X(J)+YT(II)
      IF(Z(J).GE.0.0D0.AND.YV(II).GE.0.0D0) Z(J)=(Z(J)**2+YV(II)**2)**.5
   76 CONTINUE
      WRITE(4,109) (X(J), J=1,NS)
      WRITE(4,107) (Z(J), J=1,NS)
   47 CONTINUE
   25 CONTINUE
C --- Compute Lx, store in D (rows 1-10); fx, store in D (rows 11-20);
C     Vx, store in E (rows 1-10); population Lx, store in E (row 11);
C     population fx, (row 12); and population Vx, store in E (row 13).
C --- Compute fx in terms of total number of newborns (weighted for
C     reproductive value, NBE=0, or unweighted, NBE=1).
C --- P+F is stored in F; (B ISN'T USED ANYMORE)
      CALL STORE(70,70,70,70,NS,NS,S,F)
      WRITE(4,135)
      IF(ISCRN.NE.0) WRITE(*,135)
C --- Scale the left eigenvector (reproductive value)
      VB=EIGL(1)*BJ(1)
      DO 45 J=2,NS
      VB=VB+EIGL(J)*BJ(J)
   45 CONTINUE
      DO 55 J=1,NS
      EIGL(J)=EIGL(J)/VB
   55 CONTINUE
      E(11,1)=1.0D0
      E(12,1)=0.0D0
      E(13,1)=1.0D0
      E(11,2)=0.0D0
      E(12,2)=0.0D0
      E(13,2)=0.0D0
      DO 26 J=1,NNBC
      LX(1)=1
      D(1,N(J))=1.0D0
      LF(1)=11
      D(11,N(J))=G(N(J))
      E(1,N(J))=EIGL(N(J))
      LX(2)=2
      D(2,N(J))=S(1,N(J))
      LF(2)=12
      D(12,N(J))=S(1,N(J))*G(1)
      E(2,N(J))=S(1,N(J))*EIGL(1)
      DO 81 I=2,NS
      D(2,N(J))=D(2,N(J))+S(I,N(J))
      D(12,N(J))=D(12,N(J))+S(I,N(J))*G(I)
      E(2,N(J))=E(2,N(J))+S(I,N(J))*EIGL(I)
   81 CONTINUE
      E(13,2)=E(13,2)+E(2,N(J))*BJ(N(J))
      D(12,N(J))=D(12,N(J))/D(2,N(J))
      E(2,N(J))=E(2,N(J))/D(2,N(J))
      E(12,1)=E(12,1)+D(11,N(J))*BJ(J)
      E(11,2)=E(11,2)+D(2,N(J))*BJ(J)
      E(12,2)=E(12,2)+D(12,N(J))*BJ(J)
   26 CONTINUE
      E(13,2)=E(13,2)/E(11,2)
      DO 39 L=3,10
      E(11,L)=0.0D0
      E(12,L)=0.0D0
      E(13,L)=0.0D0
      LX(L)=L
      LF(L)=L+10
      CALL SMMULT(F,S,S,B,70,NS) 
      DO 39 J=1,NNBC
      D(L,N(J))=S(1,N(J))
      D(LF(L),N(J))=S(1,N(J))*G(1)
      E(L,N(J))=S(1,N(J))*EIGL(1)
      DO 75 I=2,NS
      D(L,N(J))=D(L,N(J))+S(I,N(J))
      D(LF(L),N(J))=D(LF(L),N(J))+S(I,N(J))*G(I)
      E(L,N(J))=E(L,N(J))+S(I,N(J))*EIGL(I)
   75 CONTINUE
      E(11,L)=E(11,L)+D(L,N(J))*BJ(J)
      E(13,L)=E(13,L)+E(L,N(J))*BJ(N(J))
      D(LF(L),N(J))=D(LF(L),N(J))/D(L,N(J))
      E(L,N(J))=E(L,N(J))/D(L,N(J))
      E(12,L)=E(12,L)+D(LF(L),N(J))*BJ(J)
   39 CONTINUE
      J=1
   77 J1=J+1
      IF(J1.GT.NNBC)  J1=J
      WRITE(4,143) (N(JJ),LX(1),D(1,N(JJ)),D(LF(1),N(JJ)),E(1,N(JJ)),
     $JJ=J,J1)
      DO 15 L=2,10
      WRITE(4,144) (LX(L),D(L,N(JJ)),D(LF(L),N(JJ)),E(L,N(JJ)), JJ=J,J1)
   15 CONTINUE
      J=J+2
      IF(J.LE.NNBC) GO TO 77
      T=E(11,10)
      IF(NNBC.EQ.1) GO TO 73
      WRITE(4,*) 'Population   lx          fx        Vx/V1'
      DO 82 L=1,10
      E(13,L)=E(13,L)/E(11,L)
   82 CONTINUE
      WRITE(4,145) (LX(L),E(11,L),E(12,L),E(13,L), L=1,10)
   73 MAX10=MAX10-1
      DO 38 L=1,MAX10
      IF(T.LT.TLX) GO TO 38
      DO 54 LL=1,10
      LX(LL)=L*10+LL
      CALL SMMULT(F,S,S,B,70,NS) 
      E(11,LL)=0.0D0
      E(12,LL)=0.0D0
      E(13,LL)=0.0D0
      DO 54 J=1,NNBC
      D(LL,N(J))=S(1,N(J))
      D(LF(LL),N(J))=S(1,N(J))*G(1)
      E(LL,N(J))=S(1,N(J))*EIGL(1)
      DO 89 I=2,NS
      D(LL,N(J))=D(LL,N(J))+S(I,N(J))
      D(LF(LL),N(J))=D(LF(LL),N(J))+S(I,N(J))*G(I)
      E(LL,N(J))=E(LL,N(J))+S(I,N(J))*EIGL(I)
   89 CONTINUE
      E(11,LL)=E(11,LL)+D(LL,N(J))*BJ(J)
      E(13,LL)=E(13,LL)+E(LL,N(J))*BJ(J)
      D(LF(LL),N(J))=D(LF(LL),N(J))/D(LL,N(J))
      E(LL,N(J))=E(LL,N(J))/D(LL,N(J))
      E(12,LL)=E(12,LL)+D(LF(LL),N(J))*BJ(J) 
   54 CONTINUE
      IF(NNBC.GT.1) WRITE(4,135)
      J=1
   53 J1=J+1
      IF(J1.GT.NNBC)  J1=J
      WRITE(4,143) (N(JJ),LX(1),D(1,N(JJ)),D(LF(1),N(JJ)),E(1,N(JJ)),
     $JJ=J,J1)
      DO 37 LLL=2,10
      WRITE(4,144) (LX(LLL),D(LLL,N(J)),D(LF(LLL),N(J)),E(LLL,N(J)),
     $JJ=J,J1)
   37 CONTINUE
      J=J+2
      IF(J.LE.NNBC) GO TO 53
      T=E(11,10)
      IF(NNBC.EQ.1) GO TO 38
      WRITE(4,*) 'Population   lx          fx        Vx/V1'
      DO 84 LLL=1,10
      E(13,LLL)=E(13,LLL)/E(11,LLL)
   84 CONTINUE
      WRITE(4,145) (LX(LLL),E(11,LLL),E(12,LLL),E(13,LLL), LLL=1,10)
   38 CONTINUE
C --- Compute Net Reproductive Rate, R0
C --- Store (I-C)=(I-(P+F)) in the first part of C
      IF(INDC.EQ.1) WRITE(*,142)
      IF(INDC.EQ.1) WRITE(4,142)
      IF(INDC.EQ.1) GO TO 88
      CALL STORE(70,140,70,70,NS,NS,C,F)
C --- Invert (stored in the last half of C)
      CALL GAUSS(C,NS,70,140,.1E-15)
C --- Square the inverse and store in the first half of C
      DO 12 I=1,NS
      DO 12 J=1,NS
      JJ=J+NS
      S(I,J)=C(I,JJ)
   12 CONTINUE
      CALL SMMULT(S,S,D,B,70,NS)
      CALL STORE(70,140,70,70,NS,NS,C,D)
C --- Cube the inverse, store in D
      CALL SMMULT(S,D,D,B,70,NS)
      CALL STORE(70,70,70,70,NS,NS,S,F)
      DO 90 I=1,NS
      S(I,I)=S(I,I)+1.0D0
   90 CONTINUE
C --- Premultiply inverse cubed times I+P+F
      CALL SMMULT(S,D,S,B,70,NS)
      Y=0.0D0
      T=0.0D0
      Z(1)=0.0D0
C --- Store RO(J) in YV(J)
C --- Compute the means and variances for MU1
C --- Compute MU1 for each class (STORE IN YT(N(J)))
C --- Multiply by G(I), store results in X
      IF(ISCRN.EQ.1) WRITE(*,*) 'Computing R0 and Mu1'
      DO 17 J=1,NNBC
      JJ=N(J)+NS
      YV(N(J))=C(1,JJ)*G(1)
      X(N(J))=S(1,N(J))*G(1)
      YT(N(J))=C(1,N(J))*G(1)
      DO 23 I=2,NS
      YV(N(J))=YV(N(J))+C(I,JJ)*G(I)
      X(N(J))=X(N(J))+S(I,N(J))*G(I)
      YT(N(J))=YT(N(J))+C(I,N(J))*G(I)
   23 CONTINUE
C --- Compute R0 for population,MU1 numerator for population
      Y=Y+YV(N(J))*BJ(N(J))
      Z(1)=Z(1)+X(N(J))*BJ(N(J))
C --- Compute MU1 for population (STORE IN T)
      T=T+YT(N(J))*BJ(N(J))
      YT(N(J))=YT(N(J))/YV(N(J))
      X(N(J))=(X(N(J))/YV(N(J))-YT(N(J))**2)
      IF(X(N(J)).LT.0.0D0.AND.X(N(J)).GT.-.5D-7) X(N(J))=0.0D0
      IF(X(N(J)).GT.0.0D0) X(N(J))=X(N(J))**.5
      WRITE(4,120) N(J),YV(N(J))
      WRITE(4,123) YT(N(J)),X(N(J))
   17 CONTINUE
      IF(NNBC.EQ.1) GO TO 85
      WRITE(4,121) Y
      IF(ISCRN.EQ.2) WRITE(*,121) Y
C --- Finish computing population MU1 (STORE IN T)
      T=T/Y
C --- Divide all by appropriate R0'S
      Z(1)=((Z(1)/Y)-T**2)
      IF(Z(1).LT.0.0D0.AND.Z(1).GT.-.5D-7) Z(1)=0.0D0
      IF(Z(1).GT.0.0D0) Z(1)=Z(1)**.5
      WRITE(4,123) T,Z(1)
C --- Compute population T
      Y=DLOG(Y)/DLOG(DOMEIG)
      WRITE(4,122) Y
C --- Compute age of residence (NUM IN FIRST HALF OF C, DENOM
C     IN LAST HALF OF C, NUM OF VAR IN S)
C --- Store population mean numerator in YT, denominator in T, 
C     variance numerator in YV
   85 CALL INITV(YT,70,NS)
      CALL INITV(YV,70,NS)
      DO 64 I=1,NS
      T=0.0D0
      DO 61 J=1,NNBC
      JJ=N(J)+NS
C --- If the denominator is 0, so are both numerators
      IF(C(I,JJ).EQ.0.0D0) GO TO 61
      T=T+C(I,JJ)*BJ(N(J))
      YT(I)=YT(I)+C(I,N(J))*BJ(N(J))
      YV(I)=YV(I)+S(I,N(J))*BJ(N(J)) 
      C(I,N(J))=C(I,N(J))/C(I,JJ)
      S(I,N(J))=S(I,N(J))/C(I,JJ)-C(I,N(J))**2
   61 CONTINUE
      YT(I)=YT(I)/T
      YV(I)=((YV(I)/T)**.5-YT(I))*((YV(I)/T)**.5+YT(I))
      IF(YV(I).LT.0.0D0.AND.YV(I).GT.-.5D-7) YV(I)=0.0D0
      IF(YV(I).GT.0.0D0) YV(I)=YV(I)**.5
   64 CONTINUE
C --- Finalize population variance.
      DO 63 I=1,NS
      DO 74 J=1,NNBC
      IF(S(I,N(J)).GE.0.0D0) S(I,N(J))=S(I,N(J))**.5
      IF(S(I,N(J)).LT.0.0D0.AND.S(I,N(J)).GT.-.5D-7) S(I,N(J))=0.0D0
   74 CONTINUE
      WRITE(4,133) I
      WRITE(4,131) (N(J),C(I,N(J)),S(I,N(J)), J=1,NNBC) 
      IF(NNBC.GT.1) WRITE(4,132) YT(I),YV(I)
   63 CONTINUE
C --- Compute mean time to first passage
   88 DO 46 J=1,NS
      IF(P(J,J).EQ.1.0D0) X(J)=1.0D3
      IF(P(J,J).EQ.1.0D0) GO TO 46
      X(J)=1.0D0/(1.0D0-P(J,J))
   46 CONTINUE
      WRITE(4,110) (J,X(J), J=1,NS)
C --- Recall F=P+F
C --- Compute age at maturity for each newborn class
C --- Compute I-Q, store in first half of C
      DO 58 M=1,NS
      DO 58 I=1,NS
      IF(G(M).GT.0.0D0) F(I,M)=0.0D0
      C(I,M)=F(I,M)
   58 CONTINUE 
C --- Invert
      CALL GAUSS(C,NS,70,140,.1E-15)
C --- Square (infirst half of C)
      DO 97 I=1,NS
      DO 97 J=1,NS
      JJ=J+NS
      S(I,J)=C(I,JJ)
   97 CONTINUE
      CALL SMMULT(S,S,D,B,70,NS)
      CALL STORE(70,140,70,70,NS,NS,C,D)
C --- Cube, store in D, premultiply by I + Q, store in S
      CALL SMMULT(S,D,D,B,70,NS)
      DO 52 I=1,NS
      CALL STORE(70,70,70,70,NS,NS,S,F)
      S(I,I)=S(I,I)+1.0D0
   52 CONTINUE
      CALL SMMULT(S,D,S,B,70,NS)
C --- Numerator in X, denominator in Z, variance numerator in YV
      DO 96 J=1,NNBC
      JJ=N(J)+NS
      X(N(J))=0.0D0
      Z(N(J))=0.0D0
      YV(N(J))=0.0D0
      DO 96 I=1,NS
      IF(G(I).EQ.0.0D0) GO TO 96
      X(N(J))=C(I,N(J))+X(N(J))
      Z(N(J))=C(I,JJ)+Z(N(J))
      YV(N(J))=S(I,N(J))+YV(N(J))
   96 CONTINUE
      T=0.0D0
      Y=0.0D0
      DO 51 J=1,NNBC
      T=T+X(N(J))*BJ(J)
      Y=Y+Z(N(J))*BJ(J)
      X(N(J))=X(N(J))/Z(N(J))
      YV(N(J))=(YV(N(J))/Z(N(J))-(X(N(J)))**2)
      IF(YV(N(J)).GE.0.0D0) YV(N(J))=YV(N(J))**.5 
C --- Probability of reaching maturity = Z
   51 CONTINUE
      WRITE(4,108) (N(J),Z(N(J)),X(N(J)),YV(N(J)), J=1,NNBC)
C --- Compute population age at maturity
      T=T/Y
      IF(NNBC.GT.1) WRITE(4,124) Y,T
  100 FORMAT(3(I2,1X),2(I1,1X),I3,1X,D11.5,1X,I1,1X,A32)
  101 FORMAT(4(I3,19(I4),/))
  102 FORMAT(/,'Dominant Eigenvalue= ',D11.5,/,'Eigenvector= ',
     $6(D10.5,1X),11(/,13X,6(D10.5,1X)))
  103 FORMAT('Mean time to I',2X,70(D11.5,1X))
  104 FORMAT(3('(',D11.5,' + ',D11.5,' I)',1X))
  105 FORMAT(/,'Target stage I=',I2)
  106 FORMAT(6(D11.5,1X))
  107 FORMAT('Standard dev. ',2X,6(D11.5,1X)/(14X,6(D11.5,1X)))
  108 FORMAT(/,('Newborn type ',I2,' Pr(maturity) = ',D11.5,
     $' Age(maturity) = ',D11.5,/,48X,'St. dev. = ',D11.5))
  109 FORMAT('Total lifespan',2X,16(D11.5,1X))
  110 FORMAT(/,'Mean time to first passage',/,5('st.',4X,'time',5X),/,
     $(5(I2,1X,D11.5,2X)))
  111 FORMAT(/,('Mean age in stage ',I2,' is ',D13.6,' with st. dev. ',
     $D13.6))
  112 FORMAT('Welcome to STAGECOACH, ver. 2.2 by M. E. Cochran.',/,
     $'(1992. ESA Supplementary Publication Service Document No. 9202.'
     $,/,'The Ecological Society of America, 328 East State Street,',
     $/,'Ithaca, NY  14850-4318 USA.)',/,'For further details, see: ',
     $/,'Cochran and Ellner, 1992. Ecol. Monog. 62:345-364.',/,/)
  113 FORMAT(/,3(I2,1X),2(I1,1X),A32)
  114 FORMAT(2X,I2,2X,D11.5,2X,D11.5,2X,I4,2X,D11.5,2X,I4)
  115 FORMAT(1X,I2,1X,5(D11.5,2x)/(4X,5(D11.5,2x)))
  116 FORMAT(/,' St.  Mean lifespan   St dev')
  117 FORMAT(I2,1X,I2,1X,D16.8)
  118 FORMAT(16(I2,1X))
  119 FORMAT(1X,I2,1X,6(D11.5,2x))
  120 FORMAT(/,'Net reproductive rate for stage ',I2,' = ',
     $D11.5)
  121 FORMAT(/,'population net reproductive rate = ',D11.5)
  122 FORMAT('Generation time T = ',D11.5)
  123 FORMAT('Generation time MU1 = ',D11.5,3X,'St.dev. = ',D11.5)
  124 FORMAT(/,'pop. Pr(maturity) = ',D11.5,5X,'pop. Age(maturity) = ',
     $D11.5)
  125 FORMAT(/,'population generation time ABAR = ',D11.5)
  126 FORMAT('F(2,1)=',D16.10)
  127 FORMAT(/,'Birth matrix')
  128 FORMAT(/,'Fission matrix')
  129 FORMAT(/,'Transition matrix')
  130 FORMAT(/,'Warning:  THE C MATRIX IS NOT CONVERGENT')
  131 FORMAT('Newborn type ',I2,2X,'Mean age = ',D11.5,3X,
     $'St.dev. = ',D11.5)
  132 FORMAT(10X,'Population mean = ',D11.5,3X,'St.dev. = ',D11.5)
  133 FORMAT(/,'Residence ages in state ',I2)
  135 FORMAT(/,2('type',2X,'X',7X,'lx',9x,'fx',9x,'Vx/V1',4x))
  136 FORMAT(5(I4,1X,D11.5,1X)/(5(I4,1X,D11.5,1X)))
  137 FORMAT(/,'Stable Age Distribution')
  142 FORMAT('Since C is not convergent, R0 and Mu1 do not converge.')
  143 FORMAT((2(1x,I2,1X,I4,3(1X,D11.5))))
  144 FORMAT((2(4X,I4,3(1X,D11.5))))
  145 FORMAT((4X,I4,3(1X,D11.5)))
    1 CONTINUE
      STOP
      END


      SUBROUTINE   EIGVEC(NM,N,A,WR,WI,Z,MODE)
      INTEGER      NM,N,LOW,IGH,IERR,MODE
      DOUBLE PRECISION         A(NM,N),Z(NM,N), WR(N),WI(N)
      DOUBLE PRECISION         SCALE(100 )
COMMENT*XXXX;
C**PURPOSE     FINDS ALL N EIGENVALUES AND NORMALIZED EIGENVECTORS OF
C   THE GENERAL DOUBLE PRECISION MATRIX A.
COMMENT*YYYY;
C INPUT       NM - ROW DIMENSION OF THE TW0-DIM. ARRAYS A AND Z AS
C                  DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT
C             N  - ORDER OF MATRIX
C             A  - CONTAINS THE GENERAL DOUBLE PRECISION INPUT MATRIX.
C             MODE - DETERMINES TYPE OF NORMALIZATION FOR THE EIGEN-
C                    VECTORS.
C OUTPUT      A  - THE INPUT MATRIX IS DESTROYED BY THE COMPUTATION
C             WR - THE REAL PARTS OF THE EIGENVALUES ARE STORED IN THE
C                  ARRAY WR
C             WI - THE IMAGINARY PARTS OF THE EIGENVALUES ARE STORED
C                  IN THE ARRAY WI.  COMPLEX CONJUGATE EIGENVALUES
C                  ARE STORED CONSECUTIVELY.
C             Z  - TW0-DIMENSIONAL ARRAY CONTAINING THE NORMALIZED
C                  EIGENVECTORS OF A IN ITS FIRST N COLUMNS: IF THE
C                  K-TH EIGENVALUE IS REAL, COLUMN K CONTAINS ITS
C                  EIGENVECTOR; IF THE K-TH EIGENVALUE IS COMPLEX WITH
C                  POSITIVE IMAGINARY PART, COLUMN K CONTAINS THE REAL
C                  PART AND COLUMN K+1 THE IMAGINARY PART OF THE KTH
C                  EIGENVECTOR. THE COMPLEX CONJUGATE EIGENVECTOR IS
C                  NOT STORED.
C
C SUBROUTINES REQUIRED:  BALANC, ORTHES, ORTRAN, HQR2, BALBAK, NRMLZ
C
C RESTRICTION  THE ORDER N OF THE MATRIX MUST NOT EXCEED 100.  FOR
C   LARGER MATRICES, THE DIMENSION OF THE WORKING STORAGE ARRAY
C   SCALE MUST BE INCREASED.  THE ARRAYS A AND Z  MUST HAVE THE SAME
C   ACTUAL ROW DIMENSION AND MUST OCCUPY SEPARATE STORAGE.
C
C METHOD      A IS FIRST PREPARED BY SUBROUTINE BALANC, WHICH DETECTS
C   ANY ISOLATED EIGENVALUES AND REDUCES THE NORM OF A WITHOUT INTRO-
C   DUCING ANY ROUNDOFF ERRORS.  ORTHES IS THEN CALLED TO REDUCE A TO
C   UPPER HESSENBERG FORM USING ELEMENTARY ORTHOGONAL (HOUSEHOLDER)
C   TRANSFORMATIONS.  THE MATRIX DEFINING THE SIMILARITY REDUCTION OF
C   A TO UPPER HESSENBERG FORM IS THEN FORMED AND STORED IN THE ARRAY
C   U BY SUBROUTINE ORTRANS USING THE INFORMATION LEFT BELOW THE SUB-
C   DIAGONAL OF A BY ORTHES, COMPLETING THE PRELIMINARY REDUCTIONS AND
C   TRANSFORMATIONS.  SUBROUTINE HQR2 IS EMPLOYED TO OBTAIN THE EIGEN-
C   VALUES OF A AND THE EIGENVECTORS OF THE BALANCED MATRIX.  HQR2
C   UTILIZES THE FRANCIS DOUBLE-STEP QR ALGORITHM TO COMPUTE THE
C   EIGENVALUES AND A METHOD OF BACK SUBSTITUTION FOR THE EIGENVECTORS.
C   THE EIGENVECTORS OF THE BALANCED MATRIX ARE NEXT TRANSFORMED INTO
C   EIGENVECTORS OF THE ORIGINAL MATRIX BY SUBROUTINE BALBAK.  FINALLY,
C   THE EIGENVECTORS OF A ARE NORMALIZED SO THAT A REAL EIGENVECTOR
C   HAS ITS COMPONENT OF MAXIMUM MAGNITUDE EQUAL TO 1.0, AND A COMPLEX
C   EIGENVECTOR HAS ITS MAXIMUM MODULUS COMPONENT EQUAL TO 1.0 + 0*I.
C
C REMARKS     BALANCING IS STRONGLY RECOMMENDED SINCE IT GENERALLY
C   IMPROVES ACCURACY AT A SMALL ADDITIONAL COST IN TIME.  HOWEVER, IT
C   IS NOT NECESSARY TO BALANCE ORTHOGONAL OR NORMAL MATRICES.  IF
C   BALANC IS NOT USED, LOW MUST BE SET TO 1 AND UPP TO N.
C             IF THE ORIGINAL MATRIX IS ALREADY OF UPPER HESSENBERG
C   FORM (ZEROES BELOW THE FIRST SUB-DIAGONAL), THE CALLS TO ORTHES
C   AND ORTRANS MAY BE REMOVED (SUBROUTINE BALANC WILL PRESERVE THE
C   HESSENBERG FORM).  IN THIS CASE U MUST BE SET TO THE IDENTITY
C   MATRIX BEFORE HQR2 IS ENTERED.
C             IF UNORMALIZED EIGENVECTORS ARE SUFFICIENT, THE CALL TO
C   SUBROUTINE NORMALIZE MAY BE REMOVED.
C             IF THE K-TH EIGENVALUE HAS NOT BEEN DETERMINED AFTER
C  30 ITERATIONS, THE VALUE OF K WILL BE PRINTED AND THE SUBROUTINE
C   EXITED NORMALLY WITH THE EIGENVALUES OF INDICES K+1, K+2, ..., N
C   STORED CORRECTLY.  FAILURES WILL BE EXCEEDINGLY RARE.
C
C REFERENCE   FOR MORE INFORMATION REFER TO G. PETERS AND J. H.
C   WILKINSON, "EIGENVECTORS OF REAL AND COMPLEX MATRICES BY LR AND
C   QR TRIANGULARIZATIONS," NUMERISCHE MATHEMATIK 16, 1970, PP 181-204.
C
C ACKNOWLEDGEMNT  THESE FORTRAN TRANSLATIONS OF BALANC, ORTHES,
C   ORTRAN, HQR2 AND BALBAK ARE PART OF THE EISPACK SERIES OF
C   SUBROUTINES DEVELOPED AND MADE AVAILABLE BY  NATS PROJECT OF THE
C   THE APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORIES,
C   DIVISION, ARGONNE NATIONAL LABORATORIES, ARGONNE, ILLINOIS.
C
COMMENT*ZZZZ;
      CALL BALANC(NM,N,A,LOW,IGH,SCALE)
      CALL ORTHES(NM,N,LOW,IGH,A,WR)
C     ARRAY WR IS USED TEMPORARILY TO PASS INFO FROM ORTHES TO ORTRAN
      CALL ORTRAN(NM,N,LOW,IGH,A,WR,Z)
      CALL HQR2(NM,N,LOW,IGH,A,WR,WI,Z,IERR)
      CALL BALBAK(NM,N,LOW,IGH,SCALE,N,Z)
      CALL NRMLZ(NM,N,WI,Z,MODE)
      IF (IERR .NE. 0) WRITE(*,1) IERR
    1 FORMAT(' ',I3,'-TH EIGENVALUE WAS NOT FOUND AFTER 3O ITERATIONS')
      RETURN
      END
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      DOUBLE PRECISION A(NM,N),SCALE(N)
      DOUBLE PRECISION C,F,G,R,S,B2,RADIX
      DOUBLE PRECISION DABS
      LOGICAL NOCONV
COMMENT*XXXX;
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     BALANC   BALANCES A DOUBLE PRECISION MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
COMMENT*YYYY;
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT-
C
C        A CONTAINS THE BALANCED MATRIX,
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     -----------------------------------------------------------------
C
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
C
C                **********
COMMENT*ZZZZ;
      RADIX = 2.
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     ********** IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE **********
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN **********
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT **********
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
      DO 180 I = K, L
  180 SCALE(I) = 1.0
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0
         R = 0.0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + DABS(A(J,I))
            R = R + DABS(A(I,J))
  200    CONTINUE
C
         G = R / RADIX
         F = 1.0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     ********** NOW BALANCE **********
  240    IF ((C + R) / F .GE. 0.95 * S) GO TO 270
         G = 1.0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
C
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
C**********************************************************************
C
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      DOUBLE PRECISION A(NM,N),ORT(IGH)
      DOUBLE PRECISION F,G,H,SCALE
      DOUBLE PRECISION DSQRT,DABS,DSIGN
COMMENT*XXXX;
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A DOUBLE PRECISION GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
COMMENT*YYYY;
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N,
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT-
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX,
C
C          ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     -----------------------------------------------------------------
C
COMMENT*ZZZZ;
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0
         ORT(M) = 0.0
         SCALE = 0.0
C     ********** SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) **********
         DO 90 I = M, IGH
   90    SCALE = SCALE + DABS(A(I,M-1))
C
         IF (SCALE .EQ. 0.0) GO TO 180
         MP = M + IGH
C     ********** FOR I=IGH STEP -1 UNTIL M DO -- **********
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
C
         G = -DSIGN(DSQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
C     ********** FORM (I-(U*UT)/H) * A **********
         DO 130 J = M, N
            F = 0.0
C     ********** FOR I=IGH STEP -1 UNTIL M DO -- **********
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
C
            F = F / H
C
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
C
  130    CONTINUE
C     ********** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) **********
         DO 160 I = 1, IGH
            F = 0.0
C     ********** FOR J=IGH STEP -1 UNTIL M DO -- **********
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
C
            F = F / H
C
          DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
C
  160    CONTINUE
C
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
C
  200 RETURN
      END
C***********************************************************************
C
      SUBROUTINE ORTRAN(NM,N,LOW,IGH,A,ORT,Z)
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      DOUBLE PRECISION A(NM,IGH),ORT(IGH),Z(NM,N)
      DOUBLE PRECISION G,H
COMMENT*XXXX;
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE ORTRANS,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     ORTRAN   ACCUMULATES THE ORTHOGONAL SIMILARITY
C     TRANSFORMATIONS USED IN THE REDUCTION OF A DOUBLE PRECISION GENERAL
C     MATRIX TO UPPER HESSENBERG FORM BY  ORTHES.
COMMENT*YYYY;
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N,
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES
C          IN ITS STRICT LOWER TRIANGLE,
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  ORTHES,
C
C        ORT HAS BEEN ALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     -----------------------------------------------------------------
C
COMMENT*ZZZZ;
C     ********** INITIALIZE Z TO IDENTITY MATRIX **********
      DO 80 I = 1, N
C
         DO 60 J = 1, N
   60    Z(I,J) = 0.0
C
         Z(I,I) = 1.0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     ********** FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********
      DO 140 MM = 1, KL
         MP = IGH - MM
C     ********** H BELOW IS NEGATIVE OF H FORMED IN ORTHES **********
         H = A(MP,MP-1) * ORT(MP)
         IF (H .EQ. 0.0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = MP, IGH
            G = 0.0
C
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
C
            G = G / H
C
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
C***********************************************************************
C
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITS,LOW,MP2,ENM2,IERR
      DOUBLE PRECISION H(NM,N),WR(N),WI(N),Z(NM,N)
      DOUBLE PRECISION P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
     # ,TEMP,TA,TB
      DOUBLE PRECISION DSQRT,DABS,DSIGN
      INTEGER MIN0
      LOGICAL NOTLAS
COMMENT*XXXX;
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     HQR2   FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A DOUBLE PRECISION UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C    EIGENVECTORS OF A DOUBLE PRECISION GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
COMMENT*YYYY;
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N,
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        H HAS BEEN DESTROYED,
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N,
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL
C     PROCEDURE CDIV BY COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     -----------------------------------------------------------------
C
COMMENT*ZZZZ;
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP = 2D+00**(-56)
C
      IERR = 0
C     ********** STORE ROOTS ISOLATED BY BALANC **********
      DO 50 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0
   50 CONTINUE
C
      EN = IGH
      T = 0.0
C     ********** SEARCH FOR NEXT EIGENVALUES **********
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- **********
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         IF (DABS(H(L,L-1)) .LE. MACHEP * (DABS(H(L-1,L-1))
     X      + DABS(H(L,L)))) GO TO 100
   80 CONTINUE
C     ********** FORM SHIFT **********
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     ********** FORM EXCEPTIONAL SHIFT **********
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
      X = 0.75 * S
      Y = X
      W = -0.4375 * S * S
  130 ITS = ITS + 1
C     ********** LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- **********
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         IF (DABS(H(M,M-1)) * (DABS(Q) + DABS(R)) .LE. MACHEP * DABS(P)
     X    * (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0
       IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0
  160 CONTINUE
C     ********** DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN **********
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0
         IF (NOTLAS) R = H(K+2,K-1)
         X = DABS(P) + DABS(Q) + DABS(R)
         IF (X .EQ. 0.0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     ********** ROW MODIFICATION **********
         DO 210 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C
         J = MIN0(EN,K+3)
C     ********** COLUMN MODIFICATION **********
         DO 230 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C     ********** ACCUMULATE TRANSFORMATIONS **********
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            IF (.NOT. NOTLAS) GO TO 240
            P = P + ZZ * Z(I,K+2)
          Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K) = Z(I,K) - P
  250    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     ********** ONE ROOT FOUND **********
  270 H(EN,EN) = X  +  T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0
      EN = NA
      GO TO 60
C     ********** TWO ROOTS FOUND **********
  280 P = (Y - X) / 2.0
      Q = P * P + W
      ZZ = DSQRT(DABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0) GO TO 320
C     ********** REAL PAIR **********
      ZZ = P + DSIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0
      WI(EN) = 0.0
      X = H(EN,NA)
      R = DSQRT(X*X+ZZ*ZZ)
      P = X / R
      Q = ZZ / R
C     ********** ROW MODIFICATION **********
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     ********** COLUMN MODIFICATION **********
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     ********** ACCUMULATE TRANSFORMATIONS **********
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     ********** COMPLEX PAIR **********
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
330   EN = ENM2
      GO TO 60
C     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM **********
  340 NORM = 0.0
      K = 1
C
      DO 360 I = 1, N
C
         DO 350 J = K, N
  350    NORM = NORM + DABS(H(I,J))
C
         K = I
  360 CONTINUE
C
      IF (NORM .EQ. 0.0) GO TO 1001
C     ********** FOR EN=N STEP -1 UNTIL 1 DO -- **********
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     ********** REAL VECTOR **********
  600    M = EN
         H(EN,EN) = 1.0
         IF (NA .EQ. 0) GO TO 800
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = H(I,EN)
            IF (M .GT. NA) GO TO 620
C
            DO 610 J = M, NA
  610       R = R + H(I,J) * H(J,EN)
C
  620       IF (WI(I) .GE. 0.0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0) GO TO 640
            T = W
            IF (W .EQ. 0.0) T = MACHEP * NORM
            H(I,EN) = -R / T
            GO TO 700
C     ********** SOLVE REAL EQUATIONS **********
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (DABS(X) .LE. DABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 700
  650       H(I+1,EN) = (-S - Y * T) / ZZ
  700    CONTINUE
C     ********** END REAL VECTOR **********
         GO TO 800
C     ********** COMPLEX VECTOR **********
  710    M = NA
C     ********** LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR **********
         IF (DABS(H(EN,NA)) .LE. DABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    TEMP=(H(NA,NA)-P)**2+Q**2
         TA=-H(NA,EN)*Q/TEMP
         TB=(P-H(NA,NA))*H(NA,EN)/TEMP
         H(NA,EN)=TB
         H(NA,NA)=TA
  730    H(EN,NA) = 0.0
         H(EN,EN) = 1.0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C
         DO 790 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0
            SA = H(I,EN)
C
            DO 760 J = M, NA
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 790
  770       M = I
            IF (WI(I) .NE. 0.0) GO TO 780
            TEMP=W**2+Q**2
            H(I,NA)=(-RA*W-SA*Q)/TEMP
            H(I,EN)=(-SA*W+RA*Q)/TEMP
            GO TO 790
C     ********** SOLVE COMPLEX EQUATIONS **********
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0 * Q
            IF (VR .EQ. 0.0 .AND. VI .EQ. 0.0) VR = MACHEP * NORM
     X       * (DABS(W) + DABS(Q) + DABS(X) + DABS(Y) + DABS(ZZ))
            TA=X*R-ZZ*RA+Q*SA
            TB=X*S-ZZ*SA-Q*RA
            TEMP=VR**2+VI**2
            H(I,NA)=(TA*VR+TB*VI)/TEMP
            H(I,EN)=(TB*VR-TA*VI)/TEMP
            IF (DABS(X) .LE. DABS(ZZ) + DABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       TA=-R-Y*H(I,NA)
            TB=-S-Y*H(I,EN)
            H(I+1,NA)=(TA*ZZ+TB*Q)/TEMP
            H(I+1,EN)=(TB*ZZ-TA*Q)/TEMP
  790    CONTINUE
C     ********** END COMPLEX VECTOR **********
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS **********
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- **********
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZ = 0.0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = EN
 1001 RETURN
C     ********** LAST CARD OF HQR2 **********
      END
C***********************************************************************
C
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      DOUBLE PRECISION SCALE(N),Z(NM,M)
      DOUBLE PRECISION S
C**PURPOSE.     TRANSLATION OF THE ALGOL PROCEDURE BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     BALBAK   FORMS THE EIGENVECTORS OF A DOUBLE PRECISION GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  BALANC.
COMMENT*YYYY;
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC,
C
C      SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  BALANC,
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     -----------------------------------------------------------------
C
COMMENT*ZZZZ;
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0/SCALE(I). **********
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
C
  110 CONTINUE
C     ********- FOR I=LOW-1 STEP -1 UNTIL 1,
C               IGH+1 STEP 1 UNTIL N DO -- **********
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE   NRMLZ(NM,N,WI,Z,MODE)
      INTEGER      I,J,ISAVE,N,NM,MODE
      DOUBLE PRECISION         S,S1,RECIPR,RECIPI
      DOUBLE PRECISION         Z(NM,N), WI(N)
       DOUBLE PRECISION          DABS
COMMENT*XXXX
C IF MODE=0
C     NORMALIZES A REAL EIGENVECTOR SO THAT ITS COMPONENT OF MAXIMUM
C     MAGNITUDE IS 1.0 AND A COMPLEX EIGENVECTOR SO THAT ITS COMPONENT
C     OF MAXIMUM MAGNITUDE IS 1.0 + 0.0*I
C IF MODE=1
C    NORMALIZES SO THAT MAXIMUM COMPONENT OF EACH NONZERO PART OF
C    THE EIGENVECTORS (REAL AND COMPLEX) IS 1.0
C MODE>1
C    MODE-1 TH COMPONENT IS 1.0+0.0*I
COMMENT*YYYY;
COMMENT*ZZZZ;
      IF(MODE.EQ.1)GO TO 1000
      ISKIP=0
      DO 500 J=1,N
      IF(ISKIP.EQ.1)GO TO 450
         S = 0.0
         IF (WI(J) .NE. 0) GO TO 250
C        REAL EIGENVECTOR
        IF(MODE.EQ.0)GO TO 90
        ISAVE=MODE-1
        S=DABS(Z(ISAVE,J))
        GO TO 110
 90      DO 100 I=1,N
            S1 = DABS(Z(I,J))
            IF(S1 .LE. S) GO TO 100
            S = S1
            ISAVE = I
 100     CONTINUE
 110     IF (S .NE. 0) S = 1./Z(ISAVE,J)
         DO 200 I=1,N
 200     Z(I,J) = Z(I,J)*S
         GO TO 500
C        COMPLEX EIGENVECTOR STORED IN COLUMN J AND J+1
 250    IF(MODE.EQ.0)GO TO 260
        ISAVE=MODE-1
        S=Z(ISAVE,J)**2+Z(ISAVE,J+1)**2
        GO TO 310
  260     DO 300 I=1,N
            S1 = Z(I,J)**2 + Z(I,J+1)**2
            IF (S1 .LE. S) GO TO 300
            S = S1
            ISAVE = I
 300     CONTINUE
 310     IF (S .EQ. 0) GO TO 500
         RECIPR = Z(ISAVE,J)/S
         RECIPI =-Z(ISAVE,J+1)/S
         DO 400 I=1,N
            S1       = Z(I,J)*RECIPR - Z(I,J+1)*RECIPI
            Z(I,J+1) = Z(I,J)*RECIPI + Z(I,J+1)*RECIPR
            Z(I,J)   = S1
 400     CONTINUE
         Z(ISAVE,J) = 1.0
         Z(ISAVE,J+1) = 0.0
         GO TO 460
 450     ISKIP=0
         GO TO 500
 460     ISKIP=1
 500  CONTINUE
      RETURN
 1000 CONTINUE
      DO 1500 J=1,N
         S = 0.0
         DO 1100 I=1,N
            S1 = DABS(Z(I,J))
            IF(S1 .LE. S) GO TO 1100
            S = S1
            ISAVE = I
 1100    CONTINUE
         IF (S .NE. 0) S = 1./Z(ISAVE,J)
         DO 1200 I=1,N
 1200    Z(I,J) = Z(I,J)*S
 1500  CONTINUE
      RETURN
      END
      SUBROUTINE GAUSS(AA,NEQ,NDIM,MDIM,ERROR)
C --- This subroutine inverts a square matrix (I-AA) using Gaussian 
C --- Elimination with partial pivoting on an augmented matrix consisting 
C --- of the matrix to be inverted followed by the identity matrix.
      DOUBLE PRECISION AA
      DIMENSION AA(NDIM,MDIM)
C     
C --- AA	Coefficient matrix augmented with the identity matrix
C --- NEQ       Number of equations
C --- NC	Total number of columns in the augmented matrix (2*NEQ)
C --- NDIM	first dimension of matrix AA in the calling program
C --- MDIM      second dimension of matrix AA in the calling program
C --- THE INVERSE IS RETURNED IN THE AUGMENTATION COLUMNS OF AA.
C
C --- Form the matrix (I-AA) and augment it with the identity matrix
      DO 7 I=1,NEQ
      DO 7 J=1,NEQ
      JJ=J+NEQ
      AA(I,JJ)=0.0D0
      IF(I.EQ.J) AA(I,JJ)=1.0D0
      AA(I,J)=AA(I,JJ)-AA(I,J)
    7 CONTINUE
C --- Begin the reduction
      NC=2*NEQ
      NEQ1=NEQ-1
      DO 1 I=1,NEQ1
C --- Find the row number of the pivot row, then interchange rows to
C     put the pivot element on the diagonal.
      IPVT=I
      IP1=I+1
      DO 2 J=IP1,NEQ
      IF(ABS(AA(IPVT,I)).LT.ABS(AA(J,I))) IPVT=J
    2 CONTINUE
C --- If the pivot element is too small, print a message and return.
      IF (ABS(AA(IPVT,I)).LT.ERROR) GO TO 99
C --- Unless the pivot element is already on the diagonal, interchange.
      IF (IPVT.EQ.I) GO TO 3
      DO 4 JCOL=I,NC
      SAVE=AA(I,JCOL)
      AA(I,JCOL)=AA(IPVT,JCOL)
      AA(IPVT,JCOL)=SAVE
    4 CONTINUE
C --- Reduce all elementsbelow the diagonal in the I-th row. If a zero 
C     is already present, do not reduce that row.
    3 CONTINUE 
      DO 5 JROW=IP1,NEQ
      IF(AA(JROW,I).EQ.0) GO TO 5
      RATIO=AA(JROW,I)/AA(I,I)
      DO 6 KCOL=IP1,NC
      AA(JROW,KCOL)=AA(JROW,KCOL)-RATIO*AA(I,KCOL)
    6 CONTINUE
    5 CONTINUE
    1 CONTINUE
C --- Check AA(NEQ,NEQ) for size.
      IF (ABS(AA(NEQ,NEQ)).LT.ERROR) GO TO 99
C --- Back sustitute.
      NC1=NEQ+1
      DO 10 KCOL=NC1,NC
      AA(NEQ,KCOL)=AA(NEQ,KCOL)/AA(NEQ,NEQ)
      DO 11 J=2,NEQ
      NVAR=NC1-J
      L=NVAR+1
      VALUE=AA(NVAR,KCOL)
      DO 12 K=L,NEQ
      VALUE=VALUE-AA(NVAR,K)*AA(K,KCOL)
   12 CONTINUE
      AA(NVAR,KCOL)=VALUE/AA(NVAR,NVAR)
   11 CONTINUE
   10 CONTINUE
      RETURN
   99 WRITE(*,*) 'NEAR ZERO PIVOT ENCOUNTERED'
      WRITE(4,*) 'NEAR ZERO PIVOT ENCOUNTERED'
      WRITE(4,5000) IPVT,I
 5000 FORMAT('IPVT = ',I2,' I = ',I2)
      RETURN
      END
      SUBROUTINE SMMULT(A,B,C,D,M,N)
C --- This subroutine multiplies two square matrices
C
C --- A = first matrix
C --- B = second matrix
C --- C = AB
C --- D = dummy array
C --- M = maximal dimension of A,B,and C in calling program 
C --- N = actual dimension of A, B, and C in calling program
C
      DOUBLE PRECISION A(M,M),B(M,M),C(M,M),D(M,M)
      DO 1 I=1,N
      DO 1 J=1,N
      D(I,J)=A(I,1)*B(1,J)
      DO 1 K=2,N
      D(I,J)=D(I,J)+A(I,K)*B(K,J)
    1 CONTINUE
      CALL STORE(M,M,M,M,N,N,C,D)
      RETURN
      END
      SUBROUTINE INITV(X,N,M)
      DOUBLE PRECISION X(N)
      DO 1 I=1,M
      X(I)=0.0D0
    1 CONTINUE
      RETURN
      END
      SUBROUTINE INITM(A,N,M,K,L)
      DOUBLE PRECISION A(N,M)
      DO 1 I=1,K
      DO 1 J=1,L
      A(I,J)=0.0D0
    1 CONTINUE
      RETURN
      END
      SUBROUTINE COLSUM(A,N,M,NI,MI,X,I)
      DOUBLE PRECISION A(N,M), X(I)
      DO 1 J=1,NI
      X(J)=A(1,J)
      DO 1 K=2,MI
      X(J)=X(J)+A(K,J)
    1 CONTINUE
      RETURN
      END
      SUBROUTINE PROB1(M,M1,N1,NN,N2,PP)
C --- THIS SUBROUTINE COMPUTES PROBABILITIES FROM N1-1 CLASSES TO
C --- N1 CLASSES
      DOUBLE PRECISION PP
      DIMENSION NN(M,M1),PP(M,M),N2(M)
      NI=N1-1
      DO 2 I=1,NI
      N2(I)=NN(I,N1)
      DO 3 J=1,NI
      N2(I)=N2(I)+NN(I,J)
    3 CONTINUE
      DO 5 J=1,NI
      IF(N2(I).EQ.0.OR.NN(I,J).EQ.0) GO TO 4
      PP(J,I)=FLOAT(NN(I,J))/FLOAT(N2(I))
      GO TO 5
    4 PP(J,I)=0.0
    5 CONTINUE
    2 CONTINUE
      RETURN
      END

      SUBROUTINE STORE(MMAXA,NMAXA,MMAXB,NMAXB,M,N,A,B)
      DOUBLE PRECISION A(MMAXA,NMAXA),B(MMAXB,NMAXB)
      DO 1 I=1,M
      DO 1 J=1,N
      A(I,J)=B(I,J)
    1 CONTINUE
      RETURN
      END

      SUBROUTINE DOMINANT(M,XREAL,XIMAG,XNORM,N,IND)
C --- This subroutine finds the dominant eigenvalue, where 
C --- XREAL(I) = real part of eigenvalue I
C --- XIMAG(I) = imanginary part of eigenvalue I
      DOUBLE PRECISION XREAL(M),XIMAG(M),XNORM(M)
      DO 1 I=1,N
      XNORM(I)=(XREAL(I)**2+XIMAG(I)**2)**.5
    1 CONTINUE
C --- Find the dominant eigenvalue of C
      IND=1
      DO 2 I=2,N
      IF(XNORM(I).GT.XNORM(IND)) IND=I
    2 CONTINUE
      RETURN
      END
      SUBROUTINE NORMLZ(M,N,V,I,E)
C
C --- M = maximum dimension of the eigenvector matrix in the calling program
C --- N = actual dimension of the eigenvector
C --- V = eigenvector matrix
C --- I = column of v containing the dominant eigenvector
C --- E = eigenvector
C
      DOUBLE PRECISION E,X,V
      DIMENSION E(M),V(M,M)
      DO 1 J=1,N
      E(J)=V(J,I)
    1 CONTINUE
      X=E(1)
      DO 2 J=2,N
      X=X+E(J)
    2 CONTINUE
      DO 3 J=1,N
      E(J)=E(J)/X
    3 CONTINUE
      RETURN
      END

