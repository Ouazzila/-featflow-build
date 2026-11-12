************************************************************************
* bnstok.f  -  Stage 2: (generalized) Stokes problem with Bernstein
*              Qk/P(k-1)disc, compared with the original pp2d pair
*              Q1~/Q0 (FEAT E030/E031), same solver for both.
*
*   [ a M + nu A ] u + B p = f ,   B^T u = 0 ,   u = g on the boundary
*
*   B(i,q) = -int psi_q div phi_i  (see BNASM).  a = 0: steady Stokes;
*   a = 1/(theta K), nu = alpha/Re: matrix of a pp2d theta-substep.
*
* Solver (Glowinski, Handbook IX, Ch. IV par. 21): conjugate gradient
* (Uzawa-CG) on the Schur complement S = B^T (aM+nuA)^(-1) B with the
* Cahouet-Chabard preconditioner
*      S^(-1) ~ nu Mp^(-1) + a (B^T Ml^(-1) B)^(-1)
* (Mp pressure mass matrix, element blocks; Ml lumped velocity mass =
* row sums of M). The inner velocity solves use SSOR-PCG, the
* pressure Poisson part Jacobi-CG; all on the same CSR data for both
* element families (COMMON /BNFAM/, see bern.f).
* Pressure is fixed by mean zero (all-Dirichlet problem): dual
* vectors (residuals) are projected to sum 0, primal vectors
* (pressures) to integral 0.
************************************************************************
C
************************************************************************
*   Exact solution (divergence-free, stream function cos3x sin4y/...)
*     u1 = 4 cos3x cos4y ,  u2 = 3 sin3x sin4y ,  -Lap u = 25 u
*     p  = cos2x e^y
*     f  = (a + 25 nu) u + grad p     with a, nu from /BNSTKP/
************************************************************************
      DOUBLE PRECISION FUNCTION SKU1(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKU1=4D0*COS(3D0*X)*COS(4D0*Y)
      END
      DOUBLE PRECISION FUNCTION SKU1X(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKU1X=-12D0*SIN(3D0*X)*COS(4D0*Y)
      END
      DOUBLE PRECISION FUNCTION SKU1Y(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKU1Y=-16D0*COS(3D0*X)*SIN(4D0*Y)
      END
      DOUBLE PRECISION FUNCTION SKU2(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKU2=3D0*SIN(3D0*X)*SIN(4D0*Y)
      END
      DOUBLE PRECISION FUNCTION SKU2X(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKU2X=9D0*COS(3D0*X)*SIN(4D0*Y)
      END
      DOUBLE PRECISION FUNCTION SKU2Y(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKU2Y=12D0*SIN(3D0*X)*COS(4D0*Y)
      END
      DOUBLE PRECISION FUNCTION SKP(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SKP=COS(2D0*X)*EXP(Y)
      END
      DOUBLE PRECISION FUNCTION SKF1(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BNSTKP/ DSA,DSNU
      SKF1=(DSA+25D0*DSNU)*4D0*COS(3D0*X)*COS(4D0*Y)
     *    -2D0*SIN(2D0*X)*EXP(Y)
      END
      DOUBLE PRECISION FUNCTION SKF2(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BNSTKP/ DSA,DSNU
      SKF2=(DSA+25D0*DSNU)*3D0*SIN(3D0*X)*SIN(4D0*Y)
     *    +COS(2D0*X)*EXP(Y)
      END
C
************************************************************************
*   BNBDV: Dirichlet coefficients of the Bernstein Qk velocity
*     on each boundary edge (KADJ=0): vertex coefficients = U(vertex)
*     (B_0, B_K interpolate the end points), edge coefficients by the
*     1D L2 projection of U - c_0 B_0 - c_K B_K onto B_1..B_(K-1).
*     Edge position IS counted from KVERT(IME) as in BNDFMP.
************************************************************************
      SUBROUTINE BNBDV(K,DCORVG,KVERT,KADJ,KDFG,NEL,U,DU)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNG=NNBK+2,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KADJ(NNVE,*)
      DIMENSION KDFG((K+1)**2,*),DU(*)
      DIMENSION XG(NNG),WG(NNG),DB(0:NNBK),DDB(0:NNBK)
      DIMENSION DG(NNBK,NNBK),DR(NNBK),LPOS(0:NNBK)
      EXTERNAL U
C
      K1=K+1
      NG=K+2
      CALL BNGAUS(NG,XG,WG)
      DO 100 IEL=1,NEL
      DO 100 IME=1,4
       IF (KADJ(IME,IEL).NE.0) GOTO 100
       IV1=KVERT(IME,IEL)
       IV2=KVERT(MOD(IME,4)+1,IEL)
       X1=DCORVG(1,IV1)
       Y1=DCORVG(2,IV1)
       X2=DCORVG(1,IV2)
       Y2=DCORVG(2,IV2)
       DO 10 IS=0,K
        IF (IME.EQ.1) LPOS(IS)=IS+1
        IF (IME.EQ.2) LPOS(IS)=K+IS*K1+1
        IF (IME.EQ.3) LPOS(IS)=(K-IS)+K*K1+1
        IF (IME.EQ.4) LPOS(IS)=(K-IS)*K1+1
10     CONTINUE
       C0=U(X1,Y1)
       CK=U(X2,Y2)
       DU(KDFG(LPOS(0),IEL))=C0
       DU(KDFG(LPOS(K),IEL))=CK
C *** 1D projection of the remainder (Gram matrix on [-1,1])
       DO 20 IR=1,K-1
        DR(IR)=0D0
        DO 20 IC=1,K-1
20     DG(IR,IC)=0D0
       DO 30 IQ=1,NG
        S=XG(IQ)
        T=0.5D0*(S+1D0)
        CALL BNB1D(K,S,DB,DDB)
        RES=U(X1+T*(X2-X1),Y1+T*(Y2-Y1))-C0*DB(0)-CK*DB(K)
        DO 25 IR=1,K-1
         DR(IR)=DR(IR)+WG(IQ)*RES*DB(IR)
         DO 25 IC=1,K-1
25      DG(IR,IC)=DG(IR,IC)+WG(IQ)*DB(IR)*DB(IC)
30     CONTINUE
       NM=K-1
       CALL BNGSOL(DG,NNBK,NM,DR)
       DO 40 IS=1,K-1
40     DU(KDFG(LPOS(IS),IEL))=DR(IS)
100   CONTINUE
C
      END
C
************************************************************************
*   BQBDV: Dirichlet values of Q1~ on boundary edges:
*          E031 (IFAMQ=31): value in the edge midpoint,
*          E030 (IFAMQ=30): integral mean over the edge (4-pt Gauss)
************************************************************************
      SUBROUTINE BQBDV(IFAMQ,DCORVG,KVERT,KMID,KADJ,NEL,NVT,U,DU)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*)
      DIMENSION DU(*),XG(4),WG(4)
      EXTERNAL U
C
      CALL BNGAUS(4,XG,WG)
      DO 100 IEL=1,NEL
      DO 100 IME=1,4
       IF (KADJ(IME,IEL).NE.0) GOTO 100
       IV1=KVERT(IME,IEL)
       IV2=KVERT(MOD(IME,4)+1,IEL)
       X1=DCORVG(1,IV1)
       Y1=DCORVG(2,IV1)
       X2=DCORVG(1,IV2)
       Y2=DCORVG(2,IV2)
       IMT=KMID(IME,IEL)-NVT
       IF (IFAMQ.EQ.31) THEN
        DU(IMT)=U(0.5D0*(X1+X2),0.5D0*(Y1+Y2))
       ELSE
        S=0D0
        DO 10 IQ=1,4
         T=0.5D0*(XG(IQ)+1D0)
         S=S+0.5D0*WG(IQ)*U(X1+T*(X2-X1),Y1+T*(Y2-Y1))
10      CONTINUE
        DU(IMT)=S
       ENDIF
100   CONTINUE
C
      END
C
************************************************************************
*   BNGSOL: Gauss elimination (no pivoting, SPD matrices) for the
*           N x N system DG x = DR (leading dimension LDG), x -> DR
************************************************************************
      SUBROUTINE BNGSOL(DG,LDG,N,DR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DG(LDG,*),DR(*)
C
      DO 30 IC=1,N
      DO 30 IR=IC+1,N
       F=DG(IR,IC)/DG(IC,IC)
       DO 25 J=IC,N
25     DG(IR,J)=DG(IR,J)-F*DG(IC,J)
       DR(IR)=DR(IR)-F*DR(IC)
30    CONTINUE
      DO 35 IR=N,1,-1
       S=DR(IR)
       DO 34 J=IR+1,N
34     S=S-DG(IR,J)*DR(J)
       DR(IR)=S/DG(IR,IR)
35    CONTINUE
C
      END
C
************************************************************************
*   BNMPI: inverse element pressure mass blocks DMPI(NPL,NPL,IEL)
*          and DWP(q) = int psi_q  (Mp * 1)
************************************************************************
      SUBROUTINE BNMPI(NEL,DCORVG,KVERT,DMPI,DWP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),DMPI(NPL,NPL,*),DWP(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      DIMENSION DG(NNPL,NNPL),DR(NNPL)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
C
      CALL BNGAUS(NG,XG,WG)
      DO 50 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNGSET(IEL)
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
       DO 20 IQ=1,NPL
        S=0D0
        DO 10 IG=1,NG*NG
10      S=S+DWT(IG)*DPSI(IQ,IG)
        DWP((IEL-1)*NPL+IQ)=S
20     CONTINUE
C *** inverse column by column
       DO 40 JC=1,NPL
        DO 30 I=1,NPL
         DR(I)=0D0
         DO 30 J=1,NPL
          S=0D0
          DO 25 IG=1,NG*NG
25        S=S+DWT(IG)*DPSI(I,IG)*DPSI(J,IG)
          DG(I,J)=S
30      CONTINUE
        DR(JC)=1D0
        CALL BNGSOL(DG,NNPL,NPL,DR)
        DO 35 I=1,NPL
35      DMPI(I,JC,IEL)=DR(I)
40     CONTINUE
50    CONTINUE
      CALL BNGSET(0)
C
      END
C
************************************************************************
*   BNPERR: L2 error of the pressure, both means removed
************************************************************************
      SUBROUTINE BNPERR(NEL,DCORVG,KVERT,KDFP,DP,P,ERR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFP(NPL,*),DP(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      EXTERNAL P
C
      CALL BNGAUS(NG,XG,WG)
      DO 100 IPASS=1,2
       IF (IPASS.EQ.1) THEN
        AREA=0D0
        PM=0D0
        PHM=0D0
       ELSE
        PM=PM/AREA
        PHM=PHM/AREA
        ERR=0D0
       ENDIF
       DO 50 IEL=1,NEL
        DO 5 IV=1,4
         DXV(IV)=DCORVG(1,KVERT(IV,IEL))
         DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5       CONTINUE
        CALL BNGSET(IEL)
        CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
        DO 40 IG=1,NG*NG
         PH=0D0
         DO 30 IQ=1,NPL
30       PH=PH+DP(KDFP(IQ,IEL))*DPSI(IQ,IG)
         PE=P(DXQ(IG),DYQ(IG))
         IF (IPASS.EQ.1) THEN
          AREA=AREA+DWT(IG)
          PM =PM +DWT(IG)*PE
          PHM=PHM+DWT(IG)*PH
         ELSE
          ERR=ERR+DWT(IG)*((PH-PHM)-(PE-PM))**2
         ENDIF
40      CONTINUE
50     CONTINUE
100   CONTINUE
      CALL BNGSET(0)
      ERR=SQRT(ERR)
C
      END
C
************************************************************************
*   BNSSOR: z = (D+U)^(-1) D (D+L)^(-1) r on the free dofs (KBD=0),
*           z = 0 on Dirichlet dofs (symmetric CSR, diagonal first)
************************************************************************
      SUBROUTINE BNSSOR(DA,KLD,KCOL,NEQ,KBD,DR,DZ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),KLD(*),KCOL(*),KBD(*),DR(*),DZ(*)
C
      DO 20 I=1,NEQ
       IF (KBD(I).EQ.1) THEN
        DZ(I)=0D0
        GOTO 20
       ENDIF
       S=DR(I)
       DO 10 IA=KLD(I)+1,KLD(I+1)-1
        J=KCOL(IA)
        IF ((J.LT.I).AND.(KBD(J).EQ.0)) S=S-DA(IA)*DZ(J)
10     CONTINUE
       DZ(I)=S/DA(KLD(I))
20    CONTINUE
      DO 40 I=NEQ,1,-1
       IF (KBD(I).EQ.1) GOTO 40
       S=DA(KLD(I))*DZ(I)
       DO 30 IA=KLD(I)+1,KLD(I+1)-1
        J=KCOL(IA)
        IF ((J.GT.I).AND.(KBD(J).EQ.0)) S=S-DA(IA)*DZ(J)
30     CONTINUE
       DZ(I)=S/DA(KLD(I))
40    CONTINUE
C
      END
C
************************************************************************
*   BNVSOL: SSOR-PCG for the velocity matrix DA, Dirichlet dofs
*           (KBD=1) keep their value in DX. |r| <= EPS*|r0|.
*           ITE = number of iterations (-ITMAX: no convergence)
************************************************************************
      SUBROUTINE BNVSOL(DA,KLD,KCOL,NEQ,DX,DB,KBD,DR,DP,DQ,DZ,
     *                  EPS,ITMAX,ITE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),KLD(*),KCOL(*),DX(*),DB(*),KBD(*)
      DIMENSION DR(*),DP(*),DQ(*),DZ(*)
C
      CALL BNAX(DA,KLD,KCOL,NEQ,DX,DR)
      RN0=0D0
      DO 10 I=1,NEQ
       DR(I)=DB(I)-DR(I)
       IF (KBD(I).EQ.1) DR(I)=0D0
       RN0=RN0+DR(I)**2
10    CONTINUE
      RN0=SQRT(RN0)
      ITE=0
      IF (RN0.LE.1D-300) RETURN
      CALL BNSSOR(DA,KLD,KCOL,NEQ,KBD,DR,DZ)
      RZ=0D0
      DO 15 I=1,NEQ
       DP(I)=DZ(I)
       RZ=RZ+DR(I)*DZ(I)
15    CONTINUE
C
      DO 50 ITE=1,ITMAX
       CALL BNAX(DA,KLD,KCOL,NEQ,DP,DQ)
       PQ=0D0
       DO 20 I=1,NEQ
        IF (KBD(I).EQ.1) DQ(I)=0D0
        PQ=PQ+DP(I)*DQ(I)
20     CONTINUE
       ALPHA=RZ/PQ
       RN=0D0
       DO 30 I=1,NEQ
        DX(I)=DX(I)+ALPHA*DP(I)
        DR(I)=DR(I)-ALPHA*DQ(I)
        RN=RN+DR(I)**2
30     CONTINUE
       IF (SQRT(RN).LE.EPS*RN0) RETURN
       CALL BNSSOR(DA,KLD,KCOL,NEQ,KBD,DR,DZ)
       RZN=0D0
       DO 35 I=1,NEQ
35     RZN=RZN+DR(I)*DZ(I)
       GAM=RZN/RZ
       RZ=RZN
       DO 40 I=1,NEQ
40     DP(I)=DZ(I)+GAM*DP(I)
50    CONTINUE
      ITE=-ITMAX
C
      END
C
************************************************************************
*   BNBX : DY1 = B1 p, DY2 = B2 p                      (velocity space)
*   BNBTX: DR  = B1^T u1 + B2^T u2  (= -int psi div u_h) (pressure dual)
************************************************************************
      SUBROUTINE BNBX(DB1,DB2,KLDB,KCOLB,NU,DP,DY1,DY2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DB1(*),DB2(*),KLDB(*),KCOLB(*),DP(*),DY1(*),DY2(*)
      DO 20 I=1,NU
       S1=0D0
       S2=0D0
       DO 10 IA=KLDB(I),KLDB(I+1)-1
        S1=S1+DB1(IA)*DP(KCOLB(IA))
        S2=S2+DB2(IA)*DP(KCOLB(IA))
10     CONTINUE
       DY1(I)=S1
       DY2(I)=S2
20    CONTINUE
      END
C
      SUBROUTINE BNBTX(DB1,DB2,KLDB,KCOLB,NU,NP,DU1,DU2,DR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DB1(*),DB2(*),KLDB(*),KCOLB(*),DU1(*),DU2(*),DR(*)
      DO 5 IQ=1,NP
5     DR(IQ)=0D0
      DO 20 I=1,NU
       DO 10 IA=KLDB(I),KLDB(I+1)-1
        IQ=KCOLB(IA)
        DR(IQ)=DR(IQ)+DB1(IA)*DU1(I)+DB2(IA)*DU2(I)
10     CONTINUE
20    CONTINUE
      END
C
************************************************************************
*   BNPJD: dual vector  r <- r - (sum r / |Omega|) * DWP   (sum r = 0)
*   BNPJP: primal vector p <- p - (DWP.p / |Omega|)        (int p = 0)
************************************************************************
      SUBROUTINE BNPJD(DR,DWP,NP,AREA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DR(*),DWP(*)
      S=0D0
      DO 10 IQ=1,NP
10    S=S+DR(IQ)
      S=S/AREA
      DO 20 IQ=1,NP
20    DR(IQ)=DR(IQ)-S*DWP(IQ)
      END
C
      SUBROUTINE BNPJP(DP,DWP,NP,AREA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DP(*),DWP(*)
      S=0D0
      DO 10 IQ=1,NP
10    S=S+DWP(IQ)*DP(IQ)
      S=S/AREA
      DO 20 IQ=1,NP
20    DP(IQ)=DP(IQ)-S
      END
C
************************************************************************
*   BNCX: DY = C x,  C = B^T Ml^(-1) B with Dirichlet rows of B
*         removed (the discrete projection matrix of pp2d, projma.f)
*   BNCDG: diagonal of C
************************************************************************
      SUBROUTINE BNCX(DB1,DB2,KLDB,KCOLB,NU,NP,DML,KBD,DX,DY,DT1,DT2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DB1(*),DB2(*),KLDB(*),KCOLB(*),DML(*),KBD(*)
      DIMENSION DX(*),DY(*),DT1(*),DT2(*)
      CALL BNBX(DB1,DB2,KLDB,KCOLB,NU,DX,DT1,DT2)
      DO 10 I=1,NU
       IF (KBD(I).EQ.1) THEN
        DT1(I)=0D0
        DT2(I)=0D0
       ELSE
        DT1(I)=DT1(I)/DML(I)
        DT2(I)=DT2(I)/DML(I)
       ENDIF
10    CONTINUE
      CALL BNBTX(DB1,DB2,KLDB,KCOLB,NU,NP,DT1,DT2,DY)
      END
C
      SUBROUTINE BNCDG(DB1,DB2,KLDB,KCOLB,NU,NP,DML,KBD,DCD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DB1(*),DB2(*),KLDB(*),KCOLB(*),DML(*),KBD(*),DCD(*)
      DO 5 IQ=1,NP
5     DCD(IQ)=0D0
      DO 20 I=1,NU
       IF (KBD(I).EQ.1) GOTO 20
       DO 10 IA=KLDB(I),KLDB(I+1)-1
        IQ=KCOLB(IA)
        DCD(IQ)=DCD(IQ)+(DB1(IA)**2+DB2(IA)**2)/DML(I)
10     CONTINUE
20    CONTINUE
      END
C
************************************************************************
*   BNPREC: Cahouet-Chabard preconditioner, z = nu Mp^(-1) r
*           + a C^(-1) r (Jacobi-CG on C, |res| <= EPSC*|r|);
*           z is returned with integral 0. ITC: CG iterations.
*           DPA(NP,4) auxiliary, DT(NU,2) auxiliary.
************************************************************************
      SUBROUTINE BNPREC(DB1,DB2,KLDB,KCOLB,NU,NP,NEL,DML,KBD,DMPI,
     *                  DWP,DCD,AREA,DSA,DSNU,DR,DZ,DPA,DT,EPSC,ITC)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DB1(*),DB2(*),KLDB(*),KCOLB(*),DML(*),KBD(*)
      DIMENSION DMPI(NPL,NPL,*),DWP(*),DCD(*),DR(*),DZ(*)
      DIMENSION DPA(NP,4),DT(NU,2)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
C
C *** z = nu Mp^(-1) r (element blocks)
      DO 20 IEL=1,NEL
       IOFF=(IEL-1)*NPL
       DO 20 I=1,NPL
        S=0D0
        DO 10 J=1,NPL
10      S=S+DMPI(I,J,IEL)*DR(IOFF+J)
        DZ(IOFF+I)=DSNU*S
20    CONTINUE
      ITC=0
      IF (DSA.EQ.0D0) GOTO 99
C
C *** x = C^(-1) r by Jacobi-CG, x in DPA(.,1), r in DPA(.,2),
C     direction DPA(.,3), C*direction / preconditioned r in DPA(.,4).
C     C is singular (constant pressures): the residual is kept dual
C     mean-free and the preconditioned residual primal mean-free, so
C     x never picks up a (large) constant from rounding.
      RN0=0D0
      DO 30 IQ=1,NP
       DPA(IQ,1)=0D0
       DPA(IQ,2)=DR(IQ)
       RN0=RN0+DR(IQ)**2
30    CONTINUE
      RN0=SQRT(RN0)
      IF (RN0.LE.1D-300) GOTO 90
      CALL BNPJD(DPA(1,2),DWP,NP,AREA)
      RZ=0D0
      DO 32 IQ=1,NP
32    DPA(IQ,4)=DPA(IQ,2)/DCD(IQ)
      CALL BNPJP(DPA(1,4),DWP,NP,AREA)
      DO 34 IQ=1,NP
       DPA(IQ,3)=DPA(IQ,4)
       RZ=RZ+DPA(IQ,2)*DPA(IQ,4)
34    CONTINUE
      DO 60 ITC=1,10000
       CALL BNCX(DB1,DB2,KLDB,KCOLB,NU,NP,DML,KBD,DPA(1,3),DPA(1,4),
     *           DT(1,1),DT(1,2))
       PQ=0D0
       DO 40 IQ=1,NP
40     PQ=PQ+DPA(IQ,3)*DPA(IQ,4)
       ALPHA=RZ/PQ
       DO 45 IQ=1,NP
        DPA(IQ,1)=DPA(IQ,1)+ALPHA*DPA(IQ,3)
        DPA(IQ,2)=DPA(IQ,2)-ALPHA*DPA(IQ,4)
45     CONTINUE
       CALL BNPJD(DPA(1,2),DWP,NP,AREA)
       RN=0D0
       DO 46 IQ=1,NP
        RN=RN+DPA(IQ,2)**2
        DPA(IQ,4)=DPA(IQ,2)/DCD(IQ)
46     CONTINUE
       IF (SQRT(RN).LE.EPSC*RN0) GOTO 90
       CALL BNPJP(DPA(1,4),DWP,NP,AREA)
       RZN=0D0
       DO 48 IQ=1,NP
48     RZN=RZN+DPA(IQ,2)*DPA(IQ,4)
       GAM=RZN/RZ
       RZ=RZN
       DO 50 IQ=1,NP
50     DPA(IQ,3)=DPA(IQ,4)+GAM*DPA(IQ,3)
60    CONTINUE
      ITC=-10000
90    DO 95 IQ=1,NP
95    DZ(IQ)=DZ(IQ)+DSA*DPA(IQ,1)
99    CALL BNPJP(DZ,DWP,NP,AREA)
C
      END
C
************************************************************************
*   BNUZAW: Uzawa-CG for  Ah u + B p = f, B^T u = 0  (Ah = aM + nu A)
*     (flexible PCG, Polak-Ribiere direction update)
*     input : U1,U2 with the Dirichlet values (0 elsewhere), F1,F2
*     output: U1,U2,P (P with integral 0)
*     stop  : |B^T u|_2 <= EPS*|B^T u_0|_2 (dual residual, mean removed)
*     ITUZ Uzawa steps, ITV velocity PCG steps (total), ITC Poisson
*     CG steps in the preconditioner (total), RES0/RES residuals
*     DAUX(NU,8): W1 W2 Y1 Y2 + 4 PCG vectors, DPAUX(NP,8)
************************************************************************
      SUBROUTINE BNUZAW(NU,NP,NEL,DAH,KLDA,KCOLA,DB1,DB2,KLDB,KCOLB,
     *                  KBD,DML,DMPI,DWP,DCD,AREA,DSA,DSNU,F1,F2,
     *                  U1,U2,P,DAUX,DPAUX,EPS,ITMAX,ITUZ,ITV,ITC,
     *                  RES0,RES)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DAH(*),KLDA(*),KCOLA(*),DB1(*),DB2(*),KLDB(*),KCOLB(*)
      DIMENSION KBD(*),DML(*),DMPI(*),DWP(*),DCD(*),F1(*),F2(*)
      DIMENSION U1(*),U2(*),P(*),DAUX(NU,8),DPAUX(NP,8)
      PARAMETER (EPSV=1D-12,EPSC=1D-12,ITVMAX=20000)
C
      ITV=0
      ITC=0
      ITUZ=0
      DO 5 IQ=1,NP
5     P(IQ)=0D0
C
C *** u_0 = Ah^(-1) f (with boundary values),  r_0 = B^T u_0
      CALL BNVSOL(DAH,KLDA,KCOLA,NU,U1,F1,KBD,DAUX(1,5),DAUX(1,6),
     *            DAUX(1,7),DAUX(1,8),EPSV,ITVMAX,ITE)
      ITV=ITV+ABS(ITE)
      CALL BNVSOL(DAH,KLDA,KCOLA,NU,U2,F2,KBD,DAUX(1,5),DAUX(1,6),
     *            DAUX(1,7),DAUX(1,8),EPSV,ITVMAX,ITE)
      ITV=ITV+ABS(ITE)
C     DPAUX: 1 r, 2 z, 3 d, 4 S d, 5..8 preconditioner
      CALL BNBTX(DB1,DB2,KLDB,KCOLB,NU,NP,U1,U2,DPAUX(1,1))
      CALL BNPJD(DPAUX(1,1),DWP,NP,AREA)
      RES0=0D0
      DO 10 IQ=1,NP
10    RES0=RES0+DPAUX(IQ,1)**2
      RES0=SQRT(RES0)
      RES=RES0
      IF (RES0.LE.1D-300) RETURN
      CALL BNPREC(DB1,DB2,KLDB,KCOLB,NU,NP,NEL,DML,KBD,DMPI,DWP,DCD,
     *            AREA,DSA,DSNU,DPAUX(1,1),DPAUX(1,2),DPAUX(1,5),
     *            DAUX(1,3),EPSC,ITE)
      ITC=ITC+ABS(ITE)
      RZ=0D0
      DO 15 IQ=1,NP
       DPAUX(IQ,3)=DPAUX(IQ,2)
       RZ=RZ+DPAUX(IQ,1)*DPAUX(IQ,2)
15    CONTINUE
C
      DO 100 ITUZ=1,ITMAX
C ***  w = Ah^(-1) B d (homogeneous boundary values), S d = B^T w
       CALL BNBX(DB1,DB2,KLDB,KCOLB,NU,DPAUX(1,3),DAUX(1,3),DAUX(1,4))
       DO 20 I=1,NU
        DAUX(I,1)=0D0
        DAUX(I,2)=0D0
20     CONTINUE
       CALL BNVSOL(DAH,KLDA,KCOLA,NU,DAUX(1,1),DAUX(1,3),KBD,
     *             DAUX(1,5),DAUX(1,6),DAUX(1,7),DAUX(1,8),
     *             EPSV,ITVMAX,ITE)
       ITV=ITV+ABS(ITE)
       CALL BNVSOL(DAH,KLDA,KCOLA,NU,DAUX(1,2),DAUX(1,4),KBD,
     *             DAUX(1,5),DAUX(1,6),DAUX(1,7),DAUX(1,8),
     *             EPSV,ITVMAX,ITE)
       ITV=ITV+ABS(ITE)
       CALL BNBTX(DB1,DB2,KLDB,KCOLB,NU,NP,DAUX(1,1),DAUX(1,2),
     *            DPAUX(1,4))
       CALL BNPJD(DPAUX(1,4),DWP,NP,AREA)
C ***  flexible CG: alpha = (r,d)/(d,Sd)  (the preconditioner contains
C      an inner CG solve and is not exactly a fixed linear operator)
       DSD=0D0
       RD=0D0
       DO 30 IQ=1,NP
        DSD=DSD+DPAUX(IQ,3)*DPAUX(IQ,4)
        RD=RD+DPAUX(IQ,1)*DPAUX(IQ,3)
30     CONTINUE
       ALPHA=RD/DSD
C ***  p += alpha d,  u -= alpha w,  r -= alpha S d
       RES=0D0
       DO 40 IQ=1,NP
        P(IQ)=P(IQ)+ALPHA*DPAUX(IQ,3)
        DPAUX(IQ,1)=DPAUX(IQ,1)-ALPHA*DPAUX(IQ,4)
        RES=RES+DPAUX(IQ,1)**2
40     CONTINUE
       RES=SQRT(RES)
       DO 45 I=1,NU
        U1(I)=U1(I)-ALPHA*DAUX(I,1)
        U2(I)=U2(I)-ALPHA*DAUX(I,2)
45     CONTINUE
       IF (RES.LE.EPS*RES0) GOTO 200
C ***  z = Cahouet-Chabard(r), new direction
       CALL BNPREC(DB1,DB2,KLDB,KCOLB,NU,NP,NEL,DML,KBD,DMPI,DWP,DCD,
     *             AREA,DSA,DSNU,DPAUX(1,1),DPAUX(1,2),DPAUX(1,5),
     *             DAUX(1,3),EPSC,ITE)
       ITC=ITC+ABS(ITE)
C ***  Polak-Ribiere: gam = (z_new, r_new-r_old)/(z_old,r_old),
C      r_new - r_old = -alpha S d
       RZN=0D0
       ZSD=0D0
       DO 50 IQ=1,NP
        RZN=RZN+DPAUX(IQ,1)*DPAUX(IQ,2)
        ZSD=ZSD+DPAUX(IQ,2)*DPAUX(IQ,4)
50     CONTINUE
       GAM=-ALPHA*ZSD/RZ
       RZ=RZN
       DO 60 IQ=1,NP
60     DPAUX(IQ,3)=DPAUX(IQ,2)+GAM*DPAUX(IQ,3)
100   CONTINUE
      ITUZ=-ITMAX
200   CALL BNPJP(P,DWP,NP,AREA)
C
      END
C
************************************************************************
*   BNSTKT: Stokes comparison test on the levels NLMIN..NLMAX
*     IBERN = k >= 2 : Bernstein Qk / P(k-1)disc
*     IBERN = -30    : original pair, E030 (Q1~ edge means) / Q0
*     IBERN = -31    : original pair, E031 (Q1~ midpoints)  / Q0
*   Case S: steady Stokes, a = 0, nu = 1
*   Case G: theta-substep of pp2d, a = 1/(theta K), K = 0.01,
*           nu = alpha/Re, Re = 1000  (theta = 1-1/sqrt2, alpha=2-sqrt2)
************************************************************************
      SUBROUTINE BNSTKT(IBERN,MFILE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      CHARACTER SUB*6,FMT*15,CPARAM*120,CFAM*16
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNBK=8)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      COMMON /BNSTKP/ DSA,DSNU
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      DIMENSION ERR(NNLEV,2,4)
      EXTERNAL SKU1,SKU1X,SKU1Y,SKU2,SKU2X,SKU2Y,SKP,SKF1,SKF2
      SAVE
C
      SUB='BNSTKT'
      IF (IBERN.GE.2.AND.IBERN.LE.NNBK) THEN
       IFAM=0
       KORD=IBERN
       NLDF=(KORD+1)**2
       NPL =KORD*(KORD+1)/2
       NG  =KORD+2
       WRITE(CFAM,'(A,I1,A,I1,A)') 'Bern. Q',KORD,'/P',KORD-1,'disc'
      ELSE IF ((IBERN.EQ.-30).OR.(IBERN.EQ.-31).OR.
     *         (IBERN.EQ.-130).OR.(IBERN.EQ.-131)) THEN
       IFAM=-IBERN
       KORD=1
       NLDF=4
       NPL =1
       NG  =3
       IF (IFAM.LT.100) THEN
        WRITE(CFAM,'(A,I2,A)') 'Q1~ (E0',IFAM,')/Q0'
       ELSE
        WRITE(CFAM,'(A,I2,A)') 'Q1~ (EM',IFAM-100,')/Q0'
       ENDIF
      ELSE
       WRITE(MTERM,*) 'BNSTKT: IBERN must be 2..',NNBK,
     *                ' or -30,-31,-130,-131'
       RETURN
      ENDIF
      THETA=1D0-SQRT(0.5D0)
      FALPHA=(1D0-2D0*THETA)/(1D0-THETA)
C
      DO 200 ILEV=NLMIN,NLMAX
       CALL SETLEV(1)
       IF (IFAM.EQ.0) THEN
        NU=NVT+(KORD-1)*NMT+(KORD-1)**2*NEL
       ELSE
        NU=NMT
       ENDIF
       NP=NPL*NEL
C
C *** dofs, boundary marker, matrices
       CALL ZNEW(NLDF*NEL,3,LDFG,'KDFG  ')
       CALL ZNEW(NPL *NEL,3,LDFP,'KDFP  ')
       CALL ZNEW(NU      ,3,LBD ,'KBD   ')
       IF (IER.NE.0) STOP 'BNSTKT: not enough workspace'
       IF (IFAM.EQ.0) THEN
        CALL BNDFMP(KORD,KWORK(L(LVERT)),KWORK(L(LMID)),NEL,NVT,NMT,
     *              KWORK(L(LDFG)))
        CALL BNBDMK(KORD,KWORK(L(LDFG)),KWORK(L(LADJ)),NEL,NU,
     *              KWORK(L(LBD)))
       ELSE
        CALL BQDFMP(KWORK(L(LMID)),KWORK(L(LADJ)),NEL,NVT,NMT,
     *              KWORK(L(LDFG)),KWORK(L(LBD)))
       ENDIF
       CALL BNDFPP(NPL,NEL,KWORK(L(LDFP)))
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NLDF,KWORK(L(LDFG)),NU,
     *            .TRUE.,LLDA,LCOLA,NA)
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NPL,KWORK(L(LDFP)),NU,
     *            .FALSE.,LLDB,LCOLB,NB)
       CALL ZNEW(NA,1,LM  ,'DMBN  ')
       CALL ZNEW(NA,1,LA  ,'DABN  ')
       CALL ZNEW(NA,1,LAH ,'DAHBN ')
       CALL ZNEW(NB,1,LB1 ,'DB1BN ')
       CALL ZNEW(NB,1,LB2 ,'DB2BN ')
       CALL ZNEW(NU,1,LML ,'DMLBN ')
       CALL ZNEW(NPL*NPL*NEL,1,LMPI,'DMPIBN')
       CALL ZNEW(NP,1,LWP ,'DWPBN ')
       CALL ZNEW(NP,1,LCD ,'DCDBN ')
       CALL ZNEW(NU,1,LUD1,'DUD1BN')
       CALL ZNEW(NU,1,LUD2,'DUD2BN')
       CALL ZNEW(NU,1,LF1 ,'DF1BN ')
       CALL ZNEW(NU,1,LF2 ,'DF2BN ')
       CALL ZNEW(NU,1,LU1 ,'DU1BN ')
       CALL ZNEW(NU,1,LU2 ,'DU2BN ')
       CALL ZNEW(NP,1,LP  ,'DPBN  ')
       CALL ZNEW(8*NU,1,LAUX ,'DAUXBN')
       CALL ZNEW(8*NP,1,LPAUX,'DPAUXB')
       IF (IER.NE.0) STOP 'BNSTKT: not enough workspace'
       CALL BNASM(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),KWORK(L(LDFP)),
     *            KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *            KWORK(L(LLDB)),KWORK(L(LCOLB)),
     *            DWORK(L(LM)),DWORK(L(LA)),DWORK(L(LB1)),DWORK(L(LB2)),
     *            DAREA)
C *** lumped mass Ml = row sums of M (> 0 for Bernstein and Q1~)
       DO 20 IR=1,NU
        S=0D0
        DO 15 IA=KWORK(L(LLDA)+IR-1),KWORK(L(LLDA)+IR)-1
15      S=S+DWORK(L(LM)+IA-1)
        DWORK(L(LML)+IR-1)=S
20     CONTINUE
       CALL BNMPI(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            DWORK(L(LMPI)),DWORK(L(LWP)))
       AREA=0D0
       DO 25 IQ=1,NP
25     AREA=AREA+DWORK(L(LWP)+IQ-1)
       CALL BNCDG(DWORK(L(LB1)),DWORK(L(LB2)),KWORK(L(LLDB)),
     *            KWORK(L(LCOLB)),NU,NP,DWORK(L(LML)),KWORK(L(LBD)),
     *            DWORK(L(LCD)))
C *** Dirichlet values of the exact velocity
       IF (IFAM.EQ.0) THEN
        CALL BNBDV(KORD,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LADJ)),KWORK(L(LDFG)),NEL,SKU1,
     *             DWORK(L(LUD1)))
        CALL BNBDV(KORD,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LADJ)),KWORK(L(LDFG)),NEL,SKU2,
     *             DWORK(L(LUD2)))
       ELSE
        CALL BQBDV(MOD(IFAM,100),DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LMID)),KWORK(L(LADJ)),NEL,NVT,SKU1,
     *             DWORK(L(LUD1)))
        CALL BQBDV(MOD(IFAM,100),DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LMID)),KWORK(L(LADJ)),NEL,NVT,SKU2,
     *             DWORK(L(LUD2)))
       ENDIF
C
       DO 100 ICASE=1,2
        IF (ICASE.EQ.1) THEN
         DSA =0D0
         DSNU=1D0
        ELSE
         DSA =1D0/(THETA*0.01D0)
         DSNU=FALPHA*1D-3
        ENDIF
        DO 30 IA=1,NA
30      DWORK(L(LAH)+IA-1)=DSA*DWORK(L(LM)+IA-1)
     *                    +DSNU*DWORK(L(LA)+IA-1)
        CALL LCL1(DWORK(L(LF1)),NU)
        CALL LCL1(DWORK(L(LF2)),NU)
        CALL BNRHS(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),SKF1,DWORK(L(LF1)))
        CALL BNRHS(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),SKF2,DWORK(L(LF2)))
        CALL LCP1(DWORK(L(LUD1)),DWORK(L(LU1)),NU)
        CALL LCP1(DWORK(L(LUD2)),DWORK(L(LU2)),NU)
C
        CALL ZTIME(T0)
        CALL BNUZAW(NU,NP,NEL,DWORK(L(LAH)),KWORK(L(LLDA)),
     *       KWORK(L(LCOLA)),DWORK(L(LB1)),DWORK(L(LB2)),
     *       KWORK(L(LLDB)),KWORK(L(LCOLB)),KWORK(L(LBD)),
     *       DWORK(L(LML)),DWORK(L(LMPI)),DWORK(L(LWP)),DWORK(L(LCD)),
     *       AREA,DSA,DSNU,DWORK(L(LF1)),DWORK(L(LF2)),
     *       DWORK(L(LU1)),DWORK(L(LU2)),DWORK(L(LP)),
     *       DWORK(L(LAUX)),DWORK(L(LPAUX)),1D-10,500,
     *       ITUZ,ITV,ITC,RES0,RES)
        CALL ZTIME(T1)
C
        CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),DWORK(L(LU1)),SKU1,SKU1X,SKU1Y,
     *             E1L2,E1H1)
        CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),DWORK(L(LU2)),SKU2,SKU2X,SKU2Y,
     *             E2L2,E2H1)
        CALL BNPERR(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *              KWORK(L(LDFP)),DWORK(L(LP)),SKP,EPL2)
        ERR(ILEV,ICASE,1)=SQRT(E1L2**2+E2L2**2)
        ERR(ILEV,ICASE,2)=SQRT(E1H1**2+E2H1**2)
        ERR(ILEV,ICASE,3)=EPL2
C *** best approximation |p-Pi_h p|L2 (L2 projection, uses DPAUX)
        CALL BNPL2(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),SKP,
     *             DWORK(L(LPAUX)),EPBA)
        ERR(ILEV,ICASE,4)=EPBA
C
        DO 90 IU=1,2
         IF (IU.EQ.1) MU=MTERM
         IF (IU.EQ.2) MU=MFILE
         IF ((ILEV.EQ.NLMIN).AND.(ICASE.EQ.1))
     *    WRITE(MU,1000) CFAM,THETA,FALPHA
         IF (ICASE.EQ.1) WRITE(MU,1010) ILEV,NEL,NU,NP,NA,NB
         WRITE(MU,1020) ICASE,ITUZ,DBLE(ITV)/DBLE(2*ABS(ITUZ)+2),
     *    DBLE(ITC)/DBLE(ABS(ITUZ)+1),RES/RES0,T1-T0,
     *    (ERR(ILEV,ICASE,J),J=1,4)
         IF (ILEV.GT.NLMIN) WRITE(MU,1030)
     *    (LOG(ERR(ILEV-1,ICASE,J)/ERR(ILEV,ICASE,J))/LOG(2D0),J=1,4)
90      CONTINUE
100    CONTINUE
C
       CALL ZDISP(0,LPAUX,'DPAUXB')
       CALL ZDISP(0,LAUX ,'DAUXBN')
       CALL ZDISP(0,LP   ,'DPBN  ')
       CALL ZDISP(0,LU2  ,'DU2BN ')
       CALL ZDISP(0,LU1  ,'DU1BN ')
       CALL ZDISP(0,LF2  ,'DF2BN ')
       CALL ZDISP(0,LF1  ,'DF1BN ')
       CALL ZDISP(0,LUD2 ,'DUD2BN')
       CALL ZDISP(0,LUD1 ,'DUD1BN')
       CALL ZDISP(0,LCD  ,'DCDBN ')
       CALL ZDISP(0,LWP  ,'DWPBN ')
       CALL ZDISP(0,LMPI ,'DMPIBN')
       CALL ZDISP(0,LML  ,'DMLBN ')
       CALL ZDISP(0,LB2  ,'DB2BN ')
       CALL ZDISP(0,LB1  ,'DB1BN ')
       CALL ZDISP(0,LAH  ,'DAHBN ')
       CALL ZDISP(0,LA   ,'DABN  ')
       CALL ZDISP(0,LM   ,'DMBN  ')
       CALL ZDISP(0,LCOLB,'KCOLBN')
       CALL ZDISP(0,LLDB ,'KLDBN ')
       CALL ZDISP(0,LCOLA,'KCOLBN')
       CALL ZDISP(0,LLDA ,'KLDBN ')
       CALL ZDISP(0,LBD  ,'KBD   ')
       CALL ZDISP(0,LDFP ,'KDFP  ')
       CALL ZDISP(0,LDFG ,'KDFG  ')
200   CONTINUE
C
1000  FORMAT(/' STOKES TEST  ',A16,
     * '   case 1: a=0, nu=1;  case 2: a=1/(theta K), K=0.01,',
     * ' nu=alpha/1000'/'   theta=',F8.5,' alpha=',F8.5/
     * ' case  UZ  vel.PCG/solve  CG(C)/prec  |div|/|div0|',
     * '  CPU   |u-uh|L2   |u-uh|H1   |p-ph|L2 |p-Pp|L2'/' ',
     * 110('-'))
1010  FORMAT(' LEV',I2,' NEL=',I7,' NU=',I8,' (x2)  NP=',I8,
     * '  nnz(A)=',I9,' nnz(B)=',I9)
1020  FORMAT(I4,I5,F12.1,F12.1,D14.3,F8.2,4D11.3)
1030  FORMAT('   rates',47X,4F11.2)
C
      END
