************************************************************************
* bnafc.f  -  AFC (algebraic flux correction) for Bernstein Qk,
*             scalar steady convection-diffusion test (IBNTST=4)
*
*   -eps Lap u + b.grad u = f   on the c2d channel with cylinder,
*   b = potential flow past the cylinder (centre (0.5,0.2), R=0.05),
*   div b = 0, b.n = 0 on the cylinder.
*
* Schemes (all with the same Bernstein Qk space):
*   GAL  Galerkin:            G u = f,  G = eps*A + K
*   LOW  discrete upwinding:  L u = f,  L = G - D, D symmetric,
*        d_ij = max(0,g_ij,g_ji) (i/=j), d_ii = -sum_j d_ij
*        -> L has zero row sums (A*1 = K*1 = 0) and nonpositive
*           off-diagonals (M-matrix with Dirichlet rows)
*   AFC  L u = f + sum_j alpha_ij f_ij,  f_ij = d_ij (u_i - u_j),
*        Zalesak-type limiter with local bounds
*        u_i^max/min = max/min of u over the stencil of i,
*        Q_i^+- = CQ * (sum_j d_ij) * (u_i^max/min - u_i),
*        R_i^+- = min(1, Q_i^+-/P_i^+-), alpha_ij = min(R_i^+,R_j^-)
*        for f_ij>0 (min(R_i^-,R_j^+) otherwise); fixed-point
*        iteration L u^(m+1) = f + F(u^m) with relaxation OMEGA.
*        Since the Bernstein basis is nonnegative with partition of
*        unity, bounds on the coefficients bound the function.
*
* Test A (smooth): exact u = sin(3x)cos(2y)+xy, Dirichlet on the whole
*   boundary (L2-projection coefficients), errors L2/H1 and rates.
* Test B (layer):  f = 0, Dirichlet only at the inflow x=0 with
*   u = 1 for y>0.3, 0 otherwise (coefficients from the value at the
*   Greville point), Neumann elsewhere. Exact solution in [0,1];
*   reported: min/max of the coefficients and of u at the cubature
*   points, smearing  W = int 4 v(1-v), v = min(1,max(0,u)).
************************************************************************
C
************************************************************************
*   BNBVEL: potential flow past the cylinder, U=1
************************************************************************
      SUBROUTINE BNBVEL(X,Y,W1,W2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (XC=0.5D0,YC=0.2D0,RC=0.05D0)
      XX=X-XC
      YY=Y-YC
      R2=XX*XX+YY*YY
      IF (R2.LT.1D-12) R2=1D-12
      W1=1D0-RC*RC*(XX*XX-YY*YY)/(R2*R2)
      W2=   -RC*RC*2D0*XX*YY/(R2*R2)
      END
C
************************************************************************
*   Right hand side of test A: f = -eps Lap u + b.grad u,
*   u = BNUEX (bern.f), -Lap u = BNFEX
************************************************************************
      DOUBLE PRECISION FUNCTION AFFEX(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BNAFCP/ EPSD
      DOUBLE PRECISION BNFEX,BNUEXX,BNUEXY
      EXTERNAL BNFEX,BNUEXX,BNUEXY
      CALL BNBVEL(X,Y,W1,W2)
      AFFEX=EPSD*BNFEX(X,Y)+W1*BNUEXX(X,Y)+W2*BNUEXY(X,Y)
      END
C
************************************************************************
*   BNCONV: convection matrix DK(ij) += int phi_i b.grad phi_j
************************************************************************
      SUBROUTINE BNCONV(NEL,DCORVG,KVERT,KDFG,KLDA,KCOLA,DK)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG(NLDF,*)
      DIMENSION KLDA(*),KCOLA(*),DK(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      DIMENSION DEK(NNLD,NNLD)
      INTEGER*8 BNFIND
C
      CALL BNGAUS(NG,XG,WG)
      DO 100 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
       DO 10 J=1,NLDF
       DO 10 I=1,NLDF
10     DEK(I,J)=0D0
       DO 30 IG=1,NG*NG
        CALL BNBVEL(DXQ(IG),DYQ(IG),W1,W2)
        DO 20 J=1,NLDF
         C=DWT(IG)*(W1*DPX(J,IG)+W2*DPY(J,IG))
         DO 20 I=1,NLDF
20       DEK(I,J)=DEK(I,J)+C*DPHI(I,IG)
30     CONTINUE
       DO 50 I=1,NLDF
        IR=KDFG(I,IEL)
        DO 40 J=1,NLDF
         IA=BNFIND(KLDA,KCOLA,IR,KDFG(J,IEL))
         IF (IA.EQ.0) STOP 'BNCONV: entry not in structure'
         DK(IA)=DK(IA)+DEK(I,J)
40      CONTINUE
50     CONTINUE
100   CONTINUE
      END
C
************************************************************************
*   BNDFXY: Greville point of every global dof (local (i,j) ->
*           reference (-1+2i/k,-1+2j/k), bilinear map) = Bernstein
*           coefficients of x, y; curved elements (IISO, bniso.f):
*           the coefficients of the geometry map (x, y in Qk exactly)
************************************************************************
      SUBROUTINE BNDFXY(K,NEL,DCORVG,KVERT,KDFG,DXD,DYD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG((K+1)**2,*)
      DIMENSION DXD(*),DYD(*),DXV(4),DYV(4)
      COMMON /BNCUR/  ICUR,MGC,CGX(16),CGY(16)
      DO 20 IEL=1,NEL
       CALL BNGSET(IEL)
       IF ((ICUR.EQ.1).AND.(MGC.EQ.K)) THEN
        DO 3 L=1,(K+1)**2
         DXD(KDFG(L,IEL))=CGX(L)
         DYD(KDFG(L,IEL))=CGY(L)
3       CONTINUE
        GOTO 20
       ENDIF
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       DO 10 J=0,K
       DO 10 I=0,K
        XI =-1D0+2D0*DBLE(I)/DBLE(K)
        ETA=-1D0+2D0*DBLE(J)/DBLE(K)
        CALL BNMAP(XI,ETA,DXV,DYV,X,Y,D11,D12,D21,D22,DET)
        IG=KDFG(I+J*(K+1)+1,IEL)
        DXD(IG)=X
        DYD(IG)=Y
10     CONTINUE
20    CONTINUE
      CALL BNGSET(0)
      END
C
************************************************************************
*   BNTRIX: KT(ia) = position of the transposed entry (j,i)
************************************************************************
      SUBROUTINE BNTRIX(KLD,KCOL,N,KT)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KLD(*),KCOL(*),KT(*)
      INTEGER*8 BNFIND
      DO 20 IR=1,N
      DO 20 IA=KLD(IR),KLD(IR+1)-1
       KT(IA)=BNFIND(KLD,KCOL,KCOL(IA),IR)
       IF (KT(IA).EQ.0) STOP 'BNTRIX: structure not symmetric'
20    CONTINUE
      END
C
************************************************************************
*   BNDISC: discrete upwinding. DD(ia) = d_ij >= 0 for off-diagonal
*           entries (0 on the diagonal), DL = DG - D (D has zero row
*           sums), DQS(i) = sum_j d_ij
************************************************************************
      SUBROUTINE BNDISC(DG,KLD,KCOL,KT,N,DD,DL,DQS)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DG(*),KLD(*),KCOL(*),KT(*),DD(*),DL(*),DQS(*)
      DO 20 IR=1,N
       S=0D0
       DD(KLD(IR))=0D0
       DO 10 IA=KLD(IR)+1,KLD(IR+1)-1
        D=MAX(0D0,DG(IA),DG(KT(IA)))
        DD(IA)=D
        DL(IA)=DG(IA)-D
        S=S+D
10     CONTINUE
       DL(KLD(IR))=DG(KLD(IR))+S
       DQS(IR)=S
20    CONTINUE
      END
C
************************************************************************
*   BNDIRR: Dirichlet rows (KBD=1) of a matrix -> identity rows
************************************************************************
      SUBROUTINE BNDIRR(DA,KLD,N,KBD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),KLD(*),KBD(*)
      DO 20 IR=1,N
       IF (KBD(IR).NE.1) GOTO 20
       DA(KLD(IR))=1D0
       DO 10 IA=KLD(IR)+1,KLD(IR+1)-1
10     DA(IA)=0D0
20    CONTINUE
      END
C
************************************************************************
*   BNLIM: limited antidiffusion DF(i) = sum_j alpha_ij f_ij,
*          f_ij = d_ij (u_i - u_j); returns ALAVG = mean alpha over
*          the fluxes with d_ij>0, NLIM = number with alpha<0.999
*          DRP,DRM: work (R^+, R^-)
************************************************************************
      SUBROUTINE BNLIM(DU,DD,DQS,KLD,KCOL,N,KBD,CQ,DF,DRP,DRM,
     *                 ALAVG,NLIM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DU(*),DD(*),DQS(*),KLD(*),KCOL(*),KBD(*),DF(*)
      DIMENSION DRP(*),DRM(*)
C
      DO 20 IR=1,N
       UI=DU(IR)
       UMX=UI
       UMN=UI
       PP=0D0
       PM=0D0
       DO 10 IA=KLD(IR)+1,KLD(IR+1)-1
        UJ=DU(KCOL(IA))
        UMX=MAX(UMX,UJ)
        UMN=MIN(UMN,UJ)
        F=DD(IA)*(UI-UJ)
        PP=PP+MAX(0D0,F)
        PM=PM+MIN(0D0,F)
10     CONTINUE
       QP=CQ*DQS(IR)*(UMX-UI)
       QM=CQ*DQS(IR)*(UMN-UI)
       DRP(IR)=1D0
       DRM(IR)=1D0
C *** R = min(1,Q/P) without overflow for tiny P
       IF (QP.LT.PP) DRP(IR)=QP/PP
       IF (QM.GT.PM) DRM(IR)=QM/PM
       IF (KBD(IR).EQ.1) THEN
        DRP(IR)=1D0
        DRM(IR)=1D0
       ENDIF
20    CONTINUE
C
      ALSUM=0D0
      NFL=0
      NLIM=0
      DO 40 IR=1,N
       S=0D0
       DO 30 IA=KLD(IR)+1,KLD(IR+1)-1
        IF (DD(IA).LE.0D0) GOTO 30
        JC=KCOL(IA)
        F=DD(IA)*(DU(IR)-DU(JC))
        IF (F.GT.0D0) THEN
         ALPHA=MIN(DRP(IR),DRM(JC))
        ELSE
         ALPHA=MIN(DRM(IR),DRP(JC))
        ENDIF
        S=S+ALPHA*F
        ALSUM=ALSUM+ALPHA
        NFL=NFL+1
        IF (ALPHA.LT.0.999D0) NLIM=NLIM+1
30     CONTINUE
       DF(IR)=S
       IF (KBD(IR).EQ.1) DF(IR)=0D0
40    CONTINUE
      ALAVG=1D0
      IF (NFL.GT.0) ALAVG=ALSUM/DBLE(NFL)
      END
C
************************************************************************
*   BNSSRP: z = M^(-1) r, M = SSOR(1) of the CSR matrix DP (diagonal
*           first, columns in any order)
************************************************************************
      SUBROUTINE BNSSRP(DP,KLD,KCOL,N,DR,DZ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DP(*),KLD(*),KCOL(*),DR(*),DZ(*)
      DO 20 IR=1,N
       S=DR(IR)
       DO 10 IA=KLD(IR)+1,KLD(IR+1)-1
10     IF (KCOL(IA).LT.IR) S=S-DP(IA)*DZ(KCOL(IA))
       DZ(IR)=S/DP(KLD(IR))
20    CONTINUE
      DO 40 IR=N,1,-1
       S=0D0
       DO 30 IA=KLD(IR)+1,KLD(IR+1)-1
30     IF (KCOL(IA).GT.IR) S=S+DP(IA)*DZ(KCOL(IA))
       DZ(IR)=DZ(IR)-S/DP(KLD(IR))
40    CONTINUE
      END
C
************************************************************************
*   BNBICG: right-preconditioned BiCGStab for DA x = DB, preconditioner
*           SSOR(DP) (DP may differ from DA). DW(N,8) work.
*           ITE<0: no convergence (|r| <= EPS |r0|) within ITMAX
************************************************************************
      SUBROUTINE BNBICG(DA,DP,KLD,KCOL,N,DX,DB,DW,ITMAX,EPS,ITE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),DP(*),KLD(*),KCOL(*),DX(*),DB(*),DW(N,8)
C     columns: 1 r, 2 r0hat, 3 p, 4 v, 5 phat, 6 s, 7 shat, 8 t
C
      CALL BNAX(DA,KLD,KCOL,N,DX,DW(1,1))
      RN0=0D0
      DO 10 I=1,N
       DW(I,1)=DB(I)-DW(I,1)
       DW(I,2)=DW(I,1)
       DW(I,3)=0D0
       DW(I,4)=0D0
       RN0=RN0+DW(I,1)**2
10    CONTINUE
      RN0=SQRT(RN0)
      ITE=0
      IF (RN0.EQ.0D0) RETURN
      RHO=1D0
      ALPHA=1D0
      OM=1D0
      DO 100 ITE=1,ITMAX
       RHO1=0D0
       DO 15 I=1,N
15     RHO1=RHO1+DW(I,2)*DW(I,1)
       IF (RHO1.EQ.0D0) GOTO 200
       GBET=(RHO1/RHO)*(ALPHA/OM)
       RHO=RHO1
       DO 20 I=1,N
20     DW(I,3)=DW(I,1)+GBET*(DW(I,3)-OM*DW(I,4))
       CALL BNSSRP(DP,KLD,KCOL,N,DW(1,3),DW(1,5))
       CALL BNAX(DA,KLD,KCOL,N,DW(1,5),DW(1,4))
       S=0D0
       DO 25 I=1,N
25     S=S+DW(I,2)*DW(I,4)
       IF (S.EQ.0D0) GOTO 200
       ALPHA=RHO/S
       SN=0D0
       DO 30 I=1,N
        DW(I,6)=DW(I,1)-ALPHA*DW(I,4)
        SN=SN+DW(I,6)**2
30     CONTINUE
       IF (SQRT(SN).LE.EPS*RN0) THEN
        DO 35 I=1,N
35      DX(I)=DX(I)+ALPHA*DW(I,5)
        RETURN
       ENDIF
       CALL BNSSRP(DP,KLD,KCOL,N,DW(1,6),DW(1,7))
       CALL BNAX(DA,KLD,KCOL,N,DW(1,7),DW(1,8))
       TT=0D0
       TS=0D0
       DO 40 I=1,N
        TT=TT+DW(I,8)**2
        TS=TS+DW(I,8)*DW(I,6)
40     CONTINUE
       IF (TT.EQ.0D0) GOTO 200
       OM=TS/TT
       RN=0D0
       DO 50 I=1,N
        DX(I)=DX(I)+ALPHA*DW(I,5)+OM*DW(I,7)
        DW(I,1)=DW(I,6)-OM*DW(I,8)
        RN=RN+DW(I,1)**2
50     CONTINUE
       IF (SQRT(RN).LE.EPS*RN0) RETURN
       IF (OM.EQ.0D0) GOTO 200
100   CONTINUE
200   ITE=-ABS(ITE)
      END
C
************************************************************************
*   BNAFCN: AFC nonlinear iteration for u = g(u), g(u) = L^(-1)(f+F(u))
*     MAA = 0: relaxed fixed point  u <- u + OMEGA (g(u) - u)
*     MAA > 0: Anderson acceleration with history depth MAA and
*              damping AABET:  f_k = g(u_k) - u_k,
*              gamma = argmin |f_k - dF gamma|  (DGELS, dF = last MAA
*              differences of f, dX = differences of u),
*              u <- u + AABET f_k - (dX + AABET dF) gamma;
*              history cleared if the residual exceeds 10x its minimum
*   Stops at |f + F(u) - L u| <= EPSNL |f + F(u^0) - L u^0| or after
*   NITMX steps. DL has identity Dirichlet rows, DRHS holds f with the
*   Dirichlet values; AB/IPIV: banded LU of DL (BNBFAC). DW(N,11) and
*   DAW(N*(3*MAA+2)+LWK) work, LWK=N+64*(MAA+1). Returns NIT (<0: no
*   conv.), ALAVG, NLIM of the last limiter call, LINTOT = number of
*   solves, RNL = final relative residual, NRST = Anderson restarts.
************************************************************************
      SUBROUTINE BNAFCN(DL,KLD,KCOL,N,KBD,DD,DQS,CQ,DU,DRHS,DW,
     *                  AB,LDAB,KBW,IPIV,KINV,MAA,AABET,DAW,NRST,
     *                  OMEGA,EPSNL,NITMX,NIT,ALAVG,NLIM,LINTOT,RNL)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DL(*),KLD(*),KCOL(*),KBD(*),DD(*),DQS(*)
      DIMENSION DU(*),DRHS(*),DW(N,11),AB(*),IPIV(*),KINV(*),DAW(*)
C     DW(.,1) solve work, 6/7 R+/R-, 9 F, 10 rhs / f_k, 11 g(u)
C     DAW: dF(N,MAA) | dX(N,MAA) | copy of dF | f_old | x_old | LS work
      IFH=1
      IXH=IFH+N*MAA
      IAC=IXH+N*MAA
      IFO=IAC+N*MAA
      IXO=IFO+N
      IWK=IXO+N
      LWK=N+64*(MAA+1)
      NHIST=0
      BOLD=.FALSE.
      NRST=0
      LINTOT=0
      RNL0=-1D0
      RNLMN=1D300
      DO 100 NIT=0,NITMX
       CALL BNLIM(DU,DD,DQS,KLD,KCOL,N,KBD,CQ,DW(1,9),DW(1,6),
     *            DW(1,7),ALAVG,NLIM)
       CALL BNAX(DL,KLD,KCOL,N,DU,DW(1,11))
       RNL=0D0
       DO 10 I=1,N
        DW(I,10)=DRHS(I)+DW(I,9)
        RNL=RNL+(DW(I,10)-DW(I,11))**2
10     CONTINUE
       RNL=SQRT(RNL)
       IF (RNL0.LT.0D0) RNL0=MAX(RNL,1D-300)
       IF (RNL.LE.EPSNL*RNL0) THEN
        RNL=RNL/RNL0
        RETURN
       ENDIF
       IF (NIT.EQ.NITMX) GOTO 200
       CALL BNBSOL(AB,LDAB,N,KBW,IPIV,KINV,DW(1,10),DW(1,11),DW)
       LINTOT=LINTOT+1
C
       IF (MAA.LE.0) THEN
        DO 30 I=1,N
30      DU(I)=DU(I)+OMEGA*(DW(I,11)-DU(I))
        GOTO 100
       ENDIF
C
C *** Anderson: restart if the residual grew strongly
       IF (RNL.GT.1D1*RNLMN) THEN
        NHIST=0
        BOLD=.FALSE.
        NRST=NRST+1
       ENDIF
       RNLMN=MIN(RNLMN,RNL)
C *** f_k = g(u_k) - u_k in DW(.,10)
       DO 40 I=1,N
40     DW(I,10)=DW(I,11)-DU(I)
C *** new history column (ring buffer)
       IF (BOLD) THEN
        JC=MOD(NHIST,MAA)
        DO 45 I=1,N
         DAW(IFH+JC*N+I-1)=DW(I,10)-DAW(IFO+I-1)
         DAW(IXH+JC*N+I-1)=DU(I)   -DAW(IXO+I-1)
45      CONTINUE
        NHIST=NHIST+1
       ENDIF
       DO 50 I=1,N
        DAW(IFO+I-1)=DW(I,10)
        DAW(IXO+I-1)=DU(I)
50     CONTINUE
       BOLD=.TRUE.
       MK=MIN(NHIST,MAA)
C *** least squares  min |f_k - dF gamma|, gamma -> DW(1..MK,11)
       IF (MK.GT.0) THEN
        DO 55 I=1,N*MK
55      DAW(IAC+I-1)=DAW(IFH+I-1)
        DO 56 I=1,N
56      DW(I,11)=DW(I,10)
        CALL DGELS('N',N,MK,1,DAW(IAC),N,DW(1,11),N,DAW(IWK),LWK,
     *             INFO)
        IF (INFO.NE.0) THEN
         NHIST=0
         NRST=NRST+1
         MK=0
        ENDIF
       ENDIF
C *** update u <- u + beta f_k - sum_j gamma_j (dX_j + beta dF_j)
       DO 60 I=1,N
60     DU(I)=DU(I)+AABET*DW(I,10)
       DO 70 J=1,MK
        GAM=DW(J,11)
        DO 65 I=1,N
65      DU(I)=DU(I)-GAM*(DAW(IXH+(J-1)*N+I-1)
     *                  +AABET*DAW(IFH+(J-1)*N+I-1))
70     CONTINUE
100   CONTINUE
200   NIT=-NITMX
      RNL=RNL/RNL0
      END
C
************************************************************************
*   BNRNG: min/max of u at the cubature points and smearing
*          W = int 4 v(1-v), v = min(1,max(0,u))
************************************************************************
      SUBROUTINE BNRNG(NEL,DCORVG,KVERT,KDFG,DC,UMN,UMX,WSM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG(NLDF,*),DC(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      CALL BNGAUS(NG,XG,WG)
      UMN= 1D30
      UMX=-1D30
      WSM=0D0
      DO 30 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
       DO 20 IG=1,NG*NG
        U=0D0
        DO 10 I=1,NLDF
10      U=U+DC(KDFG(I,IEL))*DPHI(I,IG)
        UMN=MIN(UMN,U)
        UMX=MAX(UMX,U)
        W=MIN(1D0,MAX(0D0,U))
        WSM=WSM+DWT(IG)*4D0*W*(1D0-W)
20     CONTINUE
30    CONTINUE
      END
C
************************************************************************
*   BNAFCT: driver, levels NLMIN..NLMAX, eps = 1e-3 and 1e-6
************************************************************************
      SUBROUTINE BNAFCT(K,MFILE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNBK=8,NNEPS=2)
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
      COMMON /BNAFCP/ EPSD
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      DIMENSION EPSV(NNEPS)
      DIMENSION EAL2(3,NNLEV,NNEPS),EAH1(3,NNLEV,NNEPS)
      PARAMETER (NNGAM=8)
      DIMENSION GAMV(NNGAM)
      DIMENSION EJL2(NNGAM,NNLEV,NNEPS),EJH1(NNGAM,NNLEV,NNEPS)
      CHARACTER*4 CSCH(3)
      CHARACTER*60 CGMV
      EXTERNAL BNUEX,BNUEXX,BNUEXY,AFFEX
      DOUBLE PRECISION BNUEX,BNUEXX,BNUEXY,AFFEX
      DATA EPSV /1D-3,1D-6/
      DATA CSCH /'GAL ','LOW ','AFC '/
      DATA GAMV /1D-3,1D-2,1D-1,1D0,4*0D0/
      SAVE
C
      SUB='BNAFCT'
      IF ((K.LT.2).OR.(K.GT.NNBK)) THEN
       WRITE(MTERM,*) 'BNAFCT: order K must be in 2..',NNBK
       RETURN
      ENDIF
      IFAM=0
      KORD=K
      NLDF=(K+1)**2
      NPL =K*(K+1)/2
      NG  =K+2
      CQ   =1D0
      OMEGA=1D0
      EPSNL=1D-8
      NITMX=500
      MAA  =0
      AABET=1D0
      IAFC =1
      IGMV =0
      NGAMR=4
C *** optional file afc.par: CQ OMEGA EPSNL NITMX [MAA AABET [IAFC
C     [IGMV [NGAMR GAMV(1..NGAMR)]]]]; IAFC=0: skip the AFC solve
C     (jump runs); IGMV>0: GMV files of all solutions on level IGMV in
C     #gmv/; NGAMR<=NNGAM jump parameters (default 1e-3,1e-2,0.1,1)
      OPEN(77,FILE='afc.par',STATUS='OLD',ERR=1)
      READ(77,*,END=1,ERR=1) CQ,OMEGA,EPSNL,NITMX,MAA,AABET,IAFC,IGMV,
     *                       NGAMR,(GAMV(I),I=1,MIN(NGAMR,NNGAM))
1     CLOSE(77)
      NGAMR=MAX(1,MIN(NGAMR,NNGAM))
C
      DO 2 IU=1,2
       IF (IU.EQ.1) MU=MTERM
       IF (IU.EQ.2) MU=MFILE
       WRITE(MU,1000) K,CQ,OMEGA,EPSNL,NITMX,MAA,AABET
       WRITE(MU,1001) (GAMV(I),I=1,NGAMR)
2     CONTINUE
C
      DO 100 ILEV=NLMIN,NLMAX
       CALL SETLEV(1)
       NDOF=NVT+(K-1)*NMT+(K-1)**2*NEL
C
       CALL ZNEW(NLDF*NEL,3,LDFG,'KDFG  ')
       CALL ZNEW(NPL *NEL,3,LDFP,'KDFP  ')
       CALL ZNEW(NDOF    ,3,LBDA,'KBDA  ')
       CALL ZNEW(NDOF    ,3,LBDB,'KBDB  ')
       CALL ZNEW(NDOF    ,3,LNOB,'KNOB  ')
       CALL ZNEW(NDOF    ,1,LXD ,'DXD   ')
       CALL ZNEW(NDOF    ,1,LYD ,'DYD   ')
       IF (IER.NE.0) STOP 'BNAFCT: not enough workspace'
       CALL BNDFMP(K,KWORK(L(LVERT)),KWORK(L(LMID)),NEL,NVT,NMT,
     *             KWORK(L(LDFG)))
       CALL BNDFPP(NPL,NEL,KWORK(L(LDFP)))
       CALL BNBDMK(K,KWORK(L(LDFG)),KWORK(L(LADJ)),NEL,NDOF,
     *             KWORK(L(LBDA)))
       CALL BNDFXY(K,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),DWORK(L(LXD)),DWORK(L(LYD)))
C *** test B: Dirichlet only at the inflow x=0
       NINF=0
       DO 4 I=1,NDOF
        KWORK(L(LNOB)+I-1)=0
        KWORK(L(LBDB)+I-1)=0
        IF ((KWORK(L(LBDA)+I-1).EQ.1).AND.
     *      (DWORK(L(LXD)+I-1).LT.1D-8)) THEN
         KWORK(L(LBDB)+I-1)=1
         NINF=NINF+1
        ENDIF
4      CONTINUE
C
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NLDF,KWORK(L(LDFG)),NDOF,
     *            .TRUE.,LLDA,LCOLA,NA)
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NPL,KWORK(L(LDFP)),NDOF,
     *            .FALSE.,LLDB,LCOLB,NB)
       CALL ZNEW(NA,1,LM ,'DMAF  ')
       CALL ZNEW(NA,1,LA ,'DAAF  ')
       CALL ZNEW(NB,1,LB1,'DB1AF ')
       CALL ZNEW(NB,1,LB2,'DB2AF ')
       CALL ZNEW(NA,1,LK ,'DKAF  ')
       CALL ZNEW(NA,1,LG ,'DGAF  ')
       CALL ZNEW(NA,1,LLL,'DLAF  ')
       CALL ZNEW(NA,1,LDD,'DDAF  ')
       CALL ZNEW(NA,3,LKT,'KTAF  ')
       CALL ZNEW(NDOF,1,LQS ,'DQSAF ')
       CALL ZNEW(NDOF,1,LCP ,'DCPAF ')
       CALL ZNEW(NDOF,1,LRHS,'DRHSAF')
       CALL ZNEW(NDOF,1,LU  ,'DUAF  ')
       CALL ZNEW(11*NDOF,1,LW,'DWAF  ')
       CALL ZNEW(NDOF,3,LKP ,'KPAF  ')
       CALL ZNEW(NDOF,3,LKIN,'KINAF ')
       CALL ZNEW(NDOF,3,LIPV,'IPIVAF')
       IF (IER.NE.0) STOP 'BNAFCT: not enough workspace'
       CALL BNPERM(DWORK(L(LXD)),NDOF,KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *             KWORK(L(LKP)),KWORK(L(LKIN)),KBW)
C *** jump stabilization: interior edges, extended structure, A, K, J
       CALL ZNEW(5*4*NEL,3,LEDG,'KEDGJ ')
       IF (IER.NE.0) STOP 'BNAFCT: not enough workspace (edges)'
       CALL BNJEDG(KWORK(L(LADJ)),KWORK(L(LVERT)),NEL,KWORK(L(LEDG)),
     *             NEDG)
       NPAT=NEL+NEDG
       CALL ZNEW(2*NLDF*NPAT,3,LPAT,'KPATJ ')
       IF (IER.NE.0) STOP 'BNAFCT: not enough workspace (patches)'
       CALL BNJPAT(NLDF,KWORK(L(LDFG)),NEL,KWORK(L(LEDG)),NEDG,
     *             KWORK(L(LPAT)))
       CALL BNSTR(NPAT,2*NLDF,KWORK(L(LPAT)),2*NLDF,KWORK(L(LPAT)),
     *            NDOF,.TRUE.,LLDJ,LCOLJ,NAJ)
       CALL ZNEW(NAJ,1,LMJ ,'DMJ   ')
       CALL ZNEW(NAJ,1,LAJ ,'DAJ   ')
       CALL ZNEW(NAJ,1,LKJ ,'DKJ   ')
       CALL ZNEW(NAJ,1,LJJ ,'DJJ   ')
       CALL ZNEW(NAJ,1,LGJ ,'DGJ   ')
       CALL ZNEW(NB ,1,LB1J,'DB1J  ')
       CALL ZNEW(NB ,1,LB2J,'DB2J  ')
       IF (IER.NE.0) STOP 'BNAFCT: not enough workspace (jump)'
       CALL BNASM(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),KWORK(L(LDFP)),
     *            KWORK(L(LLDJ)),KWORK(L(LCOLJ)),
     *            KWORK(L(LLDB)),KWORK(L(LCOLB)),
     *            DWORK(L(LMJ)),DWORK(L(LAJ)),DWORK(L(LB1J)),
     *            DWORK(L(LB2J)),DAREAJ)
       CALL BNCONV(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),KWORK(L(LLDJ)),KWORK(L(LCOLJ)),
     *             DWORK(L(LKJ)))
       CALL BNJASM(K,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),KWORK(L(LEDG)),NEDG,
     *             KWORK(L(LLDJ)),KWORK(L(LCOLJ)),DWORK(L(LJJ)))
C *** check: J x = J y = 0 (Greville coordinates = coefficients of the
C     linear functions x, y on bilinear elements), scaled by max|J|
       DJMX=0D0
       DO 8 IA=1,NAJ
8      DJMX=MAX(DJMX,ABS(DWORK(L(LJJ)+IA-1)))
       DJLIN=0D0
       DO 9 IXY=1,2
        IF (IXY.EQ.1) LCRD=LXD
        IF (IXY.EQ.2) LCRD=LYD
        CALL BNAX(DWORK(L(LJJ)),KWORK(L(LLDJ)),KWORK(L(LCOLJ)),NDOF,
     *            DWORK(L(LCRD)),DWORK(L(LW)))
        DO 9 I=1,NDOF
         DJLIN=MAX(DJLIN,ABS(DWORK(L(LW)+I-1)))
9      CONTINUE
       DJLIN=DJLIN/MAX(DJMX,1D-300)
       CALL BNPERM(DWORK(L(LXD)),NDOF,KWORK(L(LLDJ)),KWORK(L(LCOLJ)),
     *             KWORK(L(LKP)),KWORK(L(LKIN)),KBWJ)
       LDAB=3*KBW+1
       LDABJ=3*KBWJ+1
       CALL ZNEW(MAX(LDAB,LDABJ)*NDOF,1,LAB,'DABAF ')
       CALL ZNEW(NDOF*(3*MAX(MAA,1)+2)+NDOF+64*(MAA+1),1,LAAW,'DAAWAF')
       IF (IER.NE.0) STOP 'BNAFCT: not enough workspace (band)'
       CALL BNASM(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),KWORK(L(LDFP)),
     *            KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *            KWORK(L(LLDB)),KWORK(L(LCOLB)),
     *            DWORK(L(LM)),DWORK(L(LA)),DWORK(L(LB1)),DWORK(L(LB2)),
     *            DAREA)
       CALL BNCONV(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFG)),KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *             DWORK(L(LK)))
       CALL BNTRIX(KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,KWORK(L(LKT)))
C *** max |K*1| (should vanish)
       CALL LCL1(DWORK(L(LU)),NDOF)
       DO 6 I=1,NDOF
6      DWORK(L(LU)+I-1)=1D0
       CALL BNAX(DWORK(L(LK)),KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,
     *           DWORK(L(LU)),DWORK(L(LCP)))
       DK1=0D0
       DO 7 I=1,NDOF
7      DK1=MAX(DK1,ABS(DWORK(L(LCP)+I-1)))
C *** L2 projection of the exact u (Dirichlet data of test A)
       CALL LCL1(DWORK(L(LRHS)),NDOF)
       CALL LCL1(DWORK(L(LCP)),NDOF)
       CALL BNRHS(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),BNUEX,DWORK(L(LRHS)))
       CALL BNCG(DWORK(L(LM)),KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,
     *           DWORK(L(LCP)),DWORK(L(LRHS)),KWORK(L(LNOB)),
     *           DWORK(L(LW)),DWORK(L(LW)+NDOF),DWORK(L(LW)+2*NDOF),
     *           10000,1D-13,ITP)
C *** best approximation: error of the L2 projection (test A)
       CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),DWORK(L(LCP)),BNUEX,BNUEXX,BNUEXY,
     *            EPL2,EPH1)
C
       DO 3 IU=1,2
        IF (IU.EQ.1) MU=MTERM
        IF (IU.EQ.2) MU=MFILE
        WRITE(MU,1010) ILEV,NEL,NDOF,NA,NINF,DK1,KBW
        WRITE(MU,1011) NEDG,NAJ,KBWJ,DJLIN
        WRITE(MU,1012) EPL2,EPH1
3      CONTINUE
C
       DO 90 IE=1,NNEPS
        EPSD=EPSV(IE)
        DO 80 ITST=1,2
C *** Galerkin matrix G = eps A + K, discrete upwinding
        DO 10 IA=1,NA
10      DWORK(L(LG)+IA-1)=EPSD*DWORK(L(LA)+IA-1)+DWORK(L(LK)+IA-1)
        CALL BNDISC(DWORK(L(LG)),KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *              KWORK(L(LKT)),NDOF,DWORK(L(LDD)),DWORK(L(LLL)),
     *              DWORK(L(LQS)))
        IF (ITST.EQ.1) LBD=LBDA
        IF (ITST.EQ.2) LBD=LBDB
        CALL BNDIRR(DWORK(L(LG)),KWORK(L(LLDA)),NDOF,KWORK(L(LBD)))
        CALL BNDIRR(DWORK(L(LLL)),KWORK(L(LLDA)),NDOF,KWORK(L(LBD)))
C *** right hand side with Dirichlet values
        CALL LCL1(DWORK(L(LRHS)),NDOF)
        IF (ITST.EQ.1)
     *   CALL BNRHS(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *              KWORK(L(LDFG)),AFFEX,DWORK(L(LRHS)))
        DO 15 I=1,NDOF
         IF (KWORK(L(LBD)+I-1).NE.1) GOTO 15
         IF (ITST.EQ.1) THEN
          DWORK(L(LRHS)+I-1)=DWORK(L(LCP)+I-1)
         ELSE
          G=0D0
          IF (DWORK(L(LYD)+I-1).GT.0.3D0+1D-12) G=1D0
          DWORK(L(LRHS)+I-1)=G
         ENDIF
15      CONTINUE
C
        DO 70 IS=1,3
         IF ((IS.EQ.3).AND.(IAFC.EQ.0)) GOTO 70
         CALL ZTIME(T0)
         CALL LCL1(DWORK(L(LU)),NDOF)
         DO 20 I=1,NDOF
20       IF (KWORK(L(LBD)+I-1).EQ.1)
     *    DWORK(L(LU)+I-1)=DWORK(L(LRHS)+I-1)
         NIT=0
         ITL=0
         ALAVG=1D0
         NLIM=0
         RNL=0D0
C *** direct banded LU: Galerkin G, low order L (factor of L reused
C     by AFC), checked by the residual of the solve
         IF (IS.LE.2) THEN
          IF (IS.EQ.1) LMAT=LG
          IF (IS.EQ.2) LMAT=LLL
          CALL BNBFAC(DWORK(L(LMAT)),KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *                NDOF,KWORK(L(LKIN)),KBW,DWORK(L(LAB)),LDAB,
     *                KWORK(L(LIPV)),INFO)
          IF (INFO.NE.0) WRITE(MTERM,*) 'BNAFCT: DGBTRF INFO=',INFO
         ENDIF
         CALL BNBSOL(DWORK(L(LAB)),LDAB,NDOF,KBW,KWORK(L(LIPV)),
     *               KWORK(L(LKIN)),DWORK(L(LRHS)),DWORK(L(LU)),
     *               DWORK(L(LW)))
         ITL=1
         IF (IS.EQ.3) THEN
C *** AFC, started from the low-order solution
          CALL BNAFCN(DWORK(L(LLL)),KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *                NDOF,KWORK(L(LBD)),DWORK(L(LDD)),DWORK(L(LQS)),
     *                CQ,DWORK(L(LU)),DWORK(L(LRHS)),DWORK(L(LW)),
     *                DWORK(L(LAB)),LDAB,KBW,KWORK(L(LIPV)),
     *                KWORK(L(LKIN)),MAA,AABET,DWORK(L(LAAW)),NRST,
     *                OMEGA,EPSNL,NITMX,NIT,ALAVG,NLIM,ITL,RNL)
          ITL=ITL+1
         ENDIF
C *** residual of the discrete problem solved (G or L u = f, AFC:
C     L u - f - F(u)), relative to |f|
         IF (IS.EQ.1) LMAT=LG
         IF (IS.GE.2) LMAT=LLL
         CALL BNAX(DWORK(L(LMAT)),KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,
     *             DWORK(L(LU)),DWORK(L(LW)+NDOF))
         IF (IS.EQ.3)
     *    CALL BNLIM(DWORK(L(LU)),DWORK(L(LDD)),DWORK(L(LQS)),
     *               KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,KWORK(L(LBD)),
     *               CQ,DWORK(L(LW)+2*NDOF),DWORK(L(LW)+3*NDOF),
     *               DWORK(L(LW)+4*NDOF),ALAVG,NLIM)
         RF=0D0
         RR=0D0
         DO 25 I=1,NDOF
          R=DWORK(L(LRHS)+I-1)-DWORK(L(LW)+NDOF+I-1)
          IF (IS.EQ.3) R=R+DWORK(L(LW)+2*NDOF+I-1)
          RR=RR+R*R
          RF=RF+DWORK(L(LRHS)+I-1)**2
25       CONTINUE
         RES=SQRT(RR/MAX(RF,1D-300))
         CALL ZTIME(T1)
C *** coefficient range
         CMN= 1D30
         CMX=-1D30
         DO 30 I=1,NDOF
          CMN=MIN(CMN,DWORK(L(LU)+I-1))
          CMX=MAX(CMX,DWORK(L(LU)+I-1))
30       CONTINUE
         DO 60 IU=1,2
          IF (IU.EQ.1) MU=MTERM
          IF (IU.EQ.2) MU=MFILE
          IF (ITST.EQ.1) THEN
           CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *                KWORK(L(LDFG)),DWORK(L(LU)),BNUEX,BNUEXX,BNUEXY,
     *                EAL2(IS,ILEV,IE),EAH1(IS,ILEV,IE))
           RL2=0D0
           RH1=0D0
           IF (ILEV.GT.NLMIN) THEN
            RL2=LOG(EAL2(IS,ILEV-1,IE)/EAL2(IS,ILEV,IE))/LOG(2D0)
            RH1=LOG(EAH1(IS,ILEV-1,IE)/EAH1(IS,ILEV,IE))/LOG(2D0)
           ENDIF
           WRITE(MU,1020) EPSD,CSCH(IS),EAL2(IS,ILEV,IE),RL2,
     *      EAH1(IS,ILEV,IE),RH1,RES,NIT,ALAVG,NLIM,T1-T0
          ELSE
           CALL BNRNG(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *                KWORK(L(LDFG)),DWORK(L(LU)),UMN,UMX,WSM)
           WRITE(MU,1030) EPSD,CSCH(IS),CMN,CMX,UMN,UMX,WSM,
     *      RES,NIT,ALAVG,NLIM,T1-T0
          ENDIF
60       CONTINUE
         IF (ILEV.EQ.IGMV) THEN
          WRITE(CGMV,'(A,I1,A,I1,A,A3,A,I1,A)') '#gmv/t',ITST,'_e',
     *     NINT(-LOG10(EPSD)),'_',CSCH(IS)(1:3),'_L',ILEV,'.gmv'
          CALL BNGMVW(K,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *                KWORK(L(LDFG)),DWORK(L(LU)),8,ITST,CGMV)
         ENDIF
70      CONTINUE
C
C *** jump stabilization (eps A + K + gamma J) u = f, direct solve
        DO 76 IGM=1,NGAMR
         GAMJ=GAMV(IGM)
         CALL ZTIME(T0)
         DO 71 IA=1,NAJ
71       DWORK(L(LGJ)+IA-1)=EPSD*DWORK(L(LAJ)+IA-1)
     *     +DWORK(L(LKJ)+IA-1)+GAMJ*DWORK(L(LJJ)+IA-1)
         CALL BNDIRR(DWORK(L(LGJ)),KWORK(L(LLDJ)),NDOF,KWORK(L(LBD)))
         CALL BNBFAC(DWORK(L(LGJ)),KWORK(L(LLDJ)),KWORK(L(LCOLJ)),
     *               NDOF,KWORK(L(LKIN)),KBWJ,DWORK(L(LAB)),LDABJ,
     *               KWORK(L(LIPV)),INFO)
         IF (INFO.NE.0) WRITE(MTERM,*) 'BNAFCT: DGBTRF jump INFO=',INFO
         CALL BNBSOL(DWORK(L(LAB)),LDABJ,NDOF,KBWJ,KWORK(L(LIPV)),
     *               KWORK(L(LKIN)),DWORK(L(LRHS)),DWORK(L(LU)),
     *               DWORK(L(LW)))
         CALL ZTIME(T1)
         CALL BNAX(DWORK(L(LGJ)),KWORK(L(LLDJ)),KWORK(L(LCOLJ)),NDOF,
     *             DWORK(L(LU)),DWORK(L(LW)+NDOF))
         RF=0D0
         RR=0D0
         DO 72 I=1,NDOF
          R=DWORK(L(LRHS)+I-1)-DWORK(L(LW)+NDOF+I-1)
          RR=RR+R*R
          RF=RF+DWORK(L(LRHS)+I-1)**2
72       CONTINUE
         RES=SQRT(RR/MAX(RF,1D-300))
         CMN= 1D30
         CMX=-1D30
         DO 73 I=1,NDOF
          CMN=MIN(CMN,DWORK(L(LU)+I-1))
          CMX=MAX(CMX,DWORK(L(LU)+I-1))
73       CONTINUE
         DO 75 IU=1,2
          IF (IU.EQ.1) MU=MTERM
          IF (IU.EQ.2) MU=MFILE
          IF (ITST.EQ.1) THEN
           CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *                KWORK(L(LDFG)),DWORK(L(LU)),BNUEX,BNUEXX,BNUEXY,
     *                EJL2(IGM,ILEV,IE),EJH1(IGM,ILEV,IE))
           RL2=0D0
           RH1=0D0
           IF (ILEV.GT.NLMIN) THEN
            RL2=LOG(EJL2(IGM,ILEV-1,IE)/EJL2(IGM,ILEV,IE))/LOG(2D0)
            RH1=LOG(EJH1(IGM,ILEV-1,IE)/EJH1(IGM,ILEV,IE))/LOG(2D0)
           ENDIF
           WRITE(MU,1022) EPSD,GAMJ,EJL2(IGM,ILEV,IE),RL2,
     *      EJH1(IGM,ILEV,IE),RH1,RES,T1-T0
          ELSE
           CALL BNRNG(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *                KWORK(L(LDFG)),DWORK(L(LU)),UMN,UMX,WSM)
           WRITE(MU,1032) EPSD,GAMJ,CMN,CMX,UMN,UMX,WSM,RES,T1-T0
          ENDIF
75       CONTINUE
         IF (ILEV.EQ.IGMV) THEN
          WRITE(CGMV,'(A,I1,A,I1,A,I1,A,I1,A)') '#gmv/t',ITST,'_e',
     *     NINT(-LOG10(EPSD)),'_JMP',IGM,'_L',ILEV,'.gmv'
          CALL BNGMVW(K,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *                KWORK(L(LDFG)),DWORK(L(LU)),8,ITST,CGMV)
         ENDIF
76      CONTINUE
80      CONTINUE
90     CONTINUE
C
       CALL ZDISP(0,LAAW,'DAAWAF')
       CALL ZDISP(0,LB2J,'DB2J  ')
       CALL ZDISP(0,LB1J,'DB1J  ')
       CALL ZDISP(0,LGJ ,'DGJ   ')
       CALL ZDISP(0,LJJ ,'DJJ   ')
       CALL ZDISP(0,LKJ ,'DKJ   ')
       CALL ZDISP(0,LAJ ,'DAJ   ')
       CALL ZDISP(0,LMJ ,'DMJ   ')
       CALL ZDISP(0,LCOLJ,'KCOLBN')
       CALL ZDISP(0,LLDJ ,'KLDBN ')
       CALL ZDISP(0,LPAT,'KPATJ ')
       CALL ZDISP(0,LEDG,'KEDGJ ')
       CALL ZDISP(0,LAB ,'DABAF ')
       CALL ZDISP(0,LIPV,'IPIVAF')
       CALL ZDISP(0,LKIN,'KINAF ')
       CALL ZDISP(0,LKP ,'KPAF  ')
       CALL ZDISP(0,LW  ,'DWAF  ')
       CALL ZDISP(0,LU  ,'DUAF  ')
       CALL ZDISP(0,LRHS,'DRHSAF')
       CALL ZDISP(0,LCP ,'DCPAF ')
       CALL ZDISP(0,LQS ,'DQSAF ')
       CALL ZDISP(0,LKT ,'KTAF  ')
       CALL ZDISP(0,LDD ,'DDAF  ')
       CALL ZDISP(0,LLL ,'DLAF  ')
       CALL ZDISP(0,LG  ,'DGAF  ')
       CALL ZDISP(0,LK  ,'DKAF  ')
       CALL ZDISP(0,LB2 ,'DB2AF ')
       CALL ZDISP(0,LB1 ,'DB1AF ')
       CALL ZDISP(0,LA  ,'DAAF  ')
       CALL ZDISP(0,LM  ,'DMAF  ')
       CALL ZDISP(0,LCOLB,'KCOLBN')
       CALL ZDISP(0,LLDB ,'KLDBN ')
       CALL ZDISP(0,LCOLA,'KCOLBN')
       CALL ZDISP(0,LLDA ,'KLDBN ')
       CALL ZDISP(0,LYD ,'DYD   ')
       CALL ZDISP(0,LXD ,'DXD   ')
       CALL ZDISP(0,LNOB,'KNOB  ')
       CALL ZDISP(0,LBDB,'KBDB  ')
       CALL ZDISP(0,LBDA,'KBDA  ')
       CALL ZDISP(0,LDFP,'KDFP  ')
       CALL ZDISP(0,LDFG,'KDFG  ')
100   CONTINUE
C
1000  FORMAT(/' AFC TEST  Bernstein Q',I1,
     * '  -eps Lap u + b.grad u = f, b = potential flow'/
     * ' CQ=',F5.2,' OMEGA=',F5.2,' EPSNL=',D8.1,' NITMX=',I5,
     * ' Anderson MAA=',I3,' beta=',F5.2/
     * ' A: smooth u=sin3x cos2y+xy, Dirichlet everywhere'/
     * ' B: inflow step u=1 (y>0.3), Neumann elsewhere, u in [0,1]'/
     * ' ',100('-'))
1001  FORMAT(' jump gamma values:',8D10.2)
1010  FORMAT(/' LEV',I2,' NEL=',I7,' NDOF=',I8,' NA=',I9,
     * ' inflow dofs=',I5,' max|K*1|=',D9.2,' bandwidth=',I6)
1011  FORMAT('   jump: interior edges=',I8,' NA(ext)=',I9,
     * ' bandwidth=',I6,' max|J*x|,|J*y| / max|J| =',D9.2)
1012  FORMAT('   L2 projection of u (best L2 approximation): L2=',D10.3,
     * ' H1=',D10.3)
1022  FORMAT(' A eps=',D7.1,' JMP g=',D7.1,' L2=',D10.3,' (',F5.2,
     * ') H1=',D10.3,' (',F5.2,') res=',D8.1,' t=',F7.2)
1032  FORMAT(' B eps=',D7.1,' JMP g=',D7.1,' coef[',F8.4,',',F8.4,
     * '] u[',F8.4,',',F8.4,'] W=',D10.3,' res=',D8.1,' t=',F7.2)
1020  FORMAT(' A eps=',D7.1,1X,A4,' L2=',D10.3,' (',F5.2,') H1=',D10.3,
     * ' (',F5.2,') res=',D8.1,' nl=',I5,' alpha=',F6.3,' nlim=',I7,
     * ' t=',F7.2)
1030  FORMAT(' B eps=',D7.1,1X,A4,' coef[',F8.4,',',F8.4,'] u[',F8.4,
     * ',',F8.4,'] W=',D10.3,' res=',D8.1,' nl=',I5,' alpha=',F6.3,
     * ' nlim=',I7,' t=',F7.2)
      END
C
************************************************************************
*   BNPERM: ordering of the dofs by increasing x (Greville points) for
*           the banded LU; KP(new)=old, KINV(old)=new; KBW = half
*           bandwidth of the permuted matrix structure
************************************************************************
      SUBROUTINE BNPERM(DXD,N,KLD,KCOL,KP,KINV,KBW)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DXD(*),KLD(*),KCOL(*),KP(*),KINV(*)
      DO 5 I=1,N
5     KP(I)=I
C *** shell sort of KP by DXD
      IGAP=N/2
10    IF (IGAP.LT.1) GOTO 40
      DO 30 I=IGAP+1,N
       ITMP=KP(I)
       X=DXD(ITMP)
       J=I
20     IF (J.GT.IGAP) THEN
        IF (DXD(KP(J-IGAP)).GT.X) THEN
         KP(J)=KP(J-IGAP)
         J=J-IGAP
         GOTO 20
        ENDIF
       ENDIF
       KP(J)=ITMP
30    CONTINUE
      IGAP=IGAP/2
      GOTO 10
40    DO 45 I=1,N
45    KINV(KP(I))=I
      KBW=0
      DO 60 IR=1,N
      DO 60 IA=KLD(IR),KLD(IR+1)-1
       KBW=MAX(KBW,ABS(KINV(IR)-KINV(KCOL(IA))))
60    CONTINUE
      END
C
************************************************************************
*   BNBFAC: banded LU (LAPACK DGBTRF) of the CSR matrix DA in the
*           ordering KINV, KL=KU=KBW, LDAB=3*KBW+1
************************************************************************
      SUBROUTINE BNBFAC(DA,KLD,KCOL,N,KINV,KBW,AB,LDAB,IPIV,INFO)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),KLD(*),KCOL(*),KINV(*),AB(LDAB,*),IPIV(*)
      DO 10 J=1,N
      DO 10 I=1,LDAB
10    AB(I,J)=0D0
      DO 20 IR=1,N
      DO 20 IA=KLD(IR),KLD(IR+1)-1
       IRN=KINV(IR)
       JCN=KINV(KCOL(IA))
       AB(2*KBW+1+IRN-JCN,JCN)=DA(IA)
20    CONTINUE
      CALL DGBTRF(N,N,KBW,KBW,AB,LDAB,IPIV,INFO)
      END
C
************************************************************************
*   BNBSOL: DX = A^(-1) DB with the factors of BNBFAC, DY work (N)
************************************************************************
      SUBROUTINE BNBSOL(AB,LDAB,N,KBW,IPIV,KINV,DB,DX,DY)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION AB(LDAB,*),IPIV(*),KINV(*),DB(*),DX(*),DY(*)
      DO 10 I=1,N
10    DY(KINV(I))=DB(I)
      CALL DGBTRS('N',N,KBW,KBW,1,AB,LDAB,IPIV,DY,N,INFO)
      IF (INFO.NE.0) STOP 'BNBSOL: DGBTRS failed'
      DO 20 I=1,N
20    DX(I)=DY(KINV(I))
      END
C
************************************************************************
*   Edge-oriented jump stabilization for Bernstein Qk (scalar test)
*     J(u,v) = sum_E h_E^2 int_E [grad u].[grad v] ds  (interior edges)
*   assembled on an extended structure (all dofs of the two elements
*   at an edge are coupled), scheme JMP: (eps A + K + gamma J) u = f.
************************************************************************
*   BNJEDG: interior edges KEDG(1..5,IE) = IEL1, IME1, IEL2, IME2,
*           ISGN (-1: the edge runs in opposite direction in IEL2)
************************************************************************
      SUBROUTINE BNJEDG(KADJ,KVERT,NEL,KEDG,NEDG)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION KADJ(NNVE,*),KVERT(NNVE,*),KEDG(5,*)
      NEDG=0
      DO 20 IEL=1,NEL
      DO 20 IME=1,4
       JEL=KADJ(IME,IEL)
       IF (JEL.LE.IEL) GOTO 20
       JME=0
       DO 10 J=1,4
10     IF (KADJ(J,JEL).EQ.IEL) JME=J
       IF (JME.EQ.0) STOP 'BNJEDG: neighbour relation broken'
       NEDG=NEDG+1
       KEDG(1,NEDG)=IEL
       KEDG(2,NEDG)=IME
       KEDG(3,NEDG)=JEL
       KEDG(4,NEDG)=JME
       IF (KVERT(IME,IEL).EQ.KVERT(MOD(JME,4)+1,JEL)) THEN
        KEDG(5,NEDG)=-1
       ELSEIF (KVERT(IME,IEL).EQ.KVERT(JME,JEL)) THEN
        KEDG(5,NEDG)=1
       ELSE
        STOP 'BNJEDG: edge vertices do not match'
       ENDIF
20    CONTINUE
      END
C
************************************************************************
*   BNJPAT: patches for the extended structure: element dofs (twice)
*           and, per interior edge, the dofs of both elements
************************************************************************
      SUBROUTINE BNJPAT(NLDF,KDFG,NEL,KEDG,NEDG,KPAT)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KDFG(NLDF,*),KEDG(5,*),KPAT(2*NLDF,*)
      DO 10 IEL=1,NEL
      DO 10 I=1,NLDF
       KPAT(I,IEL)=KDFG(I,IEL)
       KPAT(NLDF+I,IEL)=KDFG(I,IEL)
10    CONTINUE
      DO 20 IE=1,NEDG
      DO 20 I=1,NLDF
       KPAT(I,NEL+IE)=KDFG(I,KEDG(1,IE))
       KPAT(NLDF+I,NEL+IE)=KDFG(I,KEDG(3,IE))
20    CONTINUE
      END
C
************************************************************************
*   BNGRDP: physical gradients of the Qk Bernstein basis at the
*           reference point (XI,ETA) of the element (DXV,DYV)
************************************************************************
      SUBROUTINE BNGRDP(K,XI,ETA,DXV,DYV,DPX,DPY)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8)
      DIMENSION DXV(4),DYV(4),DPX(*),DPY(*)
      DIMENSION DV1(0:NNBK),D1(0:NNBK),DV2(0:NNBK),D2(0:NNBK)
      CALL BNMAP(XI,ETA,DXV,DYV,X,Y,DJ11,DJ12,DJ21,DJ22,DET)
      CALL BNB1D(K,XI ,DV1,D1)
      CALL BNB1D(K,ETA,DV2,D2)
      DO 20 J=0,K
      DO 20 I=0,K
       L=I+J*(K+1)+1
       DPXI=D1(I)*DV2(J)
       DPET=DV1(I)*D2(J)
       DPX(L)=( DJ22*DPXI-DJ21*DPET)/DET
       DPY(L)=(-DJ12*DPXI+DJ11*DPET)/DET
20    CONTINUE
      END
C
************************************************************************
*   BNJASM: DJ += sum_E h_E^2 int_E [grad phi_j].[grad phi_i] ds,
*           Gauss with K+1 points per edge, h_E = edge length
************************************************************************
      SUBROUTINE BNJASM(K,NEL,DCORVG,KVERT,KDFG,KEDG,NEDG,KLD,KCOL,DJ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNG=NNBK+2,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG((K+1)**2,*),KEDG(5,*)
      DIMENSION KLD(*),KCOL(*),DJ(*)
      DIMENSION XG(NNG),WG(NNG),DXV1(4),DYV1(4),DXV2(4),DYV2(4)
      DIMENSION G1X(NNLD),G1Y(NNLD),G2X(NNLD),G2Y(NNLD)
      DIMENSION XR(4),YR(4)
      INTEGER*8 BNFIND
      DATA XR /-1D0, 1D0, 1D0,-1D0/
      DATA YR /-1D0,-1D0, 1D0, 1D0/
C
      NLDF=(K+1)**2
      NGE=K+1
      CALL BNGAUS(NGE,XG,WG)
      DO 100 IE=1,NEDG
       IEL1=KEDG(1,IE)
       IME1=KEDG(2,IE)
       IEL2=KEDG(3,IE)
       IME2=KEDG(4,IE)
       ISGN=KEDG(5,IE)
       DO 5 IV=1,4
        DXV1(IV)=DCORVG(1,KVERT(IV,IEL1))
        DYV1(IV)=DCORVG(2,KVERT(IV,IEL1))
        DXV2(IV)=DCORVG(1,KVERT(IV,IEL2))
        DYV2(IV)=DCORVG(2,KVERT(IV,IEL2))
5      CONTINUE
       IA1=IME1
       IB1=MOD(IME1,4)+1
       IA2=IME2
       IB2=MOD(IME2,4)+1
       HE=SQRT((DXV1(IB1)-DXV1(IA1))**2+(DYV1(IB1)-DYV1(IA1))**2)
       DO 50 IG=1,NGE
        S1=XG(IG)
        S2=DBLE(ISGN)*S1
        XI1=0.5D0*((1D0-S1)*XR(IA1)+(1D0+S1)*XR(IB1))
        ET1=0.5D0*((1D0-S1)*YR(IA1)+(1D0+S1)*YR(IB1))
        XI2=0.5D0*((1D0-S2)*XR(IA2)+(1D0+S2)*XR(IB2))
        ET2=0.5D0*((1D0-S2)*YR(IA2)+(1D0+S2)*YR(IB2))
C ***   curved elements (IISO): their map in BNMAP; the interior edges
C       stay straight (blend = 0 at the vertices), HE and ds unchanged
        CALL BNGSET(IEL1)
        CALL BNGRDP(K,XI1,ET1,DXV1,DYV1,G1X,G1Y)
        CALL BNGSET(IEL2)
        CALL BNGRDP(K,XI2,ET2,DXV2,DYV2,G2X,G2Y)
        W=HE*HE*WG(IG)*0.5D0*HE
C *** [grad phi] = grad phi|_1 (dofs of IEL1) - grad phi|_2 (IEL2)
        DO 40 I=1,NLDF
         IR1=KDFG(I,IEL1)
         IR2=KDFG(I,IEL2)
         DO 30 J=1,NLDF
          JC1=KDFG(J,IEL1)
          JC2=KDFG(J,IEL2)
          IA=BNFIND(KLD,KCOL,IR1,JC1)
          DJ(IA)=DJ(IA)+W*(G1X(I)*G1X(J)+G1Y(I)*G1Y(J))
          IA=BNFIND(KLD,KCOL,IR1,JC2)
          DJ(IA)=DJ(IA)-W*(G1X(I)*G2X(J)+G1Y(I)*G2Y(J))
          IA=BNFIND(KLD,KCOL,IR2,JC1)
          DJ(IA)=DJ(IA)-W*(G2X(I)*G1X(J)+G2Y(I)*G1Y(J))
          IA=BNFIND(KLD,KCOL,IR2,JC2)
          DJ(IA)=DJ(IA)+W*(G2X(I)*G2X(J)+G2Y(I)*G2Y(J))
30       CONTINUE
40      CONTINUE
50     CONTINUE
100   CONTINUE
      CALL BNGSET(0)
      END
C
************************************************************************
*   BNGMVW: GMV file of a Bernstein Qk function: every element is cut
*           into NS x NS quads, u_h evaluated at the (NS+1)^2 points
*           (nodes are not shared, the output is exact per element);
*           variables uh, and for ITST=1 uex and err = uh - uex
************************************************************************
      SUBROUTINE BNGMVW(K,NEL,DCORVG,KVERT,KDFG,DC,NS,ITST,CFILE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNVE=4)
      CHARACTER*(*) CFILE
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG((K+1)**2,*),DC(*)
      DIMENSION DXV(4),DYV(4),DV1(0:NNBK),D1(0:NNBK)
      DIMENSION DV2(0:NNBK),D2(0:NNBK)
      DOUBLE PRECISION BNUEX
      EXTERNAL BNUEX
      MG=91
      NP1=NS+1
      NNOD=NEL*NP1*NP1
      OPEN(MG,FILE=CFILE,STATUS='UNKNOWN')
      WRITE(MG,'(A)') 'gmvinput ascii'
      WRITE(MG,'(A,I10)') 'nodes ',NNOD
      DO 30 ICO=1,3
      DO 30 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       DO 20 JS=0,NS
       DO 20 IS=0,NS
        XI =-1D0+2D0*DBLE(IS)/DBLE(NS)
        ETA=-1D0+2D0*DBLE(JS)/DBLE(NS)
        CALL BNMAP(XI,ETA,DXV,DYV,X,Y,D11,D12,D21,D22,DET)
        IF (ICO.EQ.1) WRITE(MG,'(E16.8)') X
        IF (ICO.EQ.2) WRITE(MG,'(E16.8)') Y
        IF (ICO.EQ.3) WRITE(MG,'(E16.8)') 0D0
20     CONTINUE
30    CONTINUE
      WRITE(MG,'(A,I10)') 'cells ',NEL*NS*NS
      DO 40 IEL=1,NEL
      DO 40 JS=0,NS-1
      DO 40 IS=0,NS-1
       N1=(IEL-1)*NP1*NP1+JS*NP1+IS+1
       WRITE(MG,'(A)') 'quad 4'
       WRITE(MG,'(4I10)') N1,N1+1,N1+1+NP1,N1+NP1
40    CONTINUE
      WRITE(MG,'(A)') 'variable'
      NVAR=1
      IF (ITST.EQ.1) NVAR=3
      DO 70 IVAR=1,NVAR
       IF (IVAR.EQ.1) WRITE(MG,'(A)') 'uh 1'
       IF (IVAR.EQ.2) WRITE(MG,'(A)') 'uex 1'
       IF (IVAR.EQ.3) WRITE(MG,'(A)') 'err 1'
       DO 60 IEL=1,NEL
        DO 45 IV=1,4
         DXV(IV)=DCORVG(1,KVERT(IV,IEL))
         DYV(IV)=DCORVG(2,KVERT(IV,IEL))
45      CONTINUE
        DO 55 JS=0,NS
        DO 55 IS=0,NS
         XI =-1D0+2D0*DBLE(IS)/DBLE(NS)
         ETA=-1D0+2D0*DBLE(JS)/DBLE(NS)
         CALL BNMAP(XI,ETA,DXV,DYV,X,Y,D11,D12,D21,D22,DET)
         CALL BNB1D(K,XI ,DV1,D1)
         CALL BNB1D(K,ETA,DV2,D2)
         U=0D0
         DO 50 J=0,K
         DO 50 I=0,K
50       U=U+DC(KDFG(I+J*(K+1)+1,IEL))*DV1(I)*DV2(J)
         IF (IVAR.EQ.2) U=BNUEX(X,Y)
         IF (IVAR.EQ.3) U=U-BNUEX(X,Y)
         IF (ABS(U).LT.1D-99) U=0D0
         WRITE(MG,'(E16.8)') U
55      CONTINUE
60     CONTINUE
70    CONTINUE
      WRITE(MG,'(A)') 'endvars'
      WRITE(MG,'(A)') 'endgmv'
      CLOSE(MG)
      END
C
************************************************************************
*   BNCHKX: check of BNDFXY: max over the elements and 3x3 points of
*     |grad(sum_l DXD(l) phi_l) - (1,0)| and the same for y -> (0,1);
*     separately for curved (IISO) and bilinear elements
************************************************************************
      SUBROUTINE BNCHKX(K,NEL,DCORVG,KVERT,KDFG,DXD,DYD,ERRS,ERRC)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG((K+1)**2,*)
      DIMENSION DXD(*),DYD(*),DXV(4),DYV(4),GX(NNLD),GY(NNLD)
      COMMON /BNCUR/  ICUR,MGC,CGX(16),CGY(16)
      ERRS=0D0
      ERRC=0D0
      DO 20 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNGSET(IEL)
       DO 10 IA=-1,1
       DO 10 IB=-1,1
        CALL BNGRDP(K,0.9D0*IA,0.9D0*IB,DXV,DYV,GX,GY)
        XX=0D0
        XY=0D0
        YX=0D0
        YY=0D0
        DO 8 L=1,(K+1)**2
         XX=XX+DXD(KDFG(L,IEL))*GX(L)
         XY=XY+DXD(KDFG(L,IEL))*GY(L)
         YX=YX+DYD(KDFG(L,IEL))*GX(L)
         YY=YY+DYD(KDFG(L,IEL))*GY(L)
8       CONTINUE
        E=MAX(ABS(XX-1D0),ABS(XY),ABS(YX),ABS(YY-1D0))
        IF (ICUR.EQ.1) THEN
         ERRC=MAX(ERRC,E)
        ELSE
         ERRS=MAX(ERRS,E)
        ENDIF
10     CONTINUE
20    CONTINUE
      CALL BNGSET(0)
      END
