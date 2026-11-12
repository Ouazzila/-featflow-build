************************************************************************
* bern.f  -  Bernstein Qk / P(k-1)disc mixed finite elements
*
* Stage 1: element kernel (Bernstein bases, Gauss cubature, bilinear
* map), global dof numbering, assembly of M, A (Laplacian), B1, B2,
* and the test driver BNTEST (switch IBERN=k in pp2d.dat).
*
* Velocity  Qk  : phi_ij(xi,eta) = B_i^k(xi) B_j^k(eta), i,j=0..k,
*                 B_i^k(x) = C(k,i) t^i (1-t)^(k-i), t=(x+1)/2,
*                 local index l = i + j*(k+1) + 1
* Pressure  P(k-1)disc : Bernstein basis of P_(k-1) w.r.t. the
*                 reference triangle (-1,-1),(3,-1),(-1,3),
*                 lam2=(xi+1)/4, lam3=(eta+1)/4, lam1=1-lam2-lam3
* Global velocity dofs: NVT vertices, (k-1) per edge, (k-1)^2 per
* element (see NOTES.md). Matrices are CSR, square ones with the
* diagonal entry first in every row (FEAT storage technique 7).
* Own assembly (no /ELEM/, NDFGL): FEAT limits elements to 21 basis
* functions and knows no dofs beyond edge midpoints.
************************************************************************
C
************************************************************************
*   BNB1D: 1D Bernstein polynomials of degree K on [-1,1] in X,
*          values DB(0:K) and derivatives d/dx DDB(0:K)
*          (de Casteljau recursion, d/dt B_i^K = K (B_(i-1)^(K-1)
*           - B_i^(K-1)), dt/dx = 1/2)
************************************************************************
      SUBROUTINE BNB1D(K,X,DB,DDB)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8)
      DIMENSION DB(0:*),DDB(0:*),DW(0:NNBK)
C
      T=0.5D0*(X+1D0)
      S=1D0-T
      DW(0)=1D0
      DDB(0)=0D0
      DO 20 N=1,K
       IF (N.EQ.K) THEN
        DO 10 I=0,K
         DM1=0D0
         DI =0D0
         IF (I.GE.1)   DM1=DW(I-1)
         IF (I.LE.K-1) DI =DW(I)
         DDB(I)=0.5D0*DBLE(K)*(DM1-DI)
10      CONTINUE
       ENDIF
       DW(N)=T*DW(N-1)
       DO 15 I=N-1,1,-1
        DW(I)=S*DW(I)+T*DW(I-1)
15     CONTINUE
       DW(0)=S*DW(0)
20    CONTINUE
      DO 30 I=0,K
30    DB(I)=DW(I)
C
      END
C
************************************************************************
*   BNPB: Bernstein basis of P_N (N=k-1) on the reference square,
*         barycentric coordinates of the triangle (-1,-1),(3,-1),(-1,3)
*         DP(l), l=1..(N+1)(N+2)/2, ordering: a2=0..N, a3=0..N-a2
************************************************************************
      SUBROUTINE BNPB(N,XI,ETA,DP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DP(*)
      DOUBLE PRECISION BNFAC
C
      DL2=0.25D0*(XI +1D0)
      DL3=0.25D0*(ETA+1D0)
      DL1=1D0-DL2-DL3
      L=0
      DO 10 IA2=0,N
      DO 10 IA3=0,N-IA2
       IA1=N-IA2-IA3
       L=L+1
       DP(L)=BNFAC(N)/(BNFAC(IA1)*BNFAC(IA2)*BNFAC(IA3))
     *       *DL1**IA1*DL2**IA2*DL3**IA3
10    CONTINUE
C
      END
C
      DOUBLE PRECISION FUNCTION BNFAC(N)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      BNFAC=1D0
      DO 10 I=2,N
10    BNFAC=BNFAC*DBLE(I)
      END
C
************************************************************************
*   BNGAUS: NG-point Gauss-Legendre rule on [-1,1] (Newton iteration)
************************************************************************
      SUBROUTINE BNGAUS(NG,XG,WG)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION XG(*),WG(*)
C
      PI=4D0*ATAN(1D0)
      DO 30 I=1,NG
       Z=COS(PI*(DBLE(I)-0.25D0)/(DBLE(NG)+0.5D0))
       DO 20 IT=1,100
        P1=1D0
        P2=0D0
        DO 10 J=1,NG
         P3=P2
         P2=P1
         P1=(DBLE(2*J-1)*Z*P2-DBLE(J-1)*P3)/DBLE(J)
10      CONTINUE
        PP=DBLE(NG)*(Z*P1-P2)/(Z*Z-1D0)
        Z1=Z
        Z=Z1-P1/PP
        IF (ABS(Z-Z1).LE.1D-15) GOTO 25
20     CONTINUE
25     XG(I)=-Z
       WG(I)=2D0/((1D0-Z*Z)*PP*PP)
30    CONTINUE
C
      END
C
************************************************************************
*   BNEVQK: values/physical gradients of the Qk Bernstein basis and
*           values of the P(k-1) basis in all NG*NG Gauss points of
*           the element with vertices (DXV,DYV) (bilinear map).
*           DWT = Gauss weight * det(J), DXQ/DYQ = physical points.
*           Reference vertices (-1,-1),(1,-1),(1,1),(-1,1) = KVERT(1..4)
************************************************************************
      SUBROUTINE BNEVQK(K,NG,XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,
     *                  DWT,DXQ,DYQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2)
      DIMENSION XG(*),WG(*),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,*),DPX(NNLD,*),DPY(NNLD,*),DPSI(NNPL,*)
      DIMENSION DWT(*),DXQ(*),DYQ(*)
      DIMENSION DV1(0:NNBK),D1(0:NNBK),DV2(0:NNBK),D2(0:NNBK)
C
      IG=0
      DO 40 JQ=1,NG
      DO 40 IQ=1,NG
       IG=IG+1
       XI =XG(IQ)
       ETA=XG(JQ)
       CALL BNMAP(XI,ETA,DXV,DYV,X,Y,DJ11,DJ12,DJ21,DJ22,DET)
       DXQ(IG)=X
       DYQ(IG)=Y
       DWT(IG)=WG(IQ)*WG(JQ)*DET
C *** tensor-product Bernstein basis, grad = J^(-T) * reference grad
       CALL BNB1D(K,XI ,DV1,D1)
       CALL BNB1D(K,ETA,DV2,D2)
       DO 20 J=0,K
       DO 20 I=0,K
        L=I+J*(K+1)+1
        DPXI=D1(I)*DV2(J)
        DPET=DV1(I)*D2(J)
        DPHI(L,IG)=DV1(I)*DV2(J)
        DPX(L,IG)=( DJ22*DPXI-DJ21*DPET)/DET
        DPY(L,IG)=(-DJ12*DPXI+DJ11*DPET)/DET
20     CONTINUE
       CALL BNPB(K-1,XI,ETA,DPSI(1,IG))
40    CONTINUE
C
      END
C
************************************************************************
*   BNMAP: bilinear map of the reference square onto the element with
*          vertices (DXV,DYV), point (X,Y) and Jacobian
*          DJ11=dx/dxi DJ12=dx/deta DJ21=dy/dxi DJ22=dy/deta, DET=det J
************************************************************************
      SUBROUTINE BNMAP(XI,ETA,DXV,DYV,X,Y,DJ11,DJ12,DJ21,DJ22,DET)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DXV(4),DYV(4),DN(4),DNX(4),DNE(4)
C *** curved element (isoparametric boundary, bniso.f, set by BNGSET)
      COMMON /BNCUR/  ICUR,MGC,CGX(16),CGY(16)
C
      IF (ICUR.EQ.1) THEN
       CALL BNCMAP(MGC,CGX,CGY,XI,ETA,X,Y,DJ11,DJ12,DJ21,DJ22,DET)
       RETURN
      ENDIF
      DN(1)=0.25D0*(1D0-XI)*(1D0-ETA)
      DN(2)=0.25D0*(1D0+XI)*(1D0-ETA)
      DN(3)=0.25D0*(1D0+XI)*(1D0+ETA)
      DN(4)=0.25D0*(1D0-XI)*(1D0+ETA)
      DNX(1)=-0.25D0*(1D0-ETA)
      DNX(2)= 0.25D0*(1D0-ETA)
      DNX(3)= 0.25D0*(1D0+ETA)
      DNX(4)=-0.25D0*(1D0+ETA)
      DNE(1)=-0.25D0*(1D0-XI)
      DNE(2)=-0.25D0*(1D0+XI)
      DNE(3)= 0.25D0*(1D0+XI)
      DNE(4)= 0.25D0*(1D0-XI)
      X=0D0
      Y=0D0
      DJ11=0D0
      DJ12=0D0
      DJ21=0D0
      DJ22=0D0
      DO 10 IV=1,4
       X=X+DN(IV)*DXV(IV)
       Y=Y+DN(IV)*DYV(IV)
       DJ11=DJ11+DNX(IV)*DXV(IV)
       DJ12=DJ12+DNE(IV)*DXV(IV)
       DJ21=DJ21+DNX(IV)*DYV(IV)
       DJ22=DJ22+DNE(IV)*DYV(IV)
10    CONTINUE
      DET=DJ11*DJ22-DJ12*DJ21
C
      END
C
************************************************************************
*   BNEVAL: element evaluation for the active family (COMMON /BNFAM/)
*     IFAM = 0  : Bernstein Qk / P(k-1)disc, k = KORD     (BNEVQK)
*     IFAM = 30 : FEAT E030 (Q1~, edge means) / Q0         (BQEVAL)
*     IFAM = 31 : FEAT E031 (Q1~, edge midpoints) / Q0     (BQEVAL)
*     IFAM = 130: EM30 nonparametric Q1~ (edge means) / Q0 (BQEVAL)
*     IFAM = 131: EM31 nonparametric Q1~ (midpoints) / Q0  (BQEVAL)
*   NG*NG Gauss points, arrays as in BNEVQK
************************************************************************
      SUBROUTINE BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION XG(*),WG(*),DXV(4),DYV(4)
      DIMENSION DPHI(*),DPX(*),DPY(*),DPSI(*),DWT(*),DXQ(*),DYQ(*)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
C
      IF (IFAM.EQ.0) THEN
       CALL BNEVQK(KORD,NG,XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,
     *             DWT,DXQ,DYQ)
      ELSE
       CALL BQEVAL(IFAM,NG,XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,
     *             DWT,DXQ,DYQ)
      ENDIF
C
      END
C
************************************************************************
*   BQEVAL: the original pp2d pair Q1~/Q0 evaluated with the FEAT
*           library elements E030/E031 (parametric, as IELT=1/0 in
*           pp2d): /ELEM/ gets DJAC, DETJ, BDER, the element returns
*           DBAS(i,1..3) = value, d/dx, d/dy; local dof i = edge i.
*           IFAMQ = 130/131: nonparametric EM30/EM31 (util.f, as
*           IELT=3/2, pp2d default IELT=3), called as in ABM7
*           (xmrout.f): DX,DY,KVE set, ELE(0,0,-2) once per element,
*           then ELE(x,y,-3) in the physical cubature points.
************************************************************************
      SUBROUTINE BQEVAL(IFAMQ,NG,XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,
     *                  DWT,DXQ,DYQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4)
      DIMENSION XG(*),WG(*),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,*),DPX(NNLD,*),DPY(NNLD,*),DPSI(NNPL,*)
      DIMENSION DWT(*),DXQ(*),DYQ(*)
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE /ELEM/
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      EXTERNAL E030,E031,EM30,EM31
C
      DO 5 I=1,NNDER
5     BDER(I)=I.LE.3
      IF (IFAMQ.GE.100) THEN
C *** nonparametric: basis coefficients from the element geometry
       DO 7 IVE=1,4
        DX(IVE)=DXV(IVE)
        DY(IVE)=DYV(IVE)
7      CONTINUE
       IF (IFAMQ.EQ.130) CALL EM30(0D0,0D0,-2)
       IF (IFAMQ.EQ.131) CALL EM31(0D0,0D0,-2)
      ENDIF
      IG=0
      DO 40 JQ=1,NG
      DO 40 IQ=1,NG
       IG=IG+1
       XI =XG(IQ)
       ETA=XG(JQ)
       CALL BNMAP(XI,ETA,DXV,DYV,X,Y,DJ11,DJ12,DJ21,DJ22,DET)
       DXQ(IG)=X
       DYQ(IG)=Y
       DWT(IG)=WG(IQ)*WG(JQ)*DET
       DJAC(1,1)=DJ11
       DJAC(1,2)=DJ12
       DJAC(2,1)=DJ21
       DJAC(2,2)=DJ22
       DETJ=DET
       IPAR=0
       IF (IFAMQ.EQ.30) THEN
        CALL E030(XI,ETA,IPAR)
       ELSE IF (IFAMQ.EQ.31) THEN
        CALL E031(XI,ETA,IPAR)
       ELSE IF (IFAMQ.EQ.130) THEN
        CALL EM30(X,Y,-3)
       ELSE
        CALL EM31(X,Y,-3)
       ENDIF
       DO 20 L=1,4
        DPHI(L,IG)=DBAS(L,1)
        DPX (L,IG)=DBAS(L,2)
        DPY (L,IG)=DBAS(L,3)
20     CONTINUE
       DPSI(1,IG)=1D0
40    CONTINUE
C
      END
C
************************************************************************
*   BQDFMP: dofs of Q1~/Q0: KDFG(i,IEL) = KMID(i,IEL)-NVT (edge i),
*           boundary marker KBD for edges without neighbour
************************************************************************
      SUBROUTINE BQDFMP(KMID,KADJ,NEL,NVT,NMT,KDFG,KBD)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION KMID(NNVE,*),KADJ(NNVE,*),KDFG(4,*),KBD(*)
C
      DO 5 IMT=1,NMT
5     KBD(IMT)=0
      DO 10 IEL=1,NEL
      DO 10 IME=1,4
       KDFG(IME,IEL)=KMID(IME,IEL)-NVT
       IF (KADJ(IME,IEL).EQ.0) KBD(KMID(IME,IEL)-NVT)=1
10    CONTINUE
C
      END
C
************************************************************************
*   BNDFMP: global velocity dofs KDFG(l,IEL), l = i+j(K+1)+1
*     vertices: KVERT; edge m = KVERT(m)->KVERT(m+1), global edge
*     KMID(m)-NVT, dof s=1..K-1 counted from the lower vertex number;
*     interior dofs after all edge dofs.
************************************************************************
      SUBROUTINE BNDFMP(K,KVERT,KMID,NEL,NVT,NMT,KDFG)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KDFG((K+1)**2,*)
C
      K1=K+1
      DO 30 IEL=1,NEL
C *** vertices (0,0),(K,0),(K,K),(0,K)
       KDFG(1        ,IEL)=KVERT(1,IEL)
       KDFG(K+1      ,IEL)=KVERT(2,IEL)
       KDFG(K+K*K1+1 ,IEL)=KVERT(3,IEL)
       KDFG(K*K1+1   ,IEL)=KVERT(4,IEL)
C *** edges, position IS from KVERT(IME) towards KVERT(IME+1)
       DO 20 IME=1,4
        IMT=KMID(IME,IEL)-NVT
        IV1=KVERT(IME,IEL)
        IV2=KVERT(MOD(IME,4)+1,IEL)
        DO 10 IS=1,K-1
         IF (IV1.LT.IV2) THEN
          IG=NVT+(IMT-1)*(K-1)+IS
         ELSE
          IG=NVT+(IMT-1)*(K-1)+K-IS
         ENDIF
         IF (IME.EQ.1) L=IS+1
         IF (IME.EQ.2) L=K+IS*K1+1
         IF (IME.EQ.3) L=(K-IS)+K*K1+1
         IF (IME.EQ.4) L=(K-IS)*K1+1
         KDFG(L,IEL)=IG
10      CONTINUE
20     CONTINUE
C *** interior
       IOFF=NVT+NMT*(K-1)+(IEL-1)*(K-1)**2
       DO 25 J=1,K-1
       DO 25 I=1,K-1
        KDFG(I+J*K1+1,IEL)=IOFF+(I-1)+(J-1)*(K-1)+1
25     CONTINUE
30    CONTINUE
C
      END
C
************************************************************************
*   BNDFPP: global pressure dofs KDFP(l,IEL) = (IEL-1)*NPL+l
************************************************************************
      SUBROUTINE BNDFPP(NPL,NEL,KDFP)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KDFP(NPL,*)
      DO 10 IEL=1,NEL
      DO 10 L=1,NPL
10    KDFP(L,IEL)=(IEL-1)*NPL+L
      END
C
************************************************************************
*   BNBDMK: mark boundary velocity dofs (KBD=1): vertices and edge
*           dofs of edges without neighbour (KADJ=0)
************************************************************************
      SUBROUTINE BNBDMK(K,KDFG,KADJ,NEL,NDOF,KBD)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION KDFG((K+1)**2,*),KADJ(NNVE,*),KBD(*)
C
      K1=K+1
      DO 5 I=1,NDOF
5     KBD(I)=0
      DO 20 IEL=1,NEL
      DO 20 IME=1,4
       IF (KADJ(IME,IEL).NE.0) GOTO 20
       DO 10 IS=0,K
        IF (IME.EQ.1) L=IS+1
        IF (IME.EQ.2) L=K+IS*K1+1
        IF (IME.EQ.3) L=(K-IS)+K*K1+1
        IF (IME.EQ.4) L=(K-IS)*K1+1
        KBD(KDFG(L,IEL))=1
10     CONTINUE
20    CONTINUE
C
      END
C
************************************************************************
*   BNSTL: matrix structure as linked lists; rows KDR(1..NLR,IEL),
*          columns KDC(1..NLC,IEL). BSQ: square, diagonal inserted
*          first. Returns NA (number of entries) or NA=-1 on overflow.
*   BNSTC: compress the lists into KLD(NROW+1), KCOL(NA)
************************************************************************
      SUBROUTINE BNSTL(NEL,NLR,KDR,NLC,KDC,NROW,BSQ,
     *                 KHEAD,KNEXT,KCOLT,NMAX,NA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KDR(NLR,*),KDC(NLC,*),KHEAD(*),KNEXT(*),KCOLT(*)
C
      NA=0
      DO 5 IR=1,NROW
       KHEAD(IR)=0
       IF (BSQ) THEN
        NA=NA+1
        KCOLT(NA)=IR
        KNEXT(NA)=0
        KHEAD(IR)=NA
       ENDIF
5     CONTINUE
C
      DO 30 IEL=1,NEL
      DO 30 I=1,NLR
       IR=KDR(I,IEL)
       DO 20 J=1,NLC
        JC=KDC(J,IEL)
        IP=KHEAD(IR)
        ILAST=0
10      IF (IP.NE.0) THEN
         IF (KCOLT(IP).EQ.JC) GOTO 20
         ILAST=IP
         IP=KNEXT(IP)
         GOTO 10
        ENDIF
        IF (NA.GE.NMAX) THEN
         NA=-1
         RETURN
        ENDIF
        NA=NA+1
        KCOLT(NA)=JC
        KNEXT(NA)=0
        IF (ILAST.EQ.0) THEN
         KHEAD(IR)=NA
        ELSE
         KNEXT(ILAST)=NA
        ENDIF
20     CONTINUE
30    CONTINUE
C
      END
C
      SUBROUTINE BNSTC(NROW,KHEAD,KNEXT,KCOLT,KLD,KCOL)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KHEAD(*),KNEXT(*),KCOLT(*),KLD(*),KCOL(*)
C
      IA=0
      DO 20 IR=1,NROW
       KLD(IR)=IA+1
       IP=KHEAD(IR)
10     IF (IP.NE.0) THEN
        IA=IA+1
        KCOL(IA)=KCOLT(IP)
        IP=KNEXT(IP)
        GOTO 10
       ENDIF
20    CONTINUE
      KLD(NROW+1)=IA+1
C
      END
C
************************************************************************
*   BNSTR: builds a CSR structure with FEAT workspace arrays
*          output: handles LLD, LCOL (KWORK), number of entries NA
************************************************************************
      SUBROUTINE BNSTR(NEL,NLR,KDR,NLC,KDC,NROW,BSQ,LLD,LCOL,NA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNWORK=1)
      DIMENSION KDR(*),KDC(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      NMAX=NEL*NLR*NLC+NROW
      CALL ZNEW(NROW,3,LHEAD,'KHEAD ')
      CALL ZNEW(NMAX,3,LNEXT,'KNEXT ')
      CALL ZNEW(NMAX,3,LCOLT,'KCOLT ')
      IF (IER.NE.0) STOP 'BNSTR: not enough workspace'
      CALL BNSTL(NEL,NLR,KDR,NLC,KDC,NROW,BSQ,KWORK(L(LHEAD)),
     *           KWORK(L(LNEXT)),KWORK(L(LCOLT)),NMAX,NA)
      IF (NA.LT.0) STOP 'BNSTR: structure overflow'
      CALL ZNEW(NROW+1,3,LLD ,'KLDBN ')
      CALL ZNEW(NA    ,3,LCOL,'KCOLBN')
      IF (IER.NE.0) STOP 'BNSTR: not enough workspace'
      CALL BNSTC(NROW,KWORK(L(LHEAD)),KWORK(L(LNEXT)),
     *           KWORK(L(LCOLT)),KWORK(L(LLD)),KWORK(L(LCOL)))
      CALL ZDISP(0,LCOLT,'KCOLT ')
      CALL ZDISP(0,LNEXT,'KNEXT ')
      CALL ZDISP(0,LHEAD,'KHEAD ')
C
      END
C
************************************************************************
*   BNFIND: position of entry (IR,JC) in a CSR matrix (0 if missing)
************************************************************************
      INTEGER*8 FUNCTION BNFIND(KLD,KCOL,IR,JC)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KLD(*),KCOL(*)
      DO 10 IA=KLD(IR),KLD(IR+1)-1
       IF (KCOL(IA).EQ.JC) THEN
        BNFIND=IA
        RETURN
       ENDIF
10    CONTINUE
      BNFIND=0
      END
C
************************************************************************
*   BNASM: assembly of the Bernstein matrices
*     DM (i,j) =  int phi_j phi_i                      (mass)
*     DA (i,j) =  int grad phi_j . grad phi_i          (Laplacian)
*     DB1(i,q) = -int psi_q d/dx phi_i                 (gradient B)
*     DB2(i,q) = -int psi_q d/dy phi_i
*   DAREA = |Omega_h| (sum of the cubature weights)
************************************************************************
      SUBROUTINE BNASM(NEL,DCORVG,KVERT,KDFG,KDFP,KLDA,KCOLA,
     *                 KLDB,KCOLB,DM,DA,DB1,DB2,DAREA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG(NLDF,*)
      DIMENSION KDFP(NPL,*),KLDA(*),KCOLA(*),KLDB(*),KCOLB(*)
      DIMENSION DM(*),DA(*),DB1(*),DB2(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      DIMENSION DEM(NNLD,NNLD),DEA(NNLD,NNLD)
      DIMENSION DEB1(NNLD,NNPL),DEB2(NNLD,NNPL)
      INTEGER*8 BNFIND
C
      NGP =NG*NG
      CALL BNGAUS(NG,XG,WG)
      DAREA=0D0
C
      DO 100 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNGSET(IEL)
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
C
C *** element matrices
       DO 10 J=1,NLDF
       DO 10 I=1,NLDF
        DEM(I,J)=0D0
        DEA(I,J)=0D0
10     CONTINUE
       DO 11 IQ=1,NPL
       DO 11 I=1,NLDF
        DEB1(I,IQ)=0D0
        DEB2(I,IQ)=0D0
11     CONTINUE
       DO 30 IG=1,NGP
        W=DWT(IG)
        DAREA=DAREA+W
        DO 20 I=1,NLDF
         DO 15 J=1,NLDF
          DEM(I,J)=DEM(I,J)+W*DPHI(I,IG)*DPHI(J,IG)
          DEA(I,J)=DEA(I,J)+W*(DPX(I,IG)*DPX(J,IG)
     *                        +DPY(I,IG)*DPY(J,IG))
15       CONTINUE
         DO 16 IQ=1,NPL
          DEB1(I,IQ)=DEB1(I,IQ)-W*DPSI(IQ,IG)*DPX(I,IG)
          DEB2(I,IQ)=DEB2(I,IQ)-W*DPSI(IQ,IG)*DPY(I,IG)
16       CONTINUE
20      CONTINUE
30     CONTINUE
C
C *** scatter into the global matrices
       DO 50 I=1,NLDF
        IR=KDFG(I,IEL)
        DO 40 J=1,NLDF
         IA=BNFIND(KLDA,KCOLA,IR,KDFG(J,IEL))
         IF (IA.EQ.0) STOP 'BNASM: entry not in structure A'
         DM(IA)=DM(IA)+DEM(I,J)
         DA(IA)=DA(IA)+DEA(I,J)
40      CONTINUE
        DO 45 IQ=1,NPL
         IA=BNFIND(KLDB,KCOLB,IR,KDFP(IQ,IEL))
         IF (IA.EQ.0) STOP 'BNASM: entry not in structure B'
         DB1(IA)=DB1(IA)+DEB1(I,IQ)
         DB2(IA)=DB2(IA)+DEB2(I,IQ)
45      CONTINUE
50     CONTINUE
100   CONTINUE
      CALL BNGSET(0)
C
      END
C
************************************************************************
*   BNRHS: DB(i) = int F(x,y) phi_i  (F external function)
************************************************************************
      SUBROUTINE BNRHS(NEL,DCORVG,KVERT,KDFG,F,DB)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG(NLDF,*),DB(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      EXTERNAL F
C
      CALL BNGAUS(NG,XG,WG)
      DO 30 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
       DO 20 IG=1,NG*NG
        FW=DWT(IG)*F(DXQ(IG),DYQ(IG))
        DO 10 I=1,NLDF
         DB(KDFG(I,IEL))=DB(KDFG(I,IEL))+FW*DPHI(I,IG)
10      CONTINUE
20     CONTINUE
30    CONTINUE
C
      END
C
************************************************************************
*   BNL2E: L2 error |u_h - U| and H1 semi-norm error |grad(u_h - U)|
*          of the Bernstein function with coefficients DC
*          (U, UX, UY external functions)
************************************************************************
      SUBROUTINE BNL2E(NEL,DCORVG,KVERT,KDFG,DC,U,UX,UY,ERRL2,ERRH1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KDFG(NLDF,*),DC(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      EXTERNAL U,UX,UY
C
      CALL BNGAUS(NG,XG,WG)
      ERRL2=0D0
      ERRH1=0D0
      DO 30 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
       DO 20 IG=1,NG*NG
        UH=0D0
        UHX=0D0
        UHY=0D0
        DO 10 I=1,NLDF
         C=DC(KDFG(I,IEL))
         UH =UH +C*DPHI(I,IG)
         UHX=UHX+C*DPX(I,IG)
         UHY=UHY+C*DPY(I,IG)
10      CONTINUE
        X=DXQ(IG)
        Y=DYQ(IG)
        ERRL2=ERRL2+DWT(IG)*(UH-U(X,Y))**2
        ERRH1=ERRH1+DWT(IG)*((UHX-UX(X,Y))**2+(UHY-UY(X,Y))**2)
20     CONTINUE
30    CONTINUE
      ERRL2=SQRT(ERRL2)
      ERRH1=SQRT(ERRH1)
C
      END
C
************************************************************************
*   BNPL2: element-wise L2 projection of P(x,y) onto P(k-1)disc
*          (local Gram system, Gauss elimination) -> coefficients DPC
*          and L2 error ERRL2.  Tests the pressure basis.
************************************************************************
      SUBROUTINE BNPL2(NEL,DCORVG,KVERT,P,DPC,ERRL2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      PARAMETER (NNBK=8,NNLD=(NNBK+1)**2,NNPL=NNBK*(NNBK+1)/2,
     *           NNG=NNBK+2,NNGP=NNG*NNG,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),DPC(*)
      DIMENSION XG(NNG),WG(NNG),DXV(4),DYV(4)
      DIMENSION DPHI(NNLD,NNGP),DPX(NNLD,NNGP),DPY(NNLD,NNGP)
      DIMENSION DPSI(NNPL,NNGP),DWT(NNGP),DXQ(NNGP),DYQ(NNGP)
      DIMENSION DG(NNPL,NNPL),DR(NNPL)
      EXTERNAL P
C
      CALL BNGAUS(NG,XG,WG)
      ERRL2=0D0
      DO 50 IEL=1,NEL
       DO 5 IV=1,4
        DXV(IV)=DCORVG(1,KVERT(IV,IEL))
        DYV(IV)=DCORVG(2,KVERT(IV,IEL))
5      CONTINUE
       CALL BNEVAL(XG,WG,DXV,DYV,DPHI,DPX,DPY,DPSI,DWT,DXQ,DYQ)
       DO 10 I=1,NPL
        DR(I)=0D0
        DO 10 J=1,NPL
10     DG(I,J)=0D0
       DO 20 IG=1,NG*NG
        PW=DWT(IG)*P(DXQ(IG),DYQ(IG))
        DO 15 I=1,NPL
         DR(I)=DR(I)+PW*DPSI(I,IG)
         DO 15 J=1,NPL
15      DG(I,J)=DG(I,J)+DWT(IG)*DPSI(I,IG)*DPSI(J,IG)
20     CONTINUE
C *** Gauss elimination without pivoting (Gram matrix is SPD)
       DO 30 IC=1,NPL
        DO 30 IR=IC+1,NPL
         F=DG(IR,IC)/DG(IC,IC)
         DO 25 J=IC,NPL
25       DG(IR,J)=DG(IR,J)-F*DG(IC,J)
         DR(IR)=DR(IR)-F*DR(IC)
30     CONTINUE
       DO 35 IR=NPL,1,-1
        S=DR(IR)
        DO 34 J=IR+1,NPL
34      S=S-DG(IR,J)*DR(J)
        DR(IR)=S/DG(IR,IR)
35     CONTINUE
       DO 40 I=1,NPL
40     DPC((IEL-1)*NPL+I)=DR(I)
       DO 45 IG=1,NG*NG
        PH=0D0
        DO 42 I=1,NPL
42      PH=PH+DR(I)*DPSI(I,IG)
        ERRL2=ERRL2+DWT(IG)*(PH-P(DXQ(IG),DYQ(IG)))**2
45     CONTINUE
50    CONTINUE
      ERRL2=SQRT(ERRL2)
C
      END
C
************************************************************************
*   BNAX: DY = A*DX for a CSR matrix (NROW rows)
************************************************************************
      SUBROUTINE BNAX(DA,KLD,KCOL,NROW,DX,DY)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),KLD(*),KCOL(*),DX(*),DY(*)
      DO 20 IR=1,NROW
       S=0D0
       DO 10 IA=KLD(IR),KLD(IR+1)-1
10     S=S+DA(IA)*DX(KCOL(IA))
       DY(IR)=S
20    CONTINUE
      END
C
************************************************************************
*   BNCG: Jacobi-preconditioned CG for a symmetric CSR matrix with the
*         diagonal first in each row; dofs with KBD=1 are kept fixed
*         (Dirichlet). DR,DP,DQ: auxiliary vectors.
*         Stops at |r|<=EPS*|r0|; ITE = iterations (-ITE if no conv.)
************************************************************************
      SUBROUTINE BNCG(DA,KLD,KCOL,NEQ,DX,DB,KBD,DR,DP,DQ,
     *                ITMAX,EPS,ITE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DA(*),KLD(*),KCOL(*),DX(*),DB(*),KBD(*)
      DIMENSION DR(*),DP(*),DQ(*)
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
      IF (RN0.EQ.0D0) RETURN
      RHO0=0D0
      DO 15 I=1,NEQ
       DP(I)=DR(I)/DA(KLD(I))
       RHO0=RHO0+DR(I)*DP(I)
15    CONTINUE
C
      DO 50 ITE=1,ITMAX
       CALL BNAX(DA,KLD,KCOL,NEQ,DP,DQ)
       PQ=0D0
       DO 20 I=1,NEQ
        IF (KBD(I).EQ.1) DQ(I)=0D0
        PQ=PQ+DP(I)*DQ(I)
20     CONTINUE
       ALPHA=RHO0/PQ
       RN=0D0
       RHO1=0D0
       DO 30 I=1,NEQ
        DX(I)=DX(I)+ALPHA*DP(I)
        DR(I)=DR(I)-ALPHA*DQ(I)
        RN=RN+DR(I)**2
        RHO1=RHO1+DR(I)**2/DA(KLD(I))
30     CONTINUE
       IF (SQRT(RN).LE.EPS*RN0) RETURN
       GAM=RHO1/RHO0
       RHO0=RHO1
       DO 40 I=1,NEQ
40     DP(I)=DR(I)/DA(KLD(I))+GAM*DP(I)
50    CONTINUE
      ITE=-ITMAX
C
      END
C
************************************************************************
*   Exact functions for BNTEST
*   u = sin(3x)cos(2y)+xy, -Lap u = 13 sin(3x)cos(2y), p = cos(2x)e^y
************************************************************************
      DOUBLE PRECISION FUNCTION BNUEX(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BNUEX=SIN(3D0*X)*COS(2D0*Y)+X*Y
      END
      DOUBLE PRECISION FUNCTION BNUEXX(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BNUEXX=3D0*COS(3D0*X)*COS(2D0*Y)+Y
      END
      DOUBLE PRECISION FUNCTION BNUEXY(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BNUEXY=-2D0*SIN(3D0*X)*SIN(2D0*Y)+X
      END
      DOUBLE PRECISION FUNCTION BNFEX(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BNFEX=13D0*SIN(3D0*X)*COS(2D0*Y)
      END
      DOUBLE PRECISION FUNCTION BNPEX(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BNPEX=COS(2D0*X)*EXP(Y)
      END
C
************************************************************************
*   BNTEST: stage-1 tests of the Bernstein Qk/P(k-1)disc pair on all
*           levels NLMIN..NLMAX (called from pp2d if IBERN=K>0)
*     1  number of dofs, every global dof used
*     2  sum(M) = sum(Ml) = |Omega_h|, min(Ml) > 0 (lumped mass)
*     3  A*1 = 0
*     4  B1*1 = B2*1 = 0 in interior dofs (partition of unity of psi)
*     5  L2 projection onto Qk: L2 error O(h^(k+1))
*     6  Poisson -Lap u = f, Dirichlet data = projection coefficients:
*        L2 error O(h^(k+1)), H1 error O(h^k)
*     7  L2 projection onto P(k-1)disc: L2 error O(h^k)
************************************************************************
      SUBROUTINE BNTEST(K,MFILE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      CHARACTER SUB*6,FMT*15,CPARAM*120
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
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      DIMENSION E5(NNLEV),E6(NNLEV),E6H(NNLEV),E7(NNLEV)
      EXTERNAL BNUEX,BNUEXX,BNUEXY,BNFEX,BNPEX
      DOUBLE PRECISION BNUEX,BNUEXX,BNUEXY,BNFEX,BNPEX
      SAVE
C
      SUB='BNTEST'
      IF ((K.LT.2).OR.(K.GT.NNBK)) THEN
       WRITE(MTERM,*) 'BNTEST: order K must be in 2..',NNBK
       RETURN
      ENDIF
      IFAM=0
      KORD=K
      NLDF=(K+1)**2
      NPL =K*(K+1)/2
      NG  =K+2
C
      WRITE(MTERM,1000) K,K-1
      WRITE(MFILE,1000) K,K-1
C
      DO 100 ILEV=NLMIN,NLMAX
       CALL SETLEV(1)
       NDOF =NVT+(K-1)*NMT+(K-1)**2*NEL
       NPDOF=NPL*NEL
C
C *** dof maps and boundary marker
       CALL ZNEW(NLDF*NEL,3,LDFG,'KDFG  ')
       CALL ZNEW(NPL *NEL,3,LDFP,'KDFP  ')
       CALL ZNEW(NDOF    ,3,LBD ,'KBD   ')
       IF (IER.NE.0) STOP 'BNTEST: not enough workspace'
       CALL BNDFMP(K,KWORK(L(LVERT)),KWORK(L(LMID)),NEL,NVT,NMT,
     *             KWORK(L(LDFG)))
       CALL BNDFPP(NPL,NEL,KWORK(L(LDFP)))
C *** test 1: every dof 1..NDOF used (KBD as counter)
       DO 5 I=1,NDOF
5      KWORK(L(LBD)+I-1)=0
       IBAD=0
       DO 6 I=1,NLDF*NEL
        IG=KWORK(L(LDFG)+I-1)
        IF ((IG.LT.1).OR.(IG.GT.NDOF)) THEN
         IBAD=IBAD+1
        ELSE
         KWORK(L(LBD)+IG-1)=1
        ENDIF
6      CONTINUE
       NUNUSD=0
       DO 7 I=1,NDOF
7      IF (KWORK(L(LBD)+I-1).EQ.0) NUNUSD=NUNUSD+1
       CALL BNBDMK(K,KWORK(L(LDFG)),KWORK(L(LADJ)),NEL,NDOF,
     *             KWORK(L(LBD)))
       NBDOF=0
       DO 8 I=1,NDOF
8      NBDOF=NBDOF+KWORK(L(LBD)+I-1)
C
C *** matrix structures and assembly
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NLDF,KWORK(L(LDFG)),NDOF,
     *            .TRUE.,LLDA,LCOLA,NA)
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NPL,KWORK(L(LDFP)),NDOF,
     *            .FALSE.,LLDB,LCOLB,NB)
       CALL ZNEW(NA,1,LM ,'DMBN  ')
       CALL ZNEW(NA,1,LA ,'DABN  ')
       CALL ZNEW(NB,1,LB1,'DB1BN ')
       CALL ZNEW(NB,1,LB2,'DB2BN ')
       IF (IER.NE.0) STOP 'BNTEST: not enough workspace'
       CALL BNASM(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),KWORK(L(LDFP)),
     *            KWORK(L(LLDA)),KWORK(L(LCOLA)),
     *            KWORK(L(LLDB)),KWORK(L(LCOLB)),
     *            DWORK(L(LM)),DWORK(L(LA)),DWORK(L(LB1)),DWORK(L(LB2)),
     *            DAREA)
C
C *** test 2/3: sum(M), lumped mass, A*1
       SUMM=0D0
       DMLMIN=1D30
       DMLMAX=-1D30
       DA1MAX=0D0
       DO 20 IR=1,NDOF
        DML=0D0
        DA1=0D0
        DO 15 IA=KWORK(L(LLDA)+IR-1),KWORK(L(LLDA)+IR)-1
         DML=DML+DWORK(L(LM)+IA-1)
         DA1=DA1+DWORK(L(LA)+IA-1)
15      CONTINUE
        SUMM=SUMM+DML
        DMLMIN=MIN(DMLMIN,DML)
        DMLMAX=MAX(DMLMAX,DML)
        DA1MAX=MAX(DA1MAX,ABS(DA1))
20     CONTINUE
C *** test 4: B*1 in interior dofs
       DB1MAX=0D0
       DO 30 IR=1,NDOF
        IF (KWORK(L(LBD)+IR-1).EQ.1) GOTO 30
        S1=0D0
        S2=0D0
        DO 25 IA=KWORK(L(LLDB)+IR-1),KWORK(L(LLDB)+IR)-1
         S1=S1+DWORK(L(LB1)+IA-1)
         S2=S2+DWORK(L(LB2)+IA-1)
25      CONTINUE
        DB1MAX=MAX(DB1MAX,ABS(S1),ABS(S2))
30     CONTINUE
C
C *** test 5: L2 projection M c = (u,phi)
       CALL ZNEW(NDOF,1,LC  ,'DCBN  ')
       CALL ZNEW(NDOF,1,LCP ,'DCPBN ')
       CALL ZNEW(NDOF,1,LRHS,'DRHSBN')
       CALL ZNEW(NDOF,1,LR  ,'DRBN  ')
       CALL ZNEW(NDOF,1,LP  ,'DPBN  ')
       CALL ZNEW(NDOF,1,LQ  ,'DQBN  ')
       CALL ZNEW(NDOF,3,LNOB,'KNOBN ')
       CALL ZNEW(NPDOF,1,LPC,'DPCBN ')
       IF (IER.NE.0) STOP 'BNTEST: not enough workspace'
       CALL BNRHS(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),BNUEX,DWORK(L(LRHS)))
       CALL BNCG(DWORK(L(LM)),KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,
     *           DWORK(L(LCP)),DWORK(L(LRHS)),KWORK(L(LNOB)),
     *           DWORK(L(LR)),DWORK(L(LP)),DWORK(L(LQ)),
     *           10000,1D-13,ITE5)
       CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),DWORK(L(LCP)),BNUEX,BNUEXX,BNUEXY,
     *            E5(ILEV),DUMMY)
C
C *** test 6: Poisson, boundary dofs from the projection
       CALL LCL1(DWORK(L(LRHS)),NDOF)
       CALL BNRHS(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),BNFEX,DWORK(L(LRHS)))
       DO 40 I=1,NDOF
        IF (KWORK(L(LBD)+I-1).EQ.1) THEN
         DWORK(L(LC)+I-1)=DWORK(L(LCP)+I-1)
        ELSE
         DWORK(L(LC)+I-1)=0D0
        ENDIF
40     CONTINUE
       CALL BNCG(DWORK(L(LA)),KWORK(L(LLDA)),KWORK(L(LCOLA)),NDOF,
     *           DWORK(L(LC)),DWORK(L(LRHS)),KWORK(L(LBD)),
     *           DWORK(L(LR)),DWORK(L(LP)),DWORK(L(LQ)),
     *           100000,1D-13,ITE6)
       CALL BNL2E(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LDFG)),DWORK(L(LC)),BNUEX,BNUEXX,BNUEXY,
     *            E6(ILEV),E6H(ILEV))
C
C *** test 7: pressure projection
       CALL BNPL2(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),BNPEX,
     *            DWORK(L(LPC)),E7(ILEV))
C
       DO 50 IU=1,2
        IF (IU.EQ.1) MU=MTERM
        IF (IU.EQ.2) MU=MFILE
        WRITE(MU,1010) ILEV,NEL,NDOF,NBDOF,NPDOF,NA,NB,NUNUSD,IBAD
        WRITE(MU,1020) DAREA,SUMM,DMLMIN,DMLMAX,DA1MAX,DB1MAX
        WRITE(MU,1030) E5(ILEV),ITE5,E6(ILEV),E6H(ILEV),ITE6,E7(ILEV)
        IF (ILEV.GT.NLMIN) WRITE(MU,1040)
     *   LOG(E5 (ILEV-1)/E5 (ILEV))/LOG(2D0),
     *   LOG(E6 (ILEV-1)/E6 (ILEV))/LOG(2D0),
     *   LOG(E6H(ILEV-1)/E6H(ILEV))/LOG(2D0),
     *   LOG(E7 (ILEV-1)/E7 (ILEV))/LOG(2D0)
50     CONTINUE
C
       CALL ZDISP(0,LPC ,'DPCBN ')
       CALL ZDISP(0,LNOB,'KNOBN ')
       CALL ZDISP(0,LQ  ,'DQBN  ')
       CALL ZDISP(0,LP  ,'DPBN  ')
       CALL ZDISP(0,LR  ,'DRBN  ')
       CALL ZDISP(0,LRHS,'DRHSBN')
       CALL ZDISP(0,LCP ,'DCPBN ')
       CALL ZDISP(0,LC  ,'DCBN  ')
       CALL ZDISP(0,LB2 ,'DB2BN ')
       CALL ZDISP(0,LB1 ,'DB1BN ')
       CALL ZDISP(0,LA  ,'DABN  ')
       CALL ZDISP(0,LM  ,'DMBN  ')
       CALL ZDISP(0,LCOLB,'KCOLBN')
       CALL ZDISP(0,LLDB ,'KLDBN ')
       CALL ZDISP(0,LCOLA,'KCOLBN')
       CALL ZDISP(0,LLDA ,'KLDBN ')
       CALL ZDISP(0,LBD ,'KBD   ')
       CALL ZDISP(0,LDFP,'KDFP  ')
       CALL ZDISP(0,LDFG,'KDFG  ')
100   CONTINUE
C
1000  FORMAT(/' BERNSTEIN TEST  Q',I1,' / P',I1,'disc'/
     * ' ',70('-'))
1010  FORMAT(' LEV',I2,' NEL=',I7,' NDOF=',I8,' NBD=',I6,' NPDOF=',I8/
     * '       NA=',I9,' NB=',I9,' unused dofs=',I4,' bad=',I4)
1020  FORMAT('   |Omega|=',D14.7,' sum(M)=',D14.7,
     * ' Ml min/max=',2D11.3/
     * '   max|A*1|=',D10.3,' max|B*1| (interior)=',D10.3)
1030  FORMAT('   L2proj err=',D11.4,' (CG',I6,')',
     * '  Poisson L2/H1=',2D11.4,' (CG',I7,')'/
     * '   P(k-1) proj err=',D11.4)
1040  FORMAT('   rates:  L2proj ',F5.2,'  Poisson L2 ',F5.2,
     * '  H1 ',F5.2,'  pressure ',F5.2)
C
      END
