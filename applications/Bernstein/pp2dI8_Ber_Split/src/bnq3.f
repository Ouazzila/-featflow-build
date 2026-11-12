************************************************************************
*   bnq3.f: Bernstein Q3/P2disc inside pp2d (IELT=5), and generic
*   (k=2,3) versions of the k-dependent parts of bnpp.f
*
*   Velocity: Bernstein Q3 as a FEAT element (EB30), numbered by NDFGL
*   like E014 (type 14): vertices, 2 dofs per edge (NVT+2(IMT-1)+s,
*   s=1 next to the lower vertex number), 4 per element
*   (NVT+2NMT+4(IEL-1)+1..4). This is the numbering of BNDFMP (bern.f)
*   for k=3 if the interior dofs 13..16 are (i,j)=(1,1),(2,1),(1,2),
*   (2,2) (E014 orders them counter-clockwise; they are element-local,
*   so only EB30 and BNQKDF have to agree).
*   Pressure: P2disc, Bernstein basis of the reference triangle
*   (-1,-1),(3,-1),(-1,3) (BNPB with N=2), 6 dofs per element.
*
*   The degree k is IBNPP (COMMON /BNLEV/): 0 = off, 2 = Q2, 3 = Q3.
*   The Q2 code of bnpp.f is kept unchanged (bit-identical results);
*   the generic routines here are used for k=3 and, for k=2, only to
*   check them against the Q2 matrices (BNGCHK).
************************************************************************
*
************************************************************************
*   EBK: element dispatcher for FEAT (XMAP7, XMAB07, XVB0, SUPWDG,
*        ELPQU, BNFVOL): EB20 (IBNPP=2) or EB30 (IBNPP=3)
************************************************************************
      SUBROUTINE EBK(XI1,XI2,IPAR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNLEV=9)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      SAVE /BNLEV/
C
      IF (IBNPP.EQ.3) THEN
       CALL EB30(XI1,XI2,IPAR)
      ELSE
       CALL EB20(XI1,XI2,IPAR)
      ENDIF
      END
C
************************************************************************
*   EB30: Bernstein Q3 on the reference square [-1,1]^2, parametric
*         (bilinear map), FEAT /ELEM/ interface as E014:
*     IPAR = -1 : return element type 14 (Q3 numbering in NDFGL)
*     IPAR = -2 : evaluate on all cubature points of /CUB/
*     IPAR = -3 : use the values of cubature point ICUBP
*     IPAR =  0 : evaluate in (XI1,XI2)
*   Local order: v1..v4, edge m (v_m -> v_m+1): 2m+3 next to v_m,
*   2m+4 next to v_m+1, interior 13..16 = (1,1),(2,1),(1,2),(2,2)
*   phi = B_i^3(xi) B_j^3(eta) (BNLORD gives i,j)
************************************************************************
      SUBROUTINE EB30(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=16)
      DIMENSION DHELP(NDOFL,3,NNCUBP)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
      SAVE DHELP
C
      SUB='EB30'
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2)
10    ICUBP0=1
      CALL EBKA(3,XI1,XI2,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element (numbering as Q3, E014)
1     IPAR=14
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
2     DO 11 ICUBP0=1,NCUBP
      CALL EBKA(3,DXI(ICUBP0,1),DXI(ICUBP0,2),DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combinations of the values in DHELP
C *** for the actual element
3     ICUBP0=ICUBP
C
4     IF (BDER(4).OR.BDER(5).OR.BDER(6)) THEN
       CALL WERR(-131,'EB30  ')
       GOTO 99999
      ENDIF
      IF (BDER(1)) THEN
       DO 1301 IDFL=1,NDOFL
1301   DBAS(IDFL,1)=DHELP(IDFL,1,ICUBP0)
      ENDIF
      XJ1=1D0/DETJ
      IF (BDER(2)) THEN
       DO 1303 IDFL=1,NDOFL
1303   DBAS(IDFL,2)=(DHELP(IDFL,2,ICUBP0)*DJAC(2,2)
     *              -DHELP(IDFL,3,ICUBP0)*DJAC(2,1))*XJ1
      ENDIF
      IF (BDER(3)) THEN
       DO 1305 IDFL=1,NDOFL
1305   DBAS(IDFL,3)=(-DHELP(IDFL,2,ICUBP0)*DJAC(1,2)
     *               +DHELP(IDFL,3,ICUBP0)*DJAC(1,1))*XJ1
      ENDIF
C
99999 END
C
************************************************************************
*   EBKA: values and reference derivatives of the (K+1)^2 Bernstein
*         Qk basis functions (element-local order BNLORD) in (X1,X2),
*         DHELP(NL,3,ICUBP0), NL=(K+1)^2 is the leading dimension
************************************************************************
      SUBROUTINE EBKA(K,X1,X2,DHELP,ICUBP0)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8,NNL=16)
      DIMENSION DHELP((K+1)**2,3,*)
      DIMENSION DB(0:NNBK),DDB(0:NNBK),DC(0:NNBK),DDC(0:NNBK)
      DIMENSION KI(NNL),KJ(NNL)
C
      CALL BNLORD(K,KI,KJ)
      CALL BNB1D(K,X1,DB,DDB)
      CALL BNB1D(K,X2,DC,DDC)
      DO 10 IDFL=1,(K+1)**2
       I=KI(IDFL)
       J=KJ(IDFL)
       DHELP(IDFL,1,ICUBP0)= DB(I)* DC(J)
       DHELP(IDFL,2,ICUBP0)=DDB(I)* DC(J)
       DHELP(IDFL,3,ICUBP0)= DB(I)*DDC(J)
10    CONTINUE
      END
C
************************************************************************
*   BNLORD: (i,j) of the element-local dofs of EB20 (K=2) / EB30 (K=3)
************************************************************************
      SUBROUTINE BNLORD(K,KI,KJ)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KI(*),KJ(*),KI2(9),KJ2(9),KI3(16),KJ3(16)
      DATA KI2/0,2,2,0, 1,2,1,0, 1/
      DATA KJ2/0,0,2,2, 0,1,2,1, 1/
      DATA KI3/0,3,3,0, 1,2, 3,3, 2,1, 0,0, 1,2,1,2/
      DATA KJ3/0,0,3,3, 0,0, 1,2, 3,3, 2,1, 1,1,2,2/
      IF (K.EQ.2) THEN
       DO 10 I=1,9
        KI(I)=KI2(I)
10      KJ(I)=KJ2(I)
      ELSE IF (K.EQ.3) THEN
       DO 20 I=1,16
        KI(I)=KI3(I)
20      KJ(I)=KJ3(I)
      ELSE
       STOP 'BNLORD: only K=2,3'
      ENDIF
      END
C
************************************************************************
*   BNQKDF: global dofs of element IEL in element-local order
*     (K=2: BNQ2DF; K=3: as NDFGL type 14 with EB30 order)
************************************************************************
      SUBROUTINE BNQKDF(K,IEL,KVERT,KMID,NVT,NMT,KDF)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KVERT(4,*),KMID(4,*),KDF(*)
      IF (K.EQ.2) THEN
       CALL BNQ2DF(IEL,KVERT,KMID,NVT,NMT,KDF)
       RETURN
      ENDIF
      DO 10 I=1,4
       KDF(I)=KVERT(I,IEL)
       IMT=KMID(I,IEL)-NVT
       IV1=KVERT(I,IEL)
       IV2=KVERT(MOD(I,4)+1,IEL)
       IB=NVT+2*(IMT-1)
       IF (IV1.LT.IV2) THEN
        KDF(2*I+3)=IB+1
        KDF(2*I+4)=IB+2
       ELSE
        KDF(2*I+3)=IB+2
        KDF(2*I+4)=IB+1
       ENDIF
10    CONTINUE
      IOFF=NVT+2*NMT+4*(IEL-1)
      DO 20 I=1,4
20    KDF(12+I)=IOFF+I
      END
C
************************************************************************
*   COMMON /BNTRG/ (generic transfer, BNGTM):
*     GTU(lf,m,k) : Bernstein coefficient of fine local dof lf (child k)
*                   of the coarse basis function m (element-local order)
*     GSP(lf,q,k) : P(k-1)disc coefficient lf on child k of coarse psi_q
*     GLI(i,a)    : 1D Lagrange (nodes -1+2a/k) -> Bernstein coeff. i
*     KNC(a,b), GWN(lf,a,b): child containing the coarse Lagrange node
*                   (a,b) and the child basis values there
*     KGT         : degree the arrays were built for (0 = none)
************************************************************************
*
************************************************************************
*   BNGTM: build /BNTRG/ for degree K (exact, the child maps are
*     affine in the parent reference coordinates)
************************************************************************
      SUBROUTINE BNGTM(K)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNL=16,NNP=6,NNBK=8)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION CX(4),CY(4),QX(4),QY(4),DHELP(NNL,3,1)
      DIMENSION DV(0:3,0:3,NNL),DB(0:NNBK),DDB(0:NNBK)
      DIMENSION DG(4,4),DR(4),GP(NNP,NNP),GPI(NNP,NNP),DPV(NNP)
      DIMENSION PS(NNP),PT(NNP),GQ(NNP,NNP),DRP(NNP)
      DIMENSION KI(NNL),KJ(NNL)
      DATA CX/-1D0,1D0,1D0,-1D0/, CY/-1D0,-1D0,1D0,1D0/
      SAVE
C
      IF (KGT.EQ.K) RETURN
      IF ((K.LT.2).OR.(K.GT.3)) STOP 'BNGTM: only K=2,3'
      NL=(K+1)**2
      N=K-1
      NPL=(N+1)*(N+2)/2
      CALL BNLORD(K,KI,KJ)
C
C *** 1D: GLI = inverse of V(a,i) = B_i(s_a), s_a = -1+2a/K
      DO 10 IA=0,K
       DO 5 IC=0,K
        DO 4 JR=0,K
         S=-1D0+2D0*DBLE(JR)/DBLE(K)
         CALL BNB1D(K,S,DB,DDB)
         DG(JR+1,IC+1)=DB(IC)
4       CONTINUE
5      CONTINUE
       DO 6 JR=0,K
6      DR(JR+1)=0D0
       DR(IA+1)=1D0
       CALL BNGSOL(DG,4,K+1,DR)
       DO 7 IC=0,K
7      GLI(IC,IA)=DR(IC+1)
10    CONTINUE
C
C *** pressure lattice points of P_N on the reference triangle
      IP=0
      DO 12 IA2=0,N
      DO 12 IA3=0,N-IA2
       IP=IP+1
       IF (N.EQ.0) THEN
        PS(IP)=0D0
        PT(IP)=0D0
       ELSE
        PS(IP)=-1D0+4D0*DBLE(IA2)/DBLE(N)
        PT(IP)=-1D0+4D0*DBLE(IA3)/DBLE(N)
       ENDIF
12    CONTINUE
      DO 14 IP=1,NPL
       CALL BNPB(N,PS(IP),PT(IP),DPV)
       DO 13 LF=1,NPL
13     GP(IP,LF)=DPV(LF)
14    CONTINUE
C *** GPI = GP^-1 (column by column)
      DO 16 JC=1,NPL
       DO 15 I=1,NPL
        DRP(I)=0D0
        DO 15 J=1,NPL
15      GQ(I,J)=GP(I,J)
       DRP(JC)=1D0
       CALL BNGSOL(GQ,NNP,NPL,DRP)
       DO 16 I=1,NPL
16     GPI(I,JC)=DRP(I)
C
      DO 100 KC=1,4
       KP=MOD(KC,4)+1
       KM=MOD(KC+2,4)+1
       QX(1)=CX(KC)
       QY(1)=CY(KC)
       QX(2)=0.5D0*(CX(KC)+CX(KP))
       QY(2)=0.5D0*(CY(KC)+CY(KP))
       QX(3)=0D0
       QY(3)=0D0
       QX(4)=0.5D0*(CX(KC)+CX(KM))
       QY(4)=0.5D0*(CY(KC)+CY(KM))
C *** velocity: coarse basis values at the child Lagrange nodes
       DO 20 IA=0,K
       DO 20 IB=0,K
        SA=-1D0+2D0*DBLE(IA)/DBLE(K)
        SB=-1D0+2D0*DBLE(IB)/DBLE(K)
        CALL BNCHLD(QX,QY,SA,SB,X,Y)
        CALL EBKA(K,X,Y,DHELP,1)
        DO 18 M=1,NL
18      DV(IA,IB,M)=DHELP(M,1,1)
20     CONTINUE
       DO 30 M=1,NL
       DO 30 LF=1,NL
        I=KI(LF)
        J=KJ(LF)
        S=0D0
        DO 25 IA=0,K
        DO 25 IB=0,K
25      S=S+GLI(I,IA)*GLI(J,IB)*DV(IA,IB,M)
        GTU(LF,M,KC)=S
30     CONTINUE
C *** pressure: coarse psi_q at the child lattice points -> child coeff.
       DO 50 IQ=1,NPL
        DO 45 IP=1,NPL
         CALL BNCHLD(QX,QY,PS(IP),PT(IP),X,Y)
         CALL BNPB(N,X,Y,DPV)
         DRP(IP)=DPV(IQ)
45      CONTINUE
        DO 48 LF=1,NPL
         S=0D0
         DO 47 IP=1,NPL
47       S=S+GPI(LF,IP)*DRP(IP)
         GSP(LF,IQ,KC)=S
48      CONTINUE
50     CONTINUE
100   CONTINUE
C
C *** coarse Lagrange nodes -> child (first child that contains it)
C     and child basis values (solution restriction BNRSLG)
      DO 200 IA=0,K
      DO 200 IB=0,K
       XA=-1D0+2D0*DBLE(IA)/DBLE(K)
       YB=-1D0+2D0*DBLE(IB)/DBLE(K)
       KNC(IA,IB)=0
       DO 150 KC=1,4
        KP=MOD(KC,4)+1
        KM=MOD(KC+2,4)+1
        QX(1)=CX(KC)
        QY(1)=CY(KC)
        QX(2)=0.5D0*(CX(KC)+CX(KP))
        QY(2)=0.5D0*(CY(KC)+CY(KP))
        QX(4)=0.5D0*(CX(KC)+CX(KM))
        QY(4)=0.5D0*(CY(KC)+CY(KM))
C ***   affine: (x,y) = C + s*ES + t*ET, C = centre of the child
        XC=0.5D0*QX(1)
        YC=0.5D0*QY(1)
        ESX=0.5D0*(QX(2)-QX(1))
        ESY=0.5D0*(QY(2)-QY(1))
        ETX=0.5D0*(QX(4)-QX(1))
        ETY=0.5D0*(QY(4)-QY(1))
        DET=ESX*ETY-ESY*ETX
        S=((XA-XC)*ETY-(YB-YC)*ETX)/DET
        T=(ESX*(YB-YC)-ESY*(XA-XC))/DET
        IF ((ABS(S).LE.1D0+1D-12).AND.(ABS(T).LE.1D0+1D-12)) THEN
         KNC(IA,IB)=KC
         CALL EBKA(K,S,T,DHELP,1)
         DO 140 LF=1,NL
140      GWN(LF,IA,IB)=DHELP(LF,1,1)
         GOTO 200
        ENDIF
150    CONTINUE
       STOP 'BNGTM: coarse node in no child'
200   CONTINUE
C
      KGT=K
      END
C
************************************************************************
*   BNGCHK: K=2 only: compare the generic matrices with the Q2 ones of
*     BNTMAT (DTU, DSP), result to MFILE/MTERM
************************************************************************
      SUBROUTINE BNGCHK(MFILE,MSHOW)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNL=16,NNP=6)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /BNTRF/  DTU(9,9,4),DSP(3,3,4),ITRF
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      SAVE
C
      CALL BNTMAT
      CALL BNGTM(2)
      EU=0D0
      EP=0D0
      DO 10 KC=1,4
       DO 5 I=1,9
       DO 5 J=1,9
5      EU=MAX(EU,ABS(GTU(I,J,KC)-DTU(I,J,KC)))
       DO 8 I=1,3
       DO 8 J=1,3
8      EP=MAX(EP,ABS(GSP(I,J,KC)-DSP(I,J,KC)))
10    CONTINUE
      IF (MSHOW.GE.2) WRITE(MTERM,1000) EU,EP
      IF (MSHOW.GE.0) WRITE(MFILE,1000) EU,EP
1000  FORMAT(' BNGCHK: generic vs Q2 transfer: max|GTU-DTU|',D11.3,
     *       '  max|GSP-DSP|',D11.3)
      END
C
************************************************************************
*   BNPRLG/BNRSTG: generic velocity prolongation / restriction
*     (BNPRLU/BNRSTU for K=3), owner arrays KLBOWN from BNOWN
************************************************************************
      SUBROUTINE BNPRLG(DUC,DUF,ILEVF,K)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNL=16,NNP=6)
      DIMENSION DUC(*),DUF(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION KWORK(1),KDF(NNL)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      CALL BNGTM(K)
      NL=(K+1)**2
      ILEVC=ILEVF-1
      KOWN=L(KLBOWN(ILEVF))-1
      KV1=L(KLVERT(ILEVC))
      KM1=L(KLMID(ILEVC))
      NVT1=KNVT(ILEVC)
      NMT1=KNMT(ILEVC)
      IELO=0
      DO 10 J=1,KNU(ILEVF)
       IEL1=KWORK(KOWN+3*J-2)
       KC  =KWORK(KOWN+3*J-1)
       LF  =KWORK(KOWN+3*J)
       IF (IEL1.NE.IELO) THEN
        CALL BNQKDF(K,IEL1,KWORK(KV1),KWORK(KM1),NVT1,NMT1,KDF)
        IELO=IEL1
       ENDIF
       S=0D0
       DO 5 MC=1,NL
5      S=S+GTU(LF,MC,KC)*DUC(KDF(MC))
       DUF(J)=S
10    CONTINUE
      END
C
      SUBROUTINE BNRSTG(DDF,DDC,ILEVC,K)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNL=16,NNP=6)
      DIMENSION DDF(*),DDC(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION KWORK(1),KDF(NNL)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      CALL BNGTM(K)
      NL=(K+1)**2
      ILEVF=ILEVC+1
      KOWN=L(KLBOWN(ILEVF))-1
      KV1=L(KLVERT(ILEVC))
      KM1=L(KLMID(ILEVC))
      NVT1=KNVT(ILEVC)
      NMT1=KNMT(ILEVC)
      DO 1 I=1,KNU(ILEVC)
1     DDC(I)=0D0
      IELO=0
      DO 10 J=1,KNU(ILEVF)
       IEL1=KWORK(KOWN+3*J-2)
       KC  =KWORK(KOWN+3*J-1)
       LF  =KWORK(KOWN+3*J)
       IF (IEL1.NE.IELO) THEN
        CALL BNQKDF(K,IEL1,KWORK(KV1),KWORK(KM1),NVT1,NMT1,KDF)
        IELO=IEL1
       ENDIF
       D=DDF(J)
       DO 5 MC=1,NL
5      DDC(KDF(MC))=DDC(KDF(MC))+GTU(LF,MC,KC)*D
10    CONTINUE
      END
C
************************************************************************
*   BNPRPG/BNRSPG: generic P(k-1)disc prolongation / restriction
************************************************************************
      SUBROUTINE BNPRPG(DPC,DPF,ILEVF,K)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNL=16,NNP=6)
      DIMENSION DPC(*),DPF(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION KWORK(1),KCH(4)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      CALL BNGTM(K)
      NPL=K*(K+1)/2
      KA2=L(KLADJ(ILEVF))-1
      DO 20 IEL1=1,KNEL(ILEVF-1)
       KCH(1)=IEL1
       DO 5 KC=2,4
5      KCH(KC)=KWORK(KA2+4*(KCH(KC-1)-1)+2)
       IC0=NPL*(IEL1-1)
       DO 10 KC=1,4
        IF0=NPL*(KCH(KC)-1)
        DO 10 LF=1,NPL
         S=0D0
         DO 8 IQ=1,NPL
8        S=S+GSP(LF,IQ,KC)*DPC(IC0+IQ)
         DPF(IF0+LF)=S
10     CONTINUE
20    CONTINUE
      END
C
      SUBROUTINE BNRSPG(DDF,DDC,ILEVC,K)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNL=16,NNP=6)
      DIMENSION DDF(*),DDC(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION KWORK(1),KCH(4)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      CALL BNGTM(K)
      NPL=K*(K+1)/2
      KA2=L(KLADJ(ILEVC+1))-1
      DO 20 IEL1=1,KNEL(ILEVC)
       KCH(1)=IEL1
       DO 5 KC=2,4
5      KCH(KC)=KWORK(KA2+4*(KCH(KC-1)-1)+2)
       IC0=NPL*(IEL1-1)
       DO 8 IQ=1,NPL
8      DDC(IC0+IQ)=0D0
       DO 10 KC=1,4
        IF0=NPL*(KCH(KC)-1)
        DO 10 IQ=1,NPL
        DO 10 LF=1,NPL
10      DDC(IC0+IQ)=DDC(IC0+IQ)+GSP(LF,IQ,KC)*DDF(IF0+LF)
20    CONTINUE
      END
C
************************************************************************
*   BNRSLG: generic solution restriction (BNRSOL for K=3): coarse Qk
*     interpolation of the fine solution. Coarse Lagrange node values
*     from the child that contains the node (GWN, KNC), converted to
*     Bernstein coefficients with GLI x GLI. Shared coarse dofs get the
*     same value from every element (continuity).
************************************************************************
      SUBROUTINE BNRSLG(DU1C,DU2C,DU1F,DU2F,ILEVC,K)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNL=16,NNP=6)
      DIMENSION DU1C(*),DU2C(*),DU1F(*),DU2F(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION KWORK(1),KDF(NNL),KDFF(NNL,4),KCH(4),KI(NNL),KJ(NNL)
C *** V* is not in the implicit DOUBLE PRECISION ranges: declare
      DOUBLE PRECISION V1(0:3,0:3),V2(0:3,0:3)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      CALL BNGTM(K)
      CALL BNLORD(K,KI,KJ)
      NL=(K+1)**2
      ILEVF=ILEVC+1
      KV1=L(KLVERT(ILEVC))
      KM1=L(KLMID(ILEVC))
      KV2=L(KLVERT(ILEVF))
      KM2=L(KLMID(ILEVF))
      KA2=L(KLADJ(ILEVF))-1
      NVT1=KNVT(ILEVC)
      NMT1=KNMT(ILEVC)
      NVT2=KNVT(ILEVF)
      NMT2=KNMT(ILEVF)
      DO 50 IEL1=1,KNEL(ILEVC)
       CALL BNQKDF(K,IEL1,KWORK(KV1),KWORK(KM1),NVT1,NMT1,KDF)
       KCH(1)=IEL1
       DO 5 KC=2,4
5      KCH(KC)=KWORK(KA2+4*(KCH(KC-1)-1)+2)
       DO 6 KC=1,4
6      CALL BNQKDF(K,KCH(KC),KWORK(KV2),KWORK(KM2),NVT2,NMT2,
     *             KDFF(1,KC))
C *** values in the coarse Lagrange nodes
       DO 20 IA=0,K
       DO 20 IB=0,K
        KC=KNC(IA,IB)
        S1=0D0
        S2=0D0
        DO 15 LF=1,NL
         JF=KDFF(LF,KC)
         S1=S1+GWN(LF,IA,IB)*DU1F(JF)
         S2=S2+GWN(LF,IA,IB)*DU2F(JF)
15      CONTINUE
        V1(IA,IB)=S1
        V2(IA,IB)=S2
20     CONTINUE
C *** Lagrange -> Bernstein
       DO 40 MC=1,NL
        I=KI(MC)
        J=KJ(MC)
        S1=0D0
        S2=0D0
        DO 30 IA=0,K
        DO 30 IB=0,K
         W=GLI(I,IA)*GLI(J,IB)
         S1=S1+W*V1(IA,IB)
         S2=S2+W*V2(IA,IB)
30      CONTINUE
        DU1C(KDF(MC))=S1
        DU2C(KDF(MC))=S2
40     CONTINUE
50    CONTINUE
      END
C
************************************************************************
*   BNPERP: L2 error of the P(k-1)disc pressure on the current level
*     (ELPQP for IELT>=4), both means removed (BNPERR, bnstok.f):
*     |p_h - p|, |p| (both mean-free) and the relative error to MFILE
************************************************************************
      SUBROUTINE BNPERP(DP,PE,MFILE)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNWORK=1)
      DIMENSION DP(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      EXTERNAL PE
      SAVE
C
      CALL ZNEW(NPL*NEL,3,LDFP,'KDFPER')
      CALL ZNEW(NPL*NEL,1,LZER,'DZERPE')
      IF (IER.NE.0) RETURN
      CALL BNDFPP(NPL,NEL,KWORK(L(LDFP)))
      CALL BNPERR(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LDFP)),
     *            DP,PE,ERR)
      CALL BNPERR(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LDFP)),
     *            DWORK(L(LZER)),PE,PNRM)
      CALL ZDISP(0,LZER,'DZERPE')
      CALL ZDISP(0,LDFP,'KDFPER')
      PREL=0D0
      IF (PNRM.GT.0D0) PREL=ERR/PNRM
      WRITE(MTERM,1000) KORD-1,ERR,PNRM,PREL
      WRITE(MFILE,1000) KORD-1,ERR,PNRM,PREL
1000  FORMAT(' *BNPERP* P',I1,'disc: L2 |p-ph|',D14.6,'  |p|',D14.6,
     *       '  rel.',D14.6,' (means removed)')
      END
