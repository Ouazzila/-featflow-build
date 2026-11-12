************************************************************************
*   bniso.f: isoparametric (curved) boundary elements for the Bernstein
*   pairs in pp2d (IELT=4,5), switch IISO=1 (optional pp2d.dat line)
*
*   Geometry of degree m = k (Q2 map for Q2/P1disc, Q3 map for Q3/P2disc)
*   on the elements that have an edge on a curved boundary; all other
*   elements keep the bilinear map (unchanged code path).
*   Geometry nodes: Lagrange nodes (xi_a, eta_b) = -1+2a/m, -1+2b/m.
*   Boundary edge from vertex V1 (parameter t1) with midpoint parameter
*   tm (FEAT DVBDP/DMBDP): node at s (0<s<1, from V1) on the curve at
*   t = t1 + 2 s (tm - t1) (t1 shifted by TMAX across the start of the
*   boundary component). Deviation from the chord d(s) = P(t) - chord(s).
*   Interior nodes by Gordon-Hall blending: x = bilinear + sum over the
*   curved edges of d(s) * (linear blend towards the opposite edge).
*   The map is stored as Bernstein Qm coefficients (lexicographic
*   l = i + j(m+1) + 1, as bern.f), obtained from the node positions
*   with the 1D Lagrange->Bernstein matrix GLI (BNGTM, bnq3.f).
*
*   COMMON /BNISO/ IISO      : 0 = off, 1 = isoparametric boundary
*                  MGEO      : degree m of the geometry (= k)
*                  KLCRV(L)  : per element: 0 = bilinear, ic > 0 = curved
*                              element number ic in KLGEO
*                  KLGEO(L)  : 2*(m+1)^2 coefficients per curved element
*                  NCRV(L)   : number of curved elements
*   COMMON /BNCUR/ geometry of the current element for the own assembly
*                  (BNASM, BNMPI, BNPERR via BNEVQK/BNMAP): BNGSET(IEL)
*
*   Users of the map: BNMAP (bern.f: M, ST, B, Mp, pressure error), and
*   BNGJAC in the FEAT-type loops SUPWDG (convection), ELPQU (velocity
*   error), BNFVOL (forces). Not changed: right hand side f (FEAT XVB0,
*   bilinear; f=0 in DFG and Kovasznay), jump term (IJUMP=0 with IISO).
************************************************************************
*
************************************************************************
*   BNISOL: curved element geometry of level ILEVI (current level set
*     by SETLEV, /TRIAD/ and /TRIAA/ of ILEVI)
************************************************************************
      SUBROUTINE BNISOL(ILEVI,MFILE,MSHOW)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1,NNL=16,NNP=6)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /BNISO/  IISO,MGEO,KLCRV(NNLEV),KLGEO(NNLEV),NCRV(NNLEV)
      COMMON /BNTRG/  GTU(NNL,NNL,4),GSP(NNP,NNP,4),GLI(0:3,0:3),
     *                GWN(NNL,0:3,0:3),KNC(0:3,0:3),KGT
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
C *** V* is not in the implicit DOUBLE PRECISION ranges: declare
      DOUBLE PRECISION VX(0:3,0:3),VY(0:3,0:3)
      DIMENSION XV(4),YV(4)
      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      SAVE
C
      MG=MGEO
      NLG=(MG+1)**2
      CALL BNGTM(MG)
      CALL ZNEW(NEL,3,KLCRV(ILEVI),'KCRVIS')
C *** deviations d(s_is) of the curved edges: DEV(2,3,4,NEL)
      CALL ZNEW(24*NEL,1,LDEV,'DDEVIS')
      IF (IER.NE.0) STOP 'BNISOL: not enough workspace'
      KCRV=L(KLCRV(ILEVI))-1
      KDEV=L(LDEV)-1
      KCV =L(LCORVG)-1
      KVT =L(LVERT)-1
      DEVMX=0D0
      DO 30 IVBD=1,NVBD
       CALL GETMBD(IMID,IV1,IV2,IVBD,KWORK(L(LVBD)),KWORK(L(LEBD)),
     *             KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),INPR)
       IF (IMID.EQ.0) GOTO 30
       IEL=KWORK(L(LEBD)+IVBD-1)
       II=0
       DO 5 I=1,4
5      IF (KWORK(KVT+4*(IEL-1)+I).EQ.IV1) II=I
       IF (II.EQ.0) STOP 'BNISOL: boundary vertex not in element'
       T1=DWORK(L(LVBDP)+IVBD-1)
       TM=DWORK(L(LMBDP)+IVBD-1)
       TMX=TMAX(INPR)
       IF (TM-T1.GT.0.5D0*TMX) T1=T1+TMX
       IF (T1-TM.GT.0.5D0*TMX) T1=T1-TMX
       X1=DWORK(KCV+2*IV1-1)
       Y1=DWORK(KCV+2*IV1)
       X2=DWORK(KCV+2*IV2-1)
       Y2=DWORK(KCV+2*IV2)
       HE=SQRT((X2-X1)**2+(Y2-Y1)**2)
       DMX=0D0
       DO 10 IS=1,MG-1
        S=DBLE(IS)/DBLE(MG)
        T=T1+2D0*S*(TM-T1)
        T=MOD(T,TMX)
        IF (T.LT.0D0) T=T+TMX
        DX=PARX(T,INPR)-(X1+S*(X2-X1))
        DY=PARY(T,INPR)-(Y1+S*(Y2-Y1))
        IP=KDEV+24*(IEL-1)+6*(II-1)+2*(IS-1)
        DWORK(IP+1)=DX
        DWORK(IP+2)=DY
        DMX=MAX(DMX,SQRT(DX*DX+DY*DY))
10     CONTINUE
       IF (DMX.GT.1D-12*HE) THEN
        KWORK(KCRV+IEL)=-1
        DEVMX=MAX(DEVMX,DMX/HE)
       ELSE
        DO 12 IS=1,6
12      DWORK(KDEV+24*(IEL-1)+6*(II-1)+IS)=0D0
       ENDIF
30    CONTINUE
C
C *** number the curved elements, allocate the coefficients
      NC=0
      DO 40 IEL=1,NEL
       IF (KWORK(KCRV+IEL).NE.0) THEN
        NC=NC+1
        KWORK(KCRV+IEL)=NC
       ENDIF
40    CONTINUE
      NCRV(ILEVI)=NC
      KLGEO(ILEVI)=0
      IF (NC.GT.0) THEN
       CALL ZNEW(2*NLG*NC,1,KLGEO(ILEVI),'DGEOIS')
       IF (IER.NE.0) STOP 'BNISOL: not enough workspace'
      ENDIF
C
      DETMN=1D30
      DO 100 IEL=1,NEL
       IC=KWORK(KCRV+IEL)
       IF (IC.EQ.0) GOTO 100
       DO 45 IV=1,4
        JV=KWORK(KVT+4*(IEL-1)+IV)
        XV(IV)=DWORK(KCV+2*JV-1)
        YV(IV)=DWORK(KCV+2*JV)
45     CONTINUE
       KD=KDEV+24*(IEL-1)
C *** node positions: bilinear + blended edge deviations
       DO 60 IA=0,MG
       DO 60 IB=0,MG
        XI =-1D0+2D0*DBLE(IA)/DBLE(MG)
        ETA=-1D0+2D0*DBLE(IB)/DBLE(MG)
        F1=0.25D0*(1D0-XI)*(1D0-ETA)
        F2=0.25D0*(1D0+XI)*(1D0-ETA)
        F3=0.25D0*(1D0+XI)*(1D0+ETA)
        F4=0.25D0*(1D0-XI)*(1D0+ETA)
        X=F1*XV(1)+F2*XV(2)+F3*XV(3)+F4*XV(4)
        Y=F1*YV(1)+F2*YV(2)+F3*YV(3)+F4*YV(4)
C ***   edge 1 (eta=-1, v1->v2, s=a/m), edge 2 (xi=1, v2->v3, s=b/m),
C       edge 3 (eta=1, v3->v4, s=(m-a)/m), edge 4 (xi=-1, v4->v1,
C       s=(m-b)/m)
        IF ((IA.GT.0).AND.(IA.LT.MG)) THEN
         W=0.5D0*(1D0-ETA)
         X=X+W*DWORK(KD+0 +2*(IA-1)+1)
         Y=Y+W*DWORK(KD+0 +2*(IA-1)+2)
         W=0.5D0*(1D0+ETA)
         X=X+W*DWORK(KD+12+2*(MG-IA-1)+1)
         Y=Y+W*DWORK(KD+12+2*(MG-IA-1)+2)
        ENDIF
        IF ((IB.GT.0).AND.(IB.LT.MG)) THEN
         W=0.5D0*(1D0+XI)
         X=X+W*DWORK(KD+6 +2*(IB-1)+1)
         Y=Y+W*DWORK(KD+6 +2*(IB-1)+2)
         W=0.5D0*(1D0-XI)
         X=X+W*DWORK(KD+18+2*(MG-IB-1)+1)
         Y=Y+W*DWORK(KD+18+2*(MG-IB-1)+2)
        ENDIF
        VX(IA,IB)=X
        VY(IA,IB)=Y
60     CONTINUE
C *** Lagrange values -> Bernstein coefficients
       KG=L(KLGEO(ILEVI))+2*NLG*(IC-1)-1
       DO 80 J=0,MG
       DO 80 I=0,MG
        SX=0D0
        SY=0D0
        DO 70 IA=0,MG
        DO 70 IB=0,MG
         W=GLI(I,IA)*GLI(J,IB)
         SX=SX+W*VX(IA,IB)
         SY=SY+W*VY(IA,IB)
70      CONTINUE
        LC=I+J*(MG+1)+1
        DWORK(KG+LC)=SX
        DWORK(KG+NLG+LC)=SY
80     CONTINUE
C *** check: det J > 0 in the 6x6 Gauss points
       DO 90 IA=1,6
       DO 90 IB=1,6
        XI =-1D0+2D0*(DBLE(IA)-0.5D0)/6D0
        ETA=-1D0+2D0*(DBLE(IB)-0.5D0)/6D0
        CALL BNCMAP(MG,DWORK(KG+1),DWORK(KG+NLG+1),XI,ETA,X,Y,
     *              D11,D12,D21,D22,DET)
        DETMN=MIN(DETMN,DET)
90     CONTINUE
100   CONTINUE
      CALL ZDISP(0,LDEV,'DDEVIS')
C
      IF (DETMN.LE.0D0) THEN
       WRITE(MTERM,*) 'BNISOL: non-positive det J of a curved element'
       STOP
      ENDIF
      IF (MSHOW.GE.2) WRITE(MTERM,1000) ILEVI,MG,NC,DEVMX,DETMN
      IF (MSHOW.GE.0) WRITE(MFILE,1000) ILEVI,MG,NC,DEVMX,DETMN
1000  FORMAT(' BNISOL: level',I2,' geometry degree',I2,' curved el.',
     *       I7,' max |d|/h_E',D11.3,' min det J',D11.3)
      END
C
************************************************************************
*   BNCMAP: map of degree MG with Bernstein coefficients CX, CY,
*     point (X,Y) and Jacobian (DJ11=dx/dxi, DJ12=dx/deta, DJ21=dy/dxi,
*     DJ22=dy/deta), DET = det J
************************************************************************
      SUBROUTINE BNCMAP(MG,CX,CY,XI,ETA,X,Y,DJ11,DJ12,DJ21,DJ22,DET)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBK=8)
      DIMENSION CX(*),CY(*)
      DIMENSION DV1(0:NNBK),D1(0:NNBK),DV2(0:NNBK),D2(0:NNBK)
C
      CALL BNB1D(MG,XI ,DV1,D1)
      CALL BNB1D(MG,ETA,DV2,D2)
      X=0D0
      Y=0D0
      DJ11=0D0
      DJ12=0D0
      DJ21=0D0
      DJ22=0D0
      DO 10 J=0,MG
      DO 10 I=0,MG
       LC=I+J*(MG+1)+1
       F =DV1(I)*DV2(J)
       FX=D1(I)*DV2(J)
       FE=DV1(I)*D2(J)
       X=X+F*CX(LC)
       Y=Y+F*CY(LC)
       DJ11=DJ11+FX*CX(LC)
       DJ12=DJ12+FE*CX(LC)
       DJ21=DJ21+FX*CY(LC)
       DJ22=DJ22+FE*CY(LC)
10    CONTINUE
      DET=DJ11*DJ22-DJ12*DJ21
      END
C
************************************************************************
*   BNGSET: geometry of element IEL of the current level ILEV into
*     /BNCUR/ for BNMAP (own assembly); IEL=0 or IISO=0 -> bilinear
************************************************************************
      SUBROUTINE BNGSET(IEL)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /BNISO/  IISO,MGEO,KLCRV(NNLEV),KLGEO(NNLEV),NCRV(NNLEV)
      COMMON /BNCUR/  ICUR,MGC,CGX(16),CGY(16)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      ICUR=0
      IF ((IISO.LE.0).OR.(IEL.LE.0)) RETURN
      IF (KLCRV(ILEV).EQ.0) RETURN
      IC=KWORK(L(KLCRV(ILEV))+IEL-1)
      IF (IC.LE.0) RETURN
      MGC=MGEO
      NLG=(MGEO+1)**2
      KG=L(KLGEO(ILEV))+2*NLG*(IC-1)-1
      DO 10 LC=1,NLG
       CGX(LC)=DWORK(KG+LC)
       CGY(LC)=DWORK(KG+NLG+LC)
10    CONTINUE
      ICUR=1
      END
C
************************************************************************
*   BNGJAC: FEAT-type loops (SUPWDG, ELPQU, BNFVOL): for a curved
*     element IEL of the current level replace DJAC, DETJ, the weight
*     OM = W*DETJ and the physical point (XX,YY); BCRV = .TRUE. then.
*     Straight elements: nothing changed, BCRV = .FALSE.
************************************************************************
      SUBROUTINE BNGJAC(IEL,XI1,XI2,W,DJAC,DETJ,OM,XX,YY,BCRV)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION DJAC(2,2)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /BNISO/  IISO,MGEO,KLCRV(NNLEV),KLGEO(NNLEV),NCRV(NNLEV)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      BCRV=.FALSE.
      IF (IISO.LE.0) RETURN
      IF (KLCRV(ILEV).EQ.0) RETURN
      IC=KWORK(L(KLCRV(ILEV))+IEL-1)
      IF (IC.LE.0) RETURN
      NLG=(MGEO+1)**2
      KG=L(KLGEO(ILEV))+2*NLG*(IC-1)
      CALL BNCMAP(MGEO,DWORK(KG),DWORK(KG+NLG),XI1,XI2,XX,YY,
     *            D11,D12,D21,D22,DETJ)
      DJAC(1,1)=D11
      DJAC(1,2)=D12
      DJAC(2,1)=D21
      DJAC(2,2)=D22
      OM=W*DETJ
      BCRV=.TRUE.
      END
