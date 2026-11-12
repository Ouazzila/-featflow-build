************************************************************************
*   bnpp.f: Bernstein Q2/P1disc inside pp2d (stage 3, IELT=4)
*   Q3/P2disc (IELT=5): k-dependent parts branch to bnq3.f (IBNPP=3)
*
*   Velocity: Bernstein Q2 as a FEAT element (EB20), numbered by NDFGL
*   like E013 (type 13): vertices, KMID (=NVT+IMT), NVT+NMT+IEL.
*   Pressure: P1disc, Bernstein basis w.r.t. the reference triangle
*   (-1,-1),(3,-1),(-1,3) (as BNPB with N=1), 3 dofs per element,
*   global index 3*(IEL-1)+l.
************************************************************************
*
************************************************************************
*   EB20: Bernstein Q2 on the reference square [-1,1]^2, parametric
*         (bilinear map), FEAT /ELEM/ interface as E013:
*     IPAR = -1 : return element type 13 (Q2 numbering in NDFGL)
*     IPAR = -2 : evaluate on all cubature points of /CUB/
*     IPAR = -3 : use the values of cubature point ICUBP
*     IPAR =  0 : evaluate in (XI1,XI2)
*   Local order (E013): v1..v4 = (-1,-1),(1,-1),(1,1),(-1,1),
*   m1..m4 = edge v_i -> v_(i+1), 9 = centre.
*   phi = B_i(xi) B_j(eta), B_0=(1-x)^2/4, B_1=(1-x^2)/2, B_2=(1+x)^2/4
************************************************************************
      SUBROUTINE EB20(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=9)
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
      SUB='EB20'
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
C *** Evaluate at point (XI1,XI2)
10    ICUBP0=1
      CALL EB20A(XI1,XI2,DHELP,ICUBP0)
      GOTO 4
C
C *** Set number of element (numbering as Q2, E013)
1     IPAR=13
      GOTO 99999
C
C *** Evaluate basis functions on the reference element
C *** for each cubature point
2     DO 11 ICUBP0=1,NCUBP
      CALL EB20A(DXI(ICUBP0,1),DXI(ICUBP0,2),DHELP,ICUBP0)
11    CONTINUE
      GOTO 99999
C
C *** Form linear combinations of the values in DHELP
C *** for the actual element
3     ICUBP0=ICUBP
C
4     IF (BDER(4).OR.BDER(5).OR.BDER(6)) THEN
       CALL WERR(-131,'EB20  ')
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
*   EB20A: values and reference derivatives of the 9 Bernstein Q2
*          basis functions in (X1,X2), stored in DHELP(.,.,ICUBP0)
************************************************************************
      SUBROUTINE EB20A(X1,X2,DHELP,ICUBP0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NDOFL=9)
      DIMENSION DHELP(NDOFL,3,*)
      DIMENSION DB(0:2),DDB(0:2),DC(0:2),DDC(0:2)
      DIMENSION KI(NDOFL),KJ(NDOFL)
C *** (i,j) of B_i(xi) B_j(eta) for the local dofs v1..v4,m1..m4,c
      DATA KI/0,2,2,0, 1,2,1,0, 1/
      DATA KJ/0,0,2,2, 0,1,2,1, 1/
C
      DB(0) =0.25D0*(1D0-X1)**2
      DB(1) =0.5D0 *(1D0-X1*X1)
      DB(2) =0.25D0*(1D0+X1)**2
      DDB(0)=-0.5D0*(1D0-X1)
      DDB(1)=-X1
      DDB(2)= 0.5D0*(1D0+X1)
      DC(0) =0.25D0*(1D0-X2)**2
      DC(1) =0.5D0 *(1D0-X2*X2)
      DC(2) =0.25D0*(1D0+X2)**2
      DDC(0)=-0.5D0*(1D0-X2)
      DDC(1)=-X2
      DDC(2)= 0.5D0*(1D0+X2)
      DO 10 IDFL=1,NDOFL
       I=KI(IDFL)
       J=KJ(IDFL)
       DHELP(IDFL,1,ICUBP0)= DB(I)* DC(J)
       DHELP(IDFL,2,ICUBP0)=DDB(I)* DC(J)
       DHELP(IDFL,3,ICUBP0)= DB(I)*DDC(J)
10    CONTINUE
C
      END
C
************************************************************************
*   COMMON /BNLEV/ (stage 3, IELT=4), all handles per level:
*     IBNPP     : degree k of the active Bernstein pair, 0 = off,
*                 2 = Q2/P1disc (IELT=4), 3 = Q3/P2disc (IELT=5)
*     KLBOWN(L) : owner of each fine velocity dof (L>NLMIN): coarse
*                 element, child position k, local dof l (3*NU)
*     KLBWP(L)  : int psi_q for all pressure dofs (NP)
*     KLBMPI(L) : inverse element pressure mass blocks (9*NEL)
*     KLBBV(L)  : Dirichlet list: vertex 1, vertex 2, INPR (3*NMBD)
*     DBNAR(L)  : |Omega_h|
************************************************************************
*
************************************************************************
*   BNBDLS: Dirichlet dof list for Bernstein Q2 on the current level
*     (after BDRNEU, which marks Neumann edges by KNPR(NVT+IMID)=0):
*     KMBD = all vertex dofs of Dirichlet boundary edges, followed by
*     the edge dofs NVT+IMID; DDBD = parameter of the edge midpoint
*     (0 for vertices). KLMBD/KLDBD are reallocated, NMBD returned.
************************************************************************
      SUBROUTINE BNBDLS(ILEVB,NMBD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      SAVE
C
      CALL ZDISP(0,KLDBD(ILEVB),'DDBD  ')
      CALL ZDISP(0,KLMBD(ILEVB),'KMBD  ')
      KB=MAX(IBNPP,2)
      CALL ZNEW(KB*NVBD,3,KLMBD(ILEVB),'KMBD  ')
      CALL ZNEW(KB*NVBD,1,KLDBD(ILEVB),'DDBD  ')
      CALL ZNEW(3*KB*NVBD,3,KLBBV(ILEVB),'KBBV  ')
      CALL ZNEW(NVT   ,3,LMARK       ,'KMARK ')
      IF (IER.NE.0) STOP 'BNBDLS: not enough workspace'
C
      NMBD=0
      DO 30 IPASS=1,2
      DO 20 IVBD=1,NVBD
       CALL GETMBD(IMID,IV1,IV2,IVBD,KWORK(L(LVBD)),KWORK(L(LEBD)),
     *             KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),INPR)
       IF (IMID.EQ.0) GOTO 20
C ***  Neumann edge (BDRNEU)
       IF (KWORK(L(LNPR)+NVT+IMID-1).EQ.0) GOTO 20
       IF (IPASS.EQ.1) THEN
C ***   pass 1: vertices of Dirichlet edges, each once
        DO 10 IV=1,2
         IF (IV.EQ.1) IVT=IV1
         IF (IV.EQ.2) IVT=IV2
         IF (KWORK(L(LMARK)+IVT-1).NE.0) GOTO 10
         KWORK(L(LMARK)+IVT-1)=1
         NMBD=NMBD+1
         KWORK(L(KLMBD(ILEVB))+NMBD-1)=IVT
         DWORK(L(KLDBD(ILEVB))+NMBD-1)=0D0
         KWORK(L(KLBBV(ILEVB))+3*NMBD-3)=IVT
         KWORK(L(KLBBV(ILEVB))+3*NMBD-2)=0
         KWORK(L(KLBBV(ILEVB))+3*NMBD-1)=INPR
10      CONTINUE
       ELSE IF (IBNPP.EQ.3) THEN
C ***   pass 2, Q3: two edge dofs, NVT+2(IMID-1)+1 next to the lower
C       vertex number; KBBV = (vertex next to the dof, other vertex,
C       INPR), values by BNBSET from the straight edge points
        IVL=MIN(IV1,IV2)
        IVH=MAX(IV1,IV2)
C ***   DDBD = boundary parameter of the curve point at 1/3 of the edge
C       from the vertex next to the dof (used with IISO>0, BNBSET)
        T1=DWORK(L(LVBDP)+IVBD-1)
        TM=DWORK(L(LMBDP)+IVBD-1)
        TMX=TMAX(INPR)
        IF (TM-T1.GT.0.5D0*TMX) T1=T1+TMX
        IF (T1-TM.GT.0.5D0*TMX) T1=T1-TMX
        DO 15 IS=1,2
         NMBD=NMBD+1
         KWORK(L(KLMBD(ILEVB))+NMBD-1)=NVT+2*(IMID-1)+IS
         S=1D0/3D0
         IF (((IS.EQ.1).AND.(IVL.NE.IV1)).OR.
     *       ((IS.EQ.2).AND.(IVH.NE.IV1))) S=2D0/3D0
         T=MOD(T1+2D0*S*(TM-T1),TMX)
         IF (T.LT.0D0) T=T+TMX
         DWORK(L(KLDBD(ILEVB))+NMBD-1)=T
         IF (IS.EQ.1) THEN
          KWORK(L(KLBBV(ILEVB))+3*NMBD-3)=IVL
          KWORK(L(KLBBV(ILEVB))+3*NMBD-2)=IVH
         ELSE
          KWORK(L(KLBBV(ILEVB))+3*NMBD-3)=IVH
          KWORK(L(KLBBV(ILEVB))+3*NMBD-2)=IVL
         ENDIF
         KWORK(L(KLBBV(ILEVB))+3*NMBD-1)=INPR
15      CONTINUE
       ELSE
C ***   pass 2: edge dofs
        NMBD=NMBD+1
        KWORK(L(KLMBD(ILEVB))+NMBD-1)=NVT+IMID
        DWORK(L(KLDBD(ILEVB))+NMBD-1)=DWORK(L(LMBDP)+IVBD-1)
        KWORK(L(KLBBV(ILEVB))+3*NMBD-3)=IV1
        KWORK(L(KLBBV(ILEVB))+3*NMBD-2)=IV2
        KWORK(L(KLBBV(ILEVB))+3*NMBD-1)=INPR
       ENDIF
20    CONTINUE
30    CONTINUE
      CALL ZDISP(0,LMARK,'KMARK ')
C
      END
C
************************************************************************
*   BNBSET: BDRSET for Bernstein Q2 (current level ILEV, list from
*     BNBDLS): vertex dof = UE(vertex),
*     edge dof = 2*UE(midpoint) - (UE(v1)+UE(v2))/2
*     (Q2 interpolation written in the Bernstein basis)
*     Q3 (IBNPP=3): edge dof next to vertex a (other vertex b):
*     c = (-5 v0 + 18 v1 - 9 v2 + 2 v3)/6, v0=UE(a), v3=UE(b),
*     v1, v2 = UE at 1/3, 2/3 of the straight edge from a to b
*     (cubic interpolation in the Bernstein basis)
************************************************************************
      SUBROUTINE BNBSET(DU1,DU2,DF1,DF2,KMBD,DDBD,NMBD,
     *                  PARX,PARY,UE,ICOMP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION DU1(*),DU2(*),DF1(*),DF2(*),KMBD(*),DDBD(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      COMMON /BNISO/  IISO,MGEO,KLCRV(NNLEV),KLGEO(NNLEV),NCRV(NNLEV)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      EXTERNAL PARX,PARY,UE
      SAVE
C
      KBV=L(KLBBV(ILEV))-1
      KCV=L(LCORVG)-1
      DO 10 IMBD=1,NMBD
       IDOF=KMBD(IMBD)
       IV1 =KWORK(KBV+3*IMBD-2)
       IV2 =KWORK(KBV+3*IMBD-1)
       INPR=KWORK(KBV+3*IMBD)
       X1=DWORK(KCV+2*IV1-1)
       Y1=DWORK(KCV+2*IV1)
       DO 5 IC=1,ICOMP
        EV1=UE(X1,Y1,IC)
        IF (IV2.EQ.0) THEN
         U=EV1
        ELSE
         X2=DWORK(KCV+2*IV2-1)
         Y2=DWORK(KCV+2*IV2)
         EV2=UE(X2,Y2,IC)
         IF ((IBNPP.EQ.3).AND.(IISO.GT.0)) THEN
C ***     curve points: own parameter (1/3 from IV1) and the one of the
C         partner dof of the same edge (1/3 from IV2)
          IPART=IMBD+1
          IF (IV1.GT.IV2) IPART=IMBD-1
          EVA=UE(PARX(DDBD(IMBD),INPR),PARY(DDBD(IMBD),INPR),IC)
          EVB=UE(PARX(DDBD(IPART),INPR),PARY(DDBD(IPART),INPR),IC)
          U=(-5D0*EV1+18D0*EVA-9D0*EVB+2D0*EV2)/6D0
         ELSE IF (IBNPP.EQ.3) THEN
          EVA=UE(X1+(X2-X1)/3D0,Y1+(Y2-Y1)/3D0,IC)
          EVB=UE(X1+2D0*(X2-X1)/3D0,Y1+2D0*(Y2-Y1)/3D0,IC)
          U=(-5D0*EV1+18D0*EVA-9D0*EVB+2D0*EV2)/6D0
         ELSE
          PX=PARX(DDBD(IMBD),INPR)
          PY=PARY(DDBD(IMBD),INPR)
          EVM=UE(PX,PY,IC)
          U=2D0*EVM-0.5D0*(EV1+EV2)
         ENDIF
        ENDIF
        IF (IC.EQ.1) THEN
         DU1(IDOF)=U
         DF1(IDOF)=U
        ELSE
         DU2(IDOF)=U
         DF2(IDOF)=U
        ENDIF
5      CONTINUE
10    CONTINUE
C
      END
C
************************************************************************
*   BNMATS: matrices of Bernstein Q2/P1disc on all levels, into the
*     pp2d handles (/MGFLD/, /MGDIM/, /MGPROJ/); the A structure
*     (KLCOLA, KLLDA, KNA) must exist (XMAP7 with EB20).
*     M, ST (=NY*Laplace), B1, B2 by BNASM (4x4 Gauss), Ml = row sums
*     of M, C = B^T Ml^-1 B over the non-Dirichlet velocity dofs,
*     pressure mass data (BNMPI), MG owner arrays (BNOWN).
*     IMASS=0: M is not kept. IPRECA=0: ST, M stored REAL*4.
************************************************************************
      SUBROUTINE BNMATS(MFILE,MSHOW)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      COMMON /BNISO/  IISO,MGEO,KLCRV(NNLEV),KLGEO(NNLEV),NCRV(NNLEV)
      INTEGER  VIPARM
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C
      IFAM=0
      KORD=MAX(IBNPP,2)
      NLDF=(KORD+1)**2
      NPL =KORD*(KORD+1)/2
      NG  =KORD+2
C
      DO 100 ILEV=NLMIN,NLMAX
       CALL SETLEV(1)
C *** isoparametric boundary elements (bniso.f), before BNASM
       IF (IISO.GT.0) CALL BNISOL(ILEV,MFILE,MSHOW)
       NU=KNU(ILEV)
       NP=KNP(ILEV)
       NA=KNA(ILEV)
C *** element dof maps (Bernstein numbering = FEAT Q2 (type 13) / Q3
C     (type 14, with the EB30 interior order) numbering)
       CALL ZNEW(NLDF*NEL,3,LDFG,'KDFGBN')
       CALL ZNEW(NPL*NEL,3,LDFP,'KDFPBN')
       IF (IER.NE.0) STOP 'BNMATS: not enough workspace'
       CALL BNDFMP(KORD,KWORK(L(LVERT)),KWORK(L(LMID)),NEL,NVT,NMT,
     *             KWORK(L(LDFG)))
       CALL BNDFPP(NPL,NEL,KWORK(L(LDFP)))
C *** B structure (rows = velocity dofs, technique 9)
       CALL BNSTR(NEL,NLDF,KWORK(L(LDFG)),NPL,KWORK(L(LDFP)),NU,
     *            .FALSE.,KLLDB(ILEV),KLCOLB(ILEV),NB)
       KNB(ILEV)=NB
C *** M, ST, B1, B2
       CALL ZNEW(NA,1,KLMASS(ILEV),'DMBN  ')
       CALL ZNEW(NA,1,KLST(ILEV)  ,'DSTBN ')
       CALL ZNEW(NB,1,KLB1(ILEV)  ,'DB1   ')
       CALL ZNEW(NB,1,KLB2(ILEV)  ,'DB2   ')
       CALL ZNEW(NU,1,KLM(ILEV)   ,'DMASS ')
       IF (IER.NE.0) STOP 'BNMATS: not enough workspace'
       CALL BNASM(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LDFG)),
     *            KWORK(L(LDFP)),KWORK(L(KLLDA(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDB(ILEV))),
     *            KWORK(L(KLCOLB(ILEV))),DWORK(L(KLMASS(ILEV))),
     *            DWORK(L(KLST(ILEV))),DWORK(L(KLB1(ILEV))),
     *            DWORK(L(KLB2(ILEV))),DAREA)
       DBNAR(ILEV)=DAREA
       CALL LSC1(DWORK(L(KLST(ILEV))),NA,NY)
C *** lumped mass = row sums of the Bernstein mass matrix (>0)
       KLDA=L(KLLDA(ILEV))-1
       DO 20 IU=1,NU
        S=0D0
        DO 15 IA=KWORK(KLDA+IU),KWORK(KLDA+IU+1)-1
15      S=S+DWORK(L(KLMASS(ILEV))+IA-1)
        DWORK(L(KLM(ILEV))+IU-1)=S
20     CONTINUE
C *** pressure mass data
       CALL ZNEW(NP   ,1,KLBWP(ILEV) ,'DWPBN ')
       CALL ZNEW(NPL*NPL*NEL,1,KLBMPI(ILEV),'DMPIBN')
       IF (IER.NE.0) STOP 'BNMATS: not enough workspace'
       CALL BNMPI(NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            DWORK(L(KLBMPI(ILEV))),DWORK(L(KLBWP(ILEV))))
C *** C = B^T Ml^-1 B
       CALL BNCMAT(ILEV,NU,NP)
C *** MG owner arrays for the velocity prolongation
       IF (ILEV.GT.NLMIN) CALL BNOWN(ILEV,MFILE,MSHOW)
C
       CALL ZDISP(0,LDFP,'KDFPBN')
       CALL ZDISP(0,LDFG,'KDFGBN')
C
       IF (IMASS.EQ.0) CALL ZDISP(0,KLMASS(ILEV),'DMBN  ')
       IF (IPRECA.EQ.0) THEN
        CALL ZCTYPE(2,KLST(ILEV),'DSTBN ')
        IF (IMASS.EQ.1) CALL ZCTYPE(2,KLMASS(ILEV),'DMBN  ')
       ENDIF
C
C *** ILU(C) for the pressure smoother/solver, as in INIT1
       IF ((ISMP.EQ.4).OR.(ISLP.EQ.3).OR.(ISLP.EQ.4)) THEN
        CALL ZNEW(KNC(ILEV),-1,KLCILU(ILEV),'DDCILU')
        IF (IER.NE.0) STOP 'BNMATS: not enough workspace'
        CALL ZCPY(KLC(ILEV),'CC    ',KLCILU(ILEV),'DDCILU')
        CALL IFD17(DWORK(L(KLCILU(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),NP,1,0D0,1D-12)
        IF (IER.NE.0) STOP 'BNMATS: ILU(C) failed'
        CALL ZCTYPE(2,KLCILU(ILEV),'DDCILU')
       ENDIF
C
       IF ((ILEV.EQ.NLMIN).AND.(MSHOW.GE.0))
     *  WRITE(MFILE,'(A,I2,A,I3,A,I3,A,I2)') ' BNMATS: k=',KORD,
     *  ' velocity dofs/el.',NLDF,' pressure dofs/el.',NPL,' Gauss',NG
       IF (MSHOW.GE.2) WRITE(MTERM,1000) ILEV,NU,NP,NA,NB,KNC(ILEV),
     *                                   KNMBD(ILEV),DAREA
       IF (MSHOW.GE.0) WRITE(MFILE,1000) ILEV,NU,NP,NA,NB,KNC(ILEV),
     *                                   KNMBD(ILEV),DAREA
100   CONTINUE
C
C
C *** k=2: generic transfer matrices (bnq3.f) against the Q2 ones
      IF (KORD.EQ.2) CALL BNGCHK(MFILE,MSHOW)
C *** check of the grid transfer (nested spaces): restriction of
C     int phi_j (= lumped mass) and of int psi_q must give the coarse
C     vectors (exact up to the elements at the curved boundary)
      DO 200 ILEV=NLMIN,NLMAX-1
       NUC=KNU(ILEV)
       NPC=KNP(ILEV)
       CALL ZNEW(NUC,1,LTMP,'DTMPBN')
       IF (IER.NE.0) STOP 'BNMATS: not enough workspace'
       CALL BNRSTU(DWORK(L(KLM(ILEV+1))),DWORK(L(LTMP)),ILEV)
       EU=0D0
       EUM=0D0
       DO 210 I=1,NUC
        EU=MAX(EU,ABS(DWORK(L(LTMP)+I-1)-DWORK(L(KLM(ILEV))+I-1)))
210    EUM=MAX(EUM,ABS(DWORK(L(KLM(ILEV))+I-1)))
       CALL ZDISP(0,LTMP,'DTMPBN')
       CALL ZNEW(NPC,1,LTMP,'DTMPBN')
       CALL BNRSTP(DWORK(L(KLBWP(ILEV+1))),DWORK(L(LTMP)),ILEV)
       EP=0D0
       EPM=0D0
       DO 220 I=1,NPC
        EP=MAX(EP,ABS(DWORK(L(LTMP)+I-1)-DWORK(L(KLBWP(ILEV))+I-1)))
220    EPM=MAX(EPM,ABS(DWORK(L(KLBWP(ILEV))+I-1)))
       CALL ZDISP(0,LTMP,'DTMPBN')
       IF (MSHOW.GE.2) WRITE(MTERM,1010) ILEV,EU/EUM,EP/EPM
       IF (MSHOW.GE.0) WRITE(MFILE,1010) ILEV,EU/EUM,EP/EPM
200   CONTINUE
C
1010  FORMAT(' BNMATS: transfer check level',I2,': |R int phi - int',
     *       ' phi|/max',D11.3,'  |R int psi - int psi|/max',D11.3)
1000  FORMAT(' Bernstein Qk/P(k-1)disc: ILEV=',I2,' NU=',I8,' NP=',I8,
     *       ' NA=',I9,' NB=',I9,' NC=',I9,' NDIR=',I6,' |Omega|=',
     *       F16.12)
      END
C
************************************************************************
*   BNCMAT: C = B^T Ml^-1 B (level ILEVC), sum over the velocity dofs
*     that are not Dirichlet dofs (list KLMBD), CSR with the diagonal
*     first (technique 7) into KLC, KLCOLC, KLLDC, KNC
************************************************************************
      SUBROUTINE BNCMAT(ILEVC,NU,NP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      NB=KNB(ILEVC)
      CALL ZNEW(NU  ,3,LBD  ,'KBDC  ')
      CALL ZNEW(NP+1,3,LLDT ,'KLDT  ')
      CALL ZNEW(NB  ,3,LROWT,'KROWT ')
      CALL ZNEW(NB  ,3,LPOST,'KPOST ')
      CALL ZNEW(NP  ,3,LMRK ,'KMRK  ')
      CALL ZNEW(NP+1,3,KLLDC(ILEVC),'KLDC  ')
      IF (IER.NE.0) STOP 'BNCMAT: not enough workspace'
C *** Dirichlet marker
      DO 5 I=1,KNMBD(ILEVC)
       IDOF=KWORK(L(KLMBD(ILEVC))+I-1)
       IF (IDOF.GT.0) KWORK(L(LBD)+IDOF-1)=1
5     CONTINUE
C *** transpose structure of B: velocity rows of each pressure column
      CALL BNTRSP(NU,NP,KWORK(L(KLLDB(ILEVC))),KWORK(L(KLCOLB(ILEVC))),
     *            KWORK(L(LLDT)),KWORK(L(LROWT)),KWORK(L(LPOST)))
C *** count, allocate, fill
      CALL BNCSTR(NP,KWORK(L(KLLDB(ILEVC))),KWORK(L(KLCOLB(ILEVC))),
     *            KWORK(L(LLDT)),KWORK(L(LROWT)),KWORK(L(LBD)),
     *            KWORK(L(LMRK)),KWORK(L(KLLDC(ILEVC))),KWORK(1),
     *            DWORK(1),DWORK(1),DWORK(1),DWORK(1),KWORK(1),0,NC)
      KNC(ILEVC)=NC
      CALL ZNEW(NC,3,KLCOLC(ILEVC),'KCOLC ')
      CALL ZNEW(NC,1,KLC(ILEVC)   ,'CC    ')
      IF (IER.NE.0) STOP 'BNCMAT: not enough workspace'
      CALL BNCSTR(NP,KWORK(L(KLLDB(ILEVC))),KWORK(L(KLCOLB(ILEVC))),
     *            KWORK(L(LLDT)),KWORK(L(LROWT)),KWORK(L(LBD)),
     *            KWORK(L(LMRK)),KWORK(L(KLLDC(ILEVC))),
     *            KWORK(L(KLCOLC(ILEVC))),DWORK(L(KLC(ILEVC))),
     *            DWORK(L(KLB1(ILEVC))),DWORK(L(KLB2(ILEVC))),
     *            DWORK(L(KLM(ILEVC))),KWORK(L(LPOST)),1,NC)
C *** technique 7 as XMAP7: diagonal first, then ascending columns
C     (needed by the ILU and SOR routines of FEAT)
      CALL BNSRT7(NP,KWORK(L(KLLDC(ILEVC))),KWORK(L(KLCOLC(ILEVC))),
     *            DWORK(L(KLC(ILEVC))))
      CALL ZDISP(0,LMRK ,'KMRK  ')
      CALL ZDISP(0,LPOST,'KPOST ')
      CALL ZDISP(0,LROWT,'KROWT ')
      CALL ZDISP(0,LLDT ,'KLDT  ')
      CALL ZDISP(0,LBD  ,'KBDC  ')
C
      END
C
************************************************************************
*   BNTRSP: transpose of a CSR structure (NU rows, NP columns):
*     column q has the rows KROWT(KLDT(q)..KLDT(q+1)-1), KPOST = the
*     positions of these entries in the original arrays
************************************************************************
      SUBROUTINE BNTRSP(NU,NP,KLD,KCOL,KLDT,KROWT,KPOST)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KLD(*),KCOL(*),KLDT(*),KROWT(*),KPOST(*)
C
      DO 5 IQ=1,NP+1
5     KLDT(IQ)=0
      DO 10 IA=1,KLD(NU+1)-1
10    KLDT(KCOL(IA)+1)=KLDT(KCOL(IA)+1)+1
      KLDT(1)=1
      DO 15 IQ=1,NP
15    KLDT(IQ+1)=KLDT(IQ+1)+KLDT(IQ)
      DO 20 IU=1,NU
      DO 20 IA=KLD(IU),KLD(IU+1)-1
       IQ=KCOL(IA)
       KROWT(KLDT(IQ))=IU
       KPOST(KLDT(IQ))=IA
       KLDT(IQ)=KLDT(IQ)+1
20    CONTINUE
      DO 25 IQ=NP,1,-1
25    KLDT(IQ+1)=KLDT(IQ)
      KLDT(1)=1
C
      END
C
************************************************************************
*   BNCSTR: IFILL=0: row lengths of C -> KLDC, NC (KMRK = marker)
*           IFILL=1: columns (diagonal first) and values of C,
*                    KMRK = position of column r in the current row
************************************************************************
      SUBROUTINE BNCSTR(NP,KLDB,KCOLB,KLDT,KROWT,KBD,KMRK,KLDC,KCOLC,
     *                  DC,DB1,DB2,DML,KPOST,IFILL,NC)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KLDB(*),KCOLB(*),KLDT(*),KROWT(*),KBD(*),KMRK(*)
      DIMENSION KLDC(*),KCOLC(*),DC(*),DB1(*),DB2(*),DML(*),KPOST(*)
C
      DO 5 IQ=1,NP
5     KMRK(IQ)=0
      NC=0
      DO 50 IQ=1,NP
       N=1
       IF (IFILL.EQ.0) THEN
        KMRK(IQ)=IQ
       ELSE
        IROW0=KLDC(IQ)
        KCOLC(IROW0)=IQ
        DC(IROW0)=0D0
        KMRK(IQ)=IROW0
       ENDIF
       DO 40 IT=KLDT(IQ),KLDT(IQ+1)-1
        IU=KROWT(IT)
        IF (KBD(IU).NE.0) GOTO 40
        IF (IFILL.EQ.1) THEN
         IAQ=KPOST(IT)
         F1=DB1(IAQ)/DML(IU)
         F2=DB2(IAQ)/DML(IU)
        ENDIF
        DO 30 IA=KLDB(IU),KLDB(IU+1)-1
         IR=KCOLB(IA)
         IF (IFILL.EQ.0) THEN
          IF (KMRK(IR).NE.IQ) THEN
           KMRK(IR)=IQ
           N=N+1
          ENDIF
         ELSE
C ***     KMRK(IR) points into the current row only if >= IROW0
          IP=KMRK(IR)
          IF (IP.LT.IROW0) THEN
           IP=IROW0+N
           N=N+1
           KMRK(IR)=IP
           KCOLC(IP)=IR
           DC(IP)=0D0
          ENDIF
          DC(IP)=DC(IP)+F1*DB1(IA)+F2*DB2(IA)
         ENDIF
30      CONTINUE
40     CONTINUE
       IF (IFILL.EQ.0) THEN
        KLDC(IQ)=NC+1
        NC=NC+N
       ENDIF
50    CONTINUE
      IF (IFILL.EQ.0) KLDC(NP+1)=NC+1
      IF (IFILL.EQ.1) NC=KLDC(NP+1)-1
C
      END
C
************************************************************************
*   BNTMAT: transfer matrices for one coarse element and its 4
*     children (child k: vertex 1 = coarse vertex k, vertex 2 =
*     midpoint of coarse edge k, vertex 3 = centre, vertex 4 =
*     midpoint of coarse edge k-1, see MP030), in reference coords:
*     DTU(l,m,k): Bernstein Q2 coefficient of fine local dof l (child
*                 k) of the coarse basis function m (EB20 order)
*     DSP(l,q,k): P1disc coefficient l on child k of coarse psi_q
*     Exact (affine relation of the reference coordinates).
************************************************************************
      SUBROUTINE BNTMAT
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      COMMON /BNTRF/  DTU(9,9,4),DSP(3,3,4),ITRF
      DIMENSION CX(4),CY(4),QX(4),QY(4),DHELP(9,3,1)
      DIMENSION DLI(0:2,0:2),SN(0:2),DV(0:2,0:2,9)
      DIMENSION DG(3,3),DGI(3,3),DPV(3),PS(3),PT(3)
      DIMENSION KI(9),KJ(9)
      DATA CX/-1D0,1D0,1D0,-1D0/, CY/-1D0,-1D0,1D0,1D0/
      DATA SN/-1D0,0D0,1D0/
      DATA PS/-1D0,1D0,-1D0/, PT/-1D0,-1D0,1D0/
      DATA KI/0,2,2,0, 1,2,1,0, 1/
      DATA KJ/0,0,2,2, 0,1,2,1, 1/
      SAVE
C
      IF (ITRF.EQ.1) RETURN
C *** Lagrange values at 0,1/2,1 -> Bernstein coefficients (Q2, 1D)
      DO 2 I=0,2
      DO 2 J=0,2
2     DLI(I,J)=0D0
      DLI(0,0)= 1D0
      DLI(1,0)=-0.5D0
      DLI(1,1)= 2D0
      DLI(1,2)=-0.5D0
      DLI(2,2)= 1D0
C
      DO 100 K=1,4
       KP=MOD(K,4)+1
       KM=MOD(K+2,4)+1
       QX(1)=CX(K)
       QY(1)=CY(K)
       QX(2)=0.5D0*(CX(K)+CX(KP))
       QY(2)=0.5D0*(CY(K)+CY(KP))
       QX(3)=0D0
       QY(3)=0D0
       QX(4)=0.5D0*(CX(K)+CX(KM))
       QY(4)=0.5D0*(CY(K)+CY(KM))
C *** velocity: coarse basis values at the 3x3 child nodes
       DO 20 IA=0,2
       DO 20 IB=0,2
        CALL BNCHLD(QX,QY,SN(IA),SN(IB),X,Y)
        CALL EB20A(X,Y,DHELP,1)
        DO 15 M=1,9
15      DV(IA,IB,M)=DHELP(M,1,1)
20     CONTINUE
       DO 30 M=1,9
       DO 30 LF=1,9
        I=KI(LF)
        J=KJ(LF)
        S=0D0
        DO 25 IA=0,2
        DO 25 IB=0,2
25      S=S+DLI(I,IA)*DLI(J,IB)*DV(IA,IB,M)
        DTU(LF,M,K)=S
30     CONTINUE
C *** pressure: values of child and coarse basis at 3 child points
       DO 40 IP=1,3
        CALL BNPB(1,PS(IP),PT(IP),DPV)
        DO 35 LF=1,3
35      DG(IP,LF)=DPV(LF)
40     CONTINUE
       CALL BNINV3(DG,DGI)
       DO 50 IQ=1,3
        DO 45 IP=1,3
         CALL BNCHLD(QX,QY,PS(IP),PT(IP),X,Y)
         CALL BNPB(1,X,Y,DPV)
         DG(IP,1)=DPV(IQ)
45      CONTINUE
        DO 48 LF=1,3
48      DSP(LF,IQ,K)=DGI(LF,1)*DG(1,1)+DGI(LF,2)*DG(2,1)
     *              +DGI(LF,3)*DG(3,1)
50     CONTINUE
100   CONTINUE
      ITRF=1
C
      END
C
************************************************************************
*   BNCHLD: point (S,T) of the child reference square -> parent
*     reference coordinates (X,Y), bilinear through QX,QY
************************************************************************
      SUBROUTINE BNCHLD(QX,QY,S,T,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QX(4),QY(4)
      F1=0.25D0*(1D0-S)*(1D0-T)
      F2=0.25D0*(1D0+S)*(1D0-T)
      F3=0.25D0*(1D0+S)*(1D0+T)
      F4=0.25D0*(1D0-S)*(1D0+T)
      X=F1*QX(1)+F2*QX(2)+F3*QX(3)+F4*QX(4)
      Y=F1*QY(1)+F2*QY(2)+F3*QY(3)+F4*QY(4)
      END
C
************************************************************************
*   BNINV3: inverse of a 3x3 matrix (adjugate / determinant)
************************************************************************
      SUBROUTINE BNINV3(A,AI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(3,3),AI(3,3)
      AI(1,1)= A(2,2)*A(3,3)-A(2,3)*A(3,2)
      AI(1,2)=-(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      AI(1,3)= A(1,2)*A(2,3)-A(1,3)*A(2,2)
      AI(2,1)=-(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      AI(2,2)= A(1,1)*A(3,3)-A(1,3)*A(3,1)
      AI(2,3)=-(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      AI(3,1)= A(2,1)*A(3,2)-A(2,2)*A(3,1)
      AI(3,2)=-(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      AI(3,3)= A(1,1)*A(2,2)-A(1,2)*A(2,1)
      DET=A(1,1)*AI(1,1)+A(1,2)*AI(2,1)+A(1,3)*AI(3,1)
      DO 10 I=1,3
      DO 10 J=1,3
10    AI(I,J)=AI(I,J)/DET
      END
C
************************************************************************
*   BNQ2DF: global Q2 dofs of element IEL (EB20 order), FEAT numbering
************************************************************************
      SUBROUTINE BNQ2DF(IEL,KVERT,KMID,NVT,NMT,KDF)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KVERT(4,*),KMID(4,*),KDF(9)
      DO 10 I=1,4
       KDF(I)  =KVERT(I,IEL)
       KDF(4+I)=KMID(I,IEL)
10    CONTINUE
      KDF(9)=NVT+NMT+IEL
      END
C
************************************************************************
*   BNOWN: owner of every velocity dof of level ILEVF (coarse element,
*     child position, local dof), first come first served; check of
*     the child convention: max distance of the fine vertices from the
*     image of the predicted parent reference points (interior
*     elements: round-off, boundary elements: curved boundary)
************************************************************************
      SUBROUTINE BNOWN(ILEVF,MFILE,MSHOW)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
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
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      DIMENSION KDF(16),CX(4),CY(4),KCH(4)
      DATA CX/-1D0,1D0,1D0,-1D0/, CY/-1D0,-1D0,1D0,1D0/
      SAVE
C
      KB=MAX(IBNPP,2)
      NLB=(KB+1)**2
      IF (KB.EQ.2) THEN
       CALL BNTMAT
      ELSE
       CALL BNGTM(KB)
      ENDIF
      ILEVC=ILEVF-1
      NU2=KNU(ILEVF)
      CALL ZNEW(3*NU2,3,KLBOWN(ILEVF),'KBOWN ')
      IF (IER.NE.0) STOP 'BNOWN: not enough workspace'
      KOWN=L(KLBOWN(ILEVF))-1
      KV1=L(KLVERT(ILEVC))
      KV2=L(KLVERT(ILEVF))
      KM2=L(KLMID(ILEVF))
      KA2=L(KLADJ(ILEVF))-1
      KC1=L(KLCVG(ILEVC))-1
      KC2=L(KLCVG(ILEVF))-1
      NVT2=KNVT(ILEVF)
      NMT2=KNMT(ILEVF)
      KM1=L(KLMID(ILEVC))
      NBADN=0
      DEVMX=0D0
      HMX=0D0
      DO 50 IEL1=1,KNEL(ILEVC)
       KCH(1)=IEL1
       DO 5 K=2,4
5      KCH(K)=KWORK(KA2+4*(KCH(K-1)-1)+2)
       DO 40 K=1,4
        IELH=KCH(K)
        CALL BNQKDF(KB,IELH,KWORK(KV2),KWORK(KM2),NVT2,NMT2,KDF)
        DO 10 LF=1,NLB
         J=KDF(LF)
         IF (KWORK(KOWN+3*J-2).EQ.0) THEN
          KWORK(KOWN+3*J-2)=IEL1
          KWORK(KOWN+3*J-1)=K
          KWORK(KOWN+3*J)  =LF
         ENDIF
10      CONTINUE
C ***   FEAT numbering used by BNRSOL: fine vertices of child k are
C       coarse vertex k, coarse KMID(k) (=NVT1+IMT), NVT1+NMT1+IEL1
        IVF1=KWORK(KV2+4*(IELH-1))
        IVF2=KWORK(KV2+4*(IELH-1)+1)
        IVF3=KWORK(KV2+4*(IELH-1)+2)
        IF ((IVF1.NE.KWORK(KV1+4*(IEL1-1)+K-1)).OR.
     *      (IVF2.NE.KWORK(KM1+4*(IEL1-1)+K-1)).OR.
     *      (IVF3.NE.KNVT(ILEVC)+KNMT(ILEVC)+IEL1)) NBADN=NBADN+1
C ***   check: fine vertices = coarse map of the predicted points
        KP=MOD(K,4)+1
        KM=MOD(K+2,4)+1
        DO 30 IV=1,4
         IF (IV.EQ.1) THEN
          XI=CX(K)
          ETA=CY(K)
         ELSE IF (IV.EQ.2) THEN
          XI=0.5D0*(CX(K)+CX(KP))
          ETA=0.5D0*(CY(K)+CY(KP))
         ELSE IF (IV.EQ.3) THEN
          XI=0D0
          ETA=0D0
         ELSE
          XI=0.5D0*(CX(K)+CX(KM))
          ETA=0.5D0*(CY(K)+CY(KM))
         ENDIF
         X=0D0
         Y=0D0
         DO 20 IC=1,4
          F=0.25D0*(1D0+CX(IC)*XI)*(1D0+CY(IC)*ETA)
          IVC=KWORK(KV1+4*(IEL1-1)+IC-1)
          X=X+F*DWORK(KC1+2*IVC-1)
          Y=Y+F*DWORK(KC1+2*IVC)
          IF (IC.GT.1) THEN
           IVP=KWORK(KV1+4*(IEL1-1)+IC-2)
           HMX=MAX(HMX,SQRT((DWORK(KC1+2*IVC-1)-DWORK(KC1+2*IVP-1))**2
     *                     +(DWORK(KC1+2*IVC)-DWORK(KC1+2*IVP))**2))
          ENDIF
20       CONTINUE
         IVF=KWORK(KV2+4*(IELH-1)+IV-1)
         DEV=SQRT((X-DWORK(KC2+2*IVF-1))**2+(Y-DWORK(KC2+2*IVF))**2)
         DEVMX=MAX(DEVMX,DEV)
30      CONTINUE
40     CONTINUE
50    CONTINUE
C
      DO 60 J=1,NU2
       IF (KWORK(KOWN+3*J-2).EQ.0) STOP 'BNOWN: dof without owner'
60    CONTINUE
      IF (NBADN.GT.0) THEN
       WRITE(MTERM,*) 'BNOWN: fine vertex numbering differs from the',
     *                ' coarse Q2 numbering in',NBADN,' children'
       STOP
      ENDIF
      IF (MSHOW.GE.2) WRITE(MTERM,1000) ILEVF,DEVMX,HMX
      IF (MSHOW.GE.0) WRITE(MFILE,1000) ILEVF,DEVMX,HMX
1000  FORMAT(' BNOWN: level',I2,' max. deviation of the fine vertices',
     *       ' from the parent map',D11.3,' (coarse h_max',D11.3,')')
C
      END
C
************************************************************************
*   BNPRLU: prolongation DUF = P DUC of one velocity component
*           (fine level ILEVF, coarse ILEVF-1), exact for Q2
*   BNRSTU: restriction DDC = P^T DDF (fine level ILEVC+1)
************************************************************************
      SUBROUTINE BNPRLU(DUC,DUF,ILEVF)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
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
      COMMON /BNTRF/  DTU(9,9,4),DSP(3,3,4),ITRF
      DIMENSION KWORK(1),KDF(9)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF (IBNPP.EQ.3) THEN
       CALL BNPRLG(DUC,DUF,ILEVF,3)
       RETURN
      ENDIF
      ILEVC=ILEVF-1
      KOWN=L(KLBOWN(ILEVF))-1
      KV1=L(KLVERT(ILEVC))
      KM1=L(KLMID(ILEVC))
      NVT1=KNVT(ILEVC)
      NMT1=KNMT(ILEVC)
      IELO=0
      DO 10 J=1,KNU(ILEVF)
       IEL1=KWORK(KOWN+3*J-2)
       K   =KWORK(KOWN+3*J-1)
       LF  =KWORK(KOWN+3*J)
       IF (IEL1.NE.IELO) THEN
        CALL BNQ2DF(IEL1,KWORK(KV1),KWORK(KM1),NVT1,NMT1,KDF)
        IELO=IEL1
       ENDIF
       S=0D0
       DO 5 MC=1,9
5      S=S+DTU(LF,MC,K)*DUC(KDF(MC))
       DUF(J)=S
10    CONTINUE
C
      END
C
      SUBROUTINE BNRSTU(DDF,DDC,ILEVC)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
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
      COMMON /BNTRF/  DTU(9,9,4),DSP(3,3,4),ITRF
      DIMENSION KWORK(1),KDF(9)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF (IBNPP.EQ.3) THEN
       CALL BNRSTG(DDF,DDC,ILEVC,3)
       RETURN
      ENDIF
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
       K   =KWORK(KOWN+3*J-1)
       LF  =KWORK(KOWN+3*J)
       IF (IEL1.NE.IELO) THEN
        CALL BNQ2DF(IEL1,KWORK(KV1),KWORK(KM1),NVT1,NMT1,KDF)
        IELO=IEL1
       ENDIF
       D=DDF(J)
       DO 5 MC=1,9
5      DDC(KDF(MC))=DDC(KDF(MC))+DTU(LF,MC,K)*D
10    CONTINUE
C
      END
C
************************************************************************
*   BNPRLP: prolongation of the P1disc pressure (fine level ILEVF)
*   BNRSTP: restriction = transpose (coarse level ILEVC)
************************************************************************
      SUBROUTINE BNPRLP(DPC,DPF,ILEVF)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION DPC(*),DPF(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /BNTRF/  DTU(9,9,4),DSP(3,3,4),ITRF
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      DIMENSION KWORK(1),KCH(4)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF (IBNPP.EQ.3) THEN
       CALL BNPRPG(DPC,DPF,ILEVF,3)
       RETURN
      ENDIF
      CALL BNTMAT
      KA2=L(KLADJ(ILEVF))-1
      DO 20 IEL1=1,KNEL(ILEVF-1)
       KCH(1)=IEL1
       DO 5 K=2,4
5      KCH(K)=KWORK(KA2+4*(KCH(K-1)-1)+2)
       IC0=3*(IEL1-1)
       DO 10 K=1,4
        IF0=3*(KCH(K)-1)
        DO 10 LF=1,3
10      DPF(IF0+LF)=DSP(LF,1,K)*DPC(IC0+1)+DSP(LF,2,K)*DPC(IC0+2)
     *             +DSP(LF,3,K)*DPC(IC0+3)
20    CONTINUE
C
      END
C
      SUBROUTINE BNRSTP(DDF,DDC,ILEVC)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION DDF(*),DDC(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /BNTRF/  DTU(9,9,4),DSP(3,3,4),ITRF
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      DIMENSION KWORK(1),KCH(4)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF (IBNPP.EQ.3) THEN
       CALL BNRSPG(DDF,DDC,ILEVC,3)
       RETURN
      ENDIF
      CALL BNTMAT
      KA2=L(KLADJ(ILEVC+1))-1
      DO 20 IEL1=1,KNEL(ILEVC)
       KCH(1)=IEL1
       DO 5 K=2,4
5      KCH(K)=KWORK(KA2+4*(KCH(K-1)-1)+2)
       IC0=3*(IEL1-1)
       DO 8 IQ=1,3
8      DDC(IC0+IQ)=0D0
       DO 10 K=1,4
        IF0=3*(KCH(K)-1)
        DO 10 IQ=1,3
        DO 10 LF=1,3
10      DDC(IC0+IQ)=DDC(IC0+IQ)+DSP(LF,IQ,K)*DDF(IF0+LF)
20    CONTINUE
C
      END
C
************************************************************************
*   BNTOL0: P1disc pressure into L2_0 on level ILEV (TOL20A for IELT=4)
*     p := p - (int p)/|Omega| (Bernstein: constant c = all coeff. c)
************************************************************************
      SUBROUTINE BNTOL0(P)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION P(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      SAVE
C
      KWP=L(KLBWP(ILEV))-1
      PINT=0D0
      DO 10 IQ=1,KNP(ILEV)
10    PINT=PINT+DWORK(KWP+IQ)*P(IQ)
      C=PINT/DBNAR(ILEV)
      DO 20 IQ=1,KNP(ILEV)
20    P(IQ)=P(IQ)-C
C
      END
C
************************************************************************
*   BNMPAP: P := P + A1 * Mp^-1 FP (finest level), block diagonal
*           P1disc pressure mass matrix (XDFKDV for IELT=4)
************************************************************************
      SUBROUTINE BNMPAP(P,FP,A1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION P(*),FP(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      SAVE
C
      KMPI=L(KLBMPI(NLEV))-1
      DO 20 IEL=1,KNEL(NLEV)
       I0=NPL*(IEL-1)
       K0=KMPI+NPL*NPL*(IEL-1)
       DO 10 I=1,NPL
        S=0D0
        DO 5 J=1,NPL
5       S=S+DWORK(K0+I+NPL*(J-1))*FP(I0+J)
        P(I0+I)=P(I0+I)+A1*S
10     CONTINUE
20    CONTINUE
C
      END
C
************************************************************************
*   BNINTP: P1disc pressure -> vertex values (INTPV for IELT=4):
*     area-weighted mean of the element values in the vertex
************************************************************************
      SUBROUTINE BNINTP(DP,DPL,DAUX,AREA,KVERT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      REAL AREA
      PARAMETER (NNVE=4)
      DIMENSION DP(*),DPL(*),DAUX(*),KVERT(NNVE,*),AREA(*)
      DIMENSION CX(4),CY(4),DPV(6)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
      DATA CX/-1D0,1D0,1D0,-1D0/, CY/-1D0,-1D0,1D0,1D0/
      SAVE
C
      DO 1 IVT=1,NVT
      DPL (IVT)=0D0
1     DAUX(IVT)=0D0
      DO 10 IEL=1,NEL
       DAREA=DBLE(AREA(IEL))
       DO 10 IV=1,4
        CALL BNPB(KORD-1,CX(IV),CY(IV),DPV)
        IF (NPL.EQ.3) THEN
         PV=DPV(1)*DP(3*IEL-2)+DPV(2)*DP(3*IEL-1)+DPV(3)*DP(3*IEL)
        ELSE
         PV=0D0
         DO 8 LP=1,NPL
8        PV=PV+DPV(LP)*DP(NPL*(IEL-1)+LP)
        ENDIF
        IVT=KVERT(IV,IEL)
        DPL (IVT)=DPL (IVT)+DAREA*PV
        DAUX(IVT)=DAUX(IVT)+DAREA
10    CONTINUE
      DO 20 IVT=1,NVT
20    DPL(IVT)=DPL(IVT)/DAUX(IVT)
C
      END
C
************************************************************************
*   BNFVOL: drag/lift on boundary component 2 by the volume integral,
*     BDFVOL for IELT=4 (Bernstein Q2/P1disc, parametric EB20):
*     F_i = - [ (dt u_i + u.grad u_i, a) + nu (grad u_i, grad a)
*               - (p, da/dx_i) ],  c = 2 F / DPF(2)
*     (Q3/P2disc: KORD=3, 5x5 Gauss)
*     a = sum of the Q2 basis functions of the dofs on the cylinder
*     edges (vertices + edge dofs; = 1 on these edges, partition of
*     unity); 4x4 Gauss (ICUB=10). Arguments as BDFVOL.
************************************************************************
      SUBROUTINE BNFVOL(DU1,DU2,DP,DU1O,DU2O,DTK,DALPHA,
     *                  KVERT,KMID,KVBD,KEBD,KMM,DCORVG,ELE,DFWV,DAWV,
     *                  DFWP,DAWP,DFWN,DAWN)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      DIMENSION DU1(*),DU2(*),DP(*),DU1O(*),DU2O(*),DALPHA(*)
      DIMENSION DCORVG(2,*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*),KMM(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),DPV(6)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      COMMON /BNFAM/  IFAM,KORD,NLDF,NPL,NG
C *** isoparametric boundary elements (IELT>=4, IISO>0), bniso.f
      COMMON /BNISO/  IISO,MGEO,KLCRV(9),KLGEO(9),NCRV(9)
      SAVE
      EXTERNAL ELE
C
      DFWV=0D0
      DAWV=0D0
      DFWP=0D0
      DAWP=0D0
      DFWN=0D0
      DAWN=0D0
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
C *** test function a: Q2 dofs of the edges of component 2
      DO 1 IU=1,NU
1     DALPHA(IU)=0D0
      DO 10 IVBD=1,NVBD
      IF (KVBD(IVBD).EQ.KMM(1,2)) GOTO 15
10    CONTINUE
      WRITE(MTERM,*) 'ERROR 1 IN BNFVOL'
      RETURN
15    IVBD1=IVBD
      DO 20 IVBD=IVBD1,NVBD
      IELB=KEBD(IVBD)
      IVT1=KVBD(IVBD)
      IF (IVT1.EQ.KMM(2,2)) THEN
       IVT2=KVBD(IVBD1)
      ELSE
       IVT2=KVBD(IVBD+1)
      ENDIF
      DO 21 II=1,4
      JJ=MOD(II,4)+1
      IF (((KVERT(II,IELB).EQ.IVT1).AND.(KVERT(JJ,IELB).EQ.IVT2)).OR.
     *    ((KVERT(II,IELB).EQ.IVT2).AND.(KVERT(JJ,IELB).EQ.IVT1))) THEN
       IF (KORD.EQ.3) THEN
        IB3=NVT+2*(KMID(II,IELB)-NVT-1)
        DALPHA(IB3+1)=1D0
        DALPHA(IB3+2)=1D0
       ELSE
        DALPHA(KMID(II,IELB))=1D0
       ENDIF
       DALPHA(IVT1)=1D0
       DALPHA(IVT2)=1D0
      ENDIF
21    CONTINUE
      IF (IVT1.EQ.KMM(2,2)) GOTO 25
20    CONTINUE
25    CONTINUE
C *** volume integral, 4x4 Gauss formula, parametric element
      DO 30 IDER=1,NNDER
30    BDER(IDER)=.FALSE.
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      ICUB=10
      IF (KORD.EQ.3) ICUB=11
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
      DTKI=0D0
      IF (DTK.GT.0D0) DTKI=1D0/DTK
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
      BCYL=.FALSE.
      DO 101 IDFL1=1,IDFL
101   IF (DALPHA(KDFG(IDFL1)).NE.0D0) BCYL=.TRUE.
      IF (.NOT.BCYL) GOTO 100
      DO 110 IVE=1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
110   CONTINUE
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
      DO 120 ICUBP=1,NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
      IF (IISO.GT.0) CALL BNGJAC(IEL,XI1,XI2,DOMEGA(ICUBP),DJAC,DETJ,
     *                           OM,XXC,YYC,BCRV)
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
      U1V=0D0
      U1X=0D0
      U1Y=0D0
      U2V=0D0
      U2X=0D0
      U2Y=0D0
      DTU1=0D0
      DTU2=0D0
      AV=0D0
      AX=0D0
      AY=0D0
      DO 130 IDFL1=1,IDFL
      IG=KDFG(IDFL1)
      DB1=DBAS(KDFL(IDFL1),1)
      DB2=DBAS(KDFL(IDFL1),2)
      DB3=DBAS(KDFL(IDFL1),3)
      U1V=U1V+DU1(IG)*DB1
      U1X=U1X+DU1(IG)*DB2
      U1Y=U1Y+DU1(IG)*DB3
      U2V=U2V+DU2(IG)*DB1
      U2X=U2X+DU2(IG)*DB2
      U2Y=U2Y+DU2(IG)*DB3
      DTU1=DTU1+(DU1(IG)-DU1O(IG))*DB1
      DTU2=DTU2+(DU2(IG)-DU2O(IG))*DB1
      AV=AV+DALPHA(IG)*DB1
      AX=AX+DALPHA(IG)*DB2
      AY=AY+DALPHA(IG)*DB3
130   CONTINUE
C *** P1disc pressure in the cubature point
      CALL BNPB(KORD-1,XI1,XI2,DPV)
      IF (NPL.EQ.3) THEN
       DPCONT=DPV(1)*DP(3*IEL-2)+DPV(2)*DP(3*IEL-1)+DPV(3)*DP(3*IEL)
      ELSE
       DPCONT=0D0
       DO 135 LP=1,NPL
135    DPCONT=DPCONT+DPV(LP)*DP(NPL*(IEL-1)+LP)
      ENDIF
      AH1=(DTKI*DTU1+U1V*U1X+U2V*U1Y)*AV
     *   +DPF(1)*(U1X*AX+U1Y*AY)-DPCONT*AX
      AH2=(DTKI*DTU2+U1V*U2X+U2V*U2Y)*AV
     *   +DPF(1)*(U2X*AX+U2Y*AY)-DPCONT*AY
      DFWV=DFWV-AH1*OM
      DAWV=DAWV-AH2*OM
      DFWP=DFWP+DPCONT*AX*OM
      DAWP=DAWP+DPCONT*AY*OM
      DFWN=DFWN-DPF(1)*(U1X*AX+U1Y*AY)*OM
      DAWN=DAWN-DPF(1)*(U2X*AX+U2Y*AY)*OM
120   CONTINUE
100   CONTINUE
      DFWV=2D0*DFWV/DPF(2)
      DAWV=2D0*DAWV/DPF(2)
      DFWP=2D0*DFWP/DPF(2)
      DAWP=2D0*DAWP/DPF(2)
      DFWN=2D0*DFWN/DPF(2)
      DAWN=2D0*DAWN/DPF(2)
99999 END
C
************************************************************************
*   BNSRT7: sort the off-diagonal entries of each row ascending by
*           column (insertion sort), diagonal stays first
************************************************************************
      SUBROUTINE BNSRT7(NEQ,KLD,KCOL,DA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KLD(*),KCOL(*),DA(*)
C
      DO 30 IEQ=1,NEQ
       DO 20 IA=KLD(IEQ)+2,KLD(IEQ+1)-1
        JC=KCOL(IA)
        DH=DA(IA)
        IB=IA-1
10      IF (IB.GT.KLD(IEQ)) THEN
         IF (KCOL(IB).GT.JC) THEN
          KCOL(IB+1)=KCOL(IB)
          DA(IB+1)=DA(IB)
          IB=IB-1
          GOTO 10
         ENDIF
        ENDIF
        KCOL(IB+1)=JC
        DA(IB+1)=DH
20     CONTINUE
30    CONTINUE
C
      END
C
************************************************************************
*   BNRSOL: solution restriction for the coarse convection matrices
*     (RESTRU for IELT=4). The coarse Q2 nodes are fine vertices with
*     the same numbers (FEAT refinement: coarse vertices, NVT1+IMT,
*     NVT1+NMT1+IEL; checked in BNOWN), and Bernstein vertex
*     coefficients are point values: coarse nodal values = fine vertex
*     coefficients, converted to Bernstein coefficients (Q2
*     interpolation of the fine solution on the coarse grid):
*     vertex c=v, edge c=2v_m-(v_a+v_b)/2,
*     centre c=4v_c-sum(v_m)+sum(v_vertex)/4
************************************************************************
      SUBROUTINE BNRSOL(DU1C,DU2C,DU1F,DU2F,ILEVC)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION DU1C(*),DU2C(*),DU1F(*),DU2F(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /BNLEV/  IBNPP,KLBOWN(NNLEV),KLBWP(NNLEV),KLBMPI(NNLEV),
     *                KLBBV(NNLEV),DBNAR(NNLEV)
      DIMENSION KWORK(1),KDF(9)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF (IBNPP.EQ.3) THEN
       CALL BNRSLG(DU1C,DU2C,DU1F,DU2F,ILEVC,3)
       RETURN
      ENDIF
      KV1=L(KLVERT(ILEVC))
      KM1=L(KLMID(ILEVC))
      NVT1=KNVT(ILEVC)
      NMT1=KNMT(ILEVC)
      DO 30 IEL1=1,KNEL(ILEVC)
       CALL BNQ2DF(IEL1,KWORK(KV1),KWORK(KM1),NVT1,NMT1,KDF)
       DO 20 IC=1,2
        DO 10 I=1,4
         IV=KDF(I)
         IM=KDF(4+I)
         IA=KDF(I)
         IB=KDF(MOD(I,4)+1)
         IF (IC.EQ.1) THEN
          DU1C(IV)=DU1F(IV)
          DU1C(IM)=2D0*DU1F(IM)-0.5D0*(DU1F(IA)+DU1F(IB))
         ELSE
          DU2C(IV)=DU2F(IV)
          DU2C(IM)=2D0*DU2F(IM)-0.5D0*(DU2F(IA)+DU2F(IB))
         ENDIF
10      CONTINUE
        IF (IC.EQ.1) THEN
         DU1C(KDF(9))=4D0*DU1F(KDF(9))
     *        -(DU1F(KDF(5))+DU1F(KDF(6))+DU1F(KDF(7))+DU1F(KDF(8)))
     *        +0.25D0*(DU1F(KDF(1))+DU1F(KDF(2))+DU1F(KDF(3))
     *                +DU1F(KDF(4)))
        ELSE
         DU2C(KDF(9))=4D0*DU2F(KDF(9))
     *        -(DU2F(KDF(5))+DU2F(KDF(6))+DU2F(KDF(7))+DU2F(KDF(8)))
     *        +0.25D0*(DU2F(KDF(1))+DU2F(KDF(2))+DU2F(KDF(3))
     *                +DU2F(KDF(4)))
        ENDIF
20     CONTINUE
30    CONTINUE
C
      END
