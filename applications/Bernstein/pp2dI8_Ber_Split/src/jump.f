************************************************************************
*     Edge-oriented jump stabilization (IJUMP>0), ported from
*     applications/myapp/cc2d_powerlaw_new (stabil.f: JUMPC, ADJELEM,
*     BASEJUMP, BASESF; util.f: apex71/apex72; jmpsep.f)
*
*       J(u,v) = DJUMP * sum_E  sum_g  w_g h_E^3  [grad u]:[grad v](x_g)
*              = DJUMP * sum_E  2 h_E^2  int_E [grad u]:[grad v] ds
*
*     (IJUMP=2; 3-point Gauss formula on each inner edge E, weights w_g
*      on [-1,1]; IJUMP=1: nu*[u][v] with weight w_g, IJUMP=3/4: only
*      the normal/tangential derivative). J does not depend on the
*     solution and acts on both velocity components alike (blocks 11
*     and 22), so it is assembled once per level (JMPINI) and stored
*     in DOUBLE on the structure of A, which is extended for this to
*     the A*A pattern (JMPSTR: the 7 dofs of the two elements at an
*     edge are coupled).
*
*     J is added in SUPWDG/SUPWNP (they call JMPADD), i.e. in all
*     matrices, defects and right hand sides, with the weight
*       IPRECA=4: THSTEP*VISFAC, like nu*A (linear, symmetric);
*                 ISPLIT>0: theta-scheme splitting as for nu*A,
*                 Stokes  alpha*theta*K impl. + beta*theta*K expl.,
*                 Burgers beta*theta'*K impl. + alpha*theta'*K expl.
*       IPRECA<4: THSTEP*CNVFAC, like the convection (ISPLIT=0 only;
*                 for ISPLIT>0 J would be explicit in the Stokes
*                 substeps, which was unstable for K=0.05, DJUMP=0.01)
*     In the original scheme both weights are the same.
*     Only nonparametric elements (IELT=2,3) and IUPW=0.
************************************************************************
      SUBROUTINE JMPSTR(KLCOLA,KLLDA,KNA,KNEQ)
************************************************************************
*     Extends the structure of A on all levels to the A*A pattern
*     (diagonal first, then the other columns in increasing order)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION KLCOLA(*),KLLDA(*),KNA(*),KNEQ(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE
C
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      DO 10 JLEV=NLMIN,NLMAX
      NEQ=KNEQ(JLEV)
      CALL ZNEW(NEQ,-3,LAUX,'KAUXJ ')
      IF (IER.NE.0) GOTO 99999
      CALL APEX71(NEQ,KWORK(L(KLCOLA(JLEV))),KWORK(L(KLLDA(JLEV))),
     *            NA1,KWORK(L(LAUX)))
      CALL ZNEW(NA1,-3,LCOL1,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NEQ+1,-3,LLD1,'KLD   ')
      IF (IER.NE.0) GOTO 99999
      CALL APEX72(NEQ,KWORK(L(KLCOLA(JLEV))),KWORK(L(KLLDA(JLEV))),
     *            KWORK(L(LCOL1)),KWORK(L(LLD1)),NA1,KWORK(L(LAUX)))
      IF (IER.NE.0) GOTO 99999
C
      WRITE(MTERM,'(A,I3,A,I10,A,I10)') ' jump: level',JLEV,
     *  '  NA =',KNA(JLEV),'  -> extended NA =',NA1
C
      CALL ZDISP(0,LAUX,'KAUXJ ')
      CALL ZDISP(0,KLCOLA(JLEV),'KCOL  ')
      CALL ZDISP(0,KLLDA(JLEV),'KLD   ')
      IF (IER.NE.0) GOTO 99999
      KLCOLA(JLEV)=LCOL1
      KLLDA(JLEV) =LLD1
      KNA(JLEV)   =NA1
10    CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE APEX71(NEQ,KCOL,KLD,NA1,KAUX)
************************************************************************
*     Number of nonzeros NA1 of the A*A pattern: a(i,j)<>0 and
*     a(j,k)<>0 -> (i,k) (KAUX(k)=i marks k as counted in row i)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KCOL(*),KLD(*),KAUX(*)
C
      CALL LCL3(KAUX,NEQ)
      NA1=0
      DO 10 IEQ=1,NEQ
      DO 20 IA=KLD(IEQ),KLD(IEQ+1)-1
      JEQ=KCOL(IA)
      DO 30 JA=KLD(JEQ),KLD(JEQ+1)-1
      KEQ=KCOL(JA)
      IF (KAUX(KEQ).NE.IEQ) THEN
       KAUX(KEQ)=IEQ
       NA1=NA1+1
      ENDIF
30    CONTINUE
20    CONTINUE
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE APEX72(NEQ,KCOL,KLD,KCOL1,KLD1,NA1,KAUX)
************************************************************************
*     Fills the A*A pattern KCOL1,KLD1 (storage technique 7: diagonal
*     first, off-diagonal columns sorted increasingly)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION KCOL(*),KLD(*),KCOL1(*),KLD1(*),KAUX(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
      CALL LCL3(KAUX,NEQ)
      ICOL1=1
      DO 10 IEQ=1,NEQ
      KLD1(IEQ)=ICOL1
      KCOL1(ICOL1)=IEQ
      KAUX(IEQ)=IEQ
      ICOL1=ICOL1+1
      IROW1=ICOL1
C
      DO 20 IA=KLD(IEQ),KLD(IEQ+1)-1
      JEQ=KCOL(IA)
      DO 30 JA=KLD(JEQ),KLD(JEQ+1)-1
      KEQ=KCOL(JA)
      IF (KAUX(KEQ).NE.IEQ) THEN
       KAUX(KEQ)=IEQ
       KCOL1(ICOL1)=KEQ
       ICOL1=ICOL1+1
      ENDIF
30    CONTINUE
20    CONTINUE
C
C *** insertion sort of the off-diagonal columns of row IEQ
      DO 40 IA=IROW1+1,ICOL1-1
      KEQ=KCOL1(IA)
      JA=IA-1
41    IF (JA.GE.IROW1) THEN
       IF (KCOL1(JA).GT.KEQ) THEN
        KCOL1(JA+1)=KCOL1(JA)
        JA=JA-1
        GOTO 41
       ENDIF
      ENDIF
      KCOL1(JA+1)=KEQ
40    CONTINUE
C
10    CONTINUE
      KLD1(NEQ+1)=ICOL1
C
      IF (NA1.NE.ICOL1-1) THEN
       WRITE(MTERM,*) 'ERROR in APEX72: wrong number of nonzeros'
       IER=-1
      ENDIF
C
      END
C
C
C
************************************************************************
      SUBROUTINE JMPINI
************************************************************************
*     Assembles J (without THSTEP) on all levels into KLJMP(ILEV)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
      INTEGER  VIPARM
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C *** jump stabilization (GDAT), J of each level
      COMMON /NSJUMP/ DJUMP,IJUMP,KLJMP(NNLEV)
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
      SAVE
C
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL EM30,EM31
C
      DO 10 JLEV=1,NNLEV
10    KLJMP(JLEV)=0
      IF (IJUMP.EQ.0) RETURN
C
      DO 20 JLEV=NLMIN,NLMAX
      ILEV=JLEV
      ISETLV=2
      CALL SETLEV(ISETLV)
C
      CALL ZNEW(NA,1,KLJMP(JLEV),'DJUMP ')
      IF (IER.NE.0) GOTO 99999
C *** the two elements at each edge (KMEL is not built by XMSB2)
      CALL ZNEW(2*NMT,3,LMELJ,'KMELJ ')
      IF (IER.NE.0) GOTO 99999
      CALL JMPMEL(KWORK(L(LMID)),KWORK(L(LMELJ)))
C
      IF (IELT.EQ.2)
     * CALL JMPASM(DWORK(L(KLJMP(JLEV))),KWORK(KCOLA),KWORK(KLDA),
     *             KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LMELJ)),
     *             DWORK(L(LCORVG)),EM31)
      IF (IELT.EQ.3)
     * CALL JMPASM(DWORK(L(KLJMP(JLEV))),KWORK(KCOLA),KWORK(KLDA),
     *             KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LMELJ)),
     *             DWORK(L(LCORVG)),EM30)
      IF (IER.NE.0) GOTO 99999
C *** IELT=4 (Bernstein Q2): J = DJUMP*sum_E 2 h_E^2 int_E [grad u].
C     [grad v] ds by BNJASM (bnafc.f, 3 Gauss points per edge); FEAT's
C     Q2 numbering (vertices, NVT+IMT, NVT+NMT+IEL) is the Bernstein
C     numbering of BNDFMP for k=2. Check: J*x = J*y = 0 (Greville
C     points = coefficients of the linear functions x, y)
      IF (IELT.GE.4) THEN
C ***  k=IELT-2 (Q2: FEAT type 13, Q3: type 14 = BNDFMP numbering)
       KBJ=IELT-2
       NUJ=NVT+(KBJ-1)*NMT+(KBJ-1)**2*NEL
       CALL ZNEW((KBJ+1)**2*NEL,3,LDFGJ,'KDFGJ ')
       IF (IER.NE.0) GOTO 99999
       CALL ZNEW(20*NEL,3,LEDGJ,'KEDGJ ')
       IF (IER.NE.0) GOTO 99999
       CALL ZNEW(3*NUJ ,1,LXYJ ,'DXYJ  ')
       IF (IER.NE.0) GOTO 99999
       CALL BNDFMP(KBJ,KWORK(L(LVERT)),KWORK(L(LMID)),NEL,NVT,NMT,
     *             KWORK(L(LDFGJ)))
       CALL BNJEDG(KWORK(L(LADJ)),KWORK(L(LVERT)),NEL,
     *             KWORK(L(LEDGJ)),NEDGJ)
       CALL BNJASM(KBJ,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFGJ)),KWORK(L(LEDGJ)),NEDGJ,
     *             KWORK(KLDA),KWORK(KCOLA),DWORK(L(KLJMP(JLEV))))
       CALL LSC1(DWORK(L(KLJMP(JLEV))),NA,2D0*DJUMP)
       CALL BNDFXY(KBJ,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFGJ)),DWORK(L(LXYJ)),DWORK(L(LXYJ)+NUJ))
       CALL BNCHKX(KBJ,NEL,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *             KWORK(L(LDFGJ)),DWORK(L(LXYJ)),DWORK(L(LXYJ)+NUJ),
     *             ERRXS,ERRXC)
       WRITE(MFILE1,'(A,I3,A,D10.2,A,D10.2)') ' JMPINI: level',JLEV,
     *  ' |grad x-e|: bilinear el.',ERRXS,'  curved el.',ERRXC
       DJMXJ=0D0
       DO 30 IA=1,NA
30     DJMXJ=MAX(DJMXJ,ABS(DWORK(L(KLJMP(JLEV))+IA-1)))
       DJLNJ=0D0
       DO 32 IXY=0,1
        CALL BNAX(DWORK(L(KLJMP(JLEV))),KWORK(KLDA),KWORK(KCOLA),NUJ,
     *            DWORK(L(LXYJ)+IXY*NUJ),DWORK(L(LXYJ)+2*NUJ))
        DO 32 I=1,NUJ
32      DJLNJ=MAX(DJLNJ,ABS(DWORK(L(LXYJ)+2*NUJ+I-1)))
       IF (DJMXJ.GT.0D0) DJLNJ=DJLNJ/DJMXJ
       WRITE(MTERM,1001) JLEV,NEDGJ,NA,DJLNJ
       WRITE(MFILE1,1001) JLEV,NEDGJ,NA,DJLNJ
       CALL ZDISP(0,LXYJ ,'DXYJ  ')
       CALL ZDISP(0,LEDGJ,'KEDGJ ')
       CALL ZDISP(0,LDFGJ,'KDFGJ ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZDISP(0,LMELJ,'KMELJ ')
      IF (IER.NE.0) GOTO 99999
20    CONTINUE
C
      ILEV=NLMAX
      ISETLV=2
      CALL SETLEV(ISETLV)
C
1001  FORMAT(' JMPINI IELT=4: level',I3,' inner edges',I8,' NA',I9,
     *       ' max|J*x|,|J*y|/max|J| =',D9.2)
99999 END
C
C
C
************************************************************************
      SUBROUTINE JMPMEL(KMID,KMELJ)
************************************************************************
*     KMELJ(1:2,IMT) = elements at edge IMT (KMELJ(2,IMT)=0: boundary)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION KMID(NNVE,*),KMELJ(2,*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      DO 10 IEL=1,NEL
      DO 10 IVE=1,NNVE
      IMT=KMID(IVE,IEL)-NVT
      IF (KMELJ(1,IMT).EQ.0) THEN
       KMELJ(1,IMT)=IEL
      ELSE
       KMELJ(2,IMT)=IEL
      ENDIF
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE JMPASM(DJ,KCOLA,KLDA,KVERT,KMID,KMELJ,DCORVG,ELE)
************************************************************************
*     J = DJUMP * jump form, loop over all inner edges (after JUMPC
*     of cc2d_powerlaw_new, solution independent forms IJUMP=1..4)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNBAS=21,NNDER=6,NNVE=4,NNLEV=9)
      DIMENSION DJ(*),KCOLA(*),KLDA(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KMELJ(2,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DENTRY(7,7),DXG(3),DYG(3),DALPHA(3)
      DIMENSION ILEADJ(2),INUG(8),INUG7(7)
      DIMENSION DBASSF(9,7),DBASEJ(9,8)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
      COMMON /NSJUMP/ DJUMP,IJUMP,KLJMP(NNLEV)
      SAVE
      EXTERNAL ELE
C
      ICUB=3
      IDFL7=7
      DNY=NY
C
      DO 1 I=1,NNDER
1     BDER(I)=.FALSE.
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL ELE(0D0,0D0,-2)
C
C *** 3-point Gauss formula on [-1,1]
      DCONS2=SQRT(3D0/5D0)
      DALPHA(1)=8D0/9D0
      DALPHA(2)=5D0/9D0
      DALPHA(3)=DALPHA(2)
C
C *** Loop over all inner edges
      DO 100 IMT=1,NMT
C
      CALL ADJELEM(KVERT,KMID,KMELJ,IMT,NVT,ILEADJ,IVT1,IVT2,
     *             INUG,INUG7)
      IF (ILEADJ(2).EQ.0) GOTO 100
C
      DO 110 JDOFE=1,IDFL7
      DO 110 IDOFE=1,IDFL7
110   DENTRY(JDOFE,IDOFE)=0D0
C
      PX1=DCORVG(1,IVT1)
      PY1=DCORVG(2,IVT1)
      PX2=DCORVG(1,IVT2)
      PY2=DCORVG(2,IVT2)
      HLOCAL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
      DT1=(PX2-PX1)/HLOCAL
      DT2=(PY2-PY1)/HLOCAL
C
      DXG(1)=0.5D0*(PX1+PX2)
      DXG(2)=0.5D0*( DCONS2*(PX2-PX1)+(PX1+PX2))
      DXG(3)=0.5D0*(-DCONS2*(PX2-PX1)+(PX1+PX2))
      DYG(1)=0.5D0*(PY1+PY2)
      DYG(2)=0.5D0*( DCONS2*(PY2-PY1)+(PY1+PY2))
      DYG(3)=0.5D0*(-DCONS2*(PY2-PY1)+(PY1+PY2))
C
      DO 200 ICUBP=1,ICUB
      XI1=DXG(ICUBP)
      XI2=DYG(ICUBP)
      OM1=DALPHA(ICUBP)
      OM3=DALPHA(ICUBP)*HLOCAL**3
C
C *** basis functions of both elements and the jumps across IMT
      CALL BASEJUMP(KVERT,ILEADJ,DT1,DT2,DCORVG,ELE,XI1,XI2,DBASEJ)
      CALL BASESF(IMT,INUG,DBASEJ,DBASSF)
C
      DO 240 JDOFE=1,IDFL7
      DO 240 IDOFE=1,IDFL7
      IF (IJUMP.EQ.1) AH=OM1*DNY*DBASSF(1,JDOFE)*DBASSF(1,IDOFE)
      IF (IJUMP.EQ.2) AH=OM3*(DBASSF(2,JDOFE)*DBASSF(2,IDOFE)
     *                       +DBASSF(3,JDOFE)*DBASSF(3,IDOFE))
      IF (IJUMP.EQ.3) AH=OM3*DBASSF(4,JDOFE)*DBASSF(4,IDOFE)
      IF (IJUMP.EQ.4) AH=OM3*DBASSF(5,JDOFE)*DBASSF(5,IDOFE)
      DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+AH
240   CONTINUE
C
200   CONTINUE
C
C *** add DJUMP*DENTRY into J (the entries exist in the A*A pattern)
      DO 400 JDOFE=1,IDFL7
      JDFG=INUG7(JDOFE)
      DO 400 IDOFE=1,IDFL7
      IDFG=INUG7(IDOFE)
      DO 410 ICOL=KLDA(JDFG),KLDA(JDFG+1)-1
      IF (KCOLA(ICOL).EQ.IDFG) THEN
       DJ(ICOL)=DJ(ICOL)+DJUMP*DENTRY(JDOFE,IDOFE)
       GOTO 400
      ENDIF
410   CONTINUE
      WRITE(MTERM,*) 'ERROR in JMPASM: entry not in the A*A pattern',
     *               JDFG,IDFG
      IER=-1
      GOTO 99999
400   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE ADJELEM(KVERT,KMID,KMELJ,IMT,NVT,ILEADJ,IVT1,IVT2,
     *                   INUG,INUG7)
************************************************************************
*     Elements ILEADJ(1:2) at edge IMT, its vertices IVT1,IVT2,
*     INUG(1:8) = edge dofs of both elements, INUG7(1:7) = the 7
*     different dofs (element 1, then element 2 without IMT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KMELJ(2,*)
      DIMENSION ILEADJ(2),INUG(8),INUG7(7)
C
      ILEADJ(1)=KMELJ(1,IMT)
      ILEADJ(2)=KMELJ(2,IMT)
      IEL1=ILEADJ(1)
C
      DO 10 IVE=1,NNVE
      IVEN=MOD(IVE,NNVE)+1
      IF (KMID(IVE,IEL1)-NVT.EQ.IMT) THEN
       IVT1=KVERT(IVE ,IEL1)
       IVT2=KVERT(IVEN,IEL1)
      ENDIF
10    CONTINUE
C
      DO 20 IN=1,8
20    INUG(IN)=0
      DO 21 IN=1,7
21    INUG7(IN)=0
C
      DO 30 IVE=1,NNVE
      INUG(IVE) =KMID(IVE,IEL1)-NVT
      INUG7(IVE)=KMID(IVE,IEL1)-NVT
30    CONTINUE
C
      IF (ILEADJ(2).NE.0) THEN
       ICONT=0
       DO 40 IVE=1,NNVE
       IMID=KMID(IVE,ILEADJ(2))-NVT
       INUG(NNVE+IVE)=IMID
       IF (IMID.NE.IMT) THEN
        ICONT=ICONT+1
        INUG7(NNVE+ICONT)=IMID
       ENDIF
40     CONTINUE
      ENDIF
C
      END
C
C
C
************************************************************************
      SUBROUTINE BASEJUMP(KVERT,ILEADJ,DT1,DT2,DCORVG,ELE,XI1,XI2,
     *                    DBASEJ)
************************************************************************
*     Values/derivatives of the basis functions of both elements at
*     the point (XI1,XI2) on the edge (nonparametric element):
*     DBASEJ(1..3,.) = phi, phi_x, phi_y; (4) normal, (5) tangential
*     derivative; (6..9) as in cc2d_powerlaw_new (not used here)
*     columns 1..4: element ILEADJ(1), 5..8: element ILEADJ(2)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),DCORVG(2,*),ILEADJ(2),DBASEJ(9,8)
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      SAVE
      EXTERNAL ELE
C
      DO 1 I=1,8
      DO 1 J=1,9
1     DBASEJ(J,I)=0D0
C
      DN1=-DT2
      DN2= DT1
C
      DO 100 IADJ=1,2
      IF (ILEADJ(IADJ).EQ.0) GOTO 100
      IEL=ILEADJ(IADJ)
      IOFF=(IADJ-1)*NNVE
      DO 10 IVE=1,NNVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
10    CONTINUE
      CALL ELE(0D0,0D0,-2)
      CALL ELE(XI1,XI2,-3)
      DO 20 IVE=1,NNVE
      DBASEJ(1,IOFF+IVE)=DBAS(IVE,1)
      DBASEJ(2,IOFF+IVE)=DBAS(IVE,2)
      DBASEJ(3,IOFF+IVE)=DBAS(IVE,3)
      DBASEJ(4,IOFF+IVE)=DN1*DBAS(IVE,2)+DN2*DBAS(IVE,3)
      DBASEJ(5,IOFF+IVE)=DT1*DBAS(IVE,2)+DT2*DBAS(IVE,3)
      DBASEJ(6,IOFF+IVE)=DBASEJ(5,IOFF+IVE)*DN1
      DBASEJ(7,IOFF+IVE)=DBASEJ(5,IOFF+IVE)*DN2
      DBASEJ(8,IOFF+IVE)=DBASEJ(4,IOFF+IVE)*DT1
      DBASEJ(9,IOFF+IVE)=DBASEJ(4,IOFF+IVE)*DT2
20    CONTINUE
100   CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE BASESF(IMT,INUG,DBASEJ,DBASSF)
************************************************************************
*     Jumps of the 7 basis functions across edge IMT (ordering INUG7):
*     element 1 positive, element 2 negative, the common dof IMT gets
*     the difference of both
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION INUG(8),DBASEJ(9,8),DBASSF(9,7)
C
      DO 1 I=1,7
      DO 1 J=1,9
1     DBASSF(J,I)=0D0
C
      IND=0
      DO 10 I=1,4
      IF (INUG(I).EQ.IMT) IND=I
      DO 10 J=1,9
10    DBASSF(J,I)=DBASEJ(J,I)
C
      ICONT=0
      DO 20 I=5,8
      IF (INUG(I).EQ.IMT) THEN
       DO 21 J=1,9
21     DBASSF(J,IND)=DBASSF(J,IND)-DBASEJ(J,I)
      ELSE
       ICONT=ICONT+1
       DO 22 J=1,9
22     DBASSF(J,4+ICONT)=-DBASEJ(J,I)
      ENDIF
20    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE JMPADD(A,NA,KCOLA,KLDA,U1,U2,D1,D2,IDEF,ICOMP,DFAC)
************************************************************************
*     Adds DFAC*J of the current level (identified by NMT) to the
*     REAL matrix A (IDEF<2) and -DFAC*J*U to the defect D (IDEF>0),
*     as SUPWDG does for the convection
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      REAL A
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      DIMENSION A(*),KCOLA(*),KLDA(*),U1(*),U2(*),D1(*),D2(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /NSJUMP/ DJUMP,IJUMP,KLJMP(NNLEV)
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      SAVE
C
      IF ((IJUMP.EQ.0).OR.(DFAC.EQ.0D0)) RETURN
C
      DO 10 JLEV=NLMIN,NLMAX
      IF ((KNMT(JLEV).EQ.NMT).AND.(KLJMP(JLEV).NE.0)) GOTO 20
10    CONTINUE
      RETURN
C
20    KJ=L(KLJMP(JLEV))-1
C
      IF (IDEF.LT.2) THEN
       DO 30 IA=1,NA
30     A(IA)=A(IA)+REAL(DFAC*DWORK(KJ+IA))
      ENDIF
C
C *** defect over all velocity rows of the level (KNU = NMT for Q1~,
C     NVT+NMT+NEL for Q2; the loop over NMT missed the vertex/centre
C     rows of Q2, i.e. matrix and defect were inconsistent)
      IF (IDEF.GT.0) THEN
       DO 40 IEQ=1,KNU(JLEV)
       S1=0D0
       S2=0D0
       DO 41 IA=KLDA(IEQ),KLDA(IEQ+1)-1
       S1=S1+DWORK(KJ+IA)*U1(KCOLA(IA))
       S2=S2+DWORK(KJ+IA)*U2(KCOLA(IA))
41     CONTINUE
                       D1(IEQ)=D1(IEQ)-DFAC*S1
       IF (ICOMP.EQ.2) D2(IEQ)=D2(IEQ)-DFAC*S2
40     CONTINUE
      ENDIF
C
      END
