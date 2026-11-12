************************************************************************
      SUBROUTINE PROSTG (MFILE,MSHOW,ITMOD,ITNSL,TSTEPH,TMSTEP,
     *                   BSTOP,BNLEND)
************************************************************************
*
*   Purpose: - performs 1 substep of Glowinski's fractional-step-theta
*              operator splitting (called from PROSTP if ISPLIT>0)
*
*              K = macro time step = 3*TSTEP, theta' = 1-2*theta,
*              alpha = (1-2theta)/(1-theta), beta = theta/(1-theta)
*
*     ITMOD=1,3: generalized Stokes problem on a theta*K step
*       [M + alpha*theta*K*nu*A] u + theta*K*B p
*                  = [M - beta*theta*K*nu*A] u_old
*                    - theta*K*N(u_old)u_old + theta*K*f
*       B^T u = 0
*       solved by an Uzawa iteration whose preconditioner is the
*       discrete projection step of pp2d (Cahouet-Chabard:
*       S^-1 = PRDIF1/(theta*K) (-Lap_h)^-1 + PRDIF2*alpha*nu*Mp^-1)
*
*     ITMOD=2: Burgers problem on a theta'*K step (no constraint)
*       [M + beta*theta'*K*nu*A + theta'*K*N(w)] u
*                  = [M - alpha*theta'*K*nu*A] u_old
*                    - theta'*K*B p_old + theta'*K*f
*       ISPLIT=1: w = u_old                                  (1. order)
*       ISPLIT=2: w = u_old + theta'/theta*(u_old - u^n)     (2. order)
*       ISPLIT=3: w = u (fixed point iteration in NSDEF)
*
*   On entry (set by PROSTP): level NLMAX active, KF1/KF2 contain the
*   right hand side f-part, TIMENS = time of the new unknown,
*   THSTEP = -beta*theta*K (ITMOD=1,3) or -alpha*theta'*K (ITMOD=2).
*   TSTEPH = theta*K or theta'*K, TMSTEP = alpha*theta*K = beta*theta'*K
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNAB=21,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CSTEP*7
C
C *** Arrays for multigrid modul M010
      DIMENSION  KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),KNEQ(NNLEV),
     *           KIT(NNLEV),KIT0(NNLEV)
C-----------------------------------------------------------------------
C     C O M M O N S
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
      COMMON /NSCOUN/ NNONL,NMGU,NMGP
      COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32
      COMMON /NSSPLT/ EPSUZ,CNVFAC,ISPLIT,NUZMAX
      COMMON /NSREFI/ IREFIN,IFRZW,LWFRZ
C *** CNVBAS: base factor of the convection (1, or 0 = Stokes limit)
      COMMON /NSCNVB/ CNVBAS
C *** VISFAC: extra weight of nu*A in SUPWDG/SUPWNP for IPRECA=4
      COMMON /NSVISF/ VISFAC
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
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
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Parametrization of the domain
      EXTERNAL PARX,PARY
C *** Coefficient of stiffness matrix, exact solution
      EXTERNAL COEFFN,UE
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30,EB20,EBK
C *** Multigrid components
      EXTERNAL  YAXP,YPROLP,YRESTP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP
C
C *** Additional vectors (allocated once, kept for the whole run):
C     LF0  : RHS of the Stokes substep without the pressure term
C     LUN  : velocity u^n at the beginning of the macro step
C     LDIV : B^T u for the divergence output
      DATA LF0/0/,LUN/0/,LDIV/0/
C
C=======================================================================
C
      BSTOP =.FALSE.
      BNLEND=.TRUE.
C
      IF (LF0.EQ.0) THEN
       CALL ZNEW(2*NU,-1,LF0 ,'DF0GL ')
       CALL ZNEW(2*NU,-1,LUN ,'DUNGL ')
       CALL ZNEW(NP  ,-1,LDIV,'DDIVGL')
C ***  LWFRZ: frozen transport velocity w for the refined Burgers solve
       IF (IREFIN.EQ.1) CALL ZNEW(2*NU,-1,LWFRZ,'DWFRZ ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C *** save the parameters which are changed temporarily
      ISTOKH=ISTOK
      INLMIH=INLMIN
      INLMAH=INLMAX
      ITEXLH=ITEXL
      ITEXL =0
C
      NNONL0=NNONL
      NMGU0 =NMGU
      NMGP0 =NMGP
      IUZ   =0
      DIVUZ =0D0
      RELP2 =0D0
      RELPM =0D0
C
      BSTOK=(ITMOD.NE.2)
      IF (BSTOK) THEN
       CSTEP='Stokes '
      ELSE
       CSTEP='Burgers'
      ENDIF
C
C *** store u^n at the beginning of the macro step (for ISPLIT=2)
      IF (ITMOD.EQ.1) CALL LCP1(DWORK(KU1),DWORK(L(LUN)),2*NU)
C
C=======================================================================
C *** Right hand side
C     Stokes : F = f-part + [M - beta*theta*K*nu*A] u_old
C     Burgers: F = f-part + [M - alpha*theta'*K*nu*A] u_old
C                         - theta'*K*B p_old
C     (THSTEP carries the weight -beta*theta*K or -alpha*theta'*K)
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF (BSTOK) THEN
       IPRJ=0
      ELSE
       IPRJ=1
      ENDIF
      CALL XMADF1(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *            KF1,KF2,KU1,KU2,KP,NA,NU,THSTEP,TSTEPH,IPRJ)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C-----------------------------------------------------------------------
C *** Stokes: explicit convection  F = F - theta*K*N(u_old)u_old
C     (skipped for ISTOK=1, i.e. a real Stokes calculation)
C *** IPRECA=4: XMADF1 gives only F = f-part + Ml u_old, M (IMASS=1)
C     and nu*A are assembled in SUPWDG/SUPWNP as well:
C     Stokes : F = F + M u_old - theta*K*(beta*nu*A + N(u_old)) u_old
C              (THSTEP=theta*K, VISFAC=beta)
C     Burgers: F = F + M u_old - alpha*theta'*K*nu*A u_old
C              (THSTEP=alpha*theta'*K, CNVFAC=0)
C     The jump term J (IJUMP>0, needs IPRECA=4) goes with nu*A.
C-----------------------------------------------------------------------
C
      IF ((BSTOK.AND.(ISTOK.NE.1)).OR.(IPRECA.EQ.4)) THEN
       CALL ZTIME(TTT0)
       THSTEH=THSTEP
       IF (BSTOK) THEN
        THSTEP=TSTEPH
        IF (IPRECA.EQ.4) THEN
         VISFAC=FBETA
         IF (ISTOK.EQ.1) CNVFAC=0D0
        ENDIF
       ELSE
        THSTEP=-THSTEH
        CNVFAC=0D0
       ENDIF
       IF (IUPW.EQ.1) THEN
        CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *              1D0,0D0,DWORK(KU1),DWORK(KU2),
     *              DWORK(KF1),DWORK(KF2),VWORK(KA1),KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),2,2)
       ELSE
        IF (IELT.EQ.0)
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E031,COEFFN,2,2,-1D0)
        IF (IELT.EQ.1)
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E030,COEFFN,2,2,-1D0)
        IF (IELT.EQ.2)
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM31,COEFFN,2,2,-1D0)
        IF (IELT.EQ.3)
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM30,COEFFN,2,2,-1D0)
        IF (IELT.GE.4)
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EBK,COEFFN,2,2,-1D0)
       ENDIF
       THSTEP=THSTEH
       VISFAC=1D0
       CNVFAC=CNVBAS
       CALL ZTIME(TTT1)
       TTUPW=TTUPW+TTT1-TTT0
      ENDIF
C
C-----------------------------------------------------------------------
C *** Burgers, ISPLIT=2: transport velocity (and start vector)
C       w = u^(n+theta) + theta'/theta*(u^(n+theta) - u^n)
C     must be done before BDRSET, so that w satisfies the new
C     Dirichlet values
C-----------------------------------------------------------------------
C
      IF ((.NOT.BSTOK).AND.(ISPLIT.EQ.2)) THEN
       CALL LLC1(DWORK(L(LUN)),DWORK(KU1),2*NU,-THETAP/THETA,
     *           1D0+THETAP/THETA)
      ENDIF
C
C *** Dirichlet values at the time of the new unknown into U and F
      CALL ZTIME(TTT0)
      CALL BDRSET (DWORK(KU1),DWORK(KU2),DWORK(KF1),
     *             DWORK(KF2),KWORK(L(KLMBD(ILEV))),
     *             DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
     *             KNMBD(ILEV),NVT,PARX,PARY,UE,2)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
C *** matrix weight alpha*theta*K (= beta*theta'*K) for NSDEF
      THSTEP=TMSTEP
C
      IF (BSTOK) GOTO 300
C
C=======================================================================
C *** Substep 2: Burgers problem, no pressure update
C     [M + beta*theta'*K*nu*A + theta'*K*N(w)] u = F
C     NSDEF weights nu*A and N with the same THSTEP=beta*theta'*K,
C     so the convection gets the extra factor CNVFAC=CNVBAS/beta.
C     (Glowinski: full weight on the convection in this substep; with
C      weight beta*theta'*K the scheme would not be consistent, since
C      theta*K+beta*theta'*K+theta*K < K)
C=======================================================================
C
      IF (ISPLIT.LE.2) THEN
       IF ((IREFIN.EQ.1).AND.(INLMAX.GT.1)) THEN
C *** linearized with refinement: w = U frozen (LWFRZ), defect
C     correction in NSDEF until EPSUD (at most INLMAX steps)
        CALL LCP1(DWORK(KU1),DWORK(L(LWFRZ)),2*NU)
        IFRZW =1
        INLMIN=1
       ELSE
C *** linearized: one linear solve with the transport velocity w = U
        INLMIN=1
        INLMAX=1
       ENDIF
      ENDIF
C
      CNVFAC=CNVBAS/FBETA
      CALL  NSDEF (MFILE,MSHOW,BSTOP,BNLEND)
      CNVFAC=CNVBAS
      IF (IER.NE.0) GOTO 99998
      IF (BSTOP) GOTO 99997
      GOTO 500
C
C=======================================================================
C *** Substeps 1 and 3: generalized Stokes problem, Uzawa iteration
C     [M + alpha*theta*K*nu*A] u~ = F0 - theta*K*B p_k
C     p_k+1 = p_k + S^-1 B^T u~ ,  u = u~ - theta*K*Ml^-1 B (p_k+1-p_k)
C     until |B^T u~|/|u~| <= EPSUZ  (at most NUZMAX steps)
C=======================================================================
C
300   CALL LCP1(DWORK(KF1),DWORK(L(LF0)),2*NU)
C *** linear problem: no convection, no forced 2nd defect correction
C     step (it would stagnate at round-off and trigger BSTOP in NSDEF)
      ISTOK =1
      INLMIN=1
C
      DO 200 IUZ=1,NUZMAX
C
C-----------------------------------------------------------------------
C *** F = F0 - theta*K*B p_k   (Dirichlet rows of B p set to zero)
C-----------------------------------------------------------------------
C
      CALL ZTIME(TTT0)
      IF (IPRECB.EQ.0) THEN
       CALL  LAX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(L(LD1)),1D0,0D0)
       CALL  LAX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(L(LD2)),1D0,0D0)
      ENDIF
      IF (IPRECB.EQ.1) THEN
       CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(L(LD1)),1D0,0D0)
       CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(L(LD2)),1D0,0D0)
      ENDIF
      IF (IPRECB.EQ.2) THEN
       CALL LCL1(DWORK(L(LD1)),NU)
       CALL LCL1(DWORK(L(LD2)),NU)
       CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KP),DWORK(L(LD1)),
     *             DWORK(L(LD2)),NEL,NVT,NMT,1D0,0D0)
      ENDIF
      CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),KWORK(L(KLMBD(NLEV))),
     *            KNMBD(NLEV))
      CALL LCP1(DWORK(L(LF0)),DWORK(KF1),2*NU)
      CALL LLC1(DWORK(L(LD1)),DWORK(KF1),NU,-TSTEPH,1D0)
      CALL LLC1(DWORK(L(LD2)),DWORK(KF2),NU,-TSTEPH,1D0)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C-----------------------------------------------------------------------
C *** momentum solve without convection (ISTOK=1)
C     (IPRECA=4: SUPWDG/SUPWNP are called also for ISTOK=1, since they
C      assemble nu*A, so the convection is switched off by CNVFAC=0)
C-----------------------------------------------------------------------
C
      IF (IPRECA.EQ.4) CNVFAC=0D0
      CALL  NSDEF (MFILE,MSHOW,BSTOP,BNLEND)
      CNVFAC=CNVBAS
      IF (IER.NE.0) GOTO 99998
      IF (BSTOP) GOTO 99997
C
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
C-----------------------------------------------------------------------
C *** FP = 1/(theta*K) B^T u~  and Uzawa residual |B^T u~|/|u~|
C-----------------------------------------------------------------------
C
      CALL ZTIME(TTT0)
      CALL XDFKD(KB1,KB2,KCOLB,KLDB,KFP,KU1,KU2,NU,NP,TSTEPH)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL LL21 (DWORK(KU1),NU,DNRMU1)
      CALL LL21 (DWORK(KU2),NU,DNRMU2)
      DNORMU=SQRT(DNRMU1*DNRMU1+DNRMU2*DNRMU2)
      IF (ABS(DNORMU).LT.1D-8) DNORMU=1D0
      CALL LL21 (DWORK(KFP),NP,DSXN)
      DIVUZ=TSTEPH*DSXN/DNORMU
C
C-----------------------------------------------------------------------
C *** projection step of pp2d (as in PROSTP, always incremental):
C     solve B^T Ml^-1 B q = FP
C-----------------------------------------------------------------------
C
      CALL LCP1(DWORK(KP)      ,DWORK(L(LDP))  ,NP)
      CALL LCP1(DWORK(L(LPOLD)),DWORK(KP)      ,NP)
      CALL LCP1(DWORK(L(LDP))  ,DWORK(L(LPOLD)),NP)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      DO 30 ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+2*KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+2*KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+2*KNU(ILEV)
      KNEQ (ILEV)=KNP(ILEV)
      KPRSM(ILEV)=NSMP*NSMPFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMP*NSMPFA**(NLMAX-ILEV)
30    CONTINUE
C
      ICYCLE=ICYCP
      IRELMG=1
      ITMGP =ILMINP
      IF (ILMAXP.GT.ILMINP) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXP,ITMGP,DMPPMG,DNORMU*EPSP/TSTEPH,DEFPMG,
     *            YAXP,YPROLP,YRESTP,YSMP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMGP,BMGP)
      NMGP=NMGP+ITMGP
C
      TTMGP=TTMGP+TTMG
      TTSP =TTSP+TTS
      TTEP =TTEP+TTE
      TTDP =TTDP+TTD
      TTPP =TTPP+TTP
      TTRP =TTRP+TTR
C
      IF (IER.GT.0) IER=0
      IF (IER.NE.0) GOTO 99998
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
      IF (BMGP) THEN
       BSTOP=.TRUE.
       GOTO 99997
      ENDIF
C
C-----------------------------------------------------------------------
C *** u = u~ - theta*K*Ml^-1 B q
C-----------------------------------------------------------------------
C
      CALL ZTIME(TTT0)
      CALL XDFKG(KM1,KB1,KB2,KCOLB,KLDB,L(LD1),L(LD2),KU1,KU2,KP,NU,NP,
     *           KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),TSTEPH)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C-----------------------------------------------------------------------
C *** p_k+1 = p_k + PRDIF1*q + PRDIF2*alpha*theta*K*nu*Mp^-1*FP
C     (reactive + diffusive part of the Cahouet-Chabard precond.)
C-----------------------------------------------------------------------
C
      CALL ZTIME(TTT0)
      CALL LCP1(DWORK(KFP)     ,DWORK(KAUXP)   ,NP)
C
      CALL LCP1(DWORK(L(LPOLD)),DWORK(KFP)     ,NP)
      CALL LCP1(DWORK(KP)      ,DWORK(L(LDP))  ,NP)
      CALL LLC1(DWORK(L(LPOLD)),DWORK(KP)      ,NP,1D0,PRDIF1)
      CALL LCP1(DWORK(L(LDP))  ,DWORK(L(LPOLD)),NP)
C
      APDIV=PRDIF2*TMSTEP
      IF (APDIV.NE.0D0)
     *    CALL XDFKDV(KP,KAUXP,VWORK(L(KLAREA(NLEV))),NP,APDIV)
C
C *** pressure change p_k+1 - p_k (for output)
      CALL LCP1(DWORK(KFP),DWORK(L(LDP)),NP)
      CALL LLC1(DWORK(KP) ,DWORK(L(LDP)),NP,-1D0,1D0)
      CALL LL21(DWORK(L(LDP)),NP,DSXN)
      RELP2=DSXN/(SQRT(DBLE(NP))*TSTEPH)
      CALL LLI1(DWORK(L(LDP)),NP,DSXN,INDMAX)
      RELPM=DSXN/TSTEPH
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,10003) IUZ,DIVUZ,ITMGP,RELP2,RELPM,
     *                                   RHOMGP
      IF (MSHOW.GE.1) WRITE(MFILE,10003) IUZ,DIVUZ,ITMGP,RELP2,RELPM,
     *                                   RHOMGP
C
      IF (DIVUZ.LE.EPSUZ) GOTO 500
C
200   CONTINUE
      IUZ=NUZMAX
C
C=======================================================================
C *** Output for this substep: iterations and divergence of u
C=======================================================================
C
500   CALL ZTIME(TTT0)
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
      CALL XDFKD(KB1,KB2,KCOLB,KLDB,L(LDIV),KU1,KU2,NU,NP,1D0)
      CALL LL21(DWORK(L(LDIV)),NP,DIVL2)
      CALL LLI1(DWORK(L(LDIV)),NP,DIVMX,INDMAX)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      INLS=NNONL-NNONL0
      IMGUS=NMGU-NMGU0
      IMGPS=NMGP-NMGP0
C
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,10002) ITNSL,ITNS,TIMENS,RELP2,RELPM
      IF (MSHOW.GE.2) WRITE(MTERM,10004) ITNSL,CSTEP,INLS,IMGUS,IUZ,
     *                                   IMGPS,DIVUZ,DIVL2,DIVMX
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,*)
C
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,10002) ITNSL,ITNS,TIMENS,RELP2,RELPM
      IF (MSHOW.GE.1) WRITE(MFILE,10004) ITNSL,CSTEP,INLS,IMGUS,IUZ,
     *                                   IMGPS,DIVUZ,DIVL2,DIVMX
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
      THSTEP=TSTEPH
      GOTO 99998
C
C *** unexpected STOP in a solver: report it (PROSTP returns at once
C     and, for IADTIM=0, the rest of the macro step is skipped)
99997 IF (MSHOW.GE.1) WRITE(MTERM,10005) ITNSL,CSTEP,TIMENS
      IF (MSHOW.GE.1) WRITE(MFILE,10005) ITNSL,CSTEP,TIMENS
      CNVFAC=CNVBAS
      VISFAC=1D0
C
C *** restore the temporarily changed parameters
99998 ISTOK =ISTOKH
      INLMIN=INLMIH
      INLMAX=INLMAH
      ITEXL =ITEXLH
      IFRZW =0
C
      GOTO 99999
C
C
C
   3  FORMAT(80('+'))
10002 FORMAT ('#',I4,'(',I4,')',1X,'TIME=',D10.3,1X,'REL2(P)=',
     *        D10.3,1X,'RELM(P)=',D10.3)
10003 FORMAT('  UZAWA',I4,' DIV(U~)=',D10.3,' MGP=',I4,
     *       ' REL2(DP)=',D10.3,' RELM(DP)=',D10.3,' RHOMGP=',D9.2)
10004 FORMAT('GLOW',I2,1X,A7,' NL=',I3,' MGU=',I4,' UZ=',I3,
     *       ' MGP=',I4,' DIV(U~)=',D9.2,' DIV-L2=',D9.2,
     *       ' DIV-MX=',D9.2)
10005 FORMAT('GLOW',I2,1X,A7,' *** BSTOP (solver failure) at TIME=',
     *       D10.3,', substep aborted')
C
C
C
99999 END
