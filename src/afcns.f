************************************************************************
*     Algebraic flux correction (AFC) of the convection for the
*     velocity: Q1~ (IELT=0..3) or Bernstein Q2 (IELT=4), IAFC=1 (dat
*     lines IAFC, CQAFC). IAFC=2: alpha=1 (consistency test).
*
*     SUPWDG accumulates the Galerkin convection matrix K (with its
*     weight THSTEP*THF, i.e. as added to A) in DOUBLE on the structure
*     of A (array DKAFC) and calls AFCADD, which
*       - builds the discrete diffusion d_ij = max(0,k_ij,k_ji) (i<>j)
*       - adds D (offdiag -d_ij, diag sum_j d_ij) to A (IDEF<2), so
*         that the matrices on all levels (multigrid preconditioner)
*         are the low-order operator L = K + D
*       - IDEF>0: replaces K u in the defect by L u - fbar(u), i.e.
*         D_i = D_i - sum_j (1-alpha_ij) d_ij (u_i - u_j), with the
*         Zalesak limiter of the scalar test (BNLIM, bnafc.f): local
*         bounds over the matrix neighbours, Q+- = CQAFC*sum_j d_ij*
*         (u_max/min - u_i), R+-=1 at Dirichlet dofs; each velocity
*         component is limited separately.
*     The mass matrix and nu*A stay Galerkin (consistent), the time
*     discretization (theta / fractional step) is unchanged: every
*     evaluation of the convection (implicit part, explicit right hand
*     side with u^n, OMEGA optimization) uses the same limited operator.
*     Restrictions (GDAT): IUPW=0, UPSAM=0, IJUMP=0, IREFIN=1,
*     IPRECA<4 (so that the element matrix of SUPWDG/SUPWNP is K only).
*     Splitting (ISPLIT>0, prostg.f): the explicit convection in the
*     right hand side of the Stokes substeps (weight theta*K) and the
*     Burgers substep (CNVFAC=CNVBAS/beta, transport velocity w, the
*     limiter acts on the unknown u) are limited the same way.
************************************************************************
      SUBROUTINE AFCINI
************************************************************************
*     Workspace for AFC (sizes of level NLMAX) and the transposition
*     index KT (position of a_ji for a_ij) of A on every level
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /NSAFC/  CQAFC,ALAFC,IAFC,NAFCL,NAFCC,
     *                LDKAF,LDDAF,LDQAF,LDRPAF,LDRMAF,LDFAF,LKBDAF,
     *                KLKTAF(NNLEV)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF (IAFC.LE.0) RETURN
      NAM=KNA(NLMAX)
      NUM=KNU(NLMAX)
      CALL ZNEW(NAM,1,LDKAF ,'DKAFC ')
      CALL ZNEW(NAM,1,LDDAF ,'DDAFC ')
      CALL ZNEW(NUM,1,LDQAF ,'DQAFC ')
      CALL ZNEW(NUM,1,LDRPAF,'DRPAFC')
      CALL ZNEW(NUM,1,LDRMAF,'DRMAFC')
      CALL ZNEW(NUM,1,LDFAF ,'DFAFC ')
      CALL ZNEW(NUM,3,LKBDAF,'KBDAFC')
      IF (IER.NE.0) RETURN
      DO 10 JLEV=NLMIN,NLMAX
       CALL ZNEW(KNA(JLEV),3,KLKTAF(JLEV),'KTAFC ')
       IF (IER.NE.0) RETURN
       CALL BNTRIX(KWORK(L(KLLDA(JLEV))),KWORK(L(KLCOLA(JLEV))),
     *             KNU(JLEV),KWORK(L(KLKTAF(JLEV))))
10    CONTINUE
      ALAFC=1D0
      NAFCL=0
      NAFCC=0
      END
C
C
C
************************************************************************
      SUBROUTINE AFCADD(A,NA,KCOLA,KLDA,DK,NEQ,U1,U2,D1,D2,IDEF,ICOMP)
************************************************************************
*     A (REAL): matrix of the current level ILEV, DK: accumulated
*     convection K (DOUBLE, weighted), NEQ velocity dofs per component.
*     Statistics (limited defects only): ALAFC = mean alpha of the
*     last call, NAFCL = number of alpha<0.999 there, NAFCC = calls
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      REAL A
      DIMENSION A(*),KCOLA(*),KLDA(*),DK(*),U1(*),U2(*),D1(*),D2(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /NSAFC/  CQAFC,ALAFC,IAFC,NAFCL,NAFCC,
     *                LDKAF,LDDAF,LDQAF,LDRPAF,LDRMAF,LDFAF,LKBDAF,
     *                KLKTAF(NNLEV)
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
      SAVE
C
      IF ((NEQ.NE.KNU(ILEV)).OR.(NA.NE.KNA(ILEV))) THEN
       WRITE(*,*) 'AFCADD: level sizes do not match',ILEV,NEQ,NA
       STOP
      ENDIF
      KT =L(KLKTAF(ILEV))
      KDD=L(LDDAF)
      KDQ=L(LDQAF)
C *** discrete diffusion, added to A (diagonal first in every row)
      CALL AFCDSC(DK,KLDA,KCOLA,KWORK(KT),NEQ,DWORK(KDD),DWORK(KDQ))
      IF (IDEF.LT.2) THEN
       DO 20 IR=1,NEQ
        A(KLDA(IR))=A(KLDA(IR))+REAL(DWORK(KDQ+IR-1))
        DO 10 IA=KLDA(IR)+1,KLDA(IR+1)-1
10      A(IA)=A(IA)-REAL(DWORK(KDD+IA-1))
20     CONTINUE
      ENDIF
      IF (IDEF.LE.0) RETURN
C
C *** Dirichlet flags
      KBD=L(LKBDAF)
      DO 30 I=1,NEQ
30    KWORK(KBD+I-1)=0
      KMBD=L(KLMBD(ILEV))
      DO 40 I=1,KNMBD(ILEV)
       IDF=KWORK(KMBD+I-1)
       IF ((IDF.GE.1).AND.(IDF.LE.NEQ)) KWORK(KBD+IDF-1)=1
40    CONTINUE
C
      ALS=0D0
      NLS=0
      CALL AFCDEF(U1,D1,DWORK(KDD),DWORK(KDQ),KLDA,KCOLA,NEQ,
     *            KWORK(KBD),IAFC,CQAFC,DWORK(L(LDFAF)),
     *            DWORK(L(LDRPAF)),
     *            DWORK(L(LDRMAF)),ALAVG,NLIM)
      ALS=ALS+ALAVG
      NLS=NLS+NLIM
      NC=1
      IF (ICOMP.EQ.2) THEN
       CALL AFCDEF(U2,D2,DWORK(KDD),DWORK(KDQ),KLDA,KCOLA,NEQ,
     *             KWORK(KBD),IAFC,CQAFC,DWORK(L(LDFAF)),
     *             DWORK(L(LDRPAF)),
     *            DWORK(L(LDRMAF)),ALAVG,NLIM)
       ALS=ALS+ALAVG
       NLS=NLS+NLIM
       NC=2
      ENDIF
      ALAFC=ALS/DBLE(NC)
      NAFCL=NLS
      NAFCC=NAFCC+1
      END
C
C
C
************************************************************************
      SUBROUTINE AFCDSC(DK,KLD,KCOL,KT,N,DD,DQS)
************************************************************************
*     DD(ia) = d_ij = max(0,k_ij,k_ji) (0 on the diagonal),
*     DQS(i) = sum_j d_ij
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION DK(*),KLD(*),KCOL(*),KT(*),DD(*),DQS(*)
      DO 20 IR=1,N
       S=0D0
       DD(KLD(IR))=0D0
       DO 10 IA=KLD(IR)+1,KLD(IR+1)-1
        D=MAX(0D0,DK(IA),DK(KT(IA)))
        DD(IA)=D
        S=S+D
10     CONTINUE
       DQS(IR)=S
20    CONTINUE
      END
C
C
C
************************************************************************
      SUBROUTINE AFCDEF(U,D,DD,DQS,KLD,KCOL,N,KBD,IAFC,CQ,DF,DRP,DRM,
     *                  ALAVG,NLIM)
************************************************************************
*     D_i = D_i - ( sum_j d_ij (u_i-u_j) - sum_j alpha_ij d_ij (u_i-u_j) )
*     (defect with K u replaced by L u - fbar)
*     IAFC=2: alpha=1 (test: Galerkin defect, low-order preconditioner)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      DIMENSION U(*),D(*),DD(*),DQS(*),KLD(*),KCOL(*),KBD(*)
      DIMENSION DF(*),DRP(*),DRM(*)
      IF (IAFC.EQ.2) THEN
       ALAVG=1D0
       NLIM=0
       RETURN
      ENDIF
      CALL BNLIM(U,DD,DQS,KLD,KCOL,N,KBD,CQ,DF,DRP,DRM,ALAVG,NLIM)
      DO 20 IR=1,N
       S=0D0
       UI=U(IR)
       DO 10 IA=KLD(IR)+1,KLD(IR+1)-1
10     S=S+DD(IA)*(UI-U(KCOL(IA)))
C *** BNLIM sets DF=0 in Dirichlet rows (their defect is reset anyway)
       IF (KBD(IR).NE.1) D(IR)=D(IR)-(S-DF(IR))
20    CONTINUE
      END
