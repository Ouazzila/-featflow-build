************************************************************************
      SUBROUTINE WRWGT (CFILE)
************************************************************************
*
*   Purpose: - writes the weights for mass-weighted L2 norms on the
*              finest level (formatted, for external postprocessing):
*                NU NP
*                M_l(i), i=1..NU   lumped velocity mass matrix (KM1)
*                |T_k|,  k=1..NP   element areas (KLAREA, REAL)
*              so that |u|_L2^2 = sum M_l(i)*(u1(i)^2+u2(i)^2)
*                  and |p|_L2^2 = sum |T_k|*p(k)^2
*
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER CFILE*(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
      SAVE
C
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** finest level must be active (SETLEV called by the caller)
      OPEN (UNIT=38,FILE=CFILE,FORM='FORMATTED')
      WRITE(38,'(2I12)') NU,NP
      WRITE(38,'(ES25.16)') (DWORK(KM1+I-1),I=1,NU)
      LAREA=L(KLAREA(NLEV))
      WRITE(38,'(ES25.16)') (DBLE(VWORK(LAREA+I-1)),I=1,NP)
      CLOSE(38)
C
      END
