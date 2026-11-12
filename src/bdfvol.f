************************************************************************
      SUBROUTINE BDFVOL(DU1,DU2,DP,DU1O,DU2O,DTK,DALPHA,
     *                  KVERT,KMID,KVBD,KEBD,KMM,DCORVG,ELE,DFWV,DAWV,
     *                  DFWP,DAWP,DFWN,DAWN)
************************************************************************
*    Purpose:  Calculates drag (DFWV) and lift (DAWV) on boundary
*              component 2 by the volume integral (after cc2d_powerlaw,
*              VIBF_STG/STGFORCE), here in the nonstationary form
*
*      F_x = - [ (dt u1 + u.grad u1, a) + nu (grad u1, grad a)
*                - (p, da/dx) ]
*      F_y = - [ (dt u2 + u.grad u2, a) + nu (grad u2, grad a)
*                - (p, da/dy) ]
*
*              with a = 1 in the midpoints (dofs) of the edges on
*              component 2 and a = 0 in all other dofs, so the
*              integrand lives only in the elements at the cylinder.
*              dt u = (u - u_old)/DTK  (backward difference, omitted
*              if DTK=0), rho=1 for the inertia terms, nu = DPF(1).
*              Scaling as in BDFORC: c = 2 F / DPF(2).
*
*    DU1O,DU2O  velocity of the previous call (time t - DTK)
*    DALPHA     work vector of length NU (returns the test function a)
*    DFWP,DAWP  pressure part   +2/DPF(2) (p, da/dx), (p, da/dy)
*    DFWN,DAWN  viscous part    -2/DPF(2) nu (grad u_i, grad a)
*               (inertia part (dt u + u.grad u, a) = total - both)
*    ELE        nonparametric element EM30/EM31 (as in BDFORC)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      DIMENSION DU1(*),DU2(*),DP(*),DU1O(*),DU2O(*),DALPHA(*)
      DIMENSION DCORVG(2,*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*),KMM(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      SAVE
C
      EXTERNAL ELE
C
C
C
      DFWV=0D0
      DAWV=0D0
      DFWP=0D0
      DAWP=0D0
      DFWN=0D0
      DAWN=0D0
C
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
C
C=======================================================================
C *** test function a: 1 on the edges of component 2 (same boundary
C     loop as in BDFORC: edge KVBD(IVBD)-KVBD(IVBD+1) in KEBD(IVBD))
C=======================================================================
C
      DO 1 IU=1,NU
1     DALPHA(IU)=0D0
C
      DO 10 IVBD=1,NVBD
      IF (KVBD(IVBD).EQ.KMM(1,2)) GOTO 15
10    CONTINUE
C
      WRITE(MTERM,*) 'ERROR 1 IN BDFVOL'
      RETURN
C
15    IVBD1=IVBD
C
      DO 20 IVBD=IVBD1,NVBD
      IELB=KEBD(IVBD)
      IVT1=KVBD(IVBD)
      IF (IVT1.EQ.KMM(2,2)) THEN
       IVT2=KVBD(IVBD1)
      ELSE
       IVT2=KVBD(IVBD+1)
      ENDIF
C
      DO 21 II=1,4
      JJ=MOD(II,4)+1
      IF (((KVERT(II,IELB).EQ.IVT1).AND.(KVERT(JJ,IELB).EQ.IVT2)).OR.
     *    ((KVERT(II,IELB).EQ.IVT2).AND.(KVERT(JJ,IELB).EQ.IVT1)))
     *  DALPHA(KMID(II,IELB)-NVT)=1D0
21    CONTINUE
C
      IF (IVT1.EQ.KMM(2,2)) GOTO 25
20    CONTINUE
C
25    CONTINUE
C
C=======================================================================
C *** volume integral, 3x3 Gauss formula
C=======================================================================
C
      DO 30 IDER=1,NNDER
30    BDER(IDER)=.FALSE.
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      ICUB=8
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
      DTKI=0D0
      IF (DTK.GT.0D0) DTKI=1D0/DTK
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** only the elements at the cylinder contribute
      BCYL=.FALSE.
      DO 101 IDFL1=1,IDFL
101   IF (DALPHA(KDFG(IDFL1)).NE.0D0) BCYL=.TRUE.
      IF (.NOT.BCYL) GOTO 100
C
      DO 110 IVE=1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
110   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
      CALL ELE(0D0,0D0,-2)
C
      DO 120 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *  +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
C
C *** nonparametric element: evaluated in real coordinates
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
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
C
      DPCONT=DP(IEL)
C
C *** F = - [ (dt u + u.grad u, a) + nu (grad u, grad a) - (p, div a) ]
      AH1=(DTKI*DTU1+U1V*U1X+U2V*U1Y)*AV
     *   +DPF(1)*(U1X*AX+U1Y*AY)-DPCONT*AX
      AH2=(DTKI*DTU2+U1V*U2X+U2V*U2Y)*AV
     *   +DPF(1)*(U2X*AX+U2Y*AY)-DPCONT*AY
C
      DFWV=DFWV-AH1*OM
      DAWV=DAWV-AH2*OM
C *** pressure and viscous parts separately
      DFWP=DFWP+DPCONT*AX*OM
      DAWP=DAWP+DPCONT*AY*OM
      DFWN=DFWN-DPF(1)*(U1X*AX+U1Y*AY)*OM
      DAWN=DAWN-DPF(1)*(U2X*AX+U2Y*AY)*OM
C
120   CONTINUE
C
100   CONTINUE
C
      DFWV=2D0*DFWV/DPF(2)
      DAWV=2D0*DAWV/DPF(2)
      DFWP=2D0*DFWP/DPF(2)
      DAWP=2D0*DAWP/DPF(2)
      DFWN=2D0*DFWN/DPF(2)
      DAWN=2D0*DAWN/DPF(2)
C
C
C
99999 END
