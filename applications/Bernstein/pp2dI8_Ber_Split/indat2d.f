************************************************************************
      DOUBLE PRECISION FUNCTION FDATIN(ITYP,IBLOC,X,Y,TIMENS,RE)
*
*     Prescribed data for files coeff.f and bndry.f
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      IMPLICIT INTEGER*8 (I-N)
      PARAMETER (PI=3.1415926535897931D0)
C *** time factor of the inflow (IRAMP=1,2, TRAMP, read in GDAT)
      COMMON /NSRAMP/ TRAMP,IRAMP
      SAVE /NSRAMP/
C
C
C
      FDATIN=0D0
C
C *** time factor of the inflow: 1, or sin^2(pi/2*t/TRAMP) for t<TRAMP
C     (smooth start: u(0)=0 and du/dt(0)=0 are compatible with u0=0)
      RAMP=1D0
      IF ((IRAMP.EQ.1).AND.(TIMENS.LT.TRAMP))
     *  RAMP=SIN(0.5D0*PI*TIMENS/TRAMP)**2
C *** DFG benchmark 2D-3: U(t)=1.5 sin(pi*t/8) for t in [0,8], TRAMP=8
      IF (IRAMP.EQ.2) RAMP=SIN(PI*TIMENS/TRAMP)
C
C
C=======================================================================
C *** Case 1: Velocity boundary values and/or exact solution
C=======================================================================
C
      IF (ITYP.EQ.1) THEN
C
       IF (IBLOC.EQ.1) THEN
        IF (X.EQ.0D0) FDATIN=RAMP*4D0*1.5D0/0.1681D0*Y*(0.41D0-Y)
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 2: Velocity x-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.2) THEN
C
       IF (IBLOC.EQ.1) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 3: Velocity y-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.3) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 4: Exact pressure solution
C=======================================================================
C
      IF (ITYP.EQ.4) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 5: Right hand side for momentum equation
C=======================================================================
C
      IF (ITYP.EQ.5) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 6: Right hand side for continuity equation
C=======================================================================
C
      IF (ITYP.EQ.6) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 7: Mean pressure values
C=======================================================================
C
      IF (ITYP.EQ.7) THEN
       DPAR=X
       INPR=IBLOC
C
       IF ((DPAR.GT.1D0).AND.(DPAR.LT.2D0).AND.(INPR.EQ.1)) THEN
        FDATIN=0D0
       ENDIF
C
       IF ((DPAR.GT.3D0).AND.(DPAR.LT.4D0).AND.(INPR.EQ.1)) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
*
*     Neumann-boundary part
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C
C
C=======================================================================
C *** Case 0: Set number of Neumann-boundary parts
C=======================================================================
C
      IF (INPART.EQ.0) THEN
       INPART=1      
      ENDIF
C
C
C=======================================================================
C *** Case <>0: Specify Neumann-boundary parts
C=======================================================================
C
      IF (INPART.GT.0) THEN
C
       IF (INPART.EQ.1) THEN
        DPARN1=1D0
        DPARN2=2D0
        INPRN =1
       ENDIF
C       
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE PTSDAT(TIMENS,DNU)
*
*     Data for Point-output (for fpost and bdpres and bdforc)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      SAVE
C
      EXTERNAL UE
C
C     
C
C=======================================================================
C *** Points for velocity, pressure and flux-tracing
C=======================================================================
C
      KPU(1)=264
      KPU(2)=17
C
      KPP(1)=27
      KPP(2)=30
      KPP(3)=316
      KPP(4)=17
C
      KPX(1)=31
      KPX(2)=6
C
C
C=======================================================================
C *** Parameters for 2 integral pressures in bdpres.f
C=======================================================================
C
      KPI(1)  =2
      DPI(1,1)=0D0
      DPI(2,1)=1D0
C
      KPI(2)  =2
      DPI(1,2)=0D0
      DPI(2,2)=0.5D0
C
C
C=======================================================================
C *** Parameters for lift (DFW) and drag (DAW) in bdforc.f (INPR=2)
C ***
C *** dfw=2 int_s [dpf(1) dut/dn n_y - p n_x] ds / dpf(2)
C *** daw=2 int_s [dpf(1) dut/dn n_x + p n_y] ds / dpf(2)
C ***
C=======================================================================
C
      RHO  =1.0D0
      DIST =0.1D0
      UMEAN=1.0D0
C
      DPF(1)=RHO*DNU
      DPF(2)=RHO*DIST*UMEAN**2
C
C
C
99999 END
