      PROGRAM ROTIN3 
c
c -------------------------------------------------------------------
c Program for performing rotational and instrumental convolution
c for a calculated stellar spectrum obtained by program SYNSPEC
c -------------------------------------------------------------------
c
      PARAMETER (MLAM=10000000,
     *           MCON=2000)
      PARAMETER (MCONV=500000)
      implicit real*8 (a-h,o-z)
      DIMENSION FLAM0(MLAM),FLAM1(MLAM),FLAM2(MLAM),FCN(MLAM),
     *          WLAM0(MLAM),WLAM1(MLAM),WLAM2(MLAM),
     *          FCON0(MCON),WCON0(MCON)
      real*4    DUMMY(2)
      character*20 fname7,fname17,fnout
c     EQUIVALENCE (FLAM2(1),FLAM0(1)),
c    *            (WLAM2(1),WLAM0(1)),
c    *            (FCN(1),FLAM1(1)) 
      T0=0.
c
c     ---------------------------------------------------------------
C     INPUT - from unit 5 - four lines of input
c     ---------------------------------------------------------------
c
c     1. filenames:
c
c     fname7  - name of the file containing the detailed synthetic spectrum
c               (i.e. fort.7 produced by Synspec )
c
c     fname17 - name of the file containing the continuum flux
c               (i.e. fort.17 produced by Synspec )
c
c     fnout   - name of the output file - convolved spectrum
c
c
c     2. parameters for rotational convolution 
c
c     VROT  - v sin i (in km/s)
c             if VROT=0 - rotational convolution is 
c                 a) either not calculated,
c                 b) or, if simultaneously FWHM is rather large
c                    (vrot/c*lambda < FWHM/20.),
c                    vrot is set to  FWHM/20*c/lambda;
c             if VROT >0 but the previous condition b) applies, the
c                     value of VROT is changed as  in the previous case
c             if VROT<0 - the value of abs(VROT) is used regardless of
c                     how small compared to FWHM it is
c     CHARD - characteristic scale of the variations of unconvolved
c             stellar spectrum (basically, characteristic distance
c             between two neighbouring wavelength points) - in A
c           - if =0 - program sets up default (0.01 A)
c     STEPR - wavelength step for evaluation rotational convolution;
c           - if =0, the program sets up default (the wavelength
c                    interval corresponding to the rotational velocity
c                    devided by 3.)                           
c             if <0, convolved spectrum calculated on the original
c             (detailed) SYNSPEC wavelength mesh
c
c
c     3. parameters for instrumental convolution
c
c     FWHM  - full width at half maximum for Gaussian instrumental 
c             profile
c     STEPI - wavelength step for evaluating instrumental convolution
c           - if =0, the program sets up default (FWHM/10.)
c           - if <0, convolved spectrum calculated with the previous
c                    wavelength mesh:
c                    either the original (SYNSPEC) one if vrot=0,
c                    or the one used in rotational convolution (vrot > 0)
c
c
c     4. wavelength interval and normalization of spectra
c
c     ALAM0 - initial wavelength
c     ALAM1 - final wavelength
c     IREL  - for =1 relative spectrum
c                 =0 absolute spectrum
c
      read(5,*) fname7,fname17,fnout
      read(5,*) vrot,chard,stepr
      read(5,*) fwhm,stepi
      read(5,*) alam0,alam1,irel
c
      open(unit=37,file=fname7,status='old')
      open(unit=47,file=fname17,status='old')
      open(unit=11,file=fnout)
c
c     if some of the above convolution are zero, set up defaults
c
      SLAM=0.5*(ALAM0+ALAM1)
      DROT=SLAM*ABS(VROT)/3.E5
      dins=fwhm/20.
      if(chard.le.0.) chard=0.01
      if(drot.lt.dins.and.vrot.ge.0.) then
         drot=dins
         vrot=drot*3.e5/slam
      end if
      if(vrot.lt.0.) vrot=abs(vrot)
      if(stepr.eq.0.) stepr=max(drot/5.,chard)
      if(stepi.eq.0.) stepi=max(fwhm/10.,chard)
      if(stepi.lt.stepr) stepi=stepr
         
      write(6,601) vrot,chard,stepr,fwhm,stepi,alam0,alam1,irel
  601 format(' vrot            =',f10.3/
     *       ' chard           =',f10.3/
     *       ' step(rot.conv)  =',f10.3/
     *       ' fwhm            =',f10.3/
     *       ' step(inst.conv) =',f10.3/
     *       ' lambda initial  =',f10.3/
     *       ' lambda final    =',f10.3/
     *       ' irel            =',i6/)
c
c -------------------------------------------------------------------
c    read in and transform the output from SYNSPEC
c -------------------------------------------------------------------
c
c  WLAM0 - array of wavelengths
c  FLAM0 - array of fluxes
c  WCON0 - arrray of wavelengths for continuum
c  FCON0 - array of continuum fluxes
c
      I=0
      ILAM=0
      ICON=0
      ALM=0.
   10 READ(37,*,END=20) AL,FL  
      IF(AL.LT.ALAM0) GO TO 10
      IF(AL.GT.ALAM1) GO TO 20
      ILAM=ILAM+1
      IF(ILAM.GT.MLAM) GO TO 20
      WLAM0(ILAM)=AL
      FLAM0(ILAM)=FL
      GO TO 10
   20 CONTINUE 
   30 READ(47,*,END=40) AL,FL  
      IF(AL.LT.ALAM0-10.) GO TO 30
      IF(AL.GT.ALAM1+10.) GO TO 40
      IF(AL.LE.ALM) GO TO 30
      ICON=ICON+1
      IF(ICON.GT.MCON) GO TO 40
      WCON0(ICON)=AL
      FCON0(ICON)=FL
      ALM=AL
      GO TO 30
   40 CONTINUE 
      NLAM0=ILAM
      NCON0=ICON
      IF(NLAM0.GT.MLAM) NLAM0=MLAM
      IF(NCON0.GT.MCON) NCON0=MCON
      IF(NLAM0.LE.0) GO TO 1000
      IF(NCON0.LE.0.AND.IREL.GT.0) GO TO 1000
c
c     if reguired, normalization to get a relative spectrum
c    
      if(irel.eq.1) then
         call interp(wcon0,fcon0,wlam0,fcn,ncon0,nlam0,2,0,0)
         do 50 i=1,nlam0
   50       flam0(i)=flam0(i)/fcn(i)
      end if
c
c     timing (only for Ultrix)
c
      time=etime(dummy)
      dt=time-t0
      t0=time
      write(6,602) t0,dt
  602 format(' TIME - READ SPECTRUM : ',2f10.2)
c
c -----------------------------------------------------------------
c Rotational convolution
c (no rotational convolution if vrot=0)
c -----------------------------------------------------------------
c
      if(vrot.gt.0.) then
         SLAM=0.5*(WLAM0(1)+WLAM0(NLAM0))
         XN=SLAM*VROT/3.E5/CHARD
         NROT=INT(XN)
         IF(NROT.LT.10) NROT=10
         IF(NROT.GT.MCONV) then
            write(6,*) 'nrot too large'
            stop
         end if
         if(stepr.le.0.) then
            nlam1=nlam0
            do 80 i=1,nlam1
   80          wlam1(i)=wlam0(i)
          else
            XN=(ALAM1-ALAM0)/STEPR
            NLAM1=INT(XN)+1
            DO 90 I=1,NLAM1
   90          WLAM1(I)=ALAM0+STEPR*(I-1)
         end if
         CALL ROTINS(1,FLAM0,FLAM1,WLAM0,WLAM1,NLAM0,NLAM1,
     *               NROT,VROT,0.D0)
c
       else
c
c  no rotational convolution if VROT=0
c
         nlam1=nlam0
         do 100 i=1,nlam1
            wlam1(i)=wlam0(i)
            flam1(i)=flam0(i)
  100    continue
      end if
c
c     timing (only for Ultrix)
c
      time=etime(dummy)
      dt=time-t0
      t0=time
      write(6,603) t0,dt,nrot
  603 format(' TIME - ROTATION CONV : ',2f10.2,'    integr.points:',i5)
c
c ------------------------------------------------------------------
c  Instrumental convolution
c  no instrumental convolution for FWHM=0
c ------------------------------------------------------------------
c
      if(fwhm.gt.0.) then
         IF(VROT.LE.0.) THEN
            XN=6.*FWHM/CHARD
          ELSE
            XN=MIN(1.8E6*FWHM/VROT/SLAM,6.*FWHM/CHARD)
         END IF
         NINS=INT(XN)
         IF(NINS.LT.10) NINS=10
c         write(6,*) 'number of integration points - instrum.',nins
         IF(NINS.GT.MCONV) then
            write(6,*) 'nins too large'
            stop
         end if
         if(stepi.le.0.) then
            nlam2=nlam1
            do 110 i=1,nlam2
  110          wlam2(i)=wlam1(i)
          else
            XNI=(ALAM1-ALAM0)/STEPI
            NLAM2=INT(XNI)+1
            DO 120 I=1,NLAM2
  120          WLAM2(I)=ALAM0+STEPI*(I-1)
         end if
         CALL ROTINS(2,FLAM1,FLAM2,WLAM1,WLAM2,NLAM1,NLAM2,
     *               NINS,0.D0,FWHM)
      else
         nlam2=nlam1
	 do 130 i=1,nlam2
            wlam2(i)=wlam1(i)
  130       flam2(i)=flam1(i)
      end if 
c
c     timing (only for Ultrix)
c
      time=etime(dummy)
      dt=time-t0
      t0=time
      write(6,604) t0,dt,nins
  604 format(' TIME - INSTRUM. CONV : ',2f10.2,'    integr.points:',i5)
c
c -----------------------------------------------------------------
c
c  output of final, convolved, synthetic spectrum on standard output
c  in a simple format:
c
c  values of wavelengths (in A) -   WLAM2 versus
c  corresponding values of fluxes - FLAM2
c
      do 200 i=1,nlam2
         write(11,600) wlam2(i),flam2(i)
  200 continue
  600 format(f16.9,1pe16.9)
c
c     rewind 7
c     rewind 17
c
 1000 continue
c
c     timing (only for Ultrix)
c
      time=etime(dummy)
      dt=time-t0
      t0=time
      write(6,605) t0,dt
  605 format(' TIME - FINAL         : ',2f10.2)
c
      stop
      end

c
c ********************************************************************
c

      SUBROUTINE ROTINS(MODE,HINP,HOUT,XLAM,YLAM,NLAMX,NLAMY,
     *                  NR,VROT,FWHM)
C
C ---------------------------------------------------------------
C
C Rotational and/or instrumental convolution
C
C  MODE  - for MODE = 1 - rotational convolution;
C          for MODE = 2 - instrumental convolution,
C                         with a Gaussian instrumental profile
C  HINP  - array input flux (arbitrary units)
C  XLAM  - array of input wavelengths (in angstroms)
C  HOUT  - array of output (convolved) flux
C  YLAM  - array of wavelengths in which the output flux is calculated
C  NLAMX - number of input wavelengths
C  NLAMY - number of output wavelengths
C  NR    - number of integration points for evaluating the
C          convolution integrals
C  VROT  - v sin i [in km/s]
C  FWHM  - full width at half maximum of the instrum. profile function
C ---------------------------------------------------------------------
C
      implicit real*8 (a-h,o-z)
      PARAMETER (MCONV=500000,
     *           MCONV1=MCONV+1)
      PARAMETER (MLAM=10000000)
      PARAMETER (ONE=1.,
     *           TWO=2.,
     *           HALF=0.5)
      PARAMETER (DLROT=10.)
      DIMENSION HINP(1),HOUT(1),XLAM(1),YLAM(1),G(MCONV1),
     *          IGCALC(MLAM)
C
      NR1=NR+1
      DO 10 I=1,NLAMY
         IGCALC(I)=0
   10 CONTINUE
C
      IF(MODE.EQ.1) THEN
         DLAM=YLAM(NLAMY)-YLAM(1)
         IF(DLAM.LE.DLROT) THEN
            SLAM=HALF*(YLAM(1)+YLAM(NLAMY))
          ELSE
            NCALG=DLAM/DLROT+1
            NSTEP=NLAMY/NCALG
            DO 20 I=1,NCALG
               IGCALC(I*NSTEP+1)=1
   20       CONTINUE
            SLAM=YLAM(1)
         END IF
       ELSE IF(MODE.EQ.2) THEN
      END IF
c
c  initial kernel function (rotation);
c  or the general kernel function (for instrumental)
c
      CALL KERNEL(MODE,NR,VROT,SLAM,FWHM,XLMAX,G)
      DLAM=XLMAX/NR
c
c determine the exterior intervals
c  a) beginning of the interval
c
      INTR0=0
      IEND0=0
      X0=XLAM(1)
      X1=X0+XLMAX
      HM0=HINP(1)
      DO 50 I=1,NLAMY
         IF(YLAM(I).LE.X0) THEN
            IEND0=I
            HOUT(I)=HM0
          ELSE IF(YLAM(I).LE.X1) THEN
            INTR0=I
          ELSE
            GO TO 60
         END IF
   50 CONTINUE
   60 CONTINUE
c
c  b) end of the interval
c
      INTR1=NLAMY
      IEND1=NLAMY
      X0=XLAM(NLAMX)
      IF(MODE.EQ.1.AND.YLAM(NLAMY)-YLAM(1).GT.DLROT) THEN
         XLMAX=YLAM(NLAMY)*VROT/2.997925E5
      END IF
      X1=X0-XLMAX
      HP0=HINP(NLAMX)
      DO 70 I=NLAMY,1,-1
         IF(YLAM(I).GE.X0) THEN
            IEND1=I
            HOUT(I)=HP0
          ELSE IF(YLAM(I).GE.X1) THEN
            INTR1=I
          ELSE
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
C
C ------------------------------------------------------------
C wavelength by wavelength convolution; integral calculated by
C the trapezoidal rule
C ------------------------------------------------------------
C  
C 1. points near the beginning of the interval
C
      IF(INTR0.GE.1) THEN
      K0=1
      DO 390 I=IEND0+1,INTR0
         HOUT(I)=0.
         DO 310 K=K0,NLAMX
            K0=K
            IF(XLAM(K).GT.YLAM(I)) GO TO 320
  310    CONTINUE
  320    K0=K0-1
         DO 350 J=1,NR1
            A2=(J-1)*DLAM
            ALAM=YLAM(I)+A2
            K=K0+1
  330       CONTINUE
            IF(ALAM.LT.XLAM(K)) THEN
               HPLUS=HINP(K-1)+(HINP(K)-HINP(K-1))/(XLAM(K)-XLAM(K-1))*
     *               (ALAM-XLAM(K-1))
               HOUT(I)=HOUT(I)+HPLUS*G(J)
             ELSE IF(ALAM.EQ.XLAM(K)) THEN
               HOUT(I)=HOUT(I)+HINP(K)*G(J)
             ELSE
               K=K+1
               GO TO 330
            END IF
  350    CONTINUE 
C  
         DO 380 J=1,NR1
            A2=(J-1)*DLAM
            ALAM=YLAM(I)-A2
            K=K0
  360       CONTINUE
            IF(ALAM.GT.XLAM(K)) THEN
               HMINUS=HINP(K)+(HINP(K+1)-HINP(K))/(XLAM(K+1)-XLAM(K))*
     *                (ALAM-XLAM(K))
               HOUT(I)=HOUT(I)+HMINUS*G(J)
             ELSE IF(ALAM.EQ.XLAM(K)) THEN
               HOUT(I)=HOUT(I)+HINP(K)*G(J)
             ELSE
               K=K-1
               IF(K.LT.1) THEN
                  HOUT(I)=HOUT(I)+HM0*G(J)
                  GO TO 380
               END IF
               GO TO 360
            END IF
  380    CONTINUE
         IF(K0.LE.0) K0=1
  390 CONTINUE
      END IF
C  
C 2. inner points
C
  601 format(5i5)
      if(intr0.le.0) intr0=1
      K0=1
      DO 300 I=INTR0+1,INTR1-1
C
C        re-evaluate the kernel function if necessary
c
         IF(IGCALC(I).EQ.1) THEN
            SLAM=YLAM(I)
            CALL KERNEL(MODE,NR,VROT,SLAM,FWHM,XLMAX,G)
            DLAM=XLMAX/NR
         END IF
c
c        perform the convolution integral
c
         HOUT(I)=0.
         DO 110 K=K0,NLAMX
            K0=K
            IF(XLAM(K).GT.YLAM(I)) GO TO 120
  110    CONTINUE
  120    K0=K0-1
         DO 200 J=1,NR1
            A2=(J-1)*DLAM
            ALAM=YLAM(I)+A2
            K=K0+1
  130       CONTINUE
            IF(ALAM.LT.XLAM(K)) THEN
               HPLUS=HINP(K-1)+(HINP(K)-HINP(K-1))/(XLAM(K)-XLAM(K-1))*
     *               (ALAM-XLAM(K-1))
               HOUT(I)=HOUT(I)+HPLUS*G(J)
             ELSE IF(ALAM.EQ.XLAM(K)) THEN
               HOUT(I)=HOUT(I)+HINP(K)*G(J)
             ELSE
               K=K+1
               GO TO 130
            END IF
  200    CONTINUE 
C  
         DO 220 J=1,NR1
            A2=(J-1)*DLAM
            ALAM=YLAM(I)-A2
            K=K0
  210       CONTINUE
            IF(ALAM.GT.XLAM(K)) THEN
               HMINUS=HINP(K)+(HINP(K+1)-HINP(K))/(XLAM(K+1)-XLAM(K))*
     *                (ALAM-XLAM(K))
               HOUT(I)=HOUT(I)+HMINUS*G(J)
             ELSE IF(ALAM.EQ.XLAM(K)) THEN
               HOUT(I)=HOUT(I)+HINP(K)*G(J)
             ELSE
               K=K-1
               GO TO 210
            END IF
  220    CONTINUE
         IF(K0.LE.0) K0=1
  300 CONTINUE
C  
C 3. points near the end of the interval
C
      IF(INTR1.LT.NLAMY) THEN
      IF(MODE.EQ.1.AND.DLAM.GT.DLROT) THEN
         SLAM=YLAM(NLAMY)
         CALL KERNEL(MODE,NR,VROT,SLAM,FWHM,XLMAX,G)
         DLAM=XLMAX/NR
      END IF
      K0=NLAMX
      DO 500 I=IEND1,INTR1,-1
         HOUT(I)=0.
         DO 410 K=K0,1,-1
            K0=K
            IF(XLAM(K).LT.YLAM(I)) GO TO 420
  410    CONTINUE
  420    CONTINUE
         DO 450 J=1,NR1
            A2=(J-1)*DLAM
            ALAM=YLAM(I)+A2
            K=K0+1
  430       CONTINUE
            IF(ALAM.LT.XLAM(K)) THEN
               HPLUS=HINP(K-1)+(HINP(K)-HINP(K-1))/(XLAM(K)-XLAM(K-1))*
     *               (ALAM-XLAM(K-1))
               HOUT(I)=HOUT(I)+HPLUS*G(J)
             ELSE IF(ALAM.EQ.XLAM(K)) THEN
               HOUT(I)=HOUT(I)+HINP(K)*G(J)
             ELSE
               K=K+1
               IF(K.GT.NLAMX) THEN
                  HOUT(I)=HOUT(I)+HP0*G(J)
                  GO TO 450
               END IF
               GO TO 430
            END IF
  450    CONTINUE 
C  
         DO 470 J=1,NR1
            A2=(J-1)*DLAM
            ALAM=YLAM(I)-A2
            K=K0
  460       CONTINUE
            IF(ALAM.GT.XLAM(K)) THEN
               HMINUS=HINP(K)+(HINP(K+1)-HINP(K))/(XLAM(K+1)-XLAM(K))*
     *                (ALAM-XLAM(K))
               HOUT(I)=HOUT(I)+HMINUS*G(J)
             ELSE IF(ALAM.EQ.XLAM(K)) THEN
               HOUT(I)=HOUT(I)+HINP(K)*G(J)
             ELSE
               K=K-1
               GO TO 460
            END IF
  470    CONTINUE
         IF(K0.LE.0) K0=1
  500 CONTINUE
      END IF
      RETURN
      END

C
C
C     ****************************************************************
C
C
      SUBROUTINE KERNEL(MODE,NR,VROT,SLAM,FWHM,XLMAX,G)
C     -------------------------------------------------
C
C Kernel function for the rotational and/or instrumental convolution 
C
C Input:
C
C  MODE  - for MODE = 1 - rotational convolution;
C          for MODE = 2 - instrumental convolution,
C                         with a Gaussian instrumental profile
C  NR    - number of integration points for evaluating the
C          convolution integrals
C  VROT  - v sin i [in km/s]
C  SLAM  - standard wavelength (for rotational convolution;
C          no meaning for the instrumental convolution)
C  FWHM  - full width at half maximum of the instrum. profile function
C
C  Output:
C
C   XLMAX- width of the kernel
C   G    - normalized kernel function
C
      implicit real*8 (a-h,o-z)
      PARAMETER (MCONV=500000, MCONV1=MCONV+1)
      PARAMETER (ONE=1.,
     *           TWO=2.,
     *           HALF=0.5)
      PARAMETER (EPS=0.6,
     *           GAUSLM=3.,
     *           DLROT=10.)
      DIMENSION G(MCONV1)
C
C set up integration limits for the convolution integral
C
      IF(MODE.EQ.1) THEN
         XLMAX=SLAM*VROT/2.997925E5
         E1=TWO*(ONE-EPS)
         E2=1.5707964*EPS
         E3=ONE/3.1415927/XLMAX/(ONE-EPS/3.)
       ELSE IF(MODE.EQ.2) THEN
         XLMAX=GAUSLM*FWHM
         D0=0.60056*FWHM
         D1=0.5641896/D0
      END IF
C
C evaluation of the kernel function G
C
      DLAM=XLMAX/NR
      NR1=NR+1
      DO 10 J=1,NR1
         X=(J-1)*DLAM
         IF(MODE.EQ.1) THEN
            X1=ABS(1.-(X/XLMAX)**2)
            G(J)=(E1*SQRT(X1)+E2*X1)*E3
          ELSE IF(MODE.EQ.2) THEN
            G(J)=D1*EXP(-(X/D0)**2)
          ELSE IF(MODE.EQ.3) THEN
            G(J)=(ONE-X/FWHM)/FWHM
         END IF
   10 CONTINUE
C
C renormalization in order to have   integral(G) = 1
C
      SUM=0.
      DO 20 J=2,NR
         SUM=SUM+TWO*G(J)
   20 CONTINUE
      SUM=SUM+G(1)+G(NR1)
      SUM=ONE/DLAM/SUM
C
      DO 30 J=1,NR1
         G(J)=G(J)*SUM
   30 CONTINUE
c
c  multiply by integration weights for trapezoidal integration
c
      DO 40 J=1,NR1 
         G(J)=G(J)*DLAM
   40 CONTINUE
      G(1)=HALF*G(1)
      G(NR1)=HALF*G(NR1)
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE INTERP(X,Y,XX,YY,NX,NXX,NPOL,ILOGX,ILOGY)
C     ====================================================
C
C     General interpolation procedure of the (NPOL-1)-th order
C
C     for  ILOGX = 1  logarithmic interpolation in X
C     for  ILOGY = 1  logarithmic interpolation in Y
C
C     Input:
C      X    - array of original x-coordinates
C      Y    - array of corresponding functional values Y=y(X)
C      NX   - number of elements in arrays X or Y
C      XX   - array of new x-coordinates (to which is to be 
C             interpolated
C      NXX  - number of elements in array XX
C     Output:
C      YY   - interpolated functional values YY=y(XX)
C
      implicit real*8 (a-h,o-z)
      DIMENSION X(1),Y(1),XX(1),YY(1)
      EXP10(X0)=EXP(X0*2.30258509299405D0)
      IF(NPOL.LE.0.OR.NX.LE.0) GO TO 200
      IF(ILOGX.EQ.0) GO TO 30
      DO 10 I=1,NX
   10    X(I)=LOG10(X(I))
      DO 20 I=1,NXX
   20    XX(I)=LOG10(XX(I))
   30 IF(ILOGY.EQ.0) GO TO 50
      DO 40 I=1,NX
   40    Y(I)=LOG10(Y(I))
   50 NM=(NPOL+1)/2
      NM1=NM+1
      NUP=NX+NM1-NPOL
      DO 100 ID=1,NXX
         XXX=XX(ID)
         DO 60 I=NM1,NUP
            IF(XXX.LE.X(I)) GO TO 70
   60    CONTINUE
         I=NUP
   70    J=I-NM
         JJ=J+NPOL-1
         YYY=0.
         DO 90 K=J,JJ
            T=1.
            DO 80 M=J,JJ
               IF(K.EQ.M) GO TO 80
               T=T*(XXX-X(M))/(X(K)-X(M))
   80       CONTINUE
   90    YYY=Y(K)*T+YYY
         YY(ID)=YYY
  100 CONTINUE
      IF(ILOGX.EQ.0) GO TO 130
      DO 110 I=1,NX
  110    X(I)=EXP10(X(I))
      DO 120 I=1,NXX
  120    XX(I)=EXP10(XX(I))
  130 IF(ILOGY.EQ.0) RETURN
      DO 140 I=1,NX
  140    Y(I)=EXP10(Y(I))
      DO 150 I=1,NXX
  150    YY(I)=EXP10(YY(I))
      RETURN
  200    N=NX
         IF(NXX.GE.NX) N=NXX
      DO 210 I=1,N
         XX(I)=X(I)
  210    YY(I)=Y(I)
      RETURN
      END
C
C ********************************************************************
C ********************************************************************
