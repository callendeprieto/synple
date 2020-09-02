      PROGRAM SYNSPEC
C
C =====================================================================I
C                                                                      I
C Program for evaluting synthetic spectra for a given model atmosphere I
C                                                                      I
C *****************                                                    I
C VERSION SYNSPEC53                                                    I
C *****************                                                    I
C                                                                      I
C Input: the same as input to TLUSTY or TLUSDISK - unit 5              I
C        additional 6 lines of input - unit 55 (proc. START and INIBL0)I
C        chemical composition - unit 56 (if a switch is on in unit 55) I
C        model atmosphere - unit 8  (procedures INPMOD or INKUR)       I
C        line list        - unit 19 (procedure INISET)                 I
C                                                                      I
C Output: diagnostic outprint - unit 6 (several procedures)            I
C         synthetic spectrum  - unit 7 (procedure OUTPRI)              I
C         flux in continuum   - unit 17 (procedure OUTPRI)             I
C         identification table- unit 12 (procedure INIBLA)             I
C         partial equiv.widths- unit 16 (procedure OUTPRI)             I
C         elapsed time        - unit 69 (procedure TIMING - UNIX only) I
C                                                                      I
C      -- if specific intensities are also calculated (set up by the   I
C         input on unit 55), there are two aditional output files:     I
C                                                                      I
C         specific intensities - unit 10                               I
C         specific intensities in continuum - unit 18                  I
C                                                                      I
C      -- in the iron-curtain option (IMODE=-2), there is another      I
C         output file:                                                 I
C         monochromatic opacities - unit 27                            I
C                                                                      I
C     ***  The contents of units 7 and 17 serve as an input to the     I
C          program ROTIN, which performs rotational and instrumental   I
C          ROTIN, which performs rotational and instrumental           I
C          convolutions, and sets up files for a plot.                 I
C                                                                      I
C Basic options: controlled by switch IMODE                            I
C IMODE    =  0 - normal synthetic spectrum                            I
C                 (ie. identification table + emergent flux)           I
C          =  1 - detailed profiles of a few individual lines          I
C          =  2 - emergent flux in the continuum (without the          I
C                 contribution of lines)                               I
C          = -1 - only identification table, ie. a list of lines which I
C                 contribute to opacity in a given wavelength          I
C                 region, together with their approximate equivalent   I
C                 widths. Synthetic spectrum is not calculated.        I
C          = -2 - the "iron curtain" option, ie. a monochromatic       I
c                 opacity for a homogeneous slab of a given T and n_e  I
C                                                                      I
C                                                                      I
C ==================================================================== I
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'LINDAT.FOR'
      include 'MODELP.FOR'
      include 'SYNTHP.FOR'
C
      OPEN(UNIT=12,STATUS='UNKNOWN')
      OPEN(UNIT=14,STATUS='UNKNOWN')
C
C     INITIALIZATION - INPUT OF BASIC PARAMETERS AND MODEL ATMOSPHERE
C
      CALL START
      if(ifeos.gt.0) imode=-3
      if(ibfac.gt.1) then
         LTE0=LTE
         LTE=.TRUE.
      END IF
      IF(IMODE.GE.-2.AND.IFEOS.EQ.0) THEN
         IF(INMOD.GT.0) CALL INPMOD
         IF(INMOD.EQ.0) CALL INKUR
         IF(ICHANG.NE.0) CALL CHANGE
         IF(IBFAC.GT.1) THEN
            CALL INPBF
            LTE=LTE0
         END IF
         IF(IFWIN.GT.1) CALL SETWIN
       ELSE
         CALL INGRID(0,inext,0)
      END IF
C
      CALL INIBL0
      CALL INIMOD
      CALL TINT
c
      IMODE0=IMODE
      IF(IMODE0.EQ.-4) IMODE=2
      igrd=0
    1 continue
c
      IF(IMODE0.LE.-3.and.ifeos.eq.0) CALL INIBL1(IGRD)
      IF(IFMOL.GT.0) then
         CALL MOLINI
c        call eospri
      end if
c
c     zero abundances for selected species (if required)
c
      if(imode0.le.-3) call abnchn(1)
c
      IBLANK=0
      NXTSET=0
      IF(IFMOL.GT.0.AND.IMODE.LT.2) THEN
         DO ILIST=1,NMLIST
            NXTSEM(ILIST)=0
            INACTM(ILIST)=0
            NLINMT(ILIST)=0
         END DO
      END IF
c
      if(ifeos.eq.0) then
      IF(IMODE.LT.2) CALL INILIN
C
      IF(IFMOL.GT.0.AND.IMODE.LT.2) THEN
         DO ILIST=1,NMLIST
            IF(IMODE.EQ.-3.AND.TEMP(1).LT.TMLIM(ILIST))
     *      CALL INMOLI(ILIST)
            IF(IMODE.GE.-2.and.imode.le.1) CALL INMOLI(ILIST)
         END DO
      END IF
      end if
c
    5 CONTINUE
c
C     ACTUAL CALCULATION OF THE SYNTHETIC SPECTRUM
C
      IF(IFEOS.NE.0) GO TO 30
   10 IBLANK=IBLANK+1
      IF(IFWIN.LE.0) THEN
         CALL RESOLV
         IF(IMODE0.LT.0) GO TO 20
         if(ifreq.le.10.and.inmod.le.1) then
            CALL RTECD
          else
            call RTE
         end if
       else
         CALL RESOLW
      end if
      CALL OUTPRI
   20 CONTINUE
      if((imode.ge.0.and.imode.ne.7).or.
     *   (imode.lt.0.and.iprin.ge.2)) then
         CALL IDTAB
         IF(IFMOL.GT.0) CALL IDMTAB
      end if
      IF(IBLANK.LT.NBLANK) GO TO 10
      IF(NXTSET.EQ.1.AND.IRLIST.EQ.0) THEN 
        IF(IMODE.LT.2) CALL INILIN
        GO TO 5
      END IF
      IF(IFMOL.GT.0.AND.IMODE.LT.2.AND.IRLIST.GT.0) THEN
         DO ILIST=1,NMLIST
            IF(NXTSEM(ILIST).EQ.1.and.inactm(ilist).eq.0) THEN
               CALL INMOLI(ILIST)
               iblank=0
               GO TO 5
            END IF
         END DO
      END IF
   30 CONTINUE
c
      if(imode0.lt.-2) then
         call ingrid(1,inext,igrd)
         igrd=igrd+1
c        call timing(1,igrd)
         if(inext.gt.0) go to 1
      end if
      if(imode0.le.-3.and.ifeos.eq.0) call fingrd
      call timing(2,iblank)
      END
C
C ********************************************************************
C
C
C

      SUBROUTINE START
C     ================
C
C     General input and initialization procedure
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      common/quasun/nunalp,nunbet,nungam,nunbal
C
C ------------------------------------------------
C Additional basic input parameters - from unit 55
C ------------------------------------------------
C
C     IMODE     =  0 - normal synthetic spectrum
C               =  1 - detailed profiles of a few individual lines
C               =  2 - emergent flux in the continuum (without the
C                      contribution of lines)
C               = -1 - identification table, ie. a list of lines which
C                      contribute to opacity in a given wavelength
C                      region, together with their approximate equivalent
C                      widths. Synthetic spectrum is not calculated.
C               = -2 - the "iron curtain" option, ie. a monochromatic 
C                      opacity for a homogeneous slab of a given T and n_e
C
C     IDSTD      - index of the "standard depth" (ie the depth at which
C                  the continuum optical depth is of the order of unity)
C                  (for detailed explanation see the code TLUSTY)
C
C     IPRIN      - determines the amount of output:
C                =0 - standard output:
C                     condensed output on unit 6 (basics + error messages),
C                     no output on unit 96 (depths of formation);
C                     normal output on 16 (equivalent widths);
C                     normal output on 12 (identification table)
C                >0 - more output:
C                     =1 - emergent flux on unit 6, no unit 96
C                     =2 - identification table + flux on unit 6, no unit 96
C                     =3 - as before, plus unit 96 (depths of formation);
C                     =4 - as before, plus unit 97 (contribution functions);
C                <0 - less output:
C                     =-1 - no output on unit 16
C                     =-2 - no output on units 16 and 12
C
C     INMOD      = 0  -  input model atmosphere as a Kurucz model
C                        (read by procedure INKUR)
C                = 1  -  input model atmosphere is a model calculated
C                        by the program TLUSTY
C                        (read by procedure INPMOD)
C                = 2  -  input model is a model of the vertical structure
C                        of one ring of an accretion disk
C     INTRPL     - switch indicating whether the input model has to be
C                  interpolated to the present depth scale;
C                  for details see procedure INPMOD
C     ICHANG     - switch indicating whether the populations from the
C                  input model have to be updated;
C                  for details see procedure CHANGE
C     ICHEMC     - switch indicating that new chemical composition will
C                  be read from unit 56
C     IOPHLI     - switch for treatment the Lyman line wings -see LYMLIN
C
   
      mode=0
      read(1,*,err=10,end=10) mode
   10 continue
c     IFMOL=0
      IFWIN=0
      nunalp=0
      nunbet=0
      nungam=0
      nunbal=0
      iunitm(1)=20
      nmlist=0
      NDSTEP=0
      if(ifeos.eq.0) then
      READ(55,*,END=3) IMODE,IDSTD,IPRIN
      READ(55,*,END=3) INMOD,INTRPL,ICHANG,ICHEMC
      if(mode.gt.0) inmod=2
      READ(55,*,ERR=3,END=3) IOPHLI,nunalp,nunbet,nungam,nunbal
      go to 4
    3 continue
      nunalp=0
      nunbet=0
      nungam=0
      nunbal=0
    4 continue
      end if
      IF(IMODE.LT.-90) THEN
         IMODE=-IMODE-100
         IFWIN=1
      END IF
      if(imode.gt.5) then
         imode=imode-10
         ifmol=1
         nmlist=1
         iunitm(1)=20
      end if
c     disableing an old option
      iophli=0
      call initia
c
c     if needed, read tables with data for quasimolecular satellites of
c     Lyman alpha, beta, gamma, and Balmer alpha
c
      call getlal
c
      IF(IMODE.LT.-1) THEN
         ND=1
         IDSTD=1
      END IF
      IF(INMOD.GT.0.AND.INTRPL.GT.0) READ(55,*) (DM(I),I=1,ND)
C
      return
      end
C
C
C     ****************************************************************
C
C

      SUBROUTINE INITIA
C     =================
C
C     driver for input and initializations - "new" routine
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (WI1=911.753578, WI2=227.837832)
      common/dissol/fropc(mlevel),indexp(mlevel)
      CHARACTER*10 TYPLEV(MLEVEL)
      CHARACTER*4 TYPION(MIOEX),TYPIOI
      CHARACTER*40 FIDATA(MION),FIODF1(MION),FIODF2(MION),FIBFCS(MION),
     *             FILEI
      CHARACTER*20 FINSTD
      CHARACTER*1 BLNK
      COMMON/PRINTP/TYPLEV
      COMMON/IONDAT/IATI(MION),IZI(MION),NLEVS(MION),NLLIM(MION)
      COMMON/IONFIL/FIDATA,FIODF1,FIODF2,FIBFCS
      COMMON/INUNIT/IUNIT
      COMMON/STRPAR/IMER,ITR,IC,IL,IP,NLASTE,NHOD
      common/quasex/iexpl(mlevel),iltot(mlevel)
      DIMENSION IGLE(18),IGMN(25),IGFE(26),IGNI(28)
      DATA IGLE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1/
      DATA IGMN/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,7,6/
      DATA IGFE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,30,25/
      DATA IGNI/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,28,21,10,21/
      DATA BLNK /' '/
C
      CALL READBF
C
C ------------------------------------
C Basic input parameters - atmospheres
C ------------------------------------
C
      IF(INMOD.LE.1) THEN
         READ(IBUFF,*) TEFF,GRAV
       ELSE IF(INMOD.EQ.2) THEN
C
C ------------------------------
C Basic input parameters - disks
C ------------------------------
C
         READ(IBUFF,*) DISPAR
      END IF
C
C ----------------------------
C other basic input parameters 
C ----------------------------
C
      READ(IBUFF,*) LTE,LTGREY
      READ(IBUFF,*) FINSTD
      CALL NSTPAR(FINSTD)
C
C
C ----------------------------
C Frequency points and weights
C ----------------------------
C
      READ(IBUFF,*) NFREAD
      NJREAD=NFREAD
C
      IF(NJREAD.LT.0) THEN
         NJREAD=-NJREAD
         NFREQC=NJREAD
         DO IJ=1,NJREAD
            READ(IBUFF,*) FREQEXP
         END DO   
       ELSE
         NFREQC=NJREAD
      END IF
C
      WRITE(6,601) TEFF,GRAV
C
C ----------------------------------------------------
C     turbulent velocities
C ----------------------------------------------------
C
      IF(VTB.LT.1.E3) VTB=VTB*1.E5
      DO ID=1,ND
         VTURB(ID)=VTB
      END DO
C
C ----------------------------------------------------
C Input parameters for explicit and non-explicit atoms
C ----------------------------------------------------
C
C     Input parameters are read by procedure STATE
C     (see description there)
C
      CALL STATE0(1)
      ID=1
      IF(IPRIN.GE.1) WRITE(6,607) YTOT(ID),WMY(ID),WMM(ID)
      DO I=1,MLEVEL
         ILK(I)=0
         iexpl(i)=0
         iltot(i)=0
      END DO
C
C --------------------------------------------------------------
C Input of parameters for explicit ions, levels, and transitions
C --------------------------------------------------------------
C
      ILEV=0
      IATLST=0
      ION=0
      IA=0
      IUNIT=34
      NATOM=0
      WRITE(6,613)
   10 CONTINUE
      READ(IBUFF,*,END=20,ERR=20) IATII,IZII,NLEVSI,ILASTI,ILVLIN,
     *              NONSTD,TYPIOI,FILEI
      IF(ILASTI.EQ.0) THEN
         ION=ION+1
         IATI(ION)=IATII
         IZI(ION)=IZII
         NLEVS(ION)=NLEVSI
         TYPION(ION)=TYPIOI
         FIDATA(ION)=FILEI
         NLLIM(ION)=ILVLIN
         ILIMITS(ION)=-1
         IUPSUM(ION)=0
         FIBFCS(ION)=BLNK
         MODEFF=1
         NFF=0
         IF(IATI(ION).EQ.1.AND.IZI(ION).EQ.0) THEN 
            IUPSUM(ION)=-100
            MODEFF=2
         END IF
         IF(IATI(ION).EQ.2.AND.IZI(ION).EQ.1) THEN
            MODEFF=2
         END IF
         IF(NONSTD.GE.10) THEN
            WRITE(*,*)'INITIA: QUANTUM NUMBERS AND ENERGY LIMITS WILL'
            WRITE(*,*)'        BE IGNORED FOR ION ',IATII,'    ',IZII
            ILIMITS(ION)=0
            NONSTD=NONSTD-10
         END IF
         IF(NONSTD.GT.0) THEN
            READ(IBUFF,*) IUPSUM(ION),ICUP,MODEFF,NFF
          ELSE IF(NONSTD.LT.0) THEN
            READ(IBUFF,*) ifil1,ifil2,FIODF1(ION),
     *                    FIODF2(ION),FIBFCS(ION)
            IF(FIBFCS(ION).NE.' ') THEN
               IUNIT=IUNIT+1
               INBFCS(ION)=IUNIT
            END IF
            IUPSUM(ION)=1
         END IF
C
         IF(IATI(ION).EQ.IATLST) THEN
            NFIRST(ION)=ILEV
          ELSE
            NFIRST(ION)=ILEV+1
            IATLST=IATI(ION)
            IA=IATEX(IATLST)
            N0A(IA)=NFIRST(ION)
            NATOM=MAX(NATOM,IA)
         END IF
         NLAST(ION)=NFIRST(ION)+NLEVS(ION)-1
         NNEXT(ION)=NLAST(ION)+1
         ILEV=NNEXT(ION)
         IZ(ION)=IZI(ION)+1 
         IF(NFF.GT.0) FF(ION)=EH/H*IZ(ION)*IZ(ION)/NFF/NFF
C
         N0I=NFIRST(ION)
         N1I=NLAST(ION)
         NKI=NNEXT(ION)
         IFREE(ION)=MODEFF
         DO II=N0I,N1I
            IEL(II)=ION
            IATM(II)=IA
         END DO
         ILK(NKI)=ION
         IATM(NKI)=IA
C
         IF(NUMAT(IA).EQ.1) THEN
            IATH=IA
            IF(IZ(ION).EQ.1) IELH=ION
            IF(IZ(ION).EQ.0) IELHM=ION
         END IF
         IF(NUMAT(IA).EQ.2) THEN
            IATHE=IA
            IF(IZ(ION).EQ.1) IELHE1=ION
            IF(IZ(ION).EQ.2) IELHE2=ION
         END IF
C         
         IF(IPRIN.GE.0) 
     *   WRITE(6,614) TYPION(ION),N0I,N1I,NKI,IZ(ION)
C   
       ELSE IF(ILASTI.GT.0) THEN
         ENION(ILEV)=0.
         G(ILEV)=ILASTI
         NQUANT(ILEV)=1
         TYPLEV(ILEV)=TYPIOI
         IFWOP(ILEV)=0
         IEL(ILEV)=ION
         NKA(IA)=NNEXT(ION)
         IF(ILASTI.EQ.1.AND.IATII.GT.IZII) THEN
            IF(IATII.LT.25) THEN
               G(ILEV)=IGLE(IATII-IZII)
             ELSE IF(IATII.EQ.25) THEN
               G(ILEV)=IGMN(IATII-IZII)
             ELSE IF(IATII.EQ.26) THEN
               G(ILEV)=IGFE(IATII-IZII)
             ELSE IF(IATII.EQ.28) THEN
               G(ILEV)=IGNI(IATII-IZII)
            ENDIF
         ENDIF
       ELSE
         GO TO 20
      END IF 
      GO TO 10    
   20 CONTINUE
      NION=ION
      NLEVEL=NKI
C     
      if(iath.gt.0) then
         N0H=N0A(IATH)
         N1H=NLAST(IELH)
         NKH=NNEXT(IELH)
         N0HN=NFIRST(IELH)
         N0M=0
         IF(IELHM.GT.0) THEN
            N0M=NFIRST(IELHM)
            IOPHMI=0
         end if
       else
         n0h=0
         n1h=0
         nkh=0
         n0hn=0
      end if
C
      IF(IPRIN.GE.1) WRITE(6,603) INMOD,ND,IDSTD,INTRPL,ICHANG,
     *             NATOM,NION,NLEVEL,
     *             IELH,IELHM,IATH
C
C -----------------------------------------
C Parameters for individual explicit levels
C -----------------------------------------
C
      IMER=0
      ITR=0
      IC=0
      IL=0
      IP=0
C
      DO ION=1,NION
         CALL RDATA(ION) 
         NFF=NQUANT(NLAST(ION))+1
         IF(NFF.GT.0) FF(ION)=EH/H*IZ(ION)*IZ(ION)/NFF/NFF
      END DO
C
      IF(IPRIN.GE.1) WRITE(6,615)
      DO I=1,NLEVEL
         IF(IPRIN.GE.1)
     *   WRITE(6,616) I,TYPLEV(I),TYPION(IEL(I)),ENION(I),G(I),
     *                NQUANT(I),IEL(I),ILK(I),IATM(I)
      END DO
C
C -----------------------------------------
C Input parameters for additional opacities
C -----------------------------------------
C
      IF(IPRIN.GE.0) WRITE(6,605) IOPHMI,IOPH2P,IOPHEM,IOPCH,IOPOH,
     *               IOPH2M,IOH2H2,IOH2HE,IOH2H1,IOHHE,
     *               IRSCT,IRSCH2,IRSCHE,IOPHLI
C
C
      IF(VTB.LT.1.E3) VTB=VTB*1.E5
      DO ID=1,ND
         VTURB(ID)=VTB
      END DO   
      WRITE(6,608) VTB*1.E-5
      DO I=1,ND
         VTURB(I)=VTURB(I)*VTURB(I)
      END DO   
C
  601 FORMAT(31X,'*******************************************'/
     * 31X,'I',41X,'I'/
     * 31X,'I   S Y N T H E T I C   S P E C T R U M   I'/
     * 31X,'I',41X,'I'/
     * 31X,'I',8X,'FOR MODEL ATMOSPHERE WITH',8X,'I'/
     * 31X,'I',41X,'I'/
     * 31X,'I',14X,'TEFF  =',F7.0,13X,'I'/
     * 31X,'I',14X,'LOG G =',F7.2,13X,'I'/
     * 31X,'I',41X,'I'/
     * 31X,'*******************************************')
c 602 FORMAT(//' EXPLICIT LEVELS AND OPACITY SOURCES'/
c    *         ' -----------------------------------'//
c    *'   I',2X,'IONIZ.ENERGY',5X,'G  NQUANT IEL ILK IAT IBF ',
c    *3X,'S0',5X,'ALF',4X,'BET',4X,'GAM'/)
  603 FORMAT(//' BASIC INPUT PARAMETERS'/
     *         ' ----------------------'/
     *            ' INMOD  =',I5/
     *            ' ND     =',I5/
     *            ' IDSTD  =',I5/
     *            ' INTRPL =',I5/
     *            ' ICHANG =',I5/
     *            ' NATOM  =',I5/
     *            ' NION   =',I5/
     *            ' NLEVEL =',I5/
     *            ' IELH   =',I5/
     *            ' IELHM  =',I5/
     *            ' IATH   =',I5)
  605 FORMAT(//' ADDITIONAL OPACITY SOURCES'/
     *            ' --------------------------'/
     * ' IOPHMI  (H-  OPACITY IN LTE)       =',I3/
     * ' IOPH2P  (H2+  OPACITY)             =',I3/
     * ' IOPHEM  (HE- B-F AND F-F)          =',I3/
     * ' IOPCH   (CH OPACITY)               =',I3/
     * ' IOPOH   (OH OPACITY)               =',I3/
     * ' IOPH2M  (H2- OPACITY)              =',I3/
     * ' IOH2H2  (CIA H2-H2 OPACITY         =',I3/
     * ' IOH2HE  (CIA H2-He OPACITY         =',I3/
     * ' IOH2H1  (CIA H2-H  OPACITY         =',I3/
     * ' IOHHE   (CIA H-He OPACITY          =',I3/
     * ' IRSCT   (RAYLEIGH SCAT. ON H I)    =',I3/
     * ' IRSCH2  (RAYLEIGH SCAT. ON H2      =',I3/
     * ' IRSCHE  (RAYLEIGH SCAT. ON HE I)   =',I3/
     * ' IOPHLI  (LYMAN LINES WINGS)        =',I3)
  607 FORMAT(///' YTOT   =',F11.5/' WMY    =',1PE15.5/
     *  ' WMM    =',E15.5)
  608 FORMAT(//' TURBULENT VELOCITY  -  DEPTH-INDEPENDENT   VTURB =',
     * 1PE10.3,'  KM/S'/
     *            ' ------------------'/)
c 609 FORMAT(//' TURBULENT VELOCITY  -  DEPTH-DEPENDENT'/
c    *            ' ------------------'/
c    * 1H ,1P10E10.3)
c 611 FORMAT(1H ,I3,1PE15.7,0PF9.2,I3,4I4,1X,4F7.2)
  613 FORMAT(//' EXPLICIT IONS INCLUDED'/
     *            ' ----------------------'//
     *  ' ION     N0    N1    NK    IZ'/)
  614 FORMAT(1H ,A4,4I6)
  615 FORMAT(//' EXPLICIT ENERGY LEVELS INCLUDED'/
     *            ' -------------------------------'//
     * ' NO.   LEVEL    ION   ION.EN.(ERG)        G   NQUANT',
     * '  IEL  ILK  IAT'/)
  616 FORMAT(1H ,I3,2X,A10,A4,1PE15.7,0PF10.2,4I5)
c 618 FORMAT(F11.3,F10.3)
C
      RETURN
      END
C     
C     
C ************************************************************************    
C     
C     
C  
      SUBROUTINE RDATA(ION)
C     =====================
C 
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (WI1=911.753578, WI2=227.837832)
      PARAMETER (T15=1.D-15)
      PARAMETER (ECONST=  5.03411142E15)
      CHARACTER*10 TYPLEV(MLEVEL)
      CHARACTER*40 FIDATA(MION),FIODF1(MION),FIODF2(MION),FIBFCS(MION)
      CHARACTER*1 A
      CHARACTER*1000 CADENA
      COMMON/IONDAT/IATI(MION),IZI(MION),NLEVS(MION),NLLIM(MION)
      COMMON/IONFIL/FIDATA,FIODF1,FIODF2,FIBFCS
      COMMON/TOPCS/CTOP(MFIT,MCROSS), !sigma=alog10(sigma/10^-18) of fit point
     *             XTOP(MFIT,MCROSS)  ! x = alog10(nu/nu0) of fit point 
      COMMON/PRINTP/TYPLEV
      COMMON/INUNIT/IUNIT
      COMMON/STRPAR/IMER,ITR,IC,IL,IP,NLASTE,NHOD
      common/dissol/fropc(mlevel),indexp(mlevel)
      common/quasex/iexpl(mlevel),iltot(mlevel)
      data iexp0/0/
C
c     IUNIT=IUNIT+1
      IUNIT=94
      OPEN(IUNIT,FILE=FIDATA(ION),STATUS='OLD')
C
C     read the first record - a label for the energy level input
C  
      READ(IUNIT,501) A
  501 FORMAT(A1)
C
C   -----------------------------------------------------
C   input parameters for explicit energy levels
C   -----------------------------------------------------
C  
C   If ILIMITS(ION) < 0, the program finds out whether energy and
C   quantum numbers are included in the input data files

      IF (ILIMITS(ION).LT.0) THEN
        READ(IUNIT,'(1000A)')CADENA
        BACKSPACE(IUNIT)
        CALL COUNT_WORDS(CADENA,NOW)
        IF (NOW.LT.14) THEN
          ILIMITS(ION)=0
        ELSE
          ILIMITS(ION)=1
        ENDIF
      ENDIF

C   Standard format: ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i)

      IF (ILIMITS(ION).EQ.0) THEN
C
      DO IL=1,NLEVS(ION)
         I=IL+NFIRST(ION)-1
         IE=IEL(I)
         N0I=NFIRST(IE)
         NKI=NNEXT(IE)
         ia=numat(iatm(n0i))
         if(isemex(ia).le.1) then
            iexp0=iexp0+1
            iexpl(i)=iexp0
            iltot(iexp0)=i
c      write(6,671) il,i,ia,ion,isemex(ia),iexp0,iltot(iexp0)
            if(il.eq.nlevs(ion)) then
               if(nki.eq.nka(iatm(i))) then
                  iexp0=iexp0+1
                  iexpl(nki)=iexp0
                  iltot(iexp0)=nki
c      write(6,671) il+1,nki,ia,ion,isemex(ia),iexp0,iltot(iexp0)
               end if
            end if
c 671 format('il,i,ia,ion,isem,iexp,iltot',7i4)
         end if
         IQ=I-N0I+1
         X=IQ*IQ
         ifwop(i)=0
         IZZ=IZ(IE)
         READ(IUNIT,*) 
     *   ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i)
         if(ifwop(i).lt.0.and.i.ne.nlast(ie)) 
     *   call quit('conflict in negative ifwop')
         if(ifwop(i).ge.2) ifwop(i)=0 
         IF(I.LT.NKI) THEN
            E=ENION(I)
            E0=E
            IF(E.LT.0.) THEN
               E=-E
               E0=E
            END IF
            IF(E.EQ.0.) THEN
c              if(izz.le.2) then
               if(izz.le.-2) then
               w0=wi1
               if(izz.eq.2) w0=wi2
               WL0=W0*X
               IF(WL0.GT.2000.) THEN
                  ALM=1.E8/(WL0*WL0) 
                  XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
                  WL0=WL0/(XN1*1.D-6+1.D0)
               END IF 
               E0=H*CL*1.D8/WL0       
               else
               E0=EH*IZZ*IZZ/X
               end if
            END IF
            IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018D-12*E
            IF(E.GT.100..AND.E.LT.1.D7) E0=1.9857D-16*E
            IF(E.GT.1.D7) E0=H*E
            IF(ENION(I).GE.0.) THEN
              ENION(I)=E0
            ELSE
              ENION(I)=-E0
            ENDIF
            IF(G(I).EQ.0.) G(I)=2.D0*X
            IF(NQUANT(I).EQ.0) NQUANT(I)=IQ
          ELSE
c            if(modref.ge.0) nref(iatm(i))=nka(iatm(i))
            IF(G(I).EQ.0..AND.NKI.EQ.NKA(IATM(I))) G(I)=1.
         END IF
         if(ifwop(i).lt.0) then
            enion(i)=0.
            ff(ie)=0.
            IMER=IMER+1
            IMRG(I)=IMER
            IIMER(IMER)=I
         endif
         fropc(i)=0.
      END DO

C     Upgraded format including limits for energies, and quantum numbers

      ELSE

      DO IL=1,NLEVS(ION)
         I=IL+NFIRST(ION)-1
         IE=IEL(I)
         N0I=NFIRST(IE)
         NKI=NNEXT(IE)
         ia=numat(iatm(n0i))
         if(isemex(ia).le.1) then
            iexp0=iexp0+1
            iexpl(i)=iexp0
            iltot(iexp0)=i
            if(il.eq.nlevs(ion)) then
               if(nki.eq.nka(iatm(i))) then
                  iexp0=iexp0+1
                  iexpl(nki)=iexp0
                  iltot(iexp0)=nki
               end if
            end if
         end if
         IQ=I-N0I+1
         X=IQ*IQ
         ifwop(i)=0
         IZZ=IZ(IE)
         READ(IUNIT,*)
     *   ENION(I),G(I),NQUANT(I),TYPLEV(I),ifwop(i),frdodf,imodl,
     *   ENION1(I),ENION2(I),
     *   SQUANT1(I),SQUANT2(I),
     *   LQUANT1(I),LQUANT2(I),
     *   PQUANT1(I),PQUANT2(I)
         if(ifwop(i).lt.0.and.i.ne.nlast(ie))
     *   call quit('conflict in negative ifwop')
         if(ifwop(i).ge.2) ifwop(i)=0
         IF(I.LT.NKI) THEN

C           check and, if necessary, transform ENION(I)

            E=ENION(I)
            E0=E
            IF(E.LT.0.) THEN
               E=-E
               E0=E
            END IF
            IF(E.EQ.0.) THEN
c              if(izz.le.2) then
               if(izz.le.-2) then
               w0=wi1
               if(izz.eq.2) w0=wi2
               WL0=W0*X
               IF(WL0.GT.2000.) THEN
                  ALM=1.E8/(WL0*WL0)
                  XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
                  WL0=WL0/(XN1*1.D-6+1.D0)
               END IF
               E0=H*CL*1.D8/WL0
               else
               E0=EH*IZZ*IZZ/X
               end if
            END IF
            IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018D-12*E
            IF(E.GT.100..AND.E.LT.1.D7) E0=1.9857D-16*E
            IF(E.GT.1.D7) E0=H*E
            IF(ENION(I).GE.0.) THEN
              ENION(I)=E0
            ELSE
              ENION(I)=-E0
            ENDIF

C           check and, if necessary, transform ENION1(I)

            E=ENION1(I)
            E0=E
            IF(E.LT.0.) THEN
               E=-E
               E0=E
            END IF
            IF(E.EQ.0.) THEN
c              if(izz.le.2) then
               if(izz.le.-2) then
               w0=wi1
               if(izz.eq.2) w0=wi2
               WL0=W0*X
               IF(WL0.GT.2000.) THEN
                  ALM=1.E8/(WL0*WL0)
                  XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
                  WL0=WL0/(XN1*1.D-6+1.D0)
               END IF
               E0=H*CL*1.D8/WL0
               else
               E0=EH*IZZ*IZZ/X
               end if
            END IF
            IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018D-12*E
            IF(E.GT.100..AND.E.LT.1.D7) E0=1.9857D-16*E
            IF(E.GT.1.D7) E0=H*E
            IF(ENION1(I).GE.0.) THEN
              ENION1(I)=E0
            ELSE
              ENION1(I)=-E0
            ENDIF

C           check and, if necessary, transform ENION2(I)

            E=ENION2(I)
            E0=E
            IF(E.LT.0.) THEN
               E=-E
               E0=E
            END IF
            IF(E.EQ.0.) THEN
c              if(izz.le.2) then
               if(izz.le.-2) then
               w0=wi1
               if(izz.eq.2) w0=wi2
               WL0=W0*X
               IF(WL0.GT.2000.) THEN
                  ALM=1.E8/(WL0*WL0)
                  XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
                  WL0=WL0/(XN1*1.D-6+1.D0)
               END IF
               E0=H*CL*1.D8/WL0
               else
               E0=EH*IZZ*IZZ/X
               end if
            END IF
            IF(E.GT.1.D-7.AND.E.LT.100.) E0=1.6018D-12*E
            IF(E.GT.100..AND.E.LT.1.D7) E0=1.9857D-16*E
            IF(E.GT.1.D7) E0=H*E
            IF(ENION2(I).GE.0.) THEN
              ENION2(I)=E0
            ELSE
              ENION2(I)=-E0
            ENDIF

C
C       Enforce an energy tolerance of 10% when the input files
C       do not have  any (e.g. pure levels in MODION models)
C
            IF((ENION1(I)-ENION(I))/ENION(I).LT.1e-6)
     *      ENION1(I)=ENION(I)*(1.+ERANGE)
            IF((ENION(I)-ENION2(I))/ENION(I).LT.1e-6)
     *      ENION2(I)=ENION(I)*(1.-ERANGE)

C
C       Convert ENION1,ENION2 to cm-1 from the ground level
C       so they can be directly used in NLTSET
C

            ENION1(I)=(ENION(N0I)-ENION1(I))*ECONST
            ENION2(I)=(ENION(N0I)-ENION2(I))*ECONST


            IF(G(I).EQ.0.) G(I)=2.D0*X
            IF(NQUANT(I).EQ.0) NQUANT(I)=IQ
          ELSE
c            if(modref.ge.0) nref(iatm(i))=nka(iatm(i))
            IF(G(I).EQ.0..AND.NKI.EQ.NKA(IATM(I))) G(I)=1.
         END IF
         if(ifwop(i).lt.0) then
            write(*,*)'RDATA:  IFWOP<0 and ILIMITS is not 0'
            stop
            enion(i)=0.
            ff(ie)=0.
            IMER=IMER+1
            IMRG(I)=IMER
            IIMER(IMER)=I
         endif
         fropc(i)=0.
      END DO

      END IF

c
C ----------------------------------------------------------------------
C
C   skip lines if more levels than needed, and skip the continuum transition 
C   label
C
    5 READ(IUNIT,501) A
      IF(A.NE.'*') GO TO 5
      II0=NFIRST(ION)-1
      ILLIM=NLLIM(ION)+II0
      JCORR=0
C
C   -----------------------------------------------------
C   input parameters for continuum transitions
C   -----------------------------------------------------
C
   10 CONTINUE
      READ(IUNIT,*,END=20,ERR=15) II,JJ,MODE,IFANCY,ICOLIS,
     *                            IFRQ0,IFRQ1,OSC,CPARAM
      IF(II.EQ.0) THEN
         IF(JJ.EQ.0) GO TO 30
         II0=JJ-1
         GO TO 10
      END IF
      IF(IABS(MODE).GT.100) READ(IUNIT,*) FR0INP
      if(iabs(mode).eq.2) then
         READ(IUNIT,*) kdo
         go to 10
      end if
      IF(IFANCY.GT.49.and.ifancy.lt.100) IASV=1
      if(iabs(mode).eq.3.or.iabs(mode).eq.4) go to 10
      IF(IABS(MODE).EQ.5 .OR. IABS(MODE).EQ.15) THEN
         READ(IUNIT,*) FROPCI
         if(ion.eq.ielh) then
            if(ii.eq.1.and.cutlym.ne.0) fropci=-cutlym
            if(ii.eq.2.and.cutbal.ne.0) fropci=-cutbal
         end if
         if(abs(fropci).lt.1.e10) fropci=2.997925e18/fropci
      END IF
      IF(II.EQ.1) JCORR=NLEVS(ION)+1-JJ
      II=II+II0
      JJ=JJ+II0+JCORR
      FROPC(II)=FROPCI
      N0I=NFIRST(IE)
      NKI=NNEXT(IE)
      IF(JJ.GE.NKI) THEN
      LPC=.FALSE.
      IF(IELHE2.GE.0) THEN
         IF(II.GE.NFIRST(IELHE2).AND.II.LE.NLAST(IELHE2)
     *     .AND.IFWOP(II).GE.0) LPC=.TRUE.
      END IF
      IF(II.GE.N0HN.AND.II.LE.N1H.AND.IFWOP(II).GE.0) LPC=.TRUE.
      IF(LPC) THEN
         MODE=5
         XI=NQUANT(II)
         X2=XI+3.
         if(ii.ge.8) x2=xi+2.
         IF(FROPC(II).GE.0.) THEN
            FROPC(II)=ENION(II)/6.6256E-27*(1.-XI*XI/(X2*X2))
          ELSE
            FROPC(II)=ABS(FROPC(II))
         END IF
c     write(6,671) ii,fropc(ii),enion(ii)/h,2.997925e18/fropc(ii)
c 671 format(i4,1p2e13.5,0pf10.1)
      END IF
      END IF
      IF(MODE.EQ.0) THEN
         IF(II.LT.NLAST(ION)) GO TO 10
         IF(II.EQ.NLAST(ION)) GO TO 15
      END IF
C
C   -----------------------------------------------------
C   Additional input parameters for continuum transitions
C   -----------------------------------------------------
C
C     Only for IFANCY = 2, 3, or 4
C     S0BF, ALFBF, BETBF, GAMBF  - parameters for evaluation the
C         photoionization cross-section
C
      IF(IFANCY.GE.2.AND.IFANCY.LE.4)
     *   READ(IUNIT,*) S0BF(II),ALFBF(II),BETBF(II),GAMBF(II)
C
C   -----------------------------------------------------
C   Additional input parameters for continuum transitions -TOPBASE DATA
C   -----------------------------------------------------
C
C     Only for IFANCY > 100 there are IFANCY-100 fit points
C
C     XTOP(MFIT,MCROSS) -  x = alog10(nu/nu0) of a fit point
C     CTOP(MFIT,MCROSS) -  sigma = alog10(sigma/10^-18) of a fit point
C
C     there are IFANCY-100 fit points
C
      IF(IFANCY.GT.100) THEN
         NFIT=IFANCY-100
         IF(NFIT.GT.MFIT) call quit(' nfit too large (TOPBASE fits)')
         READ(IUNIT,*) (XTOP(IFIT,II),IFIT=1,NFIT)
         READ(IUNIT,*) (CTOP(IFIT,II),IFIT=1,NFIT)
      END IF
      IBF(II)=IFANCY
      INDEXP(II)=IABS(MODE)
      IF(II.LT.NLAST(ION)) GO TO 10
   15 READ(IUNIT,501) A
      IF(A.NE.'*') GO TO 15
C
C  -----------------------------------------------------------
C  Input parameters for line transitions
C  -----------------------------------------------------------
C
   20 CONTINUE
      READ(IUNIT,*,END=30,ERR=30) II,JJ,MODE,IFANCY,ICOLIS,
     *                              IFRQ0,IFRQ1,OSC,CPARAM
      IF(IABS(MODE).GT.100) READ(IUNIT,*) FR0INP
      IF(JJ.GT.NLEVS(ION)) THEN
         IF(IABS(MODE).EQ.2) THEN
            READ(IUNIT,*) K1,K2,K3,X1,X2,X3,K4
            GO TO 20
         END IF
         IF(IABS(MODE).EQ.1) READ(IUNIT,*) LCMP
         IF(IABS(IFANCY).EQ.1) READ(IUNIT,*) GAMR,STARK1,STARK2,
     *              STARK3,VDWH
         GO TO 20
      END IF
      if(iabs(mode).eq.2) then
         READ(IUNIT,*) K1,K2,K3,X1,X2,X3,K4
         go to 20
      end if
      if(iabs(mode).eq.3.or.iabs(mode).eq.4) go to 20
      IF(MODE.EQ.0) GO TO 20
C
C  -----------------------------------------------------------
C  Additional input parameters for "clasical" line transitions
C   (i.e. those not represented by ODF's - ie ABS(MODE)=1)
C  -----------------------------------------------------------
C
      READ(IUNIT,*) LCOMP,INTMOD,NF,XMAX,TSTD
      IF(IABS(IFANCY).EQ.1) READ(IUNIT,*) GAMR,STARK1,STARK2,
     *              STARK3,VDWH
      GO TO 20
C
   30 CONTINUE
      close(iunit)
      RETURN
      END
c
c
C    *****************************************************************
c
C
      SUBROUTINE NSTPAR(FINSTD)
C     ==========================
C
C     settiing up the default values of various input flags, and
C     input of non-standard values of various input flags and parameters
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      common/hhebrd/sthe,nunhhe
      common/hydmol/anhmi,ahmol,ih2,ih2p,ihm
      common/gompar/hglim,ihgom
C 
      PARAMETER(MVAR=60)
      PARAMETER(INPFI=4)
      CHARACTER*(*) FINSTD
      CHARACTER*6 PVALUE(MVAR)
      CHARACTER*80 TEXT
      CHARACTER VARNAM(MVAR)*6
      CHARACTER*20 BLNK
      CHARACTER*6 BLNK6
C
      DATA VARNAM /'ISPLIN','IRTE  ','IBC   ','ILMCOR','ILPSCT',
     *             'IFPREC','IELCOR','ICHC  ','IRSPLT','IATREF',
     *             'POPZER','BERGFC','IHYDPR','NUNHHE','STHE  ',
     *             'IOVER ','ITLAS ','IFSUB ','NITER ','NLAMBD',
     *             'ND    ','NFREQS','IBFAC ',
     *             'INTRPL','ICHANG','IFEOS',
     *             'IOPHMI','IOPH2P','IOPHEM','IOPCH ','IOPOH ',
     *             'IOPH2M','IOH2H2','IOH2HE','IOH2H1','IOHHE ',
     *             'IRSCT ','IRSCH2','IRSCHE',
     *             'IHM   ','IH2   ','IH2P  ',
     *             'TRAD  ','WDIL  ',
     *             'HMIX0 ','VTB   ','IFMOL','TMOLIM',
     *             'MOLTAB','IRWTAB','IIRWIN',
     *             'CUTLYM','CUTBAL','IHXENB',
     *             'IHGOM ','HGLIM ',
     *             'ERANGE',
     *             'ISPICK','ILPICK','IPPICK'/

C
      DATA PVALUE /'     0','     3','     3','     1','     0',
     *             '     1','   100','     0','     1','     1',
     *             '1.D-20','  1.D0','     0','     0',' 1.e19',
     *             '     1','   100','     0','    30','     1',
     *             '    70','   120','     0',
     *             '     0','     0','     0',
     *             '     1','     1','     1','     1','     1',
     *             '     1','     1','     1','     1','     1',
     *             '     1','     1','     1',
     *             '     1','     1','     1',
     *             '    0.','    0.',
     *             '    0.','    2.','     1',' 8000.',
     *             '     0','     0','     0',
     *             '    0.','    0.','     0',
     *             '     0',' 1.e18',
     *             '  0.10',
     *             '     1','     1','     1'/
C
      DATA BLNK/'                    '/,BLNK6/'      '/
C
      IF(FINSTD.NE.BLNK) 
     *   OPEN(UNIT=INPFI,FILE=FINSTD,STATUS='UNKNOWN')
C
      INDV=-1
C
C    go through the input file line by line
c
      write(6,601)
  601 format(/' INPUT KEYWORD PARAMETERS:'/
     *        ' -------------------------')
c
   10 CONTINUE
      K0=1
      READ(INPFI,500,END=70,ERR=70) TEXT
  500 FORMAT(A)
      WRITE(6,*) TEXT
   20 CONTINUE
      CALL GETWRD(TEXT,K0,K1,K2)
      IF(K1.EQ.0) GO TO 60
      K0=K2+2
      IF(TEXT(K1:K2).EQ.'=') GO TO 20
      INDV=-INDV
      IF(INDV.EQ.1) THEN
         DO 40 I=1,MVAR
            IF(TEXT(K1:K2).EQ.VARNAM(I)(1:K2-K1+1)) GO TO 50
   40    CONTINUE
         CALL GETWRD(TEXT,K0,K1,K2)
         IF(K1.EQ.0) THEN
            K0=1
   45       READ(INPFI,500,END=70) TEXT
            CALL GETWRD(TEXT,K0,K1,K2)
            IF(K1.EQ.0) GO TO 45
         END IF
         K0=K2+2
         INDV=-INDV
         GO TO 20
   50    CONTINUE
         IVAR=I
       ELSE
         PVALUE(IVAR)=BLNK6
         PVALUE(IVAR)(6-K2+K1:6)=TEXT(K1:K2)
      END IF
      GO TO 20
   60 CONTINUE
      GO TO 10
   70 CONTINUE
C
      DO I=1,MVAR
         WRITE(84,684) PVALUE(I)
  684    FORMAT(1X,A)
      END DO
C
      CLOSE(UNIT=84)
      REWIND(84)
      READ(84,*) 
     *             ISPLIN,IRTE  ,IBC   ,ILMCOR,ILPSCT,
     *             IFPREC,IELCOR,ICHC  ,IRSPLT,IATREF,
     *             POPZER,BERGFC,IHYDPR,NUNHHE,STHE  ,
     *             IOVER ,ITLAS ,IFSUB ,NITER ,NLAMBD,
     *             ND    ,NFREQS,IBFAC ,
     *             INTRPL,ICHANG,IFEOS ,
     *             IOPHMI,IOPH2P,IOPHEM,IOPCH ,IOPOH ,
     *             IOPH2M,IOH2H2,IOH2HE,IOH2H1,IOHHE ,
     *             IRSCT ,IRSCH2,IRSCHE,
     *             IHM   ,IH2   ,IH2P  ,
     *             TRAD  ,WDIL  ,
     *             HMIX0 ,VTB   ,IFMOL ,TMOLIM,
     *             MOLTAB,IRWTAB,IIRWIN,
     *             CUTLYM,CUTBAL,IHXENB,
     *             IHGOM ,HGLIM ,
     *             ERANGE,
     *             ISPICK,ILPICK,IPPICK
C      
      LCHC=.FALSE.
      IF(ICHC.EQ.1) LCHC=.TRUE.
C
      if(imode.le.-3) then
         irsct=0
         irsche=0
         irsch2=0
      end if
C
      RETURN
      END
C
C
C    ***************************************************************
C
C
        subroutine count_words(cadena,n)
C
C       Counts the number of words separated by blanks in a string
C
        character*1000  cadena
        character*1     a,b

        n=0
        a=cadena(1:1)
        if (a.ne.' ') n=1
        do i=2,len(cadena)
         b=cadena(i:i)
         if(b.ne.' '.and.a.eq.' ') n=n+1
         a=b
        enddo
        end
C
C
C    ***************************************************************
C
C
      SUBROUTINE GETWRD(TEXT,K0,K1,K2)
C
C  FINDS NEXT WORD IN TEXT FROM INDEX K0. NEXT WORD IS TEXT(K1:K2)
C  THE NEXT WORD STARTS AT THE FIRST ALPHANUMERIC CHARACTER AT K0
C  OR AFTER. IT ENDS WITH THE LAST ALPHANUMERIC CHARACTER IN A ROW
C  FROM THE START
C
C  TAKEN FROM MULTI - M. CARLSSON (1976)
C
C     INCLUDE 'IMPLIC.FOR'
      PARAMETER (MSEPAR=7)
      CHARACTER*(*) TEXT
      CHARACTER SEPAR(MSEPAR)
      DATA SEPAR/' ','(',')','=','*','/',','/
C
      K1=0
      DO 400 I=K0,LEN(TEXT)
        IF(K1.EQ.0) THEN
          DO 100 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 200
  100     CONTINUE
          K1=I
C
C  NOT START OF WORD
C
  200     CONTINUE
        ELSE
          DO 300 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 500
  300     CONTINUE
        ENDIF
  400 CONTINUE
C
C  NO NEW WORD. RETURN K1=K2=0
C
      K1=0
      K2=0
      GOTO 999
C
C  NEW WORD IN TEXT(K1:I-1)
C
  500 CONTINUE
      K2=I-1
C
  999 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE STATE0(MODOLD)
C     =========================
C
C     Initialization of the basic parameters for the Saha equation
C
      INCLUDE 'PARAMS.FOR'
      parameter (enhe1=24.5799,enhe2=54.3999)
      character*4 DYP
      DIMENSION D(3,MATOM),XI(8,MATOM),DYP(MATOM)
C
      DATA DYP/' H  ',' He ',' Li ',' Be ',' B  ',' C  ',
     *         ' N  ',' O  ',' F  ',' Ne ',' Na ',' Mg ',
     *         ' Al ',' Si ',' P  ',' S  ',' Cl ',' Ar ',
     *         ' K  ',' Ca ',' Sc ',' Ti ',' V  ',' Cr ',
     *         ' Mn ',' Fe ',' Co ',' Ni ',' Cu ',' Zn ',
     *         ' Ga ',' Ge ',' As ',' Se ',' Br ',' Kr ',
     *         ' Rb ',' Sr ',' Y  ',' Zr ',' Nb ',' Mo ',
     *         ' Tc ',' Ru ',' Rh ',' Pd ',' Ag ',' Cd ', 
     *         ' In ',' Sn ',' Sb ',' Te ',' I  ',' Xe ',
     *         ' Cs ',' Ba ',' La ',' Ce ',' Pr ',' Nd ', 
     *         ' Pm ',' Sm ',' Eu ',' Gd ',' Tb ',' Dy ', 
     *         ' Ho ',' Er ',' Tm ',' Yb ',' Lu ',' Hf ', 
     *         ' Ta ',' W  ',' Re ',' Os ',' Ir ',' Pt ', 
     *         ' Au ',' Hg ',' Tl ',' Pb ',' Bi ',' Po ', 
     *         ' At ',' Rn ',' Fr ',' Ra ',' Ac ',' Th ', 
     *         ' Pa ',' U  ',' Np ',' Pu ',' Am ',' Cm ', 
     *         ' Bk ',' Cf ',' Es '/
C
C    Standard atomic constants for first 99 species 
C      Abundances for the first 30 from Grevesse & Sauval,
C         (1998, Space Sci. Rev. 85, 161)
C
C            Element Atomic  Solar    Std.
C                    weight abundance highest 
C 
C                                     ionization stage 
      DATA D/ 1.008, 1.0D0, 2.,
     *        4.003, 1.00D-1, 3.,
     *        6.941, 1.26D-11, 3.,
     *        9.012, 2.51D-11, 3.,
     *       10.810, 5.0D-10, 4.,
     *       12.011, 3.31D-4, 5.,
     *       14.007, 8.32D-5, 5.,
     *       16.000, 6.76D-4, 5.,
     *       18.918, 3.16D-8, 4.,
     *       20.179, 1.20D-4, 4.,
     *       22.990, 2.14D-6, 4.,
     *       24.305, 3.80D-5, 4.,
     *       26.982, 2.95D-6, 4.,
     *       28.086, 3.55D-5, 5.,
     *       30.974, 2.82D-7, 5.,
     *       32.060, 2.14D-5, 5.,
     *       35.453, 3.16D-7, 5.,
     *       39.948, 2.52D-6, 5.,
     *       39.098, 1.32D-7, 5.,
     *       40.080, 2.29D-6, 5.,
     *       44.956, 1.48D-9, 5.,
     *       47.900, 1.05D-7, 5.,
     *       50.941, 1.00D-8, 5.,
     *       51.996, 4.68D-7, 5.,
     *       54.938, 2.45D-7, 5.,
     *       55.847, 3.16D-5, 5.,
     *       58.933, 8.32D-8, 5.,
     *       58.700, 1.78D-6, 5.,
     *       63.546, 1.62D-8, 5.,
     *       65.380, 3.98D-8, 5.,
     *       69.72 ,   1.34896324e-09  ,  3.,  
     *       72.60 ,   4.26579633e-09  ,  3.,  
     *       74.92 ,   2.34422821e-10  ,  3.,  
     *       78.96 ,   2.23872066e-09  ,  3.,  
     *       79.91 ,   4.26579633e-10  ,  3.,  
     *       83.80 ,   1.69824373e-09  ,  3.,  
     *       85.48 ,   2.51188699e-10  ,  3.,  
     *       87.63 ,   8.51138173e-10  ,  3.,  
     *       88.91 ,   1.65958702e-10  ,  3.,  
     *       91.22 ,   4.07380181e-10  ,  3.,  
     *       92.91 ,   2.51188630e-11  ,  3.,   
     *       95.95 ,   9.12010923e-11  ,  3.,   
     *       99.00 ,   1.00000000e-24  ,  3.,   
     *       101.1 ,   6.60693531e-11  ,  3.,   
     *       102.9 ,   1.23026887e-11  ,  3.,   
     *       106.4 ,   5.01187291e-11  ,  3.,   
     *       107.9 ,   1.73780087e-11  ,  3.,   
     *       112.4 ,   5.75439927e-11  ,  3.,   
     *       114.8 ,   6.60693440e-12  ,  3.,   
     *       118.7 ,   1.38038460e-10  ,  3.,   
     *       121.8 ,   1.09647810e-11  ,  3.,   
     *       127.6 ,   1.73780087e-10  ,  3.,   
     *       126.9 ,   3.23593651e-11  ,  3.,   
     *       131.3 ,   1.69824373e-10  ,  3.,   
     *       132.9 ,   1.31825676e-11  ,  3.,   
     *       137.4 ,   1.62181025e-10  ,  3.,   
     *       138.9 ,   1.58489337e-11  ,  3.,   
     *       140.1 ,   4.07380293e-11  ,  3.,   
     *       140.9 ,   6.02559549e-12  ,  3.,   
     *       144.3 ,   2.95120943e-11  ,  3.,   
     *       147.0 ,   1.00000000e-24  ,  3.,   
     *       150.4 ,   9.33254366e-12  ,  3.,   
     *       152.0 ,   3.46736869e-12  ,  3.,   
     *       157.3 ,   1.17489770e-11  ,  3.,   
     *       158.9 ,   2.13796216e-12  ,  3.,   
     *       162.5 ,   1.41253747e-11  ,  3.,   
     *       164.9 ,   3.16227767e-12  ,  3.,   
     *       167.3 ,   8.91250917e-12  ,  3.,   
     *       168.9 ,   1.34896287e-12  ,  3.,   
     *       173.0 ,   8.91250917e-12  ,  3.,   
     *       175.0 ,   1.31825674e-12  ,  3.,   
     *       178.5 ,   5.37031822e-12  ,  3.,   
     *       181.0 ,   1.34896287e-12  ,  3.,   
     *       183.9 ,   4.78630102e-12  ,  3.,   
     *       186.3 ,   1.86208719e-12  ,  3.,   
     *       190.2 ,   2.39883290e-11  ,  3.,   
     *       192.2 ,   2.34422885e-11  ,  3.,   
     *       195.1 ,   4.78630036e-11  ,  3.,   
     *       197.0 ,   6.76082952e-12  ,  3.,   
     *       200.6 ,   1.23026887e-11  ,  3.,   
     *       204.4 ,   6.60693440e-12  ,  3.,   
     *       207.2 ,   1.12201834e-10  ,  3.,   
     *       209.0 ,   5.12861361e-12  ,  3.,   
     *       210.0 ,   1.00000000e-24  ,  3.,   
     *       211.0 ,   1.00000000e-24  ,  3.,   
     *       222.0 ,   1.00000000e-24  ,  3.,   
     *       223.0 ,   1.00000000e-24  ,  3.,   
     *       226.1 ,   1.00000000e-24  ,  3.,   
     *       227.1 ,   1.00000000e-24  ,  3.,   
     *       232.0 ,   1.20226443e-12  ,  3.,   
     *       231.0 ,   1.00000000e-24  ,  3.,  
     *       238.0 ,   3.23593651e-13  ,  3.,  
     *       237.0 ,   1.00000000e-24  ,  3.,  
     *       244.0 ,   1.00000000e-24  ,  3.,  
     *       243.0 ,   1.00000000e-24  ,  3.,  
     *       247.0 ,   1.00000000e-24  ,  3.,  
     *       247.0 ,   1.00000000e-24  ,  3.,  
     *       251.0 ,   1.00000000e-24  ,  3.,  
     *       254.0 ,   1.00000000e-24  ,  3./
C
C
C     Ionization potentials for first 99 species:
      DATA XI/
C
C     Element Ionization potentials (eV) 
C              I     II      III     IV       V     VI     VII    VIII
C
     *       13.595,  0.   ,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *       24.580, 54.400,  0.   ,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *        5.392, 75.619,122.451,  0.   ,  0.  ,  0.  ,  0.  ,  0.  ,
     *        9.322, 18.206,153.850,217.713,  0.  ,  0.  ,  0.  ,  0.  ,
     *        8.296, 25.149, 37.920,259.298,340.22,  0.  ,  0.  ,  0.  ,
     *       11.264, 24.376, 47.864, 64.476,391.99,489.98,  0.  ,  0.  ,
     *       14.530, 29.593, 47.426, 77.450, 97.86,551.93,667.03,  0.  ,
     *       13.614, 35.108, 54.886, 77.394,113.87,138.08,739.11,871.39,
     *       17.418, 34.980, 62.646, 87.140,114.21,157.12,185.14,953.6 ,
     *       21.559, 41.070, 63.500, 97.020,126.30,157.91,207.21,239.0 ,
     *        5.138, 47.290, 71.650, 98.880,138.37,172.09,208.44,264.16,
     *        7.664, 15.030, 80.120,102.290,141.23,186.49,224.9 ,265.96, 
     *        5.984, 18.823, 28.440,119.960,153.77,190.42,241.38,284.53, 
     *        8.151, 16.350, 33.460, 45.140,166.73,205.11,246.41,303.07, 
     *       10.484, 19.720, 30.156, 51.354, 65.01,220.41,263.31,309.26,
     *       10.357, 23.400, 35.000, 47.290, 72.50, 88.03,280.99,328.8 ,
     *       12.970, 23.800, 39.900, 53.500, 67.80, 96.7 ,114.27,348.3 ,
     *       15.755, 27.620, 40.900, 59.790, 75.00, 91.3 ,124.0 ,143.46,
     *        4.339, 31.810, 46.000, 60.900, 82.6 , 99.7 ,118.0 ,155.0 ,
     *        6.111, 11.870, 51.210, 67.700, 84.39,109.0 ,128.0 ,147.0 ,
     *        6.560, 12.890, 24.750, 73.900, 92.0 ,111.1 ,138.0 ,158.7 ,
     *        6.830, 13.630, 28.140, 43.240, 99.8 ,120.0 ,140.8 ,168.5 ,
     *        6.740, 14.200, 29.700, 48.000, 65.2 ,128.9 ,151.0 ,173.7 ,
     *        6.763, 16.490, 30.950, 49.600, 73.0 , 90.6 ,161.1 ,184.7 ,
     *        7.432, 15.640, 33.690, 53.000, 76.0 , 97.0 ,119.24,196.46,
     *        7.870, 16.183, 30.652, 54.800, 75.0 , 99.1 ,125.0 ,151.06,
     *        7.860, 17.060, 33.490, 51.300, 79.5 ,102.0 ,129.0 ,157.0 ,
     *        7.635, 18.168, 35.170, 54.900, 75.5 ,108.0 ,133.0 ,162.0 ,
     *        7.726, 20.292, 36.830, 55.200, 79.9 ,103.0 ,139.0 ,166.0 ,
     *        9.394, 17.964, 39.722, 59.400, 82.6 ,108.0 ,134.0 ,174.0 ,
     *        6.000,  20.509,   30.700, 99.99,99.99,99.99,99.99,99.99,  
     *        7.89944,15.93462, 34.058, 45.715,99.99,99.99,99.99,99.99,    
     *        9.7887, 18.5892,  28.351, 99.99,99.99,99.99,99.99,99.99,    
     *        9.750,21.500, 32.000, 99.99,99.99,99.99,99.99,99.99,    
     *       11.839,21.600, 35.900, 99.99,99.99,99.99,99.99,99.99,    
     *       13.995,24.559, 36.900, 99.99,99.99,99.99,99.99,99.99,    
     *        4.175,27.500, 40.000, 99.99,99.99,99.99,99.99,99.99,    
     *        5.692,11.026, 43.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.2171,12.2236, 20.5244,60.607,99.99,99.99,99.99,99.99,    
     *        6.63390,13.13,23.17,34.418,80.348,99.99,99.99,99.99,    
     *        6.879,14.319, 25.039, 99.99,99.99,99.99,99.99,99.99,   
     *        7.099,16.149, 27.149, 99.99,99.99,99.99,99.99,99.99,   
     *        7.280,15.259, 30.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.364,16.759, 28.460, 99.99,99.99,99.99,99.99,99.99,   
     *        7.460,18.070, 31.049, 99.99,99.99,99.99,99.99,99.99,   
     *        8.329,19.419, 32.920, 99.99,99.99,99.99,99.99,99.99,   
     *        7.574,21.480, 34.819, 99.99,99.99,99.99,99.99,99.99,   
     *        8.990,16.903, 37.470, 99.99,99.99,99.99,99.99,99.99,   
     *        5.784,18.860, 28.029, 99.99,99.99,99.99,99.99,99.99,   
     *        7.342,14.627, 30.490,72.3,99.99,99.99,99.99,99.99,   
     *        8.639,16.500, 25.299,44.2,55.7,99.99,99.99,99.99,   
     *        9.0096,18.600, 27.96, 37.4,58.7,99.99,99.99,99.99,   
     *       10.454,19.090, 32.000, 99.99,99.99,99.99,99.99,99.99,   
     *       12.12984,20.975,31.05,45.,54.14,99.99,99.99,99.99,   
     *        3.893,25.100, 35.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.210,10.000, 37.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.580,11.060, 19.169, 99.99,99.99,99.99,99.99,99.99,   
     *        5.650,10.850, 20.080, 99.99,99.99,99.99,99.99,99.99,   
     *        5.419,10.550, 23.200, 99.99,99.99,99.99,99.99,99.99,   
     *        5.490,10.730, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.550,10.899, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.629,11.069, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.680,11.250, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.159,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.849,11.519, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.930,11.670, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.020,11.800, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.099,11.930, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.180,12.050, 23.700, 99.99,99.99,99.99,99.99,99.99,   
     *        6.250,12.170, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.099,13.899, 19.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.000,14.899, 23.299, 99.99,99.99,99.99,99.99,99.99,   
     *        7.879,16.200, 24.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.86404,17.700, 25.000, 99.99,99.99,99.99,99.99,99.99,   
     *        7.870,16.600, 26.000, 99.99,99.99,99.99,99.99,99.99,   
     *        8.500,17.000, 27.000, 99.99,99.99,99.99,99.99,99.99,   
     *        9.100,20.000, 28.000, 99.99,99.99,99.99,99.99,99.99,   
     *        8.95868,18.563,33.227, 99.99,99.99,99.99,99.99,99.99,   
     *        9.220,20.500, 30.000, 99.99,99.99,99.99,99.99,99.99,   
     *       10.430,18.750, 34.200, 99.99,99.99,99.99,99.99,99.99,   
     *        6.10829,20.4283,29.852,50.72,99.99,99.99,99.99,99.99,   
     *        7.416684,15.0325,31.9373,42.33,69.,99.99,99.99,99.99,   
     *        7.285519,16.679, 25.563,45.32,56.0,88.,99.99,99.99,   
     *        8.430,19.000, 27.000, 99.99,99.99,99.99,99.99,99.99,   
     *        9.300,20.000, 29.000, 99.99,99.99,99.99,99.99,99.99,   
     *       10.745,20.000, 30.000, 99.99,99.99,99.99,99.99,99.99,   
     *        4.000,22.000, 33.000, 99.99,99.99,99.99,99.99,99.99,   
     *        5.276,10.144, 34.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.900,12.100, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,   
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,    
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99,     
     *        6.000,12.000, 20.000, 99.99,99.99,99.99,99.99,99.99/

C
C
c      DATA XIFE /8*0.,233.6,262.1/
c      DATA NTOTA /99/
C
C     An element (hydrogen through zinc) can be considered in one of
C     the three following options:
C     1. explicitly - some of energy levels of some of its ionization
C                     states are considered explicitly, ie. their
C                     populations are determined by solving statistical
C                     equilibrium
C     2. implicitly - the atom is assumed not to contribute to
C                     opacity; but is allowed to contribute to the
C                     total number of particles and to the total charge;
C                     the latter is evaluated assuming LTE ionization
C                     balance, ie. by solving a set of Saha equations
C     3. not considered at all
C
C     Input:
C
C     For each element from 1 (hydrogen) to NATOMS, the following
C     parameters:
C
C     MA     =  0  - if the element is not considered (option 3)
C            =  1  - if the element is non-explicit (option 2)
C            =  2  - if the element is explicit (option 1)
C            =  4  - if the element is semi-explicit (i.e. behaves
C                    like MA=2 for continua and MA=1 for lines
C     NA0,NAK - have the meaning only for MA=2; indicate that the
C               explicit energy levels of the present species have
C               the indices between NA0 and NAK (NAK is thus the index
C               of the highest ionization state, which is represented
C               as one-level ion).
C     ION    -  has the meaning for MA=1 only;
C               if ION=0, standard number of ionization degrees is
C                         considered
C                         (counting the neutral state also; so for 
C                         instance to treat all stages of He requires
C                         ION=3, which is a default anyhow).
C               if ION>0, then ION ionization degrees is considered
C     MODPF  -  mode of evaluation of partition functions
C            =  0  -  standard evaluation (see procedure PARTF)
C            >  0  -  partition functions evaluated from the
C                     Opacity Project ionization fraction tables
C            <  0  -  non-standard evaluation, by user supplied
C                     procedure PFSPEC
C     ABN    -  if ABN=0, solar abundance is assumed (given above;
C                         abundance here is assumed as relative
C                         to hydrogen by number
C               if ABN>0, non-solar abundance ABN is assumed; in an
C                         arbitrary scale
C               if ABN<0, non-solar abundance ABN is assumed;
C                         (-ABN times the solar value)
C     PFS    -  see above
C
      READ(IBUFF,*) NATOMS
      WRITE(6,600)
      IAT=0
      IREF=0
      IF(NATOMS.LT.0) NATOMS=-NATOMS
C
      DO 10 I=1,MATOM
         DO 10 J=1,MION0
   10       RR(I,J)=0.
      DO ID=1,ND
         YTOT(ID)=0.
         WMY(ID)=0.
      END DO
C
      DO 20 I=1,MATOM
         TYPAT(I)=DYP(I)
         LGR(I)=.TRUE.
         LRM(I)=.TRUE.
         IATEX(I)=-1
         IF(I.LE.NATOMS) THEN
            IF(MODOLD.EQ.0) THEN
               READ(IBUFF,*) MA,NA0,NAK,ION,MODPF(I),ABN,
     *                       (PFSTD(J,I),J=1,5)
               MA=IABS(MA)
             ELSE
               READ(IBUFF,*) MA,ABN,MODPF(I)
               ION=0
            END IF
          ELSE IF(MOD(IMODE,10).LE.1.and.imode.ne.-4) THEN
            MA=1
            ABN=0.
            ION=0
            MODPF(I)=0
           ELSE
            MA=0
         END IF
         AMAS(I)=D(1,I)
         ABND(I)=D(2,I)
         if(iref.gt.0) abnd(i)=d(2,i)*abnd(iref)/d(2,iref)
         IONIZ(I)=int(D(3,I))
         isemex(i)=0
C
C        increase the standard highest ionization for Teff>30,000 K
C
         IF(TEFF.GT.3.D4) THEN
           IF(I.LE.8) IONIZ(I)=I+1
           IF(I.GT.8.and.i.le.30) IONIZ(I)=9
         END IF
C
         DO J=1,9
            IF(J.LE.8) ENEV(I,J)=xi(J,I)
            if(enev(i,j).ge.enhe2) then
               inpot(i,j)=3
             else if(enev(i,j).ge.enhe1) then
               inpot(i,j)=2
             else
               inpot(i,j)=1
            end if
         END DO
         IF(MA.EQ.0) GO TO 20
         LGR(I)=.FALSE.
         IF(ABN.GT.0) ABND(I)=ABN
         IF(ABN.LT.0) ABND(I)=ABS(ABN)*D(2,I)
         IF(ION.NE.0) IONIZ(I)=ION
         IF(ABN.GT.1.E6) THEN
            READ(IBUFF,*) (ABNDD(I,ID),ID=1,ND)
          ELSE
            DO ID=1,ND
               ABNDD(I,ID)=ABND(I)
            END DO
         END IF
         IF(MA.EQ.1) THEN
            LRM(I)=.FALSE.
            IATEX(I)=0
          ELSE
            IAT=IAT+1
            IATEX(I)=IAT
            if(ma.eq.4) isemex(i)=1
            if(ma.eq.5) isemex(i)=2
            IF(IAT.EQ.IATREF) THEN
               IREF=I
               DO ID=1,ND
                  ABNREF(ID)=ABNDD(I,ID)
               END DO
            END IF
C
C           store parameters for explicit atoms
C
            DO ID=1,ND
               ABUND(IAT,ID)=ABNDD(I,ID)
            END DO
            AMASS(IAT)=AMAS(I)*HMASS
            NUMAT(IAT)=I
            IF(MODOLD.EQ.0) THEN
               N0A(IAT)=NA0
               NKA(IAT)=NAK
            END IF
         END IF
         DO ID=1,ND
            YTOT(ID)=YTOT(ID)+ABNDD(I,ID)
            WMY(ID)=WMY(ID)+ABNDD(I,ID)*AMAS(I)
         END DO
         ABN=ABND(I)/D(2,I)
         IF(MA.EQ.1) WRITE(6,601) I,TYPAT(I),ABND(I),ABN
         IF(MA.EQ.2) WRITE(6,602) I,TYPAT(I),ABND(I),ABN,IAT,NA0,NAK
   20 CONTINUE
      IF(MOD(IMODE,10).LE.1) NATOMS=MATOM
      DO ID=1,ND
         WMM(ID)=WMY(ID)*HMASS/YTOT(ID)
      END DO
      DO JJ=1,NATOMS
         DO ID=1,ND
            RELAB(JJ,ID)=1.
         END DO
      END DO
C
      IF(ICHEMC.NE.1) go to 100
C
C     abundance change with respect to the model atmosphere input
C     (unit 5);
C     this option is switched on by the parameter ICHEMC (read from
C     unit 55), if it is non-zero, an additional input from
C     unit 56 is required
C
C     unit 56 input:
C
C     NCHANG  -  number of chemical elements for which the abundances
C                are going to be changes;
C
C     then there are NCHANG records, each contains:
C
C     I       - atomic number
C     ABN     - new abundance; coded using the same conventions as in
C               the standard input
C
      READ(56,*,ERR=566,END=566) NCHANG
      WRITE(6,610)
      DO 60 II=1,NCHANG
         READ(56,*) I,ABN
         ABND(I)=D(2,I)
         IF(ABN.GT.0) ABND(I)=ABN
         IF(ABN.LT.0) ABND(I)=-ABN*D(2,I)
         if(abn.gt.1.) abnd(i)=10.**(abn-12.)
         IF(ABN.GT.1.E6) THEN
            READ(56,*) (ABNDD(I,ID),ID=1,ND)
          ELSE
            DO ID=1,ND
               ABNDD(I,ID)=ABND(I)
            END DO
         END IF
         LGR(I)=.FALSE.
         IATX=IATEX(I)
         IF(IATX.GT.0) THEN
            DO ID=1,ND
               RELAB(IATX,ID)=ABNDD(I,ID)/ABUND(IATX,ID)
               ABUND(IATX,ID)=ABNDD(I,ID)
            END DO
         END IF
         ABNR=ABND(I)/D(2,I)
         WRITE(6,601) I,TYPAT(I),ABND(I),ABNR
   60 CONTINUE
C
C     renormalize abundances to have the standard element abundance
C     equal to unity
C
  100 IF(IREF.LE.1) RETURN
      write(6,620)
      DO 30 I=1,MATOM
         IAT=IATEX(I)
         IF(IAT.LT.0) GO TO 30
         DO ID=1,ND
            ABNDD(I,ID)=ABNDD(I,ID)/ABNREF(ID)
            YTOT(ID)=YTOT(ID)+ABNDD(I,ID)
            WMY(ID)=WMY(ID)+ABNDD(I,ID)*AMAS(I)
         END DO
         ABNR=ABND(I)/D(2,I)
         IF(IAT.EQ.0) THEN
            WRITE(6,601) I,TYPAT(I),ABND(I),ABNR
          ELSE
            DO ID=1,ND
               ABUND(IAT,ID)=ABNDD(I,ID)
            END DO 
            WRITE(6,602) I,TYPAT(I),ABND(I),ABNR,IAT,N0A(IAT),NKA(IAT)
         END IF
   30 CONTINUE
      DO ID=1,ND
         WMM(ID)=WMY(ID)*HMASS/YTOT(ID)
      END DO
      RETURN
  566 WRITE(6,656)
      STOP
c
  600 FORMAT(1H0//' CHEMICAL ELEMENTS INCLUDED'/
     *            ' --------------------------'//
     * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
     * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)
  601 FORMAT(1H ,I4,3X,A5,1P2E14.2)
  602 FORMAT(1H ,I4,3X,A5,1P2E14.2,3X,
     *       'EXPLICIT: IAT=',I3,'  N0A=',I3,'  NKA=',I3)
  610 FORMAT(//'    CHEMICAL ELEMENTS INCLUDED - CHANGED (unit 56)'
     *           /'    --------------------------'//
     * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
     * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)
  620 FORMAT(1H0//'    CHEMICAL ELEMENTS INCLUDED - RENORMALIZATION'/
     *            '    --------------------------'//
     * ' NUMBER  ELEMENT           ABUNDANCE'/1H ,16X,
     * 'A=N(ELEM)/N(H)  A/A(SOLAR)'/)
  656 FORMAT(//' CHEMICAL COMPOSITION COULD NOT BE READ FROM ',
     *       'UNIT 56'//' STOP.')
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE INIMOD
C
C   SET UP COMMON/RRRVAL/  - VALUES OF  N(ION)/U(ION) FOR ALL THE ATOMS
C   AND IONS CONSIDERED
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      COMMON/HPOPST/HPOP
C
c     1. "low-temperature" ionization fractions 
c         (using Hamburg partition functions)
c
      DO 50 ID=1,ND
         IF(IFMOL.EQ.0.OR.TEMP(ID).GE.TMOLIM) THEN
            CALL STATE(ID,TEMP(ID),ELEC(ID),S1)
            HPOP=DENS(ID)/WMM(ID)/YTOT(ID)
            DO J=1,MION0
               DO I=1,MATOM
                  RRR(ID,J,I)=RR(I,J)*HPOP
               END DO
            END DO
            DO IAT=1,NATOM
               ATTOT(IAT,ID)=HPOP*ABUND(IAT,ID)
            END DO
          ELSE
            HPOP=ATTOT(1,ID)
         END IF
         IF(ID.NE.IDSTD) GO TO 50
         TSTD=TEMP(ID)
         VTS=VTURB(ID)
         DSTD=SQRT(1.4E7*TSTD+VTS)
         WRITE(6,601) ID,TEMP(ID),ELEC(ID),hpop
c        DO I=1,MATOM
         DO I=1,30    
            WRITE(6,602) TYPAT(I),(RRR(ID,J,I),J=1,MION0-1)
         END DO
c        WRITE(6,603)
c        DO I=1,MATOM
c           WRITE(6,602) TYPAT(I),(PFSTD(J,I),J=1,MION0-1)
c        END DO
   50 CONTINUE
c      write(*,*) 'in inimod'
c     id=idstd
c     do ie=1,nion
c        i=nfirst(ie)
c        j=nnext(ie)
c        ia=iatm(i)
c        write(6,702) ie,ia,i,j,abund(ia,id),popul(i,id),popul(j,id)
c     end do
c 702 format(4i4,1p3e12.5)
c
c     2. "high-temperature" ionization fractions 
c         (using the Opacity Project ionization fractions)
c
      if(teff.lt.0.) then
      CALL FRAC1
      ID=IDSTD
      HPOP=DENS(ID)/WMM(ID)/YTOT(ID)
      WRITE(6,604) ID,TEMP(ID),ELEC(ID)
      DO 60 I=1,MATOM
         WRITE(6,605) TYPAT(I),(RRR(ID,J,I)/hpop,J=1,MION)
         ioniz(i)=i+1
   60 continue
      end if
C
  601 FORMAT(/' N/U  AT THE STANDARD DEPTH  (ID =',I3, 
     *         ' ; T,Ne = ',F8.1,1P2E12.3,' )'/
     *         ' --------------------------'//)
  602 FORMAT(1H ,A4,1P8E9.2)
c 603 FORMAT(//' PARTITION FUNCTIONS AT THE STANDARD DEPTH'/
c    *        ' ------------------------------------------'//)
  604 FORMAT(/' N/U  AT THE STANDARD DEPTH  - OP DATA',
     * '  (ID =',I3,' ; T,Ne = ',F8.1,1PE12.3,' )'//)
  605 FORMAT(1H ,A4,(1P8E9.2))
      RETURN
      END
C
C
C ********************************************************************
C
      SUBROUTINE STATE(ID,TE,ANE,Q)
C
C     modified LTE Saha equations - possibly using 
C     radiation temperatures after 
C     Schaerer and Schmutz AA 288, 321, 1994
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'WINCOM.FOR'
      common/moltst/pfmol(500,mdepth),anmol(500,mdepth),
     *              pfato(100,mdepth),anato(100,mdepth),
     *              pfion(100,mdepth),anion(100,mdepth)
      dimension FFI(MION0)
C
      Q=0.
      DO 50 I=1,NATOMS
         IF(LGR(I)) GO TO 50
         ION=IONIZ(I)
         RQ=0.
         RS=1.
         T=TRAD(INPOT(I,1),ID)
         if(t.le.0.) t=te
         X=SQRT(T/ANE)
         XMX=2.145E4*SQRT(X)
         CALL PARTF(I,1,T,ANE,XMX,UM)
         PFSTD(1,I)=UM
         JMAX=1
         DO J=2,ION
            J1=J-1
            T=TRAD(INPOT(I,J),ID)
            if(t.le.0.) t=te
            TLN=LOG(T)*1.5
            TK=BOLK*T
            THL=11605./T
            X=SQRT(T/ANE)
            XMX=2.145E4*SQRT(X)
            DCH=EH/XMX/XMX/TK
            DCHT=DCH*J1
            FI=36.113+TLN-THL*ENEV(I,J1)+DCHT
            X=J
            XMAX=XMX*SQRT(X)
            CALL PARTF(I,J,T,ANE,XMAX,U)
            PFSTD(J,I)=U
            FI=EXP(FI)*U/UM/ANE
            FFI(J)=FI
            IF(FFI(J).GT.1.) JMAX=J
            UM=U
         END DO
         IF(JMAX.LT.ION) THEN
            R=1.
            RQ=JMAX-1
            DO J=JMAX+1,ION
               R=R*FFI(J)
               RR(I,J)=R/PFSTD(J,I)
               RS=RS+R
               RQ=RQ+(J-1)*R
            END DO
         END IF
         IF(JMAX.GT.1) THEN
            R=1.
            DO JJ=1,JMAX-1
               J=JMAX-JJ
               R=R/FFI(J+1)
               RR(I,J)=R/PFSTD(J,I)
               RS=RS+R
               RQ=RQ+(J-1)*R
            END DO
         END IF
         ABND(I)=ABNDD(I,ID)
         RR(I,JMAX)=ABND(I)/RS
         DO J=1,ION
            IF(J.NE.JMAX) RR(I,J)=RR(I,J)*RR(I,JMAX)
            if(rr(i,j).lt.1.e-35) rr(i,j)=0.
         END DO
         RR(I,JMAX)=RR(I,JMAX)/PFSTD(JMAX,I)
         X=RQ/RS
c        IF(LRM(I)) GO TO 50
         if(i.gt.1) Q=X*ABND(I)+Q
         anato(i,id)=rr(i,1)*pfstd(1,i)
         pfato(i,id)=pfstd(1,i)
         anion(i,id)=rr(i,2)*pfstd(2,i)
         pfion(i,id)=pfstd(2,i)
   50 CONTINUE
c
      do imol=1,500
         anmol(imol,id)=0.
         pfmol(imol,id)=0.
      end do
c
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE TINT
C
C     LOGARITHMIC INTERPOLATION COEFFICIENTS FOR INTERPOLATION OF
C     TEMP(ID) TO THE VALUES  5000,10000,20000,40000
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION TT(4)
      DATA TT /3.699, 4.000, 4.301, 4.602/
C
      DO ID=1,ND
         T=LOG10(TEMP(ID))
         J=3
         IF(T.GT.TT(3)) J=4
         JT(ID)=J
         X=(TT(J)-TT(J-1))*(TT(J)-TT(J-2))*(TT(J-1)-TT(J-2))
         TI0(ID)=(T-TT(J-2))*(T-TT(J-1))*(TT(J-1)-TT(J-2))/X
         TI1(ID)=(T-TT(J-2))*(TT(J)-T)*(TT(J)-TT(J-2))/X
         TI2(ID)=(T-TT(J-1))*(T-TT(J))*(TT(J)-TT(J-1))/X
      ENd dO
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE INIBL0
C
C     AUXILIARY INITIALIZATION PROCEDURE
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      parameter (un=1.)
      character*2 iu
      character*6 ilab
      DIMENSION CROSS(MCROSS,MFRQ),
     *          ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ),
     *          ABSOC(MFREQC),EMISC(MFREQC),SCATC(MFREQC)
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/lasers/lasdel
      common/linrej/ilne(mdepth),ilvi(mdepth)
      common/velaux/velmax,iemoff,nltoff,itrad
      common/alsave/ALAM0s,ALASTs,CUTOF0s,CUTOFSs,RELOPs,SPACEs
C
C --------------------------------------------------------------
C Parameters controlling an evaluation of the synthetic spectrum
C
C --------------------------------------------------------------
C
C     ALAM0, ALAM1 - synthetic spectrum is evaluated between wavelengths
C                    ALAM0 (initial) and ALAM1 (final), given in Anstroms
C     CUTOF0       - cutoff parameter for normal lines (given in Angstroms)
C                    ie the maximum distance from the line center, in
C                    which the opacity in the line is allowd to contribute
C                    to the total opacity (recommended 5 - 10)
C     CUTOFS  = SPACON
C     SPACON       - spacing of the continuum wavelength points 
C                    (at the midpoint of teh total interval; actual spacing
C                    is equidistant in log(lambda)
C     RELOP        - the minimum value of the ratio (opacity in the line
C                    center)/(opacity in continuum), for which is the line
C                    taken into account (usually 1d-4 to 1d-3)
C     SPACE        - the maximum distance of two neighbouring frequency
C                    points for evaluating the spectrum; in Angstroms
C
C     INLTE        = 0  -  pure LTE (no line in NLTE)
C                  ne.0 -  NLTE option, ie one or more lines treated
C                          in the exact or approximate NLTE approach
C     INLIST       = 1  -  line list is read in the "original" format
C                          (see procedure INISET)
C                  ne.1 -  line list is read in the "new" format
C                          (see procedure INISET)
C     IFHE2        gt.0 -  He II line opacity in the first four series
C                          (Lyman, Balmer, Paschen, Brackett)
C                          for lines with lambda < 3900 A
C                          is taken into account even if line list
C                          does not contain any He II lines (i.e.
C                          He II lines are treated as the hydrogen lines)
C
C     IHYDPR      = 0  - means that hydrogen lines Stark profiles
C                        are calculated by approximate formulae
C                 > 0  - hydrogen lines Stark profiles are calculated
C                        in detail, using the Schoening & Butler tables;
C                        (for 1-2 to 1-5; 2-3 to 2-10).
C                        the tables are stored in file FOR0xx.dat,
C                        where xx=IHYDPR;
C                        higher Balmer lines are calculated as before
C
C     the meaning of other parameters is quite analogous, for the
C     following lines
C
C     IHE1PR    - He I lines at 4471, 4026, 4387, and 4922 Angstroms
C                 (tables calculated by Barnard, Cooper, and Shamey)
C     IHE2PR    - for the He II lines calculated by Schoening and Butler,
C
      if(ifeos.eq.0) then
      READ(55,*) IFREQ,INLTE,ICONTL,INLIST,IFHE2
      IF(LTE) INLTE=0
      READ(55,*) IHYDPR,IHE1PR,IHE2PR
      READ(55,*) ALAM0,ALAST,CUTOF0,CUTOFS,RELOP,SPACE
      end if
C
      IF(IDSTD.EQ.0) THEN
         ID1=5
         NDSTEP=(ND-2*ID1)/2
         IDSTD=2*ND/3
       ELSE IF(IDSTD.LT.0) THEN
         ID1=1
         NDSTEP=-IDSTD
         IDSTD=2*ND/3
      END IF
      if(imode.le.-3) ndstep=1
c
      alam0s=alam0
      alasts=alast
      cutof0s=cutof0
      cutofss=cutofs
      relops=relop
      spaces=space
C
C     if ALAST.lt.0 - set up vacuum wavelengths everywhere
C
      vaclim=2000.
      if(alast.lt.0.) then
         alast=abs(alast)
         alasts=alast
         vaclim=1.e18
      end if
c
      if(inlte.lt.10) then
         lasdel=.true.
       else if(inlte.le.20) then
         inlte=inlte-10
         lasdel=.false.
       else if(inlte.le.30) then
         inlte=inlte-20
         ifreq=11
         lasdel=.true.
       else if(inlte.le.40) then
         inlte=inlte-30
         ifreq=11
         lasdel=.false.
      end if
C
      ibin(0)=mod(inlist,10)
      do ilist=1,mmlist
         tmlim(ilist)=tmolim
         ibin(ilist)=mod(inlist,10)
         ivdwli(ilist)=0
         iun=19+ilist
         write(iu,622) iun
  622    format(i2)
         amlist(ilist)='fort.' // iu
      end do
c
      if(imode.ge.-3.and.imode.le.1) then
      nmlist=0
      numlis=0
      read(55,*,err=5,end=5) nmlist,(iunitm(ilist),ilist=1,nmlist)
      do ilist=1,nmlist
         write(iu,622) iunitm(ilist)
         amlist(ilist) ='fort.' // iu
      end do
    5 continue
c
      ilist=0
      amlist(0)='fort.19'
      read(3,*,err=20,end=20) amlist(0),ibin(0)
c
      ilist=0
   10 continue
      ilist=ilist+1
      read(3,*,end=20) amlist(ilist),ibin(ilist),tmlim(ilist),
     *                 ivdwli(ilist)
      write(6,610) ilist,amlist(ilist),ibin(ilist),tmlim(ilist),
     *             ivdwli(ilist)
  610 format(i3,a40,i4,f8.1,i4)
      numlis=numlis+1
   20 continue
      if(numlis.gt.0) nmlist=numlis
      if(nmlist.gt.0.and.ifmol.eq.0) then
         write(*,*) 'NEEDS TO SET IFMOL > 0 with NMLIST>0' 
         stop
      end if
c
      ilist=0
      ilab='ATOMIC'
      write(6,623) ilist,ilab,trim(amlist(ilist)),ibin(ilist)
      ilab='MOLEC '
      do ilist=1,nmlist
         write(6,624) ilist,ilab,trim(amlist(ilist)),ibin(ilist),
     *                ivdwli(ilist),tmlim(ilist)
      end do
  623 format(/'************************'/
     *        ' LINE LISTS:'/
     *       /' ILIST',8x,'FILENAME  IBIN  IVDWLI  TMLIM'/
     *        i4,2x,a6,2x,a,2x,i4,i6,f11.1)
  624 format( i4,2x,a6,2x,a,2x,i4,i6,f11.1)
      end if 
c
C
c     VTB    - turbulent velocity (in km/s). In non-negative, this
C              value overwrites the value given by the standard input
C
      read(55,*,err=30,end=30) VTB
      if(ifwin.le.0) then
         if(vtb.ge.0.) then
            WRITE(6,608) VTB
  608       FORMAT(//' TURBULENT VELOCITY  -  CHANGED TO   VTURB =',
     *      1PE10.3,'  KM/S'/' ------------------'/)
            do id=1,nd
               vturb(id)=vtb*vtb*1.e10
            end do
         end if
      end if
C
      TSTD=TEMP(IDSTD)
      VTS=VTURB(IDSTD)
      DSTD=SQRT(1.4E7*TSTD+VTS)
   30 continue
C
C     angle points (in case the specific intensities are evaluated
C
C     NMU0      -    number of angles:
C               >0 - and if also ANG0>0, angles (mu's) equidistant 
C                    between 1 and ANG0
C               >0 - and if also ANG0<0, angles (mu's) equidistant 
C                    between 0.7 and ANG0, and sinuses equidistatnt for
C                    others
C               <0 - angles read in the next record
C     ANG0      -    minimum mu (see above)
C     IFLUX     -    mode for evaluating angle-dependent intensities and
C                    the corresponding flux:
C               =0 - no specifiec intensities are evaluated; only usual
C                    flux is stored (unit 7 and 17)
C               =1 - specific intensities are evaluated; 
C                    and stored on unit 18
C               =2 - (interesting only for the case of macroscopic
C                    velocity field); specific intensities evaluated by
C                    a simple formal solution (RESOLV)
C
      NMU0=1
      ANG0=1.
      ANGL(1)=1.
      WANGL(1)=0.
      IFLUX=0
      velmax=3.e5
      nltoff=0
      iemoff=0
      itrad=0
      do id=1,nd
        wdil(id)=un
      end do
      if(ifwin.le.0) then
      READ(55,*,end=100,err=100) NMU0,ANG0,IFLUX
C
C     determinantion of the angle points and weights
C
      IF(NMU0.LT.0) THEN
         NMU0=IABS(NMU0)
         READ(55,*) (ANGL(IMU),IMU=1,NMU0)
         DO IMU=2,NMU0-1
            WANGL(IMU)=0.5*(ANGL(IMU-1)+ANGL(IMU+1))
         END DO
         WANGL(1)=0.5*(ANGL(1)-ANGL(2))
         WANGL(NMU0)=0.5*(ANGL(NMU0-1)-ANGL(NMU0))
       ELSE
         IF(ANG0.GT.0.) THEN
            IF(NMU0.GT.1) THEN
            DMU=(1.-ANG0)/(NMU0-1)
            DO IMU=1,NMU0
               ANGL(IMU)=1.-(IMU-1)*DMU
               WANGL(IMU)=DMU
            END DO
            WANGL(1)=0.5*DMU
            WANGL(NMU0-1)=0.5*DMU
            WANGL(NMU0)=2.*DMU
            END IF 
          ELSE
            ANGH=0.70710678
            DMU=ANGH/(NMU0-1)
            DO IMU=1,NMU0
               ANGL(IMU)=(IMU-1)*DMU
               ANGL(IMU)=SQRT(1.-ANGL(IMU)**2)
               IF(IMU.GT.1.AND.IMU.LT.NMU0) 
     *            WANGL(IMU)=0.5*(ANGL(IMU-1)+ANGL(IMU+1))
            END DO
            WANGL(1)=0.5*(ANGL(1)-ANGL(2))
            WANGL(NMU0)=0.5*(ANGL(NMU0-1)-ANGL(NMU0))
            IF(ANG0.LT.0.) DMU=(ANGH+ANG0)/(NMU0-1)
            DO IMU=1,NMU0-2
               ANGL(IMU+NMU0)=ANGH-IMU*DMU
               WANGL(IMU+NMU0)=DMU
            END DO
            WANGL(NMU0)=WANGL(NMU0)+0.5*DMU
            WANGL(2*NMU0-3)=0.5*DMU
            WANGL(2*NMU0-2)=2.*DMU
            NMU0=2*NMU0-2
         END IF
      END IF
      IF(NMU0.LE.0) GO TO 100
      WRITE(6,609) NMU0,(ANGL(I),I=1,NMU0)
  609 FORMAT(//' SPECIFIC INTENSITIES COMPUTED FOR',I3,
     *         ' ANGLES  mu=cos(theta) ='/
     *         ' ---------------------------------',
     *         '------------------------'//
     *         (10F7.2))
  100 CONTINUE
      else
      itrad=1
      read(55,*,end=110,err=110) velmax,ITRAD,nltoff,iemoff
  110 write(6,602) velmax,itrad,nltoff,iemoff
      if(velmax.lt.0.) then
      velmax=3.e5
      go to 120
      end if
  602 format(//' velmax (velocity for line rejection)',
     *       ' itrad,nltoff,iemoff',f10.1,2i3)
C
C     Set up rays and weights
C
      call velset
      call radtem
      CALL SETRAY
      CALL WGTJH1
C
      end if
C
  120 CONTINUE
      velmax=velmax*1.e5
      do id=1,nd
         ilvi(id)=0
         ilne(id)=0
         if(vel(id).gt.velmax.and.iemoff.eq.0) ilvi(id)=1
         if(vel(id).gt.velmax.and.nltoff.gt.0.and.iemoff.gt.0)
     *      ilne(id)=1
      end do
C
      IF(IMODE.EQ.-1) THEN
         INLTE=0
         CUTOF0=0.
      END IF
C
C     continuum frequencies
C
      if(ifwin.le.0) then
         alam0=alam0s
         if(alam0s.eq.0.) alam0=5.e7/temp(1)/10.
         if(alam0s.lt.0.) alam0=-5.e7/temp(1)/alam0s
         alast=alasts
         if(alasts.eq.0.) alast=5.e7/temp(1)*20.
         if(alasts.lt.0.) alast=-5.e7/temp(1)*alasts
c        if(alast.gt.1.e5) alast=1.e5
         ALAMC=(ALAM0+ALAST)*0.5
         if(space.eq.0.) space=4.3e-8*sqrt(temp(idstd))*alamc
         if(space.lt.0.) space=-5.72e-8*sqrt(temp(idstd))*alamc*space
         SPACF=2.997925E18/ALAMC/ALAMC*SPACE
         WRITE(6,601) ALAM0,ALAST,CUTOF0,RELOP,SPACF,SPACE
         CUTOF0=0.1*CUTOF0
         SPACE0=SPACE*0.1
         ALAM0=1.D-1*ALAM0
         ALAST=1.D-1*ALAST
         ALAMC=ALAMC*0.1
         ALST00=ALAST
         FRLAST=2.997925D17/ALAST
         NFREQ=2
         FREQ(1)=2.997925D17/ALAM0
         FREQ(2)=FRLAST
C
      else
C
         spacon=cutofs
         IF(SPACON.EQ.0) SPACON=3.
         XFR=(ALAST-ALAM0)/SPACON 
         NFREQC=int(XFR)+1
         NFREQC=MIN(NFREQC,MFREQC)
         NFREQC=MAX(NFREQC,2)
         DLAMLO=LOG10(ALAST/ALAM0)/(NFREQC-1)
         AL0L=LOG10(ALAM0)
         alambe=alam0  
         DO IJ=1,NFREQC
            AL=AL0L+(IJ-1)*DLAMLO
            ALAM=EXP(2.3025851*AL)
            WLAMC(IJ)=ALAM
            FREQC(IJ)=2.997925E18/ALAM
         END DO
         ALAMC=(ALAM0+ALAST)*0.5
         SPACF=2.997925E18/ALAMC/ALAMC*SPACE
         WRITE(6,601) ALAM0,ALAST,CUTOF0,RELOP,SPACF,SPACE
         CUTOF0=0.1*CUTOF0
         SPACE0=SPACE*0.1
         ALAM0=1.D-1*ALAM0
         ALAST=1.D-1*ALAST
         ALAMC=ALAMC*0.1
         ALST00=ALAST
         FRLAST=2.997925D17/ALAST
         NFREQ=2
         FREQ(1)=2.997925D17/ALAM0
         FREQ(2)=FRLAST
c
      end if
c
      CALL SIGAVS
      IF(IHYDPR.NE.0) THEN
         CALL HYDINI
         CALL XENINI
      END IF
      IF(IHE1PR.GT.0) CALL HE1INI
      IF(IHE2PR.GT.0) CALL HE2INI
C
C     auxiliary quantities for dissolved fractions
C
      DO ID=1,ND
         CALL DWNFR0(ID)
         CALL WNSTOR(ID)
      END DO
C
c     pretabulate expansion coefficients for the Voigt function
c
      CALL PRETAB
c
c     calculate the characteristic standard opacity
c
      IF(IMODE.LE.2) THEN
         if(ifwin.le.0.and.ndstep.eq.0) then
c
c        old procedure 
c
            CALL CROSET(CROSS)
            DO ID=1,ND
               CALL OPAC(ID,CROSS,ABSO,EMIS,SCAT)
               ABSTD(ID)=MIN(ABSO(1),ABSO(2))
            END DO
          else
c
c       new procedure
c
            if(ifwin.le.0) then
               nfreqc=ifix(real(cutofs,4))
               if(nfreqc.eq.0) nfreqc=mfreq
               all0=log(alam0)
               all1=log(alast)
               dlc=(all1-all0)/(nfreqc-1)
               do ijc=1,nfreqc
                  wlamc(ijc)=exp(all0+(ijc-1)*dlc)
                  freqc(ijc)=2.997925e17/wlamc(ijc)
               end do
               CALL CROSEW(CROSS)
               do id=1,nd
                  CALL OPACON(ID,CROSS,ABSOC,EMISC,SCATC)
                  do ijc=1,nfreqc
                     abstdw(ijc,id)=absoc(ijc)
                  end do
               end do

c         write(*,*) 'abstdw(1,ij)',(abstdw(ij,1),ij=1,nfreqc)
c         write(*,*) 'abstdw(50,ij)',(abstdw(ij,50),ij=1,nfreqc)
c
             else
               CALL CROSEW(CROSS)
               DO ID=1,ND
                  CALL OPACW(ID,CROSS,ABSO,EMIS,ABSOC,EMISC,SCATC,0)
                  DO IJ=1,NFREQC
                     ABSTDW(IJ,ID)=ABSOC(IJ)/DENSCON(ID)
                  END DO
               END DO
            end if
         end if
      END IF
C
  601 FORMAT(//'----------------------------------------------'/
     *           ' BASIC INPUT PARAMETERS FOR SYNTHETIC SPECTRA'/
     *           ' ---------------------------------------------'/
     *           ' INITIAL LAMBDA',28X,1H=,F10.3,' ANGSTROMS'/
     *           ' FINAL   LAMBDA',28X,1H=,F10.3,' ANGSTROMS'/
     *           ' CUTOFF PARAMETER',26X,1H=,F10.3,' ANGSTROMS'/
     *           ' MINIMUM VALUE OF (LINE OPAC.)/(CONT.OPAC) =',1PE10.1/
     *           ' MAXIMUM FREQUENCY SPACING',17X,1H=,1PE10.3,'  I.E.',
     *             0PF6.3,'  ANGSTROMS'/
     *           ' ---------------------------------------------'/)
c
      write(6,612) idstd,ndstep
  612 format(/'IDSTD, NDSTEP = ',2i5/)
      RETURN
      END
C
C ***********************************************************************
C

      SUBROUTINE INIBL1(IGRD)
C     =======================
C
C     AUXILIARY INITIALIZATION PROCEDURE
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/alsave/ALAM0s,ALASTs,CUTOF0s,CUTOFSs,RELOPs,SPACEs
      common/plaopa/plalin,plcint,chcint
      common/conabs/absoc(mfreqc),emisc(mfreqc),scatc(mfreqc),
     *          plac(mfreqc)
      parameter (un=1.,bnc=1.4743e-2,hkc=4.79928e4,
     *          clc=2.997925e17)
      DIMENSION CROSS(MCROSS,MFRQ),
     *          ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ)
C
C     auxiliary quantities for dissolved fractions
C
      DO ID=1,ND
         CALL DWNFR0(ID)
         CALL WNSTOR(ID)
      END DO
      CALL TINT
c
c     reset wavelengths in case of opacity grid calculations
c
      if(igrd.ge.0) then
         alam0=alam0s
         if(alam0s.eq.0.) alam0=5.e7/temp(1)/10.
         if(alam0s.lt.0.) alam0=-5.e7/temp(1)/alam0s
         alast=alasts
         if(alasts.eq.0.) alast=5.e7/temp(1)*20.
         if(alasts.lt.0.) alast=-5.e7/temp(1)*alasts
c        if(alast.gt.1.e5) alast=1.e5
         cutof0=cutof0s
         cutofs=cutofss
         relop=relops
         if(relops.eq.0) then
            relop=1.e-15
            if(temp(1).lt.2.e6) relop=1.e-6
            if(temp(1).lt.1.e6) relop=1.e-5
            if(temp(1).lt.1.e5) relop=1.e-4
         end if
         space=spaces
         ALAMC=(ALAM0+ALAST)*0.5
         if(space.eq.0.) space=4.3e-8*sqrt(temp(idstd))*alamc
         if(space.lt.0.) space=-5.72e-8*sqrt(temp(idstd))*alamc*space
         SPACF=2.997925E18/ALAMC/ALAMC*SPACE
         CUTOF0=0.1*CUTOF0
         SPACE0=SPACE*0.1
         ALAM0=1.D-1*ALAM0
         ALAST=1.D-1*ALAST
         ALAMC=ALAMC*0.1
         ALST00=ALAST
         FRLAST=CLC/ALAST
c
         nfreqc=ifix(real(cutofs,4))
         if(nfreqc.eq.0) nfreqc=mfreq
         all0=log(alam0)
         all1=log(alast)
         dlc=(all1-all0)/(nfreqc-1)
         xcc0=hkc/temp(1)
c        write(6,641) nfreqc
c 641    format(' continuum freq ',i4)
         do ijc=1,nfreqc
            wlamc(ijc)=exp(all0+(ijc-1)*dlc)
            freqc(ijc)=clc/wlamc(ijc)
c           frc=freqc(ijc)*1.e-15
c           plac(ijc)=bnc*frc**3/(exp(xcc0*frc)-un)
         end do
         id=1
         CALL CROSEW(CROSS)
         CALL OPACON(ID,CROSS,ABSOC,EMISC,SCATC)
         wc0=(freqc(1)-freqc(2))*0.5
         wc1=(freqc(nfreqc-1)-freqc(nfreqc))*0.5
         do ijc=2,nfreqc-1
            absoc(ijc)=min(absoc(ijc),1.e30)
            write(26,642) wlamc(ijc)*10.,log(absoc(ijc)/dens(1))
         end do
  642    format(f11.3,1p5e13.5)
c
         do ijc=1,nfreqc
            abstdw(ijc,id)=absoc(ijc)
         end do
c
      end if
c
c     calculate the characteristic standard opacity
c
      IF(IMODE.LE.2.and.imode.ge.-2) THEN
         if(ifwin.le.0) then
            CALL CROSET(CROSS)
            DO ID=1,ND
               CALL OPAC(ID,CROSS,ABSO,EMIS,SCAT)
               ABSTD(ID)=MIN(ABSO(1)+SCAT(1),ABSO(2)+SCAT(2))
            END DO
          else
            CALL CROSEW(CROSS)
            DO ID=1,ND
               CALL OPACW(ID,CROSS,ABSO,EMIS,ABSOC,EMISC,SCATC,0)
               DO IJ=1,NFREQC
                  denscon(id)=1.
                  ABSTDW(IJ,ID)=ABSOC(IJ)/DENSCON(ID)
               END DO
            END DO
         end if
      END IF
C
      RETURN
      END

C
C ***********************************************************************
C
      SUBROUTINE RESOLV
C
C     driver for evaluating opacities and emissivities which then
C     enter the solution of the radiative transfer equation
C     (RTE or RTEDFE)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      DIMENSION CROSS(MCROSS,MFRQ),
     *          ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ)
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)
      COMMON/HPOPST/HPOP
C
      IHYL=-1
c     imode0=imode
c
c     if(imode.le.-3) call abnchn(1)
C
C     set up the partial line list for the current interval
C
      CALL INISET
      if(ifmol.gt.0) then
         do ilist=1,nmlist
            call molset(ilist)
         end do
      end if
C   
C     select possible hydrogen lines that may contribute to the opacity
C
      IF(IMODE.NE.-1) CALL HYLSET
C   
C     select possible He II lines that may contribute to the opacity
C
      IF(IMODE.NE.-1) CALL HE2SET
C
C     output of information about selected lines
C
      CALL INIBLA
      if(ifmol.gt.0) call iniblm
C
C     photoinization cross-sections
C
      CALL CROSET(CROSS)
C     
C     monochromatic opacity and emissivity including all contributing
C     lines and continua
C
c       write(*,*) 'resolv',nd,nfreq,3.e18/freq(1),3.e18/freq(nfreq)
      IF(IMODE.GE.-1) THEN
         DO ID=1,ND
            CALL OPAC(ID,CROSS,ABSO,EMIS,SCAT)
            ABSTD(ID)=0.5*(ABSO(1)+ABSO(2))
            DO IJ=1,NFREQ
               CH(IJ,ID)=ABSO(IJ)
               ET(IJ,ID)=EMIS(IJ)
               SC(IJ,ID)=SCAT(IJ)
c              SC(IJ,ID)=SCAT(IJ)+ELEC(ID)*SIGE
            END DO
c         write(*,*) abso(1),abso(nfreq/2),abso(nfreq)
            if(imode0.eq.-4) call ougrid(abso,scat)
         END DO
C
C        output of information about selected hydrogen lines
C
         CALL INIBLH
C
C        the "iron curtain option - output of monochromatic opacities
C        on Unit 27
C        the actual stored quantities are:
C        WLAM(IJ)      - wavelength
C        ABSO(IJ)/HPOP - opacity normalized per one hydrogen atom
C
       ELSE IF(IMODE.EQ.-2) THEN
         ID=1
         write(27,626) temp(id),dens(id),elec(id)
         CALL OPAC(ID,CROSS,ABSO,EMIS,SCAT)
         DO IJ=3,NFREQ-1
            ABSO(IJ)=(ABSO(IJ)+SCAT(IJ))/HPOP
            WRITE(27,627) WLAM(IJ),ABSO(IJ),scat(ij)
         END DO
       else
         id=1
         call opac(id,cross,abso,emis,scat)
         ch(1,id)=abso(1)
         ch(2,id)=abso(2)
c        imode=imode0
         call ougrid(abso,scat)
      END IF
  626 format(1p3e15.4)
  627 format(f15.3,1p2e15.5)
      RETURN
      END
C
C *******************************************************************
C
      SUBROUTINE RTE
C
C     solution of the radiative transfer equation by Feautrier method
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      DIMENSION D(3,3,MDEPTH),ANU(3,MDEPTH),AANU(MDEPTH),DDD(MDEPTH),
     *       AA(3,3),BB(3,3),CC(3,3),VL(3),AMU(3),WTMU(3),
     *       DT(MDEPTH),TAU(MDEPTH),
     *       RDD(MDEPTH),FKK(MDEPTH),ST0(MDEPTH),SS0(MDEPTH),
     *       RINT(MDEPTH,MMU)
      CHARACTER*4 TYPION(9) 
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      COMMON/CTRFUN/CINT1(MDEPTH),CINT2(MDEPTH),
     *   CTRI(MDEPTH),CTRR(MDEPTH),XKAR(MDEPTH),
     *   ABXLI(MFREQ),EMXLI(MFREQ),IJCTR(MFREQ)
      COMMON/REFDEP/IREFD(MFREQ)
      COMMON/CENTRL/ZND,IFZ0
      PARAMETER (UN=1.D0, HALF=0.5D0)
      PARAMETER (THIRD=UN/3., QUART=UN/4., SIXTH=UN/6.D0)
      PARAMETER (TAUREF = 0.6666666666667)
      DATA AMU/.887298334620742D0,.5D0,.112701665379258D0/,
     1     WTMU/.277777777777778D0,.444444444444444D0,.277777777777778D0
     1         /
      DATA TYPION /' I  ',' II ',' III',' IV ',' V  ',
     *             ' VI ',' VII','VIII',' IX '/
C
      NMU=3
      ND1=ND-1
C
C     Overall loop over frequencies 
C
      DO 100 IJ=1,NFREQ
      TAUMIN=CH(IJ,1)/DENS(1)*DM(1)*HALF
      TAU(1)=TAUMIN
      IREF=1
      DO 10 I=1,ND1
         DT(I)=(DM(I+1)-DM(I))*(CH(IJ,I+1)/DENS(I+1)+CH(IJ,I)/DENS(I))*
     *         HALF
         ST0(I)=ET(IJ,I)/CH(IJ,I)
         SS0(I)=-SC(IJ,I)/CH(IJ,I)
         TAU(I+1)=TAU(I)+DT(I)
         IF(TAU(I).LE.TAUREF.AND.TAU(I+1).GT.TAUREF) IREF=I
   10 CONTINUE
      IREFD(IJ)=IREF
      ST0(ND)=ET(IJ,ND)/CH(IJ,ND)
      SS0(ND)=-SC(IJ,ND)/CH(IJ,ND)
      FR=FREQ(IJ)
      BNU=BN*(FR*1.E-15)**3
      PLAND=BNU/(EXP(HK*FR/TEMP(ND  ))-UN)
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)
      DPLAN=(PLAND-DPLAN)/DT(ND1)
C
C   +++++++++++++++++++++++++++++++++++++++++
C   FIRST PART  -  VARIABLE EDDINGTON FACTORS
C   +++++++++++++++++++++++++++++++++++++++++
C
      ALB1=0.
      DO 21 I=1,NMU
C
C   ************************
C   UPPER BOUNDARY CONDITION
C   ************************
C
         ID=1
         DTP1=DT(1)
         Q0=0.
         P0=0.
C
C       allowance for non-zero optical depth at the first depth point
C
         TAMM=TAUMIN/AMU(I)
         IF(TAMM.GT.0.01) THEN
            P0=UN-EXP(-TAMM)
          ELSE
            P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM)))
         END IF
         EX=UN-P0
         Q0=Q0+P0*AMU(I)*WTMU(I)
C
         DIV=DTP1/AMU(I)*THIRD
         VL(I)=DIV*(ST0(ID)+HALF*ST0(ID+1))+ST0(ID)*P0
         DO 20 J=1,NMU
            BB(I,J)=SS0(ID)*WTMU(J)*(DIV+P0)-ALB1*WTMU(J)
   20       CC(I,J)=-HALF*DIV*SS0(ID+1)*WTMU(J)
         BB(I,I)=BB(I,I)+AMU(I)/DTP1+UN+DIV
         CC(I,I)=CC(I,I)+AMU(I)/DTP1-HALF*DIV
         ANU(I,ID)=0.
   21 CONTINUE
C
C     Matrix inversion: instead of calling MATINV, a very fast inlined 
C     routine MINV3 for a specific 3 x 3 matrix inversion
C
C     CALL MATINV(BB,NMU,3)
C
C     ******************************
      BB(2,1)=BB(2,1)/BB(1,1)
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2)
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3)
      BB(3,1)=BB(3,1)/BB(1,1)
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2)
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3)
C
      BB(3,2)=-BB(3,2)
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1)
      BB(2,1)=-BB(2,1)
C
      BB(3,3)=UN/BB(3,3)
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2)
      BB(2,2)=UN/BB(2,2)
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1)
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1)
      BB(1,1)=UN/BB(1,1)
C
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1)
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2)
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1)
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2)
      BB(3,1)=BB(3,3)*BB(3,1)
      BB(3,2)=BB(3,3)*BB(3,2)
C     ******************************
C
      DO 22 I=1,NMU
         DO 22 J=1,NMU
            S=0.
            DO 23 K=1,NMU
   23          S=S+BB(I,K)*CC(K,J)
            D(I,J,ID)=S
            ANU(I,1)=ANU(I,1)+BB(I,J)*VL(J)
   22 CONTINUE
C
C   *******************
C   NORMAL DEPTH POINTS
C   *******************
C
      DO 34 ID=2,ND1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=HALF*(DTM1+DTP1)
         AL=UN/DTM1/DT0
         GA=UN/DTP1/DT0
         BE=AL+GA
         A=(UN-HALF*AL*DTP1*DTP1)*SIXTH
         C=(UN-HALF*GA*DTM1*DTM1)*SIXTH
         B=UN-A-C
         VL0=A*ST0(ID-1)+B*ST0(ID)+C*ST0(ID+1)
         DO 30 I=1,NMU
            DO 30 J=1,NMU
               AA(I,J)=-A*SS0(ID-1)*WTMU(J)
               CC(I,J)=-C*SS0(ID+1)*WTMU(J)
               BB(I,J)=B*SS0(ID)*WTMU(J)
   30    CONTINUE
         DO 31 I=1,NMU
            DIV=AMU(I)**2
            VL(I)=VL0
            AA(I,I)=AA(I,I)+DIV*AL-A
            CC(I,I)=CC(I,I)+DIV*GA-C
            BB(I,I)=BB(I,I)+DIV*BE+B
   31    CONTINUE
         DO 32 I=1,NMU
            S1=0.
            DO 36 J=1,NMU
               S=0.
               S1=S1+AA(I,J)*ANU(J,ID-1)
               DO 35 K=1,NMU
   35             S=S+AA(I,K)*D(K,J,ID-1)
   36          BB(I,J)=BB(I,J)-S
            VL(I)=VL(I)+S1
   32    CONTINUE
C
C     Matrix inversion: instead of calling MATINV, a very fast inlined 
C     routine MINV3 for a specific 3 x 3 matrix inversion
C
C     CALL MATINV(BB,NMU,3)
C
C     ******************************
      BB(2,1)=BB(2,1)/BB(1,1)
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2)
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3)
      BB(3,1)=BB(3,1)/BB(1,1)
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2)
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3)
C
      BB(3,2)=-BB(3,2)
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1)
      BB(2,1)=-BB(2,1)
C
      BB(3,3)=UN/BB(3,3)
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2)
      BB(2,2)=UN/BB(2,2)
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1)
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1)
      BB(1,1)=UN/BB(1,1)
C
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1)
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2)
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1)
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2)
      BB(3,1)=BB(3,3)*BB(3,1)
      BB(3,2)=BB(3,3)*BB(3,2)
C     ******************************
C
         DO 33 I=1,NMU
            ANU(I,ID)=0.
            DO 33 J=1,NMU
               S=0.
               DO 37 K=1,NMU
   37             S=S+BB(I,K)*CC(K,J)
               D(I,J,ID)=S
               ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
   33    CONTINUE
   34 CONTINUE
C
C   ************
C   LOWER BOUNDARY CONDITION
C   ************
C
      ID=ND
C
C     First option:
C     b.c. is different from stellar atmospheres; expresses symmetry
C     at the central plane   I(taumax,-mu,nu)=I(taumax,+mu,nu)
C
      IF(IFZ0.EQ.0) THEN
         B=DTP1*HALF
         A=0.
         DO 190 I=1,NMU
            BI=B/AMU(I)
            AI=A/AMU(I)
            VL(I)=ST0(ID)*BI+ST0(ID-1)*AI
            DO 180 J=1,NMU
               AA(I,J)=-AI*SS0(ID-1)*WTMU(J)
  180          BB(I,J)=BI*SS0(ID)*WTMU(J)
            AA(I,I)=AA(I,I)+AMU(I)/DTP1-AI
            BB(I,I)=BB(I,I)+AMU(I)/DTP1+BI
  190    CONTINUE
         DO 220 I=1,NMU
            S1=0.
            DO 210 J=1,NMU
               S=0.
               S1=S1+AA(I,J)*ANU(J,ID-1)
               DO 200 K=1,NMU
  200             S=S+AA(I,K)*D(K,J,ID-1)
  210          BB(I,J)=BB(I,J)-S
            VL(I)=VL(I)+S1
  220    CONTINUE
C
C     Second option:
C     b.c. is the same as in stellar atmospheres - the last depth point
C     is not at the central plane
C
      ELSE
      DO 41 I=1,NMU
         AA(I,I)=AMU(I)/DTP1
         VL(I)=PLAND+AMU(I)*DPLAN+AA(I,I)*ANU(I,ID-1)
         DO 40 J=1,NMU
   40       BB(I,J)=-AA(I,I)*D(I,J,ID-1)
         BB(I,I)=BB(I,I)+AA(I,I)+UN
   41 CONTINUE
      END IF
C
C     Matrix inversion: instead of calling MATINV, a very fast inlined 
C     routine MINV3 for a specific 3 x 3 matrix inversion
C
C     CALL MATINV(BB,NMU,3)
C
C     ******************************
      BB(2,1)=BB(2,1)/BB(1,1)
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2)
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3)
      BB(3,1)=BB(3,1)/BB(1,1)
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2)
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3)
C
      BB(3,2)=-BB(3,2)
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1)
      BB(2,1)=-BB(2,1)
C
      BB(3,3)=UN/BB(3,3)
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2)
      BB(2,2)=UN/BB(2,2)
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1)
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1)
      BB(1,1)=UN/BB(1,1)
C
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1)
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2)
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1)
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2)
      BB(3,1)=BB(3,3)*BB(3,1)
      BB(3,2)=BB(3,3)*BB(3,2)
C     ******************************
C
      DO 42 I=1,NMU
         ANU(I,ID)=0.
         DO 42 J=1,NMU
            D(I,J,ID)=0.
            ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
   42 CONTINUE
C
C   ************
C   BACKSOLUTION
C   ************
C
      ID=ND
      FKK(ND)=THIRD
      AJ=0.
      AK=0.
      DO I=1,NMU
         RMU=WTMU(I)*ANU(I,ID)
         AJ=AJ+RMU
         AK=AK+RMU*AMU(I)*AMU(I)
      END DO
      RDD(ID)=AJ
c     IF(IBC.NE.0) FKK(ND)=AK/AJ
      FKK(ND)=AK/AJ
             if(fr.gt.1.34195e13.and.fr.lt.1.34316e13)
     *   write(74,674)
     *   id,ij,wlam(ij),st0(id),ss0(id),tau(id),fkk(id),aj
      DO 54 ID=ND-1,1,-1
         DO 50 I=1,NMU
            DO 50 J=1,NMU
               ANU(I,ID)=ANU(I,ID)+D(I,J,ID)*ANU(J,ID+1)
   50    CONTINUE
         AJ=0.
         AK=0.
         DO 52 I=1,NMU
            DIV=WTMU(I)*ANU(I,ID)
            AJ=AJ+DIV
            AK=AK+DIV*AMU(I)**2
   52    CONTINUE
         FKK(ID)=AK/AJ
             if(fr.gt.1.34195e13.and.fr.lt.1.34316e13)
     *   write(74,674)
     *   id,ij,wlam(ij),ch(ij,id),et(ij,id),
     *   st0(id),ss0(id),tau(id),fkk(id),aj
  674 format(2i3,f10.1,1p7e12.4)
   54 CONTINUE
C
C     surface Eddington actor
C
      AH=0.
      DO 53 I=1,NMU
   53    AH=AH+WTMU(I)*AMU(I)*ANU(I,1)
      FH=AH/AJ-HALF*ALB1
C
c     FKK(ND)=THIRD
C
C
C   +++++++++++++++++++++++++++++++++++++++++
C   SECOND PART  -  DETERMINATION OF THE MEAN INTENSITIES
C   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN EDDINGTON FACTORS
C   +++++++++++++++++++++++++++++++++++++++++
C
      DTP1=DT(1)
      DIV=DTP1*THIRD
      BBB=FKK(1)/DTP1+FH+DIV+SS0(1)*(DIV+Q0)
      CCC=FKK(2)/DTP1-HALF*DIV*(UN+SS0(2))
      VLL=DIV*(ST0(1)+HALF*ST0(2))+ST0(1)*Q0
      AANU(1)=VLL/BBB
      DDD(1)=CCC/BBB
      DO 60 ID=2,ND1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=HALF*(DTP1+DTM1)
         AL=UN/DTM1/DT0
         GA=UN/DTP1/DT0
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH
         AAA=AL*FKK(ID-1)-A*(UN+SS0(ID-1))
         CCC=GA*FKK(ID+1)-C*(UN+SS0(ID+1))
         BBB=(AL+GA)*FKK(ID)+(UN-A-C)*(UN+SS0(ID))
         VLL=A*ST0(ID-1)+C*ST0(ID+1)+(UN-A-C)*ST0(ID)
         BBB=BBB-AAA*DDD(ID-1)
         DDD(ID)=CCC/BBB
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB
   60 CONTINUE
C
C     Lower boundary condition
C     1.option -  different from stellar atmospheres
C
      IF(IFZ0.EQ.0) THEN
         B=DTP1*HALF
         BBB=FKK(ND)/DTP1+B*(UN+SS0(ND))
         AAA=FKK(ND-1)/DTP1
         VLL=B*ST0(ND)
       ELSE
C
C     Lower boundary condition
C     2.option - stellar atmospheric
C
         BBB=FKK(ND)/DTP1+HALF
         AAA=FKK(ND1)/DTP1
         VLL=HALF*PLAND+DPLAN*THIRD
      END IF
      BBB=BBB-AAA*DDD(ND1)
      RDD(ND)=(VLL+AAA*AANU(ND1))/BBB
      DO 70 IID=1,ND1
         ID=ND-IID
         RDD(ID)=AANU(ID)+DDD(ID)*RDD(ID+1)
   70 CONTINUE
      FLUX(IJ)=FH*RDD(1)
C
C     if needed (if iprin.ge.3), output of interesting physical
C     quantities at the monochromatic optical depth  tau(nu)=2/3
C
      IF(IPRIN.ge.3) THEN
      T0=LOG(TAU(IREF+1)/TAU(IREF))
      X0=LOG(TAU(IREF+1)/TAUREF)/T0
      X1=LOG(TAUREF/TAU(IREF))/T0
      DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1)
      TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1)
      STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1)
      SCREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1)
      SSREF=EXP(LOG(-SS0(IREF)*RDD(IREF))*X0+
     *           LOG(-SS0(IREF+1)*RDD(IREF+1))*X1)
      SREF=STREF+SSREF
      ALM=2.997925E18/FREQ(IJ)
      WRITE(96,636) IJ,ALM,IREF,DMREF,TREF,SCREF,STREF,SSREF,SREF
  636 FORMAT(1H ,I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3,E11.3)
      END IF
C
C   THIRD PART  -  DETERMINATION OF THE SPECIFIC INTENSITIES
C   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
C
      if(iflux.eq.0) go to 100
      DO 85 IMU=1,NMU0
      ANX=ANGL(IMU)
      DTP1=DT(1)
      DIV=DTP1*THIRD/ANX
C
      TAMM=TAUMIN/ANX
      IF(TAMM.LT.0.01) THEN
         P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM)))
       ELSE
         P0=UN-EXP(-TAMM)
      END IF
C
      BBB=ANX/DTP1+UN+DIV
      CCC=ANX/DTP1-HALF*DIV
      VLL=(DIV+P0)*(ST0(1)-SS0(1)*RDD(1))
     *    +HALF*DIV*(ST0(2)-SS0(2)*RDD(2))
      AANU(1)=VLL/BBB
      DDD(1)=CCC/BBB
      DIV=ANX*ANX
      DO 82 ID=2,ND1
         DTM1=DT(ID-1)
         DTP1=DT(ID)
         DT0=HALF*(DTP1+DTM1)
         AL=UN/DTM1/DT0
         GA=UN/DTP1/DT0
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH
         AAA=DIV*AL-A
         CCC=DIV*GA-C
         BBB=DIV*(AL+GA)+UN-A-C
         VLL=A*(ST0(ID-1)-SS0(ID-1)*RDD(ID-1))+
     *       C*(ST0(ID+1)-SS0(ID+1)*RDD(ID+1))+
     *       (UN-A-C)*(ST0(ID)-SS0(ID)*RDD(ID))
         BBB=BBB-AAA*DDD(ID-1)
         DDD(ID)=CCC/BBB
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB
   82 CONTINUE
C
C     Lower boundary condition
C     1.option -  different from stellar atmospheres
C
      IF(IFZ0.EQ.0) THEN
         B=DTP1*HALF/ANX
         BBB=ANX/DTP1+B*(UN+SS0(ND))
         AAA=ANX/DTP1
         VLL=B*ST0(ND)
       ELSE
C
C     Lower boundary condition
C     2.option - stellar atmospheric
C
         AAA=ANX/DTP1
         BBB=AAA+UN
         VLL=PLAND+ANX*DPLAN
      END IF
C
      RINT(ND,IMU)=(VLL+AAA*AANU(ND1))/(BBB-AAA*DDD(ND1))
      DO 84 IID=1,ND1
         ID=ND-IID
         RINT(ID,IMU)=AANU(ID)+DDD(ID)*RINT(ID+1,IMU)
   84 CONTINUE
   85 CONTINUE
c
      FLX=0.
      DO 86 IMU=1,NMU0
         RINT(1,IMU)=RINT(1,IMU)/HALF
         FLX=FLX+ANGL(IMU)*WANGL(IMU)*RINT(1,IMU)
   86 CONTINUE
      FLX=FLX*HALF
c     FLUX(IJ)=FLX
C
C     output of emergent specific intensities to Unit 10 
C     and 18 (continuum)
C
      IF(IJ.GT.2) THEN
         WRITE(10,641) WLAM(IJ),FLX,(RINT(1,IMU),IMU=1,NMU0)
       ELSE
         WRITE(18,641) WLAM(IJ),FLX,(RINT(1,IMU),IMU=1,NMU0)
      END IF
c
      if(iprin.eq.4) then
c
c     compute contribution function C_i (ctri) and C_r (ctrr)
c     following Magain (1986, A&A 163, 135)
c
       if(ijctr(ij).gt.0) then
         xfr0=(freq(ij)-freq(2))/(freq(1)-freq(2))
         tauc=ch(1,1)/dens(1)*dm(1)*half
         do 150 id=1,nd
           chc1=ch(1,id)
           chc2=ch(2,id)
           chcc=chc2+xfr0*(chc1-chc2)
           etc1=et(1,id)
           etc2=et(2,id)
           etcc=etc2+xfr0*(etc1-etc2)
           stcc=etcc/chcc
           cint=cint2(id)+xfr0*(cint1(id)-cint2(id))
           avx=(chc1+chc2)*0.5*relop
           call linop(id,abxli,emxli,avx)
           sli0=emxli(ij)/abxli(ij)
           abt0=ch(ij,id)
           emt0=et(ij,id)
           stt0=emt0/abt0
           Xkar(id)=abxli(ij)+chcc*stcc/cint
           ctri(id)=tauc*abt0/chc1*stt0*exp(-tau(id))
           if(tau(id).gt.70.) ctri(id)=0.
           ctrr(id)=tauc/chc1*abxli(ij)*(un-sli0/cint)
           if(id.lt.nd) then
             dtc=(ch(1,id+1)/dens(id+1)+ch(1,id)/dens(id))
             tauc=tauc+half*dtc*(dm(id+1)-dm(id))
           endif
  150    continue
         taurs=Xkar(1)/dens(1)*dm(1)*half
         xcti=ctri(1)*half*(dm(2)-dm(1))
         xctr=ctrr(1)*half*(dm(2)-dm(1))
         do 155 i=1,nd-1
           ctrr(i)=ctrr(i)*exp(-taurs)
           if(i.eq.1) xctr=xctr*exp(-taurs)
           if(i.gt.1) then
             xcti=xcti+ctri(i)*half*(dm(i+1)-dm(i-1))
             xctr=xctr+ctrr(i)*half*(dm(i+1)-dm(i-1))
           endif
           if(taurs.gt.70.) ctrr(i)=0.
           dtrs=(dm(i+1)-dm(i))*(Xkar(i+1)/dens(i+1)+Xkar(i)/dens(i))
           taurs=taurs+half*dtrs
  155    continue
         ctrr(nd)=0.
         alam=2.997925d18/freq(ij)
         il0=ijctr(ij)
         il=indlin(il0)
         iat=indat(il)/100
         ion=mod(indat(il),100)
         write(97,376) il,alam,typat(iat),typion(ion),iref,dmref,tref
  376    format(i5,f11.4,2x,2a4,i8,1pe12.4,0pf10.1)
         do 156 id=1,nd
           ctrip=ctri(id)/xcti
           ctrrp=ctrr(id)/xctr
           write(97,377) id,dm(id),tau(id),ctrip,ctrrp
  377      format(i4,1p4e12.4)
  156    continue
       else if(ij.eq.1) then
         do 161 id=1,nd
           cint1(id)=rint(id,nmu0)
  161    continue
       else if(ij.eq.2) then
         do 162 id=1,nd
           cint2(id)=rint(id,nmu0)
  162    continue
       endif
      endif
  100 CONTINUE
  641 FORMAT(1H ,f10.3,1pe15.5/(1P5E15.5))
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE OUTPRI
C
C     Output of synthetic spectrum
C
C     Output onto unit 7 serves as an input to the next program
C     ROTINS, which performs convolutions for the rotational and
C     instrumental broadening, and plots the synthetic spectrum
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (UN=1.,CAS=1./2.997925D18,EQWC=1.19917D22)
      PARAMETER (PI2=3.141592654/2.)
      DIMENSION FLX(3),REL(3),ALX(3)
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
C
      if(ifwin.le.0) then
C
C     output of synthetic spectrum on unit 7
C
      DO 10 IJ=3,NFREQ-1
         FLAM=FLUX(IJ)*FREQ(IJ)*FREQ(IJ)*CAS
         WRITE(7,701) WLAM(IJ),FLAM
   10 CONTINUE
C
C     output of the continuum flux on unit 17
C
      FLAM=FLUX(1)*FREQ(1)*FREQ(1)*CAS
      WRITE(17,701) WLAM(1),FLAM
      IF(IBLANK.EQ.NBLANK) THEN
         FLAM=FLUX(NFREQ)*FREQ(NFREQ)*FREQ(NFREQ)*CAS
         WRITE(7,701) WLAM(NFREQ),FLAM
         FLAM=FLUX(2)*FREQ(2)*FREQ(2)*CAS
         WRITE(17,701) WLAM(2),FLAM
      END IF
      else
        DO IJ=1,NFROBS
           FLAM=FLUX(IJ)*FRQOBS(IJ)*FRQOBS(IJ)*CAS*0.5
           flam=max(flam,1.e-40)
           WRITE(7,701) WLobs(IJ),FLAM
        enddo
      end if
C
C     unit 6 and 16 outputs
C
      if(iprin.le.-1) return
      if(iprin.ge.1) then
      WRITE(6,600)
      WRITE(6,601)
      end if
      K1=0
      EQW=0.
      EQWP=0.
      IF(IBLANK.EQ.1) EQWT=0.
      IF(IBLANK.EQ.1) EQWTP=0.
      XX=UN/(FREQ(2)-FREQ(1))
      XXX=UN/(FREQ(1)+FREQ(2))/(FREQ(1)+FREQ(2))
      if(ifwin.le.0) then
      DO 20 IJ=1,NFREQ
         FLAM=FLUX(IJ)*FREQ(IJ)*FREQ(IJ)*CAS
         CONT=((FREQ(IJ)-FREQ(1))*FLUX(2)+(FREQ(2)-FREQ(IJ))*FLUX(1))*XX
         RE0=FLUX(IJ)/CONT
         EQW=EQW+(UN-RE0)*W(IJ)
         REP=RE0
         IF(REP.GT.UN) REP=UN
         EQWP=EQWP+(UN-REP)*W(IJ)
         if(iprin.le.0) go to 20
         K1=K1+1
         FLX(K1)=LOG10(FLAM)
         ALX(K1)=WLAM(IJ)
         REL(K1)=RE0
         IF(K1.EQ.3.OR.IJ.EQ.NFREQ) THEN
            WRITE(6,602) (ALX(I),FLX(I),REL(I),I=1,K1)
            K1=0
         END IF
   20 CONTINUE
      else
      DO 21 IJ=1,NFROBS
         FLAM=FLUX(IJ)*FREQ(IJ)*FREQ(IJ)*CAS
         CONT=((FRQOBS(IJ)-FREQ(1))*FLUX(2)+
     *        (FREQ(2)-FRQOBS(IJ))*FLUX(1))*XX
         RE0=FLUX(IJ)/CONT
         EQW=EQW+(UN-RE0)*W(IJ)
         REP=RE0
         IF(REP.GT.UN) REP=UN
         EQWP=EQWP+(UN-REP)*W(IJ)
         if(iprin.le.0) go to 21
         K1=K1+1
         FLX(K1)=LOG10(FLAM)
         ALX(K1)=WLAM(IJ)
         REL(K1)=RE0
         IF(K1.EQ.3.OR.IJ.EQ.NFREQ) THEN
            WRITE(6,602) (ALX(I),FLX(I),REL(I),I=1,K1)
            K1=0
         END IF
   21 CONTINUE
      end if
C
C     output of partial equivalent widths on unit 16
C
      EQW=EQW*EQWC*XXX
      EQWT=EQWT+EQW
      EQWP=EQWP*EQWC*XXX
      EQWTP=EQWTP+EQWP
      if(iprin.gt.0) WRITE(6,603) EQW,EQWP,EQWT,EQWTP
      WRITE(16,616) WLAM(1),WLAM(2),EQW,EQWP,EQWT,EQWTP
C
  600 FORMAT(/' EMERGENT RADIATION'/' ------------------'/)
  601 FORMAT(3('   LAMBDA  LOG HLAM    REL')/)
  602 FORMAT(3(2X,F9.3,F8.4,F7.3))
  603 FORMAT(/,'  EQUIVALENT WIDTH THIS SET  =',2F8.1,' mA'/
     *            '  EQUIVALENT WIDTH TOTAL     =',2F8.1,' mA'//)
  616 FORMAT(2F12.3,4F12.1)
  701 FORMAT(F12.5,1PE15.5)
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE CROSET(CROSS)
C
C     SET UP ARRAY CROSS  - PHOTOIONIZATION CROSS-SECTIONS
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      DIMENSION CROSS(MCROSS,MFRQ)
      common/dissol/fropc(mlevel),indexp(mlevel)
C
      IJ0=2
      IF(NFREQ.EQ.1) IJ0=1
      IF(IMODE.EQ.2) IJ0=NFREQ
      DO IJ=1,IJ0
         DO IT=1,MCROSS
            CROSS(IT,IJ)=0.
         END DO
      END DO
      DO IT=1,NLEVEL
         IF(INDEXP(IT).NE.5) THEN
            DO IJ=1,IJ0
               FR=FREQ(IJ)
               CROSS(IT,IJ)=SIGK(FR,IT,0)
            END DO
         ELSE
            DO IJ=1,IJ0
               FR=FREQ(IJ)
               CROSS(IT,IJ)=SIGK(FR,IT,1)
               IF(FR.LT.FROPC(IT)) CROSS(IT,IJ)=0.
            END DO
         END IF
      END DO
C
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE CROSEW(CROSS)
C
C     SET UP COMMON/PHOPAR/  - PHOTOIONIZATION CROSS-SECTIONS
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      DIMENSION CROSS(MCROSS,MFRQ)
      common/dissol/fropc(mlevel),indexp(mlevel)
C
      IJ0=NFREQC
      DO IJ=1,IJ0
         DO IT=1,MCROSS
            CROSS(IT,IJ)=0.
         END DO
      END DO
      DO IT=1,NLEVEL
         IF(INDEXP(IT).NE.5) THEN
            DO IJ=1,IJ0
               FR=FREQC(IJ)
               CROSS(IT,IJ)=SIGK(FR,IT,0)
            END DO
         ELSE
            DO IJ=1,IJ0
               FR=FREQC(IJ)
               CROSS(IT,IJ)=SIGK(FR,IT,1)
               IF(FR.LT.FROPC(IT)) CROSS(IT,IJ)=0.
            END DO
         END IF
      END DO
C
      RETURN
      END
C
C ********************************************************************
C
C
 
      FUNCTION SIGK(FR,ITR,MODE)
C     ==========================
C
C     driver for evaluating the photoionization cross-sections
C
C     Input: FR  -  frequency
C            ITR -  index of the transition
c            mode - =0 cross-section equal to zero longward of edge
c            mode - >0 cross-section non-zero (extrapolated) longward of edge
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (SIH0=2.815D29, E10=2.3025851)
      parameter (wi1=911.753878, wi2=227.837832, un=1.e0)
      CHARACTER*10 TYPLEV(MLEVEL)
      COMMON/PRINTP/TYPLEV
      COMMON/TOPCS/CTOP(MFIT,MCROSS), ! sigma = alog10(sigma/10^-18) of fit point
     +             XTOP(MFIT,MCROSS)  ! x = alog10(nu/nu0) of fit point 
      common/dissol/fropc(mlevel),indexp(mlevel)
      DIMENSION XFIT(MFIT) , ! local array containing x     for OP data
     +          SFIT(MFIT)   ! local array containing sigma for OP data
C
      PEACH(X,S,A,B)  =A*X**S*(B+X*(1.-B))*1.E-18
      HENRY(X,S,A,B,C)=A*X**S*(C+X*(B-2.*C+X*(1.+C-B)))*1.E-18
C
      SIGK=0.
      II=ITR
      FR0=ENION(II)/6.6256E-27
      IF(FR0.LE.0.) RETURN
      wl0=2.997925e18/fr0
C
C     wavelength with an explicit correction to the air wavalength
C
      IF(WL0.GT.vaclim) THEN
         ALM=1.E8/(WL0*WL0)
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0/(XN1*1.D-6+UN)
         fr0=2.997925e18/wl0
      END IF 
c
      IF(mode.eq.0  .and. FR.LT.FR0) RETURN
C
C     IBF(ITR) is the switch controlling the mode of evaluation of the
C        cross-section:
C      = 0  hydrogenic cross-section, with Gaunt factor set to 1
C      = 1  hydrogenic cross-section with exact Gaunt factor
C      = 2  Peach-type expression (see function PEACH)
C      = 3  Henry-type expression (see function HENRY)
C      = 4  Butler new calculations
C      = 7  hydrogenic cross-section with Gaunt factor from K. Werner
C      = 9  Opacity project fits (routine TOPBAS - interpolations)
C      > 100 - cross-sections extracted form TOPBASE, for several points
C           In this case, IBF-100 is the number of points
C      < 0  non-standard, user supplied expression (user should update
C           subroutine SPSIGK)
C
C      for H- : for any IBF > 0  - standard expression
C      for He I:
C       for IBF = 11 or = 13  -  Opacity Project cross section
C                Seaton-Ferney's cubic fits, Hummer's procedure (HEPHOT)
C           IBF = 11  means that the multiplicity S=1 (singlet)
C           IBF = 13  means that the multiplicity S=3 (triplet)
C       for IBF = 10  - cross section, based on Opacity Project, but
C                       appropriately averaged for an averaged level
C
C
      IB=IBF(ITR)
      IQ=NQUANT(II)
      IE=IEL(II)
      IF(IB.LT.0) GO TO 60
      IF(IE.EQ.IELHM) GO TO 40
      IF(IB.EQ.0) RETURN
      IF(IE.EQ.IELHE1.AND.IB.GE.10.AND.IB.LE.13) GO TO 50
      CH=IZ(IE)*IZ(IE)
      IQ5=IQ*IQ*IQ*IQ*IQ
C
      IF(IB.EQ.0) THEN
C
C        hydrogenic expression (for IBF = 0)
C
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5
C
C        exact hydrogenic - with Gaunt factor (for IBF=1)
C
       ELSE IF(IB.EQ.1) THEN
         SIGK=SIH0/FR/FR/FR*CH*CH/IQ5
c        IF(FR.GE.FR0.OR.(IE.EQ.IELH.AND.IQ.LE.3)) 
c    *   SIGK=SIGK*GAUNT(IQ,FR/CH)
         fr0l=0.95*fr0
         if(fr.ge.fr0) then
            sigk=sigk*gaunt(iq,fr/ch)
          else if(fr.ge.fr0l) then
            gau0=gaunt(iq,fr0/ch)
            corg=(fr-fr0l)/(fr0-fr0l)*(gau0-1.)+1.
            sigk=sigk*corg
         end if
       ELSE IF(IB.EQ.2) THEN
C
C        Peach-type formula (for IBF=2)
C
         IF(GAMBF(II).GT.0) THEN
            IF(GAMBF(II).LT.1.E6) THEN
              FR0=2.997925E18/GAMBF(II)
            ELSE
              FR0=GAMBF(II)
            END IF
            IF(FR.LT.FR0) RETURN        
         END IF 
         FREL=FR0/FR
         SIGK=PEACH(FREL,S0BF(II),ALFBF(II),BETBF(II))
       ELSE IF(IB.EQ.3) THEN
C
C        Henry-type formula (for IBF=3)
C
         FREL=FR0/FR
         SIGK=HENRY(FREL,S0BF(II),ALFBF(II),BETBF(II),GAMBF(II))
C
C     Butler expression
C
       ELSE IF(IB.EQ.4) THEN
         FREL=FR0/FR
         XL=LOG(FREL)
         SL=S0BF(II)+XL*(ALFBF(II)+XL*BETBF(II))
         SIGK=EXP(SL)
C
C     exact hydrogenic - with Gaunt factor from K Werner (for IBF=7)
C
       ELSE IF(IB.EQ.7) THEN
         IQ5=IQ*IQ*IQ*IQ*IQ
         SIGK=SIH0/(FR*FR*FR)*CH*CH/IQ5*GNTK(IQ,FR/CH)
C
C     selected Opacity Project data (for IBF=9)
C     (c.-s. evaluated by routine TOPBAS which needs an input file RBF.DAT)
C
       ELSE IF(IB.EQ.9) THEN
         SIGK=TOPBAS(FR,FR0,TYPLEV(II))
C
C     other Opacity Project data (for IBF>100)
C     (c.-s. evaluated by interpolating from direct input data)
C
       ELSE IF(IB.GT.100) THEN
         NFIT=IB-100
         X = LOG10(FR/FR0)
         IF(X.LT.XTOP(1,II)) THEN
            SIGM=0.
          ELSE
            DO IFIT = 1,NFIT
               XFIT(IFIT) = XTOP(IFIT,II)
               SFIT(IFIT) = CTOP(IFIT,II)
            END DO
            SIGM  = YLINTP (X,XFIT,SFIT,NFIT,MFIT)
            SIGM  = 1.D-18*EXP(E10*SIGM)
         END IF
         SIGK=SIGM
      END IF
      if(iatm(ii).eq.iath.and.ii.gt.n0hn+2.
     *   and.ib.le.1.and.fr.lt.fr0) then
         fr1=fropc(ii)
         frdec=min(fr1*1.25,fr0)
         if(fr.gt.fr1.and.fr.lt.frdec)
     *      sigk=sigk*(fr-fr1)/(frdec-fr1)
      end if
      RETURN
C
C     special expression for H-
C
   40 SIGK=SBFHMI(FR)
      RETURN
C
C     He I cross-sections
C
   50 SIGK=SBFHE1(II,IB,FR)
      RETURN
C
C     non-standard, user supplied form of cross-section (for IBF < 0)
C
   60 CALL SPSIGK(ITR,IB,FR,SIGSP)
      SIGK=SIGSP
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION GAUNT(I,FR)
C     ====================
C
C     Hydrogenic bound-free Gaunt factor for the principal quantum
C     number I and frequency FR
C
      INCLUDE 'PARAMS.FOR'
      X=FR/2.99793E14
      GAUNT=1.
      IF(I.EQ.1) THEN
      GAUNT=1.2302628+X*(-2.9094219E-3+X*(7.3993579E-6-8.7356966E-9*X))
     *+(12.803223/X-5.5759888)/X
       ELSE IF(I.EQ.2) THEN
      GAUNT=1.1595421+X*(-2.0735860E-3+2.7033384E-6*X)+(-1.2709045+
     *(-2.0244141/X+2.1325684)/X)/X
       ELSE IF(I.EQ.3) THEN
      GAUNT=1.1450949+X*(-1.9366592E-3+2.3572356E-6*X)+(-0.55936432+
     *(-0.23387146/X+0.52471924)/X)/X
       ELSE IF(I.EQ.4) THEN
      GAUNT=1.1306695+X*(-1.3482273E-3+X*(-4.6949424E-6+2.3548636E-8*X))
     *+(-0.31190730+(0.19683564-5.4418565E-2/X)/X)/X
       ELSE IF(I.EQ.5) THEN
      GAUNT=1.1190904+X*(-1.0401085E-3+X*(-6.9943488E-6+2.8496742E-8*X))
     *+(-0.16051018+(5.5545091E-2-8.9182854E-3/X)/X)/X
       ELSE IF(I.EQ.6) THEN
      GAUNT=1.1168376+X*(-8.9466573E-4+X*(-8.8393133E-6+3.4696768E-8*X))
     *+(-0.13075417+(4.1921183E-2-5.5303574E-3/X)/X)/X
       ELSE IF(I.EQ.7) THEN
      GAUNT=1.1128632+X*(-7.4833260E-4+X*(-1.0244504E-5+3.8595771E-8*X))
     *+(-9.5441161E-2+(2.3350812E-2-2.2752881E-3/X)/X)/X
       ELSE IF(I.EQ.8) THEN
      GAUNT=1.1093137+X*(-6.2619148E-4+X*(-1.1342068E-5+4.1477731E-8*X))
     *+(-7.1010560E-2+(1.3298411E-2 -9.7200274E-4/X)/X)/X
       ELSE IF(I.EQ.9) THEN
      GAUNT=1.1078717+X*(-5.4837392E-4+X*(-1.2157943E-5+4.3796716E-8*X))
     *+(-5.6046560E-2+(8.5139736E-3-4.9576163E-4/X)/X)/X
       ELSE IF(I.EQ.10) THEN
      GAUNT=1.1052734+X*(-4.4341570E-4+X*(-1.3235905E-5+4.7003140E-8*X))
     *+(-4.7326370E-2+(6.1516856E-3-2.9467046E-4/X)/X)/X
      END IF
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION GNTK(I,FR)
C     ===================
C
C     Hydrogenic bound-free Gaunt factor for the principal quantum
C     number I and frequency FR (from Klaus Werner)
C
      INCLUDE 'PARAMS.FOR'
      GNTK=1.
      IF(I.GT.3) GO TO 16
      Y=1./FR      
      GO TO (1,2,3),I
    1 GNTK=0.9916+Y*(2.71852D13-Y*2.26846D30)
      GO TO 16
    2 GNTK=1.1050-Y*(2.37490D14-Y*4.07677D28)
      GO TO 16
    3 GNTK=1.1010-Y*(0.98632D14-Y*1.03540D28)
   16 RETURN
      END
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE SPSIGK(ITR,IB,FR,SIGSP)
C     ==================================
C
C     Non-standard evaluation of the photoionization cross-sections
C     Basically user-suppled procedure; here are some examples
C
      INCLUDE 'PARAMS.FOR'
      SIGSP=0.
      if(itr.le.0) return
C
C     Special formula for the He I ground state
C
      IF(IB.EQ.-201) SIGSP=7.3E-18*EXP(1.373-2.311E-16*FR)
C
C     Special formula for the averaged <n=2> level of He I
C
      IF(IB.EQ.-202) SIGSP=SGHE12(FR)
C
C     Carbon ground configuration levels 2p2 1D and 1S
C
      IF(IB.EQ.-602.OR.IB.EQ.-603) THEN
         CALL CARBON(IB,FR,SG)
         SIGSP=SG
      END IF
C
C     Hidalgo (Ap.J. 153, 981, 1968) photoionization data
C
      IF(IB.LE.-101.AND.IB.GE.-137) SIGSP=HIDALG(IB,FR)
C
C     Reilman and Manson (Ap.J. Suppl. 40, 815, 1979) photoionization data
C
      IF(IB.LE.-301.AND.IB.GE.-337) SIGSP=REIMAN(IB,FR)
      RETURN
      END
C
C
C
C     ****************************************************************
C
C
 
      SUBROUTINE CARBON(IB,FR,SG)
C     ===========================
C
C     Photoionization cross-section for neutral carbon 2p1D and 2p1S
C     levels (G.B.Taylor - private communication)
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION FR2(34),SG2(34),FR3(45),SG3(45)
      DATA FR2/ 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82,
     *                0.83,       0.85, 0.86, 0.87, 0.88, 0.89, 0.90,
     *          0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
     *          1.00, 1.10, 1.20, 1.30, 1.45, 1.50, 1.60, 1.80, 2./
      DATA SG2/ 12.04, 12.03, 12.09, 12.26, 12.60, 13.24, 14.36, 16.24,
     *          19.28, 23.94, 37.41, 42.88, 44.76, 43.41, 40.46, 37.19,
     *          34.26, 31.82, 29.96, 28.57, 27.68, 27.37, 27.84, 29.69,
     *          34.45, 46.35, 13.80, 11.54, 10.40,  8.96,  8.54,  7.47,
     *           6.53,  5.66/
      DATA FR3/ 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82,
     *          0.84, 0.86, 0.864,0.866,0.868,0.87, 0.874,0.876,0.88,
     *          0.882,0.884,0.886,0.888,0.89 ,0.894,0.896,0.898,0.90,
     *          0.904,0.908,0.910,0.920,0.94, 0.98, 1.00, 1.10, 1.20,
     *          1.26, 1.34, 1.36, 1.40, 1.46, 1.60, 1.70, 1.80, 2./
      DATA SG3/ 13.94, 13.29, 12.56, 11.73, 10.82, 10.18,  8.62,  7.27,
     *           5.74,  4.14,  4.61,  5.92,  6.94,  8.34, 10.21, 16.12,
     *          20.64, 34.56, 44.82, 57.71, 73.09, 89.99,106.38,127.08,
     *         128.38,124.44,117.17, 99.32, 82.95, 76.05, 52.65, 33.23,
     *          21.29, 18.69, 12.62, 11.44,  9.77,  7.53, 10.47,  9.65,
     *          10.19,  7.28,  6.70,  6.11,  4.96/
      DATA NC2,NC3/34,45/
      DATA FR0/3.28805E15/
      F=FR/FR0
      IF(IB.NE.-602) GO TO 25
      J=2
      IF(F.LE.FR2(1)) GO TO 20
      DO 10 I=2,NC2
         J=I
         IF(F.GT.FR2(I-1).AND.F.LE.FR2(I)) GO TO 20
   10 CONTINUE
   20 SG=(F-FR2(J-1))/(FR2(J)-FR2(J-1))*(SG2(J)-SG2(J-1))+SG2(J-1)
      SG=SG*1.E-18
   25 IF(IB.NE.-603) GO TO 50
      J=2
      IF(F.LE.FR3(1)) GO TO 40
      DO 30 I=2,NC3
         J=I
         IF(F.GT.FR3(I-1).AND.F.LE.FR3(I)) GO TO 40
   30 CONTINUE
   40 SG=(F-FR3(J-1))/(FR3(J)-FR3(J-1))*(SG3(J)-SG3(J-1))+SG3(J-1)
      SG=SG*1.E-18
   50 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
 
      FUNCTION SGHE12(FR)
C     ===================
C
C     Special formula for the photoionization cross-section from the
C     averaged <n=2> level of He I
C
      INCLUDE 'PARAMS.FOR'
      DATA C1/3.E0/,C2/9.E0/,C3/1.6E1/,
     * A1/6.45105E-18/,A2/3.02E-19/,A3/9.9847E-18/,A4/1.1763673E-17/,
     * A5/3.63662E-19/,A6/-2.783E2/,A7/1.488E1/,A8/-2.311E-1/,
     * E1/3.5E0/,E2/3.6E0/,E3/1.91E0/,E4/2.9E0/,E5/3.3E0/
      X=FR*1.E-15
      XX=LOG(FR)
      SGHE12=(C1*(A1/X**E1+A2/X**E2)+A3/X**E3+C2*(A4/X**E4+A5/X**E5)+
     *       C1*EXP(A6+XX*(A7+XX*A8)))/C3
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION HIDALG(IB,FR)
C     ======================
C
C     Read table of wavelengths and photo-ionization cross-sections
C     from Hidalgo (1968, Ap. J., 153, 981) for the species indicated by IB
C     (Hidalgo's number = INDEX = -IB-100).
C     Compute linearly interpolated value of the cross-section
C     at the frequency FR.
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION WL1(20),WL2(20),WLI(20),SIG0(20,24),SIGS(20)
C
      DATA WL1 /
     *  39.1, 80.9, 97.6,100.1,104.3,107.2,108.7,111.9,113.6,115.4,
     * 117.1,119.0,124.8,126.9,129.1,131.3,133.6,136.0,138.5,141.1/
      DATA WL2 /
     *  68.5, 80.9,100.1,120.9,158.8,165.7,177.3,190.6,200.7,206.2,
     * 211.9,218.0,224.5,231.3,246.3,5*0./
      DATA SIG0 /
     *120*0.,
     *.0460,.2400,.3500,.3700,.4000,.4300,.4400,.4600,.4700,.4900,
     *.5000,.5200,.5700,.6200, 6*0.,
     * 80*0.,
     *.0092,.1000,.1900,.2100,.2300,.2500,.2600,.2900,.3000,.3200,
     *.3400,.3500,.4100,.4300,.4500,.4800,.5000,.5300,.5600,.5900,
     * 20*0.,
     *.3400,.4600,.6300,.7700,.9100,1.080, 14*0.,
     * 20*0.,
     *.0064,.1100,.2200,.4100,.9400,1.000,1.300,1.600, 12*0.,
     * 80*0.,
     *.0370,.0650,.1300,.2400,.5500,.6300,.7700,.9500,1.100,1.250,
     * 10*0.,
     * 40*0.,
     *.0220,.0390,.0800,.1500,.3500,.4000,.4900,.6200,.7200,.7800,
     *.8500,.9300,1.020,
     * 7*0./
C
      INDEX=-IB-100
      NUM=20
      IF(INDEX.GE.13.AND.INDEX.LE.27) NUM=15
      DO 10 I=1,NUM
         IF(INDEX.LT.13) WLI(I)=WL1(I)
         IF(INDEX.GE.13) WLI(I)=WL2(I)
         SIGS(I)=SIG0(I,INDEX)
   10 CONTINUE
C
      WLAM=2.997925E18/FR
      IL=1
      IR=NUM
      DO 50 I=1,NUM-1
        IF(WLAM.GE.WLI(I).AND.WLAM.LE.WLI(I+1)) THEN
          IL=I
          IR=I+1
          GO TO 60
        ENDIF
 50   CONTINUE
C
C     LINEAR INTERPOLATION:
C
 60   SIGM=(SIGS(IR)-SIGS(IL))*(WLAM-WLI(IL))/(WLI(IR)-WLI(IL))
     *      + SIGS(IL)
C
C     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE:
C
       IF(WLAM.LE.WLI(1)) SIGM=SIGS(1)
       IF(WLAM.GE.WLI(NUM)) SIGM=SIGS(NUM)
C
C     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION:
C
c       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL)
C
      HIDALG=SIGM*1.E-18
      RETURN
      END
C
C
C     ****************************************************************
C
C
      FUNCTION REIMAN(IB,FR)
C     ======================
C
C     Read table of photon energies and photo-ionization cross-sections
C     from Reilman & Manson (1979, Ap. J. Suppl., 40, 815) for the species 
C     indicated by IB
C    
C     Compute linearly interpolated value of the cross-section 
C     at the frequency FR.
C
C     (At the moment, only a few transitions are considered)
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION HEV(30),F0(30),SIG0(30,2),SIGS(30)
C
      DATA HEV /
     * 130.,160.,190.,210.,240.,270.,300.,330.,360.,390.,
     * 420.,450.,480.,510.,540.,570.,600.,630.,660.,690.,
     * 720.,750.,780.,810.,840.,870.,900.,930.,960.,990./
      DATA SIG0 /
     * 3*0.,     4.422E-1, 3.478E-1, 
     * 2.794E-1, 2.286E-1, 1.899E-1, 1.598E-1, 1.360E-1,
     * 1.169E-1, 1.013E-1, 8.845E-2, 7.776E-2, 6.877E-2,
     * 6.114E-2, 5.463E-2, 4.904E-2, 4.419E-2, 3.998E-2,
     * 3.629E-2, 3.305E-2, 3.019E-2, 2.766E-2, 2.540E-2,
     * 2.339E-2, 2.158E-2, 1.996E-2, 1.850E-2, 1.718E-2,
     * 4*0.,     1.981E-1, 1.584E-1, 
     * 1.290E-1, 1.066E-1, 8.932E-2, 7.567E-2, 6.475E-2,
     * 5.589E-2, 4.862E-2, 4.259E-2, 3.754E-2, 3.329E-2,
     * 2.966E-2, 2.656E-2, 2.388E-2, 2.157E-2, 1.954E-2,
     * 1.777E-2, 1.621E-2, 1.484E-2, 1.362E-2, 1.253E-2,
     * 1.155E-2, 1.067E-2, 9.888E-3, 9.179E-3/
C
      INDEX=-IB-300
      NUM=30  
      DO 10 I=1,NUM
         F0(I)=HEV(I)*2.418573E14
         SIGS(I)=SIG0(I,INDEX)
   10 CONTINUE
C 
      IL=1
      IR=NUM
      DO 50 I=1,NUM-1
        IF(FR.GE.F0(I).AND.FR.LE.F0(I+1)) THEN
          IL=I
          IR=I+1
          GO TO 60
        ENDIF
 50   CONTINUE
C
C     LINEAR INTERPOLATION:
C
 60   SIGM=(SIGS(IR)-SIGS(IL))*(FR-F0(IL))/(F0(IR)-F0(IL))
     *      + SIGS(IL)
C
C     IF OUTSIDE WAVELENGTH RANGE SET TO FIRST(LAST) VALUE:
C
       IF(FR.LE.F0(1)) SIGM=SIGS(1)
       IF(FR.GE.F0(NUM)) SIGM=SIGS(NUM)
C
C     IF LAST NON-ZERO SIG VALUES, NO INTERPOLATION:
C
c       IF(SIGS(IR).EQ.0.) SIGM=SIGS(IL)   
C
      REIMAN=SIGM*1.E-18
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION SBFHE1(II,IB,FR)
C     =========================
C
C     Calculates photoionization cross sections of neutral helium
C     from states with n = 1, 2, 3, 4.
C
C     The levels are either non-averaged (l,s) states, or some
C     averaged levels.
C     The program allows only two standard possibilities of
C     constructing averaged levels:
C     i)  all states within given principal quantum number n (>1) are
C         lumped together
C     ii) all siglet states for given n, and all triplet states for
C         given n are lumped together separately (there are thus two
C         explicit levels for a given n)
C
C     The cross sections are calculated using appropriate averages
C     of the Opacity Project cross sections, calculated by procedure
C     HEPHOT
C
C     Input parameters:
C      II    - index of the lower level (in the numbering of explicit
C              levels)
C      IB    - photoionization switch IBF for the given transition
C            = 10  -  means that the given transition is from an
C                     averaged level
C            = 11  -  the given transition is from non-averaged
C                     singlet state
C            = 13  -  the given transition is from non-averaged
C                     triplet state
C      FR    - frequency
C
      INCLUDE 'PARAMS.FOR'
C
      NI=NQUANT(II)
      IGI=INT(G(II)+0.01)
      IS=IB-10
      sbfhe1=0.
C
C     ----------------------------------------------------------------
C     IB=11 or 13  - photoionization from an non-averaged (l,s) level
C     ----------------------------------------------------------------
C
      IF(IS.EQ.1.OR.IS.EQ.3) THEN
         IL=(IGI/IS-1)/2
         SBFHE1=HEPHOT(IS,IL,NI,FR)
      END IF
C
C     ----------------------------------------------------------------
C     IS=10 - photoionization from an averaged level
C     ----------------------------------------------------------------
C
      IF(IS.EQ.0) THEN
         IF(NI.EQ.2) THEN
C
C ********    photoionization from an averaged level with n=2
C
            IF(IGI.EQ.4) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,2,FR)+3.D0*HEPHOT(1,1,2,FR))/9.D0
            ELSE IF(IGI.EQ.12) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,2,FR)+3.D0*HEPHOT(3,1,2,FR))/9.D0
            ELSE IF(IGI.EQ.16) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,2,FR)+3.D0*(HEPHOT(1,1,2,FR)+
     *                HEPHOT(3,0,2,FR))+9.D0*HEPHOT(3,1,2,FR))/1.6D1
            ELSE
               GO TO 10
            END IF
C
C
C ********    photoionization from an averaged level with n=3
C
         ELSE IF(NI.EQ.3) THEN
            IF(IGI.EQ.9) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,3,FR)+3.D0*HEPHOT(1,1,3,FR)+
     *                5.D0*HEPHOT(1,2,3,FR))/9.D0
            ELSE IF(IGI.EQ.27) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,3,FR)+3.D0*HEPHOT(3,1,3,FR)+
     *                5.D0*HEPHOT(3,2,3,FR))/9.D0
            ELSE IF(IGI.EQ.36) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,3,FR)+3.D0*HEPHOT(1,1,3,FR)+
     *                5.D0*HEPHOT(1,2,3,FR)+
     *                3.D0*HEPHOT(3,0,3,FR)+9.D0*HEPHOT(3,1,3,FR)+
     *                15.D0*HEPHOT(3,2,3,FR))/3.6D0
            ELSE
               GO TO 10
            END IF
         ELSE IF(NI.EQ.4) THEN
C
C ********    photoionization from an averaged level with n=4
C
            IF(IGI.EQ.16) THEN
C
C      a) lower level is an averaged singlet state
C
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+
     *                 5.D0*HEPHOT(1,2,4,FR)+
     *                 7.D0*HEPHOT(1,3,4,FR))/1.6D1
            ELSE IF(IGI.EQ.48) THEN
C
C      b) lower level is an averaged triplet state
C
               SBFHE1=(HEPHOT(3,0,4,FR)+3.D0*HEPHOT(3,1,4,FR)+
     *                 5.D0*HEPHOT(3,2,4,FR)+
     *                 7.D0*HEPHOT(3,3,4,FR))/1.6D1
            ELSE IF(IGI.EQ.64) THEN
C
C      c) lower level is an average of both singlet and triplet states
C
               SBFHE1=(HEPHOT(1,0,4,FR)+3.D0*HEPHOT(1,1,4,FR)+
     *                 5.D0*HEPHOT(1,2,4,FR)+
     *                 7.D0*HEPHOT(1,3,4,FR)+
     *                 3.D0*HEPHOT(3,0,4,FR)+
     *                 9.D0*HEPHOT(3,1,4,FR)+
     *                 15.D0*HEPHOT(3,2,4,FR)+
     *                 21.D0*HEPHOT(3,3,4,FR))/6.4D1
            ELSE
               GO TO 10
            END IF
         ELSE
            GO TO 10
         END IF
      END IF
      RETURN
   10 WRITE(6,601) NI,IGI,IS
  601 FORMAT(1H0/' INCONSISTENT INPUT TO PROCEDURE SBFHE1'/
     * ' QUANTUM NUMBER =',I3,'  STATISTICAL WEIGHT',I4,'  S=',I3)
      STOP
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION HEPHOT(S,L,N,FREQ)
C     ===========================
C
C           EVALUATES HE I PHOTOIONIZATION CROSS SECTION USING SEATON
C           FERNLEY'S CUBIC FITS TO THE OPACITY PROJECT CROSS SECTIONS
C           UP TO SOME ENERGY "EFITM" IN THE RESONANCE-FREE ZONE.  BEYOND
C           THIS ENERGY LINEAR FITS TO LOG SIGMA IN LOG (E/E0) ARE USED.
C           THIS EXTRAPOLATION SHOULD BE USED UP TO THE BEGINNING OF THE
C           RESONANCE ZONE "XMAX", BUT AT PRESENT IT IS USED THROUGH IT.
C           BY CHANGING A FEW LINES THAT ARE PRESENTLY COMMENTED OUT,
C           FOR ENERGIES IN THE RESONANCE ZONE A VALUE OF 1/100 OF THE
C           THRESHOLD CROSS SECTION IS USED -- THIS IS PURELY AD HOC AND
C           ONLY A TEMPORARY MEASURE.  OBVIOUSLY ANY OTHER VALUE OR FUNCTIONAL
C           FORM CAN BE INSERTED HERE.
C
C           CALLING SEQUENCE INCLUDES:
C                S = MULTIPLICITY, EITHER 1 OR 3
C                L = ANGULAR MOMENTUM, 0, 1, OR 2;
C                    for L > 2 - hydrogenic expresion
C                FREQ = FREQUENCY
C
C           DGH JUNE 1988 JILA, slightly modified by I.H.
C
      INCLUDE 'PARAMS.FOR'
      INTEGER S,L,SS,LL
      DIMENSION COEF(4,53),IST(3,2),N0(3,2),
     *          FL0(53),A(53),B(53),XFITM(53)
c      DIMENSION XMAX(53)
C
      DATA IST/1,36,20,11,45,28/
      DATA N0/1,2,3,2,2,3/
C
      DATA FL0/
     . 2.521D-01,-5.381D-01,-9.139D-01,-1.175D+00,-1.375D+00,-1.537D+00,
     .-1.674D+00,-1.792D+00,-1.896D+00,-1.989D+00,-4.555D-01,-8.622D-01,
     .-1.137D+00,-1.345D+00,-1.512D+00,-1.653D+00,-1.774D+00,-1.880D+00,
     .-1.974D+00,-9.538D-01,-1.204D+00,-1.398D+00,-1.556D+00,-1.690D+00,
     .-1.806D+00,-1.909D+00,-2.000D+00,-9.537D-01,-1.204D+00,-1.398D+00,
     .-1.556D+00,-1.690D+00,-1.806D+00,-1.909D+00,-2.000D+00,-6.065D-01,
     .-9.578D-01,-1.207D+00,-1.400D+00,-1.558D+00,-1.692D+00,-1.808D+00,
     .-1.910D+00,-2.002D+00,-5.749D-01,-9.352D-01,-1.190D+00,-1.386D+00,
     .-1.547D+00,-1.682D+00,-1.799D+00,-1.902D+00,-1.995D+00/
C
      DATA XFITM/
     . 3.262D-01, 6.135D-01, 9.233D-01, 8.438D-01, 1.020D+00, 1.169D+00,
     . 1.298D+00, 1.411D+00, 1.512D+00, 1.602D+00, 7.228D-01, 1.076D+00,
     . 1.206D+00, 1.404D+00, 1.481D+00, 1.464D+00, 1.581D+00, 1.685D+00,
     . 1.777D+00, 9.586D-01, 1.187D+00, 1.371D+00, 1.524D+00, 1.740D+00,
     . 1.854D+00, 1.955D+00, 2.046D+00, 9.585D-01, 1.041D+00, 1.371D+00,
     . 1.608D+00, 1.739D+00, 1.768D+00, 1.869D+00, 1.803D+00, 7.360D-01,
     . 1.041D+00, 1.272D+00, 1.457D+00, 1.611D+00, 1.741D+00, 1.855D+00,
     . 1.870D+00, 1.804D+00, 9.302D-01, 1.144D+00, 1.028D+00, 1.210D+00,
     . 1.362D+00, 1.646D+00, 1.761D+00, 1.863D+00, 1.954D+00/
C
      DATA A/
     . 6.95319D-01, 1.13101D+00, 1.36313D+00, 1.51684D+00, 1.64767D+00,
     . 1.75643D+00, 1.84458D+00, 1.87243D+00, 1.85628D+00, 1.90889D+00,
     . 9.01802D-01, 1.25389D+00, 1.39033D+00, 1.55226D+00, 1.60658D+00,
     . 1.65930D+00, 1.68855D+00, 1.62477D+00, 1.66726D+00, 1.83599D+00,
     . 2.50403D+00, 3.08564D+00, 3.56545D+00, 4.25922D+00, 4.61346D+00,
     . 4.91417D+00, 5.19211D+00, 1.74181D+00, 2.25756D+00, 2.95625D+00,
     . 3.65899D+00, 4.04397D+00, 4.13410D+00, 4.43538D+00, 4.19583D+00,
     . 1.79027D+00, 2.23543D+00, 2.63942D+00, 3.02461D+00, 3.35018D+00,
     . 3.62067D+00, 3.85218D+00, 3.76689D+00, 3.49318D+00, 1.16294D+00,
     . 1.86467D+00, 2.02110D+00, 2.24231D+00, 2.44240D+00, 2.76594D+00,
     . 2.93230D+00, 3.08109D+00, 3.21069D+00/
C
      DATA B/
     .-1.29000D+00,-2.15771D+00,-2.13263D+00,-2.10272D+00,-2.10861D+00,
     .-2.11507D+00,-2.11710D+00,-2.08531D+00,-2.03296D+00,-2.03441D+00,
     .-1.85905D+00,-2.04057D+00,-2.02189D+00,-2.05930D+00,-2.03403D+00,
     .-2.02071D+00,-1.99956D+00,-1.92851D+00,-1.92905D+00,-4.58608D+00,
     .-4.40022D+00,-4.39154D+00,-4.39676D+00,-4.57631D+00,-4.57120D+00,
     .-4.56188D+00,-4.55915D+00,-4.41218D+00,-4.12940D+00,-4.24401D+00,
     .-4.40783D+00,-4.39930D+00,-4.25981D+00,-4.26804D+00,-4.00419D+00,
     .-4.47251D+00,-3.87960D+00,-3.71668D+00,-3.68461D+00,-3.67173D+00,
     .-3.65991D+00,-3.64968D+00,-3.48666D+00,-3.23985D+00,-2.95758D+00,
     .-3.07110D+00,-2.87157D+00,-2.83137D+00,-2.82132D+00,-2.91084D+00,
     .-2.91159D+00,-2.91336D+00,-2.91296D+00/
C
      DATA ((COEF(I,J),I=1,4),J=1,10)/
     . 8.734D-01,-1.545D+00,-1.093D+00, 5.918D-01, 9.771D-01,-1.567D+00,
     .-4.739D-01,-1.302D-01, 1.174D+00,-1.638D+00,-2.831D-01,-3.281D-02,
     . 1.324D+00,-1.692D+00,-2.916D-01, 9.027D-02, 1.445D+00,-1.761D+00,
     .-1.902D-01, 4.401D-02, 1.546D+00,-1.817D+00,-1.278D-01, 2.293D-02,
     . 1.635D+00,-1.864D+00,-8.252D-02, 9.854D-03, 1.712D+00,-1.903D+00,
     .-5.206D-02, 2.892D-03, 1.782D+00,-1.936D+00,-2.952D-02,-1.405D-03,
     . 1.845D+00,-1.964D+00,-1.152D-02,-4.487D-03/
      DATA ((COEF(I,J),I=1,4),J=11,19)/
     . 7.377D-01,-9.327D-01,-1.466D+00, 6.891D-01, 9.031D-01,-1.157D+00,
     .-7.151D-01, 1.832D-01, 1.031D+00,-1.313D+00,-4.517D-01, 9.207D-02,
     . 1.135D+00,-1.441D+00,-2.724D-01, 3.105D-02, 1.225D+00,-1.536D+00,
     .-1.725D-01, 7.191D-03, 1.302D+00,-1.602D+00,-1.300D-01, 7.345D-03,
     . 1.372D+00,-1.664D+00,-8.204D-02,-1.643D-03, 1.434D+00,-1.715D+00,
     .-4.646D-02,-7.456D-03, 1.491D+00,-1.760D+00,-1.838D-02,-1.152D-02/
      DATA ((COEF(I,J),I=1,4),J=20,27)/
     . 1.258D+00,-3.442D+00,-4.731D-01,-9.522D-02, 1.553D+00,-2.781D+00,
     .-6.841D-01,-4.083D-03, 1.727D+00,-2.494D+00,-5.785D-01,-6.015D-02,
     . 1.853D+00,-2.347D+00,-4.611D-01,-9.615D-02, 1.955D+00,-2.273D+00,
     .-3.457D-01,-1.245D-01, 2.041D+00,-2.226D+00,-2.669D-01,-1.344D-01,
     . 2.115D+00,-2.200D+00,-1.999D-01,-1.410D-01, 2.182D+00,-2.188D+00,
     .-1.405D-01,-1.460D-01/
      DATA ((COEF(I,J),I=1,4),J=28,35)/
     . 1.267D+00,-3.417D+00,-5.038D-01,-1.797D-02, 1.565D+00,-2.781D+00,
     .-6.497D-01,-5.979D-03, 1.741D+00,-2.479D+00,-6.099D-01,-2.227D-02,
     . 1.870D+00,-2.336D+00,-4.899D-01,-6.616D-02, 1.973D+00,-2.253D+00,
     .-3.972D-01,-8.729D-02, 2.061D+00,-2.212D+00,-3.072D-01,-1.060D-01,
     . 2.137D+00,-2.189D+00,-2.352D-01,-1.171D-01, 2.205D+00,-2.186D+00,
     .-1.621D-01,-1.296D-01/
      DATA ((COEF(I,J),I=1,4),J=36,44)/
     . 1.129D+00,-3.149D+00,-1.910D-01,-5.244D-01, 1.431D+00,-2.511D+00,
     .-3.710D-01,-1.933D-01, 1.620D+00,-2.303D+00,-3.045D-01,-1.391D-01,
     . 1.763D+00,-2.235D+00,-1.829D-01,-1.491D-01, 1.879D+00,-2.215D+00,
     .-9.003D-02,-1.537D-01, 1.978D+00,-2.213D+00,-2.066D-02,-1.541D-01,
     . 2.064D+00,-2.220D+00, 3.258D-02,-1.527D-01, 2.140D+00,-2.225D+00,
     . 6.311D-02,-1.455D-01, 2.208D+00,-2.229D+00, 7.977D-02,-1.357D-01/
      DATA ((COEF(I,J),I=1,4),J=45,53)/
     . 1.204D+00,-2.809D+00,-3.094D-01, 1.100D-01, 1.455D+00,-2.254D+00,
     .-4.795D-01, 6.872D-02, 1.619D+00,-2.109D+00,-3.357D-01,-2.532D-02,
     . 1.747D+00,-2.065D+00,-2.317D-01,-5.224D-02, 1.853D+00,-2.058D+00,
     .-1.517D-01,-6.647D-02, 1.943D+00,-2.055D+00,-1.158D-01,-6.081D-02,
     . 2.023D+00,-2.070D+00,-6.470D-02,-6.800D-02, 2.095D+00,-2.088D+00,
     .-2.357D-02,-7.250D-02, 2.160D+00,-2.107D+00, 1.065D-02,-7.542D-02/
C
      IF(L.GT.2) GO TO 20
C
C          SELECT BEGINNING AND END OF COEFFICIENTS
C
      SS=(S+1)/2
      LL=L+1
      NSL0=N0(LL,SS)
      I=IST(LL,SS)+N-NSL0
C
C          EVALUATE CROSS SECTION
C
      FL=LOG10(FREQ/3.28805E15)
      X=FL-FL0(I)
      IF(X.GE.-0.001D0) THEN
         IF(X.LT.XFITM(I)) THEN
            P=COEF(4,I)
            DO 10 K=1,3
               P=X*P+COEF(4-K,I)
   10       CONTINUE
            HEPHOT=1.D-18*1.D1**P
          ELSE
C           OTHERWISE REMOVE INSTRUCTION AND 3 FOLLOWING "C"
C         ELSE IF(X.LT.XMAX(I)) THEN
            HEPHOT=1.D-18*1.D1**(A(I)+B(I)*X)
C         ELSE
C           HEPHOT=1.D-18*1.D1**(COEF(1,I)-2.0D0)
          END IF
       ELSE
          HEPHOT=0.
      END IF
      RETURN
C
C     Hydrogenic expression for L > 2
C      [multiplied by relative population of state (s,l,n), ie.
C       by  stat.weight(s,l)/stat.weight(n)]
C
   20 GN=2.D0*N*N
      HEPHOT=2.815D29/FREQ/FREQ/FREQ/N**5*(2*L+1)*S/GN
      RETURN
      END
C
C
C     ****************************************************************
C
C
 
      FUNCTION TOPBAS(FREQ,FREQ0,TYPLV)
C     ==================================
C
C     Procedure calculates the photo-ionisation cross section SIGMA in 
C     [cm^2] at frequency FREQ. FREQ0 is the threshold frequency from
C     level I of ion KI. Threshold cross-sections will be of the order
C     of the numerical value of 10^-18.
C     Opacity-Project (OP) interpolation fit formula
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER    (E10=2.3025851)
      PARAMETER    (MMAXOP = 200,! maximum number of levels in OP data
     +              MOP    =  15 )! maximum number of fit points per level
      CHARACTER*10  IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
      CHARACTER*10  TYPLV
      COMMON /TOPB/ SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
     +              XOP(MOP,MMAXOP) ,! x = alog10(nu/nu0) of fit point 
     +              NOP(MMAXOP)     ,! number of fit points for current level
     +              NTOTOP          ,! total number of levels in OP data
     +              IDLVOP         ,! level identifyer Opacity-Project data
     +              LOPREA   ! .T. OP data read in; .F. OP data not yer read in
      DIMENSION XFIT(MOP) ,! local array containing x     for OP data
     +          SFIT(MOP)  ! local array containing sigma for OP data
C
C     Read OP data if not yet done
C
      TOPBAS=0.
      IF (.NOT.LOPREA) CALL OPDATA
      X = LOG10(FREQ/FREQ0)
      DO IOP = 1,NTOTOP
         IF (IDLVOP(IOP).EQ.TYPLV) THEN
C           level has been detected in OP-data file
            IF (NOP(IOP).LE.0) GO TO 20
            DO IFIT = 1,NOP(IOP)
               XFIT(IFIT) = XOP(IFIT,IOP)
               SFIT(IFIT) = SOP(IFIT,IOP)
            END DO
            SIGM  = YLINTP (X,XFIT,SFIT,NOP(IOP),MOP)
            SIGM  = 1.D-18*EXP(E10*SIGM)
            TOPBAS=SIGM
            GO TO 10
         END IF
      END DO
   10 RETURN
C     Level is not found ,or no data for this level, in RBF.DAT
   20 WRITE (61,100) TYPLV    
  100 FORMAT ('SIGMA.......: OP DATA NOT AVAILABLE FOR LEVEL ',A10)
      RETURN
      END
C

C     ******************************************************************
C
C
      SUBROUTINE OPDATA
C     =================
C
C     Procedure reads photo-ionization cross sections fit coefficients
C     based on Opacity-Project (OP) data from file RBF.DAT
C     Data, as stored, requires linear interpolation.
C
C     Meaning of global variables:
C        NTOTOP    = total number of levels in Opacity Project data
C        IDLVOP() = level identifyer of current level
C        NOP()     = number of fit points for current level
C        XOP(,)    = x     = alog10(nu/nu0)       of fit point
C        SOP(,)    = sigma = alog10(sigma/10^-18) of fit point 
C      
      INCLUDE 'PARAMS.FOR'
      PARAMETER    (MMAXOP = 200,! maximum number of levels in OP data
     +              MOP    =  15 )! maximum number of fit points per level
      CHARACTER*10  IDLVOP(MMAXOP) ! level identifyer Opacity-Project data
      COMMON /TOPB/ SOP(MOP,MMAXOP) ,! sigma = alog10(sigma/10^-18) of fit point
     +              XOP(MOP,MMAXOP) ,! x = alog10(nu/nu0) of fit point 
     +              NOP(MMAXOP)     ,! number of fit points for current level
     +              NTOTOP          ,! total number of levels in OP data
     +              IDLVOP         ,! level identifyer Opacity-Project data
     +              LOPREA   ! .T. OP data read in; .F. OP data not yer read in
      CHARACTER*4 IONID
C
      OPEN (UNIT=40,FILE='RBF.DAT',STATUS='OLD')
C     Skip header
      DO IREAD = 1, 21
         READ (40,*)
      END DO
      IOP = 0
C         = initialize sequential level index op Opacity Project data 
C     Read number of elements in file
      READ (40,*) NEOP
      DO IEOP = 1, NEOP
C        Skip element name header
         DO IREAD = 1, 3
            READ (40,*)
         END DO
C        Read number of ionization stages of current element in  file
         READ (40,*) NIOP
         DO IIOP = 1, NIOP
C           Read ion identifyer, atomic & electron number, # of levels 
C           for current ion
            READ (40,*) IONID, IATOM_OP, IELEC_OP, NLEVEL_OP
            DO ILOP = 1, NLEVEL_OP
C              Increase sequential level index of Opacity Project data
               IOP = IOP+1
C              Read level identifyer and number of sigma fit points
               READ (40,*) IDLVOP(IOP), NOP(IOP)
C              Read normalized log10 frequency and log10 cross section values
               DO IS = 1, NOP(IOP)
                  READ (40,*) INDEX, XOP(IS,IOP), SOP(IS,IOP)
               END DO
            END DO
         END DO
      END DO
      NTOTOP  = IOP
C             = total number of levels in Opacity Project data
      LOPREA  = .TRUE.
C             = set flag as data has been read in
C
      RETURN
      END
C
C
C
C     ******************************************************************
C
C
      FUNCTION YLINTP (XINT,X,Y,N,NTOT)
C     =================================
C
C     linear interpolation routine. Determines YINT = Y(XINT) from 
C     grid Y(X) with N points and dimension NTOT. 
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION X(NTOT),Y(NTOT)
C
C     bisection (see Numerical Recipes par 3.4 page 90)
      JL = 0
      JU = N+1
10    IF (JU-JL.GT.1) THEN
         JM = (JU+JL)/2
         IF ((X(N).GT.X(1)).EQV.(XINT.GT.X(JM))) THEN
            JL = JM
         ELSE
            JU = JM
         END IF
         GO TO 10
      END IF
      J = JL
      IF (J.EQ.N) J = J-1
      IF (J.EQ.0) J = J+1
      RC         = (Y(J+1)-Y(J))/(X(J+1)-X(J))
      YLINTP = RC*(XINT-X(J))+Y(J)
C
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE OPAC(ID,CROSS,ABSO,EMIS,SCAT)
C     ========================================
C
C     Absorption, emission, and scattering coefficients
C     at depth ID and for several frequencies (some or all)
C
C     Input: ID    - depth index
C            CROSS - two dimensional array of photoionization
C                    cross-sections
C     Output: ABSO - array of absorption coefficient
C             EMIS - array of emission coefficient
C             SCAT - array of scattering coefficient (all scattering
C                    mechanisms except electron scattering)
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      DIMENSION CROSS(MCROSS,MFRQ)
      DIMENSION ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ),
     *          ABLIN(MFREQ),EMLIN(MFREQ)
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/dissol/fropc(mlevel),indexp(mlevel)
      PARAMETER (UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.694E8)
C
      IF(IMODE.EQ.-1.AND.ID.NE.IDSTD) RETURN
      T=TEMP(ID)
      ANE=ELEC(ID)
      T1=UN/T
      HKT=HK*T1
      TK=HKT/H
      SRT=UN/SQRT(T)
      SGFF=CFF*SRT
      CON=CSB*T1*SRT
      conts=1.e-36/con
      ABLY=0.
      EMLY=0.
      SCLY=0.
      sce=ane*sige
      IJ0=2
      IF(NFREQ.EQ.1) IJ0=1
      IF(IMODE.EQ.2) IJ0=NFREQ
      M=3
      IF(ICONTL.EQ.1) M=1
C
C     Opacity and emissivity in continuum
C     **** calculated only in the first and the last frequency *****
C
      DO 200 IJ=1,IJ0
         FR=FREQ(IJ)
         FR15=FR*TEN15
         BNU=BN*FR15*FR15*FR15
         HKF=HKT*FR
         ABF=0.
         EBF=0.
         AFF=0.
         DO 100 IL=1,NION
            N0I=NFIRST(IL)
            N1I=NLAST(IL)
            NKE=NNEXT(IL)
            XN=POPUL(NKE,ID)
C
C           Bound-free contribution + possibly
c            pseudo-continuum (accounting for dissolved fraction)
C
            DO 10 II=N0I,N1I
               SG=0.
               IF(IFWOP(II).LT.0) THEN
                  SG=SGMERG(II,ID,FR)
                ELSE
                  SG=CROSS(II,IJ)
                  IF(INDEXP(II).EQ.5) THEN
                     IZZ=IZ(IEL(II))
                     FR0=ENION(II)/6.6256E-27
                     CALL DWNFR1(FR,FR0,ID,IZZ,DW1)
                     SG=SG*DW1
                  END IF
               END IF
               if(sg.le.0.) go to 10
               ABF=ABF+SG*POPUL(II,ID)
               XX=SG*XN*EXP(ENION(II)*TK)*WOP(II,ID)
               IF(XX.lt.conts) go to 10
               EBF=EBF+XX*CON*G(II)/G(NKE)
   10       CONTINUE

            IT=IFREE(IL)
            IF(IT.EQ.0) GO TO 100
C
C           Free-free contribution
C
            IE=IL
            IF(IE.EQ.IELHM) GO TO 65
            CH=IZ(IL)*IZ(IL)
            SF1=CH*XN*SGFF/(FR*FR*FR)
C
C           The following expression is the so-called modified free-free
C           opacity, ie. allowing for the photoionization from higher,
C           non-explicit, LTE energy levels of the ion IL
C
            HKFM=HKT*MIN(FF(IL),FR)
            SF2=EXP(HKFM)
            IF(IT.NE.2) GO TO 50
            SG=GFREE(T,FR/CH)
            SF2=SF2+SG-UN
   50       SFF=SF1*SF2
            GO TO 70
   65       SFF=SFFHMI(XN,FR,T)
   70       AFF=AFF+SFF
  100    CONTINUE
C
C        Additional opacities
C
         CALL OPADD(0,ID,FR,ABAD,EMAD,SCAD)   
         IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY)
C
C        Total opacity and emissivity
C
         X=EXP(-HKF)
         X1=UN-X
         BNE=BNU*X*ANE
c        ABSO(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ABAD+ABLY
         ABSO(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ABAD
         EMIS(IJ)=BNE*(AFF+EBF)+EMAD+EMLY
         SCAT(IJ)=SCAD+SCLY+sce
         IF(IJ.EQ.1) THEN
         ABLY1=ABLY
         EMLY1=EMLY
         SCLY1=SCLY
         END IF
  200 CONTINUE
      AVAB=(ABSO(1)+ABSO(2)+SCAT(1)+SCAT(2))*0.5*RELOP
      IF(NFREQ.LE.2.OR.IMODE.EQ.-1) RETURN
      IF(IMODE.EQ.2) GO TO 225
C
C     interpolated continuum opacity, emissivity, and scattering
C     for all frequencies
C
      DO IJ=3,NFREQ
         ABSO(IJ)=FRX1(IJ)*ABSO(2)+FRX2(IJ)*ABSO(1)
         EMIS(IJ)=FRX1(IJ)*EMIS(2)+FRX2(IJ)*EMIS(1)
         SCAT(IJ)=FRX1(IJ)*SCAT(2)+FRX2(IJ)*SCAT(1)
      END DO
C
C    hydrogen lines -- for IHYL = 0
C     *** calculated only for the first and the last frequency
C     and interpolated hydrogen line opacity and emissivity
C     for all frequencies
C
      IF(IHYL.EQ.0) THEN
      CALL HYDLIN(ID,1,2,ABLIN,EMLIN)
      DO IJ=M,NFREQ
         ABSO(IJ)=ABSO(IJ)+FRX1(IJ)*ABLIN(2)+FRX2(IJ)*ABLIN(1)
         EMIS(IJ)=EMIS(IJ)+FRX1(IJ)*EMLIN(2)+FRX2(IJ)*EMLIN(1)
      END DO
      END IF
C
C     **** Opacity and emissivity in lines ****
C
      CALL LINOP(ID,ABLIN,EMLIN,AVAB)
      DO IJ=3,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
C
C     **** Opacity and emissivity in molecular lines ****
C
      if(ifmol.gt.0) then
      do ilist=1,nmlist
      CALL MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST)
      DO IJ=3,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
      end do
      end if
  225 CONTINUE
C
C     **** Detailed opacity and emissivity in hydrogen lines ****
C          (for IHYL=1)
C
      IF(IHYL.GT.0.OR.IMODE.EQ.2) THEN
      CALL HYDLIN(ID,M,NFREQ,ABLIN,EMLIN)
      DO IJ=M,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
      END IF
C
C     **** Detailed opacity and emissivity in HE II lines ****
C          (for IHE2L=1)
C
      IF(IHE2L.GT.0) THEN
      CALL HE2LIN(ID,M,NFREQ,ABLIN,EMLIN)
      DO IJ=M,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
      END IF
C
C     opacity due to detailed photoinization cross-section
C     (from tables; including resonance features)
C     The two routines may be called and correspond to different formats
C     as well as difference in INPUT!
C
      CALL PHTION(ID,ABSO,EMIS,FREQ,NFREQ)
      CALL PHTX(ID,ABSO,EMIS,FREQ,0)
C
       if(imode.ge.0) then
          do ij=1,nfreq
             abso(ij)=abso(ij)+scat(ij)
          end do
       end if
C
      IF(ICONTL.EQ.1) RETURN
      ABSO(1)=ABSO(1)-ABLY1
      EMIS(1)=EMIS(1)-EMLY1
      SCAT(1)=SCAT(1)-SCLY1
      ABSO(2)=ABSO(2)-ABLY
      EMIS(2)=EMIS(2)-EMLY
      SCAT(2)=SCAT(2)-SCLY
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE OPACW(ID,CROSS,ABSO,EMIS,
     *                 ABSOC,EMISC,SCATC,MODC)
C     ========================================
C
C     Absorption, emission, and scattering coefficients
C     at depth ID and for several frequencies (some or all)
C
C     Input: ID    - depth index
C            CROSS - two dimensional array of photoionization
C                    cross-sections
C     Output: ABSO - array of absorption coefficient
C             EMIS - array of emission coefficient
C             SCAT - array of scattering coefficient (all scattering
C                    mechanisms except electron scattering)
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      DIMENSION CROSS(MCROSS,MFRQ)
      DIMENSION ABSO(MFREQ),EMIS(MFREQ),SCAT(MFREQ),
     *          ABSOC(MFREQC),EMISC(MFREQC),SCATC(MFREQC),
     *          ABLIN(MFREQ),EMLIN(MFREQ),
     *          ABL1(MFREQC),EML1(MFREQC),SCL1(MFREQC)
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/dissol/fropc(mlevel),indexp(mlevel)
      common/lasers/lasdel
      PARAMETER (UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.694E8)
C
      IF(IMODE.EQ.-1.AND.ID.NE.IDSTD) RETURN
      T=TEMP(ID)
      ANE=ELEC(ID)
      T1=UN/T
      HKT=HK*T1
      TK=HKT/H
      SRT=UN/SQRT(T)
      SGFF=CFF*SRT
      CON=CSB*T1*SRT
      conts=1.e-36/con
      ABLY=0.
      EMLY=0.
      SCLY=0.
      IJ0=2
      IF(NFREQ.EQ.1) IJ0=1
      IF(IMODE.EQ.2) IJ0=NFREQ
      M=3
C
C     Opacity and emissivity in continuum
C     **** calculated only for the continuum frequencies *****
C
      DO 200 IJ=1,NFREQC
         FR=FREQC(IJ)
         FR15=FR*TEN15
         BNU=BN*FR15*FR15*FR15
         HKF=HKT*FR
         ABF=0.
         EBF=0.
         AFF=0.
         DO 100 IL=1,NION
            N0I=NFIRST(IL)
            N1I=NLAST(IL)
            NKE=NNEXT(IL)
            XN=POPUL(NKE,ID)
C
C           Bound-free contribution + possibly
c            pseudo-continuum (accounting for dissolved fraction)
C
            DO 10 II=N0I,N1I
               SG=0.
               IF(IFWOP(II).LT.0) THEN
                  SG=SGMERG(II,ID,FR)
                ELSE
                  SG=CROSS(II,IJ)
                  IF(INDEXP(II).EQ.5) THEN
                     IZZ=IZ(IEL(II))
                     FR0=ENION(II)/6.6256E-27
                     CALL DWNFR1(FR,FR0,ID,IZZ,DW1)
                     SG=SG*DW1
                  END IF
               END IF
               ABF=ABF+SG*POPUL(II,ID)
               XX=SG*XN*EXP(ENION(II)*TK)*WOP(II,ID)
               IF(XX.lt.conts) go to 10
               EBF=EBF+XX*CON*G(II)/G(NKE)
   10       CONTINUE
            IT=IFREE(IL)
            IF(IT.EQ.0) GO TO 100
C
C           Free-free contribution
C
            IE=IL
            IF(IE.EQ.IELHM) GO TO 65
            CH=IZ(IL)*IZ(IL)
            SF1=CH*XN*SGFF/(FR*FR*FR)
C
C           The following expression is the so-called modified free-free
C           opacity, ie. allowing for the photoionization from higher,
C           non-explicit, LTE energy levels of the ion IL
C
            HKFM=HKT*MIN(FF(IL),FR)
            SF2=EXP(HKFM)
            IF(IT.NE.2) GO TO 50
            SG=GFREE(T,FR/CH)
            SF2=SF2+SG-UN
   50       SFF=SF1*SF2
            GO TO 70
   65       SFF=SFFHMI(XN,FR,T)
   70       AFF=AFF+SFF
  100    CONTINUE
C
C        Additional opacities
C
         CALL OPADD(0,ID,FR,ABAD,EMAD,SCAD)   
         IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY)
C
C        Total opacity and emissivity
C
         X=EXP(-HKF)
         X1=UN-X
         BNE=BNU*X*ANE
         ABSOC(IJ)=ABF+ANE*(X1*AFF-X*EBF)+ANE*SIGE+ABAD+ABLY
         EMISC(IJ)=BNE*(AFF+EBF)+EMAD+EMLY
         SCATC(IJ)=SCAD+SCLY
         ABL1(IJ)=ABLY
         EML1(IJ)=EMLY
         SCL1(IJ)=SCLY
  200 CONTINUE
c       
      if(modc.eq.0) return
c       
      IF(NFREQ.LE.2.OR.IMODE.EQ.-1) RETURN
C
C     interpolated continuum and hydrogen line opacity and emissivity
C     for all frequencies
C
      DO IJ=1,NFREQ
         IJC=IJCINT(IJ)
         ABSO(IJ)=FRX1(IJ)*ABSOC(IJC)+(1.-FRX1(IJ))*ABSOC(IJC+1)
         EMIS(IJ)=FRX1(IJ)*EMISC(IJC)+(1.-FRX1(IJ))*EMISC(IJC+1)
         SCAT(IJ)=FRX1(IJ)*SCATC(IJC)+(1.-FRX1(IJ))*SCATC(IJC+1)
      END DO
      IF(IMODE.EQ.2) GO TO 225
C
C     **** Opacity and emissivity in lines ****
C
      CALL LINOPW(ID,ABLIN,EMLIN)
      DO IJ=1,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
C
C     **** Opacity and emissivity in molecular lines ****
C
      if(ifmol.gt.0) then
      do ilist=1,nmlist
      CALL MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST)
      DO IJ=1,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
      end do
      end if
  225 CONTINUE
C
C     **** Detailed opacity and emissivity in hydrogen lines ****
C
      CALL HYDLIW(ID,ABLIN,EMLIN)
      DO IJ=1,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
C
C     **** Detailed opacity and emissivity in HE II lines ****
C          (for IHE2L=1)
C
      CALL HE2LIW(ID,ABLIN,EMLIN)
      DO IJ=1,NFREQ
         ABSO(IJ)=ABSO(IJ)+ABLIN(IJ)
         EMIS(IJ)=EMIS(IJ)+EMLIN(IJ)
      END DO
C
C     opacity due to detailed photoinization cross-section
C     (from tables; including resonance features)
C     The two routines may be called and correspond to different formats
C     as well as difference in INPUT!
C
      CALL PHTION(ID,ABSO,EMIS,FREQ,NFREQ)
      CALL PHTX(ID,ABSO,EMIS,FREQ,0)
C
      IF(ICONTL.EQ.1) RETURN
      DO IJ=1,NFREQC
         ABSOC(IJ)=ABSOC(IJ)-ABL1(IJ)
         EMISC(IJ)=EMISC(IJ)-EML1(IJ)
         SCATC(IJ)=SCATC(IJ)-SCL1(IJ)
      END DO
      RETURN
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE OPACON(ID,CROSS,ABSOC,EMISC,SCATC)
C     ============================================
C
C     Absorption, emission, and scattering coefficients
C     at depth ID and for several frequencies (some or all)
C
C     Input: ID    - depth index
C            CROSS - two dimensional array of photoionization
C                    cross-sections
C     Output: ABSO - array of absorption coefficient
C             EMIS - array of emission coefficient
C             SCAT - array of scattering coefficient
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      DIMENSION CROSS(MCROSS,MFRQ)
      DIMENSION ABSOC(MFREQC),EMISC(MFREQC),SCATC(MFREQC)
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/dissol/fropc(mlevel),indexp(mlevel)
      PARAMETER (UN=1.,TEN15=1.E-15,CSB=2.0706E-16,CFF=3.694E8)
C
      T=TEMP(ID)
      ANE=ELEC(ID)
      T1=UN/T
      HKT=HK*T1
      TK=HKT/H
      SRT=UN/SQRT(T)
      SGFF=CFF*SRT
      CON=CSB*T1*SRT
      ABLY=0.
      EMLY=0.
      SCLY=0.
      sce=ane*sige
C
C     Opacity and emissivity in continuum
C     **** calculated only for the continuum frequencies *****
C
      DO 200 IJ=1,NFREQC
         FR=FREQC(IJ)
         FR15=FR*TEN15
         BNU=BN*FR15*FR15*FR15
         HKF=HKT*FR
         ABF=0.
         EBF=0.
         AFF=0.
         DO 100 IL=1,NION
            N0I=NFIRST(IL)
            N1I=NLAST(IL)
            NKE=NNEXT(IL)
            XN=POPUL(NKE,ID)
C
C           Bound-free contribution + possibly
c           pseudo-continuum (accounting for dissolved fraction)
C
            DO 10 II=N0I,N1I
               SG=0.
               IF(IFWOP(II).LT.0) THEN
                  SG=SGMERG(II,ID,FR)
                ELSE
                  SG=CROSS(II,IJ)
                  if(sg.le.0.) go to 10
                  IF(INDEXP(II).EQ.5) THEN
                     IZZ=IZ(IEL(II))
                     FR0=ENION(II)/6.6256E-27
                     CALL DWNFR1(FR,FR0,ID,IZZ,DW1)
                     SG=SG*DW1
                  END IF
               END IF
               if(popul(ii,id).lt.1.e-20.or.xn.lt.1.e-20) go to 10
               ABF=ABF+SG*POPUL(II,ID)
               XX=SG*XN*EXP(ENION(II)*TK-hkf)*WOP(II,ID)
               ee=exp(enion(ii)*tk-hkf)
               EBF=EBF+XX*CON*G(II)/G(NKE)
c            if(id.eq.1.or.id.eq.50) write(*,*)'opacon',id,ij,ii,
c    *       popul(ii,id),sg,abf
   10       CONTINUE
            IT=IFREE(IL)
            IF(IT.EQ.0) GO TO 100
C
C           Free-free contribution
C
            IE=IL
            IF(IE.EQ.IELHM) GO TO 65
            CH=IZ(IL)*IZ(IL)
            SF1=CH*XN*SGFF/(FR*FR*FR)
C
C           The following expression is the so-called modified free-free
C           opacity, ie. allowing for the photoionization from higher,
C           non-explicit, LTE energy levels of the ion IL
C
            IF(IT.NE.2) GO TO 50
            SG=GFREE(T,FR/CH)
            SF2=SF2+SG-UN
   50       SFF=SF1
            GO TO 70
   65       SFF=SFFHMI(XN,FR,T)
   70       AFF=AFF+SFF
  100    CONTINUE
C
C        Additional opacities
C
         CALL OPADD(0,ID,FR,ABAD,EMAD,SCAD)
         IF(IOPHLI.NE.0) CALL LYMLIN(ID,FR,ABLY,EMLY,SCLY)
C
C        Total opacity and emissivity
C
         X=EXP(-HKF)
         X1=UN-X
         BNE=BNU*X*ANE
         ABSOC(IJ)=ABF+ANE*(X1*AFF-EBF)+ABAD+ABLY
         EMISC(IJ)=BNE*AFF+BNU*ANE*EBF+EMAD+EMLY
         SCATC(IJ)=SCAD+SCLY+sce
c            if(id.eq.1.or.id.eq.50) write(*,*)'opacon-tot',id,ij,
c    *   abf,ane,absoc(ij)

  200 CONTINUE
C
      CALL PHTION(ID,ABSOC,EMISC,FREQC,NFREQC)
      CALL PHTX(ID,ABSOC,EMISC,FREQC,1)
C
      RETURN
      END
C
C
C ********************************************************************
C
C
      FUNCTION SGMERG(II,ID,FR)
C     =========================
C     formal routine - taken from TLUSTY, but not used here
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (FRH=3.28805E15, PH2=2.815D29*2., EHB=157802.77355)
C
      sgmerg=0.
c     if(id.gt.0) return
      IE=IEL(II)
      CH=IZ(IE)*IZ(IE)
      g(ii)=gmer(imrg(ii),id)
      T1=1./TEMP(ID)
      EX=EHB*CH*T1
      II0=NQUANT(II-1)+1
      SUM=0.
      SUD=0.
      DO 10 I=II0,NLMX
         X=I
         XI=1./(X*X)
         FREDG=FRH*CH*XI
         IF(FR.LT.FREDG) GO TO 10
         EXI=EXP(EX*XI)
          S=EXI*WNHINT(I,ID)*XI/X
          SUM=SUM+S
c         SUD=SUD+S*XI
   10 CONTINUE
      SG0=PH2/(FR*FR*FR*G(II))*CH*CH
      SGMERG=SUM*SG0
c     DSG=-SUD*SG0*EX*T1
      RETURN
      END
C
C
C     ****************************************************************
C
      FUNCTION GFREE(T,FR)
C     ====================
C
C     Hydrogenic free-free Gaunt factor, for temperature T and
C     frequency FR
C
      INCLUDE 'PARAMS.FOR'
      THET=5040.4/T
      IF(THET.LT.4.E-2) THET=4.E-2
      X=FR/2.99793E14
      IF(X.GT.1) GO TO 10
      IF(X.LT.0.2) X=0.2
      GFREE=(1.0823+2.98E-2/THET)+(6.7E-3+1.12E-2/THET)/X
      RETURN
   10 C1=(3.9999187E-3-7.8622889E-5/THET)/THET+1.070192
      C2=(6.4628601E-2-6.1953813E-4/THET)/THET+2.6061249E-1
      C3=(1.3983474E-5/THET+3.7542343E-2)/THET+5.7917786E-1
      C4=3.4169006E-1+1.1852264E-2/THET
      GFREE=((C4/X-C3)/X+C2)/X+C1
      RETURN
      END
C
C ********************************************************************
C ********************************************************************
C
      FUNCTION SFFHMI_old(POPI,FR,T)
C     ==========================
C
C     Free-free cross section for H- (After Kurucz,1970,SAO 309, P.80)
C
      INCLUDE 'PARAMS.FOR'
      SFFHMI_old=(1.3727E-25+(4.3748E-10-2.5993E-7/T)/FR)*POPI/FR
      RETURN
      END
C
C 
C ********************************************************************
C
C
      SUBROUTINE LYMLIN(ID,FREQ,ABLY,EMLY,SCLY)
C     =========================================
C
C     OPACITY OF THE LYMAN LINES WINGS (ALPHA - DELTA)
C     WITH APPROXIMATE PARTIAL REDISTRIBUTION
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION SN(4),SR(4),SS(4),GS(4),FRLY(4),BNLY(4),GA(4)
      DATA FRLY / 2.4660375E15, 2.9227111E15, 3.0825469E15, 3.156528E15/
     *    ,BNLY / 5.527E-2,     4.090E-2,     2.699E-2,     1.855E-2  /,
     *     SN   / 1.308E5,      5.280E3,      5.847E2,      1.078E2   /,
     *     SR   / 1.218E-16,    9.196E-17,    1.058E-16,    1.296E-16 /,
     *     SS   / 9.478E-3,     1.600E-2,     1.441E-2,     1.547E-2  /,
     *     GS   / 7.237E-8,     5.432E-6,     5.821E-5,     4.027E-4  /,
     *     GA   / 1.000,        1.791,        2.362,        2.801 /
C
      data icomp/0/
      if(iath.le.0) return
      if(icomp.eq.0) then
         icomp=1
         read(4,*,err=10,end=10) ifstrk,ifnat,ifres,ifprd,ifsti
         go to 11
   10    continue
         ifstrk=0
         ifnat=1
         ifres=1
         ifprd=0
         ifsti=0
         if(iophli.lt.0) then
            ifstrk=1
            ifprd=1
         end if
   11    continue
      end if
c
      ABLY=0.
      EMLY=0.
      SCLY=0.
      
      if(freq.gt.3.3e15) return
      
      P=POPUL(N0HN,ID)
      T=TEMP(ID)
      ANE=ELEC(ID)
      DO 40 I=1,4
         DFR=ABS(FRLY(I)-FREQ)
         IF(DFR.LE.5.E11) DFR=1.E12
         DFR2=DFR*DFR
         DFRS=SQRT(DFR)
         COR=(2.*FREQ/(FREQ+FRLY(I)))**2
         F=1.
         IF(iabs(IOPHLI).EQ.2) F=FEAUTR(FREQ,ID)
         STARK=SS(I)*ANE*F/DFR2/DFRS
         if(ifstrk.eq.0) stark=0.
         if(ifnat.eq.0) sn(i)=0.
         if(ifres.eq.0) sr(i)=0.
         SGLY=SN(I)*(1.+SR(I)*P)*COR/DFR2+STARK
         sgly=sgly*wnhint(i+1,id)
         GAMA=1./(GA(I)+GS(I)*ANE*F/DFRS)
         if(ifprd.eq.0) gama=0.
         ABLY=ABLY+P*SGLY
         EMLY=EMLY+POPUL(N0HN+I,ID)*SGLY*BNLY(I)*(1.-GAMA)
         if(ifsti.ne.0) ably=ably-popul(n0hn+i,id)*sgly/(i+1)/(i+1)
         SCLY=SCLY+P*SGLY*GAMA
   40 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
      FUNCTION FEAUTR(FREQ,ID)
C     ========================
C
C     LYMAN-ALPHA STARK BROADENING AFTER N.FEAUTRIER
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION DL(20),F05(20),F10(20),F20(20),F40(20),X(4)
      DATA F05 / 0.0537, 0.0964, 0.1330, 0.3105, 0.4585, 0.6772, 0.8229,
     *           0.8556, 0.9250, 0.9618, 0.9733, 1.1076, 1.0644, 1.0525,
     *           0.8841, 0.8282, 0.7541, 0.7091, 0.7164, 0.7672/
      DATA F10 / 0.1986, 0.2764, 0.3959, 0.5740, 0.7385, 0.9448, 1.0292,
     *           1.0317, 0.9947, 0.8679, 0.8648, 0.9815, 1.0660, 1.0793,
     *           1.0699, 1.0357, 0.9245, 0.8603, 0.8195, 0.7928/
      DATA F20 / 0.4843, 0.5821, 0.7003, 0.8411, 0.9405, 1.0300, 1.0029,
     *           0.9753, 0.8478, 0.6851, 0.6861, 0.8554, 0.9916, 1.0264,
     *           1.0592, 1.0817, 1.0575, 1.0152, 0.9761, 0.9451/
      DATA F40 / 0.7862, 0.8566, 0.9290, 0.9915, 1.0066, 0.9878, 0.8983,
     *           0.8513, 0.6881, 0.5277, 0.5302, 0.6920, 0.8607, 0.9111,
     *           0.9651, 1.0793, 1.1108, 1.1156, 1.1003, 1.0839/
      DATA DL / -150., -120., -90., -60., -40., -20., -10., -8., -4.,
     *          -2., 2., 4., 8., 10., 20., 40., 60., 90., 120., 150./
      DLAM=2.997925E18/FREQ-1215.685
      DO 10 I=2,20
         IF(DLAM.LE.DL(I)) GO TO 20
   10 CONTINUE
      I=20
   20 J=I-1
      C=DL(J)-DL(I)
      A=(DLAM-DL(I))/C
      B=(DL(J)-DLAM)/C
      X(1)=F05(J)*A+F05(I)*B
      X(2)=F10(J)*A+F10(I)*B
      X(3)=F20(J)*A+F20(I)*B
      X(4)=F40(J)*A+F40(I)*B
      J=JT(ID)
      Y=TI0(ID)*X(J)+TI1(ID)*X(J-1)+TI2(ID)*X(J-2)
      FEAUTR=0.5*(Y+1.)
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE HYLSET
C     =================
C
C     Initialization procedure for treating the hydrogen line opacity
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      DIMENSION ALB(15)
      DATA ALB /656.28,486.13,434.05,410.17,397.01,
     *          388.91,383.54,379.79,377.06,375.02,
     *          373.44,372.19,371.20,370.39,369.72/
C
C     IHYL=-1  -  hydrogen lines are excluded a priori
C
      IHYL=-1
      if(iath.le.0) return
      IF(FREQ(2).GE.3.28805E15) RETURN
      AL0=2.997925E17/FREQ(1)
      AL1=2.997925E17/FREQ(2)
c     IF(grav.lt.6.) then
c        IF(AL0.GT.160..AND.AL1.LT.364.6) RETURN
c        IF(AL0.GT.506..AND.AL1.LT.630.) RETURN
c        IF(AL0.GT.680..AND.AL1.LT.820.3) RETURN
C         IF(AL0.GT.1500.) RETURN
c      else
          IF(AL0.GT.200..AND.AL1.LT.364.6) RETURN
          IF(AL0.GT.560..AND.AL1.LT.580.) RETURN
          IF(AL0.GT.720..AND.AL1.LT.820.3) RETURN
C         IF(AL0.GT.1500.) RETURN
c     end if
C
C     otherwise, hydrogen lines are included
C
      IHYL=0
      M20=40
      IF(AL1.LT.364.6) THEN
         ILOWH=1
         FRION=3.28805E15
         M10=int(SQRT(3.28805E15/ABS(FRION-FREQ(2))))
         IF(FRION.GT.FREQ(1)) M20=int(SQRT(3.28805E15/(FRION-FREQ(1))))
         IHYL=1
         IF(AL0.GT.123.) IHYL=0
         IF(AL0.GT.104..AND.AL1.LT.120.) IHYL=0
         IF(AL0.GT.98.5.AND.AL1.LT.102.) IHYL=0
         IF(IMODE.EQ.2.OR.IHYDPR.NE.0.OR.GRAV.GE.6.) IHYL=1
       ELSE IF(AL1.LT.820.) THEN
         ILOWH=2
         if(vaclim.lt.3600.) then
         FRION=8.2225E14
         M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2))))
         else
         FRION=8.22013E14
         M10=int(SQRT(3.28805E15/ABS(FRION-FREQ(2))))
         end if
         IF(FRION.GT.FREQ(1)) M20=int(SQRT(3.289017E15/(FRION-FREQ(1)))) 
         DO 10 I=1,15
            AL=ALB(I)
            IF(AL.LT.AL0-1..OR.AL.GT.AL1+1.) GO TO 10
            IHYL=1
            GO TO 20
   10    CONTINUE
   20    CONTINUE
         IF(IMODE.EQ.2.OR.IHYDPR.NE.0.OR.GRAV.GE.6.) IHYL=1
       ELSE
         ILOWH=3
         IHYL=1
      END IF
c      ELSE IF(AL1.LT.1458.) THEN
c        ILOWH=3
c        FRION=3.6544142E14
c        M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2))))
c        IF(FRION.GT.FREQ(1)) M20=int(SQRT(3.289017E15/(FRION-FREQ(1))))
c        IHYL=1
c        IF(AL0.GT.1310.) IHYL=0
c        IF(AL0.GT.1124..AND.AL1.LT.1250.) IHYL=0
c        IF(AL0.GT.1035..AND.AL1.LT.1060.) IHYL=0
c        IF(IMODE.EQ.2.OR.GRAV.GE.6.) IHYL=1
c      ELSE IF(AL1.LT.2278.) THEN
c        ILOWH=4
c        FRION=2.0555837E14
c        M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2))))
c        IHYL=1
c      ELSE IF(AL1.LT.3281.) THEN
c        ILOWH=5
c        FRION=1.315589E14
c        M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2))))
c        IHYL=1
c      ELSE IF(AL1.LT.4466.) THEN
c        ILOWH=6
c        FRION=9.136394E13
c        M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2))))
c        IHYL=1
c      ELSE 
c        ILOWH=7
c        FRION=6.7120228E13
c        M10=int(SQRT(3.289017E15/ABS(FRION-FREQ(2))))
c        IHYL=1
c     END IF
c     WRITE(6,601) al0,al1,ILOWH,M20+1
c 601 FORMAT(/ ' *** HYDROGEN LINES CONTRIBUTE',
c 601 FORMAT(1H0/ ' *** HYDROGEN LINES CONTRIBUTE'/
c    * '     THE NEAREST LINE ON THE SHORT-WAVELENGTH SIDE IS',
c    * 2f10.1,3x,I3,'  TO ',I3)
c
      ihyl=1
c
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE HYLSEW(IJ)
C     =====================
C
C     Initialization procedure for treating the hydrogen line opacity
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
C
C     IHYL=-1  -  hydrogen lines are excluded a priori
C
      IHYLW(IJ)=0
      if(iath.le.0) return
      FR=FREQ(IJ)
      IF(FR.GE.3.28805E15) RETURN
      AL0=2.997925E17/FR
      AL1=AL0
      IF(grav.lt.6.) then
         IF(AL0.GT.160..AND.AL1.LT.364.6) RETURN
         IF(AL0.GT.506..AND.AL1.LT.630.) RETURN
         IF(AL0.GT.680..AND.AL1.LT.820.3) RETURN
       else
         IF(AL0.GT.540..AND.AL1.LT.600.) RETURN
         IF(AL0.GT.720..AND.AL1.LT.820.3) RETURN
      end if
C
C     otherwise, hydrogen lines are included
C
      IHYLW(IJ)=1
      M20W(IJ)=40
      IF(AL1.LT.364.6) THEN
         ILOWHW(IJ)=1
         FRION=3.28805E15
       ELSE IF(AL1.LT.820.) THEN
         ILOWHW(IJ)=2
         FRION=8.2225E14
       ELSE IF(AL1.LT.1458.) THEN
         ILOWHW(IJ)=3
         FRION=3.6544142E14
       ELSE IF(AL1.LT.2278.) THEN
         ILOWHW(IJ)=4
         FRION=2.0555837E14
       ELSE IF(AL1.LT.3281.) THEN
         ILOWHW(IJ)=5
         FRION=1.315589E14
       ELSE IF(AL1.LT.4466.) THEN
         ILOWHW(IJ)=6
         FRION=9.136394E13
       ELSE 
         ILOWHW(IJ)=7
         FRION=6.7120228E13
      END IF
      IF(FRION.GT.FR) M10W(IJ)=int(SQRT(3.289017E15/ABS(FRION-FR)))
c     WRITE(6,601) ILOWH,M20+1
c 601 FORMAT(1H0/ ' *** HYDROGEN LINES CONTRIBUTE'/
c    * '     THE NEAREST LINE ON THE SHORT-WAVELENGTH SIDE IS',
c    * I3,'  TO ',I3/)
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE HYDLIN(ID,I0,I1,ABSOH,EMISH)
C     =======================================
C
C     opacity and emissivity of hydrogen lines  
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (FRH1=3.28805E15,FRH2=FRH1/4.,UN=1.,SIXTH=1./6.)
      PARAMETER (CPP=4.1412E-16,CPJ=157803.) 
      PARAMETER (C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.)
      PARAMETER (CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18)
      PARAMETER (CID1=0.01497)
      common/quasun/nunalp,nunbet,nungam,nunbal
      common/hhebrd/sthe,nunhhe
      common/gompar/hglim,ihgom
      DIMENSION PJ(40),PRF0(54),OSCH(4,22),
     *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
      dimension wlir(15),irlow(15),irupp(15)
      DATA FRH    /3.289017E15/
      data wlir/
     *    123680., 75005., 59066., 51273.,190570.,113060.,
     *     87577., 75061.,277960.,162050.,123840.,105010.,
     *    223340.,168760.,141790./
      data irlow/4*6, 4*7, 4*8, 3*9/
      data irupp/7,8,9,10,8,9,10,11,9,10,11,12,11,12,13/
      data nlinir/15/
c
      DATA INIT /0/
C
      DO IJ=I0,I1
         ABSOH(IJ)=0.
         EMISH(IJ)=0.
      END DO
c
      if(iath.le.0.or.rrr(1,1,1).eq.0.) return
      izz=1
C
      IF(INIT.EQ.0) THEN
         DO I=1,4
            DO J=I+1,22
               CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0)
               WLINE(I,J)=WL0
               OSCH(I,J)=FIJ+FIJ0
            END DO
         END DO
         INIT=1
      END IF
      DO IJ=I0,I1
         ABSO(IJ)=0.
         EMIS(IJ)=0.
      END DO
c
       if(ilowh.le.0) return
c
       T=TEMP(ID)
       T1=UN/T
       SQT=SQRT(T)
       ANE=ELEC(ID)
       ANES=EXP(SIXTH*LOG(ANE))
       TL=LOG10(T)
       ANEL=LOG10(ANE)
C
C     populations of the first 40 levels of hydrogen
C
      ANP=POPUL(NKH,ID)
      PP=CPP*ANE*ANP*T1/SQT
      NLH=N1H-N0HN+1
c     if(ifwop(n1h).lt.0) nlh=nlh-1
      nlh=nlh-1
      DO 5 IL=1,50
         X=IL*IL
         IF(IL.LE.NLH) PJ(IL)=POPUL(N0HN+IL-1,ID)
         IF(IL.GT.NLH) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhint(il,id)
    5 CONTINUE
      p2=pp*exp(cpj4*t1)*4.*wnhint(2,id)
c
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      F00=C00*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID))
C
C -------------------------------------------------------------------
C     overall loop over spectral series (only in the infrared region)
C -------------------------------------------------------------------
C
      ISERL=ILOWH
      ISERU=ILOWH
c
      if(wlam(i0).gt.14000.) iseru=4
      if(wlam(i0).gt.22700.) iseru=5
      if(wlam(i0).gt.32800.) iseru=6
      if(wlam(i0).gt.44660.) iseru=7
      if(wlam(i0).gt.60000.) iserl=4

c
      if(iserl.eq.3.and.iseru.eq.3.and.nunbal.gt.0) iserl=2
      DO IJ=I0,I1
         ABSO(IJ)=0.
         EMIS(IJ)=0.
      END DO   
C
c     ========================
c     loop over specral series
c     ========================
c
      DO 200 I=ISERL,ISERU
c
c        skip the following calculations if one uses the Gomez tables
c
         if(ihgom.gt.0.and.elec(id).gt.hglim) then
            if(i.ge.1.and.i.le.ihgom) then
               call ghydop(id,i0,i1,pj,absoh,emish)
               go to 200
            end if
         end if
c
      II=I*I
      XII=UN/II
      POPI=PJ(I)
      IF(I.EQ.1) FRH=3.28805E15
C
C     determination of which hydrogen lines contribute in a current 
C     frequency region
C
      M1=M10
      IF(I.LT.ILOWH) M1=ILOWH-1
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.3..and.M1.LE.6.AND.I.EQ.2) GO TO 10
      IF(grav.lt.3..and.M1.LE.4.AND.I.EQ.1) GO TO 10
      M1=M1-1
      M2=M20+3
      IF(M1.LT.I+1) M1=I+1
   10 CONTINUE
      if(grav.gt.3.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3 
      end if
c  new! 
         if(i.ge.3) then
            m1=i+1
            m2=i+40
         end if
         if(i.ge.4) m2=i+20
         if(i.ge.6) m2=i+10
C
C     loop over lines which contribute at given wavelength region
C
      m1=min(m1,40)
      m2=min(m2,40)
      m1=max(m1,i+1)
      m2=max(m2,i+2)
      DO 100 J=M1,M2
         IF(I.EQ.1.AND.J.LE.5.AND.IOPHLI.LT.0) GO TO 100
         ILINE=0
         JJ=J*J
         XJJ=UN/JJ
         ABTRA=PJ(I)*WNHINT(J,ID)
         EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1)
         if(i.le.2.and.j.le.i+2) then
            abtra=pj(i)
            emtra=pj(j)*wnhint(i,id)/wnhint(j,id)*
     *            ii*xjj*exp(cpj*(xii-xjj)*t1)
         end if
         IF(I.LE.4.AND.J.LE.22) ILINE=ILIN0(I,J)
c
c        quasi-molecular opacity for Lyman-alpha and beta satellites
c
         lquasi=i.eq.1.and.j.eq.2.and.nunalp.gt.0
         lquasi=lquasi.or.i.eq.1.and.j.eq.3.and.nunbet.gt.0
         lquasi=lquasi.or.i.eq.1.and.j.eq.4.and.nungam.gt.0
         lquasi=lquasi.or.i.eq.2.and.j.eq.3.and.nunbal.gt.0
         lalhhe=i.eq.1.and.j.eq.2.and.nunhhe.gt.0
         if(lquasi) then
            DO IJ=I0,I1
               call allard(wlam(ij),popi,anp,sg,i,j)
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            END DO
         end if
         ahe=0.
         if(iathe.gt.0) ahe=popul(n0a(iathe),id)
         if(lalhhe.and.ahe.gt.0.) then
            rel=1./6.2831855
            do ij=i0,i1
               call lyahhe(wlam(ij),ahe,sg0)
               sg=sg0*rel
               abso(ij)=abso(ij)+sg*abtra
               emis(ij)=emis(ij)+sg*emtra
            end do
         end if
c
c        lines with special Stark broadening tables
c
         IF(ILINE.GT.0) THEN
          FID=CID*OSCH(I,J)
c
c         switch to either original Lemke/Tremblay of Xenomorph
c
          if(ilxen(i,j).eq.0.or.anel.lt.xnemin) then
c
c         original Lemke/Tremblay
c
            NWL=NWLHYD(ILINE)
            DO  IWL=1,NWL
               PRF0(IWL)=PRFHYD(ILINE,ID,IWL)
            END DO
            DO IJ=I0,I1
               AL=ABS(WLAM(IJ)-WLINE(I,J))
               IF(AL.LT.1.E-4) AL=1.E-4
               IF(ILEMKE.EQ.1) AL=AL/F00
               AL=LOG10(AL)
               DO 30 IWL=1,NWL-1
                  IW0=IWL
                  IF(AL.LE.WLHYD(ILINE,IWL+1)) GO TO 40
   30          CONTINUE
   40          IW1=IW0+1
               PRFF=(PRF0(IW0)*(WLHYD(ILINE,IW1)-AL)+PRF0(IW1)*
     *             (AL-WLHYD(ILINE,IW0)))/
     *             (WLHYD(ILINE,IW1)-WLHYD(ILINE,IW0))
               SG=EXP(PRFF*AL10)*FID
               sg0=EXP(PRFF*AL10)
               IF(ILEMKE.EQ.1) SG=SG*WLINE(I,J)**2*CINV/F00
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            END DO
c
c         XENOMORPH data for selected lines
c
          else
             ixn=ilxen(i,j)
             nwl=nwlxen(ixn)
             fr0l=2.997925e18/wline(i,j)
             do ij=i0,i1
                al=(freq(ij)-fr0l)/f00
                if(abs(al).lt.1.e-4) al=1.e-4
                all=log10(abs(al))
                do 51 iwl=1,nwl-1
                   iw0=iwl
                   if(all.le.alxen(ixn,iwl+1)) go to 52
   51           continue
   52           iw1=iw0+1
                if(al.gt.0.) then
                   prff=(prfb(ixn,id,iw0)*(alxen(ixn,iw1)-all)+
     *                   prfb(ixn,id,iw1)*(all-alxen(ixn,iw0)))/
     *                   (alxen(ixn,iw1)-alxen(ixn,iw0))
                 else
                   prff=(prfr(ixn,id,iw0)*(alxen(ixn,iw1)-all)+
     *                   prfr(ixn,id,iw1)*(all-alxen(ixn,iw0)))/
     *                   (alxen(ixn,iw1)-alxen(ixn,iw0))
                end if
                sg=exp(prff*al10)*fid/f00
                ABSO(IJ)=ABSO(IJ)+SG*ABTRA
                EMIS(IJ)=EMIS(IJ)+SG*EMTRA
             end do
          END IF

c
c         lines without special Stark broadening tables
c
          ELSE
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            if((wl0.le.wlam(i1).and.1.25*wl0.gt.wlam(i0)). or.
     *         (wl0.ge.wlam(i0).and.0.75*wl0.lt.wlam(i1))) then
            FXK=F00*XKIJ
            FXK1=UN/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0*CINV*FXK1
            BETAD=DOP*DBETA
            FID=CID*FIJ*DBETA
c           FID0=CID1*FIJ0/DOP
            CALL DIVSTR(AD,DIV)
            fac=two
            if(lquasi) fac=un
            DO IJ=I0,I1
               fr=freq(ij)
               BETA=ABS(WLAM(IJ)-WL0)*FXK1
               IF(I.LT.5) THEN
                  SG=STARKA(BETA,AD,DIV,fac)*FID
                  if(iophli.eq.2.and.i.eq.1.and.j.eq.2) 
     *               sg=sg*feautr(fr,id)
                ELSE
                  SG=STARKIR(II,JJ,T,ANE,BETA)*FID
               END IF
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            END DO
            END IF
         END IF
  100 CONTINUE
  200 CONTINUE
C
C     far infrared hydrogen lines
C
      if(wlam(i1).gt.70000.) then
      DO I=8,13
         II=I*I
         XII=UN/II
         DO J=I+1,I+4
            JJ=J*J
            XJJ=UN/JJ
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            if((wl0.le.wlam(i1).and.1.5*wl0.gt.wlam(i0)). or.
     *         (wl0.ge.wlam(i0).and.0.5*wl0.lt.wlam(i1))) then
            FXK=F00*XKIJ
            FXK1=UN/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0*CINV*FXK1
            BETAD=DOP*DBETA
            FID=CID*FIJ*DBETA
            CALL DIVSTR(AD,DIV)
            fac=two
            DO IJ=I0,I1
               fr=freq(ij)
               BETA=ABS(WLAM(IJ)-WL0)*FXK1
               SG=STARKIR(II,JJ,T,ANE,BETA)*FID
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            END DO
            END IF
         END DO
      END DO
      END IF
c
      if(wlam(i1).gt.5.e5) then
      do ij=i0,i1
         fr=freq(ij)
         do ilir=1,nlinir
            if(wlam(ij).gt.wlir(ilir)*0.95.and.
     *         wlam(ij).lt.wlir(ilir)*1.05) then
               j=irupp(ilir)
               JJ=J*J
               i=irlow(ilir)
               II=I*I
               XII=UN/II
               XJJ=UN/JJ
               ABTRA=PJ(I)*WNHINT(J,ID)
               EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1)
               CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
               FXK=F00*XKIJ
               FXK1=UN/FXK
               DOP=DOP0/WL0
               DBETA=WL0*WL0*CINV*FXK1
               BETAD=DOP*DBETA
               FID=CID*FIJ*DBETA
               CALL DIVSTR(AD,DIV)
               fac=two
               BETA=ABS(WLAM(IJ)-WL0)*FXK1
               SG=STARKA(BETA,AD,DIV,fac)*FID
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            end if
         end do
      end do
      end if
C
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      DO IJ=I0,I1
         F=FREQ(IJ)
         F15=F*1.E-15
         XKF=EXP(-4.79928e-11*F*T1)
         XKFB=XKF*1.4743E-2*F15*F15*F15
         ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ)
         EMISH(IJ)=XKFB*EMIS(IJ)
      END DO
      RETURN
      END
C
C
C ********************************************************************
C
      SUBROUTINE HYDLIW(ID,ABSOH,EMISH)
C     =================================
C
C     opacity and emissivity of hydrogen lines  
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (FRH1=3.28805E15,FRH2=FRH1/4.,UN=1.,SIXTH=1./6.)
      PARAMETER (CPP=4.1412E-16,CPJ=157803.) 
      PARAMETER (C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.)
      PARAMETER (CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18)
      PARAMETER (CID1=0.01497)
      common/lasers/lasdel
      common/quasun/nunalp,nunbet,nungam,nunbal
      DIMENSION PJ(40),PRF0(54),OSCH(4,22),
     *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
      DATA FRH    /3.289017E15/
      DATA INIT /0/
C
      if(iath.le.0) return
      izz=1
C
      IF(INIT.EQ.0) THEN
         DO I=1,4
            DO J=I+1,22
               CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0)
               WLINE(I,J)=WL0
               OSCH(I,J)=FIJ+FIJ0
            END DO
         END DO
         INIT=1
      END IF
      DO IJ=1,NFREQ
         ABSO(IJ)=0.
         EMIS(IJ)=0.
         ABSOH(IJ)=0.
         EMISH(IJ)=0.
      END DO
       T=TEMP(ID)
       T1=UN/T
       SQT=SQRT(T)
       ANE=ELEC(ID)
       ANES=EXP(SIXTH*LOG(ANE))
C
C     populations of the first 40 levels of hydrogen
C
      ANP=POPUL(NKH,ID)
      PP=CPP*ANE*ANP*T1/SQT
      NLH=N1H-N0HN+1
      if(ifwop(n1h).lt.0) nlh=nlh-1
      DO 5 IL=1,40
         X=IL*IL
         IF(IL.LE.NLH) PJ(IL)=POPUL(N0HN+IL-1,ID)
         IF(IL.GT.NLH) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhint(il,id)
    5 CONTINUE
      p2=pp*exp(cpj4*t1)*4.*wnhint(2,id)
C
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      F00=C00*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(1.65E8*T+VTURB(ID))
C
C -------------------------------------------------------------------
C     overall loop over spectral series (only in the infrared region)
C -------------------------------------------------------------------
C
      DO 300 IJ=1,NFREQ
      IF(IHYLW(IJ).LE.0) GO TO 300
      ISERL=ILOWHW(IJ)
      ISERU=ILOWHW(IJ)
      IF(WLAM(IJ).GT.17000..AND.WLAM(IJ).LE.21000.) THEN
         ISERL=3
         ISERU=4
       ELSE IF(WLAM(IJ).GT.22700..AND.WLAM(IJ).LE.29000.) THEN
         ISERL=4
         ISERU=5
       ELSE IF(WLAM(IJ).GT.32800..AND.WLAM(IJ).LE.37000.) THEN
         ISERL=5
         ISERU=6
       ELSE IF(WLAM(IJ).GT.37000..AND.WLAM(IJ).LE.44600.) THEN
         ISERL=4
         ISERU=6
       ELSE IF(WLAM(IJ).GT.44660..AND.WLAM(IJ).LE.58300.) THEN
         ISERL=5
         ISERU=7
       ELSE IF(WLAM(IJ).GT.58300..AND.WLAM(IJ).LE.72000.) THEN
         ISERL=6
         ISERU=8
       ELSE IF(WLAM(IJ).GT.72000..AND.WLAM(IJ).LE.73800.) THEN
         ISERL=5
         ISERU=8
       ELSE IF(WLAM(IJ).GT.73800..AND.WLAM(IJ).LE.77000.) THEN
         ISERL=5
         ISERU=9
       ELSE IF(WLAM(IJ).GT.77000.) THEN
         ISERL=6
         ISERU=9
      END IF
C
      if(iserl.eq.3.and.iseru.eq.3.and.nunbal.gt.0) iserl=2
C
      ABSO(IJ)=0.
      EMIS(IJ)=0.
      DO 200 I=ISERL,ISERU
      II=I*I
      XII=UN/II
      PLTEI=PP*EXP(CPJ*T1*XII)*II
      POPI=PJ(I)
      IF(I.EQ.1) FRH=3.28805E15
C
C     determination of which hydrogen lines contribute in a current 
C     frequency region
C
      M1=M10W(IJ)
      IF(I.LT.ILOWHW(IJ)) M1=ILOWHW(IJ)-1
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.3..and.M1.LE.16.AND.I.EQ.7) GO TO 10
      IF(grav.lt.3..and.M1.LE.14.AND.I.EQ.6) GO TO 10
      IF(grav.lt.3..and.M1.LE.12.AND.I.EQ.5) GO TO 10
      IF(grav.lt.3..and.M1.LE.10.AND.I.EQ.4) GO TO 10
      IF(grav.lt.3..and.M1.LE.8.AND.I.EQ.3) GO TO 10
      IF(grav.lt.3..and.M1.LE.6.AND.I.EQ.2) GO TO 10
      IF(grav.lt.3..and.M1.LE.4.AND.I.EQ.1) GO TO 10
      M1=M1-1
      M2=M20W(IJ)+3
      IF(M1.LT.I+1) M1=I+1
   10 CONTINUE
      if(grav.gt.3.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3 
      end if
      if(grav.gt.6.) then
         m2=m2+2
         m1=m1-1
         if(m1.gt.i+6) m1=m1-1 
      end if
      IF(M1.LT.I+1) M1=I+1
c      if(m2.gt.30) then
c        m2=m20W(IJ)+8
c         m1=m1-4
c      end if
      IF(M2.GT.40) M2=40
c     if(id.eq.1) write(6,666) i,m1,m2
c 666 format(/' hydrogen lines contribute - ilow=',i2,', iup from ',i3,
c    *       ' to',i3/)
C
      A=0.
      E=0.
C
C     loop over lines which contribute at given wavelength region
C
      DO 100 J=M1,M2
         IF(I.EQ.1.AND.J.LE.5.AND.IOPHLI.LT.0) GO TO 100
         ILINE=0
         JJ=J*J
         XJJ=UN/JJ
         ABTRA=PJ(I)*WNHINT(J,ID)
         EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1)
         if(i.le.2.and.j.le.i+2) then
            abtra=pj(i)
            emtra=pj(j)*wnhint(i,id)/wnhint(j,id)*
     *            ii*xjj*exp(cpj*(xii-xjj)*t1)
         end if
         IF(I.LE.4.AND.J.LE.22) ILINE=ILIN0(I,J)
c
c        quasi-molecular opacity for Lyman-alpha and beta satellites
c
         lquasi=i.eq.1.and.j.eq.2.and.nunalp.gt.0
         lquasi=lquasi.or.i.eq.1.and.j.eq.3.and.nunbet.gt.0
         lquasi=lquasi.or.i.eq.1.and.j.eq.4.and.nungam.gt.0
         lquasi=lquasi.or.i.eq.2.and.j.eq.3.and.nunbal.gt.0
         if(lquasi) then
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            FXK=F00*XKIJ
            FXK1=UN/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0*CINV*FXK1
            BETAD=DOP*DBETA
            FID=CID*FIJ*DBETA
            CALL DIVSTR(AD,DIV)
            fr=freq(ij)
            BETA=ABS(WLAM(IJ)-WL0)*FXK1
            call allard(wlam(ij),popi,anp,sg,i,j)
            sg=sg+STARKA(BETA,AD,DIV,UN)*FID
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            go to 100
         end if
c
c        lines with special Stark broadening tables
c
         IF(ILINE.GT.0) THEN
            NWL=NWLHYD(ILINE)
            DO IWL=1,NWL
               PRF0(IWL)=PRFHYD(ILINE,ID,IWL)
            END DO
            FID=CID*OSCH(I,J)
            AL=ABS(WLAM(IJ)-WLINE(I,J))
            IF(AL.LT.1.E-4) AL=1.E-4
            IF(ILEMKE.EQ.1) AL=AL/F00
            AL=LOG10(AL)
            DO 30 IWL=1,NWL-1
               IW0=IWL
               IF(AL.LE.WLHYD(ILINE,IWL+1)) GO TO 40
   30       CONTINUE
   40       IW1=IW0+1
            PRFF=(PRF0(IW0)*(WLHYD(ILINE,IW1)-AL)+PRF0(IW1)*
     *          (AL-WLHYD(ILINE,IW0)))/
     *          (WLHYD(ILINE,IW1)-WLHYD(ILINE,IW0))
            SG=EXP(PRFF*AL10)*FID
            IF(ILEMKE.EQ.1) SG=SG*WLINE(I,J)**2*CINV/F00
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA
c
c         lines without special Stark broadening tables
c
          ELSE
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            FXK=F00*XKIJ
            FXK1=UN/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0*CINV*FXK1
            BETAD=DOP*DBETA
            FID=CID*FIJ*DBETA
            CALL DIVSTR(AD,DIV)
            fr=freq(ij)
            BETA=ABS(WLAM(IJ)-WL0)*FXK1
            SG=STARKA(BETA,AD,DIV,TWO)*FID
            if(iophli.eq.2.and.i.eq.1.and.j.eq.2) 
     *            sg=sg*feautr(fr,id)
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA
         END IF
  100 CONTINUE
  200 CONTINUE
C
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      F=FREQ(IJ)
      F15=F*1.E-15
      XKF=EXP(-4.79928e-11*F*T1)
      XKFB=XKF*1.4743E-2*F15*F15*F15
      if(abso(ij).le.0. .and. lasdel) then
         abso(ij)=0.
         emis(ij)=0.
      endif
      ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ)
      EMISH(IJ)=XKFB*EMIS(IJ)
  300 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE HE2SET
C     =================
C
C     Initialization procedure for treating the He II line opacity
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      dimension frhe(12)
      DATA FRHE /1.3158153D+16, 3.2895381D+15, 1.4624854D+15,
     *           8.2261878D+14, 5.2647201D+14, 3.6560459D+14,
     *           2.6860713D+14, 2.0565220D+14, 1.6249055D+14,
     *           1.3161730D+14, 1.0877460D+14, 9.1400851D+13/
C
C     IHE2L=-1  -  He II lines are excluded a priori
C
      IHE2L=-1
      IF(IFHE2.LE.0) RETURN
      IF(FREQ(2).GE.1.315812E16) RETURN
      AL0=2.997925E17/FREQ(1)
      AL1=2.997925E17/FREQ(2)
c      IF(AL0.GT.390.) RETURN
      if(grav.lt.6.) then
         IF(AL0.GT.31..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.26.1.AND.AL1.LT.29.8) RETURN
         IF(AL0.GT.24.8.AND.AL1.LT.25.1) RETURN
         IF(AL0.GT.122.1.AND.AL1.LT.162.9) RETURN
         IF(AL0.GT.165.1.AND.AL1.LT.204.9) RETURN
         IF(AL0.GT.109..AND.AL1.LT.120.9) RETURN
         IF(AL0.GT.103..AND.AL1.LT.107.9) RETURN
         IF(AL0.GT.99.7.AND.AL1.LT.102.) RETURN
         IF(AL0.GT.320.8.AND.AL1.LT.364.4) RETURN
         IF(AL0.GT.273.8.AND.AL1.LT.319.8) RETURN
         IF(AL0.GT.251.6.AND.AL1.LT.272.8) RETURN
         IF(AL0.GT.239.0.AND.AL1.LT.250.6) RETURN
         IF(AL0.GT.231.1.AND.AL1.LT.238.0) RETURN
         IF(AL0.GT.225.8.AND.AL1.LT.230.1) RETURN
       else if(grav.lt.7.) then
         IF(AL0.GT.33..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.124.1.AND.AL1.LT.160.9) RETURN
         IF(AL0.GT.167.1.AND.AL1.LT.202.9) RETURN
         IF(AL0.GT.111..AND.AL1.LT.118.9) RETURN
         IF(AL0.GT.322.8.AND.AL1.LT.364.4) RETURN
         IF(AL0.GT.275.8.AND.AL1.LT.317.8) RETURN
         IF(AL0.GT.253.6.AND.AL1.LT.270.8) RETURN
         IF(AL0.GT.241.0.AND.AL1.LT.248.6) RETURN
         IF(AL0.GT.233.1.AND.AL1.LT.236.0) RETURN
       else 
         IF(AL0.GT.39..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.134.1.AND.AL1.LT.150.9) RETURN
         IF(AL0.GT.177.1.AND.AL1.LT.202.9) RETURN
      end if
C
C     otherwise, He II lines are included
C
      IHE2L=1
      MHE10=60
      MHE20=60
      IF(AL1.LT.91.) THEN
         ILWHE2=1
       ELSE IF(AL0.LT.204.) THEN
         ILWHE2=2
       ELSE IF(AL0.LT.364.) THEN
         ILWHE2=3
       ELSE IF(AL0.LT.569.) THEN
         ILWHE2=4
       ELSE IF(AL0.LT.819.) THEN
         ILWHE2=5
       ELSE IF(AL0.LT.1116.) THEN
         ILWHE2=6
       ELSE IF(AL0.LT.1457.) THEN
         ILWHE2=7
       ELSE IF(AL0.LT.1844.) THEN
         ILWHE2=8
       ELSE IF(AL0.LT.2277.) THEN
         ILWHE2=9
       ELSE IF(AL0.LT.2756.) THEN
         ILWHE2=10
       ELSE IF(AL0.LT.3279.) THEN
         ILWHE2=11
       ELSE
         ILWHE2=12
      END IF
      FRION=FRHE(ILWHE2)
      FR1=FRION*ILWHE2*ILWHE2
      IF(FRION.GT.FREQ(2)) MHE10=int(SQRT(FR1/(FRION-FREQ(2))))
      IF(FRION.GT.FREQ(1)) MHE20=int(SQRT(FR1/(FRION-FREQ(1))) )
      WRITE(6,601) ILWHE2,MHE20+1
  601 FORMAT(1H0/ ' *** HE II LINES CONTRIBUTE'/
     * '     THE NEAREST LINE ON THE SHORT-WAVELENGTH SIDE IS',
     * I3,'  TO ',I3/)
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE HE2SEW(IJ)
C     =====================
C
C     Initialization procedure for treating the He II line opacity
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      dimension frhe(12)
      DATA FRHE /1.3158153D+16, 3.2895381D+15, 1.4624854D+15,
     *           8.2261878D+14, 5.2647201D+14, 3.6560459D+14,
     *           2.6860713D+14, 2.0565220D+14, 1.6249055D+14,
     *           1.3161730D+14, 1.0877460D+14, 9.1400851D+13/
C
C     IHE2L=-1  -  He II lines are excluded a priori
C
      IHE2LW(IJ)=-1
      IF(IFHE2.LE.0) RETURN
      FR=FREQ(IJ)
      AL0=2.997925E17/FR
      AL1=2.997925E17/FR
      if(grav.lt.6.) then
         IF(AL0.GT.31..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.26.1.AND.AL1.LT.29.8) RETURN
         IF(AL0.GT.24.8.AND.AL1.LT.25.1) RETURN
         IF(AL0.GT.122.1.AND.AL1.LT.162.9) RETURN
         IF(AL0.GT.165.1.AND.AL1.LT.204.9) RETURN
         IF(AL0.GT.109..AND.AL1.LT.120.9) RETURN
         IF(AL0.GT.103..AND.AL1.LT.107.9) RETURN
         IF(AL0.GT.99.7.AND.AL1.LT.102.) RETURN
         IF(AL0.GT.320.8.AND.AL1.LT.364.4) RETURN
         IF(AL0.GT.273.8.AND.AL1.LT.319.8) RETURN
         IF(AL0.GT.251.6.AND.AL1.LT.272.8) RETURN
         IF(AL0.GT.239.0.AND.AL1.LT.250.6) RETURN
         IF(AL0.GT.231.1.AND.AL1.LT.238.0) RETURN
         IF(AL0.GT.225.8.AND.AL1.LT.230.1) RETURN
       else if(grav.lt.7.) then
         IF(AL0.GT.33..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.124.1.AND.AL1.LT.160.9) RETURN
         IF(AL0.GT.167.1.AND.AL1.LT.202.9) RETURN
         IF(AL0.GT.111..AND.AL1.LT.118.9) RETURN
         IF(AL0.GT.322.8.AND.AL1.LT.364.4) RETURN
         IF(AL0.GT.275.8.AND.AL1.LT.317.8) RETURN
         IF(AL0.GT.253.6.AND.AL1.LT.270.8) RETURN
         IF(AL0.GT.241.0.AND.AL1.LT.248.6) RETURN
         IF(AL0.GT.233.1.AND.AL1.LT.236.0) RETURN
       else 
         IF(AL0.GT.39..AND.AL1.LT.91.1) RETURN
         IF(AL0.GT.134.1.AND.AL1.LT.150.9) RETURN
         IF(AL0.GT.177.1.AND.AL1.LT.202.9) RETURN
      end if
C
C     otherwise, He II lines are included
C
      IHE2LW(IJ)=1
      MHE10W(IJ)=60
      MHE20W(IJ)=60
      IF(AL1.LT.91.) THEN
         ILWHEW(IJ)=1
       ELSE IF(AL0.LT.204.) THEN
         ILWHEW(IJ)=2
       ELSE IF(AL0.LT.364.) THEN
         ILWHEW(IJ)=3
       ELSE IF(AL0.LT.569.) THEN
         ILWHEW(IJ)=4
       ELSE IF(AL0.LT.819.) THEN
         ILWHEW(IJ)=5
       ELSE IF(AL0.LT.1116.) THEN
         ILWHEW(IJ)=6
       ELSE IF(AL0.LT.1457.) THEN
         ILWHEW(IJ)=7
       ELSE IF(AL0.LT.1844.) THEN
         ILWHEW(IJ)=8
       ELSE IF(AL0.LT.2277.) THEN
         ILWHEW(IJ)=9
       ELSE IF(AL0.LT.2756.) THEN
         ILWHEW(IJ)=10
       ELSE IF(AL0.LT.3279.) THEN
         ILWHEW(IJ)=11
       ELSE
         ILWHEW(IJ)=12
      END IF
      FRION=FRHE(ILWHEW(IJ))
      FR1=FRION*ILWHEW(IJ)*ILWHEW(IJ)
      IF(FRION.GT.FR) MHE10W(IJ)=int(SQRT(FR1/(FRION-FR)))
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE HE2LIN(ID,I0,I1,ABSOH,EMISH)
C
C     opacity and emissivity of He II lines  (these which are not considered
C     explicitly)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (UN=1.,SIXTH=1./6.)
      PARAMETER (CPP=4.1412E-16,CPJ=631479.) 
      PARAMETER (C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.)
      PARAMETER (CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18)
      PARAMETER (CID1=0.01497)
      DIMENSION PJ(80),FRHE(12),OSCHE2(19),PRF0(36),
     *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
      DATA FRHE /1.3158153D+16, 3.2895381D+15, 1.4624854D+15,
     *           8.2261878D+14, 5.2647201D+14, 3.6560459D+14,
     *           2.6860713D+14, 2.0565220D+14, 1.6249055D+14,
     *           1.3161730D+14, 1.0877460D+14, 9.1400851D+13/
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,
     *        1.604E-2, 1.023E-2, 6.980E-3,
     *        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,
     *        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,
     *        1.038,    1.793E-1, 6.549E-2/
C
      I=ILWHE2
      izz=2
      DO IJ=I0,I1
         ABSO(IJ)=0.
         EMIS(IJ)=0.
         ABSOH(IJ)=0.
         EMISH(IJ)=0.
      END DO
      T=TEMP(ID)
      T1=UN/T
      SQT=SQRT(T)
      ANE=ELEC(ID)
      ANES=EXP(SIXTH*LOG(ANE))
C
C     He III populations (either LTE or NLTE, depending on input model)
C
      IF(IELHE2.GT.0) THEN
         ANP=POPUL(NNEXT(IELHE2),ID)
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1
       ELSE
         ANP=RRR(ID,3,2)
         NLHE2=0
      END IF
C
C     populations of the first 60 levels of He II
C
      PP=CPP*ANE*ANP*T1/SQT
      DO IL=1,60
         X=IL*IL
         IIL=NFIRST(IELHE2)+IL-1
         IF(IL.LE.NLHE2) PJ(IL)=POPUL(IIL,ID)
         IF(IL.GT.NLHE2) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhe2(il,id)
      END DO
C
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      F00=3.906e-11*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(4.12E7*T+VTURB(ID))
C
C -------------------------------------------------------------------
C     overall loop over spectral series (only in the infrared region)
C -------------------------------------------------------------------
C
      ISERU=ILWHE2
      IF(ILWHE2.LE.3) THEN
         ISERL=ILWHE2
       ELSE IF(ILWHE2.LE.5) THEN
         ISERL=ILWHE2-1
       ELSE IF(ILWHE2.LE.7) THEN
         ISERL=ILWHE2-2
       ELSE IF(ILWHE2.LE.9) THEN
         ISERL=ILWHE2-3
       ELSE
         ISERL=ILWHE2-4
      END IF
C
      DO IJ=I0,I1
         ABSO(IJ)=0.
         EMIS(IJ)=0.
      END DO
C
      DO 200 I=ISERL,ISERU
      II=I*I
      XII=UN/II
      POPI=PJ(I)
C
C     determination of which He II lines contribute in a current 
C     frequency region
C
      M1=MHE10
      IF(I.LT.ILWHE2.AND.FRHE(I).GT.FREQ(2)) THEN
         M1=int(SQRT(FRHE(I)*II/(FRHE(I)-FREQ(2))))
      END IF
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.6..and.M1.LE.6.AND.I.EQ.2) GO TO 10
      IF(grav.lt.6..and.M1.LE.4.AND.I.EQ.1) GO TO 10
      M1=M1-1
      M2=MHE20+3
      IF(M2.GT.60) M2=60
   10 CONTINUE
      if(grav.gt.6.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3 
      end if
      IF(M1.LT.I+1) M1=I+1
      IF(M2.GT.60) M2=60
c     A=0.
c     E=0.
C
C     loop over lines which contribute at given wavelength region
C
      DO 100 J=M1,M2
         ILINE=0   
         JJ=J*J
         XJJ=UN/JJ
         ABTRA=PJ(I)*WNHE2(J,ID)
         EMTRA=PJ(J)*WNHE2(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1)
         IF(I.LE.2) THEN
            WLIN=227.838/(XII-1./JJ)
          ELSE 
            WLIN=227.7776/(XII-1./JJ)
         END IF
         IF(I.EQ.2) THEN
            IF(J.EQ.3.AND.IHE2PR.GT.0) ILINE=1
          ELSE IF(I.EQ.3) THEN
            IF(J.EQ.4.AND.IHE2PR.GT.0) ILINE=8
            IF(J.GT.5.AND.J.LE.10.AND.IHE2PR.GT.0) ILINE=J-3
          ELSE IF(I.EQ.4) THEN
            IF(J.LE.7.AND.IHE2PR.GT.0) ILINE=J+12
            IF(J.GE.8.AND.J.LE.15.AND.IHE2PR.GT.0) ILINE=J+1
         END IF
         IF(ILINE.GT.0) THEN
            NWL=NWLHE2(ILINE)
            DO IWL=1,NWL
               PRF0(IWL)=PRFHE2(ILINE,ID,IWL)
            END DO
            FID=CID*OSCHE2(ILINE)
            DO 50 IJ=I0,I1
               AL=ABS(WLAM(IJ)-WLIN)
               IF(AL.LT.1.E-4) AL=1.E-4
               AL=LOG10(AL)
               DO IWL=1,NWL-1
                  IW0=IWL
                  IF(AL.LE.WLHE2(ILINE,IWL+1)) GO TO 40
               END DO
   40          IW1=IW0+1
               PRFF=(PRF0(IW0)*(WLHE2(ILINE,IW1)-AL)+PRF0(IW1)*
     *             (AL-WLHE2(ILINE,IW0)))/
     *             (WLHE2(ILINE,IW1)-WLHE2(ILINE,IW0))
               SG=EXP(PRFF*AL10)*FID
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
   50       CONTINUE
          ELSE
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            FXK=F00*XKIJ
            FXK1=UN/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0*CINV*FXK1
            BETAD=DOP*DBETA
            FID=CID*FIJ*DBETA
c           FID0=CID1*FIJ0/DOP
            CALL DIVHE2(AD,DIV)
            DO IJ=I0,I1
               BETA=ABS(WLAM(IJ)-WL0)*FXK1
               SG=STARKA(BETA,AD,DIV,UN)*FID
c              if(fid0.gt.0.) then
c                 xd=beta/betad
c                 if(xd.lt.5.) sg=sg+exp(-xd*xd)*fid0
c              end if
               ABSO(IJ)=ABSO(IJ)+SG*ABTRA
               EMIS(IJ)=EMIS(IJ)+SG*EMTRA
            END DO
         END IF
  100 CONTINUE
  200 CONTINUE
C
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      DO IJ=I0,I1
         F=FREQ(IJ)
         F15=F*1.E-15
         XKF=EXP(-4.79928e-11*F*T1)
         XKFB=XKF*1.4743E-2*F15*F15*F15
         ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ)
         EMISH(IJ)=XKFB*EMIS(IJ)
      END DO
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE HE2LIW(ID,ABSOH,EMISH)
C     =================================
C
C     opacity and emissivity of He II lines  (these which are not considered
C     explicitly)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (UN=1.,SIXTH=1./6.)
      PARAMETER (CPP=4.1412E-16,CPJ=631479.) 
      PARAMETER (C00=1.25E-9,CDOP=1.284523E12,CID=0.02654,TWO=2.)
      PARAMETER (CPJ4=CPJ/4.,AL10=2.3025851,CINV=UN/2.997925E18)
      PARAMETER (CID1=0.01497)
      DIMENSION PJ(80),FRHE(12),OSCHE2(19),PRF0(36),
     *          ABSO(MFREQ),EMIS(MFREQ),ABSOH(MFREQ),EMISH(MFREQ)
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
      common/lasers/lasdel
      DATA FRHE /1.3158153D+16, 3.2895381D+15, 1.4624854D+15,
     *           8.2261878D+14, 5.2647201D+14, 3.6560459D+14,
     *           2.6860713D+14, 2.0565220D+14, 1.6249055D+14,
     *           1.3161730D+14, 1.0877460D+14, 9.1400851D+13/
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,
     *        1.604E-2, 1.023E-2, 6.980E-3,
     *        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,
     *        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,
     *        1.038,    1.793E-1, 6.549E-2/
C
      I=ILWHE2
      izz=2
      DO IJ=1,NFREQ
         ABSO(IJ)=0.
         EMIS(IJ)=0.
         ABSOH(IJ)=0.
         EMISH(IJ)=0.
      END DO
      IF(IFHE2.LE.0) RETURN
      T=TEMP(ID)
      T1=UN/T
      SQT=SQRT(T)
      ANE=ELEC(ID)
      ANES=EXP(SIXTH*LOG(ANE))
C
C     He III populations (either LTE or NLTE, depending on input model)
C
      IF(IELHE2.GT.0) THEN
         ANP=POPUL(NNEXT(IELHE2),ID)
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1
       ELSE
         ANP=RRR(ID,3,2)
         NLHE2=0
      END IF
C
C     populations of the first 60 levels of He II
C
      PP=CPP*ANE*ANP*T1/SQT
      DO IL=1,60
         X=IL*IL
         IIL=NFIRST(IELHE2)+IL-1
         IF(IL.LE.NLHE2) PJ(IL)=POPUL(IIL,ID)
         IF(IL.GT.NLHE2) PJ(IL)=PP*EXP(CPJ/X*T1)*X*wnhe2(il,id)
      END DO
C
C     Frequency- and line-independent parameters for evaluating the
C     asymptotic Stark profile
C
      F00=3.906e-11*ANES*ANES*ANES*ANES
      DOP0=1.E8*SQRT(4.12E7*T+VTURB(ID))
C
C -------------------------------------------------------------------
C     overall loop over spectral series (only in the infrared region)
C -------------------------------------------------------------------
C
      DO 300 IJ=1,NFREQ
      ABSO(IJ)=0.
      EMIS(IJ)=0.
      IF(IHE2LW(IJ).le.0) GO TO 300
      I=ILWHEW(IJ)
      FR=FREQ(IJ)
      ISERU=ILWHEW(IJ)
      IF(ILWHEW(IJ).LE.3) THEN
         ISERL=ILWHEW(IJ)
       ELSE IF(ILWHEW(IJ).LE.5) THEN
         ISERL=ILWHEW(IJ)-1
       ELSE IF(ILWHEW(IJ).LE.7) THEN
         ISERL=ILWHEW(IJ)-2
       ELSE IF(ILWHEW(IJ).LE.9) THEN
         ISERL=ILWHEW(IJ)-3
       ELSE
         ISERL=ILWHEW(IJ)-4
      END IF
C
C
      DO 200 I=ISERL,ISERU
      II=I*I
      XII=UN/II
      PLTEI=PP*EXP(CPJ*T1*XII)*II
      POPI=PJ(I)
C
C     determination of which He II lines contribute in a current 
C     frequency region
C
      M1=MHE10W(IJ)
      IF(I.LT.ILWHEW(IJ).AND.FRHE(I).GT.FR) THEN
         M1=int(SQRT(FRHE(I)*II/(FRHE(I)-FR)))
      END IF
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.6..and.M1.LE.6.AND.I.EQ.2) GO TO 10
      IF(grav.lt.6..and.M1.LE.4.AND.I.EQ.1) GO TO 10
      M1=M1-1
      M2=MHE20W(IJ)+3
      IF(M2.GT.60) M2=60
   10 CONTINUE
      if(grav.gt.6.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3 
      end if
      IF(M1.LT.I+1) M1=I+1
      IF(M2.GT.60) M2=60
C
C     loop over lines which contribute at given wavelength region
C
      DO 100 J=M1,M2
         ILINE=0   
         JJ=J*J
         XJJ=UN/JJ
         ABTRA=PJ(I)*WNHE2(J,ID)
         EMTRA=PJ(J)*WNHE2(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1)
         IF(I.LE.2) THEN
            WLIN=227.838/(XII-1./JJ)
          ELSE 
            WLIN=227.7776/(XII-1./JJ)
         END IF
         IF(I.EQ.2) THEN
            IF(J.EQ.3.AND.IHE2PR.GT.0) ILINE=1
          ELSE IF(I.EQ.3) THEN
            IF(J.EQ.4.AND.IHE2PR.GT.0) ILINE=8
            IF(J.GT.5.AND.J.LE.10.AND.IHE2PR.GT.0) ILINE=J-3
          ELSE IF(I.EQ.4) THEN
            IF(J.LE.7.AND.IHE2PR.GT.0) ILINE=J+12
            IF(J.GE.8.AND.J.LE.15.AND.IHE2PR.GT.0) ILINE=J+1
         END IF
         IF(ILINE.GT.0) THEN
            NWL=NWLHE2(ILINE)
            DO IWL=1,NWL
               PRF0(IWL)=PRFHE2(ILINE,ID,IWL)
            END DO
            FID=CID*OSCHE2(ILINE)
            AL=ABS(WLAM(IJ)-WLIN)
            IF(AL.LT.1.E-4) AL=1.E-4
            AL=LOG10(AL)
            DO IWL=1,NWL-1
               IW0=IWL
               IF(AL.LE.WLHE2(ILINE,IWL+1)) GO TO 40
            END DO
   40       IW1=IW0+1
            PRFF=(PRF0(IW0)*(WLHE2(ILINE,IW1)-AL)+PRF0(IW1)*
     *          (AL-WLHE2(ILINE,IW0)))/
     *          (WLHE2(ILINE,IW1)-WLHE2(ILINE,IW0))
            SG=EXP(PRFF*AL10)*FID
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA
          ELSE
            CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
            FXK=F00*XKIJ
            FXK1=UN/FXK
            DOP=DOP0/WL0
            DBETA=WL0*WL0*CINV*FXK1
            BETAD=DOP*DBETA
            FID=CID*FIJ*DBETA
            CALL DIVHE2(AD,DIV)
            BETA=ABS(WLAM(IJ)-WL0)*FXK1
            SG=STARKA(BETA,AD,DIV,UN)*FID
            ABSO(IJ)=ABSO(IJ)+SG*ABTRA
            EMIS(IJ)=EMIS(IJ)+SG*EMTRA
         END IF
  100 CONTINUE
  200 CONTINUE
C
C     ----------------------------
C     total opacity and emissivity
C     ----------------------------
C
      F=FREQ(IJ)
      F15=F*1.E-15
      XKF=EXP(-4.79928e-11*F*T1)
      XKFB=XKF*1.4743E-2*F15*F15*F15
      ABSOH(IJ)=ABSO(IJ)-XKF*EMIS(IJ)
      EMISH(IJ)=XKFB*EMIS(IJ)
  300 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE STARK0(I,J,IZZ,XKIJ,WL0,FIJ,FIJ0)
C
C     Auxiliary procedure for evaluating the approximate Stark profile
C     of hydrogen lines - sets up necessary frequency independent
C     parameters
C
C     Input:  I     - principal quantum number of the lower level
C             J     - principal quantum number of the upper level
C             IZZ   - ionic charge (IZZ=1 for hydrogen, etc.)
C     Output: XKIJ  - coefficients K(i,j) for the Hotzmark profile;
C                     exact up to j=6, asymptotic for higher j
C             WL0   - wavelength of the line i-j
C             FIJ   - Stark f-value for the line i-j
C             FIJ0  - f-value for the undisplaced component of the line
C
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (RYD1=911.763811,RYD2=911.495745,CXKIJ=5.5E-5)
      PARAMETER (WI1=911.753578, WI2=227.837832)
      PARAMETER (UN=1.,TEN=10.,TWEN=20.,HUND=100.)
      DIMENSION FSTARK(10,4),XKIJT(5,4),FOSC0(10,4),FADD(5,5)
      DATA XKIJT/3.56E-4,5.23E-4,1.09E-3,1.49E-3,2.25E-3,.0125,.0177,
     * .028,.0348,.0493,.124,.171,.223,.261,.342,.683,.866,1.02,1.19,
     * 1.46/
      DATA FSTARK/  .1387,    .0791,   .02126,   .01394,   .00642,
     *           4.814E-3, 2.779E-3, 2.216E-3, 1.443E-3, 1.201E-3,
     *              .3921,    .1193,   .03766,   .02209,   .01139,
     *           8.036E-3, 5.007E-3,  3.85E-3, 2.658E-3, 2.151E-3,
     *              .6103,    .1506,   .04931,   .02768,   .01485,
     *             .01023, 6.588E-3, 4.996E-3, 3.524E-3, 2.838E-3,
     *              .8163,    .1788,   .05985,   .03189,   .01762,
     *             .01196, 7.825E-3, 5.882E-3, 4.233E-3, 3.375E-3/
      DATA FOSC0 / 0.27746,  0., 0.00773,  0., 0.00134, 0., 
     *             0.000404, 0., 0.000162, 0.,
     *             0.24869,  0., 0.00701,  0., 0.00131, 0.,
     *             0.000422, 0., 0.000177, 0.,
     *             0.23175,  0., 0.00653,  0., 0.00118, 0.,
     *             0.000392, 0., 0.000169, 0.,
     *             0.22148,  0.0005, 0.00563, 0.0004, 0.00108, 0.,
     *             0.000362, 0., 0.000159, 0./ 
      DATA FADD /  1.231, 0.2069, 7.448E-2, 3.645E-2, 2.104E-2,
     *             1.424, 0.2340, 8.315E-2, 4.038E-2, 2.320E-2,
     *             1.616, 0.2609, 9.163E-2, 4.416E-2, 2.525E-2,
     *             1.807, 0.2876, 1.000E-1, 4.787E-2, 2.724E-2,
     *             1.999, 0.3143, 1.083E-1, 5.152E-2, 2.918E-2/
C
      II=I*I
      JJ=J*J
      JMIN=J-I
      IF(JMIN.LE.5.and.i.le.4) THEN
         XKIJ=XKIJT(JMIN,I)
       ELSE
         XKIJ=CXKIJ*(II*JJ)*(II*JJ)/(JJ-II)
      END IF
      IF(I.LE.4) THEN
         IF(JMIN.LE.10) THEN
            FIJ=FSTARK(JMIN,I)
            FIJ0=FOSC0(JMIN,I)
          ELSE 
            CFIJ=((TWEN*I+HUND)*J/(I+TEN)/(JJ-II))
            FIJ=FSTARK(10,I)*CFIJ*CFIJ*CFIJ
            FIJ0=0.
         END IF
       ELSE IF(I.LE.9) THEN
         IF(JMIN.LE.5) THEN
            FIJ=FADD(JMIN,I-4)
            FIJ0=0.
          ELSE 
            CFIJ=((TEN*I+25.)*J/(I+5.)/(JJ-II))
            FIJ=FADD(5,I-4)*CFIJ*CFIJ*CFIJ
            FIJ0=0.
         END IF
       ELSE
         CFIJ=UN*J/(JJ-II)
         FIJ=1.96*I*CFIJ*CFIJ*CFIJ
         FIJ0=0.
      END IF
C
C     wavelength with an explicit correction to the air wavalength
C
      w0=wi1
      if(izz.eq.2) w0=wi2
      WL0=W0/(UN/II-UN/JJ)
      IF(WL0.GT.vaclim) THEN
         ALM=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0/(XN1*1.D-6+UN)
      END IF        
      RETURN
      END
C
C ********************************************************************
C
      FUNCTION STARKA(BETA,A,DIV,FAC)
C
C     Approximate expressions for the hydrogen Stark profile
C
C     Input: BETA  - delta lambda in beta units,
C            BETAD - Doppler width in beta units
C            A     - auxiliary parameter
C                    A=1.5*LOG(BETAD)-1.671
C            DIV   - only for A > 1; division point between Doppler
C                    and asymptotic Stark wing, expressed in units
C                    of betad.
C                    DIV = solution of equation
C                    exp(-(beta/betad)**2)/betad/sqrt(pi)=
C                     = 1.5*FAC*beta**-5/2
C                    (ie. the point where Doppler profile is equal to
C                     the asymptotic Holtsmark)
C                    In order to save computer time, the division point
C                    DIV is calculated in advance by routine DIVSTR.
C            FAC   - factor by which the Holtsmark profile is to be
C                    multiplied to get total Stark Profile
C                    FAC should be taken to 2 for hydrogen, (and =1
C                    for He II)
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (F0=-0.5758228,F1=0.4796232,F2=0.07209481/2.,AL=1.26)
      PARAMETER (SD=0.5641895,SLO=-2.5,TRHA=1.5,BL1=1.52,BL2=8.325)
      PARAMETER (SAC=0.07966/2.)
      XD=BETA/BETAD
C
C     for a > 1 Doppler core + asymptotic Holtzmark wing with division
C               point DIV
C
      IF(A.GT.AL) THEN
         IF(XD.LE.DIV) THEN
            STARKA=SD*EXP(-XD*XD)/BETAD
          ELSE
            STARKA=TRHA*FAC*EXP(SLO*LOG(BETA))
         END IF
       ELSE
C
C     empirical formula for a < 1
C
         IF(BETA.LE.BL1) THEN
            STARKA=SAC*FAC
          ELSE IF(BETA.LT.BL2) THEN
            XL=LOG(BETA)
            FL=(F0*XL+F1)*XL
            STARKA=F2*FAC*EXP(FL)
          ELSE
            STARKA=TRHA*FAC*EXP(SLO*LOG(BETA))
         END IF
      END IF
      RETURN
      END
C
C *******************************************************************
C *******************************************************************
C
      FUNCTION STARKIR(II,JJ,T,ANE,BETA)
C     ==================================
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (PI=3.14159265,PI2=2.*PI,
     *           OS0=0.026564,RYD=3.28805E15,
     *           Y2CON=PI*PI*0.5/OS0/CL)
C
      DEL=BETA/DBETA
      HKT=HK/T
      XII=II
      XJJ=JJ
      XX=XII/XJJ
      DD=2.*XJJ*RYD/DEL
      Y1=XJJ*DEL*0.5*HKT
      Y2=Y2CON*DEL**2/ANE
      QSTAT=1.5+.5*(Y1**2-1.384)/(Y1**2+1.384)
      QIMPA=0.
      IF(Y1.GT.8..OR.Y1.GE.Y2) GO TO 10
      EXY2=0.
      IF(Y2.LE.8.) EXY2=EXPINT(Y2)
      QIMPA=1.438*SQRT(Y1*(1.-XX))*(.4*EXP(-Y1)+EXPINT(Y1)-.5*EXY2)
   10 IF(BETA.GT.20.) GO TO 20
      PROF=8./(80.+BETA**3)
      RATIO=QSTAT+QIMPA
      GO TO 30
   20 PROF=1.5/BETA/BETA/SQRT(BETA)
      DIOI=PI2*1.48E-25*DD*ANE*(SQRT(DD)*
     *     (1.3*QSTAT+.3*QIMPT)-3.9*RYD*HKT)
      RATIO=QSTAT*MIN(1.+DIOI,1.25)+QIMPA
   30 STARKIR=PROF*RATIO
      RETURN
      END


C
C *******************************************************************
C *******************************************************************
C
      SUBROUTINE DIVSTR(A,DIV)
C     ==============================
C
C     Auxiliary procedure for STARKA - determination of the division
C     point between Doppler and asymptotic Stark profiles
C
C     Input:  BETAD - Doppler width in beta units
C     Output: A     - auxiliary parameter
C                     A=1.5*LOG(BETAD)-1.671
C             DIV   - only for A > 1; division point between Doppler
C                     and asymptotic Stark wing, expressed in units
C                     of betad.
C                     DIV = solution of equation
C                     exp(-(beta/betad)**2)/betad/sqrt(pi)=3*beta**-5/2
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (UN=1.,TWO=2.,UNQ=1.25,UNH=1.5,TWH=2.5,FO=4.,FI=5.)
      PARAMETER (CA=1.671,BL=5.821,AL=1.26,CX=0.28,DX=0.0001)
C
      A=UNH*LOG(BETAD)-CA
      IF(BETAD.LT.BL) RETURN
      IF(A.GE.AL) THEN
         X=SQRT(A)*(UN+UNQ*LOG(A)/(FO*A-FI))
      ELSE
         X=SQRT(CX+A)
      ENDIF
      DO I=1,5
         XN=X*(UN-(X*X-TWH*LOG(X)-A)/(TWO*X*X-TWH))
         IF(ABS(XN-X).LE.DX) GO TO 20
         X=XN
      END DO
   20 DIV=X
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE HYDINI
C
C     Initializes necessary arrays for evaluating hydrogen line profiles
C     from the Lemke, Tremblay-Bergeron, or Schoening-Butler tables
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
c     DIMENSION WLINE(4,22)
      DIMENSION IILW(100),IIUP(100)
      CHARACTER*1 CHAR
      DATA INIT /0/
C
      IF(INIT.EQ.0) THEN
         DO I=1,4
            DO J=I+1,22
               CALL STARK0(I,J,IZZ,XK,WL0,FIJ,FIJ0)
               WLINE(I,J)=WL0
c              OSCH(I,J)=FIJ+FIJ0
            END DO
         END DO
         INIT=1
      END IF
      DO I=1,4
         DO J=1,22
            ILIN0(I,J)=0
         END DO
      END DO
C
C --------------------------------------------
C     Schoening-Butler tables - for IHYDPR < 0
C --------------------------------------------
C
      IF(IHYDPR.LT.0) THEN
      IHYDPR=67
      ILEMKE=0
      NLINE=12
c
      OPEN(UNIT=IHYDPR,FILE='./data/hydprf.dat',STATUS='OLD')
      write(6,*) ' reading Schoening-Butler tables'
C
      DO I=1,12
         READ(IHYDPR,500)
      END DO
      DO 100 ILINE=1,NLINE
C
C     read the tables, which have to be stored in file
C     unit IHYDPR (which is the input parameter in the progarm)
C
         READ(IHYDPR,501) I,J
         IF(ILINE.EQ.12) J=10
         WL0=WLINE(I,J)
         ILIN0(I,J)=ILINE
         READ(IHYDPR,*) CHAR,NWL,(WL(I,ILINE),I=1,NWL)
         READ(IHYDPR,*) CHAR,NT,(XT(I,ILINE),I=1,NT)
         READ(IHYDPR,*) CHAR,NE,(XNE(I,ILINE),I=1,NE)
         READ(IHYDPR,500)
         NWLH(ILINE)=NWL
         NWLHYD(ILINE)=NWL
         NTH(ILINE)=NT
         NEH(ILINE)=NE
C
         DO I=1,NWL
            IF(WL(I,ILINE).LT.1.E-4) WL(I,ILINE)=1.E-4
            WLHYD(ILINE,I)=LOG10(WL(I,ILINE))
         END DO
C
         DO IE=1,NE
            DO IT=1,NT
               READ(IHYDPR,500)
               READ(IHYDPR,*) (PRF(IWL,IT,IE,ILINE),IWL=1,NWL)
            END DO
         END DO
C
C        coefficient for the asymptotic profile is determined from
C        the input data
C
         XCLOG=PRF(NWL,1,1,ILINE)+2.5*LOG10(WL(NWL,ILINE))+31.5304-
     *         XNE(1,ILINE)-2.*LOG10(WL0)
         XKLOG=0.6666667*(XCLOG-0.176)
         XK=EXP(XKLOG*2.3025851)
C
         DO ID=1,ND
C
C           temperature is modified in order to account for the
C           effect of turbulent velocity on the Doppler width
C
            T=TEMP(ID)+6.06E-9*VTURB(ID)
            ANE=ELEC(ID)
            TL=LOG10(T)
            ANEL=LOG10(ANE)
            F00=1.25E-9*ANE**0.666666667
            FXK=F00*XK
            DOP=1.E8/WL0*SQRT(1.65E8*T)
            DBETA=WL0*WL0/2.997925E18/FXK
            BETAD=DBETA*DOP
C
C       interpolation to the actual values of temperature and electron
C       density. The result is stored at array PRFHYD, having indices
C       ILINE (line number: 1 for L-alpha,..., 4 for H-delta, etc.);
C                           5 for H-alpha,..., 8 for H-delta, etc.)
C       ID - depth index
C       IWL - wavelength index 
C
            DO IWL=1,NWL
               CALL INTHYD(PROF,TL,ANEL,IWL,ILINE)
               PRFHYD(ILINE,ID,IWL)=PROF
            END DO
         END DO
  100 CONTINUE
      CLOSE(IHYDPR)
C
  500 FORMAT(1X)
  501 FORMAT(12X,I1,9X,I1)
C
      IHYDPR=-IHYDPR
      RETURN
      END IF
C
C ---------------------------------
C     read Lemke or Tremblay tables
C ---------------------------------
C
      if(ihydpr.lt.20) ihydpr=ihydpr+20
      if(ihydpr.eq.21) then
         open(unit=ihydpr,file='./data/lemke.dat',status='old')
         write(6,641) ihydpr
       else if(ihydpr.eq.22) then
         open(unit=ihydpr,file='./data/tremblay.dat',status='old')
         write(6,642) ihydpr
      end if
  641 format(' -----------'/
     *       ' reading Lemke tables; ihydpr =',i3,/
     *       ' -----------')
  642 format(' -----------'/
     *       ' reading Tremblay tables; ihydpr =',i3,/
     *       ' -----------')
C
      ILEMKE=1
      READ(IHYDPR,*) NTAB
      write(6,611) ntab 
  611 format(' ntab',i4)
      DO ITAB=1,NTAB
         ILINEB=ILINE
         READ(IHYDPR,*) NLLY
         DO ILI=1,NLLY
            ILINE=ILINE+1
            READ(IHYDPR,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,
     *                     NWL,NE,NT
            WL0=WLINE(I,J)
            ILIN0(I,J)=ILINE
            NWLH(ILINE)=NWL
            NWLHYD(ILINE)=NWL
            NTH(ILINE)=NT
            NEH(ILINE)=NE
            iilw(iline)=i
            iiup(iline)=j
            DO IWL=1,NWL
               WL(IWL,ILINE)=ALMIN+(IWL-1)*DLA
               WLHYD(ILINE,IWL)=WL(IWL,ILINE)
               WL(IWL,ILINE)=EXP(2.3025851*WL(IWL,ILINE))
            END DO
            DO INE=1,NE
               XNE(INE,ILINE)=ANEMIN+(INE-1)*DLE
            END DO
            DO IT=1,NT
               XT(IT,ILINE)=TMIN+(IT-1)*DLT
            END DO 
         END DO
c
         DO ILI=1,NLLY         
            ILNE=ILINEB+ILI
            NWL=NWLH(ILNE)
            READ(IHYDPR,500)
            DO INE=1,NEH(ILNE)
               DO IT=1,NTH(ILNE)
                  READ(IHYDPR,*) QLT,(PRF(IWL,IT,INE,ILNE),IWL=1,NWL)
               END DO
            END DO
C
            i=iilw(ilne)
            j=iiup(ilne)
            DO ID=1,ND
               CALL HYDTAB(I,J,ID)
            END DO
         END DO
      END DO
      NLIHYD=ILNE
      CLOSE(IHYDPR)
C
      RETURN
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE HYDTAB(I,J,ID)
C
C     interpolated hydrogen line broadening table for line I->J and
C     for parameters (TEMP, ELEC) at depth ID
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
C
      ILINE=ILIN0(I,J)
      IF(ILINE.EQ.0) RETURN
      WL0=WLINE(I,J)
      NWL=NWLH(ILINE)
C
C     coefficient for the asymptotic profile is determined from
C     the input data
C
      if(id.eq.1) then
         XCLOG=PRF(NWL,1,1,ILINE)+2.5*WLHYD(ILINE,NWL)-0.477121
         XKLOG=0.6666667*XCLOG
         XK=EXP(XKLOG*2.3025851)
      end if
C
C     temperature is modified in order to account for the
C     effect of turbulent velocity on the Doppler width
C
      T=TEMP(ID)+6.06E-9*VTURB(ID)
      ANE=ELEC(ID)
      TL=LOG10(T)
      ANEL=LOG10(ANE)
      F00=1.25E-9*ANE**0.666666667
      FXK=F00*XK
      DOP=1.E8/WL0*SQRT(1.65E8*T)
      DBETA=WL0*WL0/2.997925E18/FXK
      BETAD=DBETA*DOP
C
C     interpolation to the actual values of temperature and electron
C     density. The result is stored at array PRFHYD, having indices
C       ILINE - line number
C       ID    - depth index
C       IWL   - wavelength index
C
      DO IWL=1,NWL
         CALL INTHYD(PROF,TL,ANEL,IWL,ILINE)
         PRFHYD(ILINE,ID,IWL)=PROF
      END DO
C
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE INTHYD(W0,X0,Z0,IWL,ILINE)
C
C     Interpolation in temperature and electron density from the
C     hydrogen odening tables to the actual valus of
C     temperature and electron density
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (TWO=2.)
      DIMENSION ZZ(3),XX(3),WX(3),WZ(3)
C
      NX=3
      NZ=3
      NT=NTH(ILINE)
      NE=NEH(ILINE)
      BETA=WL(IWL,ILINE)/FXK
      IF(ILEMKE.EQ.1) THEN
         BETA=WL(IWL,ILINE)/XK
         NX=2
         NZ=2
      END IF
C
C     for values lower than the lowest grid value of electron density
C     the profiles are determined by the approximate expression
C     (see STARKA); not by an extrapolation in the HYD tables which may
C     be very inaccurate
C
      IF(Z0.LT.XNE(1,ILINE)*0.99.OR.Z0.GT.XNE(NE,ILINE)*1.01) THEN
         CALL DIVSTR(A,DIV)
         W0=STARKA(BETA,A,DIV,TWO)*DBETA
         W0=LOG10(W0)
         GO TO 500
      END IF
C
C     Otherwise, one interpolates (or extrapolates for higher than the
C     highes grid value of electron density) in the HYD tables
C
      DO IZZ=1,NE-1
         IPZ=IZZ
         IF(Z0.LE.XNE(IZZ+1,ILINE)) GO TO 20
      END DO
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
      N1Z=N0Z+NZ-1
C
      DO 300 IZZ=N0Z,N1Z
         I0Z=IZZ-N0Z+1
         ZZ(I0Z)=XNE(IZZ,ILINE)
C
C     Likewise, the approximate expression instead of extrapolation
C     is used for higher that the highest grid value of temperature,
C     if the Doppler width expressed in beta units (BETAD) is
C     sufficiently large (> 10)
C
         IF(X0.GT.1.01*XT(NT,ILINE).AND.BETAD.GT.10.) THEN
            CALL DIVSTR(A,DIV)
            W0=STARKA(BETA,A,DIV,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500
         END IF
C
C     Otherwise, normal inter- or extrapolation
C
C     Both interpolations (in T as well as in electron density) are
C     by default the quadratic interpolations in logarithms
C
         DO IX=1,NT-1
            IPX=IX
            IF(X0.LE.XT(IX+1,ILINE)) GO TO 40
         END DO
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
         N1X=N0X+NX-1
         DO IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XT(IX,ILINE)
            WX(I0)=PRF(IWL,IX,IZZ,ILINE)
         END DO
         IF(WX(1).LT.-99..OR.WX(2).LT.-99..OR.WX(3).LT.-99.) THEN
            CALL DIVSTR(A,DIV)
            W0=STARKA(BETA,A,DIV,TWO)*DBETA
            W0=LOG10(W0)
            GO TO 500
          ELSE
            WZ(I0Z)=YINT(XX,WX,X0)
         END IF
  300 CONTINUE
      W0=YINT(ZZ,WZ,Z0)
  500 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
      FUNCTION YINT(XL,YL,XL0)
C
C     Quadratic interpolation routine
C
C     Input:  XL - array of x
C             YL - array of f(x)
C             XL0 - the point x(0) to which one interpolates
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION XL(3),YL(3)
      A0=(XL(2)-XL(1))*(XL(3)-XL(2))*(XL(3)-XL(1))
      A1=(XL0-XL(2))*(XL0-XL(3))*(XL(3)-XL(2))
      A2=(XL0-XL(1))*(XL(3)-XL0)*(XL(3)-XL(1))
      A3=(XL0-XL(1))*(XL0-XL(2))*(XL(2)-XL(1))
      YINT=(YL(1)*A1+YL(2)*A2+YL(3)*A3)/A0
      RETURN
      END
C
C ********************************************************************
C
C
 
      SUBROUTINE HE1INI
C     =================
C
C     Initializes necessary arrays for evaluating the He I line
C     absorption profiles using data calculated by Barnard, Cooper
C     and Smith JQSRT 14, 1025, 1974 (for 4471)
C     or Shamey, unpublished PhD thesis, 1969 (for other lines)
C
C     This procedure is quite analogous to HYDINI for hydrogen lines
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      COMMON/PROHE1/PRFHE1(50,4,8,3),DLMHE1(50,8,3),XNEHE1(8),
     *              NWLAM(8,4)
      COMMON/PRO447/PRF447(80,4,7),DLM447(80,7),XNE447(7)
      DATA NT /4/
C
      IH=67
c     OPEN(UNIT=IH,STATUS='OLD')
C
C        read the Barnard, Cooper, Smith tables for He I 4471 line,
C        which have to be stored in file unit IH
C
      NE=7
      DO IE=1,NE
         READ(IH,501) IL,WL0,IE1,XXNE,NWL
         NWLAM(IE,1)=NWL
         XNE447(IE)=LOG10(XXNE)
         DO I=1,NWL
            READ(IH,502) DLM447(I,IE),
     *                (PRF447(I,IT,IE),IT=1,NT)
         END DO
      END DO
C
C     read Shamey's tables for He I 4387, 4026, and 4922 lines
C     which have to be stored in file unit IH
C
      NE=8
      DO ILN=1,3
         DO IE=1,NE
            READ(IH,501) IL,WL0,IE1,XXNE,NWL
            NWLAM(IE,ILN+1)=NWL
            XNEHE1(IE)=LOG10(XXNE)
            DO I=1,NWL
               READ(IH,*) DLMHE1(I,IE,ILN),
     *                    (PRFHE1(I,IT,IE,ILN),IT=1,NT)
            END DO
         END DO
      END DO
      CLOSE(IH)
C
  501 FORMAT(/9X,I2,7X,F10.3,13X,I2,6X,E8.1,7X,I3/)
  502 FORMAT(5E10.2)
      RETURN
      END
C
C ********************************************************************
C
 
      FUNCTION WTOT(T,ANE,ID,ILINE)
C     =============================
C
C     Evaluates the total (electron + ion) impact Stark width
C     for four HeI lines
C     After Griem (1974); and Barnard, Cooper, Smith (1974) JQSRT 14,
C     1025 for the 4471 line
C
C     Input: T     - temperature
C            ANE   - electron density
C            ID    - depth index
C            ILINE - index of the line ( = 1  for 4471,
C                                        = 2  for 4387,
C                                        = 3  for 4026,
C                                        = 4  for 4922)
C     Output: WTOT - Stark width in Angstroms
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION ALPH0(4,4),W0(4,4),ALAM0(4)
      DATA ALPH0 / 0.107, 0.119, 0.134, 0.154,
     *             0.206, 0.235, 0.272, 0.317,
     *             0.172, 0.193, 0.218, 0.249,
     *             0.121, 0.136, 0.157, 0.184/
      DATA W0    / 1.460, 1.269, 1.079, 0.898,
     *             6.130, 5.150, 4.240, 3.450,
     *             4.040, 3.490, 2.960, 2.470,
     *             2.312, 1.963, 1.624, 1.315/
      DATA ALAM0 / 4471.50, 4387.93, 4026.20, 4921.93/
C
      I=JT(ID)
      ALPHA=(TI0(ID)*ALPH0(I,ILINE)+TI1(ID)*ALPH0(I-1,ILINE)+
     *      TI2(ID)*ALPH0(I-2,ILINE))*(ANE*1.E-13)**0.25
      WE=   (TI0(ID)*W0(I,ILINE)+TI1(ID)*W0(I-1,ILINE)+
     *      TI2(ID)*W0(I-2,ILINE))*ANE*1.E-16
      F0=1.884E19/ALAM0(ILINE)/ALAM0(ILINE)
      SIG=(4.32E-5*WE/SQRT(T)*F0/ANE**0.3333)**0.3333
      WTOT=WE*(1.+1.36/SIG*ALPHA**0.8889)
      RETURN
      END
 
C
C ********************************************************************
C
 
      FUNCTION EXTPRF(DLAM,IT,ILINE,ANEL,DLAST,PLAST)
C     ===============================================
C
C     Extrapolation in wavelengths in Shamey, or Barnard, Cooper,
C     Smith tables
C     Special formula suggested by Cooper
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION W0(4,4)
      DATA W0    / 1.460, 1.269, 1.079, 0.898,
     *             6.130, 5.150, 4.240, 3.450,
     *             4.040, 3.490, 2.960, 2.470,
     *             2.312, 1.963, 1.624, 1.315/
C
      WE=W0(IT,ILINE)*EXP(ANEL*2.3025851)*1.E-16
      DLASTA=ABS(DLAST)
      D52=DLASTA*DLASTA*SQRT(DLASTA)
      F=D52*(PLAST-WE/3.14159/DLAST/DLAST)
      EXTPRF=(WE/3.14159+F/SQRT(ABS(DLAM)))/DLAM/DLAM
      RETURN
      END
 
C
C ********************************************************************
C
 
      FUNCTION PHE1(ID,FREQ,ILINE)
C     ============================
C
C     Absorption profile for four lines of He I, given by
C     Barnard, Cooper, Smith (1974) JQSRT 14, 1025 for the 4471 line;
C     Shamey (1969) PhD thesis, for other lines
C
C     Input: ID    - depth index
C            FREQ  - frequency
C            ILINE - index of the line ( = 1  for 4471,
C                                        = 2  for 4387,
C                                        = 3  for 4026,
C                                        = 4  for 4922)
C
C     Output: PHE1 - profile coefficient in frequency units,
C                    normalized to sqrt(pi) [not unity]
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (NT=4)
      COMMON/PROHE1/PRFHE1(50,NT,8,3),DLMHE1(50,8,3),XNEHE1(8),
     *              NWLAM(8,NT)
      COMMON/PRO447/PRF447(80,NT,7),DLM447(80,7),XNE447(7)
      DIMENSION WLAM0(4),XT0(NT),XX(3),WX(3),YY(2),PP(2),ZZ(3),WZ(3)
      DATA WLAM0 / 4471.50, 4387.93, 4026.20, 4921.93/
      DATA XT0/ 3.699, 4.000, 4.301, 4.602/
C
C     temperature is modified in order to account for the
C     effect of turbulent velocity on the Doppler width
C
      T=TEMP(ID)+2.42E-8*VTURB(ID)
      TL=LOG10(T)
      ANE=ELEC(ID)
      ANEL=LOG10(ANE)
      ALAM=2.997925E18/FREQ
      DLAM=ALAM-WLAM0(ILINE)
      DOPL=SQRT(4.125E7*T)*WLAM0(ILINE)/2.997925E10
C
      IF(TL.GT.XT0(NT)+0.1) GO TO 5
      IF(ILINE.EQ.1.AND.ANEL.GE.XNE447(1)) GO TO 10
      IF(ILINE.NE.1.AND.ANEL.GE.XNEHE1(1)) GO TO 10
C
C     isolated line approximation for low electron densities
C
    5 A=WTOT(T,ANE,ID,ILINE)/DOPL
      V=ABS(DLAM)/DOPL
      V1=ABS(ALAM-4471.682)/DOPL
      PHE1=VOIGTK(A,V)
      IF(ILINE.EQ.1) PHE1=(8.*PHE1+VOIGTK(A,V1))/9.
      RETURN
C
C     otherwise, interpolation (or extrapolation) in tables
C
   10 NX=3
      NZ=3
      NY=2
      NE=8
      ILNE=ILINE-1
      IF(ILINE.EQ.1) NE=7
C
C     Interpolation in electron density
C
      DO JZ=1,NE-1
         IPZ=JZ
         IF(ILINE.EQ.1.AND.ANEL.LE.XNE447(JZ+1)) GO TO 30
         IF(ILINE.NE.1.AND.ANEL.LE.XNEHE1(JZ+1)) GO TO 30
      END DO
   30 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
      N1Z=N0Z+NZ-1
      DO 300 JZ=N0Z,N1Z
         I0Z=JZ-N0Z+1
         IF(ILINE.EQ.1) ZZ(I0Z)=XNE447(JZ)
         IF(ILINE.NE.1) ZZ(I0Z)=XNEHE1(JZ)
C
C        Interpolation in temperature
C
         DO IX=1,NT-1
            IPX=IX
            IF(TL.LE.XT0(IX+1)) GO TO 50
         END DO
   50    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
         N1X=N0X+NX-1
         DO 200 IX=N0X,N1X
            I0X=IX-N0X+1
            XX(I0X)=XT0(IX)
C
C           Interpolation in wavelength
C
C           1. For delta lambda beyond tabulated values - special
C              extrapolation (Cooper's suggestion)
C
            NLST=NWLAM(JZ,ILINE)
            IF(ILINE.EQ.1) THEN
               D1=DLM447(1,JZ)
               D2=DLM447(NLST,JZ)
               IF(DLAM.LT.D1) THEN
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D1,PRF447(1,IX,JZ))
                  GO TO 150
                ELSE IF(DLAM.GT.D2) THEN
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D2,
     *                PRF447(NLST,IX,JZ))
                  GO TO 150
               END IF
             ELSE
               D1=DLMHE1(1,JZ,ILNE)
               D2=DLMHE1(NLST,JZ,ILNE)
               IF(DLAM.LT.D1) THEN
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D1,
     *                PRFHE1(1,IX,JZ,ILNE))
                  GO TO 150
                ELSE IF(DLAM.GT.D2) THEN
                  PRF0=EXTPRF(DLAM,IX,ILINE,ZZ(I0Z),D2,
     *                PRFHE1(NLST,IX,JZ,ILNE))
                  GO TO 150
               END IF
            END IF
C
C           normal linear interpolation in wavelength
C           (for 4471, linear interpolation in logarithms)
C
            DO IY=1,NLST-1
               IPY=IY
               IF(ILINE.EQ.1.AND.DLAM.LE.DLM447(IY+1,JZ)) GO TO 70
               IF(ILINE.NE.1.AND.DLAM.LE.DLMHE1(IY+1,JZ,ILNE))
     *            GO TO 70
            END DO
   70       N0Y=IPY-NY/2+1
            IF(N0Y.LT.1) N0Y=1
            IF(N0Y.GT.NLST-NY+1) N0Y=NLST-NY+1
            N1Y=N0Y+NY-1
            DO IY=N0Y,N1Y
               I0=IY-N0Y+1
               IF(ILINE.EQ.1) YY(I0)=DLM447(IY,JZ)
               IF(ILINE.EQ.1) PP(I0)=LOG(PRF447(IY,IX,JZ))
               IF(ILINE.NE.1) YY(I0)=DLMHE1(IY,JZ,ILNE)
               IF(ILINE.NE.1) PP(I0)=PRFHE1(IY,IX,JZ,ILNE)
           END DO
           IF(ILINE.NE.1) THEN
              WX(I0X)=(PP(2)*(DLAM-YY(1))+PP(1)*(YY(2)-DLAM))/
     *                (YY(2)-YY(1))
            ELSE
             WX(I0X)=(PP(2)*(DLAM-YY(1))+PP(1)*(YY(2)-DLAM))/
     *                (YY(2)-YY(1))
             WX(I0X)=EXP(WX(I0X))
           END IF
           GO TO 200
  150      WX(I0X)=PRF0
  200   CONTINUE
        WZ(I0Z)=YINT(XX,WX,TL)
  300 CONTINUE
      W0=YINT(ZZ,WZ,ANEL)
      PHE1=W0*DOPL*1.772454
      RETURN
      END
 
C
C ********************************************************************
C
 
      SUBROUTINE HE2INI
C     =================
C
C     Initializes necessary arrays for evaluating the He II line
C     absorption profiles using data calculated by Schoening and
C     Butler
C
C     This procedure is quite analogous to HYDINI for hydrogen lines
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
      COMMON/HE2DAT/WL2(36,19),XT2(6),XNE2(11,19),PRF2(36,6,11),
     *              NWL2,NT2,NE2
      DATA NLINE1 /19/
C
      IH=67
      OPEN(UNIT=IH,FILE='./data/he2prf.dat',STATUS='OLD')
C
      DO ILINE=1,NLINE1
C
C     read the Schoening and Butler tables, which have to be stored
C     in file he23prf.dat
C
         READ(IH,501) ILHE2(ILINE),IUHE2(ILINE)
         IF(ILHE2(ILINE).LE.2) THEN
            WL00=227.838
          ELSE 
            WL00=227.7776
         END IF   
         WL0=WL00/(1./ILHE2(ILINE)**2-1./IUHE2(ILINE)**2)
         READ(IH,*) NWL2,(WL2(I,ILINE),I=1,NWL2)
         READ(IH,503) NT2,(XT2(I),I=1,NT2)
         READ(IH,504) NE2,(XNE2(I,ILINE),I=1,NE2)
         READ(IH,500)
         NWLHE2(ILINE)=NWL2
C
         DO I=1,NWL2
            IF(WL2(I,ILINE).LT.1.E-4) WL2(I,ILINE)=1.E-4
            WLHE2(ILINE,I)=LOG10(WL2(I,ILINE))
         END DO
C
         DO IE=1,NE2
            DO IT=1,NT2
               READ(IH,500)
               READ(IH,505) (PRF2(IWL,IT,IE),IWL=1,NWL2)
            END DO
         END DO
C
C        coefficient for the asymptotic profile is determined from
C        the input data
C
         XCLOG=PRF2(NWL2,1,1)+2.5*LOG10(WL2(NWL2,ILINE))+31.831-
     *         XNE2(1,ILINE)-2.*LOG10(WL0)
         XKLOG=0.6666667*(XCLOG-0.176)
         XK=EXP(XKLOG*2.3025851)
         DO ID=1,ND
            T=TEMP(ID)+2.42E-8*VTURB(ID)
            ANE=ELEC(ID)
            TL=LOG10(T)
            ANEL=LOG10(ANE)
            F00=1.25E-9*ANE**0.666666667
            FXK=F00*XK
            DOP=1.E8/WL0*SQRT(4.12E7*T)
            DBETA=WL0*WL0/2.997925E18/FXK
            BETAD=DBETA*DOP
C
C     interpolation to the actual values of temperature and electron
C     density. The result is stored at array PRFHE2, which has indices
C     ILINE  - index of line
C     ID     - depth index
C     IWL    - wavelength index (notice that the wavelength grid may
C              generally be different for different lines
C
            DO IWL=1,NWL2
               CALL INTHE2(PROF,TL,ANEL,IWL,ILINE)
               PRFHE2(ILINE,ID,IWL)=PROF
           END DO
        END DO
      END DO
      CLOSE(IH)
C
  500 FORMAT(1X)
  501 FORMAT(//14X,I2,9X,I2/)
c 502 FORMAT(2X,I4,1P6E10.3,4(/5X,0P6F10.4)/5X,5F10.4)
  503 FORMAT(2X,I4,F10.3,5F12.3)
  504 FORMAT(2X,I4,F10.2,5F12.2/4X,5F12.2)
  505 FORMAT(10F8.3)
      RETURN
      END
C
C ********************************************************************
C
C
 
      SUBROUTINE INTHE2(W0,X0,Z0,IWL,ILINE)
C     =====================================
C
C     Interpolation in temperature and electron density from the
C     Schoening and Butler tables for He II lines to the actual
C     actual values of temperature and electron density
C
C     This procedure is quite analogous to INTHYD for hydrogen lines
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (UN=1.)
      COMMON/HE2DAT/WL2(36,19),XT2(6),XNE2(11,19),PRF2(36,6,11),
     *              NWL2,NT2,NE2
      DIMENSION ZZ(3),XX(3),WX(3),WZ(3)
C
      NX=3
      NZ=3
C
C     for values lower than the lowest grid value of electron density
C     the profiles are determined by the approximate expression
C     (see STARKA); not by an extrapolation in the tables which may
C     be very inaccurate
C
      IF(Z0.LT.XNE2(1,ILINE)*0.99.OR.Z0.GT.XNE2(NE2,ILINE)*1.01) THEN
         CALL DIVHE2(A,DIV)
         W0=STARKA(WL2(IWL,ILINE)/FXK,A,DIV,UN)*DBETA
         W0=LOG10(W0)
         GO TO 500
      END IF
C
C     Otherwise, one interpolates (or extrapolates for higher than the
C     highes grid value of electron density) in the Schoening and
C     Butler tables
C
      DO 10 IZZ=1,NE2-1
         IPZ=IZZ
         IF(Z0.LE.XNE2(IZZ+1,ILINE)) GO TO 20
   10 CONTINUE
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE2-NZ+1) N0Z=NE2-NZ+1
      N1Z=N0Z+NZ-1
C
      DO 300 IZZ=N0Z,N1Z
         I0Z=IZZ-N0Z+1
         ZZ(I0Z)=XNE2(IZZ,iline)
C
C     Likewise, the approximate expression instead of extrapolation
C     is used for higher that the highest grid value of temperature,
C     if the Doppler width expressed in beta units (BETAD) is
C     sufficiently large (> 10)
C
         IF(X0.GT.1.01*XT2(NT2).AND.BETAD.GT.10.) THEN
            W0=STARKA(WL2(IWL,ILINE)/FXK,A,DIV,UN)*DBETA
            W0=LOG10(W0)
            GO TO 500
         END IF
C
C     Otherwise, normal inter- or extrapolation
C
C     Both interpolations (in T as well as in electron density) are
C     by default the quadratic interpolations in logarithms
C
         DO 30 IX=1,NT2-1
            IPX=IX
            IF(X0.LE.XT2(IX+1)) GO TO 40
   30    CONTINUE
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT2-NX+1) N0X=NT2-NX+1
         N1X=N0X+NX-1
         DO 200 IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XT2(IX)
            WX(I0)=PRF2(IWL,IX,IZZ)
  200       CONTINUE
         WZ(I0Z)=YINT(XX,WX,X0)
  300 CONTINUE
      W0=YINT(ZZ,WZ,Z0)
  500 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
 
      SUBROUTINE DIVHE2(A,DIV)
C     ========================
C
C     Auxiliary procedure for evaluating approximate Stark profile
C     for He II lines
C     This procedure is quite analogous to DIVSTR for hydrogen;
C     the only difference is a somewhat different definition
C     of the parameter A ,ie. A for He II is equal to A for hydrogen
C     minus ln(2)
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (UN=1.,TWO=2.,UNQ=1.25,UNH=1.5,TWH=2.5,FO=4.,FI=5.)
      PARAMETER (CA=0.978,BL=5.821,AL=1.26,CX=0.28,DX=0.0001)
C
      A=UNH*LOG(BETAD)-CA
      IF(BETAD.LT.BL) RETURN
      IF(A.GE.AL) THEN
         X=SQRT(A)*(UN+UNQ*LOG(A)/(FO*A-FI))
      ELSE
         X=SQRT(CX+A)
      ENDIF
      DO 10 I=1,5
         XN=X*(UN-(X*X-TWH*LOG(X)-A)/(TWO*X*X-TWH))
         IF(ABS(XN-X).LE.DX) GO TO 20
         X=XN
   10 CONTINUE
   20 DIV=X
      RETURN
      END
C
C ********************************************************************
C
C
 
      SUBROUTINE PHE2(ISPEC,ID,ABLIN,EMLIN)
C     =====================================
C
C     Evaluation of the opacity and emissivity in a given He II line,
C     using profile coefficients calculated by Schoening and Butler.
C
C     Input: ISPEC - line index, defined in HE2INI
C            ID    - depth index
C     Output: ABLIN - absorption coefficient
C             EMLIN - emission coefficient
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      DIMENSION ABLIN(1),EMLIN(1),OSCHE2(19),PRF0(40),WLL(40)
      COMMON/HE2PRF/PRFHE2(19,MDEPTH,36),WLHE2(19,36),NWLHE2(19),
     *              ILHE2(19),IUHE2(19)
      common/lasers/lasdel
      DATA OSCHE2/6.407E-1, 1.506E-1, 5.584E-2, 2.768E-2,
     *        1.604E-2, 1.023E-2, 6.980E-3,
     *        8.421E-1, 3.230E-2, 1.870E-2, 1.196E-2, 8.187E-3,
     *        5.886E-3, 4.393E-3, 3.375E-3, 2.656E-3,
     *        1.038,    1.793E-1, 6.549E-2/
C
C     ILINE - line index 
C
      ILINE=ISPEC-5
C
      DO 10 IWL=1,NWLHE2(ILINE)
         PRF0(IWL)=PRFHE2(ILINE,ID,IWL)
         WLL(IWL)=WLHE2(ILINE,IWL)
   10 CONTINUE
C
      I=ILHE2(ILINE)
      J=IUHE2(ILINE)
      II=I*I
      JJ=J*J
      IF(I.LE.2) THEN
         WLIN=227.838/(1./II-1./JJ)
       ELSE 
         WLIN=227.7776/(1./II-1./JJ)
      END IF   
      T=TEMP(ID)
C
C     He III population (either LTE or NLTE, depending on input model)
C
      IF(IELHE2.GT.0.and.inlte.gt.0) THEN
         PP=POPUL(NNEXT(IELHE2),ID)
         NLHE2=NLAST(IELHE2)-NFIRST(IELHE2)+1
       ELSE
         PP=RRR(ID,3,2)
         NLHE2=0
      END IF
C
C     population of the lower level of the given transition
C     (again either LTE or NLTE)
C
      PP=PP*ELEC(ID)*4.1412E-16/T/SQRT(T)*II
      IF(I.LE.NLHE2.and.inlte.gt.0) THEN
         POPI=POPUL(NFIRST(IELHE2)+I-1,ID)
       ELSE
         POPI=PP*EXP(631479./T/II)
      END IF
C
C     population of the upper level of the given transition
C     (again either LTE or NLTE)
C
      IF(J.LE.NLHE2) THEN
         POPJ=POPUL(NFIRST(IELHE2)+J-1,ID)*II/JJ
       ELSE
         POPJ=PP*EXP(631479./T/JJ)
      END IF
 
C
C     loop over frequency points - opacity and emissivity in the given line
C     absorption coefficent is found by interpolating in previously
C     calculated tables, based on calculations of Schoening and Butler
C     (see procedure HE2INI)
C
      FID=0.02654*OSCHE2(ILINE)
      DO 50 IJ=3,NFREQ
         AL=ABS(WLAM(IJ)-WLIN)
         IF(AL.LT.1.E-4) AL=1.E-4
         AL=LOG10(AL)
         DO 20 IWL=1,NWLHE2(ILINE)-1
            IW0=IWL
            IF(AL.LE.WLL(IWL+1)) GO TO 30
   20    CONTINUE
   30    IW1=IW0+1
         PRH=(PRF0(IW0)*(WLL(IW1)-AL)+PRF0(IW1)*(AL-WLL(IW0)))/
     *       (WLL(IW1)-WLL(IW0))
         SG=EXP(PRH*2.3025851)*FID
         if((popi-popj).le.0. .and. lasdel) goto 50
         ABLIN(IJ)=ABLIN(IJ)+SG*(POPI-POPJ)
         EMLIN(IJ)=EMLIN(IJ)+SG*POPJ*1.4747E-2*(FREQ(IJ)*1.E-15)**3
   50 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
 
      FUNCTION ISPEC(IAT,ION,ALAM)
C     ============================
C
C     Auxiliary procedure for INISET
C
C     Input:  IAT  - atomic number
C             ION  - ion (=1 for neutrals, =2 for once ionized, etc.)
C             ALAM - wavelength in nanometers
C     Output: ISPEC - parameter specifying whether the given line
C                     is taken with a special (pretabulated) absorption
C                     profile - only for hydrogen and helium
C                   = 0  - profile is taken as an ordinary Voigt profile
C                   > 0  - special profile
C
      INCLUDE 'PARAMS.FOR'
C
      ISPEC=0
      IF(IAT.GT.2) RETURN
C
      IF(IAT.EQ.1) THEN
         ISPEC=1
         RETURN
       ELSE
         IF(ION.EQ.1) THEN
            IF(ABS(ALAM-447.1).LT.0.5.AND.IHE1PR.GT.0) ISPEC=2
            IF(ABS(ALAM-438.8).LT.0.2.AND.IHE1PR.GT.0) ISPEC=3
            IF(ABS(ALAM-402.6).LT.0.2.AND.IHE1PR.GT.0) ISPEC=4
            IF(ABS(ALAM-492.2).LT.0.2.AND.IHE1PR.GT.0) ISPEC=5
          ELSE
C
            IF(ALAM.LT.163..OR.ALAM.GT.1012.7) RETURN
            IF(ALAM.LT.321.) THEN
               IF(ABS(ALAM-164.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=6
               IF(ABS(ALAM-320.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=7
               IF(ABS(ALAM-273.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=8
               IF(ABS(ALAM-251.1).LT.0.2.AND.IHE2PR.GT.0) ISPEC=9
               IF(ABS(ALAM-238.5).LT.0.2.AND.IHE2PR.GT.0) ISPEC=10
               IF(ABS(ALAM-230.6).LT.0.2.AND.IHE2PR.GT.0) ISPEC=11
               IF(ABS(ALAM-225.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=12
             ELSE IF(ALAM.LT.541.) THEN
               IF(ALAM.LT.392.3) RETURN
               IF(ABS(ALAM-468.6).LT.0.2.AND.IHE2PR.GT.0) ISPEC=13
               IF(ABS(ALAM-485.9).LT.0.2.AND.IHE2PR.GT.0) ISPEC=14
               IF(ABS(ALAM-454.2).LT.0.2.AND.IHE2PR.GT.0) ISPEC=15
               IF(ABS(ALAM-433.9).LT.0.2.AND.IHE2PR.GT.0) ISPEC=16
               IF(ABS(ALAM-420.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=17
               IF(ABS(ALAM-410.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=18
               IF(ABS(ALAM-402.6).LT.0.2.AND.IHE2PR.GT.0) ISPEC=19
               IF(ABS(ALAM-396.8).LT.0.2.AND.IHE2PR.GT.0) ISPEC=20
               IF(ABS(ALAM-392.3).LT.0.2.AND.IHE2PR.GT.0) ISPEC=21
             ELSE
               IF(ABS(ALAM-1012.4).LT.0.2.AND.IHE2PR.GT.0) ISPEC=22
               IF(ABS(ALAM-656.0).LT.0.2.AND.IHE2PR.GT.0) ISPEC=23
               IF(ABS(ALAM-541.2).LT.0.2.AND.IHE2PR.GT.0) ISPEC=24
            END IF
         END IF
      END IF
      RETURN
      END
C
C
C     ******************************************************************
C
C

      SUBROUTINE HESET(IL,ALM,EXCL,EXCU,ION,IPRF0,ILWN,IUPN)
C     ======================================================
C     
C     Auxiliary procedure for INISET - set up quantities:
C     IPRF0      - index for the procedure evaluating standard absorption
C                  profile coefficient for He I lines - see GAMHE
C     ILWN,IUPN  - only in NLTE option is switched on;
C                  indices of the lower and upper level associated with
C                  the given line
C
C     Input: IL - line index
C            ALM - line wavelength in nm
C            EXCL - excitation potential of the lower level (in cm**-1)
C            EXCU - excitation potential of the upper level (in cm**-1)
C            ION  - ionisation degree (1=neutrals, 2=once ionized, etc.)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION JU(24),NU(24),IT(24)
      DATA IT/1,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0/
      DATA NU/6,6,9,3,8,4,7,5,6,6,5,4,4,4,3,4,3,3,5,5,7,8,10,2/
      DATA JU/15,3,5,9,5,3,5,3,5,1,1,15,3,5,3,1,15,5,15,5,1,1,1,9/
C
C ******* He I  ***********
C
      IF(ION.NE.1) GO TO 20
C
C     switch IPRF0 - see GAMHE
C
      IL1=IL
      ALAM=ALM*10.
      IPRF=0
      IF(ABS(ALAM-3819.60).LT.1.) IPRF=1
      IF(ABS(ALAM-3867.50).LT.1.) IPRF=2
      IF(ABS(ALAM-3871.79).LT.1.) IPRF=3
      IF(ABS(ALAM-3888.65).LT.1.) IPRF=4
      IF(ABS(ALAM-3926.53).LT.1.) IPRF=5
      IF(ABS(ALAM-3964.73).LT.1.) IPRF=6
      IF(ABS(ALAM-4009.27).LT.1.) IPRF=7
      IF(ABS(ALAM-4120.80).LT.1.) IPRF=8
      IF(ABS(ALAM-4143.76).LT.1.) IPRF=9
      IF(ABS(ALAM-4168.97).LT.1.) IPRF=10
      IF(ABS(ALAM-4437.55).LT.1.) IPRF=11
      IF(ABS(ALAM-4471.50).LT.1.) IPRF=12
      IF(ABS(ALAM-4713.20).LT.1.) IPRF=13
      IF(ABS(ALAM-4921.93).LT.1.) IPRF=14
      IF(ABS(ALAM-5015.68).LT.1.) IPRF=15
      IF(ABS(ALAM-5047.74).LT.1.) IPRF=16
      IF(ABS(ALAM-5875.70).LT.1.) IPRF=17
      IF(ABS(ALAM-6678.15).LT.1.) IPRF=18
      IF(ABS(ALAM-4026.20).LT.1.) IPRF=19
      IF(ABS(ALAM-4387.93).LT.1.) IPRF=20
      IF(ABS(ALAM-4023.97).LT.1.) IPRF=21
      IF(ABS(ALAM-3935.91).LT.1.) IPRF=22
      IF(ABS(ALAM-3833.55).LT.1.) IPRF=23
      IF(ABS(ALAM-10830.0).LT.1.) IPRF=24
      IF(IPRF.GT.0.AND.IPRF.LE.20) IPRF0=IPRF
C
C     Indices of NLTE levels associated with the given line
C
      IF(INLTE.gt.5.OR.IELHE1.EQ.0) RETURN
      N0I=NFIRST(IELHE1)
      N1I=NLAST(IELHE1)
      HC=CL*H
      EION=ENION(N0I)/HC
      ILW=0
      IUN=0
      NQL=0
      IF(IPRF.GT.0) NQL=NU(IPRF)
      DO 10 I=N0I,N1I
         NQ=NQUANT(I)
         EX=EION-ENION(I)/HC
         IF(ABS(EXCL-EX).LT.100.) THEN
            ILW=I
            IGL=INT(G(I)+0.001)
         END IF
         IF(NQ.EQ.NQL) THEN 
            IG=INT(G(I)+0.001)
            IF(IT(IPRF).EQ.0) THEN
               IF(NQ.EQ.2.AND.IG.EQ.JU(IPRF)) IUN=I
               IF(NQ.EQ.3) THEN
                  IF(IG.EQ.JU(IPRF)) THEN
                     IF(IG.EQ.1.OR.IG.EQ.5) IUN=I
                     IF(IG.EQ.3.AND.IGL.EQ.1) IUN=I
                   ELSE
                     IF(IG.EQ.9) IUN=I
                  END IF
               END IF
               IF(NQ.EQ.4) THEN
                  IF(IG.EQ.JU(IPRF)) THEN
                     IF(IG.EQ.1.OR.IG.EQ.5.OR.IG.EQ.7) IUN=I
                     IF(IG.EQ.3.AND.IGL.EQ.1) IUN=I
                   ELSE
                     IF(IG.EQ.16) IUN=I
                  END IF
               END IF
               IF(IG.EQ.25.OR.IG.EQ.36) IUN=I
               IF(IG.EQ.49.OR.IG.EQ.64.OR.IG.EQ.81) IUN=I
               IF(IG.EQ.100.OR.IG.EQ.121.OR.IG.EQ.144) IUN=I
             ELSE
               IF(NQ.EQ.3) THEN
                  IF(IG.EQ.JU(IPRF)) THEN
                     IF(IG.EQ.9.OR.IG.EQ.15) IUN=I
                     IF(IG.EQ.3.AND.IGL.EQ.9) IUN=I
                   ELSE
                     IF(IG.EQ.27) IUN=I
                  END IF
               END IF
               IF(NQ.EQ.4) THEN
                  IF(IG.EQ.JU(IPRF)) THEN
                     IF(IG.EQ.9.OR.IG.EQ.15.OR.IG.EQ.21) IUN=I
                     IF(IG.EQ.3.AND.IGL.EQ.9) IUN=I
                   ELSE
                     IF(IG.EQ.48) IUN=I
                  END IF
               END IF
               IF(IG.EQ.75) IUN=I
               IF(IG.EQ.108.OR.IG.EQ.147.OR.IG.EQ.192) IUN=I
               IF(IG.EQ.243.OR.IG.EQ.300.OR.IG.EQ.363) IUN=I
            END IF
            IF(NQ.EQ.2.AND.IG.EQ.16) IUN=I
            IF(NQ.EQ.3.AND.IG.EQ.36) IUN=I
            IF(NQ.EQ.4.AND.IG.EQ.64) IUN=I
            IF(NQ.EQ.5.AND.IG.EQ.100) IUN=I
            IF(NQ.EQ.6.AND.IG.EQ.144) IUN=I
            IF(NQ.EQ.7.AND.IG.EQ.196) IUN=I
            IF(NQ.EQ.8.AND.IG.EQ.256) IUN=I
            IF(NQ.EQ.9.AND.IG.EQ.324) IUN=I
            IF(NQ.EQ.10.AND.IG.EQ.400) IUN=I
         END IF 
   10 CONTINUE
c     print *, 'il,iprof,ilw,iupn',il,iprf,ilw,iun
      ILWN=ILW
      IUPN=IUN
C
C ******* He II ***********
C
   20 IF(ION.NE.2.OR.IELHE2.LE.0) RETURN 
      N0I=NFIRST(IELHE2)
      NLHE2=NLAST(IELHE2)-N0I+1
      XL=SQRT(1./(1.-EXCL/438916.146))
      ILW=INT(XL)
      IF((FLOAT(ILW)-XL).LT.0.) ILW=ILW+1
      XU=SQRT(1./(1.-EXCU/438916.146))
      IUN=INT(XU)
      IF((FLOAT(IUN)-XU).LT.0.) IUN=IUN+1
      IF(ILW.LE.NLHE2) ILWN=ILW+N0I-1
      IF(IUN.LE.NLHE2) IUPN=IUN+N0I-1
      RETURN
      END
C
C
C ********************************************************************
C
      SUBROUTINE INISET
C     =================
C
C     SELECTION OF LINES THAT MAY CONTRIBUTE,
C     SET UP AUXILIARY FIELDS CONTAINING LINE PARAMETERS,
C     SET UP THE SET OF FREQUENCY POINTS
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'WINCOM.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      COMMON/CTRFUN/CINT1(MDEPTH),CINT2(MDEPTH),
     *   CTRI(MDEPTH),CTRR(MDEPTH),XKAR(MDEPTH),
     *   ABXLI(MFREQ),EMXLI(MFREQ),IJCTR(MFREQ)
      SAVE ILLAST
C
      DATA CNM,CAS /2.997925D17,2.997925D18/
c      DATA C1,C2,C3 /2.3025851, 4.2014672, 1.4387886/
C
      DO 10 I=1,MFRQ
         W(I)=0.
         IJCTR(I)=0
   10 CONTINUE
C
      IL0=0
      IPRSET=0
      NLIN=0
      IREADP=1
      IRLIST=0
      IF(IBLANK.LE.1.OR.IMODE.EQ.1.OR.IMODE.EQ.-1) IREADP=0
      IF(IBLANK.LE.1) APREV=0.
      FRMIN=CNM/ALAM0
      FRM=FRMIN
      if(ifwin.le.0) then
         ij0=3
       else
         ij0=1
      end if
      IJ=IJ0
      FREQ(IJ0)=FRM
      SPACE=SPACE0
      IF(ALAMC.GT.0.) SPACE=SPACE0*ALAM0/ALAMC
      IF(SPACE0.LT.0.) SPACE=-SPACE0
      IF(IMODE.EQ.2) THEN
         NFRP=NFREQS+1
         W0=SPACE
         GO TO 105
      END IF
C
      ISTR=0
      IJMAX=0
      IMOD1L=0
      if(ifwin.le.0) then
      CUTOFF=CUTOF0
      DOPSTD=1.E7/ALAM0*DSTD
      DISTAN=0.15*DOPSTD
      SPAC=3.E16/ALAM0/ALAM0*SPACE
      DISTA0=0.14*SPAC
      ASTD=1.0
      AVAB=ABSTD(IDSTD)*RELOP
      end if
      FRLI0=FRMIN
      IF(IBLANK.GE.2.AND.IMODE.EQ.-1) IL0=ILLAST
C
   20 CONTINUE 
C
C     set up indices of lines
C     IL0 - is the current index of line in the numbering of all lines
C
      IF(IREADP.EQ.1) THEN
         IPRSET=IPRSET+1
         IL0=INDLIP(IPRSET)
         IF(FREQ0(IL0).LT.FRMIN) THEN
            IREADP=0
            IL0=INDLIP(IPRSET-1)+1
         END IF
       ELSE
         IL0=IL0+1
      END IF
      IF(IL0.GT.NLIN0) GO TO 210
      FRLIM=FRLI0
      FR0=FREQ0(IL0)
      ALAM=CNM/FR0
C  
      if(ifwin.gt.0) then
      IF(ALAMC.GT.0.) SPACE=SPACE0*ALAM/ALAMC
      IF(SPACE0.LT.0.) SPACE=-SPACE0
      CUTOFF=CUTOF0*ALAM/ALAMC
      DOPSTD=1.E7/ALAM*DSTD
      DISTAN=0.15*DOPSTD
      SPAC=SPACE
      IF(MOD(IFREQ,10).GT.0) SPAC=3.E16/ALAM/ALAM*SPACE
      DISTA0=0.14*SPAC
      end if
C
C     set up a different starting wavelength for IMODE=1
C
      IF(IMODE.NE.1) GO TO 45
      IF(ISTR.EQ.1.OR.IJ.NE.3) GO TO 45
      IF(ALAM.LT.ALAM0+2.*CUTOFF) GO TO 45
      ALAM0=ALAM-CUTOFF+0.0001
      FRMIN=CNM/ALAM0
      FRM=FRMIN
      IJ=IJ0
      FREQ(IJ0)=FRM
   45 CONTINUE
      IF(ALAM.LT.ALAM0-CUTOFF) GO TO 20
      IF(IJ.LT.NFREQS+1) GO TO 50
      IF(ALAM.GT.ALAM1+CUTOFF) GO TO 210
C
C     SECOND SELECTION : FOR LINE STRENGHTS
C
   50 CONTINUE
      ISTR=0
      IF(IMODE.GE.1) THEN
         ISTR=1
       ELSE
         EXT=EXTIN(IL0)
         FRLI0=FR0-EXT-SPAC
         IF(FRLI0.GT.FRLIM) FRLI0=FRLIM
         frmiv=frmin
         if(ifwin.gt.0) frmiv=frmiv*(1.+vinf/2.997925e10)
         IF(ALAM.LT.ALAM0.AND.FR0-FRMIv.GT.EXT+SPAC) GO TO 20
         ISTR=1
         frmav=frmax
         if(ifwin.gt.0) frmav=frmav*(1.-vinf/2.997925e10)
         IF(IJ.GE.NFREQS+1.AND.FRMAv-FR0.GT.EXT+SPAC) GO TO 20
      END IF

C
      NLIN=NLIN+1
      if(nlin.gt.mlin) call quit(' too many lines in a set')
      INDLIN(NLIN)=IL0
      ALAMCU=ALAM+CUTOFF
C
C     FREQUENCY POINTS AND WEIGHTS
C
      IF(IJ.GE.NFREQS+1) GO TO 20
      IF(FR0.GT.FRMIN) GO TO 20
  100 DELT=ABS(FRM-FR0)
      IF(DELT.LT.DISTA0.AND.IMODE.NE.1) GO TO 20
      DFREL=CNM*(1.D0/FR0-1.D0/FRM)/SPACE
      NFRP=int(DFREL)+1
      IF(NFRP.LE.2) NFRP=2
      W0=CNM*(1.D0/FR0-1.D0/FRM)/NFRP
      FRM=FR0
  105 FRACT=FREQ(IJ)
      ALACT=CNM/FRACT
C
      DO 110 K=1,NFRP
         FRACT=FRACT-W0
         ALACT=ALACT+W0
         IF(IMODE.GE.1.OR.NFRP.EQ.2) GO TO 107
         IF(FRACT.LT.FRLIM.AND.FRACT.GT.FR0+EXT+SPAC) GO TO 110
  107    IJ=IJ+1
         IF(IJ.GT.NFREQS) GO TO 130
         FREQ(IJ)=CNM/ALACT
         W(IJ)=W(IJ)+(FREQ(IJ-1)-FREQ(IJ))*0.5
         W(IJ-1)=W(IJ-1)+(FREQ(IJ-1)-FREQ(IJ))*0.5
C        IF(FREQ(IJ).LT.FRLAST) GO TO 220
         IF(IMODE.EQ.1.AND.ALACT.GT.ALAMCU) GO TO 140
  110 CONTINUE
      IJCTR(IJ)=IL0
      IF(IMOD1L.EQ.1) GO TO 210
      DISTA0=DISTAN
      GO TO 20
C
  130 FRMAX=FREQ(NFREQS)
      ALAM1=CNM/FRMAX
      NFREQ=NFREQS
      IF(IMODE.EQ.2) GO TO 210
      IF(IMOD1L.EQ.1) GO TO 210
      GO TO 20
C
  140 IJMAX=IJ
      IJMAX=MIN(IJMAX,NFREQS)
      NFREQ=IJMAX
      IF(IL0.LT.NLIN0) THEN
         NBLANK=IBLANK+1
       ELSE
         NBLANK=IBLANK
      END IF
      GO TO 240
C
  210 NBLANK=IBLANK+1
      IF(IJ.GE.NFREQS+1) GO TO 230
      IJMAX=IJ
      IJMAX=MIN(IJMAX,NFREQS)
      NFREQ=IJMAX
      IF(IMODE.NE.1) GO TO 240
      IF(IMOD1L.EQ.1) GO TO 240
C     FR0=MAX(CNM/(ALAM+CUTOFF),FRLAST*0.99999999D0)
      FR0=FRLAST*0.99999999D0
      ALAM=CNM/FR0
      IMOD1L=1
      GO TO 100
C
  230 IJMAX=NFREQS
      NFREQ=NFREQS
  240 IF(FREQ(IJMAX).LE.FRLAST) NBLANK=IBLANK
      if(alm00.gt.0.) then
         if(freq(ijmax).ge.0.999999*cnm/alm00.and.iblank.gt.1)
     *      nblank=iblank
      end if
c
c     correction for molecular lines
c
      if(nmlist.gt.0.and.ifmol.gt.0) then
      do ilist=1,nmlist
         if(alastm(ilist).gt.0..and.alastm(ilist).le.alact) then
          nblank=iblank
          irlist=1
c         write(*,*) 'iniset mol',ilist,alastm(ilist),alam
         end if
      end do
      end if
c
      if(ifwin.le.0) then
      FREQ(1)=FREQ(3)
      FREQ(2)=FREQ(IJMAX)
      W(1)=0.5*(FREQ(1)-FREQ(2))
      W(2)=W(1)
      end if
C
C     truncate the interval if the required end is reached
C
      ijmx=2
      if(ifwin.gt.0) ijmx=ijmax
      IF(FREQ(ijmx).LT.FRLAST) THEN
         FREQ(ijmx)=FRLAST
         if(ifwin.le.0) then
         W(1)=0.5*(FREQ(1)-FREQ(2))
         W(2)=W(1)
         end if
         DO 245 IJ=IJ0,NFREQ
            IF(FREQ(IJ).LT.FRLAST) GO TO 247
            IJMAX=IJ
  245    CONTINUE
  247    NFREQ=IJMAX+1
         FREQ(NFREQ)=FRLAST
         W(NFREQ)=0.5*(FREQ(NFREQ-1)-FREQ(NFREQ))
         W(NFREQ-1)=W(NFREQ)+0.5*(FREQ(NFREQ-2)-FREQ(NFREQ-1))
      END IF
C
C     frequency interpolation coefficients
C
      IF(IMODE.NE.-1) THEN
      if(ifwin.le.0) then
      XX=FREQ(2)-FREQ(1)
      DO IJ=1,NFREQ
         WLAM(IJ)=2.997925E18/FREQ(IJ)
         FRX1(IJ)=(FREQ(IJ)-FREQ(1))/XX
         FRX2(IJ)=(FREQ(2)-FREQ(IJ))/XX
      END DO
      else
      DO IJ=1,NFREQ
         WLAM(IJ)=CAS/FREQ(IJ)
         frqobs(ij)=freq(ij)
         wlobs(ij)=wlam(ij)
         fr=freq(ij)
         BNUE(IJ)=BN*fr*fr*fr
         DO IJCI=1,NFREQC-1
            IF(WLAM(IJ).LE.WLAMC(IJCI)) GO TO 248
         END DO
  248    CONTINUE
         IJC=IJCI
         IJCINT(IJ)=MAX(IJC-1,1)
         IJCI=IJCINT(IJ)
         FRX1(IJ)=(FREQ(IJ)-FREQC(IJCI+1))/
     *            (FREQC(IJCI)-FREQC(IJCI+1))
      END DO
      nfrobs=nfreq
      xx=freq(nfreq)-freq(1)
      end if
c
c     frequency indices of the line centers
c
      DFRCON=NFREQ-ij0
      DFRCON=-DFRCON/XX
      IFRCON=INT(DFRCON)
      DO 255 IL=1,NLIN
         fr0=freq0(indlin(il))
         XJC=3.+DFRCON*(FREQ(1)-FR0)
         IJC=INT(XJC)
         IJCNTR(IL)=IJC
         if(ijc.le.ij0.or.ijc.ge.nfreq) go to 255
         if(fr0.lt.freq(ijc)) then
            ijc0=ijc
            dfr0=freq(ijc0)-fr0
  252       ijc0=ijc0+1
            dfr=abs(freq(ijc0)-fr0)
            if(dfr.lt.dfr0) then
               ijc=ijc0
               ijc0=ijc0+1
               dfr0=dfr
               go to 252
            end if
          else if(fr0.gt.freq(ijc)) then
            ijc0=ijc
            dfr0=fr0-freq(ijc0)
  254       ijc0=ijc0-1
            dfr=abs(freq(ijc0)-fr0)
            if(dfr.lt.dfr0) then
               ijc=ijc0
               ijc0=ijc0-1
               dfr0=dfr
               go to 254
            end if
         end if
         IJCNTR(IL)=IJC
  255 continue
      END IF
C
      if(ifwin.gt.0) then
C
c     set up switches for hydrogen and He II line opacity
c
      DO  IJ=1,NFREQ
         call hylsew(ij)
         call he2sew(ij)
      end do
      end if
C
      NSP=0
      DO 260 IL=1,NLIN
         IL0=INDLIN(IL)
         ISP=ISPRF(IL0)
         IF(ISP.GT.5) THEN
            NSP=NSP+1
            ISP0(NSP)=ISP
         END IF
         INDLIP(IL)=INDLIN(IL)
  260 CONTINUE
      if(ifwin.le.0) then
      ILLAST=INDLIN(NLIN)
      else
      ILLAST=0
      IF(NLIN.GT.0) ILLAST=INDLIN(NLIN)
      end if
C
      CALL READPH
C
      IF(ALAM0.LE.APREV+0.001) NBLANK=IBLANK
      APREV=ALAM0
      ALAM0=ALAM1
      ALM00=CNM/FREQ(NFREQ)
c
c      write(6,611) iblank,nblank,irlist,aprev*10.,alam0*10.
c 611  format('inis ',2i6,i3,3f10.3)
      RETURN
      END
C
C ********************************************************************
C
C
      SUBROUTINE READPH
C     =================
C
C     Auxiliary routine for LINSET - read table of detailed
C     photoinization cross-section from unit IPHT1, 
C     and interpolate to the set of current wavelengths (WLAM)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/PHOTCS/PHOT(MFRQ,MPHOT),WPHT0,WPHT1,APHT(MPHOT),
     *              EPHT(MPHOT),GPHT(MPHOT),JPHT(MPHOT),
     *              NPHT
      DIMENSION PHT0(MPHOT),PHT1(MPHOT),IPHT(MPHOT),IEND(MPHOT),
     *          IFILE(MPHOT),NELEM(MPHOT),INDEX(MPHOT,MPHOT)
      PARAMETER (IPHT0=57)
      SAVE IPHT,IEND,NELEM,INDEX,NUMFIL
C
C     initialization - read basic information about files where the
C                      cross-sections are stored, 
C                      and basic parameters for starting levels
C
      IF(IBLANK.LE.1) THEN
         NPHT=0
         IPHT1=0
         NUMFIL=0
         DO 10 IJ=1,MFRQ
            DO 10 I=1,MPHOT
   10          PHOT(IJ,I)=0.
         READ(IPHT0,*,END=50,err=50) NPHT
         IF(NPHT.LE.0) RETURN
         npht1=npht
         READ(IPHT0,*,END=50) (IPHT(I),I=1,NPHT)
         READ(IPHT0,*,END=50) (APHT(I),I=1,NPHT)
         READ(IPHT0,*,END=50) (EPHT(I),I=1,NPHT)
         READ(IPHT0,*,END=50) (GPHT(I),I=1,NPHT)
         READ(IPHT0,*,END=50) (JPHT(I),I=1,NPHT)
C
C     determination of the number of files (NFILE) and the
C     partitioning of the individual cross-section to the corresponding
C     files
C
         NUMFIL=1
         IFILE(1)=1
         NELEM(1)=1
         INDEX(1,1)=1
         IF(NPHT.GT.1) THEN
            DO 30 I=2,NPHT
               DO 20 J=1,I-1
                  IF(IPHT(I).EQ.IPHT(J)) THEN
                     IFILE(I)=IFILE(J)
                     NELEM(IFILE(I))=NELEM(IFILE(I))+1
                     INDEX(IFILE(I),NELEM(IFILE(I)))=I
                     GO TO 30
                  END IF
   20          CONTINUE
               NUMFIL=NUMFIL+1
               IFILE(I)=NUMFIL
               NELEM(NUMFIL)=1
               INDEX(NUMFIL,1)=I
   30       CONTINUE
         END IF
         DO 40 IFIL=1,NUMFIL
            IEND(IFIL)=0
   40    CONTINUE
      END IF
   50 IF(NUMFIL.LE.0) RETURN
c
C     loop over individual files containing the photoionization data
C
      DO 300 IFIL=1,NUMFIL
         IF(IEND(IFIL).EQ.1) GO TO 200
         IF(IEND(IFIL).EQ.2) GO TO 300
         NPHT1=NELEM(IFIL)
         IPHT1=IPHT(INDEX(IFIL,1))
         IF(IBLANK.LE.1) THEN
  110       READ(IPHT1,*,END=200) WPHT1,(PHT1(I),I=1,NPHT1)
            IF(WPHT1.LT.WLAM(1)) GO TO 110
            BACKSPACE(IPHT1)
            BACKSPACE(IPHT1)
            READ(IPHT1,*,END=200) WPHT0,(PHT0(I),I=1,NPHT1)
          ELSE
            BACKSPACE(IPHT1)
            BACKSPACE(IPHT1)
            READ(IPHT1,*,END=200) WPHT0,(PHT0(I),I=1,NPHT1)
            READ(IPHT1,*,END=200) WPHT1,(PHT1(I),I=1,NPHT1)
         END IF
         DW=WPHT1-WPHT0
         A1=(WPHT1-WLAM(3))/DW
         A2=(WLAM(3)-WPHT0)/DW
         DO 130 I=1,NPHT1
            INDX=INDEX(IFIL,I)
            PHOT(1,INDX)=0.
            PHOT(2,INDX)=0.
            PHOT(3,INDX)=(A1*PHT0(I)+A2*PHT1(I))*1.E-18
            DO 130 IJ=4,MFRQ
               PHOT(IJ,INDX)=0.
  130    CONTINUE
         DO 190 IJ=4,MFRQ
            IF(WLAM(IJ).LE.WPHT1) THEN
               A1=(WPHT1-WLAM(IJ))/DW
               A2=(WLAM(IJ)-WPHT0)/DW
               DO 140 I=1,NPHT1
                  INDX=INDEX(IFIL,I)
                  PHOT(IJ,INDX)=(A1*PHT0(I)+A2*PHT1(I))*1.E-18
  140          CONTINUE
             ELSE
               WPHT0=WPHT1
               DO 150 I=1,NPHT1
  150             PHT0(I)=PHT1(I)
               IFSML=0
  160          READ(IPHT1,*,END=180) WPHT1,(PHT1(I),I=1,NPHT1)
               IF(WPHT1.LT.WLAM(IJ)) THEN
                  IFSML=1
                  GO TO 160
               END IF
               IF(IFSML.EQ.1) THEN
                  BACKSPACE(IPHT1)           
                  BACKSPACE(IPHT1)           
                  READ(IPHT1,*,END=180) WPHT0,(PHT0(I),I=1,NPHT1)
                  READ(IPHT1,*,END=180) WPHT1,(PHT1(I),I=1,NPHT1)
               END IF
               DW=WPHT1-WPHT0
               A1=(WPHT1-WLAM(IJ))/DW
               A2=(WLAM(IJ)-WPHT0)/DW
               DO 170 I=1,NPHT1
                  INDX=INDEX(IFIL,I)
                  PHOT(IJ,INDX)=(A1*PHT0(I)+A2*PHT1(I))*1.E-18
  170          CONTINUE
            END IF
            GO TO 190
  180       IEND(IFIL)=1      
            DO 185 I=1,NPHT1
               INDX=INDEX(IFIL,I)
               PHOT(IJ,INDX)=0.
  185       CONTINUE
  190    CONTINUE
         PHOT(1,INDX)=PHOT(3,INDX)
         PHOT(2,INDX)=PHOT(MFRQ,INDX)
         GO TO 300
  200    IEND(IFIL)=2
         DO 210 IJ=1,MFREQ
            DO 210 I=1,NELEM(IFIL)
               INDX=INDEX(IFIL,I)
               PHOT(IJ,INDX)=0.
  210    CONTINUE
  300 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
      SUBROUTINE INILIN
C     =================
C
C     read in the input line list,
C     selection of lines that may contribute,
C     set up auxiliary fields containing line parameters,
C
C     Input of line data - unit 19:
C
C     For each line, one (or two) records, containing:
C
C    ALAM    - wavelength (in nm)
C    ANUM    - code of the element and ion (as in Kurucz-Peytremann)
C              (eg. 2.00 = HeI; 26.00 = FeI; 26.01 = FeII; 6.03 = C IV)
C    GF      - log gf
C    EXCL    - excitation potential of the lower level (in cm*-1)
C    QL      - the J quantum number of the lower level
C    EXCU    - excitation potential of the upper level (in cm*-1)
C    QU      - the J quantum number of the upper level
C    AGAM    = 0. - radiation damping taken classical
C            > 0. - the value of Gamma(rad)
C
C     There are now two possibilities, called NEW and OLD, of the next
C     parameters:
C     a) NEW, next parameters are:
C    GS      = 0. - Stark broadening taken classical
C            > 0. - value of log gamma(Stark)
C    GW      = 0. - Van der Waals broadening taken classical
C            > 0. - value of log gamma(VdW)
C    INEXT   = 0  - no other record necessary for a given line
C            > 0  - a second record is present, see below
C
C    The following parameters may or may not be present,
C    in the same line, next to INEXT:
C    ISQL   >= 0  - value for the spin quantum number (2S+1) of lower level
C            < 0  - value for the spin number of the lower level unknown
C    ILQL   >= 0  - value for the L quantum number of lower level
C            < 0  - value for L of the lower level unknown
C    IPQL   >= 0  - value for the parity of lower level
C            < 0  - value for the parity of the lower level unknown
C    ISQU   >= 0  - value for the spin quantum number (2S+1) of upper level
C            < 0  - value for the spin number of the upper level unknown
C    ILQU   >= 0  - value for the L quantum number of upper level
C            < 0  - value for L of the upper level unknown
C    IPQU   >= 0  - value for the parity of upper level
C            < 0  - value for the parity of the upper level unknown
C    (by default, the program finds out whether these quantum numbers
C     are included, but the user can force the program to ignore them
C     if present by setting INLIST=11)
C
C    If INEXT was set to >0 then the following record includes:
C    WGR1,WGR2,WGR3,WGR4 - Stark broadening values from Griem (in Angst)
C                   for T=5000,10000,20000,40000 K, respectively;
C                   and n(el)=1e16 for neutrals, =1e17 for ions.
C    ILWN    = 0  - line taken in LTE (default)
C            > 0  - line taken in NLTE, ILWN is then index of the
C                   lower level
C            =-1  - line taken in approx. NLTE, with Doppler K2 function
C            =-2  - line taken in approx. NLTE, with Lorentz K2 function
C    IUN     = 0  - population of the upper level in LTE (default)
C            > 0  - index of the lower level
C    IPRF    = 0  - Stark broadening determined by GS
C            < 0  - Stark broadening determined by WGR1 - WGR4
C            > 0  - index for a special evaluation of the Stark
C                   broadening (in the present version inly for He I -
C                   see procedure GAMHE)
C      b) OLD, next parameters are
C     IPRF,ILWN,IUN - the same meaning as above
C     next record with WGR1-WGR4 - again the same meaning as above
C     (this record is automatically read if IPRF<0
C
C     The only differences between NEW and OLD is the occurence of
C     GS and GW in NEW, and slightly different format of reading.
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      COMMON/IPOTLS/IPOTL(mlin0)
C
      PARAMETER (C1     = 2.3025851,
     *           C2     = 4.2014672,
     *           C3     = 1.4387886,
     *           CNM    = 2.997925D17,
     *           ANUMIN = 1.9,
     *           ANUMAX = 99.31,
     *           AHE2   = 2.01,
     *           EXT0   = 3.17,
     *           UN     = 1.0,
     *           TEN    = 10.,
     *           HUND   = 1.D2,
     *           TENM4  = 1.D-4,
     *           TENM8  = 1.D-8,
     *           OP4    = 0.4,
     *           AGR0=2.4734E-22, 
     *           XEH=13.595, XET=8067.6, XNF=25.,
     *           R02=2.5, R12=45., VW0=4.5E-9)
      PARAMETER (ENHE1=198310.76, ENHE2=438908.85)
      CHARACTER*1000 CADENA

      DATA INLSET /0/
C
      if(ibin(0).eq.0) then
         open(unit=19,file=amlist(0),status='old')
       else
         open(unit=19,file=amlist(0),form='unformatted',status='old')
      end if
      if(imode.lt.-2) then
         call inilin_grid
         return
      end if
c
      if(ndstep.eq.0) then
         write(6,621) idstd,temp(idstd),dens(idstd)
       else
         write(6,622) 
         do id=1,nd,ndstep
            write(6,623) id,temp(id),dens(id)
         end do
      end if
  621 format(/' lines are rejected based on opacities at the',
     * ' standard depth:'/
     * ' ID =',i4,'  T = ',f10.1,',   DENS = ',1pe10.3/)
  622 format(/' lines are rejected based on opacities at depths:'/)
  623 format(' ID =',i4,'  T = ',f10.1,',   DENS = ',1pe10.3/)
c
      IL=0
      INNLT0=0
      IGRIE0=0
      IF(NXTSET.EQ.1) THEN
          ALAM0=ALM00
          ALAST=ALST00
          FRLAST=CNM/ALAST
          NXTSET=0
          REWIND 19
      END IF
      ALAM00=ALAM0
      ALAST=CNM/FRLAST
      ALAST0=ALAST
      DOPSTD=1.E7/ALAM0*DSTD
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD
      AVAB=ABSTD(IDSTD)*RELOP
      ASTD=1.0
c     IF(GRAV.GT.6.) ASTD=0.1
      CUTOFF=CUTOF0
      ALAST=CNM/FRLAST
      IF(INLTE.GE.1.AND.INLSET.EQ.0) THEN
         CALL NLTSET(0,IL,IAT,ION,ALAM0,EXCL,EXCU,QL,QU,
     *         ISQL,ILQL,IPQL,ISQU,ILQU,IPQU,IEVEN,INNLT0,ILMATCH)
         INLSET=1
         ILMATCH=0
         ILSEARCH=0
         ILFOUND=0
         ILFAIL=0
         ILMULT=0
      END IF
c
C
C    Check whether any ion needs to compare quantum number limits
C    
      MAXILIMITS=0
      DO I=1,NION
        IF (ILIMITS(I).EQ.1) MAXILIMITS=1
      ENDDO
      IF (MAXILIMITS.EQ.0.and.inlist.gt.0) INLIST=1
C
C     If INLIST=0 or 10, the program checks for the number of words
C     present in the first line of the file to determine if quantum
C     numbers are included. If  INLINST=11, it means quantum numbers 
C     are present

      IF (INLIST.EQ.0.or.INLIST.EQ.10) THEN
        CADENA=' '
        READ(19,'(1000a)')CADENA
        BACKSPACE(19)
        CALL COUNT_WORDS(CADENA,NOW)
        IF(NOW.LT.12) THEN
           INLIST=0
           WRITE(11,*) 'INILIN: NO quantum numbers given in linelist'
         ELSE      
           INLIST=10
           WRITE(11,*) 'INILIN: quantum numbers included in linelist'
        ENDIF
      ELSE IF (INLIST.EQ.11) THEN
C        INLIST=1
        WRITE(11,*)'INILIN: quant. num. limits are present and used'
      ELSE
        WRITE(11,*)'INILIN: if present, quant. num. limits are ignored'
      ENDIF

      rstd=1.e4
      if(relop.gt.0.) rstd=1./relop
      afac=10.
      if(iat.gt.15.and.iat.ne.26) afac=1.
      afac=afac*rstd*astd
C
C     first part of reading line list - read only lambda, and
C     skip all lines with wavelength below ALAM0-CUTOFF
C
      ALAM=0.
      IJC=2
    7 if(ibin(0).eq.0) then
         READ(19,510) ALAM
       else
         read(19) alam
      end if
  510 FORMAT(F10.4)
      IF(ALAM.LT.ALAM0-CUTOFF) GO TO 7
      BACKSPACE(19)
      GO TO 10
c
c   8 write(6,688) alam
c 688 format(' error in line list, lambda= ',f12.4/)
    8 continue
   10 ILWN=0
      IUN=0
      IPRF=0
      GS=0.
      GW=0.
      IF(IBIN(0).EQ.0) THEN
         IF(INLIST.EQ.0) THEN
            READ(19,*,END=100,err=8) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
     *                           GS,GW,INEXT
            IF(INEXT.NE.0) READ(19,*) WGR1,WGR2,WGR3,WGR4,ILWN,IUN,IPRF
          ELSE IF(INLIST.GE.10) THEN
            READ(19,*,END=100,err=8) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
     *                           GS,GW,INEXT,ISQL,ILQL,IPQL,ISQU,ILQU,IPQU
         END IF
       ELSE
         IF(INLIST.LT.10) THEN
            READ(19,END=100) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,GS,GW
c    *                       GS,GW,INEXT
          ELSE IF(INLIST.GT.10) THEN
            READ(19,END=100) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,GS,GW,
     *                       INEXT,ISQL,ILQL,IPQL,ISQU,ILQU,IPQU
         END IF
      END IF
      IF(INLIST.EQ.10.OR.INLIST.EQ.11) THEN
         IF(ISPICK.EQ.0) THEN
                ISQL=-1
                ISQU=-1
         END IF
         IF(ILPICK.EQ.0) THEN
                ILQL=-1
                ILQU=-1
         END IF
         IF(IPPICK.EQ.0) THEN
                IPQL=-1
                IPQU=-1
         END IF

         IF(INEXT.NE.0) READ(19,*) WGR1,WGR2,WGR3,WGR4,ILWN,IUN,IPRF
      END IF
C
c     change wavelength to vacuum for lambda > 2000
c
      if(alam.gt.200..and.vaclim.gt.2000.) then
         wl0=alam*10.
         ALM=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0*(XN1*1.D-6+UN)
         alam=wl0*0.1
      END IF        
C
C     first selection : for a given interval a atomic number
C
      IF(ALAM.GT.ALAST+CUTOFF) GO TO 100
      IF(ANUM.LT.ANUMIN.OR.ANUM.GT.ANUMAX) GO TO 10
      IF(ABS(ANUM-AHE2).LT.TENM4.AND.IFHE2.GT.0) GO TO 10
C
C     second selection : for line strenghts
C
      FR0=CNM/ALAM
      IAT=INT(ANUM)
      FRA=(ANUM-FLOAT(IAT)+TENM4)*HUND
      ION=INT(FRA)+1
      IF(ION.GT.IONIZ(IAT)) GO TO 10
      IEVEN=1
      EXCL=ABS(EXCL)
      EXCU=ABS(EXCU)
      IF(EXCL.GT.EXCU) THEN
         FRA=EXCL
         EXCL=EXCU
         EXCU=FRA
         FRA=QL
         QL=QU
         QU=FRA
         IEVEN=0
         IF(INLIST.GE.10) THEN
            IFRA=ISQL
            ISQL=ISQU
            ISQU=IFRA

            IFRA=ILQL
            ILQL=ILQU
            ILQU=IFRA

            IFRA=IPQL
            IPQL=IPQU
            IPQU=IFRA
         END IF
      END IF
      GFP=C1*GF-C2
      EPP=C3*EXCL
c  
      if(ndstep.eq.0.and.ifwin.eq.0) then
c
c        old procedure for rejecting lines
c
         GX=GFP-EPP/TSTD
         AB0=0.
         if(gx.gt.-30)
     *     AB0=EXP(GFP-EPP/TSTD)*RRR(IDSTD,ION,IAT)/DOPSTD/AVAB
c         write(6,651) iat,ion,alam*10.,gfp,epp,gx,ab0
c 651     format(2i4,f11.3,1p4e11.2)  
         IF(AB0.LT.UN) GO TO 10
C
       else
c
c        new procedure for rejecting lines
c
         DOPSTD=1.E7/ALAM*DSTD
         DOPLAM=ALAM*ALAM/CNM*DOPSTD
c        ijc=2
c        do ijcn=2,nfreqc
         do ijcn=ijc,nfreqc
            if(fr0.ge.freqc(ijcn)) go to 12
         end do
   12    continue
         ijc=ijcn
         if(ijc.gt.nfreqc) ijc=nfreqc
         tkm=1.65e8/amas(iat)
         DP0=3.33564E-11*FR0
         do id=1,nd,ndstep
            td=temp(id)
            gx=gfp-epp/td
            ab0=0.
            if(gx.gt.-30) then
               dops=dp0*sqrt(tkm/td+vturb(id))
               AB0=EXP(gx)*RRR(ID,ION,IAT)/(DOPS*abstdw(ijc,id)*relop)
            end if
            if(ab0.ge.un) go to 15
         end do
         GO TO 10
      end if
C
C     truncate line list if there are more lines than maximum allowable
C     (given by MLIN0 - see include file LINDAT.FOR)
C
   15 continue
      IL=IL+1
      IF(IL.GT.MLIN0) THEN
         WRITE(6,601) ALAM
         IL=MLIN0
         ALAST=CNM/FREQ0(IL)-CUTOFF
         FRLAST=CNM/ALAST
         NXTSET=1
         GO TO 100
      END IF
C
C     =============================================
C     line is selected, set up necessary parameters
C     =============================================
C
C     store parameters for selected lines
C
      FREQ0(IL)=FR0
      EXCL0(IL)=real(EPP)
      EXCU0(IL)=real(EXCU*C3)
      GF0(IL)=real(GFP)
      INDAT(IL)=100*IAT+ION
C
C     indices for corresponding excitation temperatures of the lower
C     and upper levels
C     (for winds)
C
      if(ifwin.gt.0) then
      IJCONT(IL)=IJC
      if(excl.ge.enhe2) then
         ipotl(il)=3
       else if(excl.ge.enhe1) then
         ipotl(il)=2
       else
         ipotl(il)=1
      end if
      end if
C
C     ****** line broadening parameters *****
C
C     1) natural broadening
C
      IF(AGAM.GT.0.) THEN
         GAMR0(IL)=real(EXP(C1*AGAM))
       ELSE
         GAMR0(IL)=real(AGR0*FR0*FR0)
      END IF
C
C     if Stark or Van der Waals broadenig assumed classical,
C     evaluate the effective quantum number
C
      IF(GS.EQ.0..OR.GW.EQ.0) THEN
         Z=FLOAT(ION)
         XNEFF2=Z**2*(XEH/(ENEV(IAT,ION)-EXCU/XET))
         IF(XNEFF2.LE.0..OR.XNEFF2.GT.XNF) XNEFF2=XNF
      END IF
C
C     2) Stark broadening
C
      IF(GS.NE.0.) THEN
         GS0(IL)=real(EXP(C1*GS))
       ELSE
         GS0(IL)=real(TENM8*XNEFF2*XNEFF2*SQRT(XNEFF2))
      END IF
C
C     3) Van der Waals broadening
C
      IF(GW.NE.0.) THEN
         GW0(IL)=real(EXP(C1*GW))
       ELSE
         IF(IAT.LT.21) THEN
            R2=R02*(XNEFF2/Z)**2
          ELSE IF(IAT.LT.45) then
            R2=(R12-FLOAT(IAT))/Z
          ELSE
            R2=0.5
         END IF
         GW0(IL)=real(VW0*R2**OP4)
      END IF
c
C     evaluation of EXTIN0 - the distance (in delta frequency) where
C     the line is supposed to contribute to the total opacity
C
      call profil(il,iat,idstd,agam)
      IF(IAT.LE.2) THEN
         EXT=SQRT(10.*AB0)
       ELSE IF(IAT.LE.14) THEN
         EX0=AB0*ASTD*10.
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
       ELSE
         EX0=AB0*ASTD
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
      END IF
      EXTIN0=EXT*DOPSTD
      EXTIN(IL)=real(EXTIN0)
C
C     4) parameters for a special profile evaluation:
C
C     a) special He I and He II line broadening parameters
C
      ISPRFF=0
      IF(IAT.LE.2) ISPRFF=ISPEC(IAT,ION,ALAM)
      IF(IAT.EQ.2) CALL HESET(IL,ALAM,EXCL,EXCU,ION,IPRF,ILWN,IUN)
      ISPRF(IL)=ISPRFF
      IPRF0(IL)=IPRF
C
C     b) parameters for Griem values of Stark broadening
C
      IF(IPRF.LT.0) THEN
         IGRIE0=IGRIE0+1
         IGRIEM(IL)=IGRIE0
         IF(IGRIE0.GT.MGRIEM) THEN
            WRITE(6,603) ALAM
            GO TO 20
         END IF
         WGR0(1,IGRIE0)=real(WGR1)
         WGR0(2,IGRIE0)=real(WGR2)
         WGR0(3,IGRIE0)=real(WGR3)
         WGR0(4,IGRIE0)=real(WGR4)
      END IF
   20 CONTINUE
C
C     implied NLTE option
C
      if(inlte.eq.-2.or.inlte.eq.12) then
          if(iat.le.20.and.excl.le.1000.) qu=-abs(qu)
       else if(inlte.eq.-3) then
          if(excl.le.1000.) qu=-abs(qu)
       else if(inlte.eq.-4) then
          qu=-abs(qu)
      end if
C
C     NLTE lines initialization
C
      INDNLT(IL)=0
      IF(QU.LT.0..OR.QL.LT.0.) THEN
         ILWN=-1
         QU=ABS(QU)
         QL=ABS(QL)
      END IF
      IF(ILWN.LT.0.AND.INLTE.NE.0) THEN
         INNLT0=INNLT0+1
         INDNLT(IL)=INNLT0
         IF(INNLT0.GT.MNLT) THEN
            WRITE(6,604) ALAM
            GO TO 100
         END IF
         GI=2.*QL+UN
         GJ=2.*QU+UN
         CALL NLTE(IL,ILWN,IUN,GI,GJ)
         ILOWN(IL)=ILWN
         IUPN(IL)=IUN
      END IF
      IF(ILWN.GT.0.AND.INLTE.NE.0) THEN
         INNLT0=INNLT0+1
         INDNLT(IL)=INNLT0
         IF(INNLT0.GT.MNLT) THEN
            WRITE(6,604) ALAM
            GO TO 100
         END IF
         GI=2.*QL+UN
         GJ=2.*QU+UN
         CALL NLTE(IL,ILWN,IUN,GI,GJ)
         ILOWN(IL)=ILWN
         IUPN(IL)=IUN
      END IF
      IF(ILWN.EQ.0.AND.INLTE.GE.1) THEN
         ILMATCH=-1
         CALL NLTSET(1,IL,IAT,ION,ALAM,EXCL,EXCU,QL,QU,
     *         ISQL,ILQL,IPQL,ISQU,ILQU,IPQU,IEVEN,INNLT0,ILMATCH)
C
C        Success accounting for nlte lines matched with quantum numbers and
C        energy limits
C
C        nlte lines searched  matching energies and quantum numbers
         IF(ILMATCH.GE.0) THEN
                ILSEARCH=ILSEARCH+1
C        nlte lines not found matching
                IF (ILMATCH.EQ.0) THEN
                        ILFAIL=ILFAIL+1
C        nlte lines with multiple matches
                ELSE IF (ILMATCH.EQ.2) THEN
                        ILMULT=ILMULT+1
C           nlte lines uniquely matched
                ELSE IF (ILMATCH.EQ.1) THEN
                        ILFOUND=ILFOUND+1
                ENDIF
         ENDIF

         IF(INDNLT(IL).GT.0) THEN
            IF(INDNLT(IL).GT.MNLT) THEN
               WRITE(6,604) ALAM
               GO TO 100
            END IF
            GI=2.*QL+UN
            GJ=2.*QU+UN
            ILWN=ILOWN(IL)
            IUN=IUPN(IL)
            IF(ILWN.EQ.IUN.AND.GI.EQ.GJ) THEN
               INDNLT(IL)=0
               ILOWN(IL)=0
               IUPN(IL)=0
             ELSE
               CALL NLTE(IL,ILWN,IUN,GI,GJ)
            END IF
         END IF
      END IF
      GO TO 10
C
  100 NLIN0=IL
      NNLT=INNLT0
      NGRIEM=IGRIE0
      ALM1=CNM/FREQ0(1)
      IF(ALAM0.LT.ALM1.AND.IMODE.NE.1) THEN
         ALAM0=ALM1-4.*DOPLAM
         IF(ALAM0.LT.ALAM00) ALAM0=ALAM00
      END IF
      ALM2=CNM/FREQ0(NLIN0)
      IF(NLIN0.GT.1) ALM2=CNM/FREQ0(NLIN0-1)
      IF(ALAST.GT.ALM2.AND.IMODE.NE.1) THEN
         ALAST=ALM2-4.*DOPLAM
         IF(ALAST.GT.ALAST0) ALAST=ALAST0
         FRLAST=CNM/ALAST
      END IF
      IBLANK=0
C 
      WRITE(11,*)'INILIN: NLTE matches using Energies and SLP limits --'
      WRITE(11,*)ILSEARCH,' lines searched'
      WRITE(11,*)ILFAIL,' lines unmatched -- set to LTE'
      WRITE(11,*)ILMULT,' lines with multiple matches'
      WRITE(11,*)ILFOUND,' lines uniquely matched'
      WRITE(11,*)'----------------------------------------------------'
C
      WRITE(*,*)'----------------------------------------------------'
      WRITE(6,611) NLIN0,NNLT
  611 FORMAT(/' LINES - TOTAL        :',I10
     *       /' LINES - NLTE         :',I10/)
  601 FORMAT(' **** MORE LINES THAN MLIN0, LINE LIST TRUNCATED '/
     *'       AT LAMBDA',F15.4,'  NM'/)
  603 FORMAT(' **** MORE LINES WITH GRIEM PROFILES THAN MGRIEM'/
     *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
  604 FORMAT(' **** MORE LINES IN NLTE OPTION THAN MNLT'/
     *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
      RETURN
      END
C
C ********************************************************************
C
C
      SUBROUTINE INILIN_grid
C     ======================
C
C     read in the input line list,
C     selection of lines that may contribute,
C     set up auxiliary fields containing line parameters,
C
C     Input of line data - unit 19:
C
C     For each line, one (or two) records, containing:
C
C    ALAM    - wavelength (in nm)
C    ANUM    - code of the element and ion (as in Kurucz-Peytremann)
C              (eg. 2.00 = HeI; 26.00 = FeI; 26.01 = FeII; 6.03 = C IV)
C    GF      - log gf
C    EXCL    - excitation potential of the lower level (in cm*-1)
C    QL      - the J quantum number of the lower level
C    EXCU    - excitation potential of the upper level (in cm*-1)
C    QU      - the J quantum number of the upper level
C    AGAM    = 0. - radiation damping taken classical
C            > 0. - the value of Gamma(rad)
C
C     There are now two possibilities, called NEW and OLD, of the next
C     parameters:
C     a) NEW, next parameters are:
C    GS      = 0. - Stark broadening taken classical
C            > 0. - value of log gamma(Stark)
C    GW      = 0. - Van der Waals broadening taken classical
C            > 0. - value of log gamma(VdW)
C    INEXT   = 0  - no other record necessary for a given line
C            > 0  - next record is read, which contains:
C    WGR1,WGR2,WGR3,WGR4 - Stark broadening values from Griem (in Angst)
C                   for T=5000,10000,20000,40000 K, respectively;
C                   and n(el)=1e16 for neutrals, =1e17 for ions.
C    ILWN    = 0  - line taken in LTE (default)
C            > 0  - line taken in NLTE, ILWN is then index of the
C                   lower level
C            =-1  - line taken in approx. NLTE, with Doppler K2 function
C            =-2  - line taken in approx. NLTE, with Lorentz K2 function
C    IUN     = 0  - population of the upper level in LTE (default)
C            > 0  - index of the lower level
C    IPRF    = 0  - Stark broadening determined by GS
C            < 0  - Stark broadening determined by WGR1 - WGR4
C            > 0  - index for a special evaluation of the Stark
C                   broadening (in the present version inly for He I -
C                   see procedure GAMHE)
C      b) OLD, next parameters are
C     IPRF,ILWN,IUN - the same meaning as above
C     next record with WGR1-WGR4 - again the same meaning as above
C     (this record is automatically read if IPRF<0
C
C     The only differences between NEW and OLD is the occurence of
C     GS and GW in NEW, and slightly different format of reading.
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/igrddd/igrdd,irelin
      common/plaopa/plalin,plcint,chcint
      common/conabs/absoc(mfreqc),emisc(mfreqc),scatc(mfreqc),
     *          plac(mfreqc)
C
      PARAMETER (C1     = 2.3025851,
     *           C2     = 4.2014672,
     *           C3     = 1.4387886,
     *           CNM    = 2.997925D17,
     *           ANUMIN = 1.9,
     *           ANUMAX = 99.31,
     *           AHE2   = 2.01,
     *           EXT0   = 3.17,
     *           UN     = 1.0,
     *           TEN    = 10.,
     *           HUND   = 1.D2,
     *           TENM4  = 1.D-4,
     *           TENM8  = 1.D-8,
     *           OP4    = 0.4,
     *           AGR0=2.4734E-22,
     *           XEH=13.595, XET=8067.6, XNF=25.,
     *           R02=2.5, R12=45., VW0=4.5E-9,
     *           bnc=1.4743e-2,hkc=4.79928e-11)
      PARAMETER (ENHE1=198310.76, ENHE2=438908.85)

      DATA INLSET /0/
C
      if(irelin.eq.0) return
c
      relop0=relop
      relop=1.e-3*relop
      if(relop.gt.1.e-4) relop=1.e-4
      if(relop.lt.1.e-5) relop=1.e-5
      plalin=0.
      ijcon=2
      IL=0
      INNLT0=0
      IGRIE0=0
      IF(NXTSET.EQ.1) THEN
          ALAM0=ALM00
          ALAST=ALST00
          FRLAST=CNM/ALAST
          NXTSET=0
          REWIND 19
      END IF
      ALAM00=ALAM0
      ALAST=CNM/FRLAST
      ALAST0=ALAST
      DOPSTD=1.E7/ALAM0*DSTD
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD
      AVAB=ABSTD(IDSTD)*RELOP
      id=idstd
      dstdid=sqrt(1.4e7*temp(idstd))
      ASTD=1.0
c     IF(GRAV.GT.6.) ASTD=0.1
      CUTOFF=CUTOF0
      ALAST=CNM/FRLAST
      absta=absoc(1)
      write(6,630) alam0,alast,abstd(idstd),absta
c     write(6,631) temp(1),elec(1),dens(1),absoc(1),scatc(1)
  630 format(/' read line list with alam0, alast',2f10.3,1p3e11.3/)
c
      rstd=1.e4
      if(relop.gt.0.) rstd=1./relop
      afac=10.
      if(iat.gt.15.and.iat.ne.26) afac=1.
      afac=afac*rstd*astd
C
      afac=afac*rstd*astd
      afilin=alast
C
C     first part of reading line list - read only lambda, and
C     skip all lines with wavelength below ALAM0-CUTOFF
C
      ALAM=0.
    7 continue
      if(ibin(0).eq.0) then
         read(19,510) alam
       else
         read(19) alam
      end if
  510 FORMAT(F10.4)
      IF(ALAM.LT.ALAM0-CUTOFF) GO TO 7
      BACKSPACE(19)
      GO TO 10
c
    8 continue
   10 ILWN=0
      IUN=0
      IPRF=0
      GS=0.
      GW=0.
      IF(IBIN(0).EQ.0) THEN
         READ(19,*,END=100,err=8) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
     *                        GS,GW
       else
         read(19,end=100) ALAM,ANUM,GF,EXCL,QL,EXCU,QU,AGAM,
     *                        GS,GW
      end if
c 
c     change wavelength to vacuum for lambda > 2000
c
      if(alam.gt.200..and.vaclim.gt.2000.) then
         wl0=alam*10.
         ALM=1.E8/(WL0*WL0)
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0*(XN1*1.D-6+UN)
             alam=wl0*0.1
      END IF
C
C     first selection : for a given interval a atomic number
C
      IF(ALAM.GT.ALAST+CUTOFF) GO TO 100
C
C     second selection : for line strengths
C
      FR0=CNM/ALAM
      if(inlist.ge.0) then
      IAT=ifix(real(ANUM,4))
      FRA=(ANUM-FLOAT(IAT)+TENM4)*HUND
      ION=INT(FRA)+1
      IF(ION.GT.IONIZ(IAT)) GO TO 10
      IEVEN=1
      EXCL=ABS(EXCL)
      EXCU=ABS(EXCU)
      IF(EXCL.GT.EXCU) THEN
         FRA=EXCL
         EXCL=EXCU
         EXCU=FRA
         FRA=QL
         QL=QU
         QU=FRA
         IEVEN=0
      END IF
      GFP=C1*GF-C2
      EPP=C3*EXCL
      else
      IF(ION.GT.IONIZ(IAT)) GO TO 10
      end if
C
      if(fr0.lt.freqc(ijcon)) then
         ijcon=ijcon+1
         absta=0.5*(absoc(ijcon)+scatc(ijcon)+
     *             absoc(ijcon-1)+scatc(ijcon-1))
      end if
      abstd(id)=absta
c
      dop=1.e7/alam*dstdid
      abct=exp(gfp-epp/temp(id))*rrr(id,ion,iat)
      abid=abct/dop/absta
      ext=sqrt(abid*afac)*dop
c
c     line part of the Planck mean opacity
c
c      if(alam.ge.alam0.and.alam.le.alast) then
c      if(abid.ge.relop) then
c        xx=exp(-hkc*fr0/temp(id))
c        pln=bnc*(fr0*1.e-15)**3*xx/(un-xx)
c        abct=abct*(un-xx)
c        plalin=plalin+pln*abct
c        write(16,643) iat,ion,alam*10.,abct,dop,absta,abid
c 643    format(2i4,0pf12.3,1p6e12.4)
c     end if
c
      ALAX0=12.
c
c      alax0=0
c
      if(imode.eq.-6) go to 10
      if(alam.lt.afilin) then
         if(abid.ge.relop) then
            afilin=alam
          else
            if(abid.lt.relop*1.e-6) go to 10
         end if
       else if(alam.lt.9500.) then
         if(abid.lt.relop) go to 10
       else if(alam.lt.9950.) then
         if(abid.lt.relop*1.e-9) go to 10
       else
         if(abid.lt.relop*1.e-19) go to 10
      end if
c
c     if(abid.lt.relop.and.alam.gt.alax0) go to 10
c     if(abid.lt.1.e-10*relop.and.alam.lt.alax0) go to 10
      IF(ANUM.LT.ANUMIN.OR.ANUM.GT.ANUMAX) GO TO 10
      IF(ANUM.GT.ANUMAX) GO TO 10
      IF(ABS(ANUM-AHE2).LT.TENM4.AND.IFHE2.GT.0) GO TO 10
c
      extin0=ext
C
C     truncate line list if there are more lines than maximum allowable
C     (given by MLIN0 - see include file LINDAT.FOR)
C
      IL=IL+1
      IF(IL.GT.MLIN0) THEN
         WRITE(6,601) ALAM
         IL=MLIN0
         ALAST=CNM/FREQ0(IL)-CUTOFF
         FRLAST=CNM/ALAST
         NXTSET=1
         GO TO 100
      END IF
C
C     =============================================
C     line is selected, set up necessary parameters
C     =============================================
C
C     evaluation of EXTIN0 - the distance (in delta frequency) where
C     the line is supposed to contribute to the total opacity
C
C     store parameters for selected lines
C
      FREQ0(IL)=FR0
      EXCL0(IL)=real(EPP,4) 
      EXCU0(IL)=real(EXCU*C3,4) 
      GF0(IL)=real(GFP,4) 
      EXTIN(IL)=real(EXTIN0,4) 
      INDAT(IL)=100*IAT+ION
c     write(14,614) il,iat,ion,alam*10.,extin(il),
c    *              2.997925e18/fr0*extin(il)/fr0,
c    *              abct/dop,absta,abid
c 614 format(3i4,0pf12.3,1p5e11.3)
C
C     ****** line broadening parameters *****
C
C     1) natural broadening
C
      IF(AGAM.GT.0.) THEN
         GAMR0(IL)=real(EXP(C1*AGAM),4) 
       ELSE
         GAMR0(IL)=real(AGR0*FR0*FR0,4) 
      END IF
C
C     if Stark or Van der Waals broadening assumed classical,
C     evaluate the effective quantum number
C
      IF(GS.EQ.0..OR.GW.EQ.0) THEN
         Z=FLOAT(ION)
         XNEFF2=Z**2*(XEH/(ENEV(IAT,ION)-EXCU/XET))
         IF(XNEFF2.LE.0..OR.XNEFF2.GT.XNF) XNEFF2=XNF
      END IF
C
C     2) Stark broadening
C
      IF(GS.NE.0.) THEN
         GS0(IL)=real(EXP(C1*GS),4)
       ELSE
         GS0(IL)=real(TENM8*XNEFF2*XNEFF2*SQRT(XNEFF2),4)
      END IF
C
C     3) Van der Waals broadening
C
      IF(GW.NE.0.) THEN
         GW0(IL)=real(EXP(C1*GW),4)
       ELSE
         IF(IAT.LT.21) THEN
            R2=R02*(XNEFF2/Z)**2
          ELSE IF(IAT.LT.45) then
            R2=(R12-FLOAT(IAT))/Z
          ELSE
            R2=0.5
         END IF
         GW0(IL)=real(VW0*R2**OP4,4)
      END IF
C
C     4) parameters for a special profile evaluation:
C
C     a) special He I and He II line broadening parameters
C
      ISPRFF=0
      IF(IAT.LE.2) ISPRFF=ISPEC(IAT,ION,ALAM)
      IF(IAT.EQ.2) CALL HESET(IL,ALAM,EXCL,EXCU,ION,IPRF,ILWN,IUN)
      ISPRF(IL)=ISPRFF
      IPRF0(IL)=IPRF
C
C     b) parameters for Griem values of Stark broadening
C
      IF(IPRF.LT.0) THEN
         IGRIE0=IGRIE0+1
         IGRIEM(IL)=IGRIE0
         IF(IGRIE0.GT.MGRIEM) THEN
            WRITE(6,603) ALAM
            GO TO 20
         END IF
         WGR0(1,IGRIE0)=real(WGR1,4)
         WGR0(2,IGRIE0)=real(WGR2,4)
         WGR0(3,IGRIE0)=real(WGR3,4)
         WGR0(4,IGRIE0)=real(WGR4,4)
      END IF
   20 CONTINUE
      GO TO 10
C
  100 NLIN0=IL
      NNLT=INNLT0
      NGRIEM=IGRIE0
      ALM1=CNM/FREQ0(1)
      IF(ALAM0.LT.ALM1.AND.IMODE.NE.1) THEN
         ALAM0=ALM1-4.*DOPLAM
         IF(ALAM0.LT.ALAM00) ALAM0=ALAM00
      END IF
      ALM2=CNM/FREQ0(NLIN0)
      IF(NLIN0.GT.1) ALM2=CNM/FREQ0(NLIN0-1)
      IF(ALAST.GT.ALM2.AND.IMODE.NE.1) THEN
         ALAST=ALM2-4.*DOPLAM
         IF(ALAST.GT.ALAST0) ALAST=ALAST0
         FRLAST=CNM/ALAST
      END IF
      IBLANK=0
      relop=relop0
C
      WRITE(6,611) NLIN0
  611 FORMAT(/' ATOMIC LINES        :',I10/)
c     WRITE(6,611) NLIN0,NNLT,NGRIEM
c 611 FORMAT(/' LINES - TOTAL        :',I10
c    *       /' LINES - NLTE         :',I10
c    *       /' LINES - GRIEM        :',I10/)
  601 FORMAT('0 **** MORE LINES THAN MLIN0, LINE LIST TRUNCATED '/
     *'       AT LAMBDA',F15.4,'  NM'/)
c 602 FORMAT('0 **** MORE LINES WITH SPECIAL PROFILES THAN MPRF'/
c    *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
  603 FORMAT('0 **** MORE LINES WITH GRIEM PROFILES THAN MGRIEM'/
     *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
c 604 FORMAT('0 **** MORE LINES IN NLTE OPTION THAN MNLT'/
c    *'       FOR LINES WITH LAMBDA GREATER THAN',F15.4,'  NM'/)
      RETURN
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE INIBLA
C     =================
C
C     driving procedure for treating a partial line list for the
C     current wavelength region
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
C
      PARAMETER (DP0=3.33564E-11, DP1=1.651E8, 
     *           VW1=0.42, VW2=0.3, TENM4=1.E-4)
c    *           VW1=0.42, VW2=0.45,TENM4=1.E-4)
      PARAMETER (UN=1.) 
C
      IF(NLIN.EQ.0) RETURN
      XX=FREQ(1)
      IF(NFREQ.GE.2) XX=0.5*(FREQ(1)+FREQ(2))
      if(ifwin.gt.0) XX=0.5*(FREQC(1)+FREQC(NFREQC))
      BNU=BN*(XX*1.E-15)**3
      HKF=HK*XX
      if(ifwin.gt.0) XX=un
      DO 20 ID=1,ND
         T=TEMP(ID)
         ANE=ELEC(ID)
         EXH=EXP(HKF/T)
         EXHK(ID)=UN/EXH
         PLAN(ID)=BNU/(EXH-UN)
         STIM(ID)=UN-EXHK(ID)
         if(iath.gt.0) then
            ANP=POPUL(NKH,ID)
            AH=DENS(ID)/WMM(ID)/YTOT(ID)-ANP
          else
            ah=rrr(id,1,1)
         end if
         AHE=RRR(ID,1,2)
         VDWC(ID)=(AH+VW1*AHE+0.85*ANH2(ID))*(T*TENM4)**VW2
         DO 10 IAT=1,MATOM
            IF(AMAS(IAT).GT.0.)
     *      DOPA1(IAT,ID)=UN/(XX*DP0*SQRT(DP1*T/AMAS(IAT)+VTURB(ID)))
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
      SUBROUTINE IDTAB
C     ================
C
C     output of selected line parameters (identification table)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      CHARACTER*4 TYPION(30)
      CHARACTER*4 APB,AP0,AP1,AP2,AP3,AP4,APR
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
      COMMON/REFDEP/IREFD(MFRQ)
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)
C
      PARAMETER (C1=2.3025851, C2=4.2014672, C3=1.4387886)
      DATA TYPION /' I  ',' II ',' III',' IV ',' V  ',
     *             ' VI ',' VII','VIII',' IX ',' X  ',
     *             ' XI ',' XII','XIII',' XIV',' XV ',
     *             ' XVI','XVII',' 18 ',' XIX',' XX ',
     *             ' XXI','XXII',' 23 ','XXIV','XXV ',
     *             'XXVI',' 27 ',' 28 ','XXIX',' XXX'/
      DATA APB,AP0,AP1,AP2,AP3,AP4 /'    ','   .','   *','  **',' ***',
     *                              '****'/
C
      IF(NLIN.EQ.0) GO TO 100
C
      ALM0=2.997925D18/FREQ(1)
      ALM1=2.997925D18/FREQ(2)
      if(ifwin.gt.0) ALM0=2.997925D18/FRQOBS(1)
      if(ifwin.gt.0) ALM1=2.997925D18/FRQOBS(NFREQ)
c     IF(IMODE.GE.0) WRITE(6,601) IBLANK,ALM0,ALM1
      IF(IPRIN.LE.-2) RETURN
      if(iprin.ge.2) then
      IF(IMODE.GE.0.OR.(IMODE.EQ.-1.AND.IBLANK.EQ.1)) WRITE(6,602)
      end if
C
      DO 20 IL0=1,NLIN
         IL=INDLIN(IL0) 
         ALAM=2.997925D18/FREQ0(IL)
         ID=IDSTD
         IJCN=IJCNTR(IL0)
         ID0=0
         IF(IJCN.GE.1.AND.IJCN.LE.NFREQS) ID0=IREFD(IJCN)
         IF(ID0.GT.0.and.id0.lt.nd) ID=ID0
         IAT=INDAT(IL)/100
         ION=MOD(INDAT(IL),100)
         CALL PROFIL(IL,IAT,ID,AGAM)
         ABCNT=EXP(GF0(IL)-EXCL0(IL)/TEMP(ID))*RRR(ID,ION,IAT)*
     *         STIM(ID)
         absta=min(ch(1,idstd),ch(2,idstd))
         if(ifwin.le.0) then
            DOP1=DOPA1(IAT,ID)
c           STR0=ABCNT*DOP1/ABSTD(ID)
            str0=abcnt*dop1/absta
          else 
            DOP1=DOPA1(IAT,ID)/FREQ0(IL)
            STR0=ABCNT*DOP1/ABSTDW(IJCONT(IL),ID)
         end if
         GF=(GF0(IL)+C2)/C1
         EXCL=EXCL0(IL)/C3
         IF(STR0.LE.1.2) THEN
            WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577))
          ELSE 
            WW1=SQRT(LOG(STR0))
         END IF
         IF(STR0.GT.55.) THEN
            WW2=0.5*SQRT(3.14*AGAM*STR0)
            IF(WW2.GT.WW1) WW1=WW2
         END IF
         EQW=ALAM/FREQ0(IL)*1.E3/DOP1*WW1
         STR=EQW*10.
         APR=APB
         IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0
         IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1
         IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2
         IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3
         IF(STR.GE.1.E4) APR=AP4
         if(alam.ge.alm0.and.alam.lt.alm1) then
         ill=ilown(il)
         ilu=iupn(il)
         if(ill.gt.0) ill=ill-nfirst(iel(ill))+1
         if(ilu.gt.0) ilu=ilu-nfirst(iel(ilu))+1
       
c        WRITE(12,603) IL0,IL,ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
c    *                 STR0,EQW,APR,ilown(il),iupn(il),id
         WRITE(12,603) ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
     *                 STR0,EQW,APR,ill,ilu,id
         end if
   20 CONTINUE
C
c 601 FORMAT(/' ',I4,'. SET:',
c    * ' INTERVAL  ',F9.3,' -',F9.3,' ANGSTROMS'/
c    *        ' ------------')
  602 FORMAT(/1H ,13X,
     * 'LAMBDA    ATOM    LOG GF       ELO    LINE/CONT',2X,
     * 'EQ.WIDTH'/)
c 603 FORMAT(1H ,I3,I7,F10.3,2X,2A4,F7.2,F12.3,1PE11.2,0PF8.1,1X,A4,
c    *       4i3)
  603 FORMAT(F11.3,2X,A4,A3,F7.2,F12.3,1PE11.2,0PF8.1,1X,A4,
     *       3i4)
C
  100 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
      SUBROUTINE INIBLH
C     =================
C
C     output information about hydrogen lines
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      CHARACTER*4 TYPION(30)
      CHARACTER*4 APB,AP0,AP1,AP2,AP3,AP4,APR
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
C
      PARAMETER (C1=2.3025851, C2=4.2014672, C3=1.4387886)
      PARAMETER (DP0=3.33564E-11, DP1=1.651E8, 
     *           VW1=0.42, VW2=0.45,TENM4=1.E-4)
      PARAMETER (UN=1.) 

      DATA TYPION /' I  ',' II ',' III',' IV ',' V  ',
     *             ' VI ',' VII','VIII',' IX ',' X  ',
     *             ' XI ',' XII','XIII',' XIV',' XV ',
     *             ' XVI','XVII',' 18 ',' XIX',' XX ',
     *             ' XXI','XXII',' 23 ','XXIV','XXV ',
     *             'XXVI',' 27 ',' 28 ','XXIX',' XXX'/
      DATA APB,AP0,AP1,AP2,AP3,AP4 /'    ','   .','   *','  **',' ***',
     *                              '****'/
C
      ALM0=2.997925D18/FREQ(1)
      ALM1=2.997925D18/FREQ(2)
      XX=FREQ(1)
      IF(NFREQ.GE.2) XX=0.5*(FREQ(1)+FREQ(2))
      BNU=BN*(XX*1.E-15)**3
      HKF=HK*XX
c      IF(IMODE.GE.0) WRITE(6,601) IBLANK,ALM0,ALM1
      IF(IPRIN.LE.-2) RETURN
      if(iprin.ge.2) then
c      IF(IMODE.GE.0.OR.(IMODE.EQ.-1.AND.IBLANK.EQ.1)) WRITE(6,602)
      end if
C
      IF(IHYL.LT.0) go to 200
      IAT=1
      ION=1
      IZZ=1
      ID=IDSTD
         T=TEMP(ID)
         ANE=ELEC(ID)
         EXH=EXP(HKF/T)
         EXHK(ID)=UN/EXH
         PLAN(ID)=BNU/(EXH-UN)
         STIM(ID)=UN-EXHK(ID)
         DOPA1(IAT,ID)=UN/(XX*DP0*SQRT(DP1*T/AMAS(IAT)+VTURB(ID)))
      ISERL=ILOWH
      ISERU=ILOWH
      IF(alm0.GT.17000..AND.alm1.LT.21000.) THEN
         ISERL=3
         ISERU=4
       ELSE IF(alm0.GT.22700.) THEN
         ISERL=4
         ISERU=5
         IF(alm0.GT.32800.) ISERU=6
         IF(alm0.GT.44660.) ISERU=7
      END IF
C
      DO 130 I=ISERL,ISERU
      II=I*I
      XII=UN/II
      M1=M10
      IF(I.LT.ILOWH) M1=ILOWH-1
      M2=M1+1
      IF(M1.LT.I+1) M1=I+1
      IF(grav.lt.3..and.M1.LE.6.AND.I.EQ.2) GO TO 110
      IF(grav.lt.3..and.M1.LE.4.AND.I.EQ.1) GO TO 110
      M1=M1-1
      M2=M20+3
      IF(M1.LT.I+1) M1=I+1
  110 CONTINUE
      if(grav.gt.3.) then
         m2=m2+5
         m1=m1-3
         if(m1.gt.i+6) m1=m1-3 
      end if
      if(grav.gt.6.) then
         m2=m2+2
         m1=m1-1
         if(m1.gt.i+6) m1=m1-1 
      end if
      IF(M1.LT.I+1) M1=I+1
      IF(M2.GT.20) M2=20
      ILINH=0
      DO 120 J=M2,M1,-1
         CALL STARK0(I,J,izz,XKIJ,WL0,FIJ,FIJ0)
         ALAM=WL0
C        if(alam.LT.alm0.OR.alam.GE.alm1) then
         if(alam.ge.alm0.and.alam.lt.alm1) then
         ILINH=ILINH+1
         GH=2.*II
         GF=LOG10(FIJ*GH)
         EXCL=109679.*(1.-XII)
         EXCL0H=EXCL*C3
         GF0H=GF*C1-C2
         ABCNT=EXP(GF0H-EXCL0H/TEMP(ID))*RRR(ID,ION,IAT)*
     *          DOPA1(IAT,ID)*STIM(ID)
         STR0=ABCNT/ABSTD(ID)
         IF(STR0.LE.1.2) THEN
            WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577))
          ELSE 
            WW1=SQRT(LOG(STR0))
         END IF
         IF(STR0.GT.55.) THEN
            agam=0.01
            WW2=0.5*SQRT(3.14*AGAM*STR0)
            IF(WW2.GT.WW1) WW1=WW2
         END IF
         EQW=ALAM*ALAM/3.E18*1.E3/DOPA1(IAT,ID)*WW1
         STR=EQW*10.
         APR=APB
         IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0
         IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1
         IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2
         IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3
         IF(STR.GE.1.E4) APR=AP4
         if(iprin.ge.2) then
         WRITE(6,603) i,ilinh,ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
     *                STR0,EQW,APR,i,j
         end if
         WRITE(14,603) i,ilinh,ALAM,TYPAT(IAT),TYPION(ION),GF,EXCL,
     *                 STR0,EQW,APR,i,j
         end if
  120    continue
  130 continue
  200 continue
C
c 601 FORMAT(/' ',I4,'. SET:',
c    * ' INTERVAL  ',F9.3,' -',F9.3,' ANGSTROMS'/
c    *        ' ------------')
c 602 FORMAT(/1H ,13X,
c    * 'LAMBDA    ATOM    LOG GF       ELO    LINE/CONT',2X,
c    * 'EQ.WIDTH'/)
  603 FORMAT(1H ,I3,I7,F10.3,2X,2A4,F7.2,F12.3,1PE11.2,0PF8.1,1X,A4,
     *2i3)
C
      RETURN
      END
C
C ********************************************************************
C
C

      SUBROUTINE NLTSET(MODE,IL,IAT,ION,ALAM0,EXCL,EXCU,QL,QU,
     *         ISQL,ILQL,IPQL,ISQU,ILQU,IPQU,IEVEN,INNLT0,ILMATCH)
C     ===============================================================
C
C     NLTE option -  automatic assignement of level indices
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      PARAMETER (MNION = MIOEX,
     *           MNLEV = MLEVEL,
     *           ECONST=  5.03411142E15)
      PARAMETER (INLLEV = 13)
      COMMON/NL2PAR/ELIMEV(MNION,MNLEV),ELIMOD(MNION,MNLEV),
     *              ELIML(MNION,MNLEV),
     *              ENREV(MNION,MNLEV),ENROD(MNION,MNLEV),
     *              INDEV(MNION,MNLEV),INDOD(MNION,MNLEV),
C     *              INDLV(MNION,MNLEV),
     *              INDLV(MNION,MNLEV),INDIO(MNION),
     *              NEVEN(MNION),NODD(MNION),NODD0,NLEVS(MNION),
     *              IATN(MNION),IONN(MNION),NNION
      COMMON/PRINTP/TYPLEV
      CHARACTER*10 TYPLEV(MLEVEL),typ
      character*4 typ1
      character*2 typ2
      character*2 typin(60)
      data typin /' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',
     *            '11','12','13','14','15','16','17','18','19','20',
     *            '21','22','23','24','25','26','27','28','29','30',
     *            '31','32','33','34','35','36','37','38','39','40',
     *            '41','42','43','44','45','46','47','48','49','50',
     *            '51','52','53','54','55','56','57','58','59','60'/
C
C     +++++++++++++++++++++++++++
C     MODE = 0  -  initialization
C     +++++++++++++++++++++++++++
C
      IF(MODE.EQ.0) THEN
         NNION=0
         READ(INLLEV,*,END=55,ERR=55) NNION    
         IF(NNION.LE.0) GO TO 55
         DO 50 I=1,NNION
            READ(INLLEV,*) IATN(I),IONN(I)
            READ(INLLEV,*) NEVEN(I)
            IF(NEVEN(I).GT.0) THEN
               DO J=1,NEVEN(I)
                  READ(INLLEV,*) ELIMEV(I,J)
               END DO
            READ(INLLEV,*) NODD(I)
            NODD0=NODD(I)
            IF(NODD(I).GT.0) THEN
               DO J=1,NODD(I)
                  READ(INLLEV,*) ELIMOD(I,J)
               END DO
             ELSE
               NODD(I)=NEVEN(I)
               DO 30 J=1,NODD(I)
                  ELIMOD(I,J)=ELIMEV(I,J)
   30          CONTINUE
            END IF
            INDION=0
            DO 40 IONEX=1,NION
               N0I=NFIRST(IONEX)
               IA=NUMAT(IATM(N0I))
               IF(IA.EQ.IATN(I).AND.IZ(IONEX)-1.EQ.IONN(I)) INDION=IONEX
   40       CONTINUE
            IF(INDION.LE.0) THEN
               call quit(' INCONSISTENCY IN UNIT 13 INPUT - NLTE')
            END IF
            NOFF=NFIRST(INDION)-1
c
            ine=1
            ino=1
            do ii=nfirst(indion),nlast(indion)
               TYP=TYPLEV(II)
               typ1=typ(2:5)
               typ2=typ(8:9)
               iev=0
               if(typ1.eq.'even') iev=1
               do k=1,60
                  if(typin(k).eq.typ2) ind=k
               end do
               if(iev.eq.1) then
                  indev(i,ine)=ii
               write(11,*) 'super-e ',i,ii,ine,elimev(i,ine)
                  ine=ine+1
                else
                  indod(i,ino)=ii
               write(11,*) 'super-o ',i,ii,ino,elimod(i,ino)
                  ino=ino+1
               end if
            end do
           END IF
C
   50    CONTINUE
   55    CONTINUE
C
         INDION=NNION
         DO 90 IONEX=1,NION
            N0I=NFIRST(IONEX)
            IA=NUMAT(IATM(N0I))
            if(isemex(ia).ge.1) go to 90
            IONM1=IZ(IONEX)-1
            IF(IA.EQ.1.OR.IA.EQ.2) GO TO 90
            DO 60 I=1,NNION
               IF(IA.EQ.IATN(I).AND.IONM1.EQ.IONN(I)) GO TO 90
   60       CONTINUE
            IF(NFIRST(IONEX).EQ.NLAST(IONEX)) GO TO 90
            INDION=INDION+1
            EION=ENION(NFIRST(IONEX))
            NLEVS(INDION)=NLAST(IONEX)-NFIRST(IONEX)+1
            INDIO(INDION)=IONEX
            NEVEN(INDION)=0
            IATN(INDION)=IA
            IONN(INDION)=IONM1
            DELE=0.
            DO 70 II=NFIRST(IONEX),NLAST(IONEX)
               I=II-NFIRST(IONEX)+1
               E=(EION-ENION(II))*ECONST
               IF(II.LT.NLAST(IONEX)) THEN
                  E1=(EION-ENION(II+1))*ECONST
                  DELE=0.5*(E1-E)
                  ELIML(INDION,I)=E+DELE
                ELSE
                  IF(INLTE.GE.2) THEN
                     ELIML(INDION,I)=E+DELE
                   ELSE
                     ELIML(INDION,I)=EION*ECONST
                  END IF
               END IF
               INDLV(INDION,I)=II
   70       CONTINUE
   90    CONTINUE
         NNION=INDION

C        Header for the table with the level assignments
C
         if(inlte.gt.0.and.iprin.ge.1)
     *      WRITE(11,*)'NLTSET: IAT ION      LAMBDA     EXCL '//
     *      '      EXCU    ILWN  IUN'

         RETURN
      END IF
C
C
C     ++++++++++++++++++++++++++++++++++++++++++
C     MODE > 0  -  level indices for the line IL
C     ++++++++++++++++++++++++++++++++++++++++++
C
      IF(NNION.LE.0) RETURN
      INION=0
      IONM1=ION-1
      DO 102 I=1,NNION
         IF(IAT.EQ.IATN(I).AND.IONM1.EQ.IONN(I)) INION=I
  102 CONTINUE
      if(isemex(iat).ge.1) RETURN
      IF(INION.LE.0) RETURN
      IF(NEVEN(INION).EQ.0) IEVEN=2
      IF(NEVEN(INION).LT.0) GOTO 400
C
      IF(IEVEN.EQ.1) THEN
         IND=0
         DO 110 J=1,NEVEN(INION)
            IF(EXCL.LE.ELIMEV(INION,J)) THEN
               IND=J
               GO TO 120
            END IF
  110    CONTINUE
         ILWN=0
         GO TO 145
  120    CONTINUE
         ILWN=INDEV(INION,IND)
C
         IND=0
         DO 130 J=1,NODD(INION)
            IF(EXCU.LE.ELIMOD(INION,J)) THEN
               IND=J
               GO TO 140
            END IF
  130    CONTINUE
         IUN=0
         GO TO 145
  140    CONTINUE
         IUN=INDOD(INION,IND)

 145     CONTINUE
C
       ELSE IF(IEVEN.EQ.0) THEN
         IND=0
         DO 150 J=1,NODD(INION)
            IF(EXCL.LE.ELIMOD(INION,J)) THEN
               IND=J
               GO TO 160
            END IF
  150    CONTINUE
         ILWN=0
         GO TO 200
  160    CONTINUE
         ILWN=INDOD(INION,IND)
C
         IND=0
         DO 170 J=1,NEVEN(INION)
            IF(EXCU.LE.ELIMEV(INION,J)) THEN
               IND=J
               GO TO 180
            END IF
  170    CONTINUE
         IUN=0
         GO TO 200
  180    CONTINUE
         IUN=INDEV(INION,IND)
  200    CONTINUE
c
c        transition between levels without a distinction in parity
c
       ELSE 

        IF (ILIMITS(INDIO(INION)).EQ.0.OR.INLIST.LT.10) THEN
C        level identification: using only energy limits
C
         IND=0
         DO 210 J=1,NLEVS(INION)
            IF(EXCL.LE.ELIML(INION,J)) THEN
               IND=J
               GO TO 220
            END IF
  210    CONTINUE
         ILWN=0
         IUN=0
         GO TO 300
  220    CONTINUE
         ILWN=INDLV(INION,IND)
C
         IND=0
         DO 230 J=1,NLEVS(INION)
            IF(EXCU.LE.ELIML(INION,J)) THEN
               IND=J
               GO TO 240
            END IF
  230    CONTINUE
         IUN=0
         GO TO 300
  240    CONTINUE
         IUN=INDLV(INION,IND)
  300    CONTINUE

        ELSE IF (ILIMITS(INDIO(INION)).EQ.1) THEN
C
C        level identification: using energy limits and quantum numbers
C

         IND=0
         INMATCHL=0
         DO 310 J=1,NLEVS(INION)

            IF(EXCL.GE.ENION1(INDLV(INION,J))   .AND.
     *         EXCL.LE.ENION2(INDLV(INION,J)).   AND.
     *         ((IPQL.GE.PQUANT1(INDLV(INION,J)).AND.
     *         IPQL.LE.PQUANT2(INDLV(INION,J))).OR.
     *         (IPQL.EQ.-1))                    .AND.
     *         ((ISQL.GE.SQUANT1(INDLV(INION,J)).AND.
     *         ISQL.LE.SQUANT2(INDLV(INION,J))).OR.
     *         (ISQL.EQ.-1))                    .AND.
     *         ((ILQL.GE.LQUANT1(INDLV(INION,J)).AND.
     *         ILQL.LE.LQUANT2(INDLV(INION,J))).OR.
     *         (ILQL.EQ.-1))
     *                                             ) THEN

               IND=J
               INMATCHL=INMATCHL+1
C               GO TO 320
            END IF
  310    CONTINUE
         IF (INMATCHL.GT.1)
     *       WRITE(11,'(A55,1X,F12.4)')
     *       ' NLTSET: WARNING-- multiple matches for lower level of ',
     *       ALAM0
         IF (INMATCHL.GT.0) GO TO 320
         ILWN=0
         IUN=0
         GO TO 350
  320    CONTINUE
         ILWN=INDLV(INION,IND)
C
C
         IND=0
         INMATCHU=0
         DO 330 J=1,NLEVS(INION)

            IF(EXCU.GE.ENION1(INDLV(INION,J))   .AND.
     *         EXCU.LE.ENION2(INDLV(INION,J)).   AND.
     *         ((IPQU.GE.PQUANT1(INDLV(INION,J)).AND.
     *         IPQU.LE.PQUANT2(INDLV(INION,J))).OR.
     *         (IPQU.EQ.-1))                    .AND.
     *         ((ISQU.GE.SQUANT1(INDLV(INION,J)).AND.
     *         ISQU.LE.SQUANT2(INDLV(INION,J))).OR.
     *         (ISQU.EQ.-1))                    .AND.
     *         ((ILQU.GE.LQUANT1(INDLV(INION,J)).AND.
     *         ILQU.LE.LQUANT2(INDLV(INION,J))).OR.
     *         (ILQU.EQ.-1))
     *                                             ) THEN

               IND=J
               INMATCHU=INMATCHU+1
C               GO TO 340
            END IF
  330    CONTINUE
         IF (INMATCHU.GT.1)
     *       WRITE(11,'(A55,1X,F12.4)')
     *       ' NLTSET: WARNING-- multiple matches for upper level of ',
     *       ALAM0
         IF (INMATCHU.GT.0) GO TO 340
         IUN=0
         GO TO 350
  340    CONTINUE
         IUN=INDLV(INION,IND)
  350    CONTINUE

         IF (INMATCHL.EQ.0.or.INMATCHU.EQ.0) THEN
                ILMATCH=0
           ELSE IF (INMATCHL.GT.1.or.INMATCHU.GT.1) THEN
                ILMATCH=2
           ELSE
                ILMATCH=1
         ENDIF

        ELSE

         write(11,*)('ILIMITS is neither 0 or 1')

        END IF

        if(inlte.gt.0.and.iprin.ge.1)
     *     WRITE(11,'(10x,2(i2,1x),3x,3(F10.3,1x),2(i4,1x))')IAT,ION,
     *              ALAM0,EXCL,EXCU,ILWN,IUN

      END IF
C
  400 IF(NEVEN(INION).LT.0) THEN
        NEV1=-NEVEN(INION)
        IF(IEVEN.EQ.1) THEN
          ILWN=0
          J=1
          DO WHILE (ILWN.EQ.0 .AND. J.LE.NEV1)
            IF(QL.EQ.ELIMEV(INION,J)) THEN
              DE=ENREV(INION,J)
              IF(EXCL.NE.0.) DE=(EXCL-DE)/EXCL
              IF(ABS(DE).LT.1.D-5) ILWN=INDEV(INION,J)
            END IF
            J=J+1
          END DO
          IUN=0
          J=1
          DO WHILE (IUN.EQ.0 .AND. J.LE.NODD(INION))
            IF(QU.EQ.ELIMOD(INION,J)) THEN
              DE=(EXCU-ENROD(INION,J))/EXCU
              IF(ABS(DE).LT.1.D-5) IUN=INDOD(INION,J)
            END IF
            J=J+1
          END DO
        ELSE IF(IEVEN.EQ.0) THEN
          ILWN=0
          J=1
          DO WHILE (ILWN.EQ.0 .AND. J.LE.NODD(INION))
            IF(QL.EQ.ELIMOD(INION,J)) THEN
              DE=ENROD(INION,J)
              IF(EXCL.NE.0.) DE=(EXCL-DE)/EXCL
              IF(ABS(DE).LT.1.D-5) ILWN=INDOD(INION,J)
            END IF
            J=J+1
          END DO
          IUN=0
          J=1
          DO WHILE (IUN.EQ.0 .AND. J.LE.NEV1)
            IF(QU.EQ.ELIMEV(INION,J)) THEN
              DE=(EXCU-ENREV(INION,J))/EXCU
              IF(ABS(DE).LT.1.D-5) IUN=INDEV(INION,J)
            END IF
            J=J+1
          END DO
        END IF
      END IF
C
      IF(INLTE.EQ.5) THEN
         INNLT0=INNLT0+1
         INDNLT(IL)=INNLT0
        ELSE IF(INLTE.EQ.4) THEN
         IF(ILWN.GT.0.AND.IUN.GT.0) THEN
            INDNLT(IL)=-1
         END IF
        ELSE IF(INLTE.EQ.3) THEN
         IF(ILWN.GT.0) THEN
            INDNLT(IL)=-1
         END IF
        ELSE
         INDNLT(IL)=-1
      END IF
      BNUL(IL)=real(BN*(FREQ0(IL)*1.E-15)**3)
      ILOWN(IL)=ILWN
      IUPN(IL)=IUN
      RETURN
      END
C
C ********************************************************************
C ********************************************************************
C
C

      SUBROUTINE PHTION(ID,ABSO,EMIS,FRE,NFRE)
C     ========================================
C
C     Opacity due to detailed photoionization (read from tables by
C     routine READPH)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/PHOTCS/PHOT(MFRQ,MPHOT),WPHT0,WPHT1,APHT(MPHOT),
     *              EPHT(MPHOT),GPHT(MPHOT),JPHT(MPHOT),
     *              NPHT
       DIMENSION ABSO(MFRQ),EMIS(MFRQ),PLANF(MFRQ),STIMU(MFRQ)
       DIMENSION FRE(MFRQ)
       PARAMETER (C3=1.4387886)
C
      IF(NPHT.LE.0) RETURN
      T=TEMP(ID)
      DO 10 IJ=1,NFRE
         XX=FRE(IJ)
         X15=XX*1.E-15
         BNU=BN*X15*X15*X15
         HKF=HK*XX
         EXH=EXP(HKF/T)
         PLANF(IJ)=BNU/(EXH-1.)
         STIMU(IJ)=1.-1./EXH
   10 CONTINUE
      DO 30 I=1,NPHT
         IF(JPHT(I).LE.0) THEN
            IAT=int(APHT(I))
            X=(APHT(I)-FLOAT(IAT)+1.E-4)*1.E2
            ION=INT(X)+1
            POP=RRR(ID,ION,IAT)*GPHT(I)*EXP(-EPHT(I)*C3/T)
          ELSE
            JJ=JPHT(I)
            POP=POPUL(JJ,ID)
         END IF
         DO 20 IJ=1,NFRE
            AB=PHOT(IJ,I)*POP*STIMU(IJ)
            ABSO(IJ)=ABSO(IJ)+AB
            EMIS(IJ)=EMIS(IJ)+AB*PLANF(IJ)
   20    CONTINUE
   30 CONTINUE
      RETURN
      END

C
C ********************************************************************
C ********************************************************************
C
      SUBROUTINE NLTE(IL,ILW,IUN,GI,GJ)
C     ===========================================
C
C     Control procedure for the NLTE option
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/NLTPOP/PNLT(MATOM,MION,MDEPTH)
      PARAMETER (UN     = 1., 
     *           C3     = 1.4387886,
     *           XET    = 8067.6,
     *           XET3   = XET*C3)
C
C     CALCULATION OF THE 
C     CENTRAL OPACITY (ABCENT) AND THE LINE SOURCE FUNCTION  (SLIN)
C
      if(gi.le.0..or.gj.le.0.) return
      ILNLT=INDNLT(IL)
      IF(ILNLT.LE.0) RETURN
      IAT=INDAT(IL)/100
      ION=MOD(INDAT(IL),100)
      EGF=EXP(GF0(IL))
      BNU=BN*(FREQ0(IL)*1.E-15)**3
      DP0=3.33564E-11*FREQ0(IL)
      DP1=1.651E8/AMAS(IAT)
      IF(ILW.LE.0) GO TO 100
C
C     line is a transition between explicit levels of the
C     input model
C
      NKI=NNEXT(IEL(ILW))
      DO 60 ID=1,ND
         T=TEMP(ID)
         COR=1.
         PP=PNLT(IAT,ION,ID)
         IF(ILW.GT.0) THEN
            PI=POPUL(ILW,ID)/G(ILW)
          ELSE
            PI=PP*EXP((ENEV(IAT,ION)*XET3-EXCL0(IL))/T)
         END IF
         IF(IUN.GT.0) THEN
            PJ=POPUL(IUN,ID)/G(IUN)
            cor=(excu0(il)-excl0(il)+
     *         (enion(iun)-enion(ilw))/1.38054e-16)/t
            cor=exp(cor)
          ELSE
            PJ=PP*EXP((ENEV(IAT,ION)*XET3-EXCU0(IL))/T)
         END IF
         if(pj.gt.0.) then
            X=PI/PJ*cor
          else
            x=un
         end if
         IF(X.EQ.UN) X=EXP(4.79928E-11*FREQ0(IL)/T)
         DOP=DP0*SQRT(DP1*T+VTURB(ID))
         SLIN(ILNLT,ID)=BNU/(X-UN)
         if(pi.gt.0.) ABCENT(ILNLT,ID)=PI*(UN-UN/X)*EGF/DOP
   60 CONTINUE
      RETURN
C
C     Approximate NLTE for resonance lines - second order escape
C     probablity theory form of the source function
C
C     Optical depth scale
C
  100 CONTINUE
      ALMIL=2.997925E17/FREQ0(IL)
      HKF=HK*FREQ0(IL)
      DO 110 ID=1,ND
         T=TEMP(ID)
         DOP=DP0*SQRT(DP1*T+VTURB(ID))
         X=EXP(HKF/T)
         ABCENT(ILNLT,ID)=EGF*EXP(-EXCL0(IL)/T)*RRR(ID,ION,IAT)/
     *                    DOP*(1.-1./X)
         AB=ABSTD(ID)+ABCENT(ILNLT,ID)*1.77245
         if(ifwin.gt.0) 
     *   AB=ABSTDW(IJCONT(IL),ID)+ABCENT(ILNLT,ID)*1.77245
         IF(ID.EQ.1) THEN
            ABM=AB/DENS(1)
            TAU=0.5*DM(1)*ABM
          ELSE
            AB0=AB/DENS(ID)
            TAU=TAU+0.5*(DM(ID)-DM(ID-1))*(AB0+ABM)
            ABM=AB0
         END IF
C
C        approximate epsilon after Kastner
C
         E=EPS(T,ELEC(ID),ALMIL,ION,IUN)
         XK2=XK2DOP(TAU)
         SLIN(ILNLT,ID)=SQRT(E/(E+(1.-E)*XK2))*BNU/(X-1.)
  110 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
 
      SUBROUTINE LINOP(ID,ABLIN,EMLIN,AVAB)
C     =====================================
C
C     TOTAL LINE OPACITY (ABLIN)  AND EMISSIVITY (EMLIN)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      PARAMETER (UN     = 1., 
     *           EXT0   = 3.17,
     *           TEN    = 10.,
     *           C3     = 1.4387886,
     *           XET    = 8067.6,
     *           XET3   = XET*C3)
      DIMENSION ABLIN(MFREQ),EMLIN(MFREQ),ABLINN(MFREQ)
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
      COMMON/NLTPOP/PNLT(MATOM,MION,MDEPTH)
      common/lasers/lasdel
C
      DO 10 IJ=1,NFREQ
         ABLIN(IJ)=0.
         ABLINN(IJ)=0.
         EMLIN(IJ)=0.
   10 CONTINUE
C
      IF(NLIN.EQ.0) RETURN
C
C     overall loop over contributing lines
C
      TEM1=UN/TEMP(ID)
      DO 100 I=1,NLIN
         IL=INDLIN(I)
         INNLT=INDNLT(IL)
         IAT=INDAT(IL)/100
         ION=MOD(INDAT(IL),100)
         LPR=.TRUE.
         ISP=ISPRF(IL)
         IF(ISP.GT.1.AND.ISP.LE.5) LPR=.FALSE.
         IF (ISP.GE.6) GO TO 100
         CALL PROFIL(IL,IAT,ID,AGAM)
         DOP1=DOPA1(IAT,ID)
         FR0=FREQ0(IL)
         IF(INNLT.EQ.0) THEN
            AB0=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*
     *          DOP1*STIM(ID)
          ELSE IF(INNLT.GT.0) THEN
            AB0=ABCENT(INNLT,ID)
            SL0=SLIN(INNLT,ID)
          ELSE
            ILW=ILOWN(IL)
            IUN=IUPN(IL)
            COR=1.
            PP=PNLT(IAT,ION,ID)
            IF(ILW.GT.0) THEN
               PI=POPUL(ILW,ID)/G(ILW)
             ELSE
               PI=PP*EXP((ENEV(IAT,ION)*XET3-EXCL0(IL))*TEM1)
            END IF
            IF(IUN.GT.0) THEN
               PJ=POPUL(IUN,ID)/G(IUN)
               cor=(excu0(il)-excl0(il)+
     *             (enion(iun)-enion(ilw))/1.38054e-16)*tem1
               cor=exp(cor)
             ELSE
               PJ=PP*EXP((ENEV(IAT,ION)*XET3-EXCU0(IL))*TEM1)
            END IF
            if(pj.gt.0.) then
               X=PI/PJ*cor
             else
               x=un
            end if
            IF(X.EQ.UN) X=EXP(4.79928E-11*FREQ0(IL)*TEM1)
            SL0=BNUL(IL)/(X-UN)
            ab0=0.
            if(pi.gt.0.) AB0=PI*(UN-UN/X)*EXP(GF0(IL))*DOP1
         END IF
         if(ab0.le.0.and.lasdel) go to 100
C
C        set up limiting frequencies where the line I is supposed to
C        contribute to the opacity 
C
         EX0=AB0/AVAB*AGAM 
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
         EXT=EXT/DOP1
         XIJEXT=DFRCON*EXT+1.5
c        IJ1=MAX(IJCNTR(I)-IJEXT,3)
c        IJ2=MIN(IJCNTR(I)+IJEXT,NFREQS)
         IJ1=int(MAX(float(IJCNTR(I))-XIJEXT,3.))
         IJ2=int(MIN(float(IJCNTR(I))+XIJEXT,float(NFREQS)))
         IF(IJ1.GE.NFREQ.OR.IJ2.LE.2) GO TO 100
C
         IF(INNLT.EQ.0) THEN
C
C        *********
C        LTE lines
C        *********
C
         IF(LPR) THEN
C
            DO 40 IJ=IJ1,IJ2
               XF=ABS(FREQ(IJ)-FR0)*DOP1
               ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF)
   40       CONTINUE
C
C        special expressions for 4 selected He I lines
C
         ELSE
            DO 60 IJ=3,NFREQ
               FR=FREQ(IJ)
               ABL=AB0*PHE1(ID,FR,ISP-1)
               ABLIN(IJ)=ABLIN(IJ)+ABL
   60       CONTINUE
         END IF
C
C        **********
C        NLTE LINES
C        **********
C
       ELSE
         IF(LPR) THEN
C
            DO 80 IJ=IJ1,IJ2
               XF=ABS(FREQ(IJ)-FR0)*DOP1
               ABL=AB0*VOIGTK(AGAM,XF)
               ABLINN(IJ)=ABLINN(IJ)+ABL
               EMLIN(IJ)=EMLIN(IJ)+ABL*SL0
   80       CONTINUE
C
C        again, special expressions for 4 selected He I lines
C
         ELSE
            DO 90 IJ=3,NFREQ
               FR=FREQ(IJ)
               ABL=AB0*PHE1(ID,FR,ISP-1)
               ABLINN(IJ)=ABLINN(IJ)+ABL
               EMLIN(IJ)=EMLIN(IJ)+ABL*SL0
   90       CONTINUE
         END IF
      END IF
  100 CONTINUE
C
      DO 110 IJ=3,NFREQ
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID)
         ABLIN(IJ)=ABLIN(IJ)+ABLINN(IJ)
  110 CONTINUE
C
C     special routine for selected He II lines
C
      IF(NSP.EQ.0) RETURN
      DO 120 IS=1,NSP
         ISP=ISP0(IS)
         IF(ISP.GE.6.AND.ISP.LE.24) CALL PHE2(ISP,ID,ABLIN,EMLIN)
  120 CONTINUE
C
      RETURN
      END
C
C ********************************************************************
C
C
 
      SUBROUTINE LINOPW(ID,ABLIN,EMLIN)
C     =================================
C
C     TOTAL LINE OPACITY (ABLIN)  AND EMISSIVITY (EMLIN)
C     (a variant for winds)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      INCLUDE 'WINCOM.FOR'
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      PARAMETER (UN     = 1., 
     *           EXT0   = 3.17,
     *           TEN    = 10.,
     *           C3     = 1.4387886,
     *           XET    = 8067.6,
     *           XET3   = XET*C3)
      DIMENSION ABLIN(MFREQ),EMLIN(MFREQ),ABLINN(MFREQ)
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
      COMMON/NLTPOP/PNLT(MATOM,MION,MDEPTH)
      COMMON/IPOTLS/IPOTL(mlin0)
      common/lasers/lasdel
      common/linrej/ilne(mdepth),ilvi(mdepth)
      common/velaux/velmax,iemoff,nltoff,itrad
C
      DO 10 IJ=1,NFREQ
         ABLIN(IJ)=0.
         ABLINN(IJ)=0.
         EMLIN(IJ)=0.
   10 CONTINUE
      wdil(id)=1.
      plw=plan(id)*wdil(id)
c      plw=xjcon(id)
C
      IF(NLIN.EQ.0) RETURN
C
C     overall loop over contributing lines
C
      TEM1=UN/TEMP(ID)
      HKT=HK*TEM1
      xx=freq(nopac)-freq(1)
      DFRCON=NOPAC-1
      DFRCON=-DFRCON/XX
      IFRCON=int(DFRCON)
      DO 100 I=1,NLIN
         IL=INDLIN(I)
         INNLT=INDNLT(IL)
c        
c        rejecting lines for v > velmax
c
         if(ilvi(id).gt.0) then
            if(innlt.eq.0) then
               go to 100
             else
               if(nltoff.ne.0) go to 100
            end if
         end if
c
c
c     frequency indices of the line centers
c
         if (id.eq.1) then

         fr0=freq0(il)
         XJC=3.+DFRCON*(FREQ(1)-FR0)
         IJC=int(XJC)
         IJCNTR(I)=IJC
         if(ijc.le.1.or.ijc.ge.nopac) go to 255
         if(fr0.lt.freq(ijc)) then
            ijc0=ijc
            dfr0=freq(ijc0)-fr0
  252       ijc0=ijc0+1
            dfr=abs(freq(ijc0)-fr0)
            if(dfr.lt.dfr0) then
               ijc=ijc0
               ijc0=ijc0+1
               dfr0=dfr
               go to 252
            end if
          else if(fr0.gt.freq(ijc)) then
            ijc0=ijc
            dfr0=fr0-freq(ijc0)
  254       ijc0=ijc0-1
            dfr=abs(freq(ijc0)-fr0)
            if(dfr.lt.dfr0) then
               ijc=ijc0
               ijc0=ijc0-1
               dfr0=dfr
               go to 254
            end if
         end if
         IJCNTR(I)=IJC
  255  continue
c        write(80,*) i,ijcntr(i),2.997925e18/freq0(il)
         endif
c
         IAT=INDAT(IL)/100
         ION=MOD(INDAT(IL),100)
         FR0=FREQ0(IL)
         LPR=.TRUE.
         ISP=ISPRF(IL)
         IF(ISP.GT.1.AND.ISP.LE.5) LPR=.FALSE.
         IF (ISP.GE.6) GO TO 100
         CALL PROFIL(IL,IAT,ID,AGAM)
         DOP1=DOPA1(IAT,ID)/FR0
         FR0=FREQ0(IL)
         IF(INNLT.EQ.0) THEN
            if(itrad.le.0) then
            AB0=EXP(GF0(IL)-EXCL0(IL)*TEM1)*RRR(ID,ION,IAT)*
     *          DOP1*(1.-exp(-hkt*fr0))
            else
            trl=trad(ipotl(il),id)
            xx=exp(-hkt*fr0)
            AB0=EXP(GF0(IL)-EXCL0(IL)/trl)*RRR(ID,ION,IAT)*
     *          DOP1*(1.-xx)
            if(excl0(il).gt.2000.) ab0=ab0*wdil(id)
            pla=1.4743e-2*(fr0*1.e-15)**3*xx/(1.-xx)
            sl0=pla*wdil(id)
            end if
          ELSE IF(INNLT.GT.0) THEN
            AB0=ABCENT(INNLT,ID)
            SL0=SLIN(INNLT,ID)
          ELSE
            ILW=ILOWN(IL)
            IUN=IUPN(IL)
            COR=1.
            PP=PNLT(IAT,ION,ID)
            IF(ILW.GT.0) THEN
               PI=POPUL(ILW,ID)/G(ILW)
             ELSE
               PI=PP*EXP((ENEV(IAT,ION)*XET3-EXCL0(IL))*TEM1)
            END IF
            IF(IUN.GT.0) THEN
               PJ=POPUL(IUN,ID)/G(IUN)
               cor=(excu0(il)-excl0(il)+
     *             (enion(iun)-enion(ilw))/1.38054e-16)*tem1
               cor=exp(cor)
             ELSE
               PJ=PP*EXP((ENEV(IAT,ION)*XET3-EXCU0(IL))*TEM1)
            END IF
            if(pj.gt.0.) then
               X=PI/PJ*cor
             else
               x=un
            end if
            IF(X.EQ.UN) X=EXP(4.79928E-11*FREQ0(IL)*TEM1)
            SL0=BNUL(IL)/(X-UN)
            ab0=0.
            if(pi.gt.0.) AB0=PI*(UN-UN/X)*EXP(GF0(IL))*DOP1
         END IF
         if(ab0.le.0.and.lasdel) go to 100
C
C        set up limiting frequencies where the line I is supposed to
C        contribute to the opacity 
C
c        if(ifwin.le.0) then
         avabw=abstdw(ijcont(il),id)*relop
         EX0=AB0/AVABw*AGAM 
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
         EXT=EXT/DOP1
         IJEXT=int((DFRCON*EXT)+1.5)
         IJ1=MAX(IJCNTR(I)-IJEXT,1)
         IJ2=MIN(IJCNTR(I)+IJEXT,NFREQ)
         IF(IJ1.GE.NFREQ.OR.IJ2.LE.2) GO TO 100
c        else
c        ij1=3
c        ij2=nfreq
c        end if
C
         IF(INNLT.EQ.0.and.itrad.le.0) THEN
C
C        *********
C        LTE lines
C        *********
C
         IF(LPR) THEN
C
            DO 40 IJ=IJ1,IJ2
               XF=ABS(FREQ(IJ)-FR0)*DOP1
               ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF)
   40       CONTINUE
C
C        special expressions for 4 selected He I lines
C
         ELSE
            DO 60 IJ=1,NFREQ
               FR=FREQ(IJ)
               ABL=AB0*PHE1(ID,FR,ISP-1)
               ABLIN(IJ)=ABLIN(IJ)+ABL
   60       CONTINUE
         END IF
C
C        **********
C        NLTE LINES
C        **********
C
       ELSE
         IF(LPR) THEN
C
            DO 80 IJ=IJ1,IJ2
               XF=ABS(FREQ(IJ)-FR0)*DOP1
               ABL=AB0*VOIGTK(AGAM,XF)
               ABLINN(IJ)=ABLINN(IJ)+ABL
               if(ilne(id).gt.0) go to 80
               EMLIN(IJ)=EMLIN(IJ)+ABL*SL0
   80       CONTINUE
C
C        again, special expressions for 4 selected He I lines
C
         ELSE
            DO 90 IJ=1,NFREQ
               FR=FREQ(IJ)
               ABL=AB0*PHE1(ID,FR,ISP-1)
               ABLINN(IJ)=ABLINN(IJ)+ABL
               if(ilne(id).gt.0) go to 90
               EMLIN(IJ)=EMLIN(IJ)+ABL*SL0
   90       CONTINUE
         END IF
      END IF
  100 CONTINUE
C
      if(vel(id).le.velmax) then
      DO 110 IJ=1,NFREQ
         PLA=BNUE(IJ)/(EXP(HKT*FREQ(IJ))-1.)
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*pla*wdil(id)
         ABLIN(IJ)=ABLIN(IJ)+ABLINN(IJ)
  110 CONTINUE
      end if
C
C     special routine for selected He II lines
C
      IF(NSP.EQ.0) RETURN
      DO 120 IS=1,NSP
         ISP=ISP0(IS)
         IF(ISP.GE.6.AND.ISP.LE.24) CALL PHE2(ISP,ID,ABLIN,EMLIN)
  120 CONTINUE
C
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE PROFIL(IL,IAT,ID,AGAM)
C     =================================
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
      DIMENSION WGR(4)
      PARAMETER (PI4=7.95774715E-2)
C
      IPRF=IPRF0(IL)
      T=TEMP(ID)
      ANE=ELEC(ID)
C
C     radiative broadening (classical)
C
      AGAM=GAMR0(IL)
C
C     Stark broadening - standard (given in the line list or classical)
C
      IF(IPRF.EQ.0) THEN
         AGAM=AGAM+GS0(IL)*ANE
C
C     Stark broadening - special expressions for He I
C
       ELSE IF(IPRF.GT.0) THEN
         ANP=POPUL(NKH,ID)
         CALL GAMHE(IPRF,T,ANE,ANP,ID,GAM)
         AGAM=AGAM+GAM
C
C     Stark broadening - Griem
C
       ELSE 
         DO 10 I=1,4
   10       WGR(I)=WGR0(I,IGRIEM(IL))
         FR=FREQ0(IL)
         ION=MOD(INDAT(IL),100)
         CALL GRIEM(ID,T,ANE,ION,FR,WGR,GAM)
         AGAM=AGAM+GAM
      END IF
C
C     Van Der Waals broadening
C
      AGAM=AGAM+GW0(IL)*VDWC(ID)
C
C     final Voigt parameter a
C
      DOP1=DOPA1(IAT,ID)
      if(ifwin.gt.0) DOP1=DOP1/FREQ0(IL)
      AGAM=AGAM*DOP1*PI4
C
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE GRIEM(ID,T,ANE,ION,FR,WGR,GAM)
C     =========================================
C
C     STARK DAMPING PARAMETER (GAM) CALCULATED FROM INPUT VALUES
C     OF STARK WIDTHS FOR  T=5000, 10000, 20000, 40000 K,
C     AND FOR  NE=1.E16 (FOR NEUTRALS)  OR  NE = 1.E17 (FOR IONS)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION WGR(4)
      if(t.le.0.) return
      J=JT(ID)
      GAM=(TI0(ID)*WGR(J)+TI1(ID)*WGR(J-1)+TI2(ID)*WGR(J-2))
     *    *ANE*1.E-10*FR*1.E-10*FR*4.2E-14
      IF(ION.GT.1) GAM=GAM*0.1
      IF(GAM.LT.0.) GAM=0.
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE GAMHE(IND,T,ANE,ANP,ID,GAM)
C     ======================================
C
C     NEUTRAL HELIUM STARK BROADENING PARAMETERS
C     AFTER DIMITRIJEVIC AND SAHAL-BRECHOT, 1984, J.Q.S.R.T. 31, 301
C     OR FREUDENSTEIN AND COOPER, 1978, AP.J. 224, 1079  (FOR C(IND).GT.0)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION W(5,20),V(4,20),C(20)
C
C   ELECTRONS T= 5000   10000   20000   40000     LAMBDA
C
      DATA W /  5.990,  6.650,  6.610,  6.210,    3819.60,
     *          2.950,  3.130,  3.230,  3.300,    3867.50,
     *          0.000,  0.000,  0.000,  0.000,    3871.79,
     *          0.142,  0.166,  0.182,  0.190,    3888.65,
     *          0.000,  0.000,  0.000,  0.000,    3926.53,
     *          1.540,  1.480,  1.400,  1.290,    3964.73,
     *         41.600, 50.500, 57.400, 65.800,    4009.27,
     *          1.320,  1.350,  1.380,  1.460,    4120.80,
     *          7.830,  8.750,  8.690,  8.040,    4143.76,
     *          5.830,  6.370,  6.820,  6.990,    4168.97,
     *          0.000,  0.000,  0.000,  0.000,    4437.55,
     *          1.630,  1.610,  1.490,  1.350,    4471.50,
     *          0.588,  0.620,  0.641,  0.659,    4713.20,
     *          2.600,  2.480,  2.240,  1.960,    4921.93,
     *          0.627,  0.597,  0.568,  0.532,    5015.68,
     *          1.050,  1.090,  1.110,  1.140,    5047.74,
     *          0.277,  0.298,  0.296,  0.293,    5875.70,
     *          0.714,  0.666,  0.602,  0.538,    6678.15,
     *          3.490,  3.630,  3.470,  3.190,    4026.20,
     *          4.970,  5.100,  4.810,  4.310,    4387.93/
C
C   PROTONS   T= 5000   10000   20000   40000
C
      DATA V /  1.520,  4.540,  9.140, 10.200,
     *          0.607,  0.710,  0.802,  0.901,
     *          0.000,  0.000,  0.000,  0.000,
     *          0.0396, 0.0434, 0.0476, 0.0526,
     *          0.000,  0.000,  0.000,  0.000,
     *          0.507,  0.585,  0.665,  0.762,
     *          0.930,  1.710, 13.600, 27.200,
     *          0.288,  0.325,  0.365,  0.410,
     *          1.330,  6.800, 12.900, 14.300,
     *          1.100,  1.370,  1.560,  1.760,
     *          0.000,  0.000,  0.000,  0.000,
     *          1.340,  1.690,  1.820,  1.630,
     *          0.128,  0.143,  0.161,  0.181,
     *          2.040,  2.740,  2.950,  2.740,
     *          0.187,  0.210,  0.237,  0.270,
     *          0.231,  0.260,  0.291,  0.327,
     *          0.0591, 0.0650, 0.0719, 0.0799,
     *          0.231,  0.260,  0.295,  0.339,
     *          2.180,  3.760,  4.790,  4.560,
     *          1.860,  5.320,  7.070,  7.150/
      DATA C /2*0.,1.83E-4,0.,1.13E-4,5*0.,1.6E-4,9*0./
C
      IF(W(1,IND).EQ.0.) GO TO 10
      J=JT(ID)
      GAM=((TI0(ID)*W(J,IND)+TI1(ID)*W(J-1,IND)+TI2(ID)*W(J-2,IND))
     *     *ANE
     *    +(TI0(ID)*V(J,IND)+TI1(ID)*V(J-1,IND)+TI2(ID)*V(J-2,IND))
     *     *ANP)*1.884E3/W(5,IND)**2
      IF(GAM.LT.0.) GAM=0.
      RETURN
   10 GAM=C(IND)*T**0.16667*ANE
      RETURN
      END
C
C ********************************************************************
C
      FUNCTION EPS(T,ANE,ALAM,ION,N)
C     ==============================
C
C     NLTE PARAMETER EPSILON (COLLISIONAL/SPONTANEOUS DEEXCITATION)
C     AFTER  KASTNER, 1981, J.Q.S.R.T. 26, 377
C
      INCLUDE 'PARAMS.FOR'
      DATA CK0,CK1 /7.75E-8, 2.58E-8/
      X=1.438E8/ALAM/T
      XKT=12390./ALAM
      TT=0.75*X
      T1=TT+1.
      A=4.36E7*XKT*XKT/(1.-EXP(-X))
      IF(ION.EQ.1) GO TO 10
      B=1.1+LOG(T1/TT)-0.4/T1/T1
      C=X*B*SQRT(T)/XKT/XKT*ANE
      IF(N.EQ.0) C=CK0*C
      IF(N.NE.0) C=CK1*C
      GO TO 20
   10 C=2.16/T/SQRT(T)/X**1.68*ANE
   20 EPS=C/(C+A)
      RETURN
      END
C
C ********************************************************************
C
      FUNCTION XK2DOP(TAU)
C     ====================
C
C     KERNEL FUNCTION K2  (AUXILIARY PROCEDURE TO NLTE)
C     AFTER  HUMMER,  1981, J.Q.S.R.T. 26, 187
C
      INCLUDE 'PARAMS.FOR'
      DATA PI2SQ,PISQ /2.506628275D0,  1.772453851D0/
      DATA A0,A1,A2,A3,A4 /
     *  1.D0,  -1.117897000D-1,  -1.249099917D-1,  -9.136358767D-3,
     *         -3.370280896D-4/
      DATA B0,B1,B2,B3,B4,B5 /
     *  1.D0,   1.566124168D-1,   9.013261660D-3,   1.908481163D-4,
     *         -1.547417750D-7,  -6.657439727D-9/
      DATA C0,C1,C2,C3,C4 /
     *  1.0D0,   1.915049608D01,   1.007986843D02,   1.295307533D02,
     *         -3.143372468D01/
      DATA D0,D1,D2,D3,D4,D5/
     *  1.D0,   1.968910391D01,   1.102576321D02,   1.694911399D02,
     *         -1.669969409D01,  -3.666448000D01/
      XK2DOP=1.D0
      IF(TAU.LE.0.) RETURN
      IF(TAU.GT.11.) GO TO 10
      P=A0+TAU*(A1+TAU*(A2+TAU*(A3+TAU*A4)))
      Q=B0+TAU*(B1+TAU*(B2+TAU*(B3+TAU*(B4+TAU*B5))))
      XK2DOP=TAU/PI2SQ*LOG(TAU/PISQ)+P/Q
      RETURN
   10 X=1.D0/LOG(TAU/PISQ)
      P=C0+X*(C1+X*(C2+X*(C3+X*C4)))
      Q=D0+X*(D1+X*(D2+X*(D3+X*(D4+X*D5))))
      XK2DOP=P/Q/2.D0/TAU/SQRT(LOG(TAU/PISQ))
      RETURN
      END
C
C ********************************************************************
C
      SUBROUTINE INKUR
C     ================
C
C     Input of a Kurucz model atmosphere
C
C     Input values (extracted from the Kurucz files):
C      TEF, G  - effective temperature, log g (appears only in output)
C      ND      - number of depth points
C     and for each depth:
C      DM      - m, m is the mass depth coordinate
C      T       - temperature
C      P       - gass pressure
C      ANE     - electron density
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION POP(MLEVEL),ES(MLEVEL,MLEVEL),BS(MLEVEL),POPLTE(MLEVEL)
      COMMON POP,ES,BS
C
      READ(8,501) TEF,GRAV
      READ(8,502) ND
      ND=ND-1
  501 FORMAT(4X,F8.0,9X,F8.5)
c  502 FORMAT(/////////////////////10X,I3)
  502 FORMAT(/////////////////////10X,I3/)
      WRITE(6,600) TEF,GRAV
      DO 10 ID=1,ND
         READ(8,*) DM(ID),TEMP(ID),P,ELEC(ID)
         AN=P/TEMP(ID)/BOLK
         DENS(ID)=WMM(ID)*(AN-ELEC(ID))
         WRITE(6,601) ID,DM(ID),TEMP(ID),ELEC(ID),DENS(ID)
         T=TEMP(ID)
         IF(IFMOL.GT.0.AND.T.LT.TMOLIM) THEN
c           AN=TOTN(ID)
            AEIN=ELEC(ID)
            CALL MOLEQ(ID,T,AN,AEIN,ANE,1)
          ELSE
            DO IAT=1,NATOM
               ATTOT(IAT,ID)=DENS(ID)/WMM(ID)/YTOT(ID)*ABUND(IAT,ID)
            END DO
         END IF
c        WRITE(6,601) ID,DM(ID),TEMP(ID),ELEC(ID),DENS(ID)
         CALL WNSTOR(ID)
         CALL SABOLF(ID)
         CALL RATMAT(ID,ES,BS)
         CALL LEVSOL(ES,BS,POPLTE,NLEVEL) 
         DO J=1,NLEVEL
            POPUL(J,ID)=POPLTE(J)
         END DO
   10 CONTINUE
      WRITE(77,503) ND, 3
      WRITE(77,504) (DM(ID),ID=1,ND)
      DO ID=1,ND
         WRITE(77,504) TEMP(ID),ELEC(ID),DENS(ID)
      END DO
c
      CLOSE(8)
c
  503 FORMAT(2I5)
  504 FORMAT(1P6E13.6)
  600 FORMAT(' INPUT KURUCZ MODEL FOR TEFF=',F7.0,'   LOG G =',
     *       F7.2//1H ,7X,'MASS',9X,'T',9X,'NE',9X,'DENS'/
     *       '-----------------------------------------------'/)
  601 FORMAT(1H ,I5,1PE10.3,0PF10.1,1P2E12.3)
      RETURN
      END
C
C ********************************************************************
C
C
C
      SUBROUTINE INPMOD
C     =================
C
C     Read an initial model atmosphere from unit 8
C     File 8 contains:
C      1. NDPTH -  number of depth points in which the initial model is
C                  given (if not equal to ND, routine interpolates
C                  automatically to the set DM by linear interpolation
C                  in log(DM)
C         NUMPAR - number of input model parameters in each depth
C                  = 3 for LTE model - ie. N, T, N(electron);
C                  > 3 for NLTE model)
C      2. DEPTH(ID),ID=1,NDPTH - mass-depth points for the input model
C      3. for each depth:
C                 T   - temperature
C                 ANE - electron density
C                 RHO - mass density
C                 level populations - only for NLTE input model
C                       Number of input level populations need not be
C                       equal to NLEVEL; in that case the procedure
C                       CHANGE is called from START to calculate the
C                       remaining level populations
C
C     Note: The output file 7, which is created by this program
C           (procedure OUTPUT) has the same structure as file 8
C           and may thus be used as input to another run of the
C           program
C     INTRPL - switch indicating whether (and, if so, how) interpolate
C              the initial model if the depth scales for the input model
C              and the present depth scale are different
C            = 0  -  no interpolation, i.e. scale DEPTH coincides with DM
C            > 0  -  polynomial interpolation of the (INTRPL-1)th order
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (MINPUT=MLEVEL+4)
      DIMENSION ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),POPLTE(MLEVEL),
     *          TOTN(MDEPTH),PLTE(MLEVEL,MDEPTH)
      COMMON ESEMAT,BESE,POPLTE,POPUL0(MLEVEL,MDEPTH),X(MINPUT),
     *       TEMP0(MDEPTH),ELEC0(MDEPTH),DENS0(MDEPTH),PPL0(MDEPTH),
     *       PPL(MDEPTH),DEPTH(MDEPTH),DM0(MDEPTH),DP(MDEPTH)
      COMMON/NLTPOP/PNLT(MATOM,MION,MDEPTH)
      common/quasex/iexpl(mlevel),iltot(mlevel)
C
      NUMLT=3
      IF(INMOD.EQ.2) NUMLT=4
      READ(8,*) NDPTH,NUMPAR
      READ(8,*) (DEPTH(I),I=1,NDPTH)
      ND=NDPTH
      NUMP=ABS(NUMPAR)
      DO 30 ID=1,NDPTH
         READ(8,*) (X(I),I=1,NUMP)
         TEMP(ID)=X(1)
         ELEC(ID)=X(2)
         DENS(ID)=X(3)
         TOTN(ID)=DENS(ID)/WMM(ID)+ELEC(ID)
         CALL WNSTOR(ID)
         CALL SABOLF(ID)
         IP=NUMLT
         IF(NUMPAR.LT.0) THEN
            IP=IP+1
            TOTN(ID)=X(IP)
         END IF
         IF(INMOD.EQ.2) IP=IP+1
c
c        first compute LTE level populations for all levels,
c        i.e. explicit, semi-explisit, and quasi-explicit
c
         NLEV0=NLEVEL
         TEMP(ID)=X(1)
         ELEC(ID)=X(2)
         DENS(ID)=X(3)
         t=temp(id)
         if(ifmol.gt.0.and.t.lt.tmolim) then
            ipri=1
            aein=elec(id)
            an=totn(id)
            call moleq(id,t,an,aein,ane,ipri)
          else
            if(imode.gt.-2) then
               DO IAT=1,NATOM
                  ATTOT(IAT,ID)=DENS(ID)/WMM(ID)/YTOT(ID)*ABUND(IAT,ID)
               END DO
             else
               DO IAT=1,NATOM
                  ATTOT(IAT,ID)=DENS(ID)/WMM(1)/YTOT(1)*ABUND(IAT,1)
               END DO
            end if
         end if
         CALL WNSTOR(ID)
         CALL SABOLF(ID)
         CALL RATMAT(ID,ESEMAT,BESE)
         CALL LEVSOL(ESEMAT,BESE,POPLTE,NLEV0) 
         DO I=1,NLEV0
            POPUL(I,ID)=POPLTE(I)
            PLTE(I,ID)=POPLTE(I)
c           if(id.eq.1) write(6,651) i,ip,popul(i,id),plte(i,id)
         END DO
c
c        if the input file fort.8 contains also NLTE level populations
c        of b-factors, replace the LTE populations by those
c
         IF(NUMP.GT.IP) THEN
            NLEV0=NUMP-IP    
            DO I=1,NLEV0
               j=iltot(i)
               POPUL(J,ID)=X(IP+I)*RELAB(IATM(I),ID)
c              if(id.eq.1) write(6,651) i,j,x(ip+i),popul(i,id)
c 651          format('in',2i4,1p2e12.4)
            END DO
c           DO I=1,NLEV0
c              j=iltot(i)
c              if(popul(j,id).le.0.) then
c                 IE=IEL(I)
c                 N0I=NFIRST(IE)
c                 NKI=NNEXT(IE)
c                 POPUL(J,ID)=ELEC(ID)*POPUL(iltot(NKI),ID)*SBF(I)
c              end if
c           END DO
c
c           in the case the input "NLTE populations are in fact b-factors,
c           compute the real populations
c
            if(ibfac.eq.1) then
               do i=1,nlev0
                  j=iltot(i)
                  popul(j,id)=popul(j,id)*plte(j,id)
               end do
            end if
         END IF
   30 CONTINUE
C
      close(8)
c
      write(6,600)
  600 format(/' INPUT TLUSTY MODEL'/
     *        ' ------------------'/
     *         1H ,8X,'MASS',9X,'T',9X,'NE',9X,'DENS'//)
         nd=ndpth
         DO 40 ID=1,ND
            DM(ID)=DEPTH(ID)
            write(6,601) id,dm(id),temp(id),elec(id),dens(id),
     *       popul(1,id)
  601       format(i6,1pe10.3,0pf10.1,1p4e12.3)
   40    CONTINUE
C
      DO 100 ID=1,ND
         BCON=ELEC(ID)/TEMP(ID)/SQRT(TEMP(ID))*2.0706E-16
         DO 100 IONE=1,NION
            ION=IZ(IONE)
            IAT=NUMAT(IATM(NFIRST(IONE)))
            NKI=NNEXT(IONE)
            IF(ION.GT.0) PNLT(IAT,ION,ID)=POPUL(NKI,ID)/G(NKI)*BCON
  100 CONTINUE
c
c     check abundances
c
c     CALL CHCKAB
      RETURN
      END

C
C ********************************************************************
C
C
 
      SUBROUTINE INPBF
C     ================
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (MINPUT=MLEVEL+4)
      DIMENSION DEPTH(MDEPTH),X(MINPUT,MDEPTH),XX(MDEPTH),BF(MDEPTH)
C
      OPEN(8,FILE='bfactors',STATUS='OLD')
      NUMLT=3
      IF(INMOD.EQ.2) NUMLT=4
      READ(8,*) NDPTH,NUMPAR
      READ(8,*) (DEPTH(I),I=1,NDPTH)
      IF(NUMPAR.LT.0) NUMLT=NUMLT+1
      NUMP=ABS(NUMPAR)
      DO ID=1,NDPTH
         READ(8,*) (X(I,ID),I=1,NUMP)
      END DO
      CLOSE(8)
c
c     interpolate the input b-factors to the original DM-scale;
c     compute new NLTE populations
c
      DO I=NUMLT+1,NUMP
         DO ID=1,NDPTH
            XX(ID)=X(I,ID)
         END DO
         CALL INTERP(DEPTH,XX,DM,BF,NDPTH,ND,2,1,1)
         DO ID=1,ND
            POPUL(I-NUMLT,ID)=POPUL(I-NUMLT,ID)*BF(ID)
         END DO
      END DO
C
      RETURN
      END

C
C
C     ****************************************************************
C
C
 
      SUBROUTINE LEVSOL(A,B,POPP,NLVCAL)
C     ==================================
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION A(MLEVEL,MLEVEL),B(MLEVEL),POPP(MLEVEL),
     *          AP(MLEVEL,MLEVEL),BP(MLEVEL),POPP1(MLEVEL)
C
C     new populations by inverting several partial rate matrices for the
C     individual chemical species
C
      if(nlvcal.le.0) return
         DO 50 IAT=1,NATOM
            N1=N0A(IAT)
            NK=NKA(IAT)
            IF(N1.LE.0) THEN
               DO 1 I=N0A(IAT),NKA(IAT)
                  N1=I
                  IF(I.GT.0) GO TO 2
    1          CONTINUE
    2          CONTINUE
            END IF
            IF(N1.LE.0) GO TO 50
            NLP=NK-N1+1
            DO 20 I=N1,NK
               DO 10 J=N1,NK
                  AP(I-N1+1,J-N1+1)=A(I,J)
   10          CONTINUE
               BP(I-N1+1)=B(I)
   20       CONTINUE
            CALL LINEQS(AP,BP,POPP1,NLP,MLEVEL)
            DO 30 I=N1,NK
                POPP(I)=POPP1(I-N1+1)
   30       CONTINUE
   50    CONTINUE
      RETURN
      END
 
C
C
C     ****************************************************************
C

      SUBROUTINE CHANGE
C     =================
C
C     This procedure controls an evaluation of initial level
C     populations in case where the system of explicit levels
C     (ie. the choice of explicit level, their numbering, or their
C     total number) is not consistent with that for the input level
C     populations read by procedure INPMOD.
C     Obviously, this procedure need be used only for NLTE input models.
C
C     Input from unit 5:
C     For each explicit level, II=1,NLEVEL, the following parameters:
C      IOLD   -  NE.0 - means that population of this level is
C                       contained in the set of input populations;
C                       IOLD is then its index in the "old" (i.e. input)
C                       numbering.
C                       All the subsequent parameters have no meaning
C                       in this case.
C             -  EQ.0 - means that this level has no equivalent in the
C                       set of "old" levels. Population of this level
C                       has thus to be evaluated.
C      MODE   -    indicates how the population is evaluated:
C             = 0  - population is equal to the population of the "old"
C                    level with index ISIOLD, multiplied by REL;
C             = 1  - population assumed to be LTE, with respect to the
C                    first state of the next ionization degree whose
C                    population must be contained in the set of "old"
C                    (ie. input) populations, with index NXTOLD in the
C                    "old" numbering.
C                    The population determined of this way may further
C                    be multiplied by REL.
C             = 2  - population determined assuming that the b-factor
C                    (defined as the ratio between the NLTE and
C                    LTE population) is the same as the b-factor of
C                    the level ISINEW (in the present numbering). The
C                    level ISINEW must have the equivalent in the "old"
C                    set; its index in the "old" set is ISIOLD, and the
C                    index of the first state of the next ionization
C                    degree, in the "old" numbering, is NXTSIO.
C                    The population determined of this way may further
C                    be multiplied by REL.
C             = 3  - level corresponds to an ion or atom which was not
C                    explicit in the old system; population is assumed
C                    to be LTE.
C      NXTOLD  -  see above
C      ISINEW  -  see above
C      ISIOLD  -  see above
C      NXTSIO  -  see above
C      REL     -  population multiplier - see above
C                 if REL=0, the program sets up REL=1
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION ESEMAT(MLEVEL,MLEVEL),BESE(MLEVEL),POPLTE(MLEVEL)
      COMMON ESEMAT,BESE,POPLTE,POPUL0(MLEVEL,MDEPTH),
     *       POPULL(MLEVEL,MDEPTH),POPL(MLEVEL)
C
      PARAMETER (S = 2.0706E-16)
      IFESE=0
      DO 100 II=1,NLEVEL
         READ(ICHANG,*) IOLD,MODE,NXTOLD,ISINEW,ISIOLD,NXTSIO,REL
         IF(MODE.GE.3) IFESE=IFESE+1
         IF(REL.EQ.0.) REL=1.
         DO 90 ID=1,ND
            IF(IOLD.EQ.0) GO TO 10
            POPUL0(II,ID)=POPUL(IOLD,ID)
            GO TO 90
   10       IF(MODE.NE.0) GO TO 20
            POPUL0(II,ID)=POPUL(ISIOLD,ID)*REL
            GO TO 90
   20       T=TEMP(ID)
            ANE=ELEC(ID)
            IF(MODE.GE.3) GO TO 40
            NXTNEW=NNEXT(IEL(II))
            SB=S/T/SQRT(T)*G(II)/G(NXTNEW)*EXP(ENION(II)/T/BOLK)
            IF(MODE.GT.1) GO TO 30
            POPUL0(II,ID)=SB*ANE*POPUL(NXTOLD,ID)*REL
            GO TO 90
   30       KK=ISINEW
            KNEXT=NNEXT(IEL(KK))
            SBK=S/T/SQRT(T)*G(KK)/G(KNEXT)*EXP(ENION(KK)/T/BOLK)
            POPUL0(II,ID)=SB/SBK*POPUL(NXTOLD,ID)/POPUL(NXTSIO,ID)*
     *                    POPUL(ISIOLD,ID)*REL
            GO TO 90
   40       IF(IFESE.EQ.1) THEN
               CALL SABOLF(ID)
               CALL RATMAT(ID,ESEMAT,BESE)
               CALL LINEQS(ESEMAT,BESE,POPLTE,NLEVEL,MLEVEL)
               DO 50 III=1,NLEVEL
   50             POPULL(III,ID)=POPLTE(III)
            END IF
            POPUL0(II,ID)=POPULL(II,ID)
   90    CONTINUE
  100 CONTINUE
      DO 110 I=1,NLEVEL
         DO 110 ID=1,ND
            POPUL(I,ID)=POPUL0(I,ID)
  110 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE RATMAT(ID,A,B)
C
C     LTE RATE MATRIX  (SAHA-BOLTZMANN EQS. + CHARGE CONSERVATION EQ.)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      parameter (un=1.)
      DIMENSION A(MLEVEL,MLEVEL),B(MLEVEL)
C
      ANE=ELEC(ID)
      DO I=1,NLEVEL
         B(I)=0.
         DO J=1,NLEVEL
            A(J,I)=0.
         END DO
      END DO
C
      DO IAT=1,NATOM
         N0I=N0A(IAT)
         NKI=NKA(IAT)
         N1I=NKI-1
         NREFI=NKI
         DO I=N0I,N1I
            A(I,I)=1.
            N=NNEXT(IEL(I))
            A(I,N)=-ANE*SBF(I)*WOP(I,ID)
         END DO
         DO I=N0I,NKI
            IL=ILK(I)
            A(NREFI,I)=UN
            IF(IL.NE.0) A(NREFI,I)=1.+ANE*USUM(IL)
         END DO
         B(NREFI)=ATTOT(IAT,ID)
      END DO
C
      RETURN
      END
C
C     ****************************************************************
C
      SUBROUTINE SABOLF(ID)
C     =====================
C
C     Saha-Boltzmann factors (SBF)
C     and "upper sums" - sum of Saha-Boltzmann factors for upper, LTE,
C     levels which are not included explicitly (USUM), and derivatives
C     wrt. temperature (T) and electron density (DUSUMN)
C
C     Input: ID  - depth index
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (UH=1.5)
      PARAMETER (CMAX=2.154D4,CCON=2.0706D-16,TWO=2.D0)
C
C     DCHI - approximate lowering of ionization potential for neutrals
C      Actual lowering is DCHI*effective charge, and is considered only
C      if IUPSUM(ION).GT.0
C
      T=TEMP(ID)
      SQT=SQRT(T)
      ANE=ELEC(ID)
      STANE=SQRT(T/ANE)
      XMAX=CMAX*SQRT(STANE)
      TK=BOLK*T
      CON=CCON/T/SQT
C
C     Saha-Boltzmann factors
C
      DO 50 ION=1,NION
         QZ=IZ(ION)
         CFN=CON/G(NNEXT(ION))
         DCH=0.
         IUPS=IUPSUM(ION)
         SSBF=0.
         USUM(ION)=0.
         nlst=nlast(ion)
         if(ifwop(nlst).ge.0) then
             nl1up=nquant(nlst)+1
          else
             nl1up=nquant(nlst)
         end if
         DO 10 II=NFIRST(ION),NLAST(ION)
            if(ifwop(ii).lt.0) then
               E=EH*QZ*QZ/TK
               SUM=0.
               DO 5 J=nl1up,NLMX
                  XJ=J
                  XI=J*J
                  X=E/XI
                  FI=XI*EXP(X)*WNHINT(J,ID)
                  SUM=SUM+FI
    5          CONTINUE
               g(ii)=sum*two
               gmer(imrg(ii),id)=g(ii)
            end if
            X=ENION(II)/TK
            if(x.gt.110.) x=110.
            SB=CFN*G(II)*EXP(X)
            SBF(II)=SB
            SSBF=SSBF+SB
   10    CONTINUE
C
C     Upper sums
C
         if(ifwop(nlst).lt.0) go to 50
         if(iups.eq.0) then
C
C     1. More exact approach - using (exact) partition functions
C
         IAT=NUMAT(IATM(NFIRST(ION)))
         XMX=XMAX*SQRT(QZ)
         CALL PARTF(IAT,IZ(ION),T,ANE,XMX,U)
         EE=ENION(NFIRST(ION))/TK
         if(ee.gt.110.) ee=110.
         CFE=CFN*EXP(EE)
         USUM(ION)=CFE*U-SSBF
         xx=(ssbf-sbf(nfirst(ion)))/sbf(nfirst(ion))
         IF(USUM(ION).LT.0.or.ee.ge.109.or.xx.lt.1.e-7) USUM(ION)=0.
         IF(USUM(ION).LT.0.) USUM(ION)=0.
C
C     2. Approximate approach - summation over fixed number of upper
C        levels, assumed hydrogenic (ie. their ionization energy and
C        statistical weight hydrogenic)
C
         else if(iups.gt.0) then
         SUM=0.
         DSUM=0.
         E=EH*QZ*QZ/TK
         DO 30 J=NQUANT(NLAST(ION))+1,IUPS
            XI=J*J
            X=E/XI
            FI=XI*EXP(X)
            SUM=SUM+FI
   30    CONTINUE
         USUM(ION)=SUM*CON*TWO
C
c        3. occupation probability form
c
         else 
         SUM=0.
         DSUM=0.
         E=EH*QZ*QZ/TK
         DO 40 J=NQUANT(NLAST(ION))+1,NLMX
            XJ=J
            XI=J*J
            X=E/XI
            FI=XI*EXP(X)*WNHINT(J,ID)
            SUM=SUM+FI
   40    CONTINUE
         USUM(ION)=SUM*CON*TWO
         end if
   50 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C

      FUNCTION SBFHMI_old(FR)
C     ===================
C
C     Bound-free cross-section for H- (negative hydrogen ion)
C
      INCLUDE 'PARAMS.FOR'
      SBFHMI=0.
      sbfhmi_old=0.
      FR0=1.8259E14
      IF(FR.LT.FR0) RETURN
      IF(FR.LT.2.111E14) GO TO 10
      X=2.997925E15/FR
      SBFHMI=(6.80133E-3+X*(1.78708E-1+X*(1.6479E-1+X*(-2.04842E-2+X*
     1        5.95244E-4))))*1.E-17
      sbfhmi_old=sbfhmi
      RETURN
   10 X=2.997925E15*(1./FR0-1./FR)
      SBFHMI=(2.69818E-1+X*(2.2019E-1+X*(-4.11288E-2+X*2.73236E-3)))
     1       *X*1.E-17
      sbfhmi_old=sbfhmi
      RETURN
      END
C
C
C
C     ****************************************************************
C
C

      SUBROUTINE OPADD(MODE,ID,FR,ABAD,EMAD,SCAD)
C     ===========================================
C
C     Additional opacities
C     This is basically user-supplied procedure; here are some more
C     important non-standard opacity sources, namely
C     Rayleigh scattering, H- opacity, H2+ opacity, and additional
C     opacity of He I and He II.
C     Inclusion of these opacities is contolled by switches transmitted
C     by COMMON/OPCPAR - see description in START.
C
C     Input parameters:
C     MODE  - controls the nature and the amount of calculations
C           = -1 - (OPADD called from START) evaluation of relevant
C                  depth-dependent quantities (usually photoionization
C                  cross-sections, but also possibly other), which are
C                  stored in array CROS
C           = 0  - evaluation of an additional opacity, emissivity, and
C                  scattering - for procedure OPAC0
C     ID    - depth index
C     FR    - frequency 
C
C     Output: 
C
C     ABAD  - absorption coefficient (at frequency FR and depth ID)
C     EMAD  - emission coefficient (at frequency FR and depth ID)
C     SCAD  - scattering coefficient (at frequency FR and depth ID)
C
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (FRAYH  =  2.463E15,
     *           FRAYHe =  5.150E15,
     *           FRAYH2 =  2.922E15,
     *           CLS    =  2.997925e18)
C
      AB0=0.
      AB1=0.
      ABAD=0.
      EMAD=0.
      SCAD=0.
C
      if(iath.gt.0) then
      N0HN=NFIRST(IELH)
      NKH=NKA(IATH)
C
      IF(MODE.GE.0) THEN
         T=TEMP(ID)
         ANE=ELEC(ID)
         HKT=HK/T
         T32=1./T/SQRT(T)
      END IF
      anh=dens(id)/(wmm(id)*ytot(id))
      anhe=rrr(id,1,2)
C
      IT=NLEVEL
C
C   -----------------------
C   HI  Rayleigh scattering
C   -----------------------
C
      IF(IRSCT.NE.0.AND.IOPHLI.NE.1.AND.IOPHLI.NE.2) THEN
         X=1.D0/(CLS/MIN(FR,FRAYH))**2
         SG=(5.799E-13+(1.422E-6+2.784*X)*X)*X*X
c        ABAD=POPUL(N0HN,ID)*SG 
         SCAD=POPUL(N0HN,ID)*SG
         scad=anh*sg
      END IF
      IF(IOPHMI.NE.0) THEN
C
C   ----------------------------
C   H-  bound-free and free-free
C   ----------------------------
C     Note: IOPHMI must not by taken non-zero if H- is considered
C           explicitly, because H- opacity would be taken twice
C
          SG=SBFHMI(FR)
          XHM=8762.9/T
          SB=1.0353E-16*T32*EXP(XHM)*POPUL(N0HN,ID)*ANE*SG
          SF=SFFHMI(POPUL(N0HN,ID),FR,T)*ANE
          AB0=SB+SF
      END IF
C
C   -----------------------
C   He I  Rayleigh scattering
C   -----------------------
C
      IF(IRSCHE.NE.0.AND.MODE.GE.0) THEN
         X=(CLS/MIN(FR,FRAYHe))**2
         CS=5.484E-14/X/X*(1.+(2.44E5+5.94E10/(X-2.90E5))/X)**2
         sg=anhe*cs
c        abad=abad+sg
         scad=scad+sg
      END IF
C
C   -----------------------
C   H2  Rayleigh scattering
C   -----------------------
C
      IF(IRSCH2.NE.0.AND.MODE.GE.0.AND.IFMOL.GT.0) THEN
         X=(CLS/MIN(FR,FRAYH2))**2
           X2=1./X/X
           CS=(8.14E-13+1.28E-6/X+1.61*X2)*X2
           sg=cs*anh2(id)
c          abad=abad+sg
           scad=scad+sg
        END IF
C
      IF(IOPH2P.GT.0.AND.IFMOL.GT.0) THEN
C
C   -----------------------------
C   H2+  bound-free and free-free
C   -----------------------------
C
         X=FR*1.E-15
         SG1=(-7.342E-3+(-2.409+(1.028+(-4.23E-1+
     *       (1.224E-1-1.351E-2*X)*X)*X)*X)*X)*1.602E-12/BOLK
         IT=IT+1
         X=LOG(FR)
         SG2=-3.0233E3+(3.7797E2+(-1.82496E1+(3.9207E-1-
     *        3.1672E-3*X)*X)*X)*X
         X2=-SG1/T+SG2
         SB=0.
         IF(X2.GT.-150.) SB=POPUL(N0HN,ID)*POPUL(NKH,ID)*EXP(X2)
         AB0=AB0+SB
      END IF
      end if
C
C   -----------------------------
C   He-  free-free
C   -----------------------------
C
      if(mode.ge.0.and.iophem.gt.0) then
         A=3.397D-46+(-5.216D-31+7.039D-15/FR)/FR
         B=-4.116D-42+(1.067D-26+8.135D-11/FR)/FR
         C=5.081D-37+(-8.724D-23-5.659D-8/FR)/FR
         cs=a*t+b+c/t
         sg=anhe*ane*cs
         ab0=ab0+sg
      end if
C
C   -----------------------------
C   H2-  free-free
C   -----------------------------
C
      IF(IOPH2M.NE.0.AND.MODE.GE.0.AND.IFMOL.GT.0.AND.T.LT.TMOLIM) THEN
         call h2minus(t,anh2(id),ane,fr,oph2)
         ab1=ab1+oph2
      END IF
C
C   -----------------------------
C     CH and OH continuuum opacity
C   -----------------------------
C 
      if(mode.ge.0.and.ifmol.gt.0.and.t.lt.tmolim) then
         if(iopch.gt.0) ab0=ab0+sbfch(fr,t)*anch(id)
         if(iopoh.gt.0) ab0=ab0+sbfoh(fr,t)*anoh(id)
C
C     ---------------------------
C     CIA H2-H2 opacity
C     ---------------------------
C
         if(ioh2h2.gt.0) then
            call cia_h2h2(t,anh2(id),fr,oph2)
            ab1=ab1+oph2
         end if
C
C     ---------------------------
C     CIA H2-He opacity
C     ---------------------------
C
         if(ioh2he.gt.0) then
            call cia_h2he(t,anh2(id),anhe,fr,oph2)
            ab1=ab1+oph2
         end if
C
C     ---------------------------
C     CIA H2-H opacity
C     ---------------------------
C
         if(ioh2h1.gt.0) then
            call cia_h2h(t,anh2(id),anh,fr,oph2)
            ab1=ab1+oph2
         end if
C
C     ---------------------------
C     CIA H-He opacity
C     ---------------------------
C  
         if(iohhe.gt.0) then
            call cia_hhe(t,anh,anhe,fr,oph2)
            ab1=ab1+oph2
         end if
      end if
C
C     ----------------------------------------------
C     The user may supply more opacity sources here:
C     ----------------------------------------------
C
C     Finally, actual absorption and emission coefficients

      IF(MODE.LT.0) RETURN
      X=EXP(-HKT*FR)
      X1=1.-X
      BNX=BN*(FR*1.E-15)**3*X
      ABAD=ABAD+X1*AB0+AB1
      EMAD=EMAD+BNX*(AB0+AB1/X1)
      RETURN
      END
C
C
C     ****************************************************************
C
C
      function wn(xn,a,id,z)
C     ======================
c
c     evaluation of the occupation probablities for a hydrogenic ion
c     using eqs (4.26), and (4.39) of Hummer,Mihalas Ap.J. 331, 794, 1988.
c     approximate evaluation of Q(beta) - Hummer
c
c     Input: xn  - real number corresponding to quantum number n
c            a   - correlation parameter
c            id  - depth index
c            z   - ionic charge
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      parameter (p1=0.1402,p2=0.1285,p3=1.,p4=3.15,p5=4.,un=1.)
      parameter (tkn=3.01,ckn=5.33333333,cb=8.59e14)
      parameter (f23=-2./3.)
      parameter (a0=0.529177e-8,wa0=-3.1415926538/6.*a0*a0*a0)
c
c     evaluation of k(n)
c
      if(xn.le.tkn) then
         xkn=un
       else
         xkn=ckn*xn/(xn+un)/(xn+un)
      end if
c
c     evaluation of beta
c
c     beta=cb*bergfc*z*z*z*xkn/(xn*xn*xn*xn)*exp(f23*log(elec(id)))
      beta=cb*z*z*z*xkn/(xn*xn*xn*xn)*exp(f23*log(elec(id)))
c
c     approximate expression for Q(beta)
c
      x=exp(p4*log(un+p3*a))
c     c1=p1*(x+p5*z*a*a*a)    ! previous expression -ERROR !!!!!!
      c1=p1*(x+p5*(z-un)*a*a*a)
      c2=p2*x
      f=(c1*beta*beta*beta)/(un+c2*beta*sqrt(beta))
      wp=f/(un+f)
c
c     contribution from neutral particles
c
      xn2=xn*xn+un
      xnh=0.
      xnhe1=0.
      if(ielh.gt.0) xnh=popul(nfirst(ielh),id)
      if(ielhe1.gt.0) xnhe1=popul(nfirst(ielhe1),id)
      w0=exp(wa0*xn2*xn2*xn2*(xnh+xnhe1))
      W0=1.
      wn=wp*w0
      return
      end
C
C
C ********************************************************************
C
C
      SUBROUTINE WNSTOR(ID)
C     =====================
C
C     Stores occupation probabilities for hydrogen levels
C     in common WNCOM for further use
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (UN=1.,TWO=2.,SIXTH=1./6.,CCOR=0.09)
      parameter (p1=0.1402,p2=0.1285,p3=1.,p4=3.15,p5=4.)
      parameter (tkn=3.01,ckn=5.33333333,cb=8.59d14,f23=-2./3.)
C
      ANE=ELEC(ID)
      A=CCOR*EXP(SIXTH*LOG(ANE))/SQRT(TEMP(ID))
      DO 20 I=1,NLMX
         XN=I
         WNHINT(I,ID)=wn(xn,a,id,un)
         WNHE2(I,ID)=wn(xn,a,id,two)
   20 CONTINUE
C
C     array WOP - occupation probabilities for explicit levels
C
      do 30 ii=1,nlevel
         wop(ii,id)=un
         if(ifwop(ii).le.0) go to 30
         ie=iel(ii)
         nq=nquant(ii)
         if(iz(ie).eq.1) then
            wop(ii,id)=wnhint(nq,id)
          else if(iz(ie).eq.2) then
            wop(ii,id)=wnhe2(nq,id)
          else
            z=iz(ie)
            xn=nq
            wop(ii,id)=wn(xn,a,id,z)
         end if
   30 continue
      RETURN
      END
C
C
C ********************************************************************
C
C
C
c      SUBROUTINE TIMING(MOD,ITER)
C     ===========================
C
C     Timing procedure (call machine dependent routine!!)
C
C     INCLUDE 'PARAMS.FOR'
c      CHARACTER ROUT*20
c      dimension dummy(2)
c      DATA T0/0./
C
c      TIME=etime(dummy)
c      DT=TIME-T0
c      T0=TIME
c      IF(MOD.EQ.0) THEN
c         IP=0
c         ROUT=' INIT     '
c      ELSE IF(MOD.EQ.1) THEN
c         IP=ITER-1
c         ROUT=' OPACITY  '
c      ELSE IF(MOD.EQ.2) THEN
c         IP=ITER
c         ROUT='  TRANSFER'
c      ENDIF
c      WRITE(69,600) IP,MOD,TIME,DT,ROUT
c  600 FORMAT(2I4,2F11.2,2X,A10)
c      RETURN
c      END
c
C *******************************************************************
c
      subroutine quit(text)
C     =====================
c
c     stops the program and writes a text
c
      INCLUDE 'PARAMS.FOR'
      character*(*) text
      write(6,10) text
   10 format(1x,a)
      stop
      end
c
c

C
C *******************************************************************
C


      function voigte(a,vs)
c     =====================
c
c     computes a voigt function  h = h(a,v)
c     a=gamma/(4*pi*dnud)   and  v=(nu-nu0)/dnud.  this  is  done after
c     traving (landolt-b\rnstein, p. 449).
c
      INCLUDE 'PARAMS.FOR'
      dimension ak(19),a1(5)
      data ak      /-1.12470432, -0.15516677,  3.28867591, -2.34357915,
     ,  0.42139162, -4.48480194,  9.39456063, -6.61487486,  1.98919585,
     , -0.22041650, 0.554153432, 0.278711796,-0.188325687, 0.042991293,
     ,-0.003278278, 0.979895023,-0.962846325, 0.532770573,-0.122727278/
      data sqp/1.772453851/,sq2/1.414213562/
c
      v = abs(vs)
      u = a + v
      v2 = v*v
      if (a.eq.0.0) go to 140
      if (a.gt.0.2) go to 120
      if (v.ge.5.0) go to 121
c
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 1
c
  100 quo = 1.
      if (v.lt.2.4) go to 101
      quo = 1./(v2 - 1.5)
      m = 11
      go to 102
c
  101 m = 6
      if (v.lt.1.3) m = 1
  102 do 103 i=1,5
         a1(i) = ak(m)
         m = m + 1
  103 continue
      h1 = quo*(a1(1) + v*(a1(2) + v*(a1(3) + v*(a1(4) + v*a1(5)))))
      if (k.gt.1) go to 110
c
c a le 0.2  and v lt 5.
c
      hh = h1*a + ex*(1. + a*a*(1. - 2.*v2))
      voigte=hh
      return
c
  110 pqs = 2./sqp
      h1p = h1 + pqs*ex
      h2p = pqs*h1p - 2.*v2*ex
      h3p = (pqs*(1. - ex*(1. - 2.*v2)) - 2.*v2*h1p)/3. + pqs*h2p
      h4p = (2.*v2*v2*ex - pqs*h1p)/3. + pqs*h3p
      psi = ak(16) + a*(ak(17) + a*(ak(18) + a*ak(19)))
c
c 0.2 lt a le 1.4  and  a + v le 3.2
c
      hh = psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
      voigte=hh
      return
c
  120 if (a.gt.1.4.or.u.gt.3.2) go to 130
      ex=0.
      if(v2.lt.100.)ex = exp(-v2)
      k = 2
      go to 100
c
c a le 0.2  and  v ge 5.
c
  121 hh = a*(15. + 6.*v2 + 4.*v2*v2)/(4.*v2*v2*v2*sqp)
      voigte=hh
      return
c
  130 a2 = a*a
      u = sq2*(a2 + v2)
      u2 = 1./(u*u)
c
c a gt 1.4  or  a + v gt 3.2
c
      hh = sq2/sqp*a/u*(1. + u2*(3.*v2 - a2) +
     ,        u2*u2*(15.*v2*v2 - 30.*v2*a2 + 3.*a2*a2))
      voigte=hh
      return
c
c a eq 0.
c
  140 hh=0.
      if(v2.lt.100.) hh=exp(-v2)
      voigte=hh
      return
      end
C
C
C ********************************************************************
C
C
      SUBROUTINE SIGAVS
C     =================
C
C     Read bound-free cross-sections for averaged levels
C     from the unit INSA (given by IFANCY), with increasing frequencies
C     It assumes that all continuum transitions for a given ion are
C     given in a successive order in the data (i.e. as in TLUSTY for
C     explicit levels. For other levels, additional input data in
C     unit 54 !!
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (HCCM=H*2.997925D10,BAM=1.D-18)
      DIMENSION CRD(MFCRA),XIFE(8),FRD(MFCRA)
      CHARACTER*40 FIDATA(MION),FIODF1(MION),FIODF2(MION),FIBFCS(MION)
      COMMON/IONFIL/FIDATA,FIODF1,FIODF2,FIBFCS
C
      DATA XIFE/63480.,130563.,247220.,442000.,605000.,799000.,
     &          1008000.,1218380./
C
      FR1=FREQ(1)
      FR2=FREQ(2)
      NUNIT=0
      NQHT=0
      IF(IASV.EQ.0) GOTO 100
c     WRITE(6,600)
c 600 FORMAT(///,' DETAILED PHOTOIONIZATION CROSS-SECTIONS',
c    * ' (EXPLICIT LEVELS)',/,
c    * ' ---------------------------------------',/)
      DO 10 I=1,NION
        N1=NFIRST(I)
        N2=NLAST(I)
        INSA=0
        DO 11 II=N1,N2
          NFCR(II)=2
          FRECR(II,1)=FR1
          FRECR(II,2)=FR2
          CROSR(II,1)=0.
          CROSR(II,2)=0.
          INSB=IBF(II)
          IF(INSB.LT.50.OR.INSB.GT.100) GO TO 11
          IF(INSA.EQ.0) INSA=INSB
          IF(INSA.NE.INSB)
     *     call quit(' Incoherent file units in SIGAVS')
   11   CONTINUE
        IF(INSA.EQ.0) GOTO 10
        IF(FIBFCS(I).NE.' ') THEN
           INSA=INBFCS(I)
           OPEN(INSA,FILE=FIBFCS(I),STATUS='OLD')
        END IF
        READ(INSA,*,END=500,ERR=500) IIAT,IIZ,NSUP
        ATI=IIAT+0.01*(IIZ-1)
        NBFI=NSUP
        IF(NSUP.GT.(N2-N1+1)) NBFI=(N2-N1+1)
c    *   call quit(' Too many bf-trans. in input file (SIGAVS)')
c       WRITE(6,601) ATI,INSA
        DO 12 II=1,NBFI
          READ(INSA,*,END=500,ERR=500) IILO,EELO,GGLO,NFCRR
          IK=N1+IILO-1
          IF (IK.GT.N2 .OR. IK.LT.N1)
     *      call quit(' Inconsistent level numbering in SIGAVS')
          IF(IIAT.NE.26) GOTO 13
          ECMR=XIFE(IIZ)-EELO
c         DE=ABS((ENION(IK)-HCCM*ECMR)/ENION(IK))
c         IF(DE.GT.1.D-4) call quit(' Incorrect energy level in SIGAVS')
   13     READ(INSA,*,END=500,ERR=500) FR0,CR0
          NFD=1
          FRD(NFD)=FR0
          CRD(NFD)=CR0
          LUV=.FALSE.
          DO 14 IJ=1,NFCRR-1
            READ(INSA,*,END=500,ERR=500) FRIN,CRIN
            IF(LUV) GOTO 14
            IF(FRIN.GT.FR1) THEN
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1
                FRD(NFD)=FR0
                CRD(NFD)=CR0
              ENDIF
              NFD=NFD+1
              FRD(NFD)=FRIN
              CRD(NFD)=CRIN
              LUV=.TRUE.
            ELSE IF(FRIN.GT.FR2) THEN
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1
                FRD(NFD)=FR0
                CRD(NFD)=CR0
              ENDIF
              NFD=NFD+1
              FRD(NFD)=FRIN
              CRD(NFD)=CRIN
              FR0=FRIN
              CR0=CRIN
            ELSE
              FR0=FRIN
              CR0=CRIN
            ENDIF
            IF(NFD.GT.MFCRA)
     *        call quit(' Too many frequencies in SIGAVS')
   14     CONTINUE
          CRMX(IK)=0.
          DO 15 IJ=1,NFD
            CRMX(IK)=MAX(CRMX(IK),CRD(IJ))
   15     CONTINUE
          IF(CRMX(IK).GT.0.) THEN
c           WRITE(6,601) ATI,IILO,EELO,NFD
c 601       FORMAT(F7.2,I6,F13.3,I8)
            NFCR(IK)=NFD
            DO 16 IJ=1,NFD
              FRECR(IK,IJ)=FRD(NFD-IJ+1)
              CROSR(IK,IJ)=CRD(NFD-IJ+1)*BAM
   16       CONTINUE
          ENDIF
   12   CONTINUE
   10 CONTINUE
C
  100 READ(54,*,END=540,ERR=540) NUNIT
      IF(NUNIT.LE.0) RETURN
      WRITE(6,602)
  602 FORMAT(///,' DETAILED PHOTOIONIZATION CROSS-SECTIONS',
     * ' (NON-EXPLICIT LEVELS)',/,
     * ' ---------------------------------------',/)
      DO 110 IN=1,NUNIT
        READ(54,*,END=540,ERR=540) ATIR,INSA,NQHTR
        NQHT=NQHT+NQHTR
        IF(NQHT.GT.MPHOT)
     *    call quit(' Too many BF cross-sections in SIGAVS')
        READ(INSA,*,END=501,ERR=501) IIAT,IIZ,NSUP
C
c       check the total number of superlevels
c
        IF(NQHTR.GT.NSUP) THEN
           WRITE(6,603) NQHTR,NSUP
  603      FORMAT(' NQHTR=',i4,' in Unit 54 input greater than NSUP=',
     *              i4,/' program resets NQHTR to NSUP'/)
           NQHTR=NSUP
        END IF
c
C      loop over superlevels - read cross-sections
c
        DO 120 I=1,NQHTR
          IK=NQHT-NQHTR+I
          READ(INSA,*,END=501,ERR=501) IILO,EELO,GGLO,NFCRR
          AQHT(IK)=ATIR
          EQHT(IK)=EELO
          GQHT(IK)=GGLO
          READ(INSA,*) FR0,CR0
          NFD=1
          FRD(NFD)=FR0
          CRD(NFD)=CR0
          LUV=.FALSE.
          DO 130 IJ=1,NFCRR-1
            READ(INSA,*) FRIN,CRIN
            IF(LUV) GOTO 130
            IF(FRIN.GT.FR1) THEN
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1
                FRD(NFD)=FR0
                CRD(NFD)=CR0
              ENDIF
              NFD=NFD+1
              FRD(NFD)=FRIN
              CRD(NFD)=CRIN
              LUV=.TRUE.
            ELSE IF(FRIN.GT.FR2) THEN
              IF(FR0.LE.FR2.AND.IJ.GT.1) THEN 
                NFD=NFD+1
                FRD(NFD)=FR0
                CRD(NFD)=CR0
              ENDIF
              NFD=NFD+1
              FRD(NFD)=FRIN
              CRD(NFD)=CRIN
              FR0=FRIN
              CR0=CRIN
            ELSE
              FR0=FRIN
              CR0=CRIN
            ENDIF
  130     CONTINUE
          CRMY(IK)=0.
          DO 140 IJ=1,NFD
            CRMY(IK)=MAX(CRMY(IK),CRD(IJ))
  140     CONTINUE
          IF(CRMY(IK).GT.0.) THEN
            WRITE(6,611) ATIR,IILO,EELO,NFD
  611       FORMAT(F7.2,I6,F13.3,I8)
            NFQHT(IK)=NFD
            DO 150 IJ=1,NFD
              FRECQ(IK,IJ)=FRD(NFD-IJ+1)
              QHOT(IK,IJ)=CRD(NFD-IJ+1)*BAM
  150       CONTINUE
          ENDIF
  120   CONTINUE
  110 CONTINUE
  540 RETURN
C
  500 call quit(' ERROR IN DATA FILE FOR BF SIG OF AVERAGED LEVELS (1)')
  501 call quit(' ERROR IN DATA FILE FOR BF SIG OF AVERAGED LEVELS (2)')
C
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE PHTX(ID,ABSO,EMIS,fre,icon)
C     ======================================
C
C     Opacity due to detailed photoionization (read from tables by
C     routine SIGAVS)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      DIMENSION ABSO(MFREQ),EMIS(MFREQ),PLANF(MFREQ),STIMU(MFREQ)
      dimension fre(mfreq)
      DIMENSION PHOTI(MCROSS,MFREQ)
      DIMENSION IJP(MLEVEL),IJQ(MPHOT)
      PARAMETER (C3=1.4387886)
      SAVE PHOTI,IJP,IJQ
C
      IF(IASV.EQ.0 .AND. NQHT.EQ.0) RETURN
      T=TEMP(ID)
      nfre=nfreq
      ij0=3
      if(icon.eq.1) then
         ij0=1
         nfre=nfreqc
      end if
c
      DO 10 IJ=1,NFRE
         XX=FRE(IJ)
         X15=XX*1.E-15
         BNU=BN*X15*X15*X15
         HKF=HK*XX
         EXH=EXP(HKF/T)
         PLANF(IJ)=BNU/(EXH-1.)
         STIMU(IJ)=1.-1./EXH
   10 CONTINUE
C
      IF(IASV.EQ.0) GOTO 100
      IF(ID.EQ.1) THEN
        DO 40 I=1,NLEVEL
          IF(CRMX(I).EQ.0.) GOTO 40
          IK1=MAX0(2,IJP(I))
          DO 42 IJ=3,NFRE
            DO 45 IK=IK1,NFCR(I)
              IF(FRECR(I,IK).LT.FRE(IJ)) THEN
                IK2=IK
                GOTO 46
              ENDIF
   45       CONTINUE
   46       IK1=IK2
            IF(IJ.EQ.3) IJP(I)=IK1
            DFR=(FRE(IJ)-FRECR(I,IK1))/(FRECR(I,IK1-1)-FRECR(I,IK1))
            PHOTI(I,IJ)=CROSR(I,IK1)+DFR*(CROSR(I,IK1-1)-CROSR(I,IK1))
   42     CONTINUE
          PHOTI(I,1)=PHOTI(I,3)
          PHOTI(I,2)=PHOTI(I,NFREQ)
   40   CONTINUE
      ENDIF
      DO 30 I=1,NLEVEL
        IF(CRMX(I).EQ.0.) GOTO 30
        POP=POPUL(I,ID)
        DO 20 IJ=1,NFRE
          AB=PHOTI(I,IJ)*POP*STIMU(IJ)
          ABSO(IJ)=ABSO(IJ)+AB
          EMIS(IJ)=EMIS(IJ)+AB*PLANF(IJ)
   20   CONTINUE
   30 CONTINUE
C
  100 IF(NQHT.EQ.0) RETURN
      IF(ID.EQ.1) THEN
        DO 110 I=1,NQHT
          IF(CRMY(I).EQ.0.) GOTO 110
          IK1=MAX0(2,IJQ(I))
          DO 120 IJ=3,NFRE
            DO 125 IK=IK1,NFQHT(I)
              IF(FRECQ(I,IK).LT.FRE(IJ)) THEN
                IK2=IK
                GOTO 126
              ENDIF
  125       CONTINUE
  126       IK1=IK2
            IF(IJ.EQ.3) IJQ(I)=IK1
            DFR=(FRE(IJ)-FRECQ(I,IK1))/(FRECQ(I,IK1-1)-FRECQ(I,IK1))
            PHOTI(I,IJ)=QHOT(I,IK1)+DFR*(QHOT(I,IK1-1)-QHOT(I,IK1))
  120     CONTINUE
  110   CONTINUE
      ENDIF
      DO 210 I=1,NQHT
        IF(CRMY(I).EQ.0.) GOTO 210
        IAT=int(AQHT(I))
        X=(AQHT(I)-FLOAT(IAT)+1.E-4)*100.
        ION=INT(X)+1
        POP=RRR(ID,ION,IAT)*GQHT(I)*EXP(-EQHT(I)*C3/T)
        DO 220 IJ=3,NFRE
          AB=PHOTI(I,IJ)*POP*STIMU(IJ)
          ABSO(IJ)=ABSO(IJ)+AB
          EMIS(IJ)=EMIS(IJ)+AB*PLANF(IJ)
  220   CONTINUE
  210 CONTINUE
C
      RETURN
      END
C
C
C ********************************************************************
C
      subroutine getlal
c     =================
c
c     getlal reads in the profile functions for Lyman alpha, beta, gamma,
c     and Balmer alpha, including the quasi-molecular satellites;
c     valid for first and second order in neutral and ionized H density
c     modified routine provided originally by D. Koester
c           
c
      INCLUDE 'PARAMS.FOR'
      parameter (NXMAX=1400,NNMAX=5)
      common/quasun/nunalp,nunbet,nungam,nunbal
      common /callarda/xlalp(NXMAX),plalp(NXMAX,NNMAX),stnnea,stncha,
     *     vneua,vchaa,nxalp,iwarna
      common /callardb/xlbet(NXMAX),plbet(NXMAX,NNMAX),stnneb,stnchb,
     *     vneub,vchab,nxbet,iwarnb
      common /callardg/xlgam(NXMAX),plgam(NXMAX,NNMAX),stnneg,stnchg,
     *     vneug,vchag,nxgam,iwarng
      common /callardc/xlbal(NXMAX),plbal(NXMAX,NNMAX),stnnec,stnchc,
     *     vneuc,vchac,nxbal,iwarnc
c
c     Lyman alpha
c
      nxalp=0
      if(nunalp.gt.0) then
      nunalp=67
      open(unit=nunalp,file='./data/laquasi.dat',status='old')
      read(nunalp,*) nxalp,stnnea,stncha,vneua,vchaa
      do i=1,nxalp
         read(nunalp,*) xlalp(i),(plalp(i,j),j=1,NNMAX)
      end do
      close(nunalp)
      stnnea=10.0**stnnea
      stncha=10.0**stncha
      iwarna=0
      close(nunalp)
      write(*,*)
      write(*,*) ' read quasi-molecular data for L alpha'
      end if
c
c     Lyman beta
c
      nxbet=0
      if(nunbet.gt.0) then
      nunbet=67
      open(unit=nunbet,file='./data/lbquasi.dat',status='old')
      read(nunbet,*) nxbet,stnneb,stnchb,vneub,vchab
      do i=1,nxbet
         read(nunbet,*) xlbet(i),(plbet(i,j),j=1,NNMAX)
      end do
      close(nunbet)
      stnneb=10.0**stnneb
      stnchb=10.0**stnchb
      iwarnb=0
      write(*,*) ' read quasi-molecular data for L beta'
      end if
c
c     Lyman gamma
c
      nxgam=0
      if(nungam.gt.0) then
      nungam=67
      open(unit=nunalp,file='./data/lgquasi.dat',status='old')
      read(nungam,*) nxgam,stnneg,stnchg,vneug,vchag
      do i=1,nxgam
         read(nungam,*) xlgam(i),(plgam(i,j),j=1,NNMAX)
      end do
      close(nungam)
      stnneg=10.0**stnneg
      stnchg=10.0**stnchg
      iwarng=0
      write(*,*) ' read quasi-molecular data for L gamma'
      end if
c
c     Balmer alpha
c
      nxbal=0
      if(nunbal.gt.0) then
      nunbal=67
      open(unit=nunalp,file='./data/lhquasi.dat',status='old')
      read(nunbal,*) nxbal,stnnec,stnchc,vneuc,vchac
      do i=1,nxbal
         read(nunbal,*) xlbal(i),(plbal(i,j),j=1,NNMAX)
      end do
      close(nunbal)
      stnnec=10.0**stnnec
      stnchc=10.0**stnchc
      iwarnc=0
      write(*,*) ' read quasi-molecular data for H alpha'
      end if
      write(*,*) 
      return
      end
c
C
C ********************************************************************
C
      subroutine allard(xl,hneutr,hcharg,prof,iq,jq)
c     ==============================================
c
c     quasi-molecular opacity for Lyman alpha, beta, and Balmer alpha
c     modified routine provided originally by D. Koester
c
c     Input:  xl:  wavelength in [A]
c             hneutr:  neutral H particle density [cm-3]
c             hcharg: ionized H particle density [cm-3]
c             iq:   quantum number of the lower level
c             jq:   quantum number of the upper level;
c                   =2  -  Lyman alpha
c                   =3  -  Lyman beta
c     Output: prof:  Lyman alpha line profile, normalized to 1.0e8
c             if integrated over A;
c             It then renormalized by multiplying by
c             8.853e-29*lambda_0^2*f_ij
c
      INCLUDE 'PARAMS.FOR'
      parameter (NXMAX=1400,NNMAX=5)
      parameter (xnorma=8.8528e-29*1215.6*1215.6*0.41618,
     *           xnormb=8.8528e-29*1025.73*1025.7*0.0791,
     *           xnormg=8.8528e-29*972.53*972.53*0.0290,
     *           xnormc=8.8528e-29*6562.*6562.*0.6407)
      common /callarda/xlalp(NXMAX),plalp(NXMAX,NNMAX),stnnea,stncha,
     *     vneua,vchaa,nxalp,iwarna
      common /callardb/xlbet(NXMAX),plbet(NXMAX,NNMAX),stnneb,stnchb,
     *     vneub,vchab,nxbet,iwarnb
      common /callardg/xlgam(NXMAX),plgam(NXMAX,NNMAX),stnneg,stnchg,
     *     vneug,vchag,nxgam,iwarng
      common /callardc/xlbal(NXMAX),plbal(NXMAX,NNMAX),stnnec,stnchc,
     *     vneuc,vchac,nxbal,iwarnc
c
      prof=0.
c
c     Lyman alpha
c
      if(iq.eq.1.and.jq.eq.2) then    
c     if(xl.lt.xlalp(1).or.xl.gt.xlalp(nxalp)) return
      if(xl.lt.xlalp(1)) return
      vn1=hneutr/stnnea
      vn2=hcharg/stncha
      vns=vn1*vneua+vn2*vchaa
      if(iwarna.eq.0) then
         if(vn1*vneua.gt.0.3.or.vn2*vchaa.gt.0.3) then
            write(*,*) '          warning: density too high for',
     *           ' Lyman alpha expansion'
            iwarna=1
         endif
      endif
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      if(xl.le.xlalp(nxalp)) then  
      jl=0
      ju=nxalp+1
 10   if(ju-jl.gt.1) then       
         jm=(ju+jl)/2
         if((xlalp(nxalp).gt.xlalp(1)).eqv.(xl.gt.xlalp(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 10
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxalp) j=j-1
      a1=(xl-xlalp(j))/(xlalp(j+1)-xlalp(j))
      p1=  vn1*((1.0-a1)*plalp(j,1)+a1*plalp(j+1,1))
      p11=vn11*((1.0-a1)*plalp(j,2)+a1*plalp(j+1,2))
      p2=  vn2*((1.0-a1)*plalp(j,3)+a1*plalp(j+1,3))
      p22=vn22*((1.0-a1)*plalp(j,4)+a1*plalp(j+1,4))
      p12=vn12*((1.0-a1)*plalp(j,5)+a1*plalp(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnorma
c
      else
      j=nxalp-1
c     a1=(xl-xlalp(j))/(xlalp(j+1)-xlalp(j))
      a1=1.
      p1=  vn1*((1.0-a1)*plalp(j,1)+a1*plalp(j+1,1))
      p11=vn11*((1.0-a1)*plalp(j,2)+a1*plalp(j+1,2))
      p2=  vn2*((1.0-a1)*plalp(j,3)+a1*plalp(j+1,3))
      p22=vn22*((1.0-a1)*plalp(j,4)+a1*plalp(j+1,4))
      p12=vn12*((1.0-a1)*plalp(j,5)+a1*plalp(j+1,5))
      pro0=(p1+p2+p11+p22+p12)*xnorm*xnorma
      xlas=xlalp(nxalp)
      x0=1215.67
      dxlas=xlalp(nxalp)-x0
      dx=xl-x0
      prof=pro0/(dx/dxlas)**2.5
c 
      end if
      return
      end if
c
c     Lyman beta
c
      if(iq.eq.1.and.jq.eq.3) then
      if(nxbet.eq.0) return
      if(xl.lt.xlbet(1).or.xl.gt.xlbet(nxbet)) return
      vn1=hneutr/stnneb
      vn2=hcharg/stnchb
      vns=vn1*vneub+vn2*vchab
      if(iwarnb.eq.0) then
         if(vn1*vneub.gt.0.3.or.vn2*vchab.gt.0.3) then
            write(*,*) '          warning: density too high for',
     *           ' Lyman beta expansion'
            iwarnb=1
         endif
      endif
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxbet+1
 20   if(ju-jl.gt.1) then       
         jm=(ju+jl)/2
         if((xlbet(nxbet).gt.xlbet(1)).eqv.(xl.gt.xlbet(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 20
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxbet) j=j-1
      a1=(xl-xlbet(j))/(xlbet(j+1)-xlbet(j))
      p1=  vn1*((1.0-a1)*plbet(j,1)+a1*plbet(j+1,1))
      p11=vn11*((1.0-a1)*plbet(j,2)+a1*plbet(j+1,2))
      p2=  vn2*((1.0-a1)*plbet(j,3)+a1*plbet(j+1,3))
      p22=vn22*((1.0-a1)*plbet(j,4)+a1*plbet(j+1,4))
      p12=vn12*((1.0-a1)*plbet(j,5)+a1*plbet(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnormb
      return
      end if
c
c     Lyman gamma
c
      if(iq.eq.1.and.jq.eq.4) then
      if(nxgam.eq.0) return
      if(xl.lt.xlgam(1).or.xl.gt.xlgam(nxgam)) return
      vn1=hneutr/stnneg
      vn2=hcharg/stnchg
      vns=vn1*vneug+vn2*vchag
      if(iwarng.eq.0) then
         if(vn1*vneug.gt.0.3.or.vn2*vchag.gt.0.3) then
            write(*,*) '          warning: density too high for',
     *           ' Lyman gamma expansion'
            iwarng=1
         endif
      endif
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxgam+1
 30   if(ju-jl.gt.1) then       
         jm=(ju+jl)/2
         if((xlgam(nxgam).gt.xlgam(1)).eqv.(xl.gt.xlgam(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 30
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxgam) j=j-1
      a1=(xl-xlgam(j))/(xlgam(j+1)-xlgam(j))
      p1=  vn1*((1.0-a1)*plgam(j,1)+a1*plgam(j+1,1))
      p11=vn11*((1.0-a1)*plgam(j,2)+a1*plgam(j+1,2))
      p2=  vn2*((1.0-a1)*plgam(j,3)+a1*plgam(j+1,3))
      p22=vn22*((1.0-a1)*plgam(j,4)+a1*plgam(j+1,4))
      p12=vn12*((1.0-a1)*plgam(j,5)+a1*plgam(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnormg
      return
      end if
c
c     Balmer alpha
c
      if(iq.eq.2.and.jq.eq.3) then    
      if(xl.lt.xlbal(1).or.xl.gt.xlbal(nxbal)) return
c      vn1=hneutr/stnnec
      vn1=0.
      vn2=hcharg/stnchc
      vns=vn1*vneuc+vn2*vchac
      vn11=vn1*vn1
      vn22=vn2*vn2
      vn12=vn1*vn2
      xnorm=1.0/(1.0+vns+0.5*vns*vns)
c
      jl=0
      ju=nxbal+1
 40   if(ju-jl.gt.1) then       
         jm=(ju+jl)/2
         if((xlbal(nxbal).gt.xlbal(1)).eqv.(xl.gt.xlbal(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 40
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxbal) j=j-1
      a1=(xl-xlbal(j))/(xlbal(j+1)-xlbal(j))
      p1=  vn1*((1.0-a1)*plbal(j,1)+a1*plbal(j+1,1))
      p11=vn11*((1.0-a1)*plbal(j,2)+a1*plbal(j+1,2))
      p2=  vn2*((1.0-a1)*plbal(j,3)+a1*plbal(j+1,3))
      p22=vn22*((1.0-a1)*plbal(j,4)+a1*plbal(j+1,4))
      p12=vn12*((1.0-a1)*plbal(j,5)+a1*plbal(j+1,5))
      prof=(p1+p2+p11+p22+p12)*xnorm*xnormc
      end if
c     
      return
      end
C
C
C     **************************************************************
C
C
      subroutine lyahhe(xl,ahe,prof)
c     ==============================
c
c     Lyman alpha broadening by helium - after N. Allard
c
      INCLUDE 'PARAMS.FOR'
      parameter (nxmax=1000)
c     parameter (sthe=1.e21)
      common/hhebrd/sthe,nunhhe
      common/calhhe/xlhhe(nxmax),sighhe(nxmax),nxhhe
      dimension xlhh0(nxmax),sighh0(nxmax)
      data iread/0/
c
      if(iread.eq.0) then
c        nxhhe=679
c        open(unit=67,
c    *   file='siglyhhe_21_T14500.lam',
c    *   status='old')
         it=0
         do i=1,nxmax
            read(67,*,err=5,end=5) xl,sig
            it=it+1
            if(nunhhe.eq.1) xl=1./(1.e-8*xl+1./1215.67)
            xlhh0(it)=xl
            sighh0(it)=sig
         end do
    5    nxhhe=it
         do i=1,nxhhe
            xlhhe(i)=xlhh0(nxhhe-i+1)
            sighhe(i)=sighh0(nxhhe-i+1)
         end do
c        do i=1,nxhhe
c           j=nxhhe-i+1
c           read(67,*) xlhhe(j),sighhe(j)
c        end do
         close(67)
         iread=1
      end if
c
      prof=0.
      if(xl.gt.xlhhe(nxhhe)) return
      jl=0
      ju=nxhhe+1
 10   if(ju-jl.gt.1) then
         jm=(ju+jl)/2
         if((xlhhe(nxhhe).gt.xlhhe(1)).eqv.(xl.gt.xlhhe(jm))) then
            jl=jm
         else
            ju=jm
         endif
         go to 10
      endif
      j=jl
c
      if(j.eq.0) j=1
      if(j.eq.nxhhe) j=j-1
      a1=(xl-xlhhe(j))/(xlhhe(j+1)-xlhhe(j))
      s1=(1.0-a1)*sighhe(j)+a1*sighhe(j+1)
      prof=s1*ahe/sthe*6.2831855
      return
      end
C
C
C     **************************************************************
C
C
      subroutine readbf
c     =================
c
c     auxiliary subroutine for enabling reading of input data with
c     comments
c
c     lines beginning with ! or * are understood as comments
c
      INCLUDE 'PARAMS.FOR'
      character*80 buff
   10 continue
      read(5,501,end=20) buff
      if(buff(1:1).eq.'!'.or.buff(1:1).eq.'*') go to 10
      write(ibuff,501) buff
      go to 10
  501 format(a)
   20 continue
      rewind ibuff
      return 
      end
C
C
C     *******************************************************************
C
C

      SUBROUTINE PRETAB
C     =================
C
C     pretabulate expansion coefficients for the Voigt function
C     200 steps per doppler width - up to 10 Doppler widths
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (VSTEPS=200.,MVOI=2001)
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      DIMENSION TABVI(81),TABH1(81)
      DATA TABVI/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,
     11.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,
     2 3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,
     3 5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,
     4 9.0,9.2,9.4,9.6,9.8,10.0,10.2,10.4,10.6,10.8,11.0,11.2,11.4,11.6,
     5 11.8,12.0/
      DATA TABH1/-1.12838,-1.10596,-1.04048,-.93703,-.80346,-.64945,
     1-.48552,-.32192,-.16772,-.03012,.08594,.17789,.24537,.28981,
     2.31394,.32130,.31573,.30094,.28027,.25648,.231726,.207528,.184882,
     3.164341,.146128,.130236,.116515,.104739,.094653,.086005,.078565,
     4 .072129,.066526,.061615,.057281,.053430,.049988,.046894,.044098,
     5 .041561,.039250,.035195,.031762,.028824,.026288,.024081,.022146,
     6 .020441,.018929,.017582,.016375,.015291,.014312,.013426,.012620,
     7 .0118860,.0112145,.0105990,.0100332,.0095119,.0090306,.0085852,
     8 .0081722,.0077885,.0074314,.0070985,.0067875,.0064967,.0062243,
     9 .0059688,.0057287,.0055030,.0052903,.0050898,.0049006,.0047217,
     T .0045526,.0043924,.0042405,.0040964,.0039595/
C
      N=MVOI
      DO 10 I=1,N
   10    H0TAB(I)=FLOAT(I-1)/VSTEPS
      CALL INTERP(TABVI,TABH1,H0TAB,H1TAB,81,N,2,0,0)
      DO 20 I=1,N
         VV=(FLOAT(I-1)/VSTEPS)**2
         H0TAB(I)=EXP(-VV)
         H2TAB(I)=H0TAB(I)-(VV+VV)*H0TAB(I)
   20 CONTINUE
      RETURN
      END
C
C
C     *******************************************************************
C
C

      FUNCTION VOIGTK(A,V)
C     ====================
C
C     Voigt function after Kurucz (in Computational Astrophysics)
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (MVOI=2001)
      PARAMETER (ONE=1., THREE=3., TEN=10., FIFTN=15., TWOH=200.,
     *           C14142=1.4142, C11283=1.12838, C15=1.5,C32=3.2,
     *           C05642=0.5642,C79788=0.79788,C02=0.2,C14=1.4,
     *           C37613=0.37613,C23=2./3.,
     *           CV1=-.122727278,CV2=.532770573,CV3=-.96284325, 
     *           CV4=.979895032)
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      IV=int(V*TWOH+C15)
      IF(A.LT.C02) THEN
         IF(V.LE.TEN) THEN 
            VOIGTK=(H2TAB(IV)*A+H1TAB(IV))*A+H0TAB(IV)
          ELSE
            VOIGTK=C05642*A/(V*V)
         END IF
         RETURN
      END IF
      IF(A.GT.C14) GO TO 10
      IF(A+V.GT.C32) GO TO 10
      VV=V*V
      HH1=H1TAB(IV)+H0TAB(IV)*C11283
      HH2=H2TAB(IV)+HH1*C11283-H0TAB(IV)
      HH3=(ONE-H2TAB(IV))*C37613-HH1*C23*VV+HH2*C11283
      HH4=(THREE*HH3-HH1)*C37613+H0TAB(IV)*C23*VV*VV
      VOIGTK=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0TAB(IV))*
     *       (((CV1*A+CV2)*A+CV3)*A+CV4)
      RETURN
   10 AA=A*A
      VV=V*V
      U=(AA+VV)*C14142
      UU=U*U
      VOIGTK=((((AA-TEN*VV)*AA*THREE+FIFTN*VV*VV)/UU+THREE*VV-AA)/UU+
     *       ONE)*A*C79788/U
      RETURN
      END
C
C
C     *******************************************************************
C
C
 
      SUBROUTINE RTECD
C     ================
C
C     solution of the radiative transfer equation by Feautrier method
C     for two continuum points 
C     used when one employs RTEDFE, ie. the DFE method for the
C     transfer equation for the inner frequency points
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      DIMENSION D(3,3,MDEPTH),ANU(3,MDEPTH),AANU(MDEPTH),DDD(MDEPTH),
     *       AA(3,3),BB(3,3),CC(3,3),VL(3),AMU(3),WTMU(3),
     *       DT(MDEPTH),TAU(MDEPTH),
     *       RDD(MDEPTH),FKK(MDEPTH),ST0(MDEPTH),SS0(MDEPTH),
     *       RINT(MDEPTH,MMU)
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
      COMMON/CONSCA/SCC1(mdepth),SCC2(MDEPTH)
      PARAMETER (UN=1.D0, HALF=0.5D0)
      PARAMETER (THIRD=UN/3., QUART=UN/4., SIXTH=UN/6.D0)
      PARAMETER (TAUREF = 0.6666666666667)
      DATA AMU/.887298334620742D0,.5D0,.112701665379258D0/,
     1     WTMU/.277777777777778D0,.444444444444444D0,.277777777777778D0
     1         /
C
      NMU=3
      ND1=ND-1
C
C     loop over two continuum frequencies 
C
      DO 100 IJ=1,2       
      TAUMIN=CH(IJ,1)/DENS(1)*DM(1)*HALF
      TAU(1)=TAUMIN
      DO 10 I=1,ND1
         DT(I)=(DM(I+1)-DM(I))*(CH(IJ,I+1)/DENS(I+1)+CH(IJ,I)/DENS(I))*
     *         HALF
         ST0(I)=ET(IJ,I)/CH(IJ,I)
         SS0(I)=-SC(IJ,I)/CH(IJ,I)
         TAU(I+1)=TAU(I)+DT(I)
         IF(TAU(I).LE.TAUREF.AND.TAU(I+1).GT.TAUREF) IREF=I
   10 CONTINUE
      ST0(ND)=ET(IJ,ND)/CH(IJ,ND)
      SS0(ND)=-SC(IJ,ND)/CH(IJ,ND)
      FR=FREQ(IJ)
      BNU=BN*(FR*1.E-15)**3
      PLAND=BNU/(EXP(HK*FR/TEMP(ND  ))-UN)
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)
      DPLAN=(PLAND-DPLAN)/DT(ND1)
C
C   +++++++++++++++++++++++++++++++++++++++++
C   FIRST PART  -  VARIABLE EDDINGTON FACTORS
C   +++++++++++++++++++++++++++++++++++++++++
C
C   Allowance for wind blanketing
C
      ALB1=0.
c      IF(IWINDB.NE.0) ALB1=ALBEDO(IJ)/(UN+ALBEDO(IJ))/HALF
      DO 21 I=1,NMU
C
C   ************************
C   UPPER BOUNDARY CONDITION
C   ************************
C
         ID=1
         DTP1=DT(1)
         Q0=0.
         P0=0.
C
C        allowance for non-zero optical depth at the first depth point
C
         TAMM=TAUMIN/AMU(I)
         IF(TAMM.GT.0.01) THEN
            P0=UN-EXP(-TAMM)
          ELSE
            P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM)))
         END IF
         EX=UN-P0
         Q0=Q0+P0*AMU(I)*WTMU(I)
C
         DIV=DTP1/AMU(I)*THIRD
         VL(I)=DIV*(ST0(ID)+HALF*ST0(ID+1))+ST0(ID)*P0
         DO 20 J=1,NMU
            BB(I,J)=SS0(ID)*WTMU(J)*(DIV+P0)-ALB1*WTMU(J)
   20       CC(I,J)=-HALF*DIV*SS0(ID+1)*WTMU(J)
         BB(I,I)=BB(I,I)+AMU(I)/DTP1+UN+DIV
         CC(I,I)=CC(I,I)+AMU(I)/DTP1-HALF*DIV
         ANU(I,ID)=0.
   21 CONTINUE
C
C     Matrix inversion: instead of calling MATINV, a very fast inlined 
C     routine MINV3 for a specific 3 x 3 matrix inversion
C
C     CALL MATINV(BB,NMU,3)
C
C     ******************************
      BB(2,1)=BB(2,1)/BB(1,1)
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2)
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3)
      BB(3,1)=BB(3,1)/BB(1,1)
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2)
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3)
C
      BB(3,2)=-BB(3,2)
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1)
      BB(2,1)=-BB(2,1)
C
      BB(3,3)=UN/BB(3,3)
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2)
      BB(2,2)=UN/BB(2,2)
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1)
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1)
      BB(1,1)=UN/BB(1,1)
C
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1)
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2)
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1)
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2)
      BB(3,1)=BB(3,3)*BB(3,1)
      BB(3,2)=BB(3,3)*BB(3,2)
C     ******************************
C
      DO 22 I=1,NMU
         DO 22 J=1,NMU
            S=0.
            DO 23 K=1,NMU
   23          S=S+BB(I,K)*CC(K,J)
            D(I,J,ID)=S
            ANU(I,1)=ANU(I,1)+BB(I,J)*VL(J)
   22 CONTINUE
C
C   *******************
C   NORMAL DEPTH POINTS
C   *******************
C
      DO 34 ID=2,ND1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=HALF*(DTM1+DTP1)
         AL=UN/DTM1/DT0
         GA=UN/DTP1/DT0
         BE=AL+GA
         A=(UN-HALF*AL*DTP1*DTP1)*SIXTH
         C=(UN-HALF*GA*DTM1*DTM1)*SIXTH
         B=UN-A-C
         VL0=A*ST0(ID-1)+B*ST0(ID)+C*ST0(ID+1)
         DO 30 I=1,NMU
            DO 30 J=1,NMU
               AA(I,J)=-A*SS0(ID-1)*WTMU(J)
               CC(I,J)=-C*SS0(ID+1)*WTMU(J)
               BB(I,J)=B*SS0(ID)*WTMU(J)
   30    CONTINUE
         DO 31 I=1,NMU
            DIV=AMU(I)**2
            VL(I)=VL0
            AA(I,I)=AA(I,I)+DIV*AL-A
            CC(I,I)=CC(I,I)+DIV*GA-C
            BB(I,I)=BB(I,I)+DIV*BE+B
   31    CONTINUE
         DO 32 I=1,NMU
            S1=0.
            DO 36 J=1,NMU
               S=0.
               S1=S1+AA(I,J)*ANU(J,ID-1)
               DO 35 K=1,NMU
   35             S=S+AA(I,K)*D(K,J,ID-1)
   36          BB(I,J)=BB(I,J)-S
            VL(I)=VL(I)+S1
   32    CONTINUE
C
C     Matrix inversion: instead of calling MATINV, a very fast inlined 
C     routine MINV3 for a specific 3 x 3 matrix inversion
C
C     CALL MATINV(BB,NMU,3)
C
C     ******************************
      BB(2,1)=BB(2,1)/BB(1,1)
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2)
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3)
      BB(3,1)=BB(3,1)/BB(1,1)
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2)
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3)
C
      BB(3,2)=-BB(3,2)
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1)
      BB(2,1)=-BB(2,1)
C
      BB(3,3)=UN/BB(3,3)
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2)
      BB(2,2)=UN/BB(2,2)
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1)
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1)
      BB(1,1)=UN/BB(1,1)
C
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1)
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2)
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1)
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2)
      BB(3,1)=BB(3,3)*BB(3,1)
      BB(3,2)=BB(3,3)*BB(3,2)
C     ******************************
C
         DO 33 I=1,NMU
            ANU(I,ID)=0.
            DO 33 J=1,NMU
               S=0.
               DO 37 K=1,NMU
   37             S=S+BB(I,K)*CC(K,J)
               D(I,J,ID)=S
               ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
   33    CONTINUE
   34 CONTINUE
C
C   ************
C   LOWER BOUNDARY CONDITION
C   ************
C
      ID=ND
      DO 41 I=1,NMU
         AA(I,I)=AMU(I)/DTP1
         VL(I)=PLAND+AMU(I)*DPLAN+AA(I,I)*ANU(I,ID-1)
         DO 40 J=1,NMU
   40       BB(I,J)=-AA(I,I)*D(I,J,ID-1)
         BB(I,I)=BB(I,I)+AA(I,I)+UN
   41 CONTINUE
C
C     Matrix inversion: instead of calling MATINV, a very fast inlined 
C     routine MINV3 for a specific 3 x 3 matrix inversion
C
C     CALL MATINV(BB,NMU,3)
C
C     ******************************
      BB(2,1)=BB(2,1)/BB(1,1)
      BB(2,2)=BB(2,2)-BB(2,1)*BB(1,2)
      BB(2,3)=BB(2,3)-BB(2,1)*BB(1,3)
      BB(3,1)=BB(3,1)/BB(1,1)
      BB(3,2)=(BB(3,2)-BB(3,1)*BB(1,2))/BB(2,2)
      BB(3,3)=BB(3,3)-BB(3,1)*BB(1,3)-BB(3,2)*BB(2,3)
C
      BB(3,2)=-BB(3,2)
      BB(3,1)=-BB(3,1)-BB(3,2)*BB(2,1)
      BB(2,1)=-BB(2,1)
C
      BB(3,3)=UN/BB(3,3)
      BB(2,3)=-BB(2,3)*BB(3,3)/BB(2,2)
      BB(2,2)=UN/BB(2,2)
      BB(1,3)=-(BB(1,2)*BB(2,3)+BB(1,3)*BB(3,3))/BB(1,1)
      BB(1,2)=-BB(1,2)*BB(2,2)/BB(1,1)
      BB(1,1)=UN/BB(1,1)
C
      BB(1,1)=BB(1,1)+BB(1,2)*BB(2,1)+BB(1,3)*BB(3,1)
      BB(1,2)=BB(1,2)+BB(1,3)*BB(3,2)
      BB(2,1)=BB(2,2)*BB(2,1)+BB(2,3)*BB(3,1)
      BB(2,2)=BB(2,2)+BB(2,3)*BB(3,2)
      BB(3,1)=BB(3,3)*BB(3,1)
      BB(3,2)=BB(3,3)*BB(3,2)
C     ******************************
C
      DO 42 I=1,NMU
         ANU(I,ID)=0.
         DO 42 J=1,NMU
            D(I,J,ID)=0.
            ANU(I,ID)=ANU(I,ID)+BB(I,J)*VL(J)
   42 CONTINUE
C
C   ************
C   BACKSOLUTION
C   ************
C
      DO 54 ID=ND-1,1,-1
         DO 50 I=1,NMU
            DO 50 J=1,NMU
               ANU(I,ID)=ANU(I,ID)+D(I,J,ID)*ANU(J,ID+1)
   50    CONTINUE
         AJ=0.
         AK=0.
         DO 52 I=1,NMU
            DIV=WTMU(I)*ANU(I,ID)
            AJ=AJ+DIV
            AK=AK+DIV*AMU(I)**2
   52    CONTINUE
         FKK(ID)=AK/AJ
   54 CONTINUE
C
C     surface Eddington actor
C
      AH=0.
      DO 53 I=1,NMU
   53    AH=AH+WTMU(I)*AMU(I)*ANU(I,1)
      FH=AH/AJ-HALF*ALB1
C
      FKK(ND)=THIRD
C
C
C   +++++++++++++++++++++++++++++++++++++++++
C   SECOND PART  -  DETERMINATION OF THE MEAN INTENSITIES
C   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN EDDINGTON FACTORS
C   +++++++++++++++++++++++++++++++++++++++++
C
      DTP1=DT(1)
      DIV=DTP1*THIRD
      BBB=FKK(1)/DTP1+FH+DIV+SS0(1)*(DIV+Q0)
      CCC=FKK(2)/DTP1-HALF*DIV*(UN+SS0(2))
      VLL=DIV*(ST0(1)+HALF*ST0(2))+ST0(1)*Q0
      AANU(1)=VLL/BBB
      DDD(1)=CCC/BBB
      DO 60 ID=2,ND1
         DTM1=DTP1
         DTP1=DT(ID)
         DT0=HALF*(DTP1+DTM1)
         AL=UN/DTM1/DT0
         GA=UN/DTP1/DT0
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH
         AAA=AL*FKK(ID-1)-A*(UN+SS0(ID-1))
         CCC=GA*FKK(ID+1)-C*(UN+SS0(ID+1))
         BBB=(AL+GA)*FKK(ID)+(UN-A-C)*(UN+SS0(ID))
         VLL=A*ST0(ID-1)+C*ST0(ID+1)+(UN-A-C)*ST0(ID)
         BBB=BBB-AAA*DDD(ID-1)
         DDD(ID)=CCC/BBB
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB
   60 CONTINUE
      BBB=FKK(ND)/DTP1+HALF
      AAA=FKK(ND1)/DTP1
      BBB=BBB-AAA*DDD(ND1)
      VLL=HALF*PLAND+DPLAN*THIRD
      RDD(ND)=(VLL+AAA*AANU(ND1))/BBB
      DO 70 IID=1,ND1
         ID=ND-IID
         RDD(ID)=AANU(ID)+DDD(ID)*RDD(ID+1)
   70 CONTINUE
      FLUX(IJ)=FH*RDD(1)
C
      if(ij.eq.1) then
         do 72 id=1,nd
   72       scc1(id)=-rdd(id)*ss0(id)*ch(1,id)
       else
         do 74 id=1,nd
   74       scc2(id)=-rdd(id)*ss0(id)*ch(2,id)
      end if
C
C     if needed (if iprin.ge.3), output of interesting physical
C     quantities at the monochromatic optical depth  tau(nu)=2/3
C
      IF(IPRIN.ge.3) THEN
      T0=LOG(TAU(IREF+1)/TAU(IREF))
      X0=LOG(TAU(IREF+1)/TAUREF)/T0
      X1=LOG(TAUREF/TAU(IREF))/T0
      DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1)
      TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1)
      STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1)
      SCREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1)
      SSREF=EXP(LOG(-SS0(IREF)*RDD(IREF))*X0+
     *           LOG(-SS0(IREF+1)*RDD(IREF+1))*X1)
      SREF=STREF+SSREF
      ALM=2.997925E18/FREQ(IJ)
      WRITE(96,636) IJ,ALM,IREF,DMREF,TREF,SCREF,STREF,SSREF,SREF
  636 FORMAT(1H ,I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3,E11.3)
      END IF
c
C   ********************************************************************
C
C   THIRD PART  -  DETERMINATION OF THE SPECIFIC INTENSITIES
C   RECALCULATION OF THE TRANSFER EQUATION WITH GIVEN SOURCE FUNCTION
C
      if(iflux.eq.0) go to 100
      DO 85 IMU=1,NMU0
      ANX=ANGL(IMU)
      DTP1=DT(1)
      DIV=DTP1*THIRD/ANX
C
      TAMM=TAUMIN/ANX
      IF(TAMM.LT.0.01) THEN
         P0=TAMM*(UN-HALF*TAMM*(UN-TAMM*THIRD*(UN-QUART*TAMM)))
       ELSE
         P0=UN-EXP(-TAMM)
      END IF
C
      BBB=ANX/DTP1+UN+DIV
      CCC=ANX/DTP1-HALF*DIV
      VLL=(DIV+P0)*(ST0(1)-SS0(1)*RDD(1))
     *    +HALF*DIV*(ST0(2)-SS0(2)*RDD(2))
      AANU(1)=VLL/BBB
      DDD(1)=CCC/BBB
      DIV=ANX*ANX
      DO 82 ID=2,ND1
         DTM1=DT(ID-1)
         DTP1=DT(ID)
         DT0=HALF*(DTP1+DTM1)
         AL=UN/DTM1/DT0
         GA=UN/DTP1/DT0
         A=(UN-HALF*DTP1*DTP1*AL)*SIXTH
         C=(UN-HALF*DTM1*DTM1*GA)*SIXTH
         AAA=DIV*AL-A
         CCC=DIV*GA-C
         BBB=DIV*(AL+GA)+UN-A-C
         VLL=A*(ST0(ID-1)-SS0(ID-1)*RDD(ID-1))+
     *       C*(ST0(ID+1)-SS0(ID+1)*RDD(ID+1))+
     *       (UN-A-C)*(ST0(ID)-SS0(ID)*RDD(ID))
         BBB=BBB-AAA*DDD(ID-1)
         DDD(ID)=CCC/BBB
         AANU(ID)=(VLL+AAA*AANU(ID-1))/BBB
   82 CONTINUE
C
C     Lower boundary condition
C
      AAA=ANX/DTP1
      BBB=AAA+UN
      VLL=PLAND+ANX*DPLAN
C
      RINT(ND,IMU)=(VLL+AAA*AANU(ND1))/(BBB-AAA*DDD(ND1))
      DO 84 IID=1,ND1
         ID=ND-IID
         RINT(ID,IMU)=AANU(ID)+DDD(ID)*RINT(ID+1,IMU)
   84 CONTINUE
   85 CONTINUE
c
      FLX=0.
      DO 86 IMU=1,NMU0
         RINT(1,IMU)=RINT(1,IMU)/HALF
         FLX=FLX+ANGL(IMU)*WANGL(IMU)*RINT(1,IMU)
   86 CONTINUE
      FLX=FLX*HALF
c     FLUX(IJ)=FLX
C
C     output of emergent specific intensities in continuum to Unit 18 
C
      if(iflux.ge.1) then
      WRITE(18,641) WLAM(IJ),FLX,(RINT(1,IMU),IMU=1,NMU0)
      end if
  100 CONTINUE
  641 FORMAT(1H ,f10.3,1pe15.5/(1P5E15.5))
c  
c     call rtedfe for the internal points
c
      CALL RTEDFE
c
      RETURN
      END
C
C
C     *******************************************************************
C
C
 
      SUBROUTINE RTEDFE
C     =================
C
C     Solution of the radiative transfer equation - frequency by
C     frequency - for the known source function.
C
C     The numerical method used:
c     Discontinuous Finite Element (DFE) method
c     Castor, Dykema, Klein, 1992, ApJ 387, 561.
C
C     Input through blank COMMON block:
C      CH     - two-dimensional array  absorption coefficient (frequency,
C               depth)
C      ET     - emission coefficient (frequency, depth)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      PARAMETER (ONE=1.,TWO=2.,HALF=0.5)
      PARAMETER (TAUREF = 0.6666666666667)
      DIMENSION DT(MDEPTH),ST0(MDEPTH),AB0(MDEPTH),DELDM(MDEPTH),
     *          dtau(mdepth),rip(mdepth),rim(mdepth),riup(mdepth),
     *          AMU(3),WTMU(3),RINT1(MMU),
     *          AMUI(MMU),AMUW(MMU),TAU(MDEPTH),SS0(MDEPTH)
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
      COMMON/CONSCA/SCC1(mdepth),SCC2(MDEPTH)
      COMMON/REFDEP/IREFD(MFREQ)
C
C     angle points (AMU) and angular integration weights (WTMU)
C
      DATA AMU/.887298334620742D0,.5D0,.112701665379258D0/,
     * WTMU/.277777777777778D0,.444444444444444D0,.277777777777778D0/
C
      DO 10 I=1,ND-1
         DELDM(I)=HALF*(DM(I+1)-DM(I))
   10 CONTINUE
C
c     angle points
C
      IF(IFLUX.EQ.0) THEN
         NMUS=NMU
         do 12 i=1,nmu
            amui(i)=amu(i)
            amuw(i)=amu(i)*wtmu(i)
   12    continue
       ELSE IF(IFLUX.EQ.1) THEN
         NMUS=NMU0
         do 14 i=1,nmus
            amui(i)=angl(i)
            amuw(i)=angl(i)*wangl(i)
   14    continue
      END IF
C
C     overall loop over frequencies
C
      DO 500 IJ=1,NFREQ
      FR=FREQ(IJ)
C
C     total source function
C  
      DO 20 ID=1,ND
         AB0(ID)=CH(IJ,ID)
         SCT=FRX1(IJ)*SCC2(ID)+FRX2(IJ)*SCC1(ID)
         ST0(ID)=(ET(IJ,ID)+SCT)/AB0(ID)
         SS0(ID)=-SCT/AB0(ID)
   20 CONTINUE
      AH=0.
C
C     optical depth scale
C
      TAU(1)=0.
      IREF=1
      DO 30 ID=1,ND-1
         DT(ID)=DELDM(ID)*(AB0(ID+1)/DENS(ID+1)+AB0(ID)/DENS(ID))
         TAU(ID+1)=TAU(ID)+DT(ID)
         IF(TAU(ID).LE.TAUREF.AND.TAU(ID+1).GT.TAUREF) IREF=ID
   30 CONTINUE
      IREFD(IJ)=IREF
C
C     quantities for the lower boundary condition
C
      FR15=FR*1.D-15
      BNU=BN*FR15*FR15*FR15
      PLAND=BNU/(EXP(HK*FR/TEMP(ND))-ONE)
      DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-ONE)
      DPLAN=(PLAND-DPLAN)/DT(ND-1)
c
c     loop over angle poits
c
      DO 100 I=1,NMUS
         do id=1,nd-1
            dtau(id)=dt(id)/amui(i)
         enddo
C
c           outgoing intensity
c 
            rip(nd)=PLAND+AMUI(I)*DPLAN
            id=nd-1
            dt0=dtau(id)
            dtaup1=dt0+one
            dtau2=dt0*dt0
            bb=two*dtaup1
            cc=dt0*dtaup1
            aa=dtau2+bb
            rim(id+1)=(aa*rip(id+1)-cc*st0(id+1)+dt0*st0(id))/bb
            do id=nd-1,1,-1
               dt0=dtau(id)
               dtaup1=dt0+one
               dtau2=dt0*dt0
               bb=two*dtaup1
               cc=dt0*dtaup1
               aa=one/(dtau2+bb)
               rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa
               rip(id+1)=(bb*rim(id+1)+cc*st0(id+1)-dt0*st0(id))*aa
            enddo
            do id=2,nd-1
               riup(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))/
     *                  (dtau(id-1)+dtau(id))
            enddo
            riup(1)=rim(1)
            riup(nd)=rip(nd)
c
         AH=AH+AMUW(I)*RIUP(1)
         RINT1(I)=RIUP(1)
         rint1(i)=max(rint1(i),1.e-40)
c
c     end of the loop over angle points
c
  100 CONTINUE
      FLUX(IJ)=AH*HALF
      if(iflux.ge.1) then
C
C     output of emergent specific intensities to Unit 10 (line points)
C     or 18 (two continuum points)
C
      IF(IJ.GT.2) THEN
      WRITE(10,618) WLAM(IJ),FLUX(IJ),(RINT1(IMU),IMU=1,NMUS)
      ELSE
      WRITE(18,618) WLAM(IJ),FLUX(IJ),(RINT1(IMU),IMU=1,NMUS)
      END IF
      end if
  618 FORMAT(1H ,f10.3,1pe15.5/(1P5E15.5))
C
C     if needed (if iprin.ge.3), output of interesting physical
C     quantities at the monochromatic optical depth  tau(nu)=2/3
C
      IF(IPRIN.GE.3) THEN
      T0=LOG(TAU(IREF+1)/TAU(IREF))
      X0=LOG(TAU(IREF+1)/TAUREF)/T0
      X1=LOG(TAUREF/TAU(IREF))/T0
      DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1)
      TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1)
      STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1)
      SSREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1)
      SREF=STREF+SSREF
      ALM=2.997925E18/FREQ(IJ)
      WRITE(96,636) IJ,ALM,IREF,DMREF,TREF,STREF,SSREF,SREF
  636 FORMAT(1H ,I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3)
      END IF
c
C
C     end of the loop over frequencies
C
  500 CONTINUE
      RETURN
      END
C
C
C    *******************************************************************
C
C
      SUBROUTINE PARTF(IAT,IZI,T,ANE,XMAXN,U)
C     =======================================
C
C     Partition functions 
C     The standard evaluation is for hydrogen through zinc, for
C     neutrals and first four ionization degrees.
C     Basically after Traving, Baschek, and Holweger, Abhand. Hamburg.
C     Sternwarte. Band VIII, Nr. 1 (1966)
C
C     For higher atomic numbers  modified Kurucz routine PFSAHA,
C     called PFHEAV here is used. The routine was provided by
C     Charles Proffitt.
C     
C     The routine calls special procedures for Fe and Ni; or
C     the values based on the tabulated Opacity Project ionization
C     fractions
C
C     Input:
C      IAT   - atomic number
C      IZI   - ionic charge (=1 for neutrals, =2 for once ionized, etc)
C      T     - temperature
C      ANE   - electron density
C      XMAXN - principal quantum number of the last bound level
C
C     Output:
C      U     - partition function
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (NIONS=123, NSS=222)
      PARAMETER (UN=1.D0, HALF=0.5D0, TWO=2.D0, TRHA=1.5D0,
     *           THIRD=UN/3.D0, SIXTH=UN/6.D0)
      REAL*4 AHH( 6),  ALB(12),  AB (11),  AC (19),  AN (30),  AO (49),
     *       AF (34),  ANN(23),  ANA(19),  AMG(15),  AAL(17),  ASI(23),
     *       AP (19),  AS (29),  ACL(28),  AAR(25),  AK (30),  ACA(17),
     *       ASC(24),  ATI(33),  AV (33),  ACR(29),  AMN(28),  AFE(35),
     *       ACO(29),  ANI(23),  ACU(20),  AZN(18)
      REAL*4 GHH( 6),  GLB(12),  GB (11),  GC (19),  GN (30),  GO (49),
     *       GF (34),  GNN(23),  GNA(19),  GMG(15),  GAL(17),  GSI(23),
     *       GP (19),  GS (29),  GCL(28),  GAR(25),  GK (30),  GCA(17),
     *       GSC(24),  GTI(33),  GV (33),  GCR(29),  GMN(28),  GFE(35),
     *       GCO(29),  GNI(23),  GCU(20),  GZN(18)
      REAL*4 XL1(99), XL2(123),  XL(222),
     *       CH1(66),  CH2(72),  CH3(55),  CH4(29),  CHION(222)
      REAL*4 ALF(678), GAM(678)
      INTEGER   II1(5,15),II2(5,15),INDEX0(5,30),
     *          IS1(53),IS2(70),IS(123),INDEXS(123),
     *          IM1(99),IM2(123),IM(222),INDEXM(222),
     *          IGP1(99),IGP2(123),IGPR(222),
     *          IG01(53),IG02(70),IG0(123)
      DIMENSION IGLE(28)
C
      EQUIVALENCE   ( AHH(1), ALF(  1)),( ALB(1), ALF(  7)),
     *              ( AB (1), ALF( 19)),
     *              ( AC (1), ALF( 30)),( AN (1), ALF( 49)),
     *              ( AO (1), ALF( 79)),( AF (1), ALF(128)),
     *              ( ANN(1), ALF(162)),( ANA(1), ALF(185)),
     *              ( AMG(1), ALF(204)),( AAL(1), ALF(219)),
     *              ( ASI(1), ALF(236)),( AP (1), ALF(259)),
     *              ( AS (1), ALF(278)),( ACL(1), ALF(307)),
     *              ( AAR(1), ALF(335)),( AK (1), ALF(360)),
     *              ( ACA(1), ALF(390)),( ASC(1), ALF(407)),
     *              ( ATI(1), ALF(431)),( AV (1), ALF(464)),
     *              ( ACR(1), ALF(497)),( AMN(1), ALF(526)),
     *              ( AFE(1), ALF(554)),( ACO(1), ALF(589)),
     *              ( ANI(1), ALF(618)),( ACU(1), ALF(641)),
     *              ( AZN(1), ALF(661))
      EQUIVALENCE   ( GHH(1), GAM(  1)),( GLB(1), GAM(  7)),
     *              ( GB (1), GAM( 19)),
     *              ( GC (1), GAM( 30)),( GN (1), GAM( 49)),
     *              ( GO (1), GAM( 79)),( GF (1), GAM(128)),
     *              ( GNN(1), GAM(162)),( GNA(1), GAM(185)),
     *              ( GMG(1), GAM(204)),( GAL(1), GAM(219)),
     *              ( GSI(1), GAM(236)),( GP (1), GAM(259)),
     *              ( GS (1), GAM(278)),( GCL(1), GAM(307)),
     *              ( GAR(1), GAM(335)),( GK (1), GAM(360)),
     *              ( GCA(1), GAM(390)),( GSC(1), GAM(407)),
     *              ( GTI(1), GAM(431)),( GV (1), GAM(464)),
     *              ( GCR(1), GAM(497)),( GMN(1), GAM(526)),
     *              ( GFE(1), GAM(554)),( GCO(1), GAM(589)),
     *              ( GNI(1), GAM(618)),( GCU(1), GAM(641)),
     *              ( GZN(1), GAM(661))
      EQUIVALENCE   ( CH1(1), CHION(  1)),
     *              ( CH2(1), CHION( 67)),
     *              ( CH3(1), CHION(139)),
     *              ( CH4(1), CHION(194)),
     *              ( XL1(1),    XL(  1)),
     *              ( XL2(1),    XL(100))
      EQUIVALENCE   ( IS1(1),  IS(1)),   ( IS2(1),  IS( 54)),
     *              ( IM1(1),  IM(1)),   ( IM2(1),  IM(100)),
     *              (IGP1(1),IGPR(1)),   (IGP2(1),IGPR(100)),
     *              (IG01(1), IG0(1)),   (IG02(1), IG0( 54)),
     *              (II1(1,1),INDEX0(1,1)),(II2(1,1),INDEX0(1,16))
C
      DATA IGLE/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,
     *          10,21,28,25,6,25,28,21,10,21/
C
      DATA II1      /   1,  -1,   0,   0,   0,
     *                  2,   3,  -1,   0,   0,
     *                  4,   5,  -2,  -1,   0,
     *                  6,   7,  -1,  -2,  -1,
     *                  8,   9,  10,  -1,  -2,
     *                 11,  12,  13,  14,  -1,
     *                 15,  16,  17,  18,  19,
     *                 20,  21,  22,  23,  24,
     *                 25,  26,  27,  28,  -6,
     *                 29,  30,  31,  32,  -9,
     *                 33,  34,  35,  36,  -4,
     *                 37,  38,  39,  40,  -9,
     *                 41,  42,  43,  44,  -6,
     *                 45,  46,  47,  48,  -1,
     *                 49,  50,  51,  52,  53                         /
      DATA II2      /  54,  55,  56,  57,  58,
     *                 59,  60,  61,  62,  63,
     *                 64,  65,  66,  67,  68,
     *                 69,  70,  71,  72,  73,
     *                 74,  75,  76,  77,  -9,
     *                 78,  76,  80,  81,  82,
     *                 83,  84,  85,  86,  87,
     *                 88,  89,  90,  91,  92,
     *                 93,  94,  95,  96,  97,
     *                 98,  99, 100, 101, 102,
     *                103, 104, 105, 106, 107,
     *                108, 109, 110, 111, -25,
     *                112, 113, 114, 115,  -1,
     *                116, 117, 118, 119,  -1,
     *                120, 121, 122, 123,  -1                         /
C
      DATA IG01     /   2,
     *                  1,   2,
     *                  2,   1,
     *                  1,   2,
     *                  2,   1,   2,
     *                  1,   2,   1,   2,
     *                  4,   1,   2,   1,   2,
     *                  5,   4,   1,   2,   1,
     *                  4,   5,   4,   1,
     *                  1,   4,   5,   4,
     *                  2,   1,   4,   5,
     *                  1,   2,   1,   4,
     *                  2,   1,   2,   1,
     *                  1,   2,   1,   2,
     *                  4,   1,   2,   1,   2                         /
      DATA  IG02    /   5,   4,   1,   2,   1,
     *                  4,   5,   4,   1,   2,
     *                  1,   4,   5,   4,   1,
     *                  2,   1,   4,   5,   4,
     *                  1,   2,   1,   4,
     *                  4,   3,   4,   1,   4,
     *                  5,   4,   5,   4,   1,
     *                  4,   1,   4,   5,   4,
     *                  7,   6,   1,   4,   5,
     *                  6,   7,   6,   1,   4,
     *                  9,  10,   9,   6,   1,
     *                 10,   9,  10,  20,
     *                  9,   6,   9,  28,
     *                  2,   1,   6,  21,
     *                  1,   2,   1,  10                              /
C
      DATA IS1      /   1,
     *                  1,   1,
     *                  1,   1,
     *                  2,   1,
     *                  1,   2,   1,
     *                  1,   2,   2,   1,
     *                  2,   2,   3,   2,   1,
     *                  3,   4,   3,   5,   2,
     *                  2,   3,   4,   3,
     *                  2,   2,   3,   2,
     *                  1,   2,   2,   3,
     *                  1,   1,   2,   2,
     *                  2,   2,   1,   2,
     *                  1,   2,   2,   1,
     *                  2,   1,   1,   1,   1                         /
      DATA  IS2     /   3,   2,   1,   2,   2,
     *                  2,   3,   2,   1,   1,
     *                  2,   2,   3,   1,   1,
     *                  1,   2,   3,   3,   2,
     *                  2,   1,   2,   2,
     *                  3,   1,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  2,   3,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  3,   2,   1,   1,   1,
     *                  3,   2,   2,   1,   1,
     *                  4,   2,   1,   1,
     *                  2,   2,   1,   1,
     *                  3,   2,   1,   1,
     *                  3,   3,   1,   1                              /
C
      DATA IM1      /   2,
     *                  2,   2,
     *                  2,   2,
     *                  3,   2,   3,
     *                  3,   3,   2,   3,
     *                  4,   3,   3,   3,   3,   3,
     *                  3,   3,   4,   3,   3,   4,   2,   3,   2,   3,
     *                  4,   2,   2,   4,   2,   3,   3,   4,   4,   2,
     *                  3,   4,   2,   2,   2,   3,   3,
     *                  3,   3,   4,   2,   2,
     *                  4,   2,   3,   2,   5,   2,   2,
     *                  2,   2,   3,   2,   4,   2,   2,   4,   2,
     *                  2,   2,   2,   3,   2,   4,   2,   2,
     *                  3,   3,   2,   2,   3,   2,
     *                  3,   2,   3,   2,   3,   2,   2,
     *                  5,   4,   4,   4,   3,   3,
     *                  3,   2,   4,   4,   3,   3                    /
      DATA  IM2     /   4,   2,   2,   4,   2,   5,   4,   2,   3,   1,
     *                  3,   2,   5,   2,   2,   4,   2,   4,   4,
     *                  2,   2,   3,   2,   4,   2,   2,   4,   4,
     *                  3,   2,   3,   3,   2,   3,
     *                  4,   2,   2,   4,   2,
     *                  3,   2,   3,   2,   2,   3,   2,
     *                  4,   3,   3,   5,   4,   2,   3,
     *                  6,   4,   3,   6,   3,   5,   4,   2,
     *                  5,   3,   5,   4,   4,   4,   4,   4,
     *                  3,   3,   3,   4,   4,   4,   4,   4,
     *                  3,   2,   3,   4,   4,   4,   4,   4,
     *                  4,   4,   3,   5,   3,   4,   4,   4,   4,
     *                  5,   3,   3,   3,   5,   4,   5,   1,
     *                  6,   3,   5,   3,   5,   1,
     *                  2,   3,   3,   4,   3,   4,   1,
     *                  2,   2,   2,   3,   3,   2,   3,   1          /
C
      DATA IGP1     /   2,
     *                  4,   2,
     *                  2,   4,
     *                  4,  12,   2,
     *                  2,   4,  12,   2,
     *                 12,   2,  18,   4,  12,   2,
     *                 18,  10,  12,  24,   2,  18,   6,   4,  12,   2,
     *                  8,  20,  12,  18,  10,   2,  10,  12,  24,  20,
     *                  2,  18,   6,  18,  10,   4,  12,
     *                 18,  10,   8,  20,  12,
     *                 18,  10,   2,  10,  12,  24,  20,
     *                  8,   4,  18,  10,   8,  20,  12,  18,  10,
     *                  2,   8,   4,  18,  10,   8,  20,  12,
     *                  4,   2,   8,   4,  18,  10,
     *                  2,  18,   4,  12,   2,   8,   4,
     *                 12,   2,  18,   4,  12,   2,
     *                 18,  10,  12,   2,   4,   2                    /
      DATA  IGP2    /   8,  20,  12,  18,  10,  12,   2,  18,   4,  12,
     *                 18,  10,   8,  20,  12,  18,  10,  12,   2,
     *                  8,   4,  18,  10,   8,  20,  12,  18,  12,
     *                  2,   8,   4,  18,  10,   2,
     *                  8,  20,  12,  18,  10,
     *                  4,  20,   2,   8,   4,  18,  10,
     *                 30,  42,  18,  20,   2,  12,  18,
     *                 56,  56,  28,  42,  10,  20,   2,  12,
     *                 50,  70,  56,  72,  64,  42,  20,   2,
     *                 12,  60,  40,  50,  18,  56,  42,  20,
     *                 14,  10,  50,  12,  72,  50,  56,  42,
     *                 60,  56,  40,  50,  18,  12,  72,  50,  56,
     *                 42,  70,  42,  18,  56,  24,  50,  12,
     *                 20,  56,  42,  18,  56,  50,
     *                  2,  30,  10,  20,  56,  42,  56,
     *                  4,   8,  12,   2,  30,  10,  20,  42          /
C
      DATA XL1      /11.0,
     *                8.0,12.0,
     *                6.0, 6.0,
     *                6.0, 4.0, 8.0,
     *                9.0, 6.0, 4.0, 6.0,
     *                6.0, 6.0, 5.0, 6.1, 5.0, 6.0,
     *                6.1, 4.0, 5.0, 3.9, 6.0, 5.0, 4.0, 6.0, 6.3, 6.0,
     *                8.0, 6.0, 3.4, 6.0, 5.0, 3.9, 3.9, 6.0, 4.9, 4.0,
     *                5.9, 5.0, 4.9, 4.0, 4.0, 6.0, 6.0,
     *                4.0, 4.0, 5.0, 4.0, 4.0,
     *                5.0, 4.0, 3.9, 4.0, 5.0, 5.0, 4.0,
     *                6.0, 6.0, 5.0, 4.0, 3.9, 4.0, 4.0, 5.0, 5.0,
     *                7.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0,
     *                7.0, 7.0, 5.0, 5.0, 5.0, 5.0,
     *                7.0, 4.0, 7.0, 4.0, 7.0, 5.0, 5.0,
     *                6.1, 5.9, 5.0, 5.0, 5.0, 7.0,
     *                5.0, 5.0, 5.0, 7.0, 8.6, 8.0                    /
      DATA  XL2     / 6.0, 5.0, 5.0, 5.0, 5.0, 3.5, 5.0,14.4, 5.0, 4.0,
     *                6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,
     *                6.0, 6.0, 5.1, 5.0, 5.0, 5.0, 5.0, 5.0, 4.0,
     *                7.0, 5.0, 5.0, 6.0, 6.0, 5.0,
     *                6.0, 5.0, 5.0, 3.6, 4.0,
     *                5.9, 6.0, 7.0, 5.0, 4.9, 5.0, 4.3,
     *                4.9, 4.9, 5.0, 5.0, 6.0, 4.6, 3.8,
     *                5.0, 4.7, 5.0, 5.0, 5.0, 5.0, 6.0, 4.8,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,11.2,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.2,
     *                6.0, 5.0, 6.0, 7.0, 5.0, 5.0, 5.0, 5.0,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 5.0, 3.6, 3.8,
     *                5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.0,
     *                5.4, 5.0, 9.0, 5.0, 5.0, 3.0,
     *                8.0, 6.0, 5.0, 7.0, 5.0, 5.0, 2.9,
     *                8.0, 5.0, 5.0, 8.0, 5.0, 5.0, 5.0, 2.8          /
C
C
      DATA CH1      /  13.595 ,
     *                 24.580 ,  54.403 ,
     *                  5.390 , 75.619 ,
     *                  9.320 ,  13.278 ,  18.206 ,
     *                  8.296 ,  25.149 ,  31.146 ,  37.920 ,
     *                 11.256 ,  24.376 ,  30.868 ,  47.871 ,  55.873 ,
     *                 64.476 ,
     *                 14.529 ,  16.428 ,  29.593 ,  36.693 ,
     *                 47.426 ,  55.765 ,  63.626 ,  77.450 ,  87.445 ,
     *                 97.863 ,
     *                 13.614 ,  16.938 ,  18.630 ,
     *                 35.108 ,  37.621 ,  40.461 ,  42.584 ,
     *                 54.886 ,  63.733 ,  70.556 ,
     *                 77.394 ,  87.609 ,  97.077 , 103.911 , 106.116 ,
     *                113.873 , 125.863 ,
     *                 17.418 ,  20.009 ,  34.977 ,  39.204 ,  41.368 ,
     *                 62.646 ,  65.774 ,  69.282 ,  71.882 ,
     *                 87.139 ,  97.852 , 106.089 ,
     *                 21.559 ,  21.656 ,  41.071 ,  44.274 ,
     *                 63.729 ,  68.806 ,  71.434 ,  97.162 , 100.917 /
      DATA CH2      /   5.138 ,  47.290 ,  47.459 ,  71.647 ,  75.504 ,
     *                 98.880 , 104.778 , 107.864 ,
     *                  7.644 ,  15.031 ,  80.117 ,  80.393 ,
     *                109.294 , 113.799 ,
     *                  5.984 ,  10.634 ,  18.823 ,  25.496 ,
     *                 28.441 , 119.957 , 120.383 ,
     *                  8.149 ,  16.339 ,  22.894 ,
     *                 33.459 ,  42.333 ,  45.130 ,
     *                 10.474 ,  11.585 ,  19.720 ,
     *                 30.156 ,  51.354 ,  65.007 ,
     *                 10.357 ,  12.200 ,  13.401 ,  23.405 ,  24.807 ,
     *                 35.047 ,  47.292 ,  57.681 ,  72.474 ,  85.701 ,
     *                 13.014 ,  14.458 ,  23.798 ,  26.041 ,  27.501 ,
     *                 39.904 ,  41.610 ,  53.450 ,  67.801 ,
     *                 15.755 ,  15.933 ,  27.619 ,  29.355 ,
     *                 40.899 ,  42.407 ,  45.234 ,  59.793 ,  75.002 ,
     *                  4.339 ,  31.810 ,  32.079 ,
     *                 45.738 ,  47.768 ,  50.515 ,
     *                 60.897 ,  63.890 ,  65.849 ,  82.799 ,  85.150 /
      DATA CH3      /   6.111 ,   7.808 ,  11.868 ,
     *                 51.207 ,  51.596 ,  67.181 ,  69.536 ,
     *                  6.538 ,   7.147 ,   8.042 ,
     *                 12.891 ,  24.752 ,  74.090 ,  91.847 ,
     *                  6.818 ,  6.953 ,  7.411 ,
     *                 13.635 ,  14.685 ,  28.137 ,  43.236 , 100.083 ,
     *                  6.738 ,   7.101 ,  14.205 ,  15.670 ,  16.277 ,
     *                 29.748 ,  48.464 ,  65.198 ,
     *                  6.763 ,   8.285 ,   9.221 ,
     *                 16.493 ,  18.662 ,  30.950 ,  49.580 ,  73.093 ,
     *                  7.432 ,   8.606 ,   9.240 ,  15.636 ,  18.963 ,
     *                 33.690 ,  53.001 ,  76.006 ,
     *                  7.896 ,   8.195 ,   8.927 ,  16.178 ,  18.662 ,
     *                 30.640 ,  34.607 ,  56.001 ,  79.001           /
      DATA CH4      /   7.863 ,   8.378 ,   9.160 ,   9.519 ,
     *                 17.052 ,  18.958 ,  33.491 ,  53.001 ,
     *                  7.633 ,   8.793 ,  18.147 ,  20.233 ,  35.161 ,
     *                 56.025 ,
     *                  7.724 ,  10.532 ,  10.980 ,
     *                 20.286 ,  27.985 ,  36.826 ,  61.975 ,
     *                  9.391 ,  17.503 ,  17.166 ,
     *                 17.959 ,  27.757 ,  28.310 ,  39.701 ,  65.074 /
C
      DATA AHH      /  20.4976, 747.5023,
     *                 28.1703, 527.8296,  22.2809, 987.7189          /
      DATA GHH      /  10.853 ,  13.342 ,
     *                 21.170 ,  24.125 ,  43.708 ,  53.542           /
C
      DATA ALB      /   8.4915,  97.5015,  23.3299, 192.6701,
     *                  9.1849,  32.9263, 183.8887,  19.9563,  88.0437,
     *                  6.0478,  35.9723, 233.9798                    /
      DATA GLB      /   2.022 ,   4.604 ,  62.032 ,  72.624 ,
     *                  2.735 ,   6.774 ,   8.569 ,  10.750 ,  11.672 ,
     *                  3.967 ,  12.758 ,  16.692                     /
C
      DATA AB       /   4.0086,  19.6741, 402.3110,
     *                  9.7257,  30.9262, 186.3466,  44.1629,  60.8371,
     *                  6.0084,  23.5767,  76.4149                    /
      DATA GB       /   0.002 ,   3.971 ,   7.882 ,
     *                  4.720 ,  13.477 ,  22.103 ,  23.056 ,  24.734 ,
     *                  6.000 ,  24.540 ,  32.300                     /
C
      DATA AC       /   8.0158,   5.8833,  33.7521, 595.3432,
     *                  4.0003,  17.0841,  82.9154,
     *                 15.9808,  48.2044, 435.8093,
     *                 10.0281,  15.7574, 186.2109,
     *                 15.4127,  55.9559, 243.6311,
     *                  6.0057,  23.5757,  76.4185                    /
      DATA GC       /   0.004 ,   1.359 ,   6.454 ,  10.376 ,
     *                  0.008 ,  16.546 ,  21.614 ,
     *                  5.688 ,  15.801 ,  26.269 ,
     *                  6.691 ,  25.034 ,  40.975 ,
     *                 17.604 ,  36.180 ,  47.133 ,
     *                  8.005 ,  40.804 ,  54.492                     /
C
      DATA AN       /  14.0499,  30.8008, 883.1443,
     *                 10.0000,  16.0000,  64.0000,
     *                  8.0462,   6.2669,  17.8696, 282.8084,
     *                  7.3751,  33.1390, 215.4829,
     *                  4.0003,  19.3533,  80.6462,
     *                 13.0998,  19.6425,  94.3035, 370.9539,
     *                 16.0000,  38.0000,
     *                 10.3289,  14.5021, 187.1624, 108.1615, 191.8383,
     *                  6.0044,  23.5612,  76.4344                    /
      DATA GN       /   2.554 ,   9.169 ,  13.651 ,
     *                 12.353 ,  13.784 ,  14.874 ,
     *                  0.014 ,   2.131 ,  15.745 ,  24.949 ,
     *                  6.376 ,  14.246 ,  29.465 ,
     *                  0.022 ,  31.259 ,  41.428 ,
     *                  7.212 ,  15.228 ,  34.387 ,  46.708 ,
     *                 46.475 ,  49.468 ,
     *                  8.693 ,  37.650 ,  65.479 ,  61.155 ,  79.196 ,
     *                  9.999 ,  60.991 ,  82.262                     /
C
      DATA AO       /   4.0029,   5.3656,  36.2853,1044.3447,
     *                131.0217, 868.9779,  14.8533,  93.1466,
     *                 12.7843,   5.6828,  98.0919, 829.4396,
     *                 50.9878, 199.0120,   2.0000,   6.0000,  10.0000,
     *                 10.0000,  30.0000,  50.0000,
     *                  8.0703,   5.7144,  84.1156, 529.0927,
     *                  5.6609,  28.9355, 111.3620, 494.0413,
     *                 45.5249, 134.4751,
     *                  4.0003,  21.2937,  78.7058,
     *                 12.8293,  16.2730, 123.6578, 327.2396,
     *                 48.7883, 102.2117,  20.0060, 161.9903,
     *                 28.4184,  61.5816,
     *                 10.5563,  13.2950, 188.1390,
     *                 14.6560, 129.4922, 470.8512                    /
      DATA GO       /   0.022 ,   2.019 ,   9.812 ,  13.087 ,
     *                 13.804 ,  16.061 ,  14.293 ,  16.114 ,
     *                  3.472 ,   7.437 ,  22.579 ,  32.035 ,
     *                 27.774 ,  33.678 ,  28.118 ,  31.019 ,  34.204 ,
     *                 30.892 ,  33.189 ,  36.181 ,
     *                  0.032 ,   2.760 ,  35.328 ,  48.277 ,
     *                  7.662 ,  16.786 ,  42.657 ,  54.522 ,
     *                 50.204 ,  56.044 ,
     *                  0.048 ,  50.089 ,  66.604 ,
     *                  8.954 ,  18.031 ,  57.755 ,  72.594 ,
     *                 68.388 ,  82.397 ,  31.960 ,  76.876 ,
     *                 75.686 ,  80.388 ,
     *                 10.747 ,  52.323 ,  94.976 ,
     *                 27.405 ,  86.350 , 109.917                     /
C
      DATA AF       /   2.0001,  39.9012, 122.0986,
     *                 10.0000,  30.0000,  50.0000,
     *                  4.0199,   5.5741,  22.1839, 190.2179,
     *                 53.0383, 126.9616,  31.6894,  75.3105,
     *                 13.5014,   7.9936,  55.7981, 298.7039,
     *                 26.2496,  63.7503,   2.0000,   6.0000,  10.0000,
     *                 28.7150,  71.2850,
     *                  8.0153,   6.1931,  21.7287,  48.7780, 278.2782,
     *                178.5560, 421.4435,  51.7632,  95.2368          /
      DATA GF       /   0.050 ,  13.317 ,  15.692 ,
     *                 15.361 ,  17.128 ,  18.498 ,
     *                  0.048 ,   2.735 ,  20.079 ,  30.277 ,
     *                 27.548 ,  32.532 ,  30.391 ,  34.707 ,
     *                  4.479 ,  12.072 ,  31.662 ,  51.432 ,
     *                 44.283 ,  50.964 ,  46.193 ,  50.436 ,  54.880 ,
     *                 50.816 ,  57.479 ,
     *                  0.058 ,   3.434 ,  14.892 ,  37.472 ,  69.883 ,
     *                 67.810 ,  83.105 ,  72.435 ,  79.747           /
C
      DATA ANN      /  34.5080, 365.4919,  16.5768, 183.4231,
     *                  2.0007,  89.5607, 380.4381,  26.4473,  63.5527,
     *                  4.0342,   5.6162,  11.5176,  72.8273,
     *                 48.5684, 131.4315,  31.1710,  76.8290,
     *                 14.0482,  13.3077,  52.7897, 467.8487,
     *                 54.2196, 195.7800                              /
      DATA GNN      /  17.796 ,  20.730 ,  17.879 ,  20.855 ,
     *                  0.097 ,  29.878 ,  37.221 ,  31.913 ,  37.551 ,
     *                  0.092 ,   3.424 ,  24.806 ,  46.616 ,
     *                 45.643 ,  54.147 ,  48.359 ,  57.420 ,
     *                  5.453 ,  18.560 ,  46.583 ,  80.101 ,
     *                 70.337 ,  85.789                               /
C
      DATA ANA      /  11.6348, 158.3593,
     *                 21.0453,  50.9546,  10.1389,  25.8611,
     *                  2.0019,  38.0569, 137.9398,  28.3106,  61.6893,
     *                  4.0334,   5.8560,  18.1786, 208.9142,
     *                 93.6895, 406.3095,  60.4276, 239.5719          /
      DATA GNA      /   2.400 ,   4.552 ,
     *                 34.367 ,  40.566 ,  34.676 ,  40.764 ,
     *                  0.170 ,  44.554 ,  57.142 ,  51.689 ,  60.576 ,
     *                  0.152 ,   4.260 ,  36.635 ,  83.254 ,
     *                 72.561 ,  89.475 ,  75.839 ,  92.582           /
C
      DATA AMG      /  10.7445, 291.5057,  53.7488,
     *                  6.2270,  31.1291, 132.6438,
     *                 40.4379, 159.5618,  20.3845,  79.6154,
     *                  2.0007, 106.8977, 343.1010,  10.1326, 237.8581/
      DATA GMG      /   2.805 ,   6.777 ,   9.254 ,
     *                  4.459 ,   9.789 ,  13.137 ,
     *                 57.413 ,  71.252 ,  58.010 ,  71.660 ,
     *                  0.276 ,  74.440 ,  94.447 ,  54.472 ,  95.858 /
C
      DATA AAL      /   4.0009,  11.7804, 142.2179,  13.6585,  96.3371,
     *                 10.0807,  49.5843, 285.3343,  14.6872,  59.3122,
     *                  6.3277,  29.5086, 134.1634,
     *                 46.3164, 153.6833,  22.9896,  77.0103          /
      DATA GAL      /   0.014 ,   3.841 ,   5.420 ,   3.727 ,   8.833 ,
     *                  4.749 ,  11.902 ,  16.719 ,  11.310 ,  18.268 ,
     *                  6.751 ,  16.681 ,  24.151 ,
     *                 83.551 , 104.787 ,  84.293 , 105.171           /
C
      DATA ASI      /   7.9658,   4.6762,   1.3512, 123.2267, 443.7797,
     *                  4.0000,   7.4186,  24.1754, 60.4060,
     *                 14.4695,  11.9721,  26.5062, 269.0521,
     *                  9.1793,   4.8766,  29.1442,  52.7998,
     *                 13.2674,  36.0417, 180.6910,
     *                  6.4839,  27.6851, 135.8301                    /
      DATA GSI      /   0.020 ,   0.752 ,   1.614 ,   5.831 ,   7.431 ,
     *                  0.036 ,   8.795 ,  11.208 ,  13.835 ,
     *                  5.418 ,   7.825 ,  14.440 ,  19.412 ,
     *                  6.572 ,  11.449 ,  18.424 ,  25.457 ,
     *                 15.682 ,  27.010 ,  34.599 ,
     *                  9.042 ,  24.101 ,  37.445                     /
C
      DATA AP       /  13.5211,  22.2130, 353.2583,  10.0000, 150.0000,
     *                  8.0241,   5.8085,  51.7542, 252.4002,
     *                  4.0021,  20.7985,  62.4194, 200.7786,
     *                 11.7414,  63.5124, 179.7420,
     *                  6.8835,  32.7777, 228.3366                    /
      DATA GP       /   1.514 ,   5.575 ,   9.247 ,   8.076 ,  10.735 ,
     *                  0.043 ,   1.212 ,   8.545 ,  15.525 ,
     *                  0.074 ,   7.674 ,  16.639 ,  25.118 ,
     *                  8.992 ,  24.473 ,  40.704 ,
     *                 11.464 ,  33.732 ,  55.455                     /
C
      DATA AS       /   3.9615,   5.0780,  15.0944, 362.8588,
     *                 51.5995, 268.4002,  12.0000, 276.0000,
     *                 11.4377,   5.5126, 141.0009, 254.0478,
     *                 33.0518, 126.9479,
     *                  4.0707,   4.0637,   5.7245, 144.6376, 106.4909,
     *                  4.0011,  19.2813,  27.5990,  35.1179,
     *                 94.7454, 283.2486,
     *                 10.5474,  28.7137,  65.7378,  24.0000          /
      DATA GS       /   0.053 ,   1.121 ,   5.812 ,   9.425 ,
     *                  8.936 ,  11.277 ,   9.600 ,  12.551 ,
     *                  1.892 ,   3.646 ,  13.550 ,  19.376 ,
     *                 16.253 ,  21.062 ,
     *                  0.043 ,   0.123 ,   1.590 ,  13.712 ,  22.050 ,
     *                  0.118 ,   9.545 ,  18.179 ,  31.441 ,
     *                 30.664 ,  56.150 ,
     *                 10.704 ,  27.075 ,  50.599 ,  43.034           /
C
      DATA ACL      /   2.0007,  62.5048, 669.4942,  29.0259, 130.9740,
     *                  3.9064,   0.3993,   5.3570,  60.3424, 119.9913,
     *                138.1567, 278.8418, 102.3681, 158.6314,
     *                 12.6089,   5.9527, 110.5635, 262.8715,
     *                 69.2035, 100.7960,
     *                  7.3458,   5.6638,  44.1256, 202.7846,
     *                  4.0037,  21.8663,  40.5363,  57.5919          /
      DATA GCL      /   0.110 ,   9.919 ,  12.280 ,  11.017 ,  13.532 ,
     *                  0.092 ,   0.581 ,   1.620 ,  13.121 ,  19.787 ,
     *                 16.365 ,  21.988 ,  18.065 ,  23.594 ,
     *                  2.358 ,   5.708 ,  19.084 ,  30.683 ,
     *                 24.880 ,  33.229 ,
     *                  0.102 ,   1.391 ,  14.709 ,  36.968 ,
     *                  0.185 ,  11.783 ,  25.653 ,  44.698           /
C
      DATA AAR      /  43.6623, 324.3375,  20.8298, 163.1701,
     *                  2.0026, 137.4515, 258.5445,  62.8129, 149.1867,
     *                  4.0495,  14.4466,  46.8234, 124.6651,
     *                151.9828, 268.0157, 101.1302, 150.8691,
     *                 13.3718,   8.6528,  60.4614, 285.5072,
     *                  6.7655,   4.7684,  12.8631,  54.5260          /
      DATA GAR      /  12.638 ,  14.958 ,  12.833 ,  15.139 ,
     *                  0.178 ,  17.522 ,  23.584 ,  20.464 ,  25.150 ,
     *                  0.151 ,   1.561 ,  17.399 ,  30.871 ,
     *                 24.684 ,  33.978 ,  27.091 ,  36.481 ,
     *                  2.810 ,   8.877 ,  24.351 ,  44.489 ,
     *                  0.144 ,   1.160 ,  10.210 ,  27.178           /
C
      DATA AK       /  12.9782, 148.6673,   6.3493,
     *                 66.3444, 101.6553,   4.0001,  13.4465,  46.5534,
     *                  2.0171, 116.4767, 713.4965,  63.5907, 396.4079,
     *                  2.0000,  10.0000,  30.0000,
     *                  4.0702,   5.7791,  52.6795, 327.4539,
     *                 62.8604, 357.1331,  55.9337, 196.0646,
     *                 10.9275,   5.5398,  43.2761,  76.2560,
     *                 42.0000,  18.0000                              /
      DATA GK       /   1.871 ,   3.713 ,  18.172 ,
     *                 21.185 ,  27.705 ,   2.059 ,  23.709 ,  28.542 ,
     *                  0.273 ,  26.709 ,  39.640 ,  31.220 ,  41.865 ,
     *                 29.955 ,  37.557 ,  42.862 ,
     *                  0.228 ,   2.274 ,  21.703 ,  50.191 ,
     *                 32.145 ,  49.262 ,  34.155 ,  51.718 ,
     *                  3.043 ,   5.479 ,  20.547 ,  30.680 ,
     *                 36.275 ,  47.345                               /
C
      DATA ACA      /  18.2366,  27.5012, 149.2617,  94.5242, 705.4711,
     *                 11.8706,  14.0710, 106.0547,
     *                 57.2414, 110.7567,  29.8121,  54.1874,
     *                  2.0184,  97.5784, 282.3939, 209.1871, 252.8129/
      DATA GCA      /   2.050 ,   3.349 ,   5.321 ,   4.873 ,   7.017 ,
     *                  1.769 ,   5.109 ,   9.524 ,
     *                 27.271 ,  41.561 ,  29.172 ,  42.140 ,
     *                  0.394 ,  28.930 ,  52.618 ,  38.593 ,  49.646 /
C
      DATA ASC      /   6.0014,  83.1958,  67.3666, 329.4354,
     *                 44.0793, 169.9969, 533.9195,
     *                 34.1642, 124.8475, 228.9879,
     *                 11.9979,  16.9280,  28.4778,  82.0418, 234.5360,
     *                  6.0042,   2.7101,  13.9801,  65.3039,
     *                 12.0000,  12.0000,
     *                  2.0051,   2.9621,  29.0306                    /
      DATA GSC      /   0.021 ,   2.056 ,   3.551 ,   5.465 ,
     *                  1.535 ,   3.797 ,   6.203 ,
     *                  2.389 ,   4.858 ,   7.141 ,
     *                  0.011 ,   0.430 ,   1.156 ,   3.711 ,   8.863 ,
     *                  0.025 ,   3.499 ,  10.463 ,  18.606 ,
     *                 41.779 ,  57.217 ,
     *                  0.539 ,  24.442 ,  51.079                     /
C
      DATA ATI      /   7.0887,   8.9186,  17.5633, 206.6832, 438.5735,
     *                654.1721,
     *                 38.0462,  69.6271, 364.2845, 832.0408,
     *                 98.8562,  57.9934, 442.1498,
     *                 19.7843,  32.0637,  37.0895, 110.6682, 288.4946,
     *                521.8837,
     *                 10.0000,  34.0000, 120.0000,
     *                 16.1691,  22.3550,  24.1646,  83.5128, 222.7963,
     *                  6.0020,   4.6177,  25.2636,  52.1162,
     *                 12.0000,   8.0000                              /
      DATA GTI      /   0.021 ,   0.048 ,   1.029 ,   2.183 ,   4.109 ,
     *                  5.785 ,
     *                  0.846 ,   1.792 ,   3.836 ,   5.787 ,
     *                  2.561 ,   4.869 ,   6.340 ,
     *                  0.023 ,   0.124 ,   0.774 ,   1.810 ,   4.980 ,
     *                  9.585 ,
     *                  1.082 ,   4.928 ,  11.279 ,
     *                  0.041 ,   1.375 ,   4.768 ,  10.985 ,  19.769 ,
     *                  0.048 ,  11.577 ,  24.531 ,  36.489 ,
     *                 54.436 ,  75.373                               /
C
      DATA AV       /  15.2627,  23.9869,  51.3053, 570.3384,1650.9417,
     *                162.2829, 298.8303, 908.8852,
     *                 23.6736,  37.1624,  86.8011, 300.7440, 864.5880,
     *                 57.8961,  79.4605, 214.9007, 864.7425,
     *                 61.8508,  64.0845, 192.8298, 718.2349,
     *                 23.8116,  68.2495, 135.0613, 536.7632,
     *                 15.9543,  22.5542,  71.4921, 248.9544,
     *                  6.0006,   5.8785,  50.5077,  97.6129          /
      DATA GV       /   0.026 ,   0.145 ,   0.718 ,   2.586 ,   5.458 ,
     *                  2.171 ,   4.153 ,   6.097 ,
     *                  0.009 ,   0.366 ,   1.504 ,   5.294 ,  10.126 ,
     *                  1.796 ,   2.353 ,   6.068 ,  12.269 ,
     *                  2.560 ,   3.674 ,   6.593 ,  12.880 ,
     *                  0.045 ,   1.684 ,   8.162 ,  21.262 ,
     *                  0.065 ,   1.746 ,  15.158 ,  33.141 ,
     *                  0.077 ,  21.229 ,  44.134 ,  60.203           /
C
      DATA ACR      /  30.1842,  79.2847, 149.5293,
     *                215.3696, 119.1974, 741.4321,
     *                184.9946,1352.5038, 784.4937,
     *                 46.6191, 160.1361, 488.0449, 657.1928,
     *                 47.1742, 267.0275, 441.1324, 150.6650,
     *                 24.3768, 122.8359, 285.5092, 794.1654,
     *                 24.2296,  75.0258, 172.9452, 543.6511,
     *                 15.9819,  17.6800,  95.2003, 225.0947          /
      DATA GCR      /   0.993 ,   3.070 ,   5.673 ,
     *                  3.339 ,   4.801 ,   7.198 ,
     *                  2.829 ,   4.990 ,   7.643 ,
     *                  1.645 ,   3.727 ,   7.181 ,  12.299 ,
     *                  2.902 ,   4.273 ,   8.569 ,  14.912 ,
     *                  0.047 ,   2.566 ,   9.441 ,  21.198 ,
     *                  0.078 ,   2.242 ,  15.638 ,  32.725 ,
     *                  0.103 ,   2.146 ,  26.153 ,  49.381           /
C
      DATA AMN      /  53.9107,  81.3931, 546.6945 ,
     *                144.1893, 407.8029,  45.6177, 298.4423,2410.9335,
     *                 22.6382,  93.8419, 183.9367, 907.5765,
     *                137.0409, 168.6783, 329.0287, 773.2513,
     *                 70.1925,  72.3372, 213.9512, 539.5165,
     *                 24.2373,  93.5415, 456.6167, 506.5484,
     *                 24.7687,  66.9896, 264.1853, 484.0161          /
      DATA GMN      /   2.527 ,   4.204 ,   6.602 ,
     *                  4.155 ,   7.321 ,   2.285 ,   5.631 ,   8.448 ,
     *                  1.496 ,   3.839 ,   7.751 ,  13.484 ,
     *                  3.681 ,   6.054 ,   9.934 ,  14.936 ,
     *                  3.531 ,   6.967 ,  15.222 ,  25.069 ,
     *                  0.071 ,   2.896 ,  20.725 ,  37.383 ,
     *                  0.126 ,   2.660 ,  28.528 ,  53.413           /
C
      DATA AFE      /  14.4102,   2.7050, 421.6612, 940.1484,
     *                 36.2187,  22.8883, 239.5997, 825.2919,
     *                110.0242, 992.3040, 640.6715,
     *                 17.0494,  32.3783,  34.3184, 420.9626,1067.2064,
     *                154.0059, 462.1117, 329.8618,
     *                 15.7906,  47.1186, 279.9292, 692.1005,
     *                 91.0206, 206.3082, 706.9927, 836.6689,
     *                 40.0790,  27.6965,  28.2243,  18.0001,
     *                 24.0899,  89.6340,  51.5756, 241.6980          /
      DATA GFE      /   0.066 ,   0.339 ,   2.897 ,   6.585 ,
     *                  0.923 ,   1.679 ,   4.620 ,   7.053 ,
     *                  4.249 ,   5.875 ,   7.781 ,
     *                  0.062 ,   0.283 ,   1.504 ,   5.430 ,  11.210 ,
     *                  2.792 ,   7.627 ,  13.623 ,
     *                  0.077 ,   3.723 ,  12.137 ,  23.700 ,
     *                  2.688 ,   7.595 ,  15.444 ,  25.587 ,
     *                  3.982 ,   4.677 ,   6.453 ,  23.561 ,
     *                  0.102 ,   3.354 ,  22.954 ,  33.796           /
C
      DATA ACO      /  11.9120,  20.4424,  28.3863, 132.5038, 600.7461,
     *                 33.3092, 237.4331, 977.2502,
     *                 55.5396, 318.8169, 619.6366,
     *                 32.6900,  83.8694, 107.4378,
     *                 11.2593,  38.2239,  22.9964, 261.3486, 637.1485,
     *                 23.0233,  41.6599, 264.6460, 181.6699,
     *                 16.0356,   7.8633,  70.3158, 423.3512, 742.3553,
     *                  0.                                            /
      DATA GCO      /   0.112 ,   0.341 ,   0.809 ,   3.808 ,   6.723 ,
     *                  2.057 ,   3.484 ,   7.210 ,
     *                  2.405 ,   5.133 ,   8.097 ,
     *                  2.084 ,   5.291 ,   8.426 ,
     *                  0.135 ,   0.517 ,   1.606 ,   6.772 ,  12.622 ,
     *                  2.512 ,   4.348 ,   8.253 ,  15.377 ,
     *                  0.132 ,   0.863 ,   3.086 ,  11.789 ,  23.263 ,
     *                  0.                                            /
C
      DATA ANI      /   7.1268,  12.4486,  11.9953,  10.0546, 114.1658,
     *                391.2064,
     *                 26.3908, 213.8081, 938.7927,
     *                  4.1421,  37.3781,  25.9712, 333.3397, 311.1633,
     *                 33.1031, 184.1854, 136.7072,
     *                 11.1915,   5.4174,  53.6793, 460.6781, 380.0056,
     *                  0.                                            /
      DATA GNI      /   0.026 ,   0.137 ,   0.315 ,   1.778 ,   4.029 ,
     *                  6.621 ,
     *                  2.249 ,   4.042 ,   7.621 ,
     *                  0.191 ,   1.235 ,   3.358 ,   8.429 ,  17.096 ,
     *                  3.472 ,   9.065 ,  16.556 ,
     *                  0.194 ,   1.305 ,   5.813 ,  14.172 ,  26.169 ,
     *                  0.                                            /
C
      DATA ACU      /  11.0549, 238.9423,  10.3077, 126.2990,1073.3876,
     *                 30.0000,  50.0000,  60.0000,
     *                 19.2984,  50.5974, 240.2021,1216.9016,
     *                 48.3048, 583.2011, 320.4931,
     *                  4.0155,  70.3264, 313.1213, 536.5331,
     *                  0.                                            /
      DATA GCU      /   4.212 ,   7.227 ,   1.493 ,   5.859 ,   9.709 ,
     *                  7.081 ,   9.362 ,  10.130 ,
     *                  2.865 ,   8.260 ,  14.431 ,  18.292 ,
     *                  9.650 ,  14.640 ,  24.320 ,
     *                  0.337 ,   8.520 ,  16.925 ,  28.342 ,
     *                  0.                                            /
C
      DATA AZN      /  15.9880, 484.0042,  18.5863, 123.4134,
     *                  3.0000, 189.0000,
     *                  6.1902,  38.9317, 204.8780,
     *                 10.2588,  89.3771, 370.3640,  30.0000, 128.0000,
     *                 24.6904, 106.7491, 439.5586,
     *                  0.                                            /
      DATA GZN      /   4.546 ,   8.840 ,  10.247 ,  16.620 ,
     *                 11.175 ,  16.321 ,
     *                  6.113 ,  12.964 ,  16.444 ,
     *                  7.926 ,  13.633 ,  24.353 ,  16.286 ,  24.910 ,
     *                 10.291 ,  20.689 ,  32.077 ,
     *                  0.                                            /
C
      DATA ICOMP /0/
C
c      save indexs, indexm, index0, is, im, ig0, igpr,
c     *     xl, chion, alf, gam
C
      IF(ICOMP.NE.0) GO TO 5
      IND=1
      DO 1 K=1,NIONS
         INDEXS(K)=IND
         IND=IND+IS(K)
    1 CONTINUE
      IND=1
      DO 2 K=1,NSS
         INDEXM(K)=IND
         IND=IND+IM(K)
    2 CONTINUE
      ICOMP=1
    5 CONTINUE
c
      IF((IAT.EQ.26.or.iat.eq.28)
     *  .AND.IZI.GE.4.AND.IZI.LE.9) GO TO 70
      IF(IAT.GT.30.AND.IZI.LE.3) GO TO 80
      IF(IAT.GT.8 .AND. IZI.GT.5) then
         u=igle(iat-izi+1)
         return
      end if
c
c     Irwin partition functions by default
c
      if(iirwin.gt.0.and.t.lt.16000.) then
         if(izi.le.2) then
            call mpartf(iat,izi,0,t,u0)
            u=u0
            return
          end if
       else if(iat.gt.30.and.izi.le.3) then
         go to 80
      end if
c
      IF(IZI.LE.0.OR.IZI.GT.9.OR.IAT.LE.0.OR.IAT.GT.30) GO TO 50
      MODE=MODPF(IAT)
      IF(MODE.LT.0) GO TO 50
      IF(MODE.GT.0) GO TO 60
      I0=INDEX0(IZI,IAT)
      IF(I0) 40,50,10
   10 QZ=IZI
C     MAX=XMAXN*SQRT(QZ)
      XMAX=XMAXN
      THET=5040.4/T
      A=31.321*QZ*QZ*THET
      XMAX2=XMAX*XMAX
      QAS1=XMAX*THIRD*(XMAX2+TRHA*XMAX+HALF)
      IS0=INDEXS(I0)
      ISS=IS0+IS(I0)-1
      SU1=0.
      SQA=0.
      DO 30 K=IS0,ISS
         XXL=XL(K)
         GPR=IGPR(K)
         X=CHION(K)*THET
         EX=0.
         IF(X.LT.30) EX=EXP(-X*2.30258029299405)
         QAS=(QAS1-XXL*THIRD*(XXL*XXL+TRHA*XXL+HALF)+(XMAX-XXL)*
     *       (UN+A*HALF/XXL/XMAX)*A)*GPR*EX
         SQA=SQA+QAS
         M0=INDEXM(K)
         M1=M0+IM(K)-1
         AL1=0.
         DO 20 M=M0,M1
            XG=GAM(M)*THET
            IF(XG.GT.20.) GO TO 20
            XM=EXP(-XG*2.30258029299405)*ALF(M)
            AL1=AL1+XM
   20    CONTINUE
         SU1=SU1+AL1
   30 CONTINUE
      U=IG0(I0)
      U=U+SU1+SQA
      IF(U.LT.0.) U=IG0(I0)
      RETURN
   40 U=FLOAT(-I0)
      RETURN
   50 CALL PFSPEC(IAT,IZI,T,ANE,U)
      RETURN
   60 u=igle(iat-izi+1)
C      U=PFSTD(IZI,IAT)
      RETURN
   70 if(iat.eq.26) call pffe(IZI,T,ANE,U)
      if(iat.eq.28) call pfni(izi,t,u,dut,dun)
      RETURN
   80 CALL PFHEAV(IAT,IZI,3,T,ANE,U)
      RETURN
      END
C
C ********************************************************************
C
C
      subroutine pffe(ion,t,ane,pf)
c     =============================
c
c     partition functions for Fe IV to Fe IX
c     after Fischel and Sparks, 1971, NASA SP-3066
c
c     Output:  PF   partition function
c
      INCLUDE 'PARAMS.FOR'
      dimension tt(50),pn(10),nca(9)
      dimension p4a(22),p4b(10,28), 
     *          p5a(30),p5b(10,20),
     *          p6a(37),p6b(10,13),
     *          p7a(40),p7b(10,10),
     *          p8a(41),p8b(10,9),
     *          p9a(45),p9b(10,5)
c
      parameter (xen=2.302585093,xmil=0.001,xmilen=xmil*xen)
      parameter (xbtz=1.38054d-16)
c
      data nca /3*0,22,30,37,40,41,45/
     *     nne /10/
c
      data tt /
     * 3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,
     * 20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,
     * 32.,34.,36.,38.,40.,42.,44.,46.,48.,
     * 50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,125.,150./
c
      data pn /-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7./
c
      data p4a /
     * 0.778, 0.778, 0.778, 0.779, 0.783, 0.789, 0.801, 0.818,
     * 0.842, 0.871, 0.906, 0.945, 0.987, 1.030, 1.074, 1.117,
     * 1.160, 1.201, 1.242, 1.280, 1.317, 1.353/
c
      data p4b /
     * 1.406,1.393,1.389,7*1.387,
     * 1.464,1.434,1.424,1.421,1.420,5*1.419,
     * 1.546,1.483,1.461,1.454,1.451,1.451,4*1.450,
     * 1.665,1.547,1.503,1.488,1.482,1.481,4*1.480,
     * 1.826,1.636,1.553,1.524,1.514,1.510,4*1.509,
     * 2.024,1.755,1.618,1.564,1.546,1.540,1.538,3*1.537,
     * 2.480,2.087,1.814,1.674,1.619,1.599,1.593,1.591,1.590,1.590,
     * 2.945,2.489,2.105,1.846,1.717,1.667,1.649,1.643,1.641,1.640,
     * 3.379,2.897,2.452,2.089,1.859,1.751,1.710,1.696,1.691,1.689,
     * 3.774,3.283,2.808,2.381,2.054,1.864,1.782,1.751,1.741,1.738,
     * 4.133,3.637,3.150,2.688,2.292,2.015,1.871,1.814,1.793,1.786,
     * 4.460,3.962,3.468,2.989,2.549,2.199,1.984,1.886,1.848,1.835,
     * 4.757,4.258,3.762,3.274,2.809,2.406,2.121,1.972,1.908,1.886,
     * 5.029,4.530,4.032,3.539,3.061,2.624,2.279,2.073,1.976,1.939,
     * 5.279,4.780,4.281,3.785,3.299,2.840,2.450,2.189,2.051,1.996,
     * 5.510,5.010,4.511,4.013,3.522,3.050,2.628,2.318,2.136,2.057,
     * 6.014,5.514,5.014,4.515,4.018,3.530,3.065,2.666,2.381,2.228,
     * 6.435,5.935,5.435,4.936,4.437,3.943,3.460,3.022,2.658,2.422,
     * 6.794,6.294,5.794,5.294,4.794,4.297,3.807,3.343,2.939,2.631,
     * 7.102,6.602,6.102,5.602,5.102,4.604,4.110,3.638,3.194,2.845,
     * 7.370,6.870,6.370,5.870,5.370,4.871,4.375,3.892,3.439,3.052,
     * 7.606,7.106,6.606,6.106,5.605,5.106,4.608,4.125,3.661,3.249,
     * 7.815,7.315,6.814,6.314,5.814,5.314,4.816,4.333,3.851,3.418,
     * 8.001,7.501,7.001,6.500,6.000,5.500,5.001,4.511,4.032,3.586,
     * 8.168,7.668,7.168,6.668,6.168,5.667,5.168,4.680,4.197,3.741,
     * 8.319,7.819,7.319,6.819,6.319,5.818,5.319,4.832,4.347,3.884,
     * 8.900,8.399,7.899,7.399,6.899,6.398,5.898,5.405,4.917,4.431,
     * 9.294,8.794,8.294,7.793,7.293,6.793,6.292,5.799,5.306,4.824/
c
      data p5a /
     * 1.235, 1.276, 1.301, 1.321, 1.339, 1.359, 1.381, 1.405,
     * 1.432, 1.460, 1.489, 1.518, 1.546, 1.574, 1.601, 1.627,
     * 1.652, 1.675, 1.697, 1.718, 1.738, 1.757, 1.775, 1.792,
     * 1.808, 1.823, 1.838, 1.851, 1.877, 1.900/
c
      data p5b /
     * 1.943,1.928,1.923,7*1.921,
     * 2.011,1.964,1.947,1.942,1.941,5*1.940,
     * 2.144,2.025,1.980,1.965,1.960,1.958,4*1.957,
     * 2.361,2.137,2.032,1.993,1.980,1.976,1.975,3*1.974,
     * 2.646,2.315,2.121,2.035,2.004,1.994,1.991,1.990,1.989,1.989,
     * 2.960,2.553,2.260,2.102,2.037,2.015,2.007,2.005,2.004,2.004,
     * 3.274,2.823,2.450,2.205,2.086,2.040,2.025,2.020,2.018,2.018,
     * 3.575,3.101,2.674,2.348,2.158,2.075,2.045,2.036,2.032,2.031,
     * 4.251,3.757,3.275,2.829,2.466,2.234,2.124,2.083,2.069,2.064,
     * 4.822,4.324,3.829,3.346,2.895,2.522,2.278,2.161,2.116,2.100,
     * 5.308,4.808,4.310,3.816,3.334,2.888,2.525,2.297,2.187,2.145,
     * 5.725,5.225,4.726,4.228,3.736,3.260,2.828,2.496,2.294,2.206,
     * 6.088,5.589,5.089,4.590,4.093,3.604,3.139,2.733,2.447,2.291,
     * 6.407,5.907,5.407,4.908,4.409,3.915,3.433,2.988,2.629,2.399,
     * 6.689,6.189,5.689,5.189,4.690,4.193,3.704,3.236,2.832,2.535,
     * 6.940,6.440,5.940,5.440,4.941,4.443,3.949,3.469,3.038,2.687,
     * 7.166,6.666,6.166,5.666,5.166,4.667,4.171,3.684,3.237,2.847,
     * 7.370,6.870,6.369,5.869,5.369,4.870,4.373,3.882,3.417,3.008,
     * 8.150,7.649,7.149,6.649,6.149,5.649,5.149,4.651,4.167,3.700,
     * 8.677,8.177,7.676,7.176,6.676,6.176,5.676,5.176,4.687,4.203/
c
      data p6a /
     * 1.218, 1.273, 1.309, 1.335, 1.358, 1.379, 1.400, 1.421,
     * 1.442, 1.463, 1.484, 1.504, 1.523, 1.542, 1.560, 1.577,
     * 1.594, 1.609, 1.624, 1.638, 1.652, 1.664, 1.677, 1.688,
     * 1.699, 1.709, 1.719, 1.729, 1.746, 1.762, 1.777, 1.790,
     * 1.803, 1.814, 1.825, 1.834, 1.843/
c
      data p6b /
     * 1.862,1.855,1.853,7*1.852,
     * 1.958,1.900,1.880,1.874,1.872,5*1.871,
     * 2.264,2.045,1.944,1.906,1.894,1.890,4*1.888,
     * 2.776,2.386,2.119,1.984,1.930,1.912,1.906,1.904,2*1.903,
     * 3.321,2.856,2.453,2.165,2.012,1.949,1.927,1.920,1.918,1.917,
     * 3.821,3.333,2.868,2.465,2.178,2.025,1.963,1.941,1.934,1.932,
     * 4.266,3.771,3.285,2.825,2.434,2.164,2.027,1.972,1.953,1.947,
     * 4.662,4.164,3.670,3.187,2.739,2.372,2.135,2.022,1.980,1.965,
     * 5.015,4.516,4.019,3.527,3.052,2.624,2.295,2.102,2.019,1.988,
     * 5.332,4.832,4.344,3.838,3.351,2.889,2.493,2.217,2.075,2.017,
     * 5.618,5.118,4.619,4.121,3.628,3.149,2.711,2.364,2.155,2.058,
     * 6.710,6.210,5.710,5.210,4.711,4.213,3.719,3.241,2.807,2.462,
     * 7.446,6.946,6.446,5.946,5.446,4.946,4.447,3.952,3.474,3.022/
c
      data p7a /
     * 1.074,1.130,1.167,1.194,1.215,1.234,1.250,1.266,1.280,1.293,
     * 1.306,1.318,1.329,1.340,1.350,1.360,1.369,1.378,1.386,1.394,
     * 1.401,1.408,1.415,1.421,1.427,1.433,1.439,1.444,1.454,1.463,
     * 1.471,1.479,1.486,1.492,1.498,1.504,1.509,1.514,1.525,1.534/
c
      data p7b /
     * 1.555,1.546,1.544,1.543,6*1.542,
     * 1.617,1.572,1.557,1.552,1.550,1.550,4*1.549,
     * 1.798,1.648,1.587,1.566,1.559,1.557,4*1.556,
     * 2.134,1.832,1.666,1.597,1.573,1.565,1.563,1.562,2*1.561,
     * 2.550,2.138,1.836,1.671,1.602,1.578,1.570,1.568,2*1.567,
     * 2.968,2.504,2.102,1.816,1.665,1.603,1.582,1.575,2*1.572,
     * 3.359,2.875,2.419,2.037,1.779,1.651,1.601,1.584,1.579,1.577,
     * 3.718,3.224,2.745,2.305,1.953,1.736,1.636,1.599,1.586,1.582,
     * 5.097,4.598,4.098,3.601,3.110,2.638,2.217,1.899,1.719,1.643,
     * 6.026,5.526,5.026,4.527,4.028,3.531,3.042,2.576,2.170,1.885/
c
      data p8a /
     * 0.809,0.849,0.875,0.894,0.908,0.918,0.927,0.934,0.939,0.944,
     * 0.948,0.952,0.955,0.958,0.960,0.962,0.964,0.966,0.967,0.969,
     * 0.970,0.971,0.973,0.974,0.975,0.975,0.976,0.977,0.978,0.980,
     * 0.981,0.982,0.983,0.984,0.984,0.985,0.986,0.986,0.987,0.988,
     * 0.989/
c
      data p8b /
     * 0.992,0.991,8*0.990,
     * 1.000,0.994,0.992,7*0.991,
     * 1.032,1.005,0.996,0.993,0.992,5*0.991,
     * 1.129,1.040,1.008,0.997,0.993,5*0.992,
     * 1.335,1.132,1.042,1.009,0.998,0.994,0.993,0.993,2*0.992,
     * 1.640,1.312,1.121,1.038,1.007,0.998,0.994,3*0.993,
     * 1.987,1.573,1.269,1.101,1.030,1.005,0.997,2*0.994,0.993,
     * 3.514,3.017,2.526,2.053,1.628,1.305,1.119,1.039,1.010,1.000,
     * 4.569,4.069,3.569,3.072,2.580,2.103,1.671,1.336,1.136,1.048/
c
      data p9a /39*0.000,0.001,0.002,0.005,0.008,0.014,0.021/
c
      data p9b /
     * 2*0.032,8*0.031,
     * 0.048,0.045,8*0.044,
     * 0.076,0.065,0.061,0.060,6*0.059,
     * 1.128,0.722,0.429,0.271,0.207,0.184,0.177,0.174,2*0.173,
     * 2.696,2.200,1.712,1.249,0.848,0.564,0.415,0.354,0.333,0.327/
c
      na=nca(ion)
      nb=50-na
      pne=log10(ane*xbtz*t)
      t0=xmil*t
      j=1
      if(pne.lt.pn(1)) go to 15
      if(pne.gt.pn(nne)) then
        j1=nne
        j2=nne
        goto 16
      endif
      do 10 j=1,nne-1      
        if(pne.ge.pn(j).and.pne.lt.pn(j+1)) go to 15
   10 continue
   15 j1=j
      j2=j1+1
      if(pne.lt.pn(1)) j2=1
   16 do 20 i=1,49
         if(t0.ge.tt(i).and.t0.lt.tt(i+1)) go to 25
   20 continue
   25 i1=i
      i2=i+1
      if(t0.gt.tt(50)) then
        i1=50
        i2=50
      endif
      if(i2.le.na) then
         if(ion.eq.4) then
           px1=p4a(i1)
           px2=p4a(i1)
           py1=p4a(i2)
           py2=p4a(i2)
         else if(ion.eq.5) then
           px1=p5a(i1)
           px2=p5a(i1)
           py1=p5a(i2)
           py2=p5a(i2)
         else if(ion.eq.6) then
           px1=p6a(i1)
           px2=p6a(i1)
           py1=p6a(i2)
           py2=p6a(i2)
         else if(ion.eq.7) then
           px1=p7a(i1)
           px2=p7a(i1)
           py1=p7a(i2)
           py2=p7a(i2)
         else if(ion.eq.8) then
           px1=p8a(i1)
           px2=p8a(i1)
           py1=p8a(i2)
           py2=p8a(i2)
         else if(ion.eq.9) then
           px1=p9a(i1)
           px2=p9a(i1)
           py1=p9a(i2)
           py2=p9a(i2)
         endif
       else if(i1.eq.na) then
         if(ion.eq.4) then
           px1=p4a(i1)
           px2=p4a(i1)
           py1=p4b(j1,i2-na)
           py2=p4b(j2,i2-na)
         else if(ion.eq.5) then
           px1=p5a(i1)
           px2=p5a(i1)
           py1=p5b(j1,i2-na)
           py2=p5b(j2,i2-na)
         else if(ion.eq.6) then
           px1=p6a(i1)
           px2=p6a(i1)
           py1=p6b(j1,i2-na)
           py2=p6b(j2,i2-na)
         else if(ion.eq.7) then
           px1=p7a(i1)
           px2=p7a(i1)
           py1=p7b(j1,i2-na)
           py2=p7b(j2,i2-na)
         else if(ion.eq.8) then
           px1=p8a(i1)
           px2=p8a(i1)
           py1=p8b(j1,i2-na)
           py2=p8b(j2,i2-na)
         else if(ion.eq.9) then
           px1=p9a(i1)
           px2=p9a(i1)
           py1=p9b(j1,i2-na)
           py2=p9b(j2,i2-na)
         endif
      else
         if(ion.eq.4) then
           px1=p4b(j1,i1-na)
           px2=p4b(j2,i1-na)
           py1=p4b(j1,i2-na)
           py2=p4b(j2,i2-na)
         else if(ion.eq.5) then
           px1=p5b(j1,i1-na)
           px2=p5b(j2,i1-na)
           py1=p5b(j1,i2-na)
           py2=p5b(j2,i2-na)
         else if(ion.eq.6) then
           px1=p6b(j1,i1-na)
           px2=p6b(j2,i1-na)
           py1=p6b(j1,i2-na)
           py2=p6b(j2,i2-na)
         else if(ion.eq.7) then
           px1=p7b(j1,i1-na)
           px2=p7b(j2,i1-na)
           py1=p7b(j1,i2-na)
           py2=p7b(j2,i2-na)
         else if(ion.eq.8) then
           px1=p8b(j1,i1-na)
           px2=p8b(j2,i1-na)
           py1=p8b(j1,i2-na)
           py2=p8b(j2,i2-na)
         else if(ion.eq.9) then
           px1=p9b(j1,i1-na)
           px2=p9b(j2,i1-na)
           py1=p9b(j1,i2-na)
           py2=p9b(j2,i2-na)
         endif
      end if
      dlgunx=px2-px1
      px=px1+(pne-pn(j1))*dlgunx
      dlguny=py2-py1
      py=py1+(pne-pn(j1))*dlguny
      delt=tt(i2)-tt(i1)
      if(delt.ne.0.) then
        dlgut=(py-px)/delt
        pf=px+(t0-tt(i1))*dlgut
      else
        pf=px
      endif
      pf=exp(xen*pf)
      return
      end
         
C
C ********************************************************************
C ********************************************************************
C

      SUBROUTINE MATINV(A,N,NR)
C     =========================
C
C     Matrix inversion
C      by LU decomposition
C
C      A  -  matrix of actual size (N x N) and maximum size (NR x NR)
C            to be inverted;
C      Inversion is accomplished in place and the original matrix is
C      replaced by its inverse
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION A(NR,NR)
      IF(N.EQ.1) GO TO 250
      DO 50 I=2,N
         IM1=I-1
         DO 20 J=1,IM1
            JM1=J-1
            DIV=A(J,J)
            SUM=0.
            IF(JM1.LT.1) GO TO 20
            DO 10 K=1,JM1
   10          SUM=SUM+A(I,K)*A(K,J)
   20       A(I,J)=(A(I,J)-SUM)/DIV
         DO 40 J=I,N
            SUM=0.
            DO 30 K=1,IM1
   30          SUM=SUM+A(I,K)*A(K,J)
   40       A(I,J)=A(I,J)-SUM
   50 CONTINUE
      DO 80 II=2,N
         I=N+2-II
         IM1=I-1
         IF(IM1.LT.1) GO TO 80
         DO 70 JJ=1,IM1
            J=I-JJ
            JP1=J+1
            SUM=0.
            IF(JP1.GT.IM1) GO TO 70
            DO 60 K=JP1,IM1
   60          SUM=SUM+A(I,K)*A(K,J)
   70       A(I,J)=-A(I,J)-SUM
   80 CONTINUE
      DO 110 II=1,N
         I=N+1-II
         DIV=A(I,I)
         IP1=I+1
         IF(IP1.GT.N) GO TO 110
         DO 100 JJ=IP1,N
            J=N+IP1-JJ
            SUM=0.
            DO 90 K=IP1,J
   90          SUM=SUM+A(I,K)*A(K,J)
            A(I,J)=-SUM/DIV
  100    CONTINUE
  110 A(I,I)=1.0D0/A(I,I)
C
      DO 240 I=1,N
         DO 230 J=1,N
            K0=I
            IF(J.GE.I) GO TO 220
            SUM=0.
  200       DO 210 K=K0,N
  210          SUM=SUM+A(I,K)*A(K,J)
            GO TO 230
  220       K0=J
            SUM=A(I,K0)
            IF(K0.EQ.N) GO TO 230
            K0=K0+1
            GO TO 200
  230    A(I,J)=SUM
  240 CONTINUE
      RETURN
  250 A(1,1)=1.0D0/A(1,1)
      RETURN
      END
C
C
C     ****************************************************************
C
C

      SUBROUTINE LINEQS(A,B,X,N,NR)
C     =============================
C
C     Solution of the linear system A*X=B
C     by Gaussian elimination with partial pivoting
C
C     Input: A  - matrix of the linear system, with actual size (N x N),
C                and maximum size (NR x NR)
C            B  - the rhs vector
C     Output: X - solution vector
C     Note that matrix A and vector B are destroyed here !
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION A(NR,NR),B(NR),X(NR),D(MLEVEL)
      DIMENSION IP(MLEVEL)
      DO 70 I=1,N
         DO 10 J=1,N
   10       D(J)=A(J,I)
         IM1=I-1
         IF(IM1.LT.1) GO TO 40
         DO 30 J=1,IM1
            IT=IP(J)
            A(J,I)=D(IT)
            D(IT)=D(J)
            JP1=J+1
            DO 20 K=JP1,N
   20          D(K)=D(K)-A(K,J)*A(J,I)
   30    CONTINUE
   40    AM=ABS(D(I))
         IP(I)=I
         DO 50 K=I,N
            IF(AM.GE.ABS(D(K))) GO TO 50
            IP(I)=K
            AM=ABS(D(K))
   50    CONTINUE
         IT=IP(I)
         A(I,I)=D(IT)
         D(IT)=D(I)
         IP1=I+1
         IF(IP1.GT.N) GO TO 80
         DO 60 K=IP1,N
   60       A(K,I)=D(K)/A(I,I)
   70 CONTINUE
   80 DO 100 I=1,N
         IT=IP(I)
         X(I)=B(IT)
         B(IT)=B(I)
         IP1=I+1
         IF(IP1.GT.N) GO TO 110
         DO 90 J=IP1,N
   90       B(J)=B(J)-A(J,I)*X(I)
  100 CONTINUE
  110 DO 140 I=1,N
         K=N-I+1
         SUM=0.
         KP1=K+1
         IF(KP1.GT.N) GO TO 130
         DO 120 J=KP1,N
  120       SUM=SUM+A(K,J)*X(J)
  130    X(K)=(X(K)-SUM)/A(K,K)
  140 CONTINUE
      RETURN
      END
C
C
C     ****************************************************************
C
C

      FUNCTION EXPINT(X)
C     ==================
C
C     First exponential integral function E1(X)
C
      
      INCLUDE 'PARAMS.FOR'
      IF(X.GT.1.0) GO TO 10
      EXPINT=-LOG(X)-0.57721566+X*(0.99999193+X*(-0.24991055
     1        +X*(0.05519968+X*(-0.00976004+X*0.00107857))))
      RETURN
   10 EXPINT=EXP(-X)*((0.2677734343+X*(8.6347608925+X*
     1       (18.059016973+X*(8.5733287401+X))))/
     2       (3.9584969228+X*(21.0996530827+X*
     3       (25.6329561486+X*(9.5733223454+X)))))/X
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
      INCLUDE 'PARAMS.FOR'
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
C

      subroutine intrp(wltab,absop,wlgrid,abgrd,nfr,nfgrid)
c     =====================================================
c
c     a more efficient interpolation routine - using bisection
c
      INCLUDE 'PARAMS.FOR'
      dimension wltab(1),absop(1),wlgrid(1),abgrd(1)
      dimension yint(mfgrid),jint(mfgrid)
c
c      set up interpolation coefficients for an interpolation
c      by bisection
c
         fr1=wltab(1)
         fr2=wltab(nfr)
         do ij=1,nfgrid
            xint=wlgrid(ij)
            jl=0
            ju=nfr+1
   10       continue
            if(ju-jl.gt.1) then
               jm=(ju+jl)/2
               if((fr2.gt.fr1).eqv.(xint.gt.wltab(jm))) then
                  jl=jm
                else
                  ju=jm
               end if
               go to 10
            end if
            j=jl
            if(j.eq.nfr) j=j-1
            if(j.eq.0) j=j+1
            jint(ij)=j
c           yint(ij)=un/log10(wltab(j+1)/wltab(j))
            yint(ij)=1./(wltab(j+1)-wltab(j))
         end do
c
         do ij=1,nfgrid
            j=jint(ij)
            rc=(absop(j+1)-absop(j))*yint(ij)
c           abgrd(ij)=rc*log10(wlgrid(ij)/wltab(j))+absop(j)
            abgrd(ij)=rc*(wlgrid(ij)-wltab(j))+absop(j)
         end do
      return
      end

C
C ********************************************************************
C
      SUBROUTINE PFSPEC(IAT,IZI,T,ANE,U)
C     ==================================
                                                                                
C     Non-standard evaluation of the partition function
C     user supplied procedure
C
C     Input:
C      IAT   - atomic number
C      IZI   - ionic charge (=1 for neutrals, =1 for once ionized, etc)
C      T     - temperature
C      ANE   - electron density
C      XMAX  - principal quantum number of the last bound level
C
C     Output:
C      U     - partition function
C
*
* Modified from the ATMOS related programme 5-April-1990
* as an addition to TLUSTY to allow high ionisation states
* of C, N and O
*
* M.A.Barstow - University of Leicester, Dept of Physics & Astronomy
*
      INCLUDE 'PARAMS.FOR'
      real nvii
      PARAMETER (MH=100,MHEI=100,MHEII=100,MCI=135,
     +     MCII=157,MCIII=156,MCIV=55,MCV=15,MCVI=100,MNI=228,MNII=122,
     +     MNIII=133,MNIV=73,MNV=51,MNVI=8,MNVII=100,MOI=174,MOII=191,
     +     MOIII=168,MOIV=166,MOV=115,MOVI=52,MOVII=16,MOVIII=100)
      DIMENSION GHYD(MH),SHYD(MH),ENHYD(MH),
     +     GHEL(MH),ENHEL(MH),SHEL(MH),
     +     GCI(MCI),ENCI(MCI),SCI(MCI),
     +     GCII(MCII),ENCII(MCII),SCII(MCII),
     +     GCIII(MCIII),ENCIII(MCIII),SCIII(MCIII),
     +     GCIV(MCIV),ENCIV(MCIV),SCIV(MCIV),
     +     GCV(MCV),ENCV(MCV),SCV(MCV),
     +     GNI(MNI),ENNI(MNI),SNI(MNI),
     +     GNII(MNII),ENNII(MNII),SNII(MNII),
     +     GNIII(MNIII),ENNIII(MNIII),SNIII(MNIII),
     +     GNIV(MNIV),ENNIV(MNIV),SNIV(MNIV),
     +     GNV(MNV),ENNV(MNV),SNV(MNV),
     +     GNVI(MNVI),ENNVI(MNVI),SNVI(MNVI),
     +     GOI(MOI),ENOI(MOI),SOI(MOI),
     +     GOII(MOII),ENOII(MOII),SOII(MOII),
     +     GOIII(MOIII),ENOIII(MOIII),SOIII(MOIII),
     +     GOIV(MOIV),ENOIV(MOIV),SOIV(MOIV),
     +     GOV(MOV),ENOV(MOV),SOV(MOV),
     +     GOVI(MOVI),ENOVI(MOVI),SOVI(MOVI),
     +     GOVII(MOVII),ENOVII(MOVII),SOVII(MOVII)
      INTEGER NHYD(MH),NHEL(MHEI),NCI(MCI),NCII(MCII),
     +     NCIII(MCIII),NCIV(MCIV),NCV(MCV),NNI(MNI),NNII(MNII),
     +     NNIII(MNIII),NNIV(MNIV),NNV(MNV),NNVI(MNVI),NOI(MOI),
     +     NOII(MOII),NOIII(MOIII),NOIV(MOIV),NOV(MOV),NOVI(MOVI),
     +     NOVII(MOVII)
      PARAMETER (HI=13.5878,HEI=24.587,HEII=54.416,CVI=489.84,
     +     NVII=666.83,OVIII=871.12)
      PARAMETER (ZH=1.0,ZHE=2.0,ZC=6.0,ZN=7.0,ZO=8.0)
C                           N***=QUANTUM NO. OF LEVEL
C      DATA FOR IONS        G***=STATISTICAL WEIGHT OF LEVEL
C                           EN***=ENERGY OF LEVEL
C                           S*=SCREENING NO. OF LEVEL
        DATA NHYD/ 1, 2, 3, 4, 5, 6,
     +          7, 8, 9,10,11,12,
     +          13,14,15,16,17,18,
     +          19,20,21,22,23,24,
     +          25,26,27,28,29,30,
     +          31,32,33,34,35,36,
     +          37,38,39,40,41,42,
     +          43,44,45,46,47,48,
     +          49,50,51,52,53,54,
     +          55,56,57,58,59,60,
     +          61,62,63,64,65,66,
     +          67,68,69,70,71,72,
     +          73,74,75,76,77,78,
     +          79,80,81,82,83,84,
     +          85,86,87,88,89,90,
     +          91,92,93,94,95,96,
     +          97,98,99, 100 /
        DATA GHYD/ 2.000000, 8.000000, 18.00000,
     +           32.00000, 50.00000, 72.00000,
     +           98.00000, 128.0000, 162.0000,
     +           200.0000, 242.0000, 288.0000,
     +           338.0000, 392.0000, 450.0000,
     +           512.0000, 578.0000, 648.0000,
     +           722.0000, 800.0000, 882.0000,
     +           968.0000, 1058.000, 1152.000,
     +           1250.000, 1352.000, 1458.000,
     +           1568.000, 1682.000, 1800.000,
     +           1922.000, 2048.000, 2178.000,
     +           2312.000, 2450.000, 2592.000,
     +           2738.000, 2888.000, 3042.000,
     +           3200.000, 3362.000, 3528.000,
     +           3698.000, 3872.000, 4050.000,
     +           4232.000, 4418.000, 4608.000,
     +           4802.000, 5000.000, 5202.000,
     +           5408.000, 5618.000, 5832.000,
     +           6050.000, 6272.000, 6498.000,
     +           6728.000, 6962.000, 7200.000,
     +           7442.000, 7688.000, 7938.000,
     +           8192.000, 8450.000, 8712.000,
     +           8978.000, 9248.000, 9522.000,
     +           9800.000, 10082.00, 10368.00,
     +           10658.00, 10952.00, 11250.00,
     +           11552.00, 11858.00, 12168.00,
     +           12482.00, 12800.00, 13122.00,
     +           13448.00, 13778.00, 14112.00,
     +           14450.00, 14792.00, 15138.00,
     +           15488.00, 15842.00, 16200.00,
     +           16562.00, 16928.00, 17298.00,
     +           17672.00, 18050.00, 18432.00,
     +           18818.00, 19208.00, 19602.00,
     +           20000.00/
        DATA ENHYD /0.0000000E+00,10.19085000000000,12.07804444444444,
     +           12.73856250000000,13.04428800000000,13.21036111111111,
     +           13.31049795918367,13.37549062500000,13.42004938271605,
     +           13.45192200000000,13.47550413223140,13.49344027777778,
     +           13.50739881656805,13.51847448979592,13.52740977777778,
     +           13.53472265625000,13.54078339100346,13.54586234567901,
     +           13.55016066481994,13.55383050000000,13.55698866213152,
     +           13.55972603305785,13.56211417769376,13.56421006944444,
     +           13.56605952000000,13.56769970414201,13.56916104252401,
     +           13.57046862244898,13.57164328180737,13.57270244444444,
     +           13.57366077003122,13.57453066406250,13.57532268135905,
     +           13.57604584775087,13.57670791836735,13.57731558641975,
     +           13.57787465303141,13.57839016620499,13.57886653517423,
     +           13.57930762500000,13.57971683521713,13.58009716553288,
     +           13.58045127095727,13.58078150826446,13.58108997530864,
     +           13.58137854442344,13.58164889090086,13.58190251736111,
     +           13.58214077467722,13.58236488000000,13.58257593233372,
     +           13.58277492603550,13.58296276254895,13.58314026063100,
     +           13.58330816528926,13.58346715561225,13.58361785164666,
     +           13.58376082045184,13.58389658144211,13.58402561111111,
     +           13.58414834721849,13.58426519250780,13.58437651801461,
     +           13.58448266601563,13.58458395266272,13.58468067033976,
     +           13.58477308977501,13.58486146193772,13.58494601974375,
     +           13.58502697959184,13.58510454274945,13.58517889660494,
     +           13.58525021580034,13.58531866325785,13.58538439111111,
     +           13.58544754155125,13.58550824759656,13.58556663379356,
     +           13.58562281685627,13.58567690625000,13.58572900472489,
     +           13.58577920880428,13.58582760923211,13.58587429138322,
     +           13.58591933564014,13.58596281773932,13.58600480908971,
     +           13.58604537706612,13.58608458527964,13.58612249382716,
     +           13.58615915952180,13.58619463610586,13.58622897444791,
     +           13.58626222272522,13.58629442659280,13.58632562934028,
     +           13.58635587203741,13.58638519366930,13.58641363126212,
     +           13.58644122000000/
        DATA SHYD/100*0.0D0/
      DATA NHEL/1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,
     +        5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,
     +        23,24,25,26,27,
     +          28,29,30,31,32,33,
     +          34,35,36,37,38,39,
     +          40,41,42,43,44,45,
     +          46,47,48,49,50,51,
     +          52,53,54,55,56,57,
     +          58,59,60,61,62,63,
     +          64,65,66,67,68,69,
     +          70,71,72,73,74,75,
     +          76,77,78,79,80,81/
      DATA GHEL/1.0D0,3.0D0,1.0D0,5.0D0,3.0D0,1.0D0,3.0D0,
     +          3.0D0,1.0D0,5.0D0,3.0D0,
     +         1.0D0,15.0D0,5.0D0,3.0D0,3.0D0,1.0D0,9.0D0,
     +         15.0D0,5.0D0,21.0D0,7.0D0,
     +         3.0D0,100.0D0,144.0D0,196.0D0,256.0D0,324.0D0,
     +         400.0D0,484.0D0,
     +         576.0D0,676.0D0,784.0D0,900.0D0,1024.0D0,1156.0D0,
     +         1296.0D0,1444.0D0,1600.0D0,1764.0D0,1936.0D0,
     +         2116.0D0,2304.0D0,2500.0D0,2704.0D0,3136.0D0,
     +          3136.000000000000,3364.000000000000,3600.000000000000,
     +          3844.000000000000,4096.000000000000,4356.000000000000,
     +          4624.000000000000,4900.000000000000,5184.000000000000,
     +          5476.000000000000,5776.000000000000,6084.000000000000,
     +          6400.000000000000,6724.000000000000,7056.000000000000,
     +          7396.000000000000,7744.000000000000,8100.000000000000,
     +          8464.000000000000,8836.000000000000,9216.000000000000,
     +          9604.000000000000,10000.00000000000,10404.00000000000,
     +          10816.00000000000,11236.00000000000,11664.00000000000,
     +          12100.00000000000,12544.00000000000,12996.00000000000,
     +          13456.00000000000,13924.00000000000,14400.00000000000,
     +          14884.00000000000,15376.00000000000,15876.00000000000,
     +          16384.00000000000,16900.00000000000,17424.00000000000,
     +          17956.00000000000,18496.00000000000,19044.00000000000,
     +          19600.00000000000,20164.00000000000,20736.00000000000,
     +          21316.00000000000,21904.00000000000,22500.00000000000,
     +          23104.00000000000,23716.00000000000,24336.00000000000,
     +          24964.00000000000,25600.00000000000,26244.00000000000/
      DATA ENHEL/0.0D0,19.819D0,20.615D0,20.964D0,
     +           20.964D0,20.964D0,21.218D0,
     +           22.718D0,22.920D0,23.007D0,23.007D0,
     +           23.007D0,23.073D0,23.074D0,
     +           23.087D0,23.593D0,23.673D0,23.707D0,
     +           23.736D0,23.736D0,23.737D0,
     +           23.737D0,23.742D0,24.028D0,24.201D0,
     +           24.304D0,24.371D0,24.417D0,
     +           24.449D0,24.473D0,24.491D0,24.506D0,
     +           24.517D0,24.526D0,24.534D0,
     +           24.540D0,24.545D0,24.549D0,24.553D0,
     +           24.556D0,24.559D0,24.562D0,
     +           24.564D0,24.566D0,24.568D0,24.570D0,
     +          24.57131951530612,24.57238228299643,24.57334055555556,
     +          24.57420759625390,24.57499462890625,24.57571120293848,
     +          24.57636548442907,24.57696448979592,24.57751427469136,
     +          24.57802008765522,24.57848649584488,24.57891748849441,
     +          24.57931656250000,24.57968679357525,24.58003089569161,
     +          24.58035127095727,24.58065005165289,24.58092913580247,
     +          24.58119021739130,24.58143481213219,24.58166427951389,
     +          24.58187984173261,24.58208260000000,24.58227354863514,
     +          24.58245358727811,24.58262353150587,24.58278412208505,
     +          24.58293603305785,24.58307987882653,24.58321622037550,
     +          24.58334557074911,24.58346839988509,24.58358513888889,
     +          24.58369618382155,24.58380189906348,24.58390262030738,
     +          24.58399865722656,24.58409029585799,24.58417780073462,
     +          24.58426141679661,24.58434137110727,24.58441787439614,
     +          24.58449112244898,24.58456129736163,24.58462856867284,
     +          24.58469309438919,24.58475502191381,24.58481448888889,
     +          24.58487162396122,24.58492654747850,24.58497937212360,
     +          24.58503020349303,24.58507914062500,24.58512627648224/
      DATA SHEL/0.375D0,0.622D0,0.622D0,0.842D0,
     +          0.842D0,0.842D0,0.842D0,0.747D0,
     +          0.747D0,0.912D0,0.912D0,0.912D0,
     +          0.993D0,0.993D0,0.912D0,0.810D0,
     +          0.810D0,0.937D0,0.995D0,0.995D0,
     +          1.000D0,1.000D0,0.937D0,0.949D0,
     +          0.958D0,75*1.000D0/
        DATA NCI/2,2,2,2,2,2,3,3,3,3,2,2,2,3,3,3,3,3,
     +          3,3,3,3,3,2,3,4,4,4,3,3,3,3,3,3,4,3,
     +          3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +          4,4,4,5,4,4,4,4,4,5,5,5,5,5,5,5,5,5,
     +          5,5,5,5,6,5,5,5,5,5,6,6,6,6,6,6,6,7,
     +          6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,
     +          8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,10,10,
     +          10,11,11,11,2,3,3,3,2,2/
        DATA GCI/1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,    
     +          5.0D0,1.0D0,3.0D0,5.0D0,3.0D0,    
     +          7.0D0,5.0D0,3.0D0,3.0D0,3.0D0,    
     +          5.0D0,7.0D0,3.0D0,1.0D0,3.0D0,    
     +          5.0D0,5.0D0,1.0D0,9.0D0,5.0D0,    
     +          1.0D0,3.0D0,5.0D0,5.0D0,7.0D0,    
     +          9.0D0,3.0D0,5.0D0,7.0D0,3.0D0,    
     +          3.0D0,3.0D0,5.0D0,3.0D0,1.0D0,    
     +          3.0D0,5.0D0,7.0D0,3.0D0,3.0D0,    
     +          1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,    
     +          5.0D0,5.0D0,7.0D0,9.0D0,3.0D0,    
     +          5.0D0,7.0D0,3.0D0,7.0D0,3.0D0,   
     +          5.0D0,3.0D0,1.0D0,3.0D0,3.0D0,   
     +          5.0D0,7.0D0,5.0D0,1.0D0,5.0D0,   
     +          5.0D0,7.0D0,9.0D0,3.0D0,5.0D0,   
     +          7.0D0,3.0D0,7.0D0,3.0D0,5.0D0,   
     +          3.0D0,1.0D0,5.0D0,5.0D0,7.0D0,   
     +          9.0D0,3.0D0,5.0D0,7.0D0,3.0D0,   
     +          7.0D0,5.0D0,3.0D0,1.0D0,3.0D0,   
     +          5.0D0,7.0D0,9.0D0,3.0D0,5.0D0,   
     +          7.0D0,7.0D0,3.0D0,5.0D0,3.0D0,   
     +          1.0D0,9.0D0,7.0D0,5.0D0,3.0D0,   
     +          5.0D0,7.0D0,7.0D0,5.0D0,3.0D0,   
     +          1.0D0,9.0D0,7.0D0,5.0D0,3.0D0,   
     +          5.0D0,7.0D0,7.0D0,3.0D0,5.0D0,   
     +          7.0D0,3.0D0,5.0D0,7.0D0,5.0D0,   
     +          3.0D0,5.0D0,7.0D0,3.0D0,3.0D0/
        DATA ENCI/0.0D0,2.0333605D-03,5.3933649D-03,1.263870,2.684086,
     +          4.182672,7.480511,7.482891,7.487915,7.684888,
     +          7.946046,7.946620,7.946474,8.537387,8.640516,
     +          8.643146,8.647287,8.771255,8.846707,8.848247,
     +          8.850785,9.002712,9.171972,9.330682,9.631248,
     +          9.683908,9.685375,9.689256,9.695577,9.697620,
     +          9.701885,9.708156,9.708925,9.710041,9.712769,
     +          9.714380,9.761111,9.833419,9.834406,9.834934,
     +          9.940317,9.942698,9.946449,9.988707,10.05592,
     +          10.08144,10.08328,10.08553,10.13833,10.19809,
     +          10.35278,10.38514,10.38514,10.38514,10.39370,
     +          10.39456,10.39580,10.40021,10.40845,10.41874,
     +          10.42750,10.42990,10.42990,10.52043,10.52041,
     +          10.52041,10.53705,10.58840,10.61635,10.67973,
     +          10.70230,10.70328,10.70328,10.70878,10.70878,
     +          10.71184,10.71407,10.71854,10.72362,10.72523,
     +          10.72684,10.72684,10.86509,10.87426,10.87513,
     +          10.87513,10.87997,10.87997,10.88257,10.88533,
     +          10.88679,10.88964,10.89075,10.89075,10.88980,
     +          10.97789,10.97854,10.97854,10.98597,10.98597,
     +          10.98597,10.98808,10.98913,10.98994,10.98994,
     +          10.98994,11.04474,11.04474,11.04487,11.05280,
     +          11.05280,11.05280,11.05392,11.05429,11.05429,
     +          11.05429,11.09049,11.09049,11.09049,11.09843,
     +          11.09843,11.09843,11.09880,11.13129,11.13129,
     +          11.13129,11.15477,11.15477,11.15477,12.13544,
     +          12.83767,12.84024,12.84331,13.11772,14.86312/
        DATA NCII/2,2,2,2,2,2,2,2,2,2,3,3,3,2,3,3,2,2,4,4,4,3,3,3,
     +          4,4,2,2,4,4,5,5,5,3,3,5,5,5,5,6,3,3,3,3,3,3,6,6,
     +          6,6,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +          4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     +          5,5,5,5,5,5,6,6,6,6,6,6,6/
        DATA GCII/2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          6.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +          2.0D0,2.0D0,4.0D0,4.0D0,4.0D0,
     +          6.0D0,6.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +          4.0D0,6.0D0,6.0D0,8.0D0,2.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +          4.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          10.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          4.0D0,6.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0,8.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,4.0D0,6.0D0,4.0D0,
     +          6.0D0,8.0D0,10.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,4.0D0,6.0D0,6.0D0,
     +          4.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +          6.0D0,8.0D0,10.0D0,6.0D0,8.0D0,
     +          6.0D0,8.0D0,10.0D0,12.0D0,8.0D0,
     +          10.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0,4.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,6.0D0,4.0D0,
     +          2.0D0,6.0D0,8.0D0,4.0D0,6.0D0,
     +          8.0D0,10.0D0,6.0D0,8.0D0,10.0D0,
     +          12.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,6.0D0,
     +          4.0D0,2.0D0/
        DATA ENCII/0.0D0,7.9350658D-03,5.331397,5.334075,5.337658,
     +          9.290338,9.290624,11.96386,13.71590,13.72101,
     +          14.44900,16.33194,16.33332,17.60895,18.04607,
     +          18.04625,18.65519,18.65582,19.49478,20.14995,
     +          20.15068,20.70119,20.70413,20.70971,20.84491,
     +          20.84496,20.92025,20.92256,20.95094,20.95094,
     +          21.49265,21.73314,21.73405,22.09347,22.13075,
     +          22.13075,22.13075,22.18799,22.18799,22.47211,
     +          22.52747,22.52929,22.53239,22.53689,22.56844,
     +          22.57086,22.82136,22.82136,22.85996,22.85996,
     +          22.89870,23.11398,23.11600,23.11878,23.38108,
     +          23.38522,24.12408,24.27024,24.27201,24.27444,
     +          24.27787,24.37010,24.37079,24.37187,24.37315,
     +          24.60198,24.60332,24.65351,24.65617,24.65793,
     +          24.78982,24.79512,25.06741,25.07039,25.98117,
     +          25.98415,25.98986,26.58329,26.58615,26.62689,
     +          26.62867,26.63139,26.63554,26.75178,26.82771,
     +          26.82771,26.83016,26.89454,26.89578,27.22147,
     +          27.22329,27.22585,27.22930,27.29263,27.29263,
     +          27.29378,27.29509,27.35131,27.35294,27.37703,
     +          27.37957,27.38104,27.41188,27.41302,27.41395,
     +          27.41395,27.41395,27.41409,27.46301,27.46301,
     +          27.46810,27.46936,27.47200,27.47561,27.47330,
     +          27.47864,27.48713,27.49096,27.49330,27.49330,
     +          27.48854,27.49412,27.55688,27.56022,27.99752,
     +          27.99752,27.99752,28.25640,28.25640,28.61124,
     +          28.61124,28.61124,28.61124,28.64683,28.64683,
     +          28.64683,28.66803,26.43629,28.66875,28.66875,
     +          28.66875,28.66875,28.70253,28.70253,28.70253,
     +          28.70253,28.70515,28.70515,28.70515,28.70515,
     +          29.31561,29.31561,29.31561,29.31561,29.33557,
     +          29.33557,29.33557/
        DATA NCIII/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,
     +          3,4,4,4,4,3,4,4,4,4,4,4,4,4,3,3,3,4,3,3,3,3,3,3,
     +          3,3,3,3,3,3,5,3,3,3,3,5,5,5,5,3,5,5,5,5,5,5,5,5,
     +          5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,
     +          7,8,8,8,8,9,9,9,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +          4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,
     +          6,6,6,6,6,6,7,7,7,7,7,7/
        DATA GCIII/1.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +          1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +          3.0D0,1.0D0,3.0D0,1.0D0,3.0D0,
     +          5.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +          1.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +          3.0D0,5.0D0,7.0D0,5.0D0,7.0D0,
     +          9.0D0,3.0D0,7.0D0,3.0D0,5.0D0,
     +          7.0D0,5.0D0,3.0D0,1.0D0,3.0D0,
     +          5.0D0,5.0D0,5.0D0,5.0D0,7.0D0,
     +          9.0D0,3.0D0,5.0D0,7.0D0,3.0D0,
     +          5.0D0,3.0D0,1.0D0,7.0D0,3.0D0,
     +          1.0D0,3.0D0,5.0D0,1.0D0,3.0D0,
     +          5.0D0,7.0D0,7.0D0,9.0D0,11.0D0,
     +          9.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +          9.0D0,7.0D0,3.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,7.0D0,9.0D0,11.0D0,
     +          9.0D0,5.0D0,5.0D0,7.0D0,9.0D0,
     +          7.0D0,3.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,3.0D0,5.0D0,7.0D0,1.0D0,
     +          3.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,1.0D0,3.0D0,5.0D0,5.0D0,
     +          5.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,1.0D0,3.0D0,5.0D0,
     +          5.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +          5.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +          7.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,5.0D0,3.0D0,1.0D0,
     +          3.0D0,5.0D0,7.0D0,1.0D0,3.0D0,5.0D0/
        DATA ENCIII/0.0D0,6.486296,6.489148,6.496191,12.69008,
     +          17.03237,17.03602,17.04185,18.08638,22.62984,
     +          29.52845,30.64541,32.10371,32.19328,32.19396,
     +          32.19555,33.47080,33.45866,33.47146,34.27982,
     +          38.20770,38.21183,38.22034,38.36164,38.43612,
     +          38.64882,39.39549,39.39549,39.39611,39.64054,
     +          39.84380,39.84582,39.84874,39.91699,39.91782,
     +          39.91892,39.97328,40.01022,40.05026,40.05341,
     +          40.05822,40.19756,40.57121,40.86969,40.87231,
     +          40.87686,41.24874,41.30157,41.32848,41.33158,
     +          41.33611,41.85783,41.80309,41.86202,42.14028,
     +          42.16117,42.16444,42.16623,42.32471,42.55869,
     +          42.67342,42.67342,42.67342,42.78661,42.83001,
     +          42.83001,42.83001,42.96405,42.96405,42.96416,
     +          42.96405,42.98029,42.98736,43.03527,43.03550,
     +          43.03579,43.25349,43.98952,44.27370,44.39248,
     +          44.39248,44.39248,44.46592,44.46592,44.46600,
     +          44.47219,44.47673,44.48596,44.48596,44.48596,
     +          44.52591,45.07626,45.24178,45.32720,45.32720,
     +          45.32720,45.38200,45.86543,45.92891,45.92891,
     +          45.92891,46.33929,46.33929,46.33929,46.69749,
     +          46.69749,46.69749,47.25143,47.35238,47.35238,
     +          47.35722,47.64920,47.64920,47.65379,47.81342,
     +          47.83558,48.06245,48.06245,48.06245,48.16114,
     +          48.16114,48.16114,48.20208,50.51542,50.55803,
     +          50.55803,50.55803,50.69428,50.69428,50.69428,
     +          50.77264,50.79460,50.90022,50.90022,50.90022,
     +          50.93829,50.93829,50.93829,52.24497,52.24497,
     +          52.24497,52.31775,52.31775,52.31775,52.43107,
     +          52.43107,52.43107,52.45302,52.45302,52.45302,
     +          53.23251,53.23251,53.23251,53.27802,53.27802,
     +          53.27802/
        DATA NCIV/2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,
     +          6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,
     +          8,8,8,8,8,8,8/
        DATA GCIV/2.0D0,2.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +          2.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +          6.0D0,8.0D0,8.0D0,10.0D0,2.0D0,
     +          2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +          8.0D0,8.0D0,10.0D0,10.0D0,12.0D0,
     +          2.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +          6.0D0,8.0D0,8.0D0,10.0D0,10.0D0,
     +          12.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          8.0D0,10.0D0,12.0D0,14.0D0,16.0D0/
        DATA ENCIV/0.0D0,7.995100,8.008378,37.54872,39.68134,
     +          39.68525,40.28040,40.28173,49.76113,50.62434,
     +          50.62599,50.87540,50.87595,50.88784,50.88784,
     +          55.21889,55.65134,55.65221,55.77947,55.77947,
     +          55.78577,55.78578,55.78703,55.78703,58.12002,
     +          58.36774,58.36774,58.44275,58.44275,58.44709,
     +          58.44709,58.44764,58.44764,58.44770,58.44770,
     +          59.84267,60.00038,60.00038,60.04725,60.04725,
     +          60.05156,60.05156,60.05191,60.05191,60.05194,
     +          60.05194,61.05946,61.05946,61.09294,61.09294,
     +          61.09319,61.09319,61.09319,61.09319,61.09319/
        DATA NCV/1,2,2,2,2,2,3,3,3,3,4,5,6,7,8/
        DATA GCV/1.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +          3.0D0,3.0D0,5.0D0,7.0D0,3.0D0,
     +          3.0D0,3.0D0,3.0D0,3.0D0,3.0D0/
        DATA ENCV/0.0D0,298.9618,304.4046,304.4030,304.4199,
     +          307.8855,354.2645,354.2645,354.2645,354.5177,
     +          370.9247,378.5349,382.6710,385.1917,386.6807/
        DATA NNI/2,2,2,2,2,3,3,3,3,3,2,2,2,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,3,3,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,5,4,4,4,4,4,4,4,4,
     +          4,4,4,4,4,4,4,4,4,3,3,3,3,6,6,6,6,6,5,5,5,5,5,5,
     +          5,5,5,5,5,5,5,5,5,5,5,7,7,7,7,7,6,6,6,6,6,6,6,6,
     +          6,6,6,6,6,6,6,6,6,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,
     +          7,7,7,9,9,9,9,9,8,8,8,8,8,8,8,8,8,8,8,8,8,10,10,10,
     +          10,10,9,9,9,9,9,9,9,9,9,9,9,9,9,11,11,11,11,11,10,
     +          10,10,10,10,10,10,10,10,10,10,10,10,12,12,12,12,12,
     +          11,11,11,11,11,11,11,11,11,11,11,11,11,13,13,12,12,
     +          12,12,12,12,12/
        DATA GNI/4.0D0,6.0D0,4.0D0,4.0D0,2.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          6.0D0,4.0D0,2.0D0,2.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +          6.0D0,4.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,10.0D0,6.0D0,
     +          8.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,4.0D0,6.0D0,
     +          2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          2.0D0,4.0D0,6.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,4.0D0,
     +          6.0D0,8.0D0,10.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,4.0D0,2.0D0,6.0D0,
     +          8.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +          6.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          4.0D0,6.0D0,8.0D0,10.0D0,4.0D0,
     +          2.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,2.0D0,4.0D0,6.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,4.0D0,6.0D0,8.0D0,
     +          10.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          4.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,6.0D0,8.0D0,
     +          4.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +          8.0D0,4.0D0,2.0D0,6.0D0,8.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +          2.0D0,6.0D0,8.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0,8.0D0,2.0D0,4.0D0,6.0D0,
     +          8.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,4.0D0,2.0D0,6.0D0,8.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +          6.0D0,4.0D0,6.0D0/
        DATA ENNI/0.0D0,2.383371,2.384363,3.575739,3.575739,
     +          10.32619,10.33038,10.33617,10.67904,10.69042,
     +          10.92429,10.92973,10.93217,11.60284,11.75037,
     +          11.75317,11.75780,11.76412,11.83769,11.83997,
     +          11.84472,11.99580,12.00032,12.00975,12.12207,
     +          12.12649,12.35701,12.35614,12.84713,12.85333,
     +          12.86185,12.91211,12.92268,12.97078,12.97568,
     +          12.97693,12.97929,12.98350,12.98958,12.99502,
     +          13.00392,13.00161,13.00483,13.00074,13.01686,
     +          13.01822,13.01983,13.02095,13.03344,13.03636,
     +          13.20179,13.23674,13.23917,13.24364,13.25041,
     +          13.26429,13.26623,13.27127,13.32189,13.61527,
     +          13.62076,13.62945,13.64202,13.65185,13.66270,
     +          13.66493,13.66914,13.67609,13.66580,13.67249,
     +          13.67410,13.68043,13.66588,13.66872,13.67695,
     +          13.68464,13.67869,13.68191,13.68836,13.69398,
     +          13.69673,13.70310,13.70607,13.92292,13.92614,
     +          13.95653,13.96207,13.97100,13.97749,13.98841,
     +          13.97948,13.98097,13.98543,13.99324,13.98568,
     +          13.98754,13.98803,13.99674,13.98865,13.98865,
     +          13.98865,13.99696,13.99237,13.99473,13.99944,
     +          14.00155,14.00384,14.13620,14.14326,14.15244,
     +          14.15045,14.15455,14.15417,14.15417,14.15417,
     +          14.15417,14.15690,14.15690,14.15690,14.16508,
     +          14.15827,14.16025,14.15864,14.16843,14.16313,
     +          14.17035,14.16645,14.16645,14.16831,14.23464,
     +          14.24468,14.25113,14.25212,14.25212,14.25683,
     +          14.25683,14.25683,14.25683,14.25882,14.25882,
     +          14.26043,14.26043,14.26545,14.27073,14.27109,
     +          14.27109,14.27109,14.36247,14.36247,14.31821,
     +          14.31821,14.31821,14.32329,14.32329,14.32329,
     +          14.32329,14.32403,14.32403,14.32465,14.32465,
     +          14.33234,14.33544,14.33494,14.33494,14.33494,
     +          14.36272,14.36272,14.36433,14.36433,14.36433,
     +          14.36830,14.36830,14.36830,14.36830,14.36854,
     +          14.36854,14.37016,14.37016,14.37896,14.38119,
     +          14.38107,14.38107,14.38107,14.39557,14.39557,
     +          14.39768,14.39768,14.39768,14.40152,14.40152,
     +          14.40202,14.40202,14.40264,14.40264,14.40264,
     +          14.40264,14.41206,14.41206,14.41442,14.41442,
     +          14.41442,14.42012,14.42012,14.42099,14.42099,
     +          14.42099,14.42583,14.42583,14.42682,14.42682,
     +          14.42781,14.42781,14.42781,14.42781,14.43636,
     +          14.43636,14.43698,14.43698,14.43698,14.46253,
     +          14.44021,14.44455,14.44455,14.45434,14.45434,
     +          14.45434,14.45980,14.45980/
        DATA NNII/2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,2,3,3,3,3,2,3,
     +          3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,3,4,4,4,
     +          4,4,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +          4,4,4,4,4,4,4,4,4,4,3,3,3,5,5,5,5,5,5,5,5,5,5,5,
     +          5,5,5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/
        DATA GNII/1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +          5.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +          3.0D0,1.0D0,5.0D0,1.0D0,3.0D0,
     +          5.0D0,3.0D0,3.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,3.0D0,3.0D0,1.0D0,
     +          3.0D0,5.0D0,5.0D0,1.0D0,5.0D0,
     +          7.0D0,9.0D0,5.0D0,3.0D0,5.0D0,
     +          7.0D0,5.0D0,3.0D0,1.0D0,7.0D0,
     +          3.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +          3.0D0,3.0D0,5.0D0,7.0D0,1.0D0,
     +          3.0D0,5.0D0,3.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,1.0D0,5.0D0,7.0D0,
     +          9.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +          5.0D0,3.0D0,1.0D0,7.0D0,5.0D0,
     +          7.0D0,9.0D0,7.0D0,7.0D0,9.0D0,
     +          11.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +          3.0D0,5.0D0,1.0D0,3.0D0,5.0D0,
     +          1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,5.0D0,7.0D0,9.0D0,
     +          7.0D0,7.0D0,9.0D0,11.0D0,9.0D0,
     +          1.0D0,3.0D0,5.0D0,7.0D0,9.0D0,
     +          3.0D0,5.0D0,7.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,9.0D0,11.0D0,7.0D0,
     +          5.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +          7.0D0,9.0D0/
        DATA ENNII/0.0D0,6.0876831D-03,1.6279284D-02,1.898923,4.052723,
     +          5.848106,11.43604,11.43781,11.43801,13.54146,
     +          13.54146,13.54228,17.87734,18.46259,18.46651,
     +          18.48341,18.49722,19.23384,20.40944,20.64636,
     +          20.65389,20.66582,20.67651,20.94027,21.14861,
     +          21.15298,21.16022,21.59986,22.10340,23.12481,
     +          23.13218,23.14229,23.19670,23.23962,23.24260,
     +          23.24636,23.41565,23.42207,23.42555,23.47490,
     +          23.57225,24.36823,24.37465,24.38944,24.53166,
     +          25.06612,25.13369,25.14001,25.15193,25.18946,
     +          25.19245,25.20124,25.23510,25.46049,25.53877,
     +          25.54572,25.55447,25.58160,25.99668,26.00464,
     +          26.01527,26.02787,26.06667,26.06994,26.07548,
     +          26.12440,26.13011,26.13327,26.16475,26.16510,
     +          26.16800,26.16849,26.17391,26.19663,26.19758,
     +          26.20937,26.20252,26.21087,26.21191,26.21252,
     +          26.22134,26.22182,26.25393,26.25770,26.26368,
     +          26.55921,26.56489,26.58065,26.63554,27.36569,
     +          27.36569,27.36569,27.40948,27.40948,27.40999,
     +          27.41783,27.42901,27.42963,27.43824,27.43947,
     +          27.77609,27.77805,27.78169,27.78704,27.79372,
     +          28.01910,28.02209,28.02755,28.54429,30.17253,
     +          30.17448,30.17763,30.18179,30.18682,30.34387,
     +          30.34864,30.35188,30.41607,30.41652,30.41750,
     +          30.41894,30.42068/
        DATA NNIII/2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,3,3,3,3,3,3,3,
     +          3,4,3,3,3,3,3,3,4,4,3,3,3,3,4,4,4,4,3,3,3,3,3,3,
     +          3,3,3,3,3,5,3,3,3,3,3,3,3,5,5,3,3,5,5,5,5,6,6,6,
     +          6,6,6,4,4,4,3,3,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,
     +          4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     +          4,4,4,4,4,4,4,3,3,5,5,5,5/
        DATA GNIII/2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          6.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +          4.0D0,6.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,2.0D0,4.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,4.0D0,6.0D0,
     +          6.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,10.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +          6.0D0,6.0D0,4.0D0,2.0D0,6.0D0,
     +          8.0D0,4.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0,8.0D0,8.0D0,10.0D0,4.0D0,
     +          6.0D0,6.0D0,8.0D0,8.0D0,10.0D0,
     +          2.0D0,4.0D0,6.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          8.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +          4.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +          6.0D0,8.0D0,10.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,6.0D0,
     +          4.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +          6.0D0,8.0D0,10.0D0,6.0D0,8.0D0,
     +          6.0D0,8.0D0,10.0D0,12.0D0,8.0D0,
     +          10.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0,4.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0/
        DATA ENNIII/0.0D0,2.1635452D-02,7.180255,7.098413,7.108480,
     +          12.52548,12.52643,16.24252,18.08651,18.10019,
     +          23.16076,25.17799,25.18006,27.43827,28.56680,
     +          28.56730,30.45896,30.46342,33.13367,33.13441,
     +          35.65022,35.65797,35.67233,36.84229,36.85629,
     +          38.44641,38.32793,38.33453,38.39367,38.39807,
     +          38.40689,38.41771,38.64517,38.64825,38.95919,
     +          39.34056,39.34595,39.35325,39.39646,39.40031,
     +          39.71098,39.71098,39.79651,39.80747,40.55027,
     +          40.94474,40.94909,40.95552,40.96437,41.26192,
     +          41.26358,41.26631,41.26982,41.37555,41.47835,
     +          41.48166,41.68555,41.69232,41.69667,42.12335,
     +          42.13715,42.39634,42.39655,42.48893,42.49769,
     +          42.49625,42.49625,42.54757,42.54757,43.95493,
     +          43.95493,44.00932,44.00932,44.04135,44.04135,
     +          45.69180,45.69957,45.71402,46.28896,46.29317,
     +          46.46321,46.47039,46.71232,46.71811,46.72555,
     +          46.73671,46.81577,46.81788,46.85206,46.86286,
     +          46.92110,47.02857,47.03412,47.04068,47.61238,
     +          47.61238,47.61845,47.62763,47.75000,47.75000,
     +          47.77108,47.77108,49.01428,47.77802,47.88887,
     +          47.88887,47.88887,47.97657,47.97913,47.98245,
     +          47.98245,47.98363,47.98760,48.07270,48.08297,
     +          48.11119,48.11662,48.12305,48.13089,48.12993,
     +          48.14229,48.14024,48.14488,48.15087,48.15427,
     +          48.15307,48.16119,49.16950,49.17073,50.71214,
     +          50.71214,50.71214,50.71214/
        DATA NNIV/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,3,3,3,3,4,3,3,3,3,3,4,4,4,3,3,3,3,4,4,4,4,3,
     +          3,3,4,4,4,4,3,4,5,5,5,5,5,5,5,6,6,6,4,4,4,4,5,5,4/
        DATA GNIV/1.0D0,1.0D0,3.0D0,7.0D0,3.0D0,
     +          1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +          3.0D0,1.0D0,1.0D0,3.0D0,5.0D0,
     +          3.0D0,5.0D0,7.0D0,5.0D0,1.0D0,
     +          3.0D0,5.0D0,3.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,3.0D0,1.0D0,3.0D0,
     +          5.0D0,5.0D0,5.0D0,5.0D0,7.0D0,
     +          9.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,7.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,5.0D0,3.0D0,1.0D0,
     +          5.0D0,5.0D0,7.0D0,9.0D0,3.0D0,
     +          7.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +          7.0D0,9.0D0,11.0D0,3.0D0,5.0D0,
     +          7.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +          3.0D0,5.0D0,7.0D0/
        DATA ENNIV/0.0D0,8.323934,8.331770,8.349648,16.20427,
     +          21.75491,21.76399,21.77946,23.41898,29.18244,
     +          46.76804,50.15470,50.32483,50.32679,50.33118,
     +          52.06988,52.07031,52.07132,53.20933,57.68086,
     +          57.69048,57.71067,58.64906,59.62210,60.05779,
     +          60.05779,60.07403,60.44809,61.27855,61.27855,
     +          61.29070,61.78379,61.95650,61.97423,61.97423,
     +          61.97423,62.44215,62.44215,62.44215,62.67301,
     +          62.67685,62.68218,62.77282,62.86333,63.40415,
     +          63.40415,63.40415,63.41109,63.41767,63.41767,
     +          63.80760,64.05482,64.05569,64.05706,64.39976,
     +          64.70402,68.21900,68.53058,68.53058,68.53058,
     +          68.73986,68.73986,68.73986,71.28416,71.28416,
     +          71.28416,73.28070,73.60580,73.60580,73.61063,
     +          78.63129,78.63129,78.63129/
        DATA NNV/2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,
     +          6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8/
        DATA GNV/2.0D0,2.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +          8.0D0,10.0D0,12.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +          8.0D0,10.0D0,12.0D0,14.0D0,2.0D0,
     +          2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +          8.0D0,8.0D0,10.0D0,12.0D0,14.0D0,16.0D0/
        DATA ENNV/0.0D0,9.976473,10.00851,56.55396,59.23740,
     +          59.24660,60.05890,60.06188,75.17694,76.26962,
     +          76.26962,76.61120,76.61120,83.55153,84.09893,
     +          84.09893,84.27598,84.27598,88.02306,88.33514,
     +          88.33514,88.43854,88.43742,88.44214,88.44214,
     +          88.44313,88.44313,88.44313,90.68689,90.88043,
     +          90.88043,90.94527,90.94527,90.94912,90.94912,
     +          90.94974,90.94974,90.94974,90.94974,92.40136,
     +          92.53167,92.53167,92.57358,92.57358,92.57618,
     +          92.57618,92.57668,92.57668,92.57668,92.57668,
     +          92.57668/
        DATA NNVI/1,2,2,2,2,2,3,4/
        DATA GNVI/1.0D0,3.0D0,1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,3.0D0/
        DATA ENNVI/0.0D0,419.8009,426.2953,426.2965,426.3325,
     +          425.7398,497.9737,521.5830/
        DATA NOI/2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,3,3,3,3,3,3,3,3,3,
     +          4,4,4,4,4,4,3,3,3,5,5,3,4,4,4,4,4,4,4,4,5,5,5,6,
     +          6,5,5,5,5,5,5,5,5,6,6,6,7,7,6,6,6,6,6,6,6,6,8,8,
     +          7,7,7,7,7,7,7,7,9,9,8,8,8,8,8,8,8,8,10,10,9,9,9,9,
     +          9,9,9,9,11,11,10,10,10,10,10,10,10,10,3,3,3,3,3,3,
     +          3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,4,4,4,2,2,2,3,
     +          3,3,3,3,5,4,4,4,4,4,4,4,4,4,4,4,3,6,5,5,5,5,5,5,5,
     +          5,5,5,7,6,6,6,2/
        DATA GOI/5.0D0,3.0D0,1.0D0,5.0D0,1.0D0,
     +          5.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +          5.0D0,3.0D0,1.0D0,5.0D0,3.0D0,
     +          9.0D0,7.0D0,5.0D0,5.0D0,3.0D0,
     +          1.0D0,7.0D0,5.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,5.0D0,3.0D0,1.0D0,
     +          7.0D0,5.0D0,3.0D0,5.0D0,3.0D0,
     +          5.0D0,9.0D0,7.0D0,5.0D0,3.0D0,
     +          1.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +          3.0D0,1.0D0,5.0D0,3.0D0,9.0D0,
     +          7.0D0,5.0D0,3.0D0,1.0D0,7.0D0,
     +          5.0D0,3.0D0,5.0D0,3.0D0,1.0D0,
     +          5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +          5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +          5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +          5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +          5.0D0,3.0D0,9.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,5.0D0,3.0D0,
     +          7.0D0,5.0D0,3.0D0,9.0D0,7.0D0,
     +          5.0D0,5.0D0,3.0D0,1.0D0,7.0D0,
     +          3.0D0,5.0D0,5.0D0,5.0D0,3.0D0,
     +          1.0D0,9.0D0,7.0D0,5.0D0,9.0D0,
     +          11.0D0,9.0D0,7.0D0,7.0D0,7.0D0,
     +          5.0D0,3.0D0,5.0D0,3.0D0,1.0D0,
     +          7.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +          5.0D0,9.0D0,7.0D0,5.0D0,9.0D0,
     +          11.0D0,9.0D0,7.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,1.0D0,5.0D0,9.0D0,
     +          7.0D0,5.0D0,9.0D0,11.0D0,9.0D0,
     +          7.0D0,5.0D0,3.0D0,1.0D0,5.0D0,
     +          5.0D0,3.0D0,1.0D0,3.0D0/
        DATA ENOI/0.0D0,01.9651687D-02,2.8082693D-02,1.967363,0.4206081,
     +          9.146132,9.521420,10.74028,10.74053,10.74098,
     +          10.98893,10.98886,10.98895,11.83768,11.93056,
     +          12.07869,12.07870,12.07870,12.07872,12.07872,
     +          12.07872,12.08711,12.08711,12.08711,12.28604,
     +          12.28612,12.28627,12.35891,12.35891,12.35891,
     +          12.53927,12.54078,12.54176,12.66092,12.69755,
     +          12.72854,12.75377,12.75377,12.75377,12.75377,
     +          12.75377,12.75911,12.75911,12.75911,12.87829,
     +          12.87829,12.87829,13.02082,13.03891,13.06624,
     +          13.06624,13.06624,13.06624,13.06624,13.06913,
     +          13.06913,13.06913,13.13145,13.13145,13.13145,
     +          13.21004,13.22030,13.23559,13.23559,13.23559,
     +          13.23559,13.23559,13.23740,13.23740,13.23740,
     +          13.32166,13.32807,13.33749,13.33749,13.33749,
     +          13.33749,13.33749,13.33869,13.33869,13.33869,
     +          13.39308,13.39756,13.40353,385.3597,13.40353,
     +          13.40353,13.40353,13.40488,13.40488,13.40488,
     +          13.44262,13.44449,13.44872,13.44872,13.44872,
     +          13.44872,13.44872,13.44966,13.44966,13.44966,
     +          13.47577,13.47812,13.48112,13.48112,13.48112,
     +          13.48112,13.48112,13.48148,13.48148,13.48148,
     +          14.04685,14.04687,14.04730,14.09888,14.09975,
     +          14.10046,14.12320,14.12450,14.12526,14.13382,
     +          14.37218,14.46048,15.22525,15.28698,15.29424,
     +          15.29817,15.40062,15.40062,15.40062,15.40372,
     +          15.40390,15.40622,15.40550,15.41465,15.59420,
     +          15.59514,15.59577,15.65520,15.66431,15.66970,
     +          15.78109,15.78181,15.78222,15.82895,15.94391,
     +          16.01073,16.07676,16.07676,16.07676,16.07836,
     +          16.07844,16.08080,16.08005,16.08545,16.11433,
     +          16.11550,16.11614,16.23505,16.35702,16.35702,
     +          16.35702,16.35702,16.39057,16.39063,16.39308,
     +          16.39308,16.40451,16.40451,16.40451,16.54127,
     +          16.56668,16.56668,16.56668,23.53702/
        DATA NOII/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,2,3,3,3,3,3,3,3,3,
     +          3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,3,3,3,3,3,3,3,4,4,4,4,4,3,4,4,4,4,4,4,4,4,3,
     +          3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,
     +          4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,
     +          5,5,5,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     +          5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3,3,3,4,4,4,4,
     +          4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,5,5,3,3,3,3,3,4/
        DATA GOII/4.0D0,6.0D0,4.0D0,4.0D0,2.0D0,
     +          6.0D0,4.0D0,2.0D0,6.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          2.0D0,2.0D0,2.0D0,4.0D0,6.0D0,
     +          8.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,4.0D0,6.0D0,4.0D0,4.0D0,
     +          2.0D0,2.0D0,4.0D0,2.0D0,6.0D0,
     +          8.0D0,6.0D0,4.0D0,4.0D0,6.0D0,
     +          8.0D0,10.0D0,6.0D0,4.0D0,2.0D0,
     +          2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          8.0D0,6.0D0,8.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,2.0D0,4.0D0,8.0D0,6.0D0,
     +          10.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,8.0D0,10.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +          6.0D0,4.0D0,2.0D0,4.0D0,2.0D0,
     +          6.0D0,8.0D0,2.0D0,6.0D0,4.0D0,
     +          8.0D0,5.80D0,4.0D0,2.0D0,6.0D0,
     +          8.0D0,10.0D0,12.0D0,8.0D0,10.0D0,
     +          4.0D0,6.0D0,4.0D0,6.0D0,8.0D0,
     +          10.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,6.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,2.0D0,
     +          4.0D0,6.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,8.0D0,6.0D0,4.0D0,
     +          2.0D0,6.0D0,8.0D0,8.0D0,6.0D0,
     +          4.0D0,2.0D0,6.0D0,8.0D0,10.0D0,
     +          12.0D0,8.0D0,10.0D0,4.0D0,6.0D0,
     +          4.0D0,6.0D0,8.0D0,10.0D0,6.0D0,
     +          8.0D0,4.0D0,6.0D0,8.0D0,6.0D0,
     +          8.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          8.0D0,10.0D0,6.0D0,8.0D0,6.0D0,
     +          4.0D0,2.0D0,4.0D0,6.0D0,10.0D0,
     +          12.0D0,2.0D0,4.0D0,4.0D0,6.20D0,
     +          10.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0/
        DATA ENOII/0.0D0,3.323850,3.326454,5.017305,5.017491,
     +          14.85813,14.87838,14.88860,20.58005,20.57736,
     +          22.96648,22.97954,23.001876,23.41940,23.44172,
     +          24.26523,25.28586,25.63160,25.63849,25.64984,
     +          25.66529,25.66142,25.66154,25.83188,25.83760,
     +          25.84900,26.22564,26.24928,26.30498,26.35845,
     +          26.37943,26.55392,26.56133,28.12621,28.35835,
     +          28.36128,28.51330,28.51270,28.67733,28.68403,
     +          28.69369,28.70637,28.82200,28.83108,28.83932,
     +          28.82414,28.82992,28.85285,28.85711,28.85729,
     +          28.85808,28.86334,28.88355,28.94193,28.95606,
     +          29.06249,29.06893,29.58618,29.59923,29.61924,
     +          29.79726,29.82051,30.42546,30.47162,30.47763,
     +          30.48836,30.50400,30.74951,30.77135,30.80112,
     +          30.81214,31.02747,31.02747,31.14773,31.14812,
     +          31.31967,31.31982,31.37404,31.37430,31.46620,
     +          31.46649,31.55199,31.55199,31.55199,31.56553,
     +          31.61407,31.61407,31.61407,31.61407,31.61407,
     +          31.62925,31.63375,31.63644,31.63766,31.65117,
     +          31.65364,31.67396,31.69345,31.70178,31.71699,
     +          31.70200,31.71709,31.72948,31.72935,31.70999,
     +          31.71043,31.71889,31.73747,31.71911,31.73823,
     +          31.72081,31.72752,31.75062,31.75112,31.75553,
     +          31.75715,31.75586,31.75803,31.95026,31.96318,
     +          31.98375,32.03889,32.06284,32.14771,32.14780,
     +          32.35511,32.35511,32.36540,32.38251,32.39264,
     +          32.39264,32.40412,32.44667,32.46798,32.88345,
     +          32.88345,32.88345,32.88345,32.90963,32.91418,
     +          32.91418,32.92780,32.92780,32.93536,32.94354,
     +          32.95061,32.96264,32.93858,32.94181,32.95049,
     +          32.97082,32.95073,32.97146,32.96227,32.96227,
     +          32.97119,32.97528,32.97826,32.97999,32.97863,
     +          32.97999,33.19875,33.19968,33.20123,34.06365,
     +          34.06901,34.08607,34.08607,34.17174,34.17174,
     +          34.20029,34.20029,34.20504,34.20504,34.21390,
     +          34.21390,34.21960,34.22819,34.22819,34.23350,
     +          34.23350,34.25269,34.25269,34.48530,34.48530,
     +          36.19083,36.18759,36.19109,36.19123,36.19131,37.05294/
        DATA NOIII/2,2,2,2,2,2,2,3,3,2,2,2,2,2,2,3,3,3,3,2,2,2,3,3,
     +          3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          2,3,3,3,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,3,3,3,4,4,
     +          4,4,4,3,3,3,3,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,5,5,
     +          5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          5,5,5,5,5,5,5,5,5,3,3,3,6,6,6,6,7,3,3,4,4,4,3,4,
     +          4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,5/
        DATA GOIII/1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +          5.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +          3.0D0,1.0D0,5.0D0,3.0D0,3.0D0,
     +          1.0D0,3.0D0,5.0D0,3.0D0,5.0D0,
     +          3.0D0,1.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,3.0D0,5.0D0,1.0D0,3.0D0,
     +          5.0D0,5.0D0,1.0D0,5.0D0,7.0D0,
     +          9.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +          5.0D0,3.0D0,1.0D0,7.0D0,3.0D0,
     +          3.0D0,5.0D0,7.0D0,1.0D0,1.0D0,
     +          3.0D0,5.0D0,1.0D0,3.0D0,5.0D0,
     +          3.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +          7.0D0,9.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +          1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +          3.0D0,5.0D0,7.0D0,5.0D0,5.0D0,
     +          7.0D0,9.0D0,5.0D0,5.0D0,3.0D0,
     +          1.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          3.0D0,1.0D0,7.0D0,3.0D0,1.0D0,
     +          3.0D0,5.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,3.0D0,5.0D0,7.0D0,9.0D0,
     +          11.0D0,1.0D0,3.0D0,5.0D0,7.0D0,
     +          9.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +          3.0D0,1.0D0,5.0D0,7.0D0,9.0D0,
     +          5.0D0,7.0D0,9.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,7.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,5.0D0,3.0D0,5.0D0,
     +          7.0D0,7.0D0,7.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,3.0D0,3.0D0,1.0D0,
     +          3.0D0,5.0D0,7.0D0,9.0D0,3.0D0,
     +          5.0D0,7.0D0,3.0D0,5.0D0,7.0D0,
     +          7.0D0,5.0D0,3.0D0,5.0D0,7.0D0,
     +          9.0D0,3.0D0,5.0D0,7.0D0,1.0D0,
     +          3.0D0,5.0D0,3.0D0/
        DATA ENOIII/0.0D0,1.4059945D-02,3.8038719D-02,2.513308,5.354124,
     +          7.477820,14.88140,14.88477,14.88550,17.65325,
     +          17.65339,17.65514,23.19140,24.43587,26.09378,
     +          33.13600,33.15068,33.18253,33.85794,35.18196,
     +          35.20895,35.22094,36.07438,36.43500,36.45190,
     +          36.47919,36.89279,36.98353,37.22392,37.23410,
     +          37.25028,38.01204,38.90675,40.22861,40.25288,
     +          40.27497,40.26230,40.57149,40.57759,40.58673,
     +          40.84922,40.86335,40.87098,41.14086,41.25951,
     +          41.97723,41.99266,42.14902,42.56451,43.39812,
     +          43.41013,43.43237,44.22956,44.24270,44.27655,
     +          44.46952,45.03978,45.31862,45.32294,45.33144,
     +          45.34384,45.35962,45.34443,45.43903,45.45230,
     +          45.47797,45.62070,45.69189,45.69899,45.71153,
     +          45.91510,45.92614,45.93959,45.98626,46.25228,
     +          46.44183,45.21283,46.46955,46.62690,46.78899,
     +          46.78899,46.78899,46.82767,46.91713,46.91867,
     +          46.92080,47.01923,47.02679,47.03461,47.20199,
     +          47.20199,47.20199,47.21141,47.24910,48.62968,
     +          48.62968,48.62968,48.69874,48.86141,48.86587,
     +          48.87442,48.91428,48.91908,48.92621,48.93560,
     +          48.94701,49.36293,49.36248,49.36198,49.36323,
     +          49.37332,49.40500,49.41368,49.41845,49.63815,
     +          49.65178,49.65844,49.76514,49.77709,49.79367,
     +          49.78386,49.78386,49.78386,49.81572,49.78386,
     +          49.78386,49.78386,50.01249,50.03133,50.31391,
     +          50.31750,50.32357,51.41365,51.47638,51.47638,
     +          51.47638,52.44297,52.69355,52.85969,53.12613,
     +          53.14089,53.16110,53.31682,54.18348,54.33549,
     +          54.33549,54.34320,54.35460,54.36977,54.46407,
     +          54.47044,54.48261,54.88958,54.88958,54.88958,
     +          55.81414,55.82281,55.82951,56.14741,56.14741,
     +          56.14741,56.31095,56.31095,56.31095,56.73994,
     +          56.73994,56.73994,58.73808/
        DATA NOIV/2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,
     +          3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,4,4,3,3,3,3,3,3,5,3,3,3,3,5,5,5,5,3,4,4,4,3,
     +          3,4,4,6,6,4,4,3,3,3,3,3,3,3,4,4,7,7,4,4,4,4,4,4,
     +          4,4,4,4,4,4,4,4,4,4,3,8,8,4,4,3,3,3,3,3,3,3,3,3,
     +          3,3,3,3,3,5,5,3,3,5,5,3,3,5,5,5,5,5,5,5,5,5,5,5,
     +          3,3,3,3,3,3,3,3,3,6,6,6,6,4,4,3,4,4,7,7,7,7/
        DATA GOIV/2.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          6.0D0,4.0D0,2.0D0,2.0D0,6.0D0,
     +          4.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +          2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,4.0D0,6.0D0,
     +          2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          10.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          4.0D0,6.0D0,6.0D0,4.0D0,2.0D0,
     +          4.0D0,6.0D0,6.0D0,8.0D0,4.0D0,
     +          2.0D0,2.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +          6.0D0,8.0D0,2.0D0,2.0D0,4.0D0,
     +          6.0D0,6.0D0,8.0D0,2.0D0,4.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,4.0D0,
     +          6.0D0,2.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +          2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          6.0D0,4.0D0,2.0D0,4.0D0,6.0D0,
     +          6.0D0,8.0D0,4.0D0,6.0D0,6.0D0,
     +          8.0D0,2.0D0,6.0D0,8.0D0,4.0D0,
     +          2.0D0,4.0D0,6.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,2.0D0,4.0D0,6.0D0,
     +          6.0D0,4.0D0,4.0D0,6.0D0,8.0D0,
     +          2.0D0,4.0D0,6.0D0,8.0D0,4.0D0,
     +          6.0D0,6.0D0,4.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,6.0D0,4.0D0,2.0D0,
     +          6.0D0,8.0D0,4.0D0,2.0D0,6.0D0,
     +          4.0D0,2.0D0,4.0D0,6.0D0,6.0D0,
     +          8.0D0,4.0D0,2.0D0,2.0D0,4.0D0,
     +          6.0D0,8.0D0,4.0D0,6.0D0,2.0D0,
     +          4.0D0,6.0D0,2.0D0,4.0D0,6.0D0,8.0D0/
        DATA ENOIV/0.0D0,4.7920357D-02,8.824909,8.841201,8.864076,
     +          15.73825,15.73998,20.37910,22.37705,22.40721,
     +          28.67474,31.63571,31.63934,35.83378,35.83476,
     +          44.33902,48.37428,48.38508,54.37857,54.39532,
     +          54.42593,56.14158,56.17444,57.92984,57.94415,
     +          58.03452,58.04428,58.06108,58.08709,58.79609,
     +          59.33789,59.34961,59.36561,59.84372,59.87542,
     +          60.23497,61.10992,61.36131,61.37108,61.38501,
     +          61.40412,61.93150,61.93509,61.94088,61.94888,
     +          62.18008,62.18691,62.46812,62.48219,62.49133,
     +          63.30199,63.30286,63.32506,63.35387,63.75540,
     +          63.77412,64.30924,64.30999,66.87376,67.85857,
     +          67.86167,68.16618,68.17400,68.44416,68.44416,
     +          68.50069,68.50069,68.74507,70.50282,70.51955,
     +          70.55017,70.76975,70.76975,71.12993,71.15609,
     +          71.21387,71.21387,71.31690,71.33785,71.39315,
     +          71.39737,71.48887,71.50672,71.53300,72.12492,
     +          72.12764,72.47591,72.50269,72.88482,72.88482,
     +          73.16019,73.37047,73.37047,73.37047,73.37047,
     +          73.52322,73.52322,73.52322,73.60108,73.61112,
     +          73.64819,73.65725,73.68911,73.71453,73.93237,
     +          73.95444,74.05078,74.06293,74.06293,74.10930,
     +          74.12628,74.40265,74.40438,74.76035,74.76035,
     +          74.76035,74.76035,75.18896,75.18896,75.18896,
     +          76.30446,76.30806,76.44791,77.47625,77.47625,
     +          77.92433,77.92433,78.12258,78.12258,78.19797,
     +          78.21979,78.41159,78.43242,78.59385,78.59385,
     +          78.59385,78.59385,78.63718,78.63718,78.63718,
     +          78.85769,78.88398,78.91572,78.91572,78.96023,
     +          78.97250,78.98019,80.20107,80.20107,80.72665,
     +          80.72900,81.00314,81.01343,81.37509,81.37509,
     +          81.37509,81.37509,81.42716,81.42716,81.83012,
     +          82.88895,82.88895,83.03365,83.03365,83.03365,83.03365/
        DATA NOV/2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
     +          3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,
     +          4,4,4,4,4,4,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,4,6,
     +          6,6,6,6,4,4,4,6,4,4,4,4,4,7,7,7,7,7,8,8,8,8,5,5,
     +          5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6/
        DATA GOV/1.0D0,1.0D0,3.0D0,5.0D0,3.0D0,
     +          1.0D0,3.0D0,5.0D0,5.0D0,1.0D0,
     +          3.0D0,1.0D0,3.0D0,1.0D0,3.0D0,
     +          5.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +          3.0D0,5.0D0,7.0D0,3.0D0,1.0D0,
     +          3.0D0,5.0D0,5.0D0,5.0D0,3.0D0,
     +          5.0D0,7.0D0,1.0D0,5.0D0,3.0D0,
     +          1.0D0,7.0D0,3.0D0,3.0D0,1.0D0,
     +          1.0D0,3.0D0,5.0D0,3.0D0,3.0D0,
     +          5.0D0,7.0D0,5.0D0,7.0D0,3.0D0,
     +          3.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          3.0D0,3.0D0,3.0D0,5.0D0,7.0D0,
     +          3.0D0,1.0D0,3.0D0,5.0D0,5.0D0,
     +          5.0D0,3.0D0,7.0D0,3.0D0,5.0D0,
     +          7.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          5.0D0,3.0D0,1.0D0,7.0D0,3.0D0,
     +          3.0D0,3.0D0,5.0D0,7.0D0,5.0D0,
     +          3.0D0,3.0D0,5.0D0,7.0D0,3.0D0,
     +          3.0D0,5.0D0,7.0D0,1.0D0,3.0D0,
     +          5.0D0,5.0D0,5.0D0,3.0D0,5.0D0,
     +          7.0D0,7.0D0,3.0D0,3.0D0,5.0D0,
     +          7.0D0,1.0D0,3.0D0,5.0D0,5.0D0/
        DATA ENOV/0.0D0,10.18183,10.19878,10.23674,19.68863,
     +          26.48845,26.50776,26.54108,28.73015,35.69651,
     +          67.83862,69.59028,72.01395,72.28146,72.28596,
     +          72.29554,74.50599,74.50733,74.50979,75.95557,
     +          80.97483,80.99497,81.03748,82.38657,83.40436,
     +          83.97941,84.00407,84.04314,84.82139,85.49855,
     +          85.51269,85.53633,86.12596,86.43890,87.33036,
     +          87.33829,87.35107,87.73579,87.80076,87.81837,
     +          87.82866,88.39750,89.17985,89.60004,90.71603,
     +          91.26665,91.26665,91.26888,91.48672,92.04689,
     +          92.04763,92.04937,92.59937,92.97132,98.72523,
     +          99.49233,106.7049,106.7049,106.7049,100.2237,
     +          102.1987,102.8568,103.0377,103.0583,103.0944,
     +          103.1870,103.5465,103.5465,103.5676,103.8792,
     +          103.8829,104.1001,104.2509,104.2990,104.2990,
     +          104.2990,104.3064,104.3181,104.3333,104.4087,
     +          104.5556,104.5689,104.5754,105.0316,105.0733,
     +          106.7358,106.9963,106.9963,106.9963,106.9274,
     +          108.4187,108.5325,108.5325,108.5325,111.4108,
     +          111.5461,111.5461,111.5461,111.7535,111.7535,
     +          111.7535,111.8898,111.9082,112.1444,112.1444,
     +          112.1444,112.3809,115.9379,116.0435,116.0435,
     +          116.0435,116.1501,116.1501,116.1501,116.2166/
        DATA NOVI/2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,6,6,6,6,
     +          6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8/
        DATA GOVI/2.0D0,2.0D0,4.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,2.0D0,2.0D0,
     +          4.0D0,4.0D0,6.0D0,6.0D0,8.0D0,
     +          2.0D0,2.0D0,4.0D0,4.0D0,6.0D0,
     +          2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +          8.0D0,8.0D0,10.0D0,12.0D0,2.0D0,
     +          2.0D0,4.0D0,4.0D0,6.0D0,6.0D0,
     +          8.0D0,8.0D0,10.0D0,12.0D0,14.0D0,
     +          2.0D0,2.0D0,4.0D0,6.0D0,8.0D0,
     +          8.0D0,10.0D0,12.0D0,14.0D0,16.0D0,4.0D0,6.0D0/
        DATA ENOVI/0.0D0,11.94909,12.01505,79.35559,82.58831,
     +          82.60773,83.64374,83.65008,105.7219,107.0408,
     +          107.0487,107.4805,107.4831,107.5050,107.5062,
     +          117.6237,118.2920,118.2920,118.5122,118.5122,
     +          124.3735,124.3735,124.5034,124.5034,124.5142,
     +          124.5142,124.5156,124.5156,124.5156,127.8017,
     +          128.0311,128.0311,128.1171,128.1171,128.1243,
     +          128.1243,128.1252,128.1252,128.1252,128.1252,
     +          130.2520,130.3984,130.3984,130.4674,130.4674,
     +          130.4680,130.4680,130.4680,130.4680,130.4680,
     +          130.4693,130.4693/
        DATA NOVII/1,2,2,2,2,2,3,3,3,3,3,3,3,4,5,6/
        DATA GOVII/1.0D0,3.0D0,1.0D0,3.0D0,5.0D0,
     +          3.0D0,1.0D0,3.0D0,5.0D0,7.0D0,
     +          5.0D0,3.0D0,3.0D0,3.0D0,3.0D0,3.0D0/
        DATA ENOVII/0.0D0,561.0761,568.6182,568.6255,568.6938,
     +          573.9532,664.1129,664.1129,664.1129,665.1804,
     +          665.1804,665.1804,665.6218,697.8022,712.7239,720.8449/
        DATA SCI/4.179704,4.179868,4.180140,4.284864,4.411317, 
     +          4.556712,4.417538,4.418036,4.419087,4.460873,
     +          5.012059,5.012145,5.012123,4.656621,4.682271,
     +          4.682931,4.683973,4.715525,4.735114,4.735517, 
     +          4.736181,4.776610,4.823287,5.245868,4.960446, 
     +          4.636463,4.637096,4.638772,4.981131,4.981795, 
     +          4.983181,4.985224,4.985476,4.985839,4.648973, 
     +          4.987257,5.002644,5.026932,5.027268,5.027448, 
     +          4.751991,4.753114,4.754886,4.775015,4.807735, 
     +          4.820394,4.821310,4.822433,4.849118,4.880081, 
     +          4.964531,4.983084,4.983084,4.983084,4.988046, 
     +          4.988550,4.989272,4.739796,4.996660,5.002712, 
     +          5.007890,5.009317,5.009317,4.830776,4.830763, 
     +          4.830763,4.843919,4.885500,4.908801,4.963564, 
     +          4.983777,4.984663,4.984663,4.989659,4.989659, 
     +          4.992449,4.793379,4.998577,5.003253,5.004741, 
     +          5.006231,5.006231,4.972326,4.984214,4.985345, 
     +          4.985345,4.991673,4.991673,4.995098,4.831870, 
     +          5.000666,5.004451,5.005935,5.005935,5.004665, 
     +          4.984616,4.985763,4.985763,4.999064,4.999064, 
     +          4.999064,5.002865,5.004758,5.006233,5.006233, 
     +          5.006233,4.984145,4.984145,4.984432,5.002990, 
     +          5.002990,5.002990,5.005628,5.006508,5.006508, 
     +          5.006508,4.983365,4.983365,4.983365,5.006886, 
     +          5.006886,5.006886,5.008002,5.012074,5.012074, 
     +          5.012074,5.014099,5.014099,5.014099,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000/ 
        DATA SCII/3.322208,3.322644,3.633090,3.633257,3.633480,
     +          3.893420,3.893440,4.089184,4.229172,4.229597, 
     +          3.436720,3.692591,3.692789,4.589103,3.953149, 
     +          3.953178,4.702748,4.702820,3.603431,3.770060, 
     +          3.770254,4.440431,4.441056,4.442240,3.961645, 
     +          3.961659,4.991754,4.992091,3.992479,3.992479, 
     +          3.697579,3.795690,3.796068,4.770876,4.780956, 
     +          3.968260,3.968260,3.994325,3.994325,3.754887, 
     +          4.893884,4.894430,4.895358,4.896706,4.906213, 
     +          4.906944,3.971236,3.971236,3.996578,3.996578, 
     +          5.011171,5.086056,5.086787,5.087796,5.188515, 
     +          5.190205,5.591661,5.735427,5.737656,5.740739, 
     +          5.745144,5.937513,5.941305,5.947732,5.956571, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000,6.000000,6.000000,6.000000, 
     +          6.000000,6.000000/
        DATA SCIII/2.247678,2.511178,2.511298,2.511595,2.783334,
     +          2.988424,2.988601,2.988887,3.040348,3.275480,
     +          2.516355,2.624130,2.770250,2.779440,2.779510,
     +          2.779673,2.913505,2.912204,2.913576,3.001503,
     +          3.471912,3.472452,3.473566,2.656192,3.501991,
     +          2.707107,2.843331,2.843331,2.843447,3.667003,
     +          2.928022,2.928408,2.928967,2.942070,2.942230,
     +          2.942442,2.952919,2.960061,3.725865,3.726323,
     +          3.727023,2.996535,3.802981,3.848412,3.848814,
     +          3.849514,3.907525,3.915897,3.920174,3.920667,
     +          3.921389,4.006182,3.997117,4.006877,2.756046,
     +          4.057183,4.057739,4.058045,4.085242,2.876864,
     +          2.910818,2.910818,2.910818,4.166811,2.957773,
     +          2.957773,2.957773,2.998549,2.998549,2.998583,
     +          2.998549,3.003525,3.005696,3.020441,3.020510,
     +          3.020601,3.088544,2.797248,2.916939,2.968366,
     +          2.968366,2.968366,3.000603,3.000603,3.000640,
     +          3.003372,3.005378,3.009464,3.009464,3.009464,
     +          3.027199,2.830505,2.926040,2.976524,2.976524,
     +          2.976524,3.009360,2.932986,2.982088,2.982088,
     +          2.982088,2.986294,2.986294,2.986294,4.828427,
     +          4.828427,4.828427,5.151012,5.224114,5.224114,
     +          5.227788,5.497263,5.497263,5.502663,5.756044,
     +          5.817122,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000,6.000000,6.000000,6.000000,6.000000,
     +          6.000000/
        DATA SCIV/1.644934,1.923884,1.924364,1.778341,1.948965,
     +          1.949284,1.998202,1.998312,1.838941,1.962835,
     +          1.963075,1.999589,1.999669,2.001418,2.001418,
     +          1.874532,1.972046,1.972243,2.001394,2.001394,
     +          2.002842,2.002845,2.003133,2.003133,1.897882,
     +          1.978616,1.978616,2.003384,2.003384,2.004821,
     +          2.004821,2.005002,2.005002,2.005023,2.005023,
     +          1.913888,1.984033,1.984033,2.005113,2.005113,
     +          2.007061,2.007061,2.007217,2.007217,2.007234,
     +          2.007234,1.989961,1.989961,2.009654,2.009654,
     +          2.009800,2.009800,2.009800,2.009800,2.009800/
        DATASCV/0.6309066,0.7688928,0.9242349,0.9241881,0.9246764,
     +          1.026124,1.003322,1.003322,1.003322,1.020118,
     +          1.021842,1.027042,1.033988,1.051912,1.003001/
        DATA SNI/4.931870,5.108953,5.109031,5.204087,5.204087,
     +          5.329969,5.330800,5.331948,5.401419,5.403777,
     +          5.968683,5.969459,5.969806,5.605717,5.641185,
     +          5.641868,5.642995,5.644538,5.662621,5.663186,
     +          5.664362,5.702335,5.703489,5.705897,5.734946,
     +          5.736104,5.797976,5.797737,5.588642,5.591228,
     +          5.594790,5.615995,5.620492,5.980873,5.982464,
     +          5.982872,5.983638,5.985012,5.986994,5.988774,
     +          5.991692,5.990931,5.991989,5.990646,5.995945,
     +          5.996395,5.996926,5.997295,6.001428,6.002394,
     +          5.745163,5.761658,5.762813,5.764937,5.768167,
     +          5.774817,5.775745,5.778173,5.802792,5.696104,
     +          5.699982,5.706142,5.715098,5.722151,5.983984,
     +          5.985277,5.987724,5.991767,5.985780,5.989671,
     +          5.990610,5.994303,5.985830,5.987479,5.992273,
     +          5.996771,5.993288,5.995173,5.998955,6.002261,
     +          6.003886,6.255743,6.257061,6.360914,6.362586,
     +          5.757127,5.763044,5.772635,5.779661,5.791555,
     +          5.984846,5.986194,5.990250,5.997386,5.990476,
     +          5.992170,5.992623,6.000597,5.993189,5.993189,
     +          5.993189,6.000802,5.996590,5.998751,6.003087,
     +          6.005032,6.007153,5.793718,5.804320,5.818227,
     +          5.815205,5.821445,5.989322,5.989322,5.989322,
     +          5.989322,5.992901,5.992901,5.992901,6.003714,
     +          5.994695,5.997311,5.995185,6.008173,6.001116,
     +          6.010741,6.005529,6.005529,6.008008,5.801159,
     +          5.821038,5.833978,5.835981,5.835981,5.989852,
     +          5.989852,5.989852,5.989852,5.993398,5.993398,
     +          5.996287,5.996287,6.005342,6.014956,6.015613,
     +          6.015613,6.015613,5.971639,5.971639,5.850568,
     +          5.850568,5.850568,5.990061,5.990061,5.990061,
     +          5.990061,5.991796,5.991796,5.993244,5.993244,
     +          6.011374,6.018783,6.017594,6.017594,6.017594,
     +          5.858175,5.858175,5.863377,5.863377,5.863377,
     +          5.988659,5.988659,5.988659,5.988659,5.989390,
     +          5.989390,5.994151,5.994151,6.020563,6.027375,
     +          6.026996,6.026996,6.026996,5.866343,5.866343,
     +          5.874646,5.874646,5.874646,5.990859,5.990859,
     +          5.992667,5.992667,5.994933,5.994933,5.994933,
     +          5.994933,6.030020,6.030020,6.038991,6.038991,
     +          6.038991,5.873279,5.873279,5.877365,5.877365,
     +          5.877365,5.992039,5.992039,5.996431,5.996431,
     +          6.000838,6.000838,6.000838,6.000838,6.039686,
     +          6.039686,6.042562,6.042562,6.042562,6.018730,
     +          5.886326,5.994597,5.994597,6.047575,6.047575,
     +          6.047575,6.078406,6.078406/
        DATA SNII/4.048939,4.049242,4.049750,4.145151,4.258360,
     +          4.356432,4.688145,4.688257,4.688270,4.826217,
     +          4.826217,4.826272,5.142618,4.284333,4.284811,
     +          4.286872,4.288557,5.253337,4.532960,4.564950,
     +          4.565974,4.567596,5.379367,4.605227,4.634192,
     +          4.634804,4.635817,4.698180,4.771750,4.928997,
     +          4.930175,4.931792,4.940517,4.947424,4.947905,
     +          4.948512,4.976005,4.977055,4.977624,4.985716,
     +          5.001774,4.517682,4.519204,4.522714,5.167543,
     +          4.688999,4.706267,4.707887,4.710950,4.720615,
     +          4.721386,4.723658,4.732426,4.791679,5.359476,
     +          5.360877,5.362645,4.824183,4.939473,4.941747,
     +          4.944790,4.948400,4.959554,4.960499,4.962098,
     +          4.976268,4.977930,4.978850,4.988034,4.988136,
     +          4.988983,4.989127,4.990716,4.997378,4.997656,
     +          5.001124,4.999108,5.001566,5.001872,5.002052,
     +          5.004650,5.004791,5.510713,5.511550,5.512880,
     +          4.633615,4.635821,4.641956,4.663452,4.970949,
     +          4.970949,4.970949,4.990886,4.990886,4.991119,
     +          4.994713,4.999842,5.000126,5.004091,5.004655,
     +          5.899771,5.900361,5.901458,5.903069,5.905087,
     +          5.975470,5.976437,5.978201,6.162114,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000/
        DATA SNIII/3.264886,3.265738,3.559230,3.555733,3.556163,
     +          3.795859,3.795903,3.971288,4.062202,4.062887,
     +          4.328298,4.441761,4.441879,3.362787,4.644640,
     +          4.644671,3.648880,3.649321,3.924340,3.924418,
     +          4.208216,4.209135,4.210838,4.353292,4.355043,
     +          3.749472,4.546074,4.546964,4.554955,4.555551,
     +          4.556746,4.558211,3.785649,3.786211,4.632736,
     +          4.686664,4.687436,4.688480,3.926233,3.926969,
     +          3.987034,3.987034,4.752837,4.754452,4.866728,
     +          4.928826,4.929522,4.930549,4.931964,4.980142,
     +          4.980413,4.980860,4.981436,3.664742,5.015918,
     +          5.016470,5.050785,5.051935,5.052674,5.126587,
     +          5.129026,3.959079,3.959144,5.192320,5.193925,
     +          3.989434,3.989434,4.005147,4.005147,3.968565,
     +          3.968565,3.992409,3.992409,4.006539,4.006539,
     +          5.571515,5.574721,5.580696,6.132490,6.134099,
     +          5.935632,5.939608,6.083616,6.087341,6.092149,
     +          6.099411,6.364476,6.365573,6.178216,6.185985,
     +          6.229222,6.316157,6.320952,6.326658,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000,7.000000,7.000000,
     +          7.000000,7.000000,7.000000/
        DATA SNIV/2.226836,2.490623,2.490878,2.491461,2.755431,
     +          2.952339,2.952669,2.953231,3.013266,3.231892,
     +          2.493613,2.749589,2.762857,2.763010,2.763353,
     +          2.901417,2.901452,2.901533,2.994477,3.382731,
     +          3.383611,3.385459,3.472422,3.564919,3.607152,
     +          3.607152,3.608737,3.645438,3.728391,3.728391,
     +          2.639493,3.779903,3.797702,3.799535,3.799535,
     +          3.799535,2.797721,2.797721,2.797721,3.872625,
     +          3.873032,3.873596,3.883205,2.857107,2.934635,
     +          2.934635,2.934635,3.951730,3.952443,3.952443,
     +          2.993446,3.029914,3.030044,3.030246,4.061023,
     +          3.127314,2.880353,2.950477,2.950477,2.950477,
     +          2.998267,2.998267,2.998267,2.959708,2.959708,
     +          2.959708,4.785085,4.873190,4.873190,4.874527,
     +          7.000000,7.000000,7.000000/
        DATA SNV/1.634565,1.915400,1.916327,1.771111,1.943796,
     +          1.944398,1.997854,1.998051,1.833396,1.959357,
     +          1.959357,1.999384,1.999384,1.870467,1.969524,
     +          1.969524,2.001982,2.001982,1.895972,1.977561,
     +          1.977561,2.004889,2.004593,2.005843,2.005843,
     +          2.006106,2.006106,2.006106,1.914798,1.983841,
     +          1.983841,2.007186,2.007186,2.008574,2.008574,
     +          2.008797,2.008797,2.008797,2.008797,1.929893,
     +          1.990742,1.990742,2.010469,2.010469,2.011696,
     +          2.011696,2.011930,2.011930,2.011930,2.011930,
     +          2.011930/
        DATA SNVI/0.6290283,0.7657156,0.9208641,0.9208946,0.9217644,
     +          0.9074407,1.024314,1.024867/
        DATA SOI/5.998809,6.000254,6.000874,6.149045,6.029965,
     +          6.280362,6.354168,6.620857,6.620917,6.621027,
     +          6.681873,6.681856,6.681878,6.554275,6.592577,
     +          6.991942,6.991948,6.991948,6.991953,6.991953,
     +          6.991953,6.994710,6.994710,6.994710,6.749977,
     +          6.750016,6.750087,6.784759,6.784759,6.784759,
     +          7.156593,7.157186,7.157570,6.676266,6.701954,
     +          7.234456,6.993919,6.993919,6.993919,6.993919,
     +          6.993919,6.997045,6.997045,6.997045,6.836974,
     +          6.836974,6.836974,6.746833,6.766088,6.996467,
     +          6.996467,6.996467,6.996467,6.996467,6.999115,
     +          6.999115,6.999115,6.869720,6.869720,6.869720,
     +          6.793480,6.808908,6.999084,6.999084,6.999084,
     +          6.999084,6.999084,7.001479,7.001479,7.001479,
     +          6.826997,6.839929,7.001803,7.001803,7.001803,
     +          7.001803,7.001803,7.003955,7.003955,7.003955,
     +          6.852827,6.864539,7.004705,8.000000,7.004705,
     +          7.004705,7.004705,7.007905,7.007905,7.007905,
     +          6.877349,6.883498,7.007761,7.007761,7.007761,
     +          7.007761,7.007761,7.010591,7.010591,7.010591,
     +          6.890946,6.900387,7.011453,7.011453,7.011453,
     +          7.011453,7.011453,7.012788,7.012788,7.012788,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000/
        DATA SOII/4.784610,4.940430,4.940555,5.022952,5.022962,
     +          5.557054,5.558274,5.558889,5.930025,5.929834,
     +          5.160761,5.162283,5.164577,5.214052,5.216705,
     +          6.210938,5.445367,5.490556,5.491464,5.492962,
     +          5.495002,5.494491,5.494507,5.517109,5.517870,
     +          5.519391,5.570158,5.573380,5.580988,6.392210,
     +          6.394130,5.615287,5.616316,5.844496,5.880436,
     +          5.880893,5.904250,5.904674,5.930839,5.931911,
     +          5.933457,5.935489,5.954107,5.955576,5.956912,
     +          5.954453,5.955388,5.959104,5.959794,5.959825,
     +          5.959952,5.960805,5.964087,5.973599,5.975908,
     +          5.993385,5.994448,5.442262,5.445265,5.449878,
     +          5.491284,5.496742,6.232405,5.654758,5.656267,
     +          5.658962,5.662895,5.725537,5.731195,5.738926,
     +          5.741796,6.348958,6.348958,6.373241,6.373322,
     +          6.408604,6.408635,6.419951,6.420006,6.439373,
     +          6.439435,5.943564,5.943564,5.943564,5.947441,
     +          5.961402,5.961402,5.961402,5.961402,6.471051,
     +          5.965786,5.967088,5.967867,5.968222,5.972136,
     +          5.972852,5.978758,6.488330,5.986873,5.991322,
     +          5.986938,5.991353,5.994984,5.994948,5.989273,
     +          5.989404,5.991879,5.997332,5.991945,5.997554,
     +          5.992443,5.994410,6.001196,6.001346,6.002642,
     +          6.003121,6.002740,6.003380,5.576062,5.580966,
     +          5.588796,5.609913,5.619139,6.121709,6.121740,
     +          5.734796,5.734796,5.738977,5.745944,5.750079,
     +          5.750079,5.754774,5.772264,5.781077,5.960446,
     +          5.960446,5.960446,5.960446,5.972282,5.974347,
     +          5.974347,5.980535,5.980535,5.983981,5.987715,
     +          5.990945,5.996456,5.985451,5.986923,5.990890,
     +          6.000214,5.991003,6.000511,5.996286,5.996286,
     +          6.000386,6.002267,6.003636,6.004436,6.003808,
     +          6.004436,6.864734,6.865004,6.865458,6.871479,
     +          6.874276,6.883227,6.883227,6.929313,6.929313,
     +          6.945119,6.945119,6.947771,6.947771,7.214550,
     +          7.214550,6.955942,6.960794,6.960794,6.963802,
     +          6.963802,6.974759,6.974759,6.897856,6.897856,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000/
        DATA SOIII/3.980091,3.980605,3.981483,4.073126,4.181011,
     +          4.263698,4.567496,2.851461,2.851508,4.688399,
     +          4.688406,4.688483,4.944256,5.004756,5.087306,
     +          4.201648,4.202927,4.205704,4.265077,5.589531,
     +          5.591178,5.591910,4.466920,4.500863,4.502461,
     +          4.505044,4.544430,5.702087,4.576288,4.577272,
     +          4.578837,4.653335,4.743010,4.880211,4.882788,
     +          4.885133,4.883788,4.916797,4.917453,4.918434,
     +          4.946753,4.948285,4.949112,4.978528,4.991553,
     +          5.071567,5.073311,5.091046,6.092469,5.236801,
     +          5.238239,5.240906,4.450987,4.453166,4.458786,
     +          4.490992,5.440957,5.477274,5.477840,5.478956,
     +          5.480584,5.482658,4.640882,4.657492,4.659830,
     +          4.664354,4.689624,5.526725,5.527675,5.529354,
     +          4.742366,4.744359,4.746791,4.755242,4.803841,
     +          5.629193,5.463434,5.633066,5.655169,4.904211,
     +          4.904211,4.904211,4.911571,5.696494,5.696715,
     +          5.697021,4.948280,4.949739,4.951245,4.983719,
     +          4.983719,4.983719,4.985557,4.992922,4.595488,
     +          4.595488,4.595488,4.614187,5.995187,5.995924,
     +          5.997336,6.003933,6.004729,6.005913,6.007472,
     +          6.009368,6.079757,6.079680,6.079593,6.079808,
     +          6.081549,6.087021,6.088523,6.089349,6.127790,
     +          6.130199,6.131379,6.150372,6.152512,6.155484,
     +          4.922875,4.922875,4.922875,4.932409,4.922875,
     +          4.922875,4.922875,4.991952,4.997716,6.251313,
     +          6.251994,6.253141,4.947119,4.974443,4.974443,
     +          4.974443,5.003925,6.782259,6.828280,6.541492,
     +          6.547457,6.555665,6.965416,7.060263,7.160808,
     +          7.160808,7.166230,7.174317,7.185196,7.256399,
     +          7.261456,7.271212,7.771374,7.771374,7.771374,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000/
        DATA SOIV/3.228562,3.230039,3.508826,3.509360,3.510109,
     +          3.741247,3.741307,3.904661,3.977057,3.978160,
     +          4.214302,4.331145,4.331291,4.503491,4.503533,
     +          3.322590,3.617383,3.618198,4.097019,4.098440,
     +          4.101037,4.249484,4.252384,4.410741,4.412061,
     +          4.420406,4.421309,4.422863,4.425270,4.491520,
     +          4.543003,4.544125,4.545658,4.591770,4.594849,
     +          3.506632,4.717019,4.742457,4.743450,4.744866,
     +          4.746809,4.800908,4.801279,4.801878,4.802707,
     +          4.826727,4.827440,4.856910,4.858391,4.859354,
     +          3.927959,3.928085,4.948471,4.951597,4.995503,
     +          4.997566,5.057139,5.057223,3.602068,5.487784,
     +          5.488193,5.528638,5.529685,3.943577,3.943577,
     +          3.956409,3.956409,5.607411,5.152442,5.155901,
     +          5.162243,5.906104,5.906104,5.285100,5.290775,
     +          3.955026,3.955026,5.325924,5.330538,6.007065,
     +          6.007765,6.023023,6.026014,6.030425,6.132527,
     +          6.133009,5.594399,5.600957,3.969004,3.969004,
     +          5.768017,5.824149,5.824149,5.824149,5.824149,
     +          5.865850,5.865850,5.865850,5.887424,5.890223,
     +          5.900586,5.903125,5.912084,5.919259,5.981793,
     +          5.988238,6.512458,4.040435,4.040435,6.034046,
     +          6.039135,6.592915,6.593322,6.679720,6.679720,
     +          6.679720,6.679720,6.791924,6.791924,6.791924,
     +          7.150804,7.152208,7.208681,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000,
     +          8.000000/
        DATA SOV/2.212299,2.477108,2.477559,2.478570,2.736373,
     +          2.929941,2.930501,2.931468,2.995396,3.204502,
     +          2.480139,2.586177,2.736415,2.753261,2.753545,
     +          2.754149,2.895501,2.895588,2.895747,2.990361,
     +          3.333698,3.335127,3.338143,3.434917,3.509306,
     +          3.551885,3.553720,3.556629,3.614975,3.666382,
     +          3.667461,3.669268,3.714561,3.738797,3.808601,
     +          3.809227,3.810236,3.840735,3.845908,3.847311,
     +          3.848131,3.893723,3.957266,2.655746,2.780047,
     +          2.842480,2.842480,2.842734,2.867645,2.932266,
     +          2.932352,2.932553,2.996815,3.040747,2.722722,
     +          2.858081,4.369750,4.369750,4.369750,2.990545,
     +          4.293693,4.399677,4.429359,4.432752,4.438707,
     +          4.454041,4.514207,4.514207,4.517767,4.570812,
     +          4.571450,2.913395,2.952782,2.965416,2.965416,
     +          2.965416,4.644915,4.646958,4.649636,2.994349,
     +          4.688903,4.691261,4.692408,4.774586,4.782194,
     +          2.928606,3.022013,3.022013,3.022013,2.997125,
     +          2.933282,2.986425,2.986425,2.986425,5.872365,
     +          5.931636,5.931636,5.931636,6.025977,6.025977,
     +          6.025977,6.090481,6.099398,6.217293,6.217293,
     +          6.217293,6.343698,8.000000,8.000000,8.000000,
     +          8.000000,8.000000,8.000000,8.000000,8.000000/
        DATA SOVI/1.626749,1.908750,1.910343,1.765577,1.939605,
     +          1.940666,1.997515,1.997865,1.829541,1.956604,
     +          1.957376,1.999562,1.999822,2.001964,2.002083,
     +          1.867336,1.968342,1.968342,2.001995,2.001995,
     +          1.976059,1.976059,2.004681,2.004681,2.007063,
     +          2.007063,2.007365,2.007365,2.007365,1.914098,
     +          1.982390,1.982390,2.008208,2.008208,2.010369,
     +          2.010369,2.010631,2.010631,2.010631,2.010631,
     +          1.930105,1.987143,1.987143,2.014184,2.014184,
     +          2.014431,2.014431,2.014431,2.014431,2.014431,
     +          2.014965,2.014965/
        DATA SOVII/0.6273875,0.7631111,0.9180546,0.9182081,0.9196253,
     +          1.029738,0.9543552,0.9543552,0.9543552,1.004676,
     +          1.004676,1.004676,1.025589,1.027917,1.034432,
     +          1.045346/
*
*       Find index for atom and ion, 10*IAT+IZI
*
c       IF(IAT.EQ.26.AND.IZI.GE.6.AND.IZI.LE.9) GO TO 260
        IF(IAT.GT.2.AND.IAT.LT.6)GO TO 9999
        IF(IAT.LT.1.OR.IAT.GT.8)GO TO 9999
        IND=10*IAT+IZI
        IF(IND.EQ.11) GO TO 11
        IF(IND.EQ.21) GO TO 21
        IF(IND.EQ.22) GO TO 22
        IF(IND.EQ.61) GO TO 61
        IF(IND.EQ.61) GO TO 62
        IF(IND.EQ.63) GO TO 63
        IF(IND.EQ.64) GO TO 64
        IF(IND.EQ.65) GO TO 65
        IF(IND.EQ.66) GO TO 66
        IF(IND.EQ.71) GO TO 71
        IF(IND.EQ.72) GO TO 72
        IF(IND.EQ.73) GO TO 73
        IF(IND.EQ.74) GO TO 74
        IF(IND.EQ.75) GO TO 75
        IF(IND.EQ.76) GO TO 76
        IF(IND.EQ.77) GO TO 77
        IF(IND.EQ.81) GO TO 81
        IF(IND.EQ.82) GO TO 82
        IF(IND.EQ.83) GO TO 83
        IF(IND.EQ.84) GO TO 84
        IF(IND.EQ.85) GO TO 85
        IF(IND.EQ.86) GO TO 86
        IF(IND.EQ.87) GO TO 87
        IF(IND.EQ.88) GO TO 88
* 
*       CALCULATING PARTITION FUNCTIONS FOR HYDROGEN
* 
 11     CALL PARTDV(T,ANE,ZH,MH,NHYD,GHYD,ENHYD,SHYD,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR HEI
* 
 21     CALL PARTDV(T,ANE,ZHE,MHEI,NHEL,GHEL,ENHEL,SHEL,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR HEII
* 
 22     CALL PARTDV(T,ANE,ZHE,MHEII,NHYD,GHYD,ENHYD,SHYD,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR CI
* 
 61     CALL PARTDV(T,ANE,ZC,MCI,NCI,GCI,ENCI,SCI,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR CII
* 
 62     CALL PARTDV(T,ANE,ZC,MCII,NCII,GCII,ENCII,SCII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR CIII
* 
 63     CALL PARTDV(T,ANE,ZC,MCIII,NCIII,GCIII,ENCIII,SCIII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR CIV
* 
 64     CALL PARTDV(T,ANE,ZC,MCIV,NCIV,GCIV,ENCIV,SCIV,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR CV
* 
 65     CALL PARTDV(T,ANE,ZC,MCV,NCV,GCV,ENCV,SCV,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR CVI
* 
 66     CALL PARTDV(T,ANE,ZC,MH,NHYD,GHYD,ENHYD,SHYD,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NI
* 
 71     CALL PARTDV(T,ANE,ZN,MNI,NNI,GNI,ENNI,SNI,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NII
* 
 72     CALL PARTDV(T,ANE,ZN,MNII,NNII,GNII,ENNII,SNII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NIII
* 
 73     CALL PARTDV(T,ANE,ZN,MNIII,NNIII,GNIII,ENNIII,SNIII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NIV
* 
 74     CALL PARTDV(T,ANE,ZN,MNIV,NNIV,GNIV,ENNIV,SNIV,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NV
* 
 75     CALL PARTDV(T,ANE,ZN,MNV,NNV,GNV,ENNV,SNV,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NVI
* 
 76     CALL PARTDV(T,ANE,ZN,MNVI,NNVI,GNVI,ENNVI,SNVI,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR NVII
* 
 77     CALL PARTDV(T,ANE,ZN,MH,NHYD,GHYD,ENHYD,SHYD,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OI
* 
 81     CALL PARTDV(T,ANE,ZO,MOI,NOI,GOI,ENOI,SOI,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OII
* 
 82     CALL PARTDV(T,ANE,ZO,MOII,NOII,GOII,ENOII,SOII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OIII
* 
 83     CALL PARTDV(T,ANE,ZO,MOIII,NOIII,GOIII,ENOIII,SOIII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OIV
* 
 84     CALL PARTDV(T,ANE,ZO,MOIV,NOIV,GOIV,ENOIV,SOIV,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OV
* 
 85     CALL PARTDV(T,ANE,ZO,MOV,NOV,GOV,ENOV,SOV,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OVI
* 
 86     CALL PARTDV(T,ANE,ZO,MOVI,NOVI,GOVI,ENOVI,SOVI,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OVII
* 
 87     CALL PARTDV(T,ANE,ZO,MOVII,NOVII,GOVII,ENOVII,SOVII,U)
        GO TO 8888
* 
*       CALCULATING PARTITION FUNCTIONS FOR OVIII
* 
 88     CALL PARTDV(T,ANE,ZO,MH,NHYD,GHYD,ENHYD,SHYD,U)
        GO TO 8888
C
C 
C       CALCULATING PARTITION FUNCTIONS FOR FE VI - FE IX
C 
C260    CALL PFFE(IZI,T,ANE,U)
 8888   CONTINUE
        RETURN
 9999   U=0
        WRITE(*,*)!! INVALID ATOM IN USER SUPPLIED ROUTINE PARTFUN !!
        STOP
        END
C
C     **************************************************************
C
C 
        SUBROUTINE PARTDV(TEMP,DNE,Z,NLEV,NE,GEE,ENRGY,S,U)
C       ===================================================
C 
        INCLUDE 'PARAMS.FOR'
        DIMENSION GEE(*),ENRGY(*),S(*)
        INTEGER NE(*)
        U=0.0
        ET=TEMP/11604.8
        P=(14.69D0-0.20-0.6667*LOG10(DNE))
C 
        DO 10 I=1,NLEV
           U1=FLOAT(NE(I))
           ZSTAR=Z-S(I)
           IF (ZSTAR.GT.0)THEN
               W=P+4.*LOG10(ZSTAR)-4.*LOG10(U1)
            ELSE
               W=0.0
           ENDIF
           IF (W.GT.1.) W=1.
C 
           IF ((ENRGY(I)/ET).LT.65.0) THEN
                U1=GEE(I)*W*EXP(-ENRGY(I)/ET)
            ELSE
                U1=0.0
           ENDIF
           U=U+U1
 10     CONTINUE
        RETURN
        END
C
C     **************************************************************
C
      subroutine pfni(ion,t,pf,dut,dun)
c     =================================
c
c     partition functions for Ni IV to Ni IX
c
c     this routine interpolates within a grid
c     calculated from all levels predicted by
c     Kurucz (1992), i.e. over 12,000 levels per ion.
c     the partition functions depend only on T !
c     (i.e. no level dissolution with increasing density)
c     TL  27-DEC-1994, 23-JAN-1995
c
c     Output:  PF   partition function
c              DUT  d(PF)/dT
c              DUN  d(PF)/d(ANE) (=0 in this version)
c
      implicit double precision (a-h,o-z)
c
      dimension g0(6)
      dimension p4a(190),p4b(170)
      dimension p5a(190),p5b(170)
      dimension p6a(190),p6b(170)
      dimension p7a(190),p7b(170)
      dimension p8a(190),p8b(170)
      dimension p9a(190),p9b(170)
      parameter (xen=2.302585093,xmil=0.001)
c
      data g0/28.,25.,6.,25.,28.,21./
c
      data p4a/
     .    1.447,1.464,1.482,1.501,1.518,1.535,1.551,1.567,1.582,1.596,
     .    1.610,1.623,1.636,1.648,1.659,1.671,1.681,1.692,1.702,1.711,
     .    1.721,1.730,1.739,1.748,1.757,1.765,1.774,1.782,1.791,1.799,
     .    1.808,1.816,1.824,1.833,1.841,1.850,1.859,1.868,1.877,1.886,
     .    1.895,1.905,1.914,1.924,1.934,1.945,1.955,1.966,1.977,1.989,
     .    2.000,2.012,2.025,2.037,2.050,2.063,2.077,2.091,2.105,2.119,
     .    2.134,2.149,2.164,2.179,2.195,2.211,2.227,2.243,2.260,2.276,
     .    2.293,2.310,2.327,2.344,2.362,2.379,2.397,2.414,2.432,2.449,
     .    2.467,2.484,2.502,2.519,2.537,2.554,2.571,2.588,2.606,2.623,
     .    2.640,2.657,2.674,2.690,2.707,2.723,2.740,2.756,2.772,2.788,
     .    2.804,2.819,2.835,2.850,2.866,2.881,2.896,2.911,2.925,2.940,
     .    2.954,2.969,2.983,2.997,3.010,3.024,3.038,3.051,3.064,3.077,
     .    3.090,3.103,3.116,3.128,3.141,3.153,3.165,3.177,3.189,3.201,
     .    3.213,3.224,3.235,3.247,3.258,3.269,3.280,3.291,3.301,3.312,
     .    3.322,3.332,3.343,3.353,3.363,3.373,3.382,3.392,3.402,3.411,
     .    3.421,3.430,3.439,3.448,3.457,3.466,3.475,3.484,3.492,3.501,
     .    3.509,3.518,3.526,3.534,3.542,3.550,3.558,3.566,3.574,3.582,
     .    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,
     .    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722/
      data p4b/
     .    3.589,3.597,3.604,3.612,3.619,3.626,3.634,3.641,3.648,3.655,
     .    3.662,3.669,3.676,3.682,3.689,3.696,3.702,3.709,3.715,3.722,
     .    3.728,3.734,3.740,3.747,3.753,3.759,3.765,3.771,3.777,3.782,
     .    3.788,3.794,3.800,3.805,3.811,3.816,3.822,3.827,3.833,3.838,
     .    3.843,3.849,3.854,3.859,3.864,3.869,3.874,3.879,3.884,3.889,
     .    3.894,3.899,3.904,3.909,3.913,3.918,3.923,3.927,3.932,3.936,
     .    3.941,3.945,3.950,3.954,3.959,3.963,3.967,3.972,3.976,3.980,
     .    3.984,3.988,3.993,3.997,4.001,4.005,4.009,4.013,4.017,4.021,
     .    4.024,4.028,4.032,4.036,4.040,4.043,4.047,4.051,4.055,4.058,
     .    4.062,4.065,4.069,4.072,4.076,4.079,4.083,4.086,4.090,4.093,
     .    4.097,4.100,4.103,4.107,4.110,4.113,4.116,4.120,4.123,4.126,
     .    4.129,4.132,4.135,4.138,4.141,4.144,4.148,4.151,4.154,4.157,
     .    4.159,4.162,4.165,4.168,4.171,4.174,4.177,4.180,4.182,4.185,
     .    4.188,4.191,4.193,4.196,4.199,4.202,4.204,4.207,4.210,4.212,
     .    4.215,4.217,4.220,4.223,4.225,4.228,4.230,4.233,4.235,4.238,
     .    4.240,4.243,4.245,4.247,4.250,4.252,4.255,4.257,4.259,4.262,
     .    4.264,4.266,4.268,4.271,4.273,4.275,4.278,4.280,4.282,4.284/
      data p5a/
     .    1.398,1.408,1.427,1.446,1.466,1.486,1.506,1.526,1.545,1.564,
     .    1.583,1.601,1.619,1.636,1.652,1.668,1.683,1.698,1.712,1.725,
     .    1.738,1.751,1.763,1.775,1.786,1.797,1.808,1.818,1.828,1.837,
     .    1.846,1.855,1.864,1.873,1.881,1.889,1.897,1.904,1.912,1.919,
     .    1.926,1.933,1.940,1.946,1.953,1.960,1.966,1.972,1.979,1.985,
     .    1.991,1.997,2.003,2.009,2.016,2.022,2.028,2.034,2.040,2.046,
     .    2.052,2.058,2.065,2.071,2.077,2.084,2.090,2.097,2.103,2.110,
     .    2.117,2.124,2.131,2.138,2.145,2.152,2.160,2.167,2.175,2.183,
     .    2.191,2.199,2.207,2.216,2.224,2.233,2.241,2.250,2.259,2.268,
     .    2.278,2.287,2.297,2.306,2.316,2.326,2.336,2.346,2.356,2.367,
     .    2.377,2.387,2.398,2.409,2.419,2.430,2.441,2.452,2.463,2.474,
     .    2.485,2.497,2.508,2.519,2.530,2.542,2.553,2.564,2.576,2.587,
     .    2.599,2.610,2.621,2.633,2.644,2.655,2.667,2.678,2.689,2.701,
     .    2.712,2.723,2.734,2.745,2.757,2.768,2.779,2.790,2.801,2.812,
     .    2.822,2.833,2.844,2.855,2.865,2.876,2.886,2.897,2.907,2.918,
     .    2.928,2.938,2.948,2.958,2.968,2.978,2.988,2.998,3.008,3.018,
     .    3.027,3.037,3.046,3.056,3.065,3.075,3.084,3.093,3.102,3.111,
     .    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,
     .    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279/
      data p5b/
     .    3.120,3.129,3.138,3.147,3.156,3.164,3.173,3.182,3.190,3.198,
     .    3.207,3.215,3.223,3.232,3.240,3.248,3.256,3.264,3.272,3.279,
     .    3.287,3.295,3.303,3.310,3.318,3.325,3.333,3.340,3.347,3.355,
     .    3.362,3.369,3.376,3.383,3.390,3.397,3.404,3.411,3.417,3.424,
     .    3.431,3.438,3.444,3.451,3.457,3.464,3.470,3.476,3.483,3.489,
     .    3.495,3.501,3.507,3.514,3.520,3.526,3.531,3.537,3.543,3.549,
     .    3.555,3.561,3.566,3.572,3.578,3.583,3.589,3.594,3.600,3.605,
     .    3.610,3.616,3.621,3.626,3.632,3.637,3.642,3.647,3.652,3.657,
     .    3.662,3.667,3.672,3.677,3.682,3.687,3.692,3.697,3.701,3.706,
     .    3.711,3.716,3.720,3.725,3.729,3.734,3.738,3.743,3.747,3.752,
     .    3.756,3.761,3.765,3.769,3.774,3.778,3.782,3.786,3.790,3.795,
     .    3.799,3.803,3.807,3.811,3.815,3.819,3.823,3.827,3.831,3.835,
     .    3.839,3.843,3.846,3.850,3.854,3.858,3.862,3.865,3.869,3.873,
     .    3.876,3.880,3.884,3.887,3.891,3.894,3.898,3.901,3.905,3.908,
     .    3.912,3.915,3.918,3.922,3.925,3.929,3.932,3.935,3.939,3.942,
     .    3.945,3.948,3.951,3.955,3.958,3.961,3.964,3.967,3.970,3.974,
     .    3.977,3.980,3.983,3.986,3.989,3.992,3.995,3.998,4.001,4.004/
      data p6a/
     .    0.778,0.804,0.817,0.834,0.854,0.876,0.901,0.928,0.957,0.987,
     .    1.017,1.048,1.079,1.109,1.139,1.169,1.197,1.225,1.253,1.279,
     .    1.304,1.329,1.353,1.376,1.398,1.419,1.440,1.459,1.478,1.497,
     .    1.515,1.532,1.548,1.564,1.580,1.594,1.609,1.623,1.636,1.649,
     .    1.662,1.674,1.686,1.698,1.709,1.720,1.730,1.740,1.750,1.760,
     .    1.769,1.779,1.788,1.796,1.805,1.813,1.821,1.829,1.837,1.845,
     .    1.852,1.860,1.867,1.874,1.881,1.888,1.894,1.901,1.907,1.914,
     .    1.920,1.926,1.932,1.938,1.944,1.950,1.956,1.962,1.968,1.974,
     .    1.979,1.985,1.991,1.996,2.002,2.007,2.013,2.018,2.024,2.029,
     .    2.035,2.041,2.046,2.052,2.057,2.063,2.068,2.074,2.080,2.086,
     .    2.091,2.097,2.103,2.109,2.115,2.121,2.127,2.133,2.139,2.145,
     .    2.152,2.158,2.164,2.171,2.177,2.184,2.190,2.197,2.204,2.211,
     .    2.218,2.225,2.232,2.239,2.246,2.253,2.261,2.268,2.276,2.283,
     .    2.291,2.298,2.306,2.314,2.322,2.330,2.338,2.346,2.354,2.362,
     .    2.370,2.379,2.387,2.395,2.404,2.412,2.420,2.429,2.438,2.446,
     .    2.455,2.463,2.472,2.481,2.489,2.498,2.507,2.516,2.524,2.533,
     .    2.542,2.551,2.560,2.569,2.577,2.586,2.595,2.604,2.613,2.622,
     .    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,
     .    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796/
      data p6b/
     .    2.631,2.639,2.648,2.657,2.666,2.675,2.683,2.692,2.701,2.710,
     .    2.718,2.727,2.736,2.744,2.753,2.761,2.770,2.779,2.787,2.796,
     .    2.804,2.812,2.821,2.829,2.838,2.846,2.854,2.862,2.871,2.879,
     .    2.887,2.895,2.903,2.911,2.919,2.927,2.935,2.943,2.951,2.958,
     .    2.966,2.974,2.982,2.989,2.997,3.005,3.012,3.020,3.027,3.035,
     .    3.042,3.049,3.057,3.064,3.071,3.078,3.086,3.093,3.100,3.107,
     .    3.114,3.121,3.128,3.135,3.141,3.148,3.155,3.162,3.169,3.175,
     .    3.182,3.188,3.195,3.202,3.208,3.214,3.221,3.227,3.234,3.240,
     .    3.246,3.252,3.259,3.265,3.271,3.277,3.283,3.289,3.295,3.301,
     .    3.307,3.313,3.319,3.325,3.330,3.336,3.342,3.348,3.353,3.359,
     .    3.364,3.370,3.376,3.381,3.386,3.392,3.397,3.403,3.408,3.413,
     .    3.419,3.424,3.429,3.434,3.440,3.445,3.450,3.455,3.460,3.465,
     .    3.470,3.475,3.480,3.485,3.490,3.495,3.499,3.504,3.509,3.514,
     .    3.518,3.523,3.528,3.533,3.537,3.542,3.546,3.551,3.555,3.560,
     .    3.564,3.569,3.573,3.578,3.582,3.586,3.591,3.595,3.599,3.604,
     .    3.608,3.612,3.616,3.621,3.625,3.629,3.633,3.637,3.641,3.645,
     .    3.649,3.653,3.657,3.661,3.665,3.669,3.673,3.677,3.681,3.685/
      data p7a/
     .    1.398,1.398,1.398,1.398,1.406,1.425,1.443,1.461,1.480,1.498,
     .    1.516,1.534,1.551,1.568,1.585,1.601,1.616,1.631,1.646,1.660,
     .    1.674,1.687,1.700,1.712,1.724,1.736,1.747,1.758,1.768,1.778,
     .    1.788,1.797,1.806,1.815,1.824,1.832,1.840,1.848,1.855,1.863,
     .    1.870,1.877,1.883,1.890,1.896,1.902,1.908,1.914,1.920,1.925,
     .    1.931,1.936,1.941,1.946,1.951,1.956,1.960,1.965,1.969,1.974,
     .    1.978,1.982,1.986,1.990,1.994,1.998,2.001,2.005,2.009,2.012,
     .    2.016,2.019,2.022,2.026,2.029,2.032,2.035,2.038,2.041,2.044,
     .    2.047,2.050,2.053,2.056,2.059,2.061,2.064,2.067,2.069,2.072,
     .    2.075,2.077,2.080,2.082,2.085,2.088,2.090,2.093,2.095,2.098,
     .    2.100,2.103,2.105,2.107,2.110,2.112,2.115,2.117,2.120,2.122,
     .    2.125,2.127,2.130,2.132,2.135,2.137,2.140,2.142,2.145,2.148,
     .    2.150,2.153,2.155,2.158,2.161,2.163,2.166,2.169,2.172,2.175,
     .    2.178,2.180,2.183,2.186,2.189,2.192,2.195,2.198,2.202,2.205,
     .    2.208,2.211,2.215,2.218,2.221,2.225,2.228,2.232,2.235,2.239,
     .    2.243,2.246,2.250,2.254,2.258,2.261,2.265,2.269,2.273,2.277,
     .    2.282,2.286,2.290,2.294,2.299,2.303,2.307,2.312,2.316,2.321,
     .    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,
     .    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420/
      data p7b/
     .    2.325,2.330,2.335,2.339,2.344,2.349,2.354,2.359,2.364,2.369,
     .    2.374,2.379,2.384,2.389,2.394,2.399,2.405,2.410,2.415,2.420,
     .    2.426,2.431,2.437,2.442,2.448,2.453,2.459,2.464,2.470,2.476,
     .    2.481,2.487,2.493,2.498,2.504,2.510,2.516,2.521,2.527,2.533,
     .    2.539,2.545,2.551,2.556,2.562,2.568,2.574,2.580,2.586,2.592,
     .    2.598,2.604,2.610,2.616,2.622,2.628,2.634,2.640,2.646,2.652,
     .    2.658,2.664,2.670,2.676,2.682,2.687,2.693,2.699,2.705,2.711,
     .    2.717,2.723,2.729,2.735,2.741,2.747,2.753,2.759,2.764,2.770,
     .    2.776,2.782,2.788,2.794,2.799,2.805,2.811,2.817,2.823,2.828,
     .    2.834,2.840,2.846,2.851,2.857,2.863,2.868,2.874,2.879,2.885,
     .    2.891,2.896,2.902,2.907,2.913,2.918,2.924,2.929,2.935,2.940,
     .    2.945,2.951,2.956,2.962,2.967,2.972,2.978,2.983,2.988,2.993,
     .    2.999,3.004,3.009,3.014,3.019,3.025,3.030,3.035,3.040,3.045,
     .    3.050,3.055,3.060,3.065,3.070,3.075,3.080,3.085,3.090,3.095,
     .    3.099,3.104,3.109,3.114,3.119,3.123,3.128,3.133,3.138,3.142,
     .    3.147,3.152,3.156,3.161,3.165,3.170,3.175,3.179,3.184,3.188,
     .    3.193,3.197,3.202,3.206,3.210,3.215,3.219,3.224,3.228,3.232/
      data p8a/
     .    1.447,1.447,1.447,1.447,1.447,1.447,1.459,1.475,1.489,1.504,
     .    1.518,1.531,1.544,1.556,1.568,1.580,1.591,1.602,1.612,1.622,
     .    1.631,1.640,1.649,1.658,1.666,1.674,1.682,1.689,1.696,1.703,
     .    1.710,1.716,1.722,1.728,1.734,1.740,1.745,1.751,1.756,1.761,
     .    1.766,1.770,1.775,1.779,1.784,1.788,1.792,1.796,1.800,1.804,
     .    1.807,1.811,1.814,1.818,1.821,1.824,1.827,1.831,1.834,1.836,
     .    1.839,1.842,1.845,1.848,1.850,1.853,1.855,1.858,1.860,1.863,
     .    1.865,1.867,1.870,1.872,1.874,1.876,1.878,1.880,1.882,1.884,
     .    1.886,1.888,1.890,1.892,1.894,1.896,1.898,1.900,1.902,1.903,
     .    1.905,1.907,1.909,1.911,1.912,1.914,1.916,1.917,1.919,1.921,
     .    1.923,1.924,1.926,1.928,1.929,1.931,1.933,1.934,1.936,1.938,
     .    1.939,1.941,1.943,1.945,1.946,1.948,1.950,1.951,1.953,1.955,
     .    1.957,1.959,1.960,1.962,1.964,1.966,1.968,1.970,1.971,1.973,
     .    1.975,1.977,1.979,1.981,1.983,1.985,1.987,1.989,1.991,1.993,
     .    1.995,1.998,2.000,2.002,2.004,2.006,2.009,2.011,2.013,2.015,
     .    2.018,2.020,2.023,2.025,2.027,2.030,2.032,2.035,2.037,2.040,
     .    2.043,2.045,2.048,2.051,2.053,2.056,2.059,2.062,2.064,2.067,
     .    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,
     .    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130/
      data p8b/
     .    2.070,2.073,2.076,2.079,2.082,2.085,2.088,2.091,2.094,2.097,
     .    2.100,2.103,2.107,2.110,2.113,2.116,2.120,2.123,2.126,2.130,
     .    2.133,2.137,2.140,2.143,2.147,2.151,2.154,2.158,2.161,2.165,
     .    2.168,2.172,2.176,2.180,2.183,2.187,2.191,2.195,2.198,2.202,
     .    2.206,2.210,2.214,2.218,2.222,2.226,2.230,2.233,2.237,2.241,
     .    2.245,2.250,2.254,2.258,2.262,2.266,2.270,2.274,2.278,2.282,
     .    2.286,2.291,2.295,2.299,2.303,2.307,2.312,2.316,2.320,2.324,
     .    2.329,2.333,2.337,2.341,2.346,2.350,2.354,2.359,2.363,2.367,
     .    2.371,2.376,2.380,2.384,2.389,2.393,2.397,2.402,2.406,2.410,
     .    2.415,2.419,2.423,2.428,2.432,2.436,2.441,2.445,2.449,2.454,
     .    2.458,2.462,2.467,2.471,2.475,2.480,2.484,2.488,2.493,2.497,
     .    2.501,2.506,2.510,2.514,2.519,2.523,2.527,2.531,2.536,2.540,
     .    2.544,2.548,2.553,2.557,2.561,2.565,2.570,2.574,2.578,2.582,
     .    2.586,2.591,2.595,2.599,2.603,2.607,2.611,2.616,2.620,2.624,
     .    2.628,2.632,2.636,2.640,2.644,2.648,2.652,2.656,2.661,2.665,
     .    2.669,2.673,2.677,2.681,2.685,2.689,2.693,2.696,2.700,2.704,
     .    2.708,2.712,2.716,2.720,2.724,2.728,2.732,2.736,2.739,2.743/
      data p9a/
     .    1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.322,1.325,
     .    1.334,1.342,1.351,1.358,1.366,1.373,1.380,1.386,1.392,1.398,
     .    1.404,1.409,1.415,1.420,1.425,1.429,1.434,1.438,1.442,1.446,
     .    1.450,1.454,1.457,1.461,1.464,1.467,1.470,1.473,1.476,1.479,
     .    1.482,1.485,1.487,1.490,1.492,1.495,1.497,1.499,1.501,1.503,
     .    1.505,1.507,1.509,1.511,1.513,1.515,1.517,1.519,1.520,1.522,
     .    1.524,1.525,1.527,1.528,1.530,1.531,1.533,1.534,1.535,1.537,
     .    1.538,1.539,1.541,1.542,1.543,1.545,1.546,1.547,1.548,1.549,
     .    1.551,1.552,1.553,1.554,1.555,1.556,1.558,1.559,1.560,1.561,
     .    1.562,1.563,1.565,1.566,1.567,1.568,1.569,1.570,1.571,1.573,
     .    1.574,1.575,1.576,1.577,1.579,1.580,1.581,1.582,1.584,1.585,
     .    1.586,1.588,1.589,1.590,1.592,1.593,1.594,1.596,1.597,1.599,
     .    1.600,1.602,1.603,1.605,1.606,1.608,1.609,1.611,1.612,1.614,
     .    1.616,1.617,1.619,1.621,1.622,1.624,1.626,1.628,1.630,1.631,
     .    1.633,1.635,1.637,1.639,1.641,1.643,1.645,1.647,1.649,1.651,
     .    1.653,1.655,1.657,1.659,1.661,1.664,1.666,1.668,1.670,1.673,
     .    1.675,1.677,1.679,1.682,1.684,1.686,1.689,1.691,1.694,1.696,
     .    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,
     .    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749/
      data p9b/
     .    1.699,1.701,1.704,1.706,1.709,1.711,1.714,1.716,1.719,1.722,
     .    1.724,1.727,1.729,1.732,1.735,1.738,1.740,1.743,1.746,1.749,
     .    1.751,1.754,1.757,1.760,1.763,1.765,1.768,1.771,1.774,1.777,
     .    1.780,1.783,1.786,1.789,1.792,1.795,1.798,1.801,1.804,1.807,
     .    1.810,1.813,1.816,1.819,1.822,1.825,1.828,1.831,1.834,1.837,
     .    1.840,1.843,1.847,1.850,1.853,1.856,1.859,1.862,1.865,1.869,
     .    1.872,1.875,1.878,1.881,1.884,1.888,1.891,1.894,1.897,1.901,
     .    1.904,1.907,1.910,1.913,1.917,1.920,1.923,1.926,1.930,1.933,
     .    1.936,1.939,1.943,1.946,1.949,1.952,1.956,1.959,1.962,1.965,
     .    1.969,1.972,1.975,1.978,1.982,1.985,1.988,1.992,1.995,1.998,
     .    2.001,2.005,2.008,2.011,2.014,2.018,2.021,2.024,2.027,2.031,
     .    2.034,2.037,2.040,2.044,2.047,2.050,2.053,2.057,2.060,2.063,
     .    2.066,2.070,2.073,2.076,2.079,2.083,2.086,2.089,2.092,2.095,
     .    2.099,2.102,2.105,2.108,2.111,2.115,2.118,2.121,2.124,2.127,
     .    2.131,2.134,2.137,2.140,2.143,2.146,2.149,2.153,2.156,2.159,
     .    2.162,2.165,2.168,2.171,2.175,2.178,2.181,2.184,2.187,2.190,
     .    2.193,2.196,2.199,2.202,2.205,2.208,2.212,2.215,2.218,2.221/
c
      if(t.lt.12000.) then
        pf=g0(ion-3)
        dut=0.
        dun=0.
        return
      endif
c
      it=int(t/1000)
      if(it.ge.350) it=349
      t1=1000.*it
      t2=t1+1000.
      if(ion.eq.4) then
        if(t.le.200000.) then
          xu1=p4a(it-10)
          xu2=p4a(it-9)
        else
          xu1=p4b(it-180)
          xu2=p4b(it-179)
        endif
      else if(ion.eq.5) then
        if(t.le.200000.) then
          xu1=p5a(it-10)
          xu2=p5a(it-9)
        else
          xu1=p5b(it-180)
          xu2=p5b(it-179)
        endif
      else if(ion.eq.6) then
        if(t.le.200000.) then
          xu1=p6a(it-10)
          xu2=p6a(it-9)
        else
          xu1=p6b(it-180)
          xu2=p6b(it-179)
        endif
      else if(ion.eq.7) then
        if(t.le.200000.) then
          xu1=p7a(it-10)
          xu2=p7a(it-9)
        else
          xu1=p7b(it-180)
          xu2=p7b(it-179)
        endif
      else if(ion.eq.8) then
        if(t.le.200000.) then
          xu1=p8a(it-10)
          xu2=p8a(it-9)
        else
          xu1=p8b(it-180)
          xu2=p8b(it-179)
        endif
      else if(ion.eq.9) then
        if(t.le.200000.) then
          xu1=p9a(it-10)
          xu2=p9a(it-9)
        else
          xu1=p9b(it-180)
          xu2=p9b(it-179)
        endif
      endif
c
      dxt=xmil*(xu2-xu1)
      xu=xu1+(t-t1)*dxt
      pf=exp(xen*xu)
      dut=xen*pf*dxt
      dun=0.
      return
      end
c 
c ******************************************************************
c
C
      SUBROUTINE PFHEAV(IIZ,JNION,MODE,t,ane,u)
C     =========================================
C
c     subset of kurucz's pfsaha for Z>28.
c     removed code for Z<28; crp- 28 aug, 1995
C     EDITED 27 JULY 1994 BY GMW - REPLACED PT III PF COEFF. AND IP             
C     MODE 3 RETURNS PARTITION FUNCTION   
C                                      
C      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'PARAMS.FOR' 
      REAL*8 IP
      PARAMETER (DEBCON=1./2.8965E-18,
     *           TVCON=8.6171E-5,
     *           HIONEV=13.595,
     *           ONE=1.,
     *           HALF=0.5,
     *           THIRD=1./3.,
     *           X18=1./18.,
     *           X120=1./120.,
     *           T211=2000./11.)
c                                                 
C     DIMENSION F(6),
      DIMENSION IP(6),PART(6),POTLO(6)
C     DIMENSION FSAVE(6)
      DIMENSION SCALE(4)                                                        
      DIMENSION NNN(6*218)                                                      
      DIMENSION NNN16(54),NNN17(54),NNN18(54),NNN19(54),NNN20(54)               
      DIMENSION NNN21(54),NNN22(54),NNN23(54),NNN24(54),NNN25(54)               
      DIMENSION NNN26(54),NNN27(54),NNN28(54),NNN29(54),NNN30(54)               
      DIMENSION NNN31(54),NNN32(54),NNN33(54),NNN34(54),NNN35(54)               
      DIMENSION NNN36(54),NNN37(54),NNN38(54),NNN39(54),NNN40(12)               
      EQUIVALENCE (NNN( 811-810),NNN16(1))                     
      EQUIVALENCE (NNN( 865-810),NNN17(1)),(NNN( 919-810),NNN18(1))
      EQUIVALENCE (NNN( 973-810),NNN19(1)),(NNN(1027-810),NNN20(1))
      EQUIVALENCE (NNN(1081-810),NNN21(1)),(NNN(1135-810),NNN22(1))
      EQUIVALENCE (NNN(1189-810),NNN23(1)),(NNN(1243-810),NNN24(1))
      EQUIVALENCE (NNN(1297-810),NNN25(1)),(NNN(1351-810),NNN26(1)) 
      EQUIVALENCE (NNN(1405-810),NNN27(1)),(NNN(1459-810),NNN28(1)) 
      EQUIVALENCE (NNN(1513-810),NNN29(1)),(NNN(1567-810),NNN30(1))
      EQUIVALENCE (NNN(1621-810),NNN31(1)),(NNN(1675-810),NNN32(1)) 
      EQUIVALENCE (NNN(1729-810),NNN33(1)),(NNN(1783-810),NNN34(1)) 
      EQUIVALENCE (NNN(1837-810),NNN35(1)),(NNN(1891-810),NNN36(1)) 
      EQUIVALENCE (NNN(1945-810),NNN37(1)),(NNN(1999-810),NNN38(1)) 
      EQUIVALENCE (NNN(2053-810),NNN39(1)),(NNN(2107-810),NNN40(1)) 
C      ( 1)( 2)   ( 3)( 4)   ( 5)( 6)   ( 7)( 8)   ( 9)(10)   ( IP ) G
      DATA NNN16/                                                               
     1 227027622, 306233052, 356839222, 446052912, 652382292,   763314,
     2 108416342, 222428472, 353944332, 577378932, 110314303,  1814900,
     3 198724282, 293236452, 468362702,  86511123, 136016073,  3516000,
     4 279836622, 461857562, 720693022, 124915873, 192522633,  5600000,
     5 262136422, 501167232,  87911303, 138916483, 190721673,  7900000,
     6 201620781, 231026761, 314737361, 450555381, 692386911,   772301,
     7 109415761, 247938311,  58910042, 190937022,  68311693,  2028903,
     8 897195961, 107212972, 165021182, 260230862, 356940532,  3682900,
     9 100010001, 100410231, 108712611, 167124841, 388460411,   939102/
      DATA NNN17/                                                               
     1 200020021, 201620761, 223726341, 351352061,  80812472,  1796001,
     2 100610471, 122617301, 300566361, 149924112, 332342352,  3970000,
     3 403245601, 493151431, 529654331, 559358091, 611065171,   600000,
     4  99710051, 104511541, 135016501, 208226431, 321837921,  2050900,
     5 199820071, 204521391, 229124761, 266028451, 302932131,  3070000,
     6 502665261, 755183501, 901496201, 102410942, 117912812,   787900,
     7 422848161, 512153401, 557458941, 636270361, 794489061,  1593000,
     8 100010261, 114613921, 175221251, 249828711, 324436181,  3421000,
     9 403143241, 491856701, 649173781, 840396751, 113013392,   981000/
      DATA NNN18/                                                               
     1 593676641, 884697521, 105911572, 129515012, 180322212,  1858700,
     2 484470541,  91510972, 125614082, 157017612, 199722912,  2829900,
     3 630172361, 799686381, 919797221, 102810942, 117712832,   975000,
     4 438055511, 691582151,  94510732, 121413672, 152016732,  2150000,
     5 651982921,  94610382, 113212492, 139515462, 169718482,  3200000,
     6 437347431, 498951671, 538559501,  74710812, 169126672,  1183910,
     7 705183611,  93510092, 111614162, 222932532, 427652992,  2160000,
     8 510869921,  87410312, 123116552, 236530712, 377744832,  3590000,
     9 100010001, 100010051, 105012781, 198535971,  65911422,  1399507/
      DATA NNN19/                                                               
     1 461049811, 522254261, 609088131, 168935052,  68612253,  2455908,
     2 759990901, 101911142, 129017782, 302856642,  99414333,  3690000,
     3 200020011, 200720361, 211523021, 269434141, 459163351,   417502,
     4 100010001, 100110321, 129524961,  61014202, 291753192,  2750004,
     5 473650891, 533156051,  66810932, 232950852,  99915303,  4000000,
     6 100110041, 104111741, 146019721, 281941411, 607785251,   569202,
     7 202621931, 255331271, 384347931, 624085761, 122417632,  1102600,
     8 100010001, 100110321, 129524961,  61014202, 291753192,  4300000,
     9 791587851, 100012192, 155119942, 254031782, 389946932,   637900/
      DATA NNN20/                                                               
     1 118217102, 220827002, 319036792, 416646512, 513256072,  1223000,
     2  92510012, 104710862, 112311612, 120212472, 132814282,  2050000,
     3 141320802, 291439702, 531170262,  92712273, 162521053,   684000,
     4 354454352, 724689652, 107212643, 148517093, 193321573,  1312900,
     5 209727032, 324537052, 415446282, 510255752, 604965222,  2298000,
     6 256636022, 465759302, 749693962, 116514243, 171520333,   687900,
     7 335157222,  84511463, 147718363, 221826083, 299933893,  1431900,
     8 223725352, 280830972, 340937362, 406844002, 473150632,  2503900,
     9 703972941,  82610822, 154822682, 327244912, 571469372,   709900/
      DATA NNN21/                                                               
     1  75714552, 274347322, 718897632, 123414913, 174920063,  1614900,
     2 267645462, 669890262, 115514323, 173620673, 242528083,  2714900,
     3  90613732, 184823562, 291735332, 419949102, 565764332,   728000,
     4 131318312, 227126932, 311735452, 397644072, 483852692,  1525900,
     5 204721673, 234725733, 284031463, 348738613, 426546943,  3000000,
     6 176824122, 318941082, 515263202, 761790472, 106112303,   736400,
     7 221934642, 501968372,  88911173, 136316243, 189221613,  1675900,
     8 210622722, 241025422, 267928262, 297731272, 327834282,  2846000,
     9 148520202, 255230902, 364942462, 489656082, 638872352,   746000/
      DATA NNN22/                                                               
     1 153421292, 288137912, 484660322, 720187062, 101011483,  1807000,
     2 254537212, 492362292, 770592182, 107312243, 137615273,  3104900,
     3 115919651, 320746011, 607576761,  95011642, 141817172,   832900,
     4 755087211, 105913442, 173122222, 282034722, 412247732,  1941900,
     5 180223462, 289735212, 414247632, 538460052, 662672472,  3292000,
     6 200020001, 200220141, 206422141, 257633021, 455164681,   757403,
     7 100810581, 125817401, 260641031,  66210072, 135316982,  2148000,
     8 795887491,  97711762, 156620252, 248329422, 340038582,  3481900,
     9 100010001, 100410241, 109212891, 176827421, 444268771,   899003/
      DATA NNN23/                                                               
     1 200020021, 201720921, 233329881, 451475371, 127520782,  1690301,
     2 100310281, 114815371, 246138311, 519265531, 791492761,  3747000,
     3 252431921, 368440461, 433746521, 512259221, 723389021,   578400,
     4 100110071, 104611651, 146118581, 225426511, 304734431,  1886000,
     5 200120111, 205021611, 243628031, 317035371, 390442701,  2802900,
     6 232637101, 488058571, 669074381, 816189091,  97210632,   734200,
     7 286335941, 408144471, 479351961, 571862901, 686274341,  1462700,
     8 100010251, 114013811, 175321601, 256829751, 338337901,  3049000,
     9 404043481, 494656811, 646772781, 813490751, 101411372,   863900/
      DATA NNN24/                                                               
     1 303147981, 618472951, 827392621, 103711702, 131214532,  1650000,
     2 313037601, 429347901, 536260591, 689477591, 862494881,  2529900,
     3 526258801, 657372351, 784284071, 897095741, 102711082,   900900,
     4 440855541, 686481251,  93810792, 125414792, 176321132,  1860000,
     5 349054751, 699883081,  96611302, 134216202, 197724212,  2800000,
     6 405342041, 438645621, 475751071, 587974491, 102214572,  1045404,
     7 568567471, 773485861,  94510362, 112712182, 130914002,  1909000,
     8 514269581,  86910562, 130716652, 215327742, 351843662,  3200000,
     9 100010001, 100010091, 109515351, 291060661, 119621482,  1212716/
      DATA NNN25/                                                               
     1 414844131, 465649111, 538464651,  87112232, 158019362,  2120000,
     2 615475101, 867797531, 112213462, 157618062, 203622662,  3209900,
     3 200020001, 201020501, 215623871, 283536181, 462756261,   389300,
     4 100010001, 100310371, 119016501, 269146361,  77912412,  2510000,
     5 424445601, 481750061, 516953311, 549356551, 581759791,  3500000,
     6 101210791, 135119351, 282340571, 574580391, 111015062,   521002,
     7 262638611, 504160621, 698579371,  91010692, 129115952,  1000000,
     8 100010001, 100310351, 118416321, 264945521,  76512182,  3700000,
     9  71111992, 172323592, 312540402, 510763182, 765791012,   558000/
      DATA NNN26/                                                      
     1 204529582, 383647882, 582469262, 807992692, 104911723,  1106000,
     2  94712552, 148416582, 179819212, 203621522, 227424042,  1916900,
     3 295959132, 103515693, 215527593, 335939413, 449650223,   565000,
     4  79718153, 289639443, 495159253, 686877533, 863794813,  1085000,
     5 298640242, 475053692, 596965912, 725379692, 872094692,  2008000,
     6 460693672, 158523823, 327242303, 519661563, 709379783,   541900,
     7 455480232, 114014653, 178521013, 240927073, 299232633,  1055000,
     8  46410533, 183826893, 354443773, 518459633, 674375243,  2320000,
     9 139623042, 364860002,  96114603, 209828633, 373446973,   549000/
      DATA NNN27/                                                               
     1 460493692, 158523823, 327142303, 519661563, 709279783,  1073000,
     2 455480232, 114014653, 178521013, 240927073, 299232633,  2000000,
     3 131720482, 280535692, 441254492, 676583972, 103412583,   555000,
     4 139623042, 364860002,  96114603, 209828633, 373446973,  1089900,
     5 460493682, 158523823, 327142303, 519661563, 709279783,  2000000,
     6  92915672, 222431062, 444763802,  89612173, 159520253,   562900,
     7 315059662,  97114563, 204627093, 342541693, 490556383,  1106900,
     8 269037812, 520270372,  91111273, 133915483, 172719093,  2000000,
     9 800080571, 851699301, 127617362, 240433032, 444958442,   568000/
      DATA NNN28/                                                               
     1 125416052, 211828182, 375549622, 644381732, 101112213,  1125000,
     2 800080571, 851699301, 127617362, 240433032, 444958442,  2000000,
     3 240432982, 427555202, 708489962, 112613853, 167319843,   615900,
     4 534793262, 139219123, 247730843, 371043333, 495055893,  1210000,
     5 364145232, 514756362, 604864112, 673870372, 732276072,  2000000,
     6 480767202,  89011393, 144118243, 230028753, 354142883,   584900,
     7 480767192,  89011393, 144118243, 230028753, 354142883,  1151900,
     8 480767202,  89011393, 144118243, 230028753, 354142883,  2000000,
     9 343147532, 645887152, 115314793, 183322063, 257729373,   593000/
      DATA NNN29/                                                               
     1 343147532, 645887142, 115314793, 183322063, 257729373,  1167000,
     2 343147532, 645887142, 115314793, 183322063, 257729373,  2000000,
     3 222635002, 542276772, 100312353, 145716713, 187020703,   602000,
     4 222635002, 542276772, 100312353, 145716713, 187020703,  1180000,
     5 222635002, 542276772, 100312353, 145716713, 187020703,  2000000,
     6 133715382, 209130152, 429859382,  79410293, 129815983,   609900,
     7 265934782, 497877532, 120517733, 245032063, 400448073,  1193000,
     8 265934782, 497877532, 120517733, 245032063, 400448073,  2000000,
     9 800381111,  87510702, 147621462, 310343462, 585475982,   618000/
      DATA NNN30/                                                               
     1 156718872, 279244452, 678196342, 128316243, 197823443,  1205000,
     2  93517192, 364666132, 103414613, 192624193, 293334613,  2370000,
     3 100010011, 101310651, 118613951, 169120661, 250629971,   625000,
     4 200120901, 270345231,  81714042, 223533112, 461959862,  1217000,
     5 100312561, 250851931,  91914182, 198626022, 323638692,  2000000,
     6 514664441, 759086851,  99211442, 133315612, 182721252,   609900,
     7 125924831, 438667801,  98714112, 199727872, 380850742,  1389900,
     8 323948621, 661297271, 158626482, 426865032,  93712843,  1900000,
     9 659294081, 128016962, 222528952, 372047062, 585171462,   700000/
      DATA NNN31/                                                               
     1  99117882, 274638812, 520867322,  84410313, 123314453,  1489900,
     2 187427702, 343739872, 448049452, 539358282, 625266642,  2329900,
     3  65210892, 171325762, 373552252, 705192012, 116414343,   787900,
     4 192837842, 600784802, 111113823, 165419233, 218524383,  1620000,
     5  99117872, 274638812, 520867312,  84410313, 123314453,  2400000,
     6 398981651, 130019172, 273438022, 516168382,  88411163,   797900,
     7 131429482, 523279952, 111414623, 183422233, 262130233,  1770000,
     8 192837842, 600784792, 111113823, 165419233, 218524383,  2500000,
     9 600963001,  75910412, 150121572, 301940972, 539168952,   787000/
      DATA NNN32/                                                               
     1  73710852, 190731262, 464964142,  83810503, 127315053,  1660000,
     2 131429482, 523279952, 111414623, 183422233, 262130233,  2600000,
     3 110815502, 216829732, 398752322, 672484682, 104612673,   850000,
     4 168225972, 362046562, 566766422, 757484612,  93010103,  1700000,
     5  73710852, 190731262, 464964142,  83810503, 127315053,  2700000,
     6 129117892, 239430882, 388748292, 596173252,  89510843,   910000,
     7 110815502, 216829732, 398752322, 672484682, 104612673,  2000000,
     8 168225972, 362046562, 566766422, 757484612,  93010103,  2800000,
     9 158918512, 207523002, 254328242, 316335762, 407246582,   900000/
      DATA NNN33/                                                               
     1  98115462, 224930742, 401150612, 623475412,  89910583,  1855900,
C    2 110815502, 216829732, 398752322, 672484682, 104612673,  2900000,
     2 146323292, 354651802,  74810923, 161723953, 348749363,  3322700,
     3 203222611, 265731251, 364042301, 494958601, 702084731,   922000,
     4 120521331, 357753801,  75310062, 130516572, 206925452,  2050000,
     5 651780821, 108814772, 195925252, 316338622, 460853882,  3000000,
     6 100010001, 100110111, 105211851, 152122101, 341552811,  1043002,
     7 200320211, 210023021, 268834231, 480472341, 111416912,  1875000,
     8 104012871, 186129471, 458664151,  82410072, 119013732,  3420000,
     9 200420711, 222424271, 265429161, 325637371, 442853911,   610500/
      DATA NNN34/                                                               
     1 100010021, 101910801, 121414641, 189525811, 358949721,  2041900,
     2 200020311, 216624611, 296337451, 489064791,  85711212,  2979900,
     3 103411711, 147819101, 244331781, 434862751,  93113762,   741404,
     4 204122231, 248227841, 311535621, 429153941, 651976431,  1502800,
     5 100210131, 106812201, 154522671, 381665951,  95512512,  3192900,
     6 400140351, 416944121, 474851591, 564362181, 690477231,   728700,
     7 106814451, 204427341, 350744811, 586879131, 108314772,  1667900,
     8 205523051, 264830231, 345439921, 469156001, 675281671,  2555900,
     9 500950661, 518153561, 559058941, 628968071, 748483501,   843000/
      DATA NNN35/                                                               
     1 443756241, 696282451,  95411012, 128615262, 182922012,  1900000,
     2 336953201, 682481011,  93810882, 127915272, 184622442,  2700000,
     3 402841621, 431544771, 463148311, 520059491, 734896851,   930000,
     4 576168741, 788387631,  96910642, 116012552, 135014462,  2000000,
     5 490265341, 812797201, 116614322, 179622692, 285035302,  2900000,
     6 100010001, 100010031, 102311051, 133018071, 264539391,  1074500,
     7 402841621, 431544771, 463148311, 520059491, 734996851,  2000000,
     8 576168741, 788387631,  96910642, 116012552, 135014462,  3000000,
     9 200020011, 201220591, 218124481, 296538611, 488859141,   400000/
      DATA NNN36/                                                               
     1 100010001, 100010031, 102311051, 133018071, 264539401,  2200000,
     2 421645151, 477449611, 511852711, 542455761, 572958821,  3300000,
     3 100010041, 105212131, 153220271, 270435641, 460258111,   527600,
     4 201221791, 258131471, 381645781, 546365131, 777592781,  1014400,
     5 100010001, 100010031, 102311051, 133018071, 264539391,  3400000,
     6 510064491,  82710872, 142718412, 232328712, 348341572,   690000,
     7 228951571,  88513232, 183324132, 305537492, 448152402,  1210000,
     8 723989131, 103511752, 130814352, 155416652, 177018682,  2000000,
     9 620099241, 162725772, 391457072,  80110833, 141818023,   600000/
      DATA NNN37/                                                               
     1 620099241, 162725772, 391457072,  80110833, 141818023,  1200000,
     2 620099251, 162725772, 391457072,  80110833, 141818023,  2000000,
     3 347877992, 129318323, 240730533, 380546863, 570368573,   600000,
     4 347877992, 129318323, 240730533, 380546863, 570368573,  1200000,
     5 347777992, 129318323, 240730533, 380546863, 570368573,  2000000,
     6 209530092, 450866762,  96613623, 186524763, 318839893,   600000,
     7 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000/
      DATA NNN38/                                                               
     1 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     2 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     3 209530092, 450866762,  96613623, 186524763, 318839893,   600000,
     4 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     5 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     6 209530092, 450866762,  96613623, 186524763, 318839893,   600000,
     7 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000/
      DATA NNN39/                                                               
     1 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     2 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     3 209530092, 450866762,  96613623, 186524763, 318839893,   600000,
     4 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     5 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     6 209530092, 450866762,  96613623, 186524763, 318839893,   600000,
     7 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     8 209530092, 450866762,  96613623, 186524763, 318839893,  2000000,
     9 209530092, 450866762,  96613623, 186524763, 318839893,   600000/
      DATA NNN40/                                                               
     1 209530092, 450866762,  96613623, 186524763, 318839893,  1200000,
     2 209530092, 450866762,  96613623, 186524763, 318839893,  2000000/
      DATA SCALE/.001,.01,.1,1./    
C
      if(mode.lt.0) return
      tk=1.38054d-16*t
      tv=8.6171d-5*t
C     LOWERING OF THE IONIZATION POTENTIAL IN VOLTS FOR UNIT ZEFF               
      CHARGE=ANE*2.                                                         
      DEBYE=SQRT(TK*DEBCON/CHARGE)                                       
C     DEBYE=SQRT(TK/12.5664/4.801E-10**2/CHARGE)                             
      POTLOW=MIN(1.D0,1.44E-7/DEBYE)                                            
      IF(IIZ.LE.28)then
         write(6,*) 'Error, routine PFHEAV for Z.GE.28 only'
         stop23
       endif
c removed elements with z<28
      if(iiz.eq.28) n=1                                                    
      IF(IIZ.GT.28) N=3*IIZ+54-135                                                     
      IF(IIZ.eq.28) NIONS=4
      IF(IIZ.GT.28) NIONS=3                                                       
      NION2=MIN0(JNION+2,NIONS)                                                  
      N=N-1                                                                     
C                                                                               
      DO 18 ION=1,NION2                                                         
      Z=ION                                                                     
      POTLO(ION)=POTLOW*Z                                                       
      N=N+1 
      nnn6n=nnn(6+6*(N-1))  
c     nnn6n=nnn(6,n)                                                                  
      NNN100=NNN6N/100 
      XN1= NNN100                                                     
      IP(ION)=XN1*1.e-3                                               
      IG=NNN6N-NNN100*100  
      GGG=IG                                                   
      T2000=IP(ION)*T211                                                  
      IT=MAX0(1,MIN0(9, INT(T/T2000-HALF))) 
      XIT=IT                                   
      DT=T/T2000-XIT-HALF                                                
      PMIN=ONE                                                                   
      I=(IT+1)/2                                                                
      nnnin=nnn(i+6*(N-1))  
c     nnnin=nnn(i,n)                                                                  
      K1=NNNIN/100000                                                        
      K2=NNNIN-K1*100000                                                     
      K3=K2/10                                                                  
      xk1=k1                                                
      xk3=k3
      KSCALE=K2-K3*10                                                           
      IF(MOD(IT,2).EQ.0)GO TO 12                                                
      P1=XK1*SCALE(KSCALE)                                                
      P2=XK3*SCALE(KSCALE)                                                
      IF(DT.GE.0.)GO TO 13                                                      
      IF(KSCALE.GT.1)GO TO 13                                                   
      KP1=int(P1)
      IF(KP1.NE. INT(P2+.5))GO TO 13                                            
      PMIN=KP1                                                                  
      GO TO 13                                                                  
   12 continue
      xk3=k3
      P1=XK3*SCALE(KSCALE)                                                
      nnni1n=nnn(i+1+6*(N-1))  
c     nnni1n=nnn(i+1,n)                                                                  
      K1=NNNI1N/100000                                                      
      KSCALE=MOD(NNNI1N,10) 
      xk1=k1                                                
      P2=XK1*SCALE(KSCALE)                                                
   13 PART(ION)= MAX (PMIN,P1+(P2-P1)*DT)                                       
      IF(GGG.EQ.0..OR.POTLO(ION).LT..1.OR.T.LT.T2000*4.)GO TO 18               
      IF(T.GT.(T2000*11.)) TV=(T2000*11.)*TVCON                          
      D1=.1/TV                                                                  
      D2=POTLO(ION)/TV  
      DX=SQRT(HIONEV*Z*Z/TV/D2)**3                                                       
      PART(ION)=PART(ION)+GGG*EXP(-IP(ION)/TV)*       
     *          (DX*(THIRD+(ONE-(HALF+(X18+D2*X120)*D2)*D2)*D2)-                              
     *           DX*(THIRD+(ONE-(HALF+(X18+D1*X120)*D1)*D1)*D1))                              
   18 CONTINUE                                                                  
      u=part(jnion)
      RETURN                                                                    
      END                                                                       
c 
c ******************************************************************
c
      subroutine frac1
c     ================
c
      include 'PARAMS.FOR'
      include 'MODELP.FOR'
      parameter (mtemp=100,melec=60,mion1=30)
      dimension xxt(mdepth),xxe(mdepth)
      dimension kt0(mdepth),kn0(mdepth)
      common/fracop/frac(mtemp,melec,mion1),fracm(mtemp,melec),
     *              itemp(mtemp),ntt
c
      do 10 id=1,nd
         xxt(id)=dlog10(temp(id))
         kt0(id)=2*int(20.*xxt(id))
         xxe(id)=dlog10(elec(id))
         kn0(id)=int(2.*xxe(id))
   10 continue
c
      DO 20 IAT=1,30
         iatnum=iat
         call fractn(iatnum)
         if(iatnum.le.0) goto 20
         do 30 id=1,nd
           if(kt0(id).lt.itemp(1)) then
             kt1=1
             write(6,611) id,temp(id)
  611        format(' (FRACOP) Extrapol. in T (low)',i4,f7.0)
             goto 41
           endif
           if(kt0(id).ge.itemp(ntt)) then
             kt1=ntt-1
             write(6,612) id,temp(id)
  612        format(' (FRACOP) Extrapol. in T (high)',i4,f12.0)
             goto 41
           endif
           do 40 it=1,ntt
             if(kt0(id).eq.itemp(it)) then
               kt1=it
               goto 41
             endif
   40      continue
   41      continue
           if(kn0(id).lt.1) then
             kn1=1
             goto 49
           endif
           if(kn0(id).ge.60) then
             kn1=59
             write(6,614) id,xxe(id)
  614        format(' (FRACOP) Extrapol. in Ne (high)',i4,f9.4)
             goto 49
           endif
           kn1=kn0(id)
   49      continue
           xt1=0.025*itemp(kt1)
           dxt=0.05
           at1=(xxt(id)-xt1)/dxt
           xn1=0.5*kn1
           dxn=0.5
           an1=(xxe(id)-xn1)/dxn
           do 50 ion=1,mion1
              x11=frac(kt1,kn1,ion)
              x21=frac(kt1+1,kn1,ion)
              x12=frac(kt1,kn1+1,ion)
              x22=frac(kt1+1,kn1+1,ion)
              x1221=x11*x21*x12*x22
              if(x1221.eq.0.) then
                  xx1=x11+at1*(x21-x11)
                  xx2=x12+at1*(x22-x12)
                  rrx=xx1+an1*(xx2-xx1)
              else
                  x11=dlog10(x11)
                  x21=dlog10(x21)
                  x12=dlog10(x12)
                  x22=dlog10(x22)
                  xx1=x11+at1*(x21-x11)
                  xx2=x12+at1*(x22-x12)
                  rrx=xx1+an1*(xx2-xx1)
                  rrx=exp(2.3025851*rrx)
              endif
              rrr(id,ion,iat)=rrx*abndd(iat,id)*
     *                        dens(id)/wmm(id)/ytot(id)
   50      continue
   30    continue
   20 CONTINUE
c
      return
      end
c 
c ******************************************************************
c
      subroutine fractn(iatnum)
c     =========================
c
      implicit double precision (a-h,o-z)
      parameter (mtemp=100,
     *           melec= 60,
     *           mion1=30,
     *           mdat = 17)
      parameter (inp=71)
      dimension frac0(-1:mion1),ioo(-1:mion1),idat(mion1)
      dimension gg(mion1,mdat),g0(mion1),z0(-1:mion1)
      dimension uu(mion1,mdat),u0(mion1)
      dimension u6(6),u7(7),u8(8),u10(10),u11(11)
      dimension u12(12),u13(13),u14(14),u16(16),u18(18),u20(20)
      dimension u24(24),u25(25),u26(26),u28(28)
      equivalence (u6(1),uu(1,3)),(u7(1),uu(1,4)),(u8(1),uu(1,5))
      equivalence (u10(1),uu(1,6)),(u11(1),uu(1,7)),(u12(1),uu(1,8))
      equivalence (u13(1),uu(1,9)),(u14(1),uu(1,10)),(u16(1),uu(1,11))
      equivalence (u18(1),uu(1,12)),(u20(1),uu(1,13)),(u24(1),uu(1,14))
      equivalence (u25(1),uu(1,15)),(u26(1),uu(1,16)),(u28(1),uu(1,17))
      common/fracop/frac(mtemp,melec,mion1),fracm(mtemp,melec),
     *              itemp(mtemp),ntt
      data idat   / 1, 2, 0, 0, 0, 3, 4, 5, 0, 6,
     *              7, 8, 9,10, 0,11, 0,12, 0,13,
     *              0, 0, 0,14,15,16, 0,17, 0, 0/ 
      data gg/2.,29*0.,2.,1.,28*0.,
     *        2.,1.,2.,1.,6.,9.,24*0.,2.,1.,2.,1.,6.,9.,4.,23*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,22*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,20*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,19*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,18*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,17*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,16*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,14*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,6.,1.,
     *        12*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *        6.,1.,2.,1.,10*0.,2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,
     *        6.,9.,4.,9.,6.,1.,10.,21.,28.,25.,6.,7.,6*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,7.,6.,5*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,30.,25.,4*0.,
     *        2.,1.,2.,1.,6.,9.,4.,9.,6.,1.,2.,1.,6.,9.,4.,9.,
     *           6.,1.,10.,21.,28.,25.,6.,25.,28.,21.,10.,21.,0.,0./
      data uu(1,1)/109.6787/
      data uu(1,2)/198.3108/
      data uu(2,2)/438.9089/
      data u6/90.82,196.665,386.241,520.178,3162.395,3952.061/
      data u7/117.225,238.751,382.704,624.866,789.537,4452.758,5380.089/
      data u8/109.837,283.24,443.086,624.384,918.657,1114.008,5963.135,
     *        7028.393/
      data u10/173.93,330.391,511.8,783.3,1018.,1273.8,1671.792,
     *         1928.462,9645.005,10986.876/
      data u11/41.449,381.395,577.8,797.8,1116.2,1388.5,1681.5,2130.8,
     *         2418.7,11817.061,13297.676/
      data u12/61.671,121.268,646.41,881.1,1139.4,1504.3,1814.3,2144.7,
     *         2645.2,2964.4,14210.261,15829.951/
      data u13/48.278,151.86,229.446,967.8,1239.8,1536.3,1947.3,2295.4,
     *         2663.4,3214.8,3565.6,16825.022,18584.138/
      data u14/65.748,131.838,270.139,364.093,1345.1,1653.9,1988.4,
     *         2445.3,2831.9,3237.8,3839.8,4222.4,19661.693,21560.63/
      data u16/83.558,188.2,280.9,381.541,586.2,710.184,2265.9,2647.4,
     *         3057.7,3606.1,4071.4,4554.3,5255.9,5703.6,26002.663,
     *         28182.535/
      data u18/127.11,222.848,328.6,482.4,605.1,734.04,1002.73,1157.08,
     *         3407.3,3860.9,4347.,4986.6,5533.8,6095.5,6894.2,7404.4,
     *         33237.173,35699.936/
      data u20/49.306,95.752,410.642,542.6,681.6,877.4,1026.,1187.6,
     *         1520.64,1704.047,4774.,5301.,5861.,6595.,7215.,7860.,
     *         8770.,9338.,41366.,44177.41/
      data u24/54.576,132.966,249.7,396.5,560.2,731.02,1291.9,1490.,
     *         1688.,1971.,2184.,2404.,2862.,3098.52,8151.,8850.,
     *         9560.,10480.,11260.,12070.,13180.,13882.,60344.,63675.9/
      data u25/59.959,126.145,271.55,413.,584.,771.1,961.44,1569.,
     *         1789.,2003.,2307.,2536.,2771.,3250.,3509.82,9152.,
     *         9872.,10620.,11590.,12410.,13260.,14420.,15162.,
     *         65660.,69137.4/
      data u26/63.737,130.563,247.22,442.,605.,799.,1008.,1218.38,
     *         1884.,2114.,2341.,2668.,2912.,3163.,3686.,3946.82,
     *         10180.,10985.,11850.,12708.,13620.,14510.,15797.,
     *         16500.,71203.,74829.6/
      data u28/61.6,146.542,283.8,443.,613.5,870.,1070.,1310.,1560.,
     *         1812.,2589.,2840.,3100.,3470.,3740.,4020.,4606.,
     *         4896.2,12430.,13290.,14160.,15280.,16220.,17190.,
     *         18510.,19351.,82984.,86909.4/
c
      if(idat(iatnum).eq.0) then
         write(6,600) iatnum
  600    format(' OP data for element no. ',i3,' do not exist')
         iatnum=-1
         return
      end if
c
      g0(iatnum+1)=1.
      do i=1,iatnum
        ig0=iatnum-i+1
        g0(ig0)=gg(i,idat(iatnum))
        u0(i)=uu(i,idat(iatnum))*1000.
      enddo
c
      if(iatnum.eq.1) open(inp,file='ioniz.dat',status='old')
      do 10 it=1,mtemp
         do 10 ie=1,melec
            fracm(it,ie)=0.
            do 10 ion=1,mion1
               frac(it,ie,ion)=0.
   10 continue
c
      read(inp,*)
      read(inp,*) it0,it1,itstp
      ntt=(it1-it0)/itstp+1
c
      do 100 it=1,ntt
         read(inp,*) itt,ie0,ie1,iestp
         itemp(it)=itt
         net=(ie1-ie0)/iestp+1
         t=exp(2.3025851*0.025*itt)
         safac0=sqrt(t)*t/2.07d-16
         tkcm=0.69496*t
         do 30 ie=1,net
            read(inp,601) iee,ion0,ion1,
     *                    (ioo(i),frac0(i),i=ion0,min(ion1,ion0+3))
            ane=exp(2.3025851*0.25*iee)
            safac=safac0/ane
            nio=ion1-ion0
            if(nio.ge.3) then
               nlin=nio/4
               do ilin=1,nlin
                  read(inp,602) (ioo(i),frac0(i),
     *                 i=ion0+4*ilin,min(ion1,ion0+4*ilin+3))
               end do
            end if
            ieind=iee/2
            do 20 ion=ion0,ion1
              if(ion.lt.iatnum) then
               if(ion.eq.ion0) then
                  z0(ion)=g0(iatnum-ion)
               else
                  z0(ion)=frac0(ion)/frac0(ion-1)*safac*z0(ion-1)
                  z0(ion)=z0(ion)*exp(-u0(iatnum-ion)/tkcm)
               endif
                  frac(it,ieind,iatnum-ion)=frac0(ion)/z0(ion)
              else
                  u0hm=6090.5
                  z0hm=frac0(ion)/frac0(ion-1)*safac
                  z0hm=z0hm*exp(-u0hm/tkcm)
                  fracm(it,ieind)=frac0(ion)/z0hm
              end if
   20       continue
c           write(6,603) it,ieind,t,ane
c           write(6,604) (frac(it,ieind,i),i=1,mion1)
   30    continue
  100 continue
  601 format(3i4,2x,4(i4,1x,e9.3))
  602 format(14x,4(i4,1x,e9.3))
      return      
      end
C
C
C     *******************************************************************
C
C

      SUBROUTINE DWNFR0(ID)
C     =====================
C
C     Auxiliary quantities for dissolved fractions
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (UN=1.,SIXTH=UN/6.,CCOR=0.09)
      parameter (p1=0.1402,p2=0.1285,p3=un,p4=3.15,p5=4.)
      parameter (f23=-2./3.)
C
      ANE=ELEC(ID)
      ELEC23(ID)=EXP(F23*LOG(ANE))
      ANES=EXP(SIXTH*LOG(ANE))
      ACOR=CCOR*ANES/SQRT(TEMP(ID))
      X=EXP(P4*LOG(UN+P3*ACOR))
      DWC2(ID)=P2*X
      A3=ACOR*ACOR*ACOR
      DO 10 IZZ=1,MZZ
         Z3(IZZ)=IZZ*IZZ*IZZ
         DWC1(IZZ,ID)=P1*(X+P5*(IZZ-1.)*A3)
   10 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE DWNFR1(FR,FR0,ID,IZZ,DW1)
C     ====================================
C
C     dissolved fraction for frequency FR
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      PARAMETER (UN=1.,TKN=3.01,CKN=5.33333333,CB=8.59d14)
      PARAMETER (SQFRH=5.734152D7)
      parameter (a0=0.529177e-8,wa0=-3.1415926538/6.*a0*a0*a0)
C
      IF(FR.LT.FR0) THEN
         XN=SQFRH*IZZ/SQRT(FR0-FR)
         if(xn.le.tkn) then
            xkn=un
          else
            xn1=un/(xn+un)
            xkn=ckn*xn*xn1*xn1
         end if
         BETA=CB*Z3(IZZ)*XKN/(XN*XN*XN*XN)*ELEC23(ID)
         beta=beta*bergfc
         BETA3=BETA*BETA*BETA
         BETA32=SQRT(BETA3)
         F=(DWC1(IZZ,ID)*BETA3)/(UN+DWC2(ID)*BETA32)
c
c     contribution from neutral particles
c
         xn2=xn*xn+un
         xnh=0.
         xnhe1=0.
         if(ielh.gt.0) xnh=popul(nfirst(ielh),id)
         if(ielhe1.gt.0) xnhe1=popul(nfirst(ielhe1),id)
         w0=exp(wa0*xn2*xn2*xn2*(xnh+xnhe1))
         W0=1.
c
         DW1=UN-F/(UN+F)*w0
       ELSE
         DW1=UN
      END IF
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE CHCKAB
C
C     check input abumdances of explicit atoms (unit 5) and those
C     which follow from the models atmosphere (unit 7) obtained by
C     summing all populations and upper sums
C     The program stops if it finds discrepancy more than 10 %
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      dimension sumpop(matom),sumiat(matom)
c
      IST=0
      DO ID1=1,3
      IF(ID1.EQ.1) ID=1
      IF(ID1.EQ.2) ID=46
      IF(ID1.EQ.3) ID=ND
      CALL WNSTOR(ID)
      ANE=ELEC(ID)
      CALL SABOLF(ID)
      DO IAT=1,NATOM
         SUM=0.
         sump=0.
         DO I=N0A(IAT),NKA(IAT) 
            IL=ILK(I)
            A=1.
            IF(IL.GT.0) A=1.+ANE*USUM(IL)
            SUM=SUM+A*POPUL(I,ID)
            SUMP=SUMP+POPUL(I,ID)
         END DO
         SUMIAT(IAT)=SUM
         SUMPOP(IAT)=SUMP
      END DO
      WRITE(6,600) ID
      DO IAT=1,NATOM
         X=SUMIAT(IAT)/SUMIAT(IATREF)
         WRITE(6,601) IAT,X,abund(iat,id),SUMPOP(IAT)/SUMPOP(IATREF)
         IF(X/abund(iat,id).GT.1.1.OR.X/abund(iat,id).LT.0.9) ist=ist+1
      END DO
      END DO
      IF(IST.GT.0) THEN
         WRITE(6,602) 
         STOP
      END IF
  600 FORMAT(' check of abundances (id =',i3/
     *       ' computed from model atmosphere  - input abundances'/)
  601 format(i5,1p3e20.3)
  602 format(' ERROR !!! - inconsistent abundances'/)
      RETURN
      END
C
C
C ********************************************************************
C
C     
      subroutine molini
c     =================
c
c     Initialization of the molecular equilibrium
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      common/moltst/pfmol(500,mdepth),anmol(500,mdepth),
     *              pfato(100,mdepth),anato(100,mdepth),
     *              pfion(100,mdepth),anion(100,mdepth)
      dimension hpo(mdepth)
c
      aeinit=1.0
c
      do 10 id=1,nd
         t=temp(id)
         tln=log(t)*1.5
         thl=11605./t
         t32=exp(tln)

         do i=1,MMOLEC
            rrmol(i,id)=0.
         end do
 
         hpo(id)=DENS(ID)/WMM(ID)/YTOT(ID)

         if(t.gt.tmolim) go to 10
         HPOP=DENS(ID)/WMM(ID)/YTOT(ID)

         an=dens(id)/wmm(id)+elec(id)
         aeinit=0.1*an
         if(t.lt.4000.) aeinit=0.01*an
         call moleq(id,t,an,aeinit,ane,0)
         

c     next initial guess will be the last ane determined for
c     previous depth point

         aeinit=ane
c
         if (id.eq.idstd) then
            write(6,600)
            nmol=nmolec
            if(id.eq.1) nmol=32
            do i=1,nmol
               write(6,601) i, cmol(i), rrmol(i,id), rrmol(i,id)/hpop
            end do
         end if
 600     format(/ 'Molecular number densities at the standard depth'/)
 601     format(i4,1x,A8,1x,1pe12.2,1x,e12.2)

   10 continue


c     update atomic populations once molecular densities are calculated

      if(imode.lt.-4) then
      do i=1,nlevel
         iat=numat(iatm(i))
         ion=iz(iel(i))
         ii=nfirst(iel(i))
         ener=(enion(ii)-enion(i))/bolk
         if((enion(i).eq.0).and.(ilk(i).gt.0)) then
            ener=0.
            ion=ion+1
         end if
         if(ifwop(i).ge.0) then
            do id=1,nd
               popul(i,id)=rrr(id,ion,iat)*g(i)
     *              *exp(-ener/temp(id))
               if(iat.eq.1.and.ion.eq.0) popul(i,id)=anhm(id)
            end do
         endif
      end do
      end if
c
c     do im=1,32
c        write(78,671) im,cmol(im)
c        do id=1,nd
c           write(78,672) id,temp(id),anmol(im,id),rrmol(im,id),
c    *                    pfmol(im,id),anmol(im,id)/hpo(id)
c        end do
c     end do
c 671 format(//i5,a10/)
c 672 format(i4,0pf10.1,1p4e12.4)

      return
      end
C
C
C ********************************************************************
C
C
      SUBROUTINE INMOLI(ILIST)
C     ========================
C
C     read in the input molecular line list,
C     selection of lines that may contribute,
C     set up auxiliary fields containing line parameters,
C
C     Input of line data - unit 20:
C
C     For each line, one (or two) records, containing:
C
C    ALAM    - wavelength (in nm)
C    ANUM    - code of the modelcule (as in Kurucz)
C              (eg. 101.00 = H2; 607.00 = CN)
C    GF      - log gf
C    EXCL    - excitation potential of the lower level (in cm*-1)
C    GR      - gamma(rad)
C    GS      - gamma(Stark)
C    GW      - gamma(VdW)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      COMMON/NXTINM/ALMM00,ALSM00
      common/alendm/alend(mmlist)
      PARAMETER (PI4=7.95774715E-2)
      PARAMETER (C1     = 2.3025851,
     *           C2     = 4.2014672,
     *           C3     = 1.4387886,
     *           CNM    = 2.997925D17,
     *           EXT0   = 3.17,
     *           UN     = 1.0,
     *           TEN    = 10.,
     *           HUND   = 1.D2,
     *           TENM4  = 1.D-4,
     *           TENM8  = 1.D-8,
     *           OP4    = 0.4,
     *           AGR0=2.4734E-22, 
     *           XEH=13.595, XET=8067.6, XNF=25.,
     *           R02=2.5, R12=45., VW0=4.5E-9)
C
c      DATA INLSET /0/
C
      if(imode.ne.-3.and.temp(idstd).gt.tmolim) return
      IUNIT=IUNITM(ILIST)

      if(ibin(ilist).eq.0) then
         open(unit=iunit,file=amlist(ilist),status='old')
       else
         open(unit=iunit,file=amlist(ilist),form='unformatted',
     *        status='old')
      end if
C
c     define a conversion table between Kurucz notation and Tsuji table 
c     through array MOLIND
C
      do i=1,11000
         molind(i)=0
      end do
      molind(101)=2
      molind(106)=5
      molind(107)=12
      molind(108)=4
      molind(111)=122
      molind(112)=32
      molind(114)=17
      molind(116)=16
      molind(120)=34
      molind(124)=198
      molind(126)=214
      molind(606)=8
      molind(607)=7
      molind(608)=6
      molind(614)=21
      molind(616)=20
      molind(707)=9
      molind(708)=11
      molind(714)=24
      molind(716)=23
      molind(808)=10
      molind(812)=126
      molind(813)=134
      molind(814)=25
      molind(816)=26
      molind(820)=179
      molind(822)=29
      molind(823)=30
      molind(10108)=3
c
c     iunit=19+ilist
C
      ALAST=CNM/FRLAST
      ALASTM(ILIST)=ALAST
      IL=0
      IF(NXTSEM(ILIST).EQ.1) THEN
          ALAM0=ALM00
          ALASTM(ILIST)=ALST00
          FRLASM(ILIST)=CNM/ALASTM(ILIST)
          NXTSEM(ILIST)=0
          REWIND IUNIT
      END IF
      ALMM00=ALAM0
c     ALASTM(ILIST)=CNM/FRLAST
c     FRLASM(ILIST)=CNM/ALASTM(ILIST)
      DOPSTD=1.E7/ALAM0*DSTD
      DOPLAM=ALAM0*ALAM0/CNM*DOPSTD
      AVAB=ABSTD(IDSTD)*RELOP
      ASTD=1.0
c     IF(GRAV.GT.6.) ASTD=0.1
      CUTOFF=CUTOF0
      ALAST=CNM/FRLAST
C
C     first part of reading line list - read only lambda, and
C     skip all lines with wavelength below ALAM0-CUTOFF
C
      REWIND IUNIT
      ALAM=0.
      IJC=2
c
    7 if(ibin(ilist).eq.0) then
         READ(IUNIT,510) ALAM
       else
         read(iunit) alam
      end if
  510 FORMAT(F10.4)
      IF(ALAM.LT.ALAM0-CUTOFF) GO TO 7
      BACKSPACE(IUNIT)
      GO TO 10
c
    8 continue 
   10 continue
      if(ibin(ilist).eq.0) then
         if(ivdwli(ilist).eq.0) then
           READ(IUNIT,*,END=100,err=8) ALAM,ANUM,GF,EXCL,GR,GS,GW
          else
           read(iunit,*,end=100) alam,anum,gf,excl,gr,gh2,xnh2,ghe,xnhe
         end if
       else
         if(ivdwli(ilist).eq.0) then
           READ(IUNIT,END=100) ALAM,ANUM,GF,EXCL,GR,GS,GW
          else
           read(iunit,end=100) alam,anum,gf,excl,gr,gh2,xnh2,ghe,xnhe
c             GH2=0.
c             GHE=0.
         end if
      end if
C
c     change wavelength to vacuum for lambda > 2000
c
      if(alam.gt.200..and.vaclim.gt.2000.) then
         wl0=alam*10.
         ALM=1.E8/(WL0*WL0) 
         XN1=64.328+29498.1/(146.-ALM)+255.4/(41.-ALM)
         WL0=WL0*(XN1*1.D-6+UN)
         alam=wl0*0.1
      END IF        
C
C     first selection : for a given interval
C
      IF(ALAM.GT.ALASTM(ILIST)+CUTOFF) GO TO 100
C
C     second selection : for line strengths
C
      FR0=CNM/ALAM
      icod=int(anum+tenm4)
      imol=molind(icod)
      if(imol.le.0.or.imol.gt.nmolec) go to 10
      EXCL=ABS(EXCL)
      GFP=C1*GF-C2
      EPP=C3*EXCL
      gx=gfp-epp/tstd
      ab0=0.
c
      if(ndstep.eq.0.and.ifwin.eq.0) then
c
c      old procedure for line rejection
c
         if(gx.gt.-30)
     *   AB0=EXP(GFP-EPP/TSTD)*RRMOL(IMOL,IDSTD)/DOPSTD/AVAB
         IF(AB0.LT.UN) GO TO 10
       else
c
c      new procedure for line rejection
c
         do ijcn=ijc,nfreqc
            if(fr0.ge.freqc(ijcn)) go to 12
         end do
   12    continue
         ijc=ijcn
         if(ijc.gt.nfreqc) ijc=nfreqc
c
         tkm=1.65e8/ammol(imol)
         DP0=3.33564E-11*FR0
         do id=1,nd,ndstep
            td=temp(id)
            gx=gfp-epp/td
            ab0=0.
            if(gx.gt.-30) then
               dops=dp0*sqrt(tkm*td+vturb(id))
               AB0=EXP(gx)*RRMOL(IMOL,ID)/(DOPS*abstdw(ijc,id)*relop)
            end if
            if(ab0.ge.un) go to 15
         end do
         GO TO 10
      end if
c
C     truncate line list if there are more lines than maximum allowable
C     (given by MLIN0 - see include file LINDAT.FOR)
C
   15 CONTINUE
      IL=IL+1
      IF(IL.GT.MLINM0) THEN
         WRITE(6,601) ALAM
         IL=MLINM0
         ALASTM(ILIST)=CNM/FREQM(IL,ILIST)-CUTOFF
         FRLASM(ILIST)=CNM/ALASTM(ILIST)
         NXTSEM(ILIST)=1
         GO TO 100
      END IF
C
C     =============================================
C     line is selected, set up necessary parameters
C     =============================================
C
C     evaluation of EXTIN0 - the distance (in delta frequency) where
C     the line is supposed to contribute to the total opacity
C
      EX0=AB0*ASTD*10.
      EXT=EXT0
      IF(EX0.GT.TEN) EXT=SQRT(EX0)
      EXTIN0=EXT*DOPSTD
C
C     store parameters for selected lines
C
      FREQM(IL,ILIST)=FR0
      EXCLM(IL,ILIST)=real(EPP)
      GFM(IL,ILIST)=real(GFP)
      EXTINM(IL,ILIST)=real(EXTIN0)
      INDATM(IL,ILIST)=imol
C
C     ****** line broadening parameters *****
C            assuming for Stark 1.e-8*effnsq**5/2, with effnsq=25
C
      GRM(IL,ILIST)=real(GR*PI4)
      GSM(IL,ILIST)=real(GS*PI4*3.125e-5)
      GWM(IL,ILIST)=real(GW*PI4)
      if(ivdwli(ilist).ne.0) then
          gvdwh2(il,ilist)=real(gh2)
          gexph2(il,ilist)=real(xnh2)
          gvdwhe(il,ilist)=real(ghe)
          gexphe(il,ilist)=real(xnhe)
          gsm(il,ilist)=0.
          gwm(il,ilist)=0.
      end if
C
      GO TO 10
  100 NLINM0(ILIST)=IL
      nlinmt(ilist)=nlinmt(ilist)+nlinm0(ilist)
      alend(ilist)=cnm/fr0
C
      xln=float(il)*1.e-6
      WRITE(6,611) IUNIT,trim(amlist(ilist)),XLN
  611 FORMAT(/' --------------------------------------------'/
     *' MOLECULAR LINES - FROM UNIT ',i3,
     *', FILE  ',a,':',f8.3,' M'/
     *' --------------------------------------------'/)
  601 FORMAT('0 **** MORE LINES THAN MLINM0, LINE LIST TRUNCATED '/
     *'       AT LAMBDA',F15.4,'  NM'/)
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE MOLSET(ILIST)
C     ========================
C
C     Selection of molecular lines that may contribute,
C     set up auxiliary fields containing line parameters.
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      common/alendm/alend(mmlist)
      SAVE IMLAST
C
      DATA CNM /2.997925D17/
C
      if(inactm(ilist).ne.0) return
      IL0=0
      IPRSEM(ILIST)=0
      NLINM=0
      IREADM(ILIST)=1
      IF(IBLANK.LE.1.OR.IMODE.EQ.1.OR.IMODE.EQ.-1) IREADM(ILIST)=0
      IF(IBLANK.LE.1) APREV=0.
      ALA0=CNM/FREQ(1)
      ALA1=CNM/FREQ(2)
c
c     skip if current wavelength larger than the largest wavelngth in the
c     line list
c
      if(ala0.gt.alend(ilist)) then
         inactm(ilist)=1
         return
      end if
c
      FRMINM=CNM/ALA0
      FRM=FRMINM
      SPACE=SPACE0
      IF(ALAMC.GT.0.) SPACE=SPACE0*ALA0/ALAMC
      IF(SPACE0.LT.0.) SPACE=-SPACE0

      CUTOFF=CUTOF0*0.2
      DOPSTD=1.E7/ALA0*DSTD
      DISTAN=0.15*DOPSTD
      SPAC=3.E16/ALA0/ALA0*SPACE
      DISTA0=0.14*SPAC
      IF(IBLANK.GE.2.AND.IMODE.EQ.-1) IL0=IMLAST
      FRLI0=FRMINM
      ASTD=1.0
      AVAB=ABSTD(IDSTD)*RELOP
C
   20 CONTINUE 
C
C     set up indices of lines
C     IL0 - is the current index of line in the numbering of all lines
C
      IF(IREADM(ILIST).EQ.1) THEN
         IPRSEM(ILIST)=IPRSEM(ILIST)+1
         IL0=INMLIP(IPRSEM(ILIST),ILIST)
         IF(FREQM(IL0,ILIST).LT.FRMINM) THEN
            IREADM(ILIST)=0
            IL0=INMLIP(IPRSEM(ILIST)-1,ILIST)+1
         END IF
c      write(*,*) 'molset',ilist,ireadm(ilist),iprsem(ilist),il0
       ELSE
         IL0=IL0+1
c      write(*,*) 'molset',ilist,ireadm(ilist),il0
      END IF
      IF(IL0.GT.NLINM0(ILIST)) GO TO 210
      FRLIM=FRLI0
      FR0=FREQM(IL0,ILIST)
      ALAM=CNM/FR0
C
      IF(ALAM.LT.ALA0-CUTOFF) GO TO 20
      IF(ALAM.GT.ALA1+CUTOFF) GO TO 210
C
C     SECOND SELECTION : FOR LINE STRENGHTS
C
      EXT=EXTINM(IL0,ILIST)
      FRLI0=FR0-EXT-SPAC
      IF(FRLI0.GT.FRLIM) FRLI0=FRLIM
      IF(ALAM.LT.ALA0.AND.FR0-FRMINM.GT.EXT+SPAC) GO TO 20
      IF(FREQ(NFREQS)-FR0.GT.EXT+SPAC) GO TO 20
C
      NLINM=NLINM+1
      if(nlinm.gt.mlinm) then
         write(*,*) 'nlinm,mlinm',nlinm,mlinm
         call quit('too many molecular lines in a set')
      end if
      INMLIN(NLINM,ILIST)=IL0
      GO TO 20
c
c     frequency indices of the line centers
c
  210 CONTINUE
      XX=FREQ(2)-FREQ(1)
      DFRCON=NFREQ-3
      DFRCON=-DFRCON/XX
      IFRCON=INT(DFRCON)
      DO 255 IL=1,NLINM
         fr0=freqm(inmlin(il,ilist),ILIST)
         XJC=3.+DFRCON*(FREQ(1)-FR0)
         IJC=INT(XJC)
         IJCMTR(IL,ILIST)=IJC
         if(ijc.le.3.or.ijc.ge.nfreq) go to 255
         if(fr0.lt.freq(ijc)) then
            ijc0=ijc
            dfr0=freq(ijc0)-fr0
  252       ijc0=ijc0+1
            dfr=abs(freq(ijc0)-fr0)
            if(dfr.lt.dfr0) then
               ijc=ijc0
               ijc0=ijc0+1
               dfr0=dfr
               go to 252
            end if
          else if(fr0.gt.freq(ijc)) then
            ijc0=ijc
            dfr0=fr0-freq(ijc0)
  254       ijc0=ijc0-1
            dfr=abs(freq(ijc0)-fr0)
            if(dfr.lt.dfr0) then
               ijc=ijc0
               ijc0=ijc0-1
               dfr0=dfr
               go to 254
            end if
         end if
         IJCMTR(IL,ILIST)=IJC
  255 continue
C
      DO 260 IL=1,NLINM
         INMLIP(IL,ILIST)=INMLIN(IL,ILIST)
  260 CONTINUE
      NLINML(ILIST)=NLINM
      IMLAST=INMLIN(NLINML(ILIST),ILIST)
C
      CALL INIBLM
C
c      write(6,611) inmlin(1,ilist),inmlin(nlinm,ilist),
c    * 2.997925e18/freqm(inmlin(1,ilist),ILIST),
c    * 2.997925e18/freqm(inmlin(nlinm,ilist),ILIST)
c 611  format('mols',2i7,2f10.3)
      RETURN
      END
C
C
C ********************************************************************
C
C

      SUBROUTINE INIBLM
C     =================
C
C     driving procedure for treating a partial molecular line list for the
C     current wavelength region
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
C
      PARAMETER (DP0=3.33564E-11, DP1=1.651E8, UN=1.) 
C
      XX=FREQ(1)
      IF(NFREQ.GE.2) XX=0.5*(FREQ(1)+FREQ(2))
      BNU=BN*(XX*1.E-15)**3
      HKF=HK*XX
      DO 20 ID=1,ND
         T=TEMP(ID)
         EXH=EXP(HKF/T)
         EXHK(ID)=UN/EXH
         PLAN(ID)=BNU/(EXH-UN)
         STIM(ID)=UN-EXHK(ID)
         DO 10 IMOL=1,NMOLEC
            IF(AMMOL(IMOL).GT.0.)
     *      DOPMOL(IMOL,ID)=UN/(XX*DP0*SQRT(DP1*T/AMMOL(IMOL)+
     *                      VTURB(ID)))
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
C
C ********************************************************************
C
C
      SUBROUTINE IDMTAB
C     =================
C
C     output of selected molecular line parameters (identification table)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'

      COMMON/REFDEP/IREFD(MFREQ)
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)

      CHARACTER*4 APB,AP0,AP1,AP2,AP3,AP4,APR
C

      PARAMETER (C1=2.3025851, C2=4.2014672, C3=1.4387886)
      DATA APB,AP0,AP1,AP2,AP3,AP4 /'    ','   .','   *','  **',' ***',
     *                              '****'/
C
      ALM0=2.997925D18/FREQ(1)
      ALM1=2.997925D18/FREQ(2)
      if(ifwin.gt.0) ALM1=2.997925D18/FREQ(NFREQ)
      IF(IPRIN.LE.-2) RETURN
      if(iprin.ge.3) then
      IF(IMODE.GE.0) WRITE(6,601) IBLANK,ALM0,ALM1
      IF(IMODE.GE.0.OR.(IMODE.EQ.-1.AND.IBLANK.EQ.1)) WRITE(6,602)
      end if
C
      ID=IDSTD
      DO 100 ILIST=1,NMLIST
      IF(NLINML(ILIST).EQ.0) GO TO 100
      DO 20 IL0=1,NLINML(ILIST)
         IL=INMLIN(IL0,ILIST) 
         ALAM=2.997925D18/FREQM(IL,ILIST)
c        ID=IDSTD
         IJCN=IJCMTR(IL0,ILIST)
c        IF(IJCN.GE.1.AND.IJCN.LE.NFREQS) ID=IREFD(IJCN)
         IMOL=INDATM(IL,ILIST)
         DOP1=DOPMOL(IMOL,ID)
         ANE=ELEC(ID)
         AGAM=(GRM(IL,ILIST)+GSM(IL,ILIST)*ANE+
     *         GVDW(IL,ILIST,ID))*DOP1
         ABCNT=EXP(GFM(IL,ILIST)-EXCLM(IL,ILIST)/TEMP(ID))*
     *         RRMOL(IMOL,ID)*DOP1*STIM(ID)
c         STR0=ABCNT/ABSTD(ID)
         absta=min(ch(1,id),ch(2,id))
         str0=abcnt/absta
         if(ifwin.gt.0) STR0=ABCNT/ABSTDW(IJCONT(IL),ID)
         GF=(GFM(IL,ILIST)+C2)/C1
         EXCL=EXCLM(IL,ILIST)/C3
         IF(STR0.LE.1.2) THEN
            WW1=0.886*STR0*(1.-STR0*(0.707-STR0*0.577))
          ELSE 
            WW1=SQRT(LOG(STR0))
         END IF
         IF(STR0.GT.55.) THEN
            WW2=0.5*SQRT(3.14*AGAM*STR0)
            IF(WW2.GT.WW1) WW1=WW2
         END IF
         EQW=ALAM/FREQM(IL,ILIST)*1.E3/DOP1*WW1
         STR=EQW*10.
         APR=APB
         IF(STR.GE.1.E0.AND.STR.LT.1.E1) APR=AP0
         IF(STR.GE.1.E1.AND.STR.LT.1.E2) APR=AP1
         IF(STR.GE.1.E2.AND.STR.LT.1.E3) APR=AP2
         IF(STR.GE.1.E3.AND.STR.LT.1.E4) APR=AP3
         IF(STR.GE.1.E4) APR=AP4
         if(alam.ge.alm0.and.alam.lt.alm1) then
c        WRITE(15,603) IL0,IL,ALAM,CMOL(IMOL),GF,EXCL,
c    *                 STR0,EQW,APR,i0,i0,id
         WRITE(15,603) ALAM,CMOL(IMOL),GF,EXCL,
     *                 STR0,EQW,APR,id,AGAM
         end if
   20 CONTINUE
C
  601 FORMAT(/' ',I4,'. SET (MOLECULAR LINES):',
     * ' INTERVAL  ',F9.3,' -',F9.3,' ANGSTROMS'/
     *        ' ------------')
  602 FORMAT(/1H ,13X,
     * 'LAMBDA  MOLECULE  LOG GF       ELO    LINE/CONT',2X,
     * 'EQ.WIDTH',8x,'AGAM'/)
  603 FORMAT(F11.3,2X,A4,4X,F7.2,F12.3,1PE11.2,0PF8.1,1X,A4,
     *       i4,1PE10.2)
C
  100 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
 
      SUBROUTINE MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST)
C     ===========================================
C
C     Total molecular line opacity (ABLIN) and emissivity (EMLIN)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      PARAMETER (UN     = 1., 
     *           EXT0   = 3.17,
     *           TEN    = 10.)
      DIMENSION ABLIN(MFREQ),EMLIN(MFREQ)
C
      DO 10 IJ=1,NFREQ
         ABLIN(IJ)=0.
         EMLIN(IJ)=0.
   10 CONTINUE
C
      if(temp(id).gt.tmolim) return
      IF(NLINML(ILIST).EQ.0) RETURN
      if(inactm(ilist).ne.0) return
C
C     overall loop over contributing lines
C
      TEM1=UN/TEMP(ID)
      ANE=ELEC(ID)
      DO 100 I=1,NLINML(ILIST)
         IL=INMLIN(I,ILIST)
         IMOL=INDATM(IL,ILIST)
         DOP1=DOPMOL(IMOL,ID)
         AGAM=(GRM(IL,ILIST)+GSM(IL,ILIST)*ANE+
     *         GVDW(IL,ILIST,ID))*DOP1
         FR0=FREQM(IL,ILIST)
         AB0=EXP(GFM(IL,ILIST)-EXCLM(IL,ILIST)*TEM1)*RRMOL(IMOL,ID)*
     *           DOP1*STIM(ID)
C
C        set up limiting frequencies where the line I is supposed to
C        contribute to the opacity 
C
         EX0=AB0/AVAB*AGAM 
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
         EXT=EXT/DOP1
         XIJEXT=DFRCON*EXT+1.5
         IJ1=int(MAX(float(IJCMTR(I,ILIST))-XIJEXT,3.))
         IJ2=int(MIN(float(IJCMTR(I,ILIST))+XIJEXT,float(NFREQS)))
C         IJ1=MAX(IJCMTR(I,ILIST)-IJEXT,3)
C         IJ2=MIN(IJCMTR(I,ILIST)+IJEXT,NFREQS)
         IF(IJ1.GE.NFREQ.OR.IJ2.LE.2) GO TO 100
C
         DO 40 IJ=IJ1,IJ2
            XF=ABS(FREQ(IJ)-FR0)*DOP1
            ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF)
   40    CONTINUE
  100 CONTINUE
C
      DO 110 IJ=3,NFREQ
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID)
  110 CONTINUE
C
      RETURN
      END
C
C
C ********************************************************************
C
C
      FUNCTION SBFHMI(FR)
C     ===================
C    
C     Bound-free cross-section for H- (negative hydrogen ion)
C     Taken from Kurucz ATLAS9
C
C     FROM MATHISEN (1984), AFTER WISHART(1979) AND BROAD AND REINHARDT (1976)
C
      INCLUDE 'PARAMS.FOR'
      DIMENSION WBF(85),BF(85)
      DATA WBF/  18.00,  19.60,  21.40,  23.60,  26.40,  29.80,  34.30,
     1   40.40,  49.10,  62.60, 111.30, 112.10, 112.67, 112.95, 113.05,
     2  113.10, 113.20, 113.23, 113.50, 114.40, 121.00, 139.00, 164.00,
     3  175.00, 200.00, 225.00, 250.00, 275.00, 300.00, 325.00, 350.00,
     4  375.00, 400.00, 425.00, 450.00, 475.00, 500.00, 525.00, 550.00,
     5  575.00, 600.00, 625.00, 650.00, 675.00, 700.00, 725.00, 750.00,
     6  775.00, 800.00, 825.00, 850.00, 875.00, 900.00, 925.00, 950.00,
     7  975.00,1000.00,1025.00,1050.00,1075.00,1100.00,1125.00,1150.00,
     8 1175.00,1200.00,1225.00,1250.00,1275.00,1300.00,1325.00,1350.00,
     9 1375.00,1400.00,1425.00,1450.00,1475.00,1500.00,1525.00,1550.00,
     A 1575.00,1600.00,1610.00,1620.00,1630.00,1643.91/
      DATA BF/   0.067,  0.088,  0.117,  0.155,  0.206,  0.283,  0.414,
     1   0.703,   1.24,   2.33,  11.60,  13.90,  24.30,  66.70,  95.00,
     2   56.60,  20.00,  14.60,   8.50,   7.10,   5.43,   5.91,   7.29,
     3   7.918,  9.453,  11.08,  12.75,  14.46,  16.19,  17.92,  19.65,
     4   21.35,  23.02,  24.65,  26.24,  27.77,  29.23,  30.62,  31.94,
     5   33.17,  34.32,  35.37,  36.32,  37.17,  37.91,  38.54,  39.07,
     6   39.48,  39.77,  39.95,  40.01,  39.95,  39.77,  39.48,  39.06,
     7   38.53,  37.89,  37.13,  36.25,  35.28,  34.19,  33.01,  31.72,
     8   30.34,  28.87,  27.33,  25.71,  24.02,  22.26,  20.46,  18.62,
     9   16.74,  14.85,  12.95,  11.07,  9.211,  7.407,  5.677,  4.052,
     A   2.575,  1.302, 0.8697, 0.4974, 0.1989,    0. /
C    Bell and Berrington J.Phys.B,vol. 20, 801-806,1987.
c
      HMINBF=0.
      IF(FR.GT.1.82365E14) THEN
         WAVE=2.99792458E17/FR
         HMINBF=YLINTP(WAVE,WBF,BF,85,85)*1.E-18
      END IF
      SBFHMI=HMINBF
      RETURN
      END
C
C
C ********************************************************************
C
C

      FUNCTION SFFHMI(POPI,FR,T)
C     ==========================
C    
C     Free-free cross-section for H- (negative hydrogen ion)
C     Taken from Kurucz ATLAS9
C
C     From Bell and Berrington J.Phys.B,vol. 20, 801-806,1987.
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (CONFF=5040.*1.380658E-16, CONTH=5040.)
      DIMENSION FFLOG(22,11),FFCS(11,22),FFLOG2(22)
      DIMENSION FFBEG(11,11),FFEND(11,11),FFTT(11),WFFLOG(22)
      DIMENSION THETAFF(11),WAVEK(22)
      EQUIVALENCE (FFCS(1,1),FFBEG(1,1)),(FFCS(1,12),FFEND(1,1))
C
      DATA WAVEK/.50,.40,.35,.30,.25,.20,.18,.16,.14,.12,.10,.09,.08,
     1 .07,.06,.05,.04,.03,.02,.01,.008,.006/
      DATA THETAFF/
     1  0.5,  0.6, 0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.8,  3.6/       
      DATA FFBEG/
     1.0178,.0222,.0308,.0402,.0498,.0596,.0695,.0795,.0896, .131, .172,
     2.0228,.0280,.0388,.0499,.0614,.0732,.0851,.0972, .110, .160, .211,
     3.0277,.0342,.0476,.0615,.0760,.0908, .105, .121, .136, .199, .262,
     4.0364,.0447,.0616,.0789,.0966, .114, .132, .150, .169, .243, .318,
     5.0520,.0633,.0859, .108, .131, .154, .178, .201, .225, .321, .418,
     6.0791,.0959, .129, .161, .194, .227, .260, .293, .327, .463, .602,
     7.0965, .117, .157, .195, .234, .272, .311, .351, .390, .549, .711,
     8 .121, .146, .195, .241, .288, .334, .381, .428, .475, .667, .861,
     9 .154, .188, .249, .309, .367, .424, .482, .539, .597, .830, 1.07,
     A .208, .250, .332, .409, .484, .557, .630, .702, .774, 1.06, 1.36,
     B .293, .354, .468, .576, .677, .777, .874, .969, 1.06, 1.45, 1.83/
      DATA FFEND/
     1 .358, .432, .572, .702, .825, .943, 1.06, 1.17, 1.28, 1.73, 2.17,
     2 .448, .539, .711, .871, 1.02, 1.16, 1.29, 1.43, 1.57, 2.09, 2.60,
     3 .579, .699, .924, 1.13, 1.33, 1.51, 1.69, 1.86, 2.02, 2.67, 3.31,
     4 .781, .940, 1.24, 1.52, 1.78, 2.02, 2.26, 2.48, 2.69, 3.52, 4.31,
     5 1.11, 1.34, 1.77, 2.17, 2.53, 2.87, 3.20, 3.51, 3.80, 4.92, 5.97,
     6 1.73, 2.08, 2.74, 3.37, 3.90, 4.50, 5.01, 5.50, 5.95, 7.59, 9.06,
     7 3.04, 3.65, 4.80, 5.86, 6.86, 7.79, 8.67, 9.50, 10.3, 13.2, 15.6,
     8 6.79, 8.16, 10.7, 13.1, 15.3, 17.4, 19.4, 21.2, 23.0, 29.5, 35.0,
     9 27.0, 32.4, 42.6, 51.9, 60.7, 68.9, 76.8, 84.2, 91.4, 117., 140.,
     A 42.3, 50.6, 66.4, 80.8, 94.5, 107., 120., 131., 142., 183., 219.,
     B 75.1, 90.0, 118., 144., 168., 191., 212., 234., 253., 325., 388./
      DATA ISTART/0/
C
      IF(ISTART.EQ.0) THEN
      ISTART=1
      DO 2 IWAVE=1,22
         WFFLOG(IWAVE)=LOG(91.134D0/WAVEK(IWAVE))
         DO 2 ITHETA=1,11
            FFLOG(IWAVE,ITHETA)=LOG(FFCS(ITHETA,IWAVE)*1.E-26)
    2 CONTINUE
      ENDIF
C
      WAVE=2.99792458E17/FR
      WAVELOG=LOG(WAVE) 
C
      DO 21 ITHETA=1,11
         DO IWAVE=1,22
            FFLOG2(IWAVE)=FFLOG(IWAVE,ITHETA)
         END DO
         FFTLOG=YLINTP(WAVELOG,WFFLOG,FFLOG2,22,22)
         FFTT(ITHETA)=EXP(FFTLOG)/THETAFF(ITHETA)*CONFF
   21 CONTINUE
c
      THETA=CONTH/T
      FFTH=YLINTP(THETA,THETAFF,FFTT,11,11)
      SFFHMI=FFTH*POPI/(1.-exp(-hk*fr/t))
      RETURN
      END
C
C
C
C     ******************************************************************
C
C
C =========================================================================
C *************************************************************************
C *************************************************************************
C
      subroutine mpartf(jatom,ion,indmol,t,u)
c     =======================================
c
c     yields partition functions with polynomial data from
c     ref. Irwin, A.W., 1981, ApJ Suppl. 45, 621.
c     ln u(temp)=sum(a(i)*(ln(temp))**(i-1)) 1<=a<=6
c
c     Input:
c       jatom = element number in periodic table
c       ion   = 1 for neutral, 2 for once ionized and 3 for twice ionized
c       indmol= index of a molecular specie (Tsuji index)
c       temp  = temperature
c     Output:
c       u     = partf.(linear scale) for iat,ion, or indmol, and temperature t
c
c
      implicit real*8 (a-h,o-z)
      real*8 a(6,3,92),aa(6),am(6,500)
      dimension indtsu(324),irw(500)
      save iread,a,am
c     data indtsu / 2,  5, 12, 4, 8, 7, 6,
c    *              9, 11, 10, 29, 50, 59, 46, 132, 52, 19,
c    *             13, 42, 38, 39, 37, 44, 36, 14, 118, 33,
c    *              3, 16, 57, 32, 49, 60, 54, 41, 107,  0,
c    *            148, 152, 153, 155, 0, 17, 24, 25, 28, 51,
c    *            112, 119,   0,   0,21, 15, 43, 56,  0, 64,
c    *             47,  65,   0,  61, 0, 62,118, 40, 66/
      data indtsu / 2,  5, 12, 4, 8, 7, 6,
     *              9, 11, 10, 29, 50, 59, 46, 132,  52, 19,
     *             13, 42, 38, 39, 37, 44, 36,  14, 117, 33,
     *              3, 16, 57, 32, 49, 60, 54,  41, 106,303,
     *            147, 151, 152, 154, 302, 17,  24,  25, 28, 51,
     *            111, 118, 102,   0,  21, 15,  43,  56,478, 64,
     *             47,  65, 413,  61, 190, 62 ,108,  40, 66,214,
     *             257*0./

      data iread /0/
c
c     read data if first call:
c
      if(iread.ne.1) then
        if(irwtab.eq.0) then
           open(67,file= './data/irwin_bc.dat',status='old')
         else
           open(67,file= './data/irwin_orig.dat',status='old')
        end if
        read(67,*)
        read(67,*)
        do j=1,92
          do i=1,3
            if(j.eq.1.and.i.eq.3) goto 10
            sp=float(j)+float(i-1)/100.
            read(67,*) spec,aa
            do k=1,6
                a(k,i,j)=aa(k)
            end do
   10    continue
         end do
       end do 
c
       read(67,*)
       read(67,*)
       read(67,*)
       do i=1,500
          irw(i)=0
       end do
       do i=1,324
          read(67,*,end=15) spec,aa
          indm=indtsu(i)
          if(indm.gt.0) then
             irw(indm)=i
             do j=1,6
                am(j,indm)=aa(j)
             end do
          end if
        end do
   15   continue
        close(67)
        iread=1
      endif
c
c     evaluation of the partition function
c     stop if T is out of limits of Irwin's tables
c
        if(t.lt.1000.) then
          stop 'partf; temp<1000 K'
        else if(t.gt.16000.) then
          stop 'partf; temp>16000 K'
        endif
        tl=log(t)
        u=0.
c
c     atomic species
c
      if(jatom.gt.0.and.ion.gt.0) then
        ulog=    a(1,ion,jatom)+
     *       tl*(a(2,ion,jatom)+
     *       tl*(a(3,ion,jatom)+
     *       tl*(a(4,ion,jatom)+
     *       tl*(a(5,ion,jatom)+
     *       tl*(a(6,ion,jatom))))))
        u=exp(ulog)
      end if
c
c     molecular species
c
      if(indmol.gt.0) then
        indm=indmol
        if(irw(indm).gt.0) then
           ulog=    am(1,indm)+
     *          tl*(am(2,indm)+
     *          tl*(am(3,indm)+
     *          tl*(am(4,indm)+
     *          tl*(am(5,indm)+
     *          tl*(am(6,indm))))))
           u=exp(ulog)
        end if
      end if
      return
      end
      
C
C =========================================================================
C *************************************************************************
C *************************************************************************
C
C
      subroutine moleq(id,tt,an,aein,ane,ipri)
c     ========================================
c
c     calculation of the equilibrium state of atoms and molecules
c
c     Input:  id    - depth point
c             tt    - temperature [K]
c             an    - number density
c             aein  - initial estimate of the electron density
c
c     Output: ane    - electron density
c
C     Output through common/atomol:
c             rrr(id,j,i) - N/U for the atom with atomic number i and
c                     ion j (j=1 for neutral, and j=2 for 1st ions)
c             rrmol(imol,id) - N/U for the molecule with index imol
c                     (the index is given by the ordering of
c                     in the input file  tsuji.molec
c
c
c     Input data for molecules  iven in the file    
c     tsuji.molec             
c             
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      character*128 MOLEC
      COMMON/COMFH1/C(600,5),PPMOL(600),APMLOG(600),
     *              XIP(100),CCOMP(100),UIIDUI(100),P(100),
     *              FP(100),XKP(100),EPS,SWITER,
     *              NELEM(5,600),NATO(5,600),MMAX(600),
     *              NELEMX(100),NMETAL,NIMAX
      common/moltst/pfmol(500,mdepth),anmol(500,mdepth),
     *              pfato(100,mdepth),anato(100,mdepth),
     *              pfion(100,mdepth),anion(100,mdepth)
      DIMENSION NATOMM(5),NELEMM(5),
     *          emass(100),uelem(100),ull(100),anden(800),
     *          aelem(100)
c     dimension anion(100,mdepth)
      dimension denso(mdepth),eleco(mdepth),wmmo(mdepth)
c
      data nmetal/92/
c
      data iread/1/
c
      MOLEC ='data/tsuji.molec_bc2'
      if(moltab.eq.1) MOLEC='data/tsuji.molec_orig'
      if(moltab.eq.2) MOLEC='data/tsuji.molec'
c
        ECONST=4.342945E-1
        AVO=0.602217E+24
        SPA=0.196E-01
        GRA=0.275423E+05
        AHE=0.100E+00
        tk=1./(tt*1.38054e-16)
        pgas=an/tk
        sahcon=1.87840e20*tt*sqrt(tt)
        nimax=3000
        eps=1.e-5
        switer=0.0
C
C---- data for atoms  ---------------- 
C
      if(iread.eq.1) then 
c
      do i=1,nmetal
         ia=i
         nelemx(i)=ia
         ccomp(ia)=abndd(ia,id)
         xip(ia)=enev(ia,1)
         emass(ia)=amas(ia)
      end do
c
c---- read molecular data from a table  ----------------------
c
        J=0
        OPEN(UNIT=26,FILE=MOLEC,STATUS='OLD')
   10   J=J+1
        IF(MOLTAB.LE.1)
     *  READ (26,510,end=20) CMOL(J),(C(J,K),K=1,5),MMAX(J),
     *               (NELEMM(M),NATOMM(M),M=1,4)
        IF(MOLTAB.EQ.2)
     *  READ (26,511,end=20) CMOL(J),(C(J,K),K=1,5),MMAX(J),
     *               (NELEMM(M),NATOMM(M),M=1,4)
  510   format(a8,5e13.5,9i3)
  511   FORMAT (A8,E11.5,4E12.5,I1,(I2,I3),3(I2,I2))
c
c      for now, exclude all molecules with 4 or more C atoms
c
         do m=1,4
            if(nelemm(m).eq.6.and.natomm(m).ge.5) then
               j=j-1
               go to 10
            end if
         end do
c        
        MMAXJ=MMAX(J)
        IF(MMAXJ.EQ.0) GO TO 20
        DO M=1,MMAXJ
           NELEM(M,J)=NELEMM(M)
           NATO(M,J)=NATOMM(M)
        END DO
c        write(6,680) j,cmol(j)
c 680    format(i5,a10)
        GO TO 10
   20   NMOLEC=J-1
c       nmolec=297
c       nmolec=427
        close(26)
c        
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           P(NELEMI)=1.D-70
        END DO
        iread=0
      endif
c
c---- end of reading atomic and molecular data  ----------------------
c
        p(99)= aein/tk
        pesave=p(99)
        p(99)=pesave
c
        THETA=5040./tt
        TEM=tt
        PGLOG=log10(Pgas)
        PG=Pgas
c
        CALL RUSSEL(TEM,PG)
c
        PE=P(99)
        ane=pe*tk
        PELOG=log10(PE)
        emass(99)=5.486e-4
        uelem(99)=2.
        aelem(99)=pe*tk/(2.*sahcon*emass(nelemi)**1.5)
        ull(99)=log10(aelem(99))
c
c----atoms-----------------------------------------------------------------
c
        tmass=0.
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FPLOG=log10(FP(NELEMI))
           anden(i)=(p(nelemi)+1.D-70)*tk
           tmass=tmass+anden(i)*emass(nelemi)
           call mpartf(nelemi,1,0,tt,u0)
           uelem(nelemi)=u0
           aelem(nelemi)=anden(i)/(u0*sahcon*emass(nelemi)**1.5)
           ull(nelemi)=log10(aelem(nelemi))
           rrr(id,1,nelemi)=anden(i)/u0
           anato(nelemi,id)=anden(i)
           pfato(nelemi,id)=u0
        END DO
        an1=anden(1)
c
c---- positive ions ---------------------------------------------------------
c
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           PLOG= log10(P(NELEMI)+1.0D-70)
           XKPLOG=log10(XKP(NELEMI)+1.0D-70)
           PIONL=PLOG+XKPLOG-PELOG
           anden(i+nmetal)=exp(pionl/econst)*tk
           tmass=tmass+anden(i+nmetal)*emass(nelemi)
           call mpartf(nelemi,2,0,tt,u1)
           anion(nelemi,id)=anden(i+nmetal)
           pfion(nelemi,id)=u1
           rrr(id,2,nelemi)=anden(i+nmetal)/u1
        END DO
c
c     H- (inaccurate in Tsuji's table)
c
c           anmol1=1.0353e-16/tt/sqrt(tt)*exp(8762.9/tt)*
c   *              anato(1,id)*ane
c       pfmol1=1.
c       rrmol(j,id)=anden(1)
c
c---- molecules-------------------------------------------------------------
c
        DO J=1,NMOLEC
           jm=j+2*nmetal
           PMOLL=log10(PPMOL(J)+1.0D-70)
           anden(jm)=exp(pmoll/econst)*tk
c          if(j.eq.1) anden(jm)=anmol1
           rrmol(j,id)=0.
           umoll=1.
           if(pmoll.gt.-30.) then
              umoll=log10(anden(jm))+c(j,2)*theta
              amasm=0.
              do jjj=1,mmax(j)
                 i=nelem(jjj,j)
                 amasm=amasm+NATO(jjj,j)*emass(i)
                 umoll=umoll-NATO(jjj,j)*ull(i)
              end do
              ammol(j)=amasm
              tmass=tmass+anden(jm)*amasm
              umoll=exp(umoll/econst)/(sahcon*amasm**1.5)
c
c     replace with Irwin data whenever available
c
             call mpartf(0,0,j,tt,um)
             if(um.gt.0.) umoll=um
c
c     Exomol - POKAZATEL partition function for water
c
              if(j.eq.3) call h2opf(tt,umoll)
c
c     replace the TiO partition function by Kurucz-Schwenke data
c
c             IF(J.eq.29) CALL TIOPF(TT,UMOLL)
c     H-
c 
              if(j.eq.1) umoll=1.
c
c     set up array RRR = number density/partition function
c
              rrmol(j,id)=anden(jm)/umoll
           end if
c
           anmol(j,id)=anden(jm)
           pfmol(j,id)=umoll
        END DO

        jm=2*nmetal
        anhm(id)=anden(1+jm)
        anh2(id)=anden(2+jm)
        anch(id)=anden(5+jm)
        anoh(id)=anden(4+jm)
C
C
C     save new density, molecular weight, and abundances of
c     atomic species
c
      ipri1=ipri
      denso(id)=dens(id)
      eleco(id)=elec(id)
      wmmo(id)=wmm(id)
      dens(id)=tmass*hmass
      elec(id)=pe*tk
      wmm(id)=dens(id)/(an-elec(id))
      ane=elec(id)
c
      do i=1,nmetal
         NELEMI=NELEMX(I)
         ia=iatex(nelemi)
         if(ia.gt.0) then
            attot(ia,id)=(anato(nelemi,id)+anion(nelemi,id))
         end if
      end do
c
      if(id.eq.nd) then
        write(86,610)
        do iid=1,nd
           write(86,611) iid,dm(iid),temp(iid),elec(iid),eleco(iid),
     *                  dens(iid),denso(iid),wmm(iid),wmmo(iid)
        end do
        end if
  610 format(/'  id     m         T       ne(old)   ne(new)',
     *       '  dens(old) dens(new) wmm(old)  wmm(new)'/)
  611 format(i4,1p8e10.2)
C
        RETURN
        END

      
C
C =========================================================================
C *************************************************************************
C *************************************************************************
C
C
      SUBROUTINE RUSSEL(TEM,PG)
c     =========================
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'

      COMMON/COMFH1/C(600,5),PPMOL(600),APMLOG(600),
     *              XIP(100),CCOMP(100),UIIDUI(100),P(100),
     *              FP(100),XKP(100),EPS,SWITER,
     *              NELEM(5,600),NATO(5,600),MMAX(600),
     *              NELEMX(100),NMETAL,NIMAX
        DIMENSION FX(100),DFX(100),Z(100),PREV(100),WA(50)
C
c       ECONST=4.342945E-1
        ECONST=4.3426E-1
        XKCON=6.667343E-1
        EPSDIE=5.0E-5
        T=5040.4/TEM
        PGLOG=log10(PG)
        tk=1./(tem*1.38054e-16)
C
C    HEH=helium/hydrogen ratio by number
C
        HEH=CCOMP(2)/CCOMP(1)
c       HEH=YTOT(1)-UN
C
C    evaluation of log XKP(MOL)
C
        DO J=1,NMOLEC
           APLOGJ=C(J,5)
           DO K=1,4
              KM5=5-K
              APLOGJ=APLOGJ*T + C(J,KM5)
           END DO
           APMLOG(J)=APLOGJ
        END DO
        apmlog(1)=-log10(1.0353e-16/tem/sqrt(tem)*tk*exp(8762.9/tem))
        DHH=(((0.1196952E-02*T-0.2125713E-01)*T+0.1545253E+00)*T
     *     -0.5161452E+01)*T+0.1277356E+02
        DHH=EXP(DHH/ECONST)
C
C  evaluation of the ionization constants
C
        TEM25=TEM**2*SQRT(TEM)
        DO I=1,NMETAL
           NELEMI = NELEMX(I)
*
* calculation of the partition functions following Irwin (1981)
C
           call mpartf(nelemi,1,0,tem,g0)
           call mpartf(nelemi,2,0,tem,g1)
c          uiidui(nelemi)=g1/g0*0.6665
           uiidui(nelemi)=g1/g0*xkcon
c        
           XKP(NELEMI)=UIIDUI(NELEMI)*TEM25*
     *                 EXP(-XIP(NELEMI)*T/ECONST)
        END DO
        HKP=XKP(1)
C
C   preliminary value of PH at high temperatures
C
        HKP=XKP(1)
        IF(T.LT.0.6) THEN
           PPH=SQRT(HKP*(PG/(1.0+HEH)+HKP))-HKP
           PH=PPH**2/HKP
         ELSE
           IF(PG/DHH.LE.0.1) THEN
              PH=PG/(1.0+HEH)
            ELSE
              PH=0.5 * (SQRT(DHH*(DHH+4.0 *PG/(1.0+HEH)))-DHH)
           END IF
        END IF
C
C  evaluation of the fictitious pressures of hydrogen
C     PG=PH+PHH+2.0*PPH+HEH*(PH+2.0*PHH+PPH)
C
        U=(1.0+2.0*HEH)/DHH
        Q=1.0+HEH
        R=(2.0+HEH)*SQRT(HKP)
        S=-1.0*PG
        X=SQRT(PH)
C
C       Russell iterations
C
        ITERAT=0
   10   CONTINUE
        F=((U*X**2+Q)*X+R)*X+S
        DF=2.0*(2.0*U*X**2+Q)*X+R
        XR=X-F/DF
C
        IF(ABS((X-XR)/XR).GT.EPSDIE) THEN
           ITERAT=ITERAT+1
           IF(ITERAT.GT.50) THEN
              WRITE(6,710) TEM,PG,X,XR,PH
  710    FORMAT(1H1, ' NOT CONVERGE IN RUSSEL '/// 'TEM=',F9.2,5X,'PG=',
     *        E12.5,5X,'X1=',E12.5,5X,'X2=',E12.5,5X,'PH=',E12.5/////)
            ELSE
              X=XR
              GO TO 10
           END IF
        END IF
        PH=XR**2
        PHH=PH**2/DHH
        PPH=SQRT(HKP*PH)
        FPH=PH+2.0*PHH+PPH
        P(100)=PPH
C
C   evaluation of the fictitious pressure of each element
C
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FP(NELEMI)=CCOMP(NELEMI)*FPH
        END DO
C
        PE=P(99)
C
C    Russell equations
C
        NITERR = 0
   20   CONTINUE
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           FX(NELEMI)=-FP(NELEMI)+P(NELEMI)*(1.0+XKP(NELEMI)/PE)
           DFX(NELEMI)=1.0+XKP(NELEMI)/PE
        END DO
C
        SPNION=0.0
        spnplu=0.
        DO J=1,NMOLEC
           MMAXJ=MMAX(J)
           PMOLJL=-APMLOG(J)
           DO M=1,MMAXJ
              NELEMJ=NELEM(M,J)
              NATOMJ=NATO(M,J)
              PMOLJL=PMOLJL+DFLOAT(NATOMJ)*log10(P(NELEMJ))
           END DO
C
           PMOLJ=EXP(PMOLJL/ECONST)
           DO M=1,MMAXJ
              NELEMJ=NELEM(M,J)
              NATOMJ=NATO(M,J)
              ATOMJ=DFLOAT(NATOMJ)
              IF(NELEMJ.EQ.99) then
                 if(natomj.ge.0) then
                    SPNION=SPNION+PMOLJ*NATOMJ
                  else
                    SPNPLU=SPNPLU-PMOLJ*NATOMJ
                 end if
              end if
              DO I=1,NMETAL
                 NELEMI=NELEMX(I)
                 IF(NELEMJ.EQ.NELEMI) THEN
                    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
                    DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*
     *                          PMOLJ/P(NELEMI)
                 END IF
              END DO
           END DO
           PPMOL(J)=PMOLJ
        END DO
C
C   solution of the Russell equations by Newton-Raphson method
C
c      iprin=5
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           WA(I)=log10(P(NELEMI)+1.0D-70)
        END DO
        IMAXP1=NMETAL+1
        WA(IMAXP1)=log10(PE+1.0D-70)
        DELTRS = 0.0
        DO I=1,NMETAL
           NELEMI=NELEMX(I)
           PREV(NELEMI)=P(NELEMI)-FX(NELEMI)/DFX(NELEMI)
           PREV(NELEMI)=ABS(PREV(NELEMI))
           IF(PREV(NELEMI).LT.1.0D-70) PREV(NELEMI)=1.0D-70
           Z(NELEMI)=PREV(NELEMI)/P(NELEMI)
           DELTRS=DELTRS+ABS(Z(NELEMI)-1.0)
           IF(SWITER.GT.0.0) THEN
              P(NELEMI)=(PREV(NELEMI)+P(NELEMI))*0.5
            ELSE
              P(NELEMI)=PREV(NELEMI)
           END IF
        END DO
C
C   ionization equilibrium
C
        PEREV = spnplu
        DO I=1,NMETAL
           NELEMI = NELEMX(I)
           PEREV=PEREV+XKP(NELEMI)*P(NELEMI)
        END DO
C
        PEREV=SQRT(PEREV/(1.0+SPNION/PE))
        DELTRS=DELTRS+ABS((PE-PEREV)/PE)
        if(iprin.gt.4)
     *     write(6,601) niterr,tem,pg*tk,fph*tk,pe*tk,perev*tk,
     *    (perev+pe)*0.5*tk,deltrs
        PE=(PEREV+PE)*0.5
        P(99)=PE
        IF(DELTRS.GT.EPS) THEN
           NITERR=NITERR+1
           IF(NITERR.LE.NIMAX) THEN
              GO TO 20
            ELSE
              WRITE(6,605) NIMAX
           END IF
        END IF
  605   FORMAT(1H0,'*DOES NOT CONVERGE AFTER ',I4,' ITERATIONS')
C
        if(iprin.gt.4) then
          write(6,601) niterr,tem,pg*tk,fph*tk,pe*tk,perev*tk,
     *    (perev+pe)*0.5*tk,deltrs
  601     format('russel iterations ',i4,1p7e13.4)
          write(*,*) ' '
        end if
c
        RETURN
        END

C
C
C ********************************************************************
C
C
c       
      subroutine tiopf(t,pf)
c     ======================
c
c     TiO partition function (data from Kurucz web site)
c
      INCLUDE 'PARAMS.FOR'
      dimension pf0(800)
      data pf0/       
     *    29.107,    55.425,    82.417,   111.190,   142.564,   176.916,
     *   214.340,   254.774,   298.065,   344.021,   392.431,   443.089,
     *   495.795,   550.365,   606.632,   664.449,   723.686,   784.230,
     *   845.981,   908.862,   972.800,  1037.739,  1103.636,  1170.451,
     *  1238.155,  1306.723,  1376.144,  1446.403,  1517.492,  1589.409,
     *  1662.152,  1735.724,  1810.122,  1885.352,  1961.428,  2038.351,
     *  2116.119,  2194.758,  2274.260,  2354.633,  2435.907,  2518.063,
     *  2601.125,  2685.096,  2769.992,  2855.809,  2942.560,  3030.257,
     *  3118.897,  3208.496,  3299.067,  3390.598,  3483.106,  3576.598,
     *  3671.095,  3766.569,  3863.048,  3960.522,  4059.035,  4158.545,
     *  4259.074,  4360.642,  4463.259,  4566.905,  4671.582,  4777.321,
     *  4884.105,  4991.937,  5100.852,  5210.813,  5321.838,  5433.972,
     *  5547.154,  5661.417,  5776.789,  5893.211,  6010.774,  6129.422,
     *  6249.173,  6370.026,  6491.973,  6615.042,  6739.240,  6864.542,
     *  6990.959,  7118.533,  7247.214,  7377.053,  7508.012,  7640.121,
     *  7773.370,  7907.764,  8043.309,  8180.032,  8317.835,  8456.861,
     *  8597.055,  8738.396,  8880.926,  9024.672,  9169.570,  9315.610,
     *  9462.927,  9611.339,  9760.963,  9911.798, 10063.900, 10217.148,
     * 10371.572, 10527.253, 10684.109, 10842.173, 11001.469, 11161.970,
     * 11323.751, 11486.758, 11650.978, 11816.415, 11983.159, 12151.134,
     * 12320.243, 12490.668, 12662.333, 12835.234, 13009.470, 13184.926,
     * 13361.601, 13539.660, 13718.891, 13899.456, 14081.252, 14264.326,
     * 14448.643, 14634.341, 14821.225, 15009.476, 15199.021, 15389.829,
     * 15581.955, 15775.377, 15970.188, 16166.239, 16363.513, 16562.006,
     * 16761.930, 16963.301, 17165.906, 17369.881, 17575.236, 17781.814,
     * 17989.816, 18198.996, 18409.707, 18621.680, 18835.068, 19049.715,
     * 19265.768, 19483.375, 19702.006, 19922.209, 20143.668, 20366.555,
     * 20590.742, 20816.402, 21043.338, 21271.672, 21501.369, 21732.563,
     * 21965.119, 22199.068, 22434.432, 22671.266, 22909.307, 23148.898,
     * 23389.893, 23632.322, 23875.969, 24121.160, 24367.707, 24615.848,
     * 24865.471, 25116.320, 25368.604, 25622.342, 25877.512, 26134.055,
     * 26392.404, 26651.764, 26912.826, 27175.250, 27439.197, 27704.539,
     * 27971.287, 28239.572, 28509.373, 28780.707, 29053.516, 29327.602,
     * 29603.338, 29880.539, 30159.105, 30439.322, 30721.055, 31004.254,
     * 31288.818, 31575.061, 31862.693, 32151.781, 32442.586, 32734.619,
     * 33027.777, 33323.023, 33619.535, 33917.707, 34217.711, 34518.996,
     * 34821.676, 35126.195, 35432.141, 35739.602, 36048.926, 36359.488,
     * 36672.023, 36985.633, 37300.863, 37617.965, 37936.469, 38256.309,
     * 38578.074, 38901.668, 39226.461, 39552.969, 39880.852, 40210.785,
     * 40541.852, 40874.691, 41209.359, 41545.535, 41883.602, 42222.715,
     * 42563.895, 42906.508, 43250.656, 43596.902, 43944.355, 44293.695,
     * 44644.504, 44997.621, 45351.590, 45707.242, 46065.008, 46424.367,
     * 46785.605, 47148.023, 47512.496, 47878.418, 48246.426, 48615.895,
     * 48987.336, 49360.082, 49734.758, 50111.004, 50489.383, 50868.996,
     * 51250.250, 51633.691, 52018.945, 52405.715, 52794.090, 53184.340,
     * 53576.375, 53970.605, 54366.176, 54763.148, 55162.430, 55563.215,
     * 55966.391, 56371.000, 56777.176, 57185.570, 57596.074, 58007.617,
     * 58421.418, 58837.172, 59254.539, 59673.418, 60094.066, 60517.410,
     * 60941.844, 61368.660, 61797.395, 62227.590, 62659.789, 63094.238,
     * 63529.695, 63967.488, 64407.887, 64849.496, 65292.867, 65735.922,
     * 66182.000, 66631.266, 67082.055, 67534.391, 67988.992, 68446.117,
     * 68904.789, 69365.180, 69827.914, 70292.781, 70759.352, 71228.500,
     * 71699.375, 72171.672, 72647.086, 73123.984, 73603.023, 74083.516,
     * 74566.359, 75050.555, 75537.758, 76027.258, 76518.125, 77012.008,
     * 77507.063, 78003.813, 78503.977, 79006.125, 79509.320, 80015.375,
     * 80522.461, 81031.938, 81544.164, 82058.313, 82574.352, 83093.914,
     * 83614.367, 84136.820, 84662.211, 85188.867, 85719.375, 86249.977,
     * 86783.781, 87319.219, 87857.180, 88396.797, 88939.805, 89484.266,
     * 90032.023, 90580.930, 91132.563, 91686.148, 92242.742, 92799.406,
     * 93360.016, 93923.453, 94488.313, 95055.211, 95625.297, 96197.477,
     * 96771.531, 97348.156, 97926.922, 98507.453, 99091.563, 99677.938,
     *100267.234,100856.438,101449.828,102045.750,102643.094,103244.117,
     *103846.969,104450.313,105057.641,105667.188,106279.516,106894.937,
     *107512.789,108133.117,108754.758,109377.687,110005.039,110634.602,
     *111266.141,111902.133,112537.984,113178.891,113819.766,114464.312,
     *115110.969,115760.687,116412.469,117068.055,117724.547,118384.383,
     *119047.469,119712.469,120380.187,121051.336,121724.102,122399.250,
     *123076.266,123756.977,124441.195,125126.406,125816.453,126506.766,
     *127202.367,127899.086,128598.266,129299.969,130004.969,130712.016,
     *131409.266,132117.719,132828.969,133544.016,134262.750,134986.344,
     *135712.891,136439.937,137170.969,137905.562,138641.578,139380.266,
     *140122.937,140868.641,141615.484,142366.703,143123.078,143880.000,
     *144638.484,145401.594,146168.125,146935.359,147707.484,148482.641,
     *149256.578,150037.281,150821.953,151606.750,152396.094,153188.766,
     *153983.391,154782.141,155582.203,156387.234,157192.719,158003.156,
     *158815.125,159632.437,160450.766,161274.750,162098.172,162926.000,
     *163756.609,164593.141,165430.859,166270.937,167114.750,167960.797,
     *168811.562,169663.906,170517.203,171376.531,172239.469,173105.891,
     *173975.250,174847.203,175721.453,176597.250,177480.984,178366.094,
     *179253.828,180145.734,181038.000,181936.031,182837.969,183739.922,
     *184645.937,185558.281,186470.844,187387.422,188307.234,189232.281,
     *190156.000,191088.234,192022.062,192957.250,193899.328,194842.984,
     *195788.391,196736.156,197687.828,198645.719,199603.422,200569.234,
     *201536.437,202508.641,203481.000,204459.016,205438.750,206424.312,
     *207409.953,208398.734,209393.234,210391.047,211390.984,212395.516,
     *213401.547,214420.141,215431.812,216453.453,217476.734,218501.266,
     *219530.219,220560.719,221597.891,222637.875,223677.750,224725.500,
     *225777.406,226829.297,227893.125,228954.547,230020.969,231086.453,
     *232157.469,233233.047,234315.406,235395.625,236480.953,237572.125,
     *238666.484,239765.125,240863.281,241969.750,243079.250,244191.719,
     *245304.812,246427.937,247548.234,248673.562,249804.984,250942.781,
     *252078.953,253222.812,254369.641,255519.359,256671.406,257827.906,
     *258988.859,260154.734,261322.281,262458.781,263606.437,264770.625,
     *265947.750,267125.156,268314.125,269507.687,270702.344,271905.156,
     *273110.156,274318.937,275531.687,276751.344,277970.781,279198.531,
     *280425.750,281663.250,282897.469,284138.906,285383.594,286637.031,
     *287891.156,289147.625,290413.312,291678.719,292946.031,294225.875,
     *295501.344,296782.656,298070.094,299363.875,300652.250,301953.750,
     *303260.062,304563.781,305874.375,307191.437,308517.031,309835.750,
     *311159.375,312490.937,313827.469,315166.781,316511.031,317860.406,
     *319214.969,320565.875,321929.344,323296.906,324660.219,326035.687,
     *327413.844,328794.406,330173.156,331566.156,332953.469,334356.187,
     *335757.625,337165.562,338566.094,339984.750,341402.937,342828.125,
     *344257.562,345686.750,347123.125,348564.250,350008.906,351453.219,
     *352908.062,354361.469,355828.000,357292.500,358765.719,360233.687,
     *361713.562,363200.187,364685.656,366174.500,367673.594,369174.906,
     *370678.969,372191.125,373708.937,375225.281,376743.719,378270.406,
     *379804.500,381334.250,382879.125,384420.812,385969.531,387519.812,
     *389078.937,390639.781,392213.875,393782.437,395359.156,396943.625,
     *398527.625,400110.937,401711.750,403310.344,404908.937,406513.875,
     *408125.781,409741.906,411356.875,412979.500,414613.125,416245.500,
     *417889.094,419530.000,421179.906,422831.531,424484.344,426153.187,
     *427816.406,429489.094,431161.312,432840.656,434517.000,436215.281,
     *437896.000,439602.594,441300.625,443016.156,444722.906,446445.437,
     *448164.812,449885.937,451615.094,453351.594,455090.125,456833.281,
     *458582.719,460335.344,462094.844,463857.094,465629.906,467402.781,
     *469178.406,470963.750,472745.906,474539.594,476333.312,478131.125,
     *479934.000,481740.750,483557.844,485376.625,487202.937,489033.562,
     *490868.031,492709.281,494547.375,496401.094,498249.594,500110.250,
     *501966.594,503836.062,505704.437,507580.687,509469.187,511349.781,
     *513239.000,515137.187,517038.812,518942.906,520858.156,522767.094,
     *524610.625,526433.812,528331.062,530253.437,532185.500,534127.875,
     *536073.937,538028.312,539983.375,541954.687,543916.312,545902.500,
     *547874.812,549857.125,551850.937,553839.937,555836.625,557838.500,
     *559849.937,561859.375,563880.625,565889.875,567916.000,569953.625,
     *571990.375,574034.937,576085.062,578127.375,580188.937,582251.000,
     *584328.812,586385.562,588464.062,590551.875,592644.625,594722.250,
     *596829.937,598931.375,601029.687,603142.812,605262.812,607384.625,
     *609513.125,611644.000,613775.875,615930.375,618073.750,620218.437,
     *622381.937,624524.312,626697.500,628869.000,631040.937,633223.562,
     *635409.187,637597.562,639800.187,642002.125,644212.562,646416.250,
     *648633.562,650864.187,653083.687,655315.312,657549.687,659795.500,
     *662032.250,664292.875,666542.312,668806.250,671071.312,673340.937,
     *675626.938,677898.750/
c
      it=int(t/10.)
      if(it.gt.800) it=800
      pf=pf0(it)
      return
      end

C
C
C ********************************************************************
C
C
      SUBROUTINE SETWIN
C     =================
C
C     Initialisation of an extended radial structure
C      (spherical symmetry is assumed)
C     with a continuous connection between the lower quasi-hydrostatic
C     layers and the upper, supersonic layers. The velocity structure
C     in the upper layers is a beta-type law (v=vinf*(1-r0/r)^beta).
C
C     Additional input are read at the end of Unit 8:
C      RCORE : Core radius (deepest layer, in solar radii or in cm)
C      NDRAD : Number of layers
C      NRCORE: Number of core rays
C      INRV  : Switch indicating the data to be read:
C           = 0 : Read an hydrostatic, plane-parallel model only; the
C                   routine builds the radial points, density and
C                   velocity structure;
C           < 0 : Read also an hydrostatic, plane-parallel model, but
C                   an empirical velocity law V(r) is read at each
C                   radial point (r(id) is read);
C           > 0 : Input from an extended model atmosphere; the velocity
C                   law is read; the density structure is recomputed for
C                   a possibly different mass-loss rate.
C      XMDOT  : Mass loss rate (in solar mass/yr)
C      BETAV, VINF : Parameters of the velocity law (VINF in km/s)
C      RD, VEL: Radial points, expansion velocity
C
C     Synspec version
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (RSUN=6.96D10)
      common/velaux/velmax,iemoff,nltoff,itrad
C
C     Read data for spherical atmosphere and velocity law
C
      READ(8,*,END=9,ERR=9) RCORE,NDRAD,NRCORE,INRV,NFIRY,NDF
      IF(RCORE.LT.1.E5) RCORE=RCORE*RSUN
      IF(NDRAD.GT.MDEPTH) CALL quit('NDRAD too large')
      READ(8,*) XMDOT,BETAV,VINF
      XMDOT=6.30289D25*XMDOT
      VINF=1.D5*VINF
      ND=NDRAD
      DO ID=1,ND
        READ(8,*) RD(ID),VEL(ID),VTURB(ID),DENSCON(ID)
        if(denscon(id).eq.0.) denscon(id)=1.
        vturb(id)=vturb(id)*vturb(id)
      END DO
C
C   Apply density contrast for clumping
C
      DO ID=1,ND
        ELEC(ID) = ELEC(ID) * DENSCON(ID)
        DENS(ID) = DENS(ID) * DENSCON(ID)
        DO I=1,NLEVEL
           POPUL(I,ID) = POPUL(I,ID) * DENSCON(ID)
        END DO
      END DO
C
C  Set up rays and weights
C
      itrad=1
      call radtem
      CALL SETRAY
      CALL WGTJH1
C
    9 continue
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE SETRAY
C     =================
C
C     Setup impact rays and angles
C      (assumes one impact ray tangent to every depth layer)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (PI4=4.*3.141592654)
      PARAMETER (UN=1., TWO=2., HALF=0.5)
      DIMENSION RS(MDEPF ),RDX(MDEPF )
      DIMENSION ZIU(MDEPTH),VIU(MDEPTH),ZIUF(MDEPF ),VIUF(MDEPF )
C
C     Fine radial grid
C
      if(ndf.eq.0.or.ndf.eq.nd) then
         ndf=nd
         DO ID=1,NDF
            DENSF(ID)=DENS(ID)
         END DO
      else
         XR1=LOG(DENS(1))
         XR2=LOG(DENS(ND))
         DXR=(XR2-XR1)/FLOAT(NDF-1)
         DO ID=1,NDF
            DENSF(ID)=EXP(XR1+FLOAT(ID-1)*DXR)
         END DO
      end if
C
C
C     Impact rays
C
      NREXT=ND
      DO ID=1,NREXT
        PIM(ID)=RD(ID)
        NUD(ID)=ID
      END DO
      DO IU=1,NRCORE
        PIM(NREXT+IU)=FLOAT(NRCORE-IU)/FLOAT(NRCORE)*RCORE
        NUD(NREXT+IU)=ND
      END DO
      KMU=NREXT+NRCORE
C
C     Angles
C
      DO ID=1,ND
        RD1=UN/RD(ID)
        DO IU=ID,KMU
          PRR=PIM(IU)*RD1
          BMU(IU,ID)=SQRT(UN-PRR*PRR)
        END DO
      END DO
C
C     Depth increments along each ray
C
      DELZ(1,1)=0.
      DFRQ(1,1)=0.
      DO IU=2,KMU
        NUDF(IU)=NUD(IU)
        IU1=IU
        IF(IU.GT.ND) IU1=ND
        DO ID=1,IU1-1
          DELZ(IU,ID)=BMU(IU,ID)*RD(ID)-BMU(IU,ID+1)*RD(ID+1)
          DFRQ(IU,ID)=BMU(IU,ID)*VEL(ID)/CL
          JD=2*NUD(IU)-ID
          DFRQ(IU,JD)=-DFRQ(IU,ID)
        END DO
        DELZ(IU,IU1)=DELZ(IU,IU1-1)
        DFRQ(IU,IU1)=0.
        IF(IU.GT.NREXT) DFRQ(IU,ND)=BMU(IU,ND)*VEL(ND)/CL
      END DO
C
C Finer grid along the NFIRY most external rays
C   velocity steps DVD(ID)
C
      XMD4=XMDOT/PI4
      CLV=UN/CL
      DO ID=1,ND
        DVD(ID)=SQRT(1.6D7*TEMP(ID)+VTURB(ID)) * 0.3
c        DVD(ID)=SQRT(1.6D7*TEMP(ID))
      END DO
      NUDX=ND
      DO IU=2,NFIRY
        IF(PIM(IU).GT.0.) THEN
          DO ID=1,NUD(IU)
            IID=NUD(IU)-ID+1
            ZIU(ID)=VEL(IID)
            VIU(ID)=DFRQ(IU,IID)*CL
          ENDDO
        ELSE
          DO ID=1,NUD(IU)
            IID=NUD(IU)-ID+1
            ZIU(ID)=RD(IID)
            VIU(ID)=DFRQ(IU,IID)*CL
          ENDDO
        ENDIF
        NUDF(IU)=1
        VIUF(1)=DFRQ(IU,1)*CL
        DO ID=1,NUD(IU)-1
          VZ1=DFRQ(IU,ID)*CL
          VZ2=DFRQ(IU,ID+1)*CL
          NFG=int((VZ1-VZ2)/DVD(ID))+1
          XFG=(VZ1-VZ2)/DFLOAT(NFG)
          IV0=NUDF(IU)
          DO IV=1,NFG
            VIUF(IV0+IV)=VZ1-DFLOAT(IV)*XFG
          ENDDO
          NUDF(IU)=NUDF(IU)+NFG
          IF(NUDF(IU).GT.MDEPF )
     +      CALL quit('Too many points in fine grid - SETRAY')
        END DO
        IF(NUDF(IU).GT.NUDX) NUDX=NUDF(IU)
        INRP=2
        IF(IU.GT.8) INRP=4
        CALL INTERP(VIU,ZIU,VIUF,ZIUF,NUD(IU),NUDF(IU),INRP,0,0)
        IF(PIM(IU).GT.0.) THEN
          DO ID=1,NUDF(IU)
            DMU=VIUF(ID)/ZIUF(ID)
            RS(ID)=PIM(IU)/SQRT(UN-DMU*DMU)
            DFRQF(IU,ID)=VIUF(ID)*CLV
            VELF(IU,ID)=ZIUF(ID)
            RDX(ID)=XMD4/(RS(ID)*RS(ID)*VELF(IU,ID))
            ZIUF(ID)=DMU*RS(ID)
          END DO
        ELSE
          DO ID=1,NUDF(IU)
            RS(ID)=ZIUF(ID)
            DFRQF(IU,ID)=VIUF(ID)*CLV
            VELF(IU,ID)=VIUF(ID)
            RDX(ID)=XMD4/(RS(ID)*RS(ID)*VELF(IU,ID))
          END DO
        END IF
        IF(IU.LE.NREXT) THEN
          DO ID=1,NUDF(IU)
            JD=2*NUDF(IU)-ID
            DFRQF(IU,JD)=-DFRQF(IU,ID)
          END DO
        END IF
        DO ID=1,NUDF(IU)-1
          DELZF(IU,ID)=ZIUF(ID)-ZIUF(ID+1)
        END DO
        DELZF(IU,NUDF(IU))=DELZF(IU,NUDF(IU)-1)
C
C   Assign depth index
C
        KRAY(IU,1)=2
        DRAY(IU,1)=0.
        IDK=1
        DO ID=2,NUDF(IU)
          DO WHILE (RDX(ID).GE.DENSF(IDK).and.idk.le.ndf)
            IDK=IDK+1
          END DO
c          IDK=IDK+1
          IF(IDK.GT.NDF) IDK=NDF
          KRAY(IU,ID)=IDK
          DRAY(IU,ID)=(RDX(ID)-DENSF(IDK-1))/(DENSF(IDK)-DENSF(IDK-1))
        END DO
        IF(IU.LE.NREXT) THEN
          DO ID=1,NUDF(IU)
            JD=2*NUDF(IU)-ID
            KRAY(IU,JD)=KRAY(IU,ID)
            DRAY(IU,JD)=DRAY(IU,ID)
          END DO
        END IF
      END DO
C
C    remaining rays (without finer grid)
C      
      IF(NFIRY.LT.KMU) THEN
        IU=KMU
        KRAY(IU,1)=2
        DRAY(IU,1)=0.
        IDK=1
        DO ID=2,NUDF(IU)
          DO WHILE (DENS(ID).GE.DENSF(IDK).and.idk.le.ndf)
            IDK=IDK+1
          END DO
c          IDK=IDK+1
          IF(IDK.GT.NDF) IDK=NDF
          KRAY(IU,ID)=IDK
          DRAY(IU,ID)=(DENS(ID)-DENSF(IDK-1))/(DENSF(IDK)-DENSF(IDK-1))
        END DO
        DO IU=NFIRY+1,KMU
          DO ID=1,NUDF(IU)
            KRAY(IU,ID)=KRAY(KMU,ID)
            DRAY(IU,ID)=DRAY(KMU,ID)
            DFRQF(IU,ID)=DFRQ(IU,ID)
            DELZF(IU,ID)=DELZ(IU,ID)
          ENDDO
          IF(IU.LE.NREXT) THEN
            DO ID=1,NUDF(IU)
              JD=2*NUDF(IU)-ID
              KRAY(IU,JD)=KRAY(IU,ID)
              DRAY(IU,JD)=DRAY(IU,ID)
              DFRQF(IU,JD)=-DFRQF(IU,ID)
            END DO
          END IF
        END DO
      END IF
C
      NFTOT=0
      DO IU=2,KMU
        IUD=NUDF(IU)
        IF(IU.LE.NREXT) IUD=2*NUDF(IU)-1
        NFTOT=NFTOT+IUD
      ENDDO
      write(10,*) 'NFTOT=',NFTOT
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE WGTJH1
C     =================
C
C     Angle quadrature weights
C      from Hummer, Kunasz, & Kunasz, 1973, Comp. Phys. Comm. 6, 38
C
C     The present version of this routine assumes that there are
C      impact rays tangent to every depth layers (i.e. NREXT=ND)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (UN=1., TWO=2., HALF=0.5)
      PARAMETER (SIX=6.)
      PARAMETER (C03=UN/3.,D03=2./3.,C04=UN/4.,C06=UN/6.)
      PARAMETER (C24=UN/24.,C45=UN/45.,D45=2./45.,C72=UN/72.)
      DIMENSION WAJ(MKU),WBJ(MKU),AHH(MKU,4)
      DIMENSION BMUH(MKU),BMUHP(MKU),WAH(MKU),WBH(MKU)
      DIMENSION WSD(MKU),WSU(MKU),WSL(MKU),WUU(MKU)
      DIMENSION WTD(MKU),WTU(MKU),WTL(MKU)
C
      DO 100 ID=1,ND
        DO IU=ID+1,KMU
          AHH(IU,1)=BMU(IU,ID)-BMU(IU-1,ID)
          AHH(IU,2)=AHH(IU,1)*AHH(IU,1)
          AHH(IU,3)=AHH(IU,2)*AHH(IU,1)
          AHH(IU,4)=AHH(IU,3)*AHH(IU,1)
          BMUH(IU)=BMU(IU,ID)*AHH(IU,1)
          BMUHP(IU)=BMU(IU-1,ID)*AHH(IU,1)
        END DO
C
C     Weights for J
C
        WAJ(ID)=HALF*AHH(ID+1,1)
        WAJ(KMU)=HALF*AHH(KMU,1)
        WBJ(ID)=-C24*AHH(ID+1,3)
        WBJ(KMU)=-C24*AHH(KMU,3)
        WSL(ID+1)=C06*AHH(ID+1,1)
        WSU(KMU-1)=0.
        WSD(ID)=C03*AHH(ID+1,1)
        WSD(KMU)=UN
        WTL(ID+1)=UN/AHH(ID+1,1)
        WTU(KMU-1)=0.
        WTD(ID)=-WTL(ID+1)
        WTD(KMU)=0.
        DO IU=ID+1,KMU-1
          WAJ(IU)=HALF*(AHH(IU,1)+AHH(IU+1,1))
          WBJ(IU)=-C24*(AHH(IU+1,3)+AHH(IU,3))
          AH1=SIX/(AHH(IU,1)+AHH(IU+1,1))
          WSL(IU+1)=C06*AH1*AHH(IU+1,1)
          WSU(IU-1)=UN-WSL(IU+1)
          WSD(IU)=TWO
          WTL(IU+1)=AH1/AHH(IU+1,1)
          WTU(IU-1)=AH1/AHH(IU,1)
          WTD(IU)=-SIX/AHH(IU,1)/AHH(IU+1,1)
        END DO
        NMUD=KMU-ID+1
        CALL TRIDAG(WSL,WSD,WSU,WBJ,WUU,NMUD)
        WMUJ(ID,ID)=WAJ(ID)+WTD(ID)*WUU(ID)+WTU(ID)*WUU(ID+1)
        WMUJ(KMU,ID)=WAJ(KMU)+WTL(KMU)*WUU(KMU-1)+WTD(KMU)*WUU(KMU)
        DO IU=ID+1,KMU-1
          WMUJ(IU,ID)=WAJ(IU)+WTL(IU)*WUU(IU-1)+
     *                  WTD(IU)*WUU(IU)+WTU(IU)*WUU(IU+1)
        END DO
C
C     Weights for emergent flux H
C
      IF(ID.GT.1) GO TO 100
      WAH(ID)=HALF*BMUH(ID+1)-C03*AHH(ID+1,2)
      WAH(KMU)=HALF*BMUHP(KMU)+C03*AHH(KMU,2)
      WBH(ID)=AHH(ID+1,3)*(C45*AHH(ID+1,1)-C24*BMU(ID+1,ID))
      WBH(KMU)=-AHH(KMU,3)*(C45*AHH(KMU,1)+C24*BMU(KMU-1,ID))
      WSL(ID+1)=0.
      WSD(ID)=UN
      WTL(ID+1)=0.
      WTD(ID)=0.
      DO IU=ID+1,KMU-1
        WAH(IU)=HALF*(BMUH(IU+1)+BMUHP(IU))-
     *             C03*(AHH(IU+1,2)-AHH(IU,2))
        WBH(IU)=-C24*(BMUH(IU+1)*AHH(IU+1,2)+BMUHP(IU)*AHH(IU,2))+
     *             C45*(AHH(IU+1,4)-AHH(IU,4))
      END DO
      CALL TRIDAG(WSL,WSD,WSU,WBH,WUU,NMUD)
      WMUH(ID)=WAH(ID)+WTD(ID)*WUU(ID)+WTU(ID)*WUU(ID+1)
      WMUH(KMU)=WAH(KMU)+WTL(KMU)*WUU(KMU-1)+WTD(KMU)*WUU(KMU)
      DO IU=ID+1,KMU-1
        WMUH(IU)=WAH(IU)+WTL(IU)*WUU(IU-1)+
     *             WTD(IU)*WUU(IU)+WTU(IU)*WUU(IU+1)
      END DO
C
  100 CONTINUE
C
C     Weights for H are overwritten by trapezoidal weigths
C
      id=1
      wmuh(1)=bmu(1,id)*(bmu(2,id)-bmu(1,id))*half
      wmuh(kmu)=bmu(kmu,id)*(bmu(kmu,id)-bmu(kmu-1,id))*half
      do iu=2,kmu-1
        wmuh(iu)=bmu(iu,id)*(bmu(iu+1,id)-bmu(iu-1,id))*half
      end do
c
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE TRIDAG(A,B,C,R,U,N)
C     ==============================
C
C     Solve tridiagonal system of equations
C      from Numerical Recipes (standard Gaussian elimination)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'WINCOM.FOR'
      DIMENSION A(N),B(N),C(N),R(N),U(N)
      DIMENSION GTRID(MKU)
C
      BTRID=B(1)
      U(1)=R(1)/BTRID
      DO J=2,N
        GTRID(J)=C(J-1)/BTRID
        BTRID=B(J)-A(J)*GTRID(J)
        U(J)=(R(J)-A(J)*U(J-1))/BTRID
      ENDDO
      DO J=N-1,1,-1
        U(J)=U(J)-GTRID(J+1)*U(J+1)
      ENDDO
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE RESOLW
C     =================
C
C     driver for evaluating opacities and emissivities which then
C     enter the solution of the radiative transfer equation (RTEWIN)
C     Setup opacities for a given frequency set 
C     Oversample in radial and frequency space for later interpolation
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (UN=1., TWO=2., HALF=0.5)
      DIMENSION CROSS(MCROSS,MOPAC),
     *          ABSO(MOPAC),EMIS(MOPAC),
     *          ABSOC(MFREQC),EMISC(MFREQC),SCATC(MFREQC)
      DIMENSION ABSD(MDEPTH),ASF(MDEPF),XDS(MDEPTH),XDSF(MDEPF)
      COMMON/CONOPA/CHC(MFREQC,MDEPTH),ETC(MFREQC,MDEPTH),
     *              SCC(MFREQC,MDEPTH)
      COMMON/HPOPST/HPOP
      COMMON/COPAC/AB(MOPAC,MDEPF),STH(MOPAC,MDEPF),SCH(MFREQC,MDEPF)
      COMMON/LIMPAR/ALAM0,ALAM1,FRMIN,FRLAST,FRLI0,FRLIM
      COMMON/BLAPAR/RELOP,SPACE0,CUTOF0,TSTD,DSTD,ALAMC
      COMMON/FRQSET/IFRS,NFRS
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
C
C     set up the partial line list for the current interval
C
      CALL INISET
C
C     output of information about selected lines
C
      IF(IMODE.LT.2) CALL INIBLA
C
C  Setup fine grid of frequencies
C
      CLV=UN/2.997925E10
      FQ1=FREQ(1)*(UN+VINF*CLV)
      FQ2=FREQ(NFREQ)*(UN-VINF*CLV)
      VXD=SQRT(0.3e7*TSTD)*FREQ(1)*CLV
      VXS=SPACE0*FREQ(1)*FREQ(1)*CLV*1.e-7
c     DVX=MAX(VXD,VXS)
      DVX=VXS
      NOPAC=int((FQ1-FQ2)/DVX)+1
      DVX=(FQ1-FQ2)/DFLOAT(NOPAC)
      NOPAC=NOPAC+3
      nopac=nfreq
      WRITE(6,600) NOPAC,NDF
      IF(NOPAC.GT.MOPAC) CALL quit('Too many freqs in fine grid')
      DO IJ=1,NOPAC
         FFQ(ij)=FQ1-DFLOAT(ij-1)*DVX
c         freq(ij)=ffq(ij)
c         wlam(ij)=2.997925e18/freq(ij)
         fr=freq(ij)*1.d-15
         BNUE(IJ)=BN*fr*fr*fr
         DO IJCI=IJC,NFREQC-1
            IF(WLAM(IJ).LE.WLAMC(IJCI)) GO TO 248
         END DO
  248    CONTINUE
         IJC=IJCI
         IJCINT(IJ)=MAX(IJC-1,1)
         IJCI=IJCINT(IJ)
         FRX1(IJ)=(FREQ(IJ)-FREQC(IJCI+1))/
     *            (FREQC(IJCI)-FREQC(IJCI+1))
c         write(80,681) ij,ijci,wlam(ij),wlamc(ijci),freq(ij),frx1(ij)
c  681 format(2i5,2f10.3,1p2e11.3)
      END DO
      nfreq=nopac
      DO JI=1,NOPAC-1
        FFQV(JI)=UN/(FFQ(JI)-FFQ(JI+1))
      END DO
      FFQV(NOPAC)=UN
c
c     the continuum opacities and radiation field - done only once
c
c     -----------------------------------
      if(iblank.le.1) then
C
c     determine the "core" radius and the factor that multiplies
c     H_nu at ID=1 to get physical flux there (R2F)
c
      ID0=ND
      DO WHILE(TEMP(ID0).GT.TEFF .AND. ID0.GT.1)
        ID0=ID0-1
      END DO
      ID0=ID0+1
      R2F=RD(1)*RD(1)/RD(ID0)/RD(ID0)
c
C     photoinization cross-sections
C
      CALL CROSEW(CROSS)
C
C     store opacity and emissivity in continuum
C
      DO ID=1,ND
         CALL OPACW(ID,CROSS,ABSO,EMIS,ABSOC,EMISC,SCATC,0)
         DO IJ=1,NFREQC
            CHC(IJ,ID)=ABSOC(IJ) / DENSCON(ID)
            ETC(IJ,ID)=EMISC(IJ) / DENSCON(ID)
            SCC(IJ,ID)=(SCATC(IJ)+ELEC(ID)*SIGE) / DENSCON(ID)
         END DO
      END DO
C
c     radiation field in the continuum
c
      call rtesca
      do ij=1,nfreqc
         write(17,640) wlamc(ij),fluxc(ij)*r2f
      end do
  640 FORMAT(1H ,F10.4,1PE15.5)
c
      end if
c     -----------------------------------
C
C     Store opacity and thermal source function in all frequencies 
C     and depths
C
      DO ID=1,ND
         CALL OPACW(ID,CROSS,ABSO,EMIS,ABSOC,EMISC,SCATC,1)
         DO IJ=1,NOPAC
            AB(IJ,ID)=ABSO(IJ) / DENSCON(ID)
            STH(IJ,ID)=EMIS(IJ)/ABSO(IJ)
         END DO
      END DO
C
c      do id=1,nd
c      do ij=1,nopac
c         write(92,693) id,ij,wlam(ij),ab(ij,id),sth(ij,id)
c      end do
c      end do
c  693 format(2i5,f10.3,1p2e10.3)
C
C  Interpolate to a finer radial (density) grid
C
      if(ndf.ne.nd)  then
      DO ID=1,ND
        XDS(ID)=LOG10(DENS(ID))
      END DO
      DO ID=1,NDF
        XDSF(ID)=LOG10(DENSF(ID))
      END DO
      DO IJ=1,NOPAC
        DO ID=1,ND
          ABSD(ID)=AB(IJ,ID)
        END DO
        CALL INTERP(XDS ,ABSD,XDSF ,ASF,ND,NDF,2,0,1)
        DO ID=1,NDF
          AB(IJ,ID)=ASF(ID)
        END DO
        DO ID=1,ND
          ABSD(ID)=STH(IJ,ID)
        END DO
        CALL INTERP(XDS ,ABSD,XDSF ,ASF,ND,NDF,2,0,1)
        DO ID=1,NDF
          STH(IJ,ID)=ASF(ID)
        END DO
      END DO
      DO IJ=1,NFREQC
        DO ID=1,ND
          ABSD(ID)=SCC(IJ,ID)
        END DO
        CALL INTERP(XDS ,ABSD,XDSF ,ASF,ND,NDF,2,0,1)
        DO ID=1,NDF
          SCH(IJ,ID)=ASF(ID)
        END DO
      END DO
      end if
      WRITE(6,601)
  600 FORMAT(/,' Opacity table for',i5,' frequencies and',/,
     *         '                  ',i5,' radial (density) points')
  601 FORMAT(' Done'/)
C
c      if(ndf.ne.nd)  then
c      do id=1,ndf
c      do ij=1,nopac
c         write(91,693) id,ij,ab(ij,id),sth(ij,id)
c      end do
c      end do
c  693 format(2i5,1p2e11.3)
c      end if
C
C
C     Loop on rays, solving radiative transfer equation
C
      DO IJ=1,NFREQ
        FLUX(IJ)=0.
      END DO
      DO IU=2,KMU
        CALL RTEWIN(IU)
      END DO
      DO IJ=1,NFREQ
        FLUX(IJ)=FLUX(IJ)*R2F
      END DO
C
      RETURN
      END
C
C
C     ****************************************************************
C
C
      SUBROUTINE RTESCA
C     =================
C
C     Solution of the radiative transfer equation
C      for deriving the scattering in continuum
C
C     Solution along every rays, for the spherically-symmetric case
C
C     Solution in the optical depth scale
C
C     The numerical method used:
C     Discontinuous Finite Element method
C     Castor, Dykema, Klein, 1992, ApJ 387, 561.
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (UN=1., TWO=2., HALF=0.5)
      PARAMETER (NTRALI=10,DJMAX=1.D-3)
      COMMON/RTEOPA/CH(MFREQ,MDEPTH),ET(MFREQ,MDEPTH),
     *              SC(MFREQ,MDEPTH)
      COMMON/CONOPA/CHC(MFREQC,MDEPTH),ETC(MFREQC,MDEPTH),
     *              SCC(MFREQC,MDEPTH)
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
      COMMON/CONSCV/SCCF(MFREQC,mdepf)
      DIMENSION ST0(mdepf ),RAD00(mdepf ),AB0(mdepf ),ALI1(mdepf ),
     *          rip(mdepf ),rim(mdepf ),riin(mdepf ),riup(mdepf ),
     *          aip(mdepf ),aim(mdepf ),aiin(mdepf ),aiup(mdepf )
      dimension dt(mdepf ),dtau(mdepf ),RDX(mdepf ),PTX(mdepf )
      dimension uf(mdepf ),af(mdepf ),ss0(mdepf ),scx(mdepth)
      dimension densr(mdepf),rdy(mdepf),
     *          abc0(mdepf),abc1(mdepf),stc0(mdepf),stc1(mdepf),
     *          scc0(mdepf),scc01(mdepf)
      COMMON/COPAC/AB(MOPAC,MDEPF),STH(MOPAC,MDEPF),SCH(MFREQC,MDEPF)
C
C     overall loop over continuum frequencies
C
      DO 500 IJ=1,NFREQC       
      FR=FREQC(IJ)
C
C     Initialisation of J=B
C
      if(ij.eq.1) then
      FR15=FR*1.D-15
      BNU=BN*FR15*FR15*FR15
      HKFR=HK*FR
      DO ID=1,ND
        RAD00(ID)=BNU/(EXP(HKFR/TEMP(ID))-UN)
      END DO
      end if
C
C     Loop over electron scattering
C
      itrali=0
   10 itrali=itrali+1
      fluxc(ij)=0.
C
      DO ID=1,ND
        RAD1(ID)=0.
        ALI1(ID)=0.
      END DO
C
C     Loop over impact rays
C
      if(nd.eq.ndf) then
         do id=1,nd
            densf(id)=dens(id)
            rdx(id)=rad00(id)
            abc0(id)=chc(ij,id)
            stc0(id)=etc(ij,id)/chc(ij,id)
            scc0(id)=scc(ij,id)
         end do
       else
         CALL INTERP(DENS,RAD00,DENSF,RDX,ND,NDF,4,1,0)
         do id=1,nd
            abc1(id)=chc(ij,id)
            stc1(id)=etc(ij,id)/chc(ij,id)
            scc01(ij)=scc(ij,id)
         end do
         CALL INTERP(DENS,abc1,DENSF,abc0,ND,NDF,4,1,0)
         CALL INTERP(DENS,stc1,DENSF,stc0,ND,NDF,4,1,0)
         CALL INTERP(DENS,scc01,DENSF,scc0,ND,NDF,4,1,0)
      end if
      DO 100 IU=1,KMU
        iud=nud(iu)
        IF(IU.LE.NFIRY) IUD=NUDF(IU)
        if(iud.le.1) goto 100
        DO ID=1,IUD
          KY=KRAY(IU,ID)
          YDR=DRAY(IU,ID)
          YDR1=UN-DRAY(IU,ID)
          DENSR(ID)=YDR1*DENSF(KY-1)+YDR*DENSF(KY)
          AB0(ID)=YDR1*abc0(KY-1)+YDR*abc0(KY)
          ST0(ID)=YDR1*stc0(KY-1)+YDR*stc0(KY)
          SC0=YDR1*scc0(KY-1)+YDR*scc0(KY)
          RDY(id)=YDR1*RDX(KY-1)+YDR*RDX(KY)
          SS0(ID)=SC0/AB0(ID)
          ST0(ID)=ST0(ID)+SS0(ID)*RDY(ID)
        END DO
        IF(IU.LE.NFIRY) THEN
          DO ID=1,IUD-1
            DTAU(ID)=HALF*(AB0(ID)+AB0(ID+1))*DELZF(IU,ID)
          END DO
        ELSE
          DO ID=1,IUD-1
            DT(ID)=HALF*(AB0(ID)+AB0(ID+1))
            DTAU(ID)=DT(ID)*DELZ(IU,ID)
          END DO
        END IF
C               
C        incoming intensity   (TAUMIN=0.)
C   
        rim(1)=0.
        aim(1)=0.
        do id=1,iud-1
          dt0=dtau(id)
          dtaup1=dt0+un
          dtau2=dt0*dt0
          bb=two*dtaup1
          cc=dt0*dtaup1
          aa=un/(dtau2+bb)
          rip(id)=(bb*rim(id)+cc*st0(id)-dt0*st0(id+1))*aa
          rim(id+1)=(two*rim(id)+dt0*st0(id)+cc*st0(id+1))*aa
          aip(id)=(cc+bb*aim(id))*aa
          aim(id+1)=cc*aa
        enddo
        do id=2,iud-1
          dtt=un/(dtau(id-1)+dtau(id))
          riin(id)=(rim(id)*dtau(id)+rip(id)*dtau(id-1))*dtt
          aiin(id)=(aim(id)*dtau(id)+aip(id)*dtau(id-1))*dtt
        enddo
        riin(1)=rim(1)
        riin(iud)=rim(iud)
        aiin(1)=aim(1)
        aiin(iud)=aim(iud)
        rip(iud)=rim(iud)
C
C        Outgoing intensity
C          symmetric boundary condition (rim(iud)=riin(iud))
C       or diffusion approx. for core rays
C 
        IF(IU.GT.NREXT) THEN
          PLAND=BNU/(EXP(HKFR/TEMP(ND))-UN)
          DPLAN=PLAND-BNU/(EXP(HKFR/TEMP(ND-1))-UN)
c         rim(iud)=PLAND+dplan/dtau(iud-1)
          rip(iud)=PLAND+dplan/dtau(iud-1)
          dt0=dtau(iud-1)
          dtaup1=dt0+un
          dtau2=dt0*dt0
          bb=two*dtaup1
          cc=dt0*dtaup1
          aa=dtau2+bb
          rim(iud)=(aa*rip(iud)-cc*st0(iud)+dt0*st0(iud-1))/bb
        ENDIF
        do id=iud-1,1,-1
          dt0=dtau(id)
          dtaup1=dt0+un
          dtau2=dt0*dt0
          bb=two*dtaup1
          cc=dt0*dtaup1
          aa=un/(dtau2+bb)
          rip(id+1)=(bb*rim(id+1)+cc*st0(id+1)-dt0*st0(id))*aa
          rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa
          aip(id+1)=(cc+bb*aim(id+1))*aa
          aim(id)=cc*aa
        enddo
        do id=2,iud-1
          dtt=un/(dtau(id-1)+dtau(id))
          riup(id)=(rim(id)*dtau(id-1)+rip(id)*dtau(id))*dtt
          aiup(id)=(aim(id)*dtau(id-1)+aip(id)*dtau(id))*dtt
        enddo
        riup(1)=rim(1)
        riup(iud)=rim(iud)
        aiup(1)=aim(1)
        aiup(iud)=aim(iud)
C
C      symmetrized (Feautrier) intensity  -- (riin+riup)/2 --
C        and interpolation in original radial grid
C
        do id=1,iud
          uf(id)=(riup(id)+riin(id))
          af(id)=(aiup(id)+aiin(id))
        end do
        if(iu.le.nfiry) then
          inrp=min(nud(iu),4)
          call interp(densr,uf,dens,ptx,iud,nud(iu),inrp,1,0)
          do id=1,nud(iu)
            uf(id)=ptx(id)
          end do
          call interp(densr,af,dens,ptx,iud,nud(iu),inrp,1,0)
          do id=1,nud(iu)
            af(id)=ptx(id)
          end do
          iud=nud(iu)
        end if
C
C       Contribution to J
C
        do id=1,nud(iu)
          rad1(id)=rad1(id)+wmuj(iu,id)*uf(id)
          ali1(id)=ali1(id)+wmuj(iu,id)*af(id)
        end do
        FLUXc(IJ)=FLUXc(IJ)+WMUH(IU)*RIM(1)
C
C    End loop over impact rays
C
  100 CONTINUE
C
C     solution of the transfer equation
C     Variables:
C     RAD1    - mean intensity
C
      NDX=NUDF(KMU)
      CALL INTERP(DENSR,SS0,DENS,SCX,NDX,ND,4,1,1)
      DJTOT=0.
      DO ID=1,ND
        RAD1(ID)=RAD1(ID)*HALF
        ALI1(ID)=ALI1(ID)*HALF
        SSS=SCX(ID)
c        DELTAJ=(UN+SSS*ALI1(ID))*(RAD1(ID)-RAD00(ID))
        DELTAJ=(RAD1(ID)-RAD00(ID))/(UN-SSS*ALI1(ID))
c        DELTAJ=RAD1(ID)-RAD00(ID)
        RAD00(ID)=RAD00(ID)+DELTAJ
        DJTOT=MAX(DJTOT,ABS(DELTAJ/RAD00(ID)))
      END DO
      write(6,1600) ij,2.997925e18/fr,itrali,djtot,djmax
      IF(DJTOT.GT.DJMAX.AND.ITRALI.LE.NTRALI) GO TO 10
 1600 format(' IJ,LAM,ITRALI,DJ',i5,f10.2,i5,1p2e12.3)
C
C     end loop for electron scattering
C
      CALL INTERP(DENS,RAD00,DENSF,RDX,ND,NDF,4,1,0)
      do id=1,ndf
        sccf(ij,id)=scc0(ID)*RDX(ID)
      enddo
      fluxc(ij)=fluxc(ij)*2.997925e18/wlamc(ij)**2*0.5
C
  500 CONTINUE
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE RTEWIN(IU)
C     =====================
C
C     Solution of the radiative transfer equation - frequency by
C     frequency - for the known source function.
C
C     The numerical method used:
c     Discontinuous Finite Element (DFE) method
c     Castor, Dykema, Klein, 1992, ApJ 387, 561.
C
C     Input through blank COMMON block:
C      AB     - two-dimensional array  absorption coefficient (frequency,
C               depth)
C      STH     - Thermal source function
C
C     Version including velocity field and extension
C      radiative transfer along ray IU
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'WINCOM.FOR'
      PARAMETER (UN=1., TWO=2., HALF=0.5)
      PARAMETER (TAUREF = 0.6666666666667)
      DIMENSION ST0(2*MDEPF ),TAU(2*MDEPF ),AB0(2*MDEPF ),
     *          rip(2*MDEPF ),rim(2*MDEPF )
c     dimension sc0(2*mdepf)
      dimension sctd(2*mdepf)
      COMMON/COPAC/AB(MOPAC,MDEPF),STH(MOPAC,MDEPF),SCH(MFREQC,MDEPF)
      COMMON/EMFLUX/FLUX(MFREQ),FLUXC(MFREQC)
      COMMON/CONSCV/SCCF(MFREQC,mdepf)
      COMMON/REFDEP/IREFD(MFREQ)
C
      IUD=NUDF(IU)
      IF(IU.LE.NREXT) IUD=2*NUDF(IU)-1
      IF(IUD.EQ.1) RETURN
      IF(NFREQ.GT.1) dlama0=(wlobs(nfrobs)-wlobs(1))/(nfrobs-1)
C
C     overall loop over frequencies (observer's frame)
C
      DO 500 IJ=1,NFROBS
      FR=FRQOBS(IJ)
      wl0=wlobs(ij)
C
C    Opacity and total source function 
c    interpolation in opacity table
C
      IVK=NOPAC-2
      DO ID=1,IUD
         KY=KRAY(IU,ID)
         YDR=DRAY(IU,ID)
         YDR1=UN-YDR
         dwlcom=wl0*DFRQF(IU,ID)
         wlcom=wl0+dwlcom
         if(wlcom.le.wlam(3)) then
            abd1=ab(1,ky-1)
            std1=sth(1,ky-1)
            abd0=ab(1,ky)
            std0=sth(1,ky)
            ij1=1
          else if(wlcom.ge.wlam(nfreq)) then
            abd1=ab(nfreq,ky-1)
            std1=sth(nfreq,ky-1)
            abd0=ab(nfreq,ky)
            std0=sth(nfreq,ky)
            ij1=nfreq
          else
            xijap=(wlcom-wlam(3))/dlama0
            ijap=int(xijap)
            ijap=max(ijap,1)
            ijap=min(ijap,nfreq)
            wlap=wlam(ijap)
            if(wlcom.lt.wlap) then
               ij1=ijap-1
               do iji=ijap-1,1,-1
                  if(wlcom.ge.wlam(iji)) go to 20
               end do
   20          continue
               ij1=iji
             else
               ij1=ijap+1
               do iji=ijap+1,nfreq
                  if(wlcom.lt.wlam(iji)) go to 30
               end do
   30          continue
               ij1=iji-1
            end if
            xfa=(wlam(ij1+1)-wlcom)/(wlam(ij1+1)-wlam(ij1))
            abd1=xfa*ab(ij1,ky-1)+(1.-xfa)*ab(ij1+1,ky-1)
            std1=xfa*sth(ij1,ky-1)+(1.-xfa)*sth(ij1+1,ky-1)
            abd0=xfa*ab(ij1,ky)+(1.-xfa)*ab(ij1+1,ky)
            std0=xfa*sth(ij1,ky)+(1.-xfa)*sth(ij1+1,ky)
         end if
        AB0(ID)=YDR1*Abd1+YDR*abd0
        ST0(ID)=YDR1*Std1+YDR*Std0
C
C   Add scattering
C
        IJC=IJCINT(IJ1)
        IF(IFREQ.NE.17) THEN
          SC1=YDR1*SCCF(ijc,KY-1)+YDR*SCCF(ijc,KY)
          SC2=YDR1*SCCF(ijc+1,KY-1)+YDR*SCCF(ijc+1,KY)
          SCT=FRX1(ij1)*SC1+(1.-FRX1(ij1))*SC2
          sctd(id)=sct/ab0(id)
          ST0(ID)=ST0(ID)+SCT/AB0(ID)
        END IF
      ENDDO
C
C    Optical depth scale
C
      TAU(1)=0.
      IREF=1
      IF(IU.LE.NFIRY) THEN
         DO ID=1,IUD-1
            JD=ID
            IF(ID.GT.NUDF(IU)) JD=2*NUDF(IU)-ID-1
            DT=HALF*(AB0(ID+1)+AB0(ID))*DELZF(IU,JD)
            TAU(ID+1)=TAU(ID)+DT
         END DO
       ELSE
         DO ID=1,IUD-1
            JD=ID
            IF(ID.GT.NUD(IU)) JD=2*NUD(IU)-ID-1
            DT=HALF*(AB0(ID+1)+AB0(ID))*DELZ(IU,JD)
            TAU(ID+1)=TAU(ID)+DT
         END DO
      END IF
      if(iu.eq.kmu) then
         DO ID=1,IUD-1
            IF(TAU(ID).LE.TAUREF.AND.TAU(ID+1).GT.TAUREF) IREF=ID
         END DO
         irefd(ij)=iref
      end if
C
C     Outgoing intensity
C 
      IF(IU.LE.NREXT) THEN
C
C     1. External rays
C
        ndt=iud
        rip(ndt)=0.
        dt0=tau(ndt)-tau(ndt-1)
        dtaup1=dt0+un
        dtau2=dt0*dt0
        bb=two*dtaup1
        cc=dt0*dtaup1
        aa=dtau2+bb
        rim(ndt)=(aa*rip(ndt)-cc*st0(ndt)+dt0*st0(ndt-1))/bb
        do id=1,iud-1
          jd=iud-id
          dt0=tau(jd+1)-tau(jd)
          dtaup1=dt0+un
          dtau2=dt0*dt0
          bb=two*dtaup1
          cc=dt0*dtaup1
          aa=un/(dtau2+bb)
          rim(jd)=(two*rim(jd+1)+dt0*st0(jd+1)+cc*st0(jd))*aa
        enddo
      ELSE
C
C      2. core rays
C
        NDT=IUD
        FR15=FR*1.D-15
        BNU=BN*FR15*FR15*FR15
        PLAND=BNU/(EXP(HK*FR/TEMP(ND))-UN)
        DPLAN=BNU/(EXP(HK*FR/TEMP(ND-1))-UN)
        DPLAN=(PLAND-DPLAN)/(TAU(IUD)-TAU(IUD-1))
        RIP(NDT)=PLAND+DPLAN
        dt0=tau(ndt)-tau(ndt-1)
        dtaup1=dt0+un
        dtau2=dt0*dt0
        bb=two*dtaup1
        cc=dt0*dtaup1
        aa=dtau2+bb
        rim(ndt)=(aa*rip(ndt)-cc*st0(ndt)+dt0*st0(ndt-1))/bb
        do id=iud-1,1,-1
           dt0=tau(id+1)-tau(id)
           dtaup1=dt0+un
           dtau2=dt0*dt0
           bb=two*dtaup1
           cc=dt0*dtaup1
           aa=un/(dtau2+bb)
           rim(id)=(two*rim(id+1)+dt0*st0(id+1)+cc*st0(id))*aa
        enddo
      ENDIF
      FLUX(IJ)=FLUX(IJ)+WMUH(IU)*RIM(1)
c
c      if(ij.eq.1.or.ij.eq.3.or.ij.eq.5.or.ij.eq.9.or.ij.eq.83) then
c      if(iu.eq.2.or.iu.eq.20.or.iu.eq.60.or.iu.eq.80) then
c      do id=1,iud
c         write(79,679) ij,iu,id,ab0(id),st0(id),sctd(id),
c     *                 tau(id),rim(id),
c     *                 flux(ij)
c      end do
c      end if
c      end if
c  679 format(3i5,1p6e12.4)
C
c      CFX=WMUH(IU)*RIM(1)
c      write(78,780) ij,iu,wlobs(ij),cfx,RIM(1)
c  780 format(2i4,f10.3,1p2e16.8)
C
c     if(iflux.ge.1) then
C
C     output of emergent specific intensities to Unit 10 (line points)
C     or 18 (two continuum points)
C
c     IF(IJ.GT.2) THEN
c     WRITE(10,618) WLAM(IJ),FLUX(IJ),RIM(1),IU
c     ELSE
c     WRITE(18,618) WLAM(IJ),FLUX(IJ),RIM(1),IU
c     END IF
c     end if
c 618 FORMAT(1H ,f10.3,2pe15.5,i5)
C
C     if needed (if iprin.ge.3), output of interesting physical
C     quantities at the monochromatic optical depth  tau(nu)=2/3
C
c     IF(IPRIN.GE.3) THEN
c     T0=LOG(TAU(IREF+1)/TAU(IREF))
c     X0=LOG(TAU(IREF+1)/TAUREF)/T0
c     X1=LOG(TAUREF/TAU(IREF))/T0
c     DMREF=EXP(LOG(DM(IREF))*X0+LOG(DM(IREF+1))*X1)
c     TREF=EXP(LOG(TEMP(IREF))*X0+LOG(TEMP(IREF+1))*X1)
c     STREF=EXP(LOG(ST0(IREF))*X0+LOG(ST0(IREF+1))*X1)
c     SSREF=EXP(LOG(-SS0(IREF))*X0+LOG(-SS0(IREF+1))*X1)
c     SREF=STREF+SSREF
c     ALM=2.997925E18/FREQ(IJ)
c     WRITE(36,636) IJ,ALM,IREF,DMREF,TREF,STREF,SSREF,SREF
c 636 FORMAT(1H ,I3,F10.3,I4,1PE10.3,0PF10.1,1X,1P3E10.3)
c     END IF
C
C       Contribution to J and H
C
c        do id=1,nud(iu)
c          rad1(id)=rad1(id)+wmuj(iu,id)*uf(id)
c          ali1(id)=ali1(id)+wmuj(iu,id)*af(id)
c        end do
c        FLUXc(IJ)=FLUXc(IJ)+WMUH(IU)*RIM(1)
C
C
C     end of the loop over frequencies
C
  500 CONTINUE
      RETURN
      END
C
C
C ***********************************************************************
C
C
      SUBROUTINE VELSET
C     =================
C
C     Determination of the macroscopic velocity as a function of depth
C
C     Input:  
C
C     RSTAR   - stellar radius (in solar radii or in cm)
C     RMAX    - maximum radial extent (in stellar radii)
C     AMLOSS  - mass loss rate ( in solar masses per year)
C     VELMAX  - maximum velocity (= V_infinity) - in km/s
C     BETA    - beta exponent in the beta-law for velocity
C     NDRAD   - Number of layers
C     NRCORE  - Number of core rays
C
C
c      parameter (un=1.,two=2.)
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'WINCOM.FOR'
      dimension zz(mdepth),vel0(mdepth),rrel(mdepth),
c     *          dvel0(mdepth),vel1(mdepth),hstt(mdepth),
     *          den0(mdepth),vel00(mdepth),ind(mdepth),
     *          densa(mdepth),eleca(mdepth),tempa(mdepth),
     *          rda(mdepth),rrela(mdepth),vel0a(mdepth)
c
      un=1
      two=2.
      read(55,*,err=100,end=100) rstar,rmax,amloss,vinf,beta,
     *           ndrad,nrcore,nfiry,ndf,nda
      rstr=rstar
      if(rstar.lt.1.e5) rstr=rstar*6.9598e10
      amdot=amloss*6.3029e25
      RCORE=RSTR
      XMDOT=amdot
      BETAV=beta
      con=amdot/12.566e5
      conr=con/rstr/rstr
      nrext0=ndrad-nd
      zz(nd+nrext0)=0.
      rd(nd+nrext0)=rstr
      rrel(nd+nrext0)=1.
      do iid=1,nd-1
         id=nd-iid
         zz(id+nrext0)=zz(id+1+nrext0)+2.*(dm(id+1)-dm(id))/
     *                 (dens(id+1)+dens(id))
         rd(id+nrext0)=rstr+zz(id+nrext0)
         rrel(id+nrext0)=rd(id+nrext0)/rstr
      end do
C
      do id=1+nrext0,nd+nrext0
         vel0(id)=con/rd(id)**2/dens(id-nrext0)
         vel00(id)=vel0(id)
         if(vel00(id).gt.vinf) vel00(id)=vinf
      end do
      vin=vel0(nrext0+1)
      r1=rrel(nrext0+1)
C
      if(rrel(1+nrext0).lt.rmax.and.nd.lt.ndrad) then
      rl1=1.-1./rrel(1+nrext0)
      rl2=1.-1./rmax
      drl=(rl2-rl1)/nrext0
      do id=1,nrext0
         rlo=rl2-(id-1)*drl
         rrel(id)=1./(1.-rlo)
         rd(id)=rrel(id)*rstr
      end do
      end if
c
      do id=nd+nrext0-1,nrext0+1,-1
         r0=rrel(id)
         numid=0
         do id1=nd+nrext0-1,nrext0+1,-1
            x=un-r0/rrel(id1)
            if(x.lt.1.e-6) x=1.e-6
            v2=vinf*x**beta 
            ind(id1)=0
            if(v2.ge.vel0(id1)) then
              ind(id1)=id1
              numid=numid+1
            end if
         end do
         if(numid.eq.0) go to 10
         rsum=0.
         isum=0
         do id1=nd+nrext0-1,nrext0+1,-1
            if(ind(id1).gt.0) then
               rsum=rsum+rrel(id1)
               isum=isum+id1
            endif
         end do
         rc=rsum/numid
         idc=isum/numid
         numid0=numid
         r00=r0
      end do
   10 continue
      v1=vel0(idc)
      r0=(r0+r00)*0.5
      if(r0.lt.rc) v2=vinf*(un-r0/rc)**beta 
      write(6,602) numid0,idc,rc,r0,v1,v2
  602 format('numid,idc,rc,r0,v1,v2 ',2i4,4f10.5)
c
      do id=nd+nrext0-1,1,-1
         if(rrel(id).gt.rc.and.rrel(id).gt.r0) 
     *      vel0(id)=vinf*(1.-r0/rrel(id))**beta
      end do
c
      t1=temp(1)
      erel=elec(1)/dens(1)
      do id=nd,1,-1
         temp(id+nrext0)=temp(id)
         den0(id+nrext0)=dens(id)
         elec(id+nrext0)=elec(id)
         do i=1,nlevel
            popul(i,id+nrext0)=popul(i,id)
         end do
         WMM(ID+nrext0)=WMM(id)
         WMY(ID+nrext0)=WMY(id)
         YTOT(ID+nrext0)=YTOT(id)
         do i=1,natom
            relab(i,id+nrext0)=relab(i,id)
            abund(i,id+nrext0)=abund(i,id)
         end do
         do i=1,matom
            abndd(i,id+nrext0)=abndd(i,id)
         end do
      end do
C
      do id=1,nrext0
         TEMP(ID)=T1
         WMM(ID)=WMM(NREXT0+1)
         WMY(ID)=WMY(NREXT0+1)
         YTOT(ID)=YTOT(NREXT0+1)
         do i=1,natom
            relab(i,id)=relab(i,nrext0+1)
            abund(i,id)=abund(i,nrext0+1)
         end do
         do i=1,matom
            abndd(i,id)=abndd(i,nrext0+1)
         end do
      end do
      idstd=idstd+nrext0
c
      VINF=vinf*1.e5
      write(6,600)
      do id=1,nd+nrext0
         if(vel0(id).gt.0.) dens(id)=con/rd(id)**2/vel0(id)
         VEL(ID)=vel0(id)*1.e5
c        velc(id)=vel0(id)/2.997925e5
      end do
c
      do id=nd,1,-1
         id1=id+nrext0
         elec(id1)=elec(id1)*dens(id1)/den0(id1)
         do i=1,nlevel
            popul(i,id1)=popul(i,id1)*dens(id1)/den0(id1)
         end do
      end do
c
      do id=1,nrext0
         elec(id)=elec(nrext0+1)*dens(id)/dens(nrext0+1)
         do i=1,nlevel
            popul(i,id)=popul(i,nrext0+1)*dens(id)/dens(nrext0+1)
         end do
      end do
C
      ND=NDRAD
      if(ndf.eq.0) ndf=nd
      do id=1,nd
         write(6,601) id,dm(id),temp(id),elec(id),dens(id),rd(id),
     *                rrel(id),vel0(id)
         write(96,601) id,dm(id),temp(id),elec(id),dens(id),rd(id),
     *                rrel(id),vel0(id),vel00(id)
      end do
  600 format('    ID    M      TEMP       ELEC      DENS             ',
     *       'R       Rrel     VEL'/)
  601 format(1h ,i3,1pe10.3,0pf8.0,1p3e12.3,0pf10.4,0p2f8.2)
C
C
      if(nda.gt.0) then
         XR1=LOG(DENS(1))
         XR2=LOG(DENS(ND))
         DXR=(XR2-XR1)/FLOAT(NDA-1)
         DO ID=1,NDA
            DENSA(ID)=EXP(XR1+FLOAT(ID-1)*DXR)
         END DO
         CALL INTERP(DENS,TEMP,DENSA,TEMPA,ND,NDA,3,1,1)
         CALL INTERP(DENS,ELEC,DENSA,ELECA,ND,NDA,3,1,1)
         CALL INTERP(DENS,RD,DENSA,RDA,ND,NDA,3,1,1)
         CALL INTERP(DENS,RREl,DENSA,RRELA,ND,NDA,3,1,1)
         CALL INTERP(DENS,VEL0,DENSA,VEL0A,ND,NDA,3,1,1)
         do id=1,nda
         write(6,603) id,tempa(id),eleca(id),densa(id),rda(id),
     *                rrela(id),vel0a(id)
         write(96,603) id,tempa(id),eleca(id),densa(id),rda(id),
     *                rrela(id),vel0a(id)
         end do
      end if
  603 format(1h ,i3,0pf8.0,1p3e12.3,0pf10.4,0p2f8.2)
C
  100 continue
      return
      end
C
C
C ***********************************************************************
C
C

      SUBROUTINE RADTEM
C     =================
C
C     determination of the radiation temperatures
C     after Schmutz (1991); inversion done by Newton-Raphson
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'WINCOM.FOR'
      common/velaux/velmax,iemoff,nltoff,itrad
      PARAMETER (CON=2.0706D-16, un=1.)
      parameter (nterad=3)
C
      DO ID=1,ND
         rx=RD(ND)/RD(ID)
c        WDIL(ID)=0.5*(1.-sqrt(1.-rx*rx))
         wdil(id)=un-sqrt(un-rx*rx)
      END DO
      DO ITRD=1,NTERAD
         if(itrad.eq.0) then
         do id=1,nd
            trad(itrd,id)=temp(id)
         end do
         else
         II=0
         JJ=0
         IF(ITRD.LE.NION) II=NFIRST(ITRD)
         IF(ITRD.LE.NION) JJ=NNEXT(ITRD)
         DO ID=1,ND
            TRAD(ITRD,ID)=TEMP(ID)
            IF(II.GT.0) THEN
c           IF(II.GT.100000) THEN
            AA=POPUL(JJ,ID)/POPUL(II,ID)*ELEC(ID)*CON
            AA=AA*G(II)/G(JJ)/WDIL(ID)/SQRT(TEMP(ID))
            TR=TEMP(ID)
            ITER=0
   10       ITER=ITER+1
            XX=ENION(II)/BOLK/TR
            DTR=(AA*EXP(XX)-TR)/(1.+XX)
            DTRR=DTR/TR
            TR=TR+DTR
            IF(ABS(DTRR).GT.1.E-3.AND.ITER.LT.100) GO TO 10
            TRAD(ITRD,ID)=TR
            END IF
         END DO
         end if
      END DO
      write(6,600)
      do id=1,nd
         write(6,601) id,temp(id),trad(1,id),trad(2,id),trad(3,id)
      end do
  600 format(/' radiation temperatures/')
  601 format(i5,4f10.1)
      RETURN
      END
C
C
C ***********************************************************************
C
C
      FUNCTION SBFCH(FR,T)
C     ====================
C     
C     cross-section times partition function for CH
C
C     from Kurucz ATLAS9
C
      INCLUDE 'PARAMS.FOR'
      parameter (fihu=500.,fihui=1./fihu,
     *           twhu=200.,twhui=1./twhu,
     *           tenl=2.30258509299405E0)
c
      DIMENSION CROSSCH(15,105),PARTCH(41),CROSSCHT(15)
      DIMENSION C1(150),C2(150),C3(150),C4(150),C5(150)
      DIMENSION C6(150),C7(150),C8(150),C9(150),C10(150)
      DIMENSION C11(75)
C
      EQUIVALENCE (CROSSCH(1, 1),C1(1)),(CROSSCH(1,11),C2(1))
      EQUIVALENCE (CROSSCH(1,21),C3(1)),(CROSSCH(1,31),C4(1))
      EQUIVALENCE (CROSSCH(1,41),C5(1)),(CROSSCH(1,51),C6(1))
      EQUIVALENCE (CROSSCH(1,61),C7(1)),(CROSSCH(1,71),C8(1))
      EQUIVALENCE (CROSSCH(1,81),C9(1)),(CROSSCH(1,91),C10(1))
      EQUIVALENCE (CROSSCH(1,101),C11(1))
C
      DATA C1/-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,
     1-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,-38.000,
     2        -32.727,-31.151,-30.133,-29.432,-28.925,-28.547,-28.257,
     2-28.030,-27.848,-27.701,-27.580,-27.479,-27.395,-27.322,-27.261,
     3        -31.588,-30.011,-28.993,-28.290,-27.784,-27.405,-27.115,
     3-26.887,-26.705,-26.558,-26.437,-26.336,-26.251,-26.179,-26.117,
     4        -30.407,-28.830,-27.811,-27.108,-26.601,-26.223,-25.932,
     4-25.705,-25.523,-25.376,-25.255,-25.154,-25.069,-24.997,-24.935,
     5        -29.513,-27.937,-26.920,-26.218,-25.712,-25.334,-25.043,
     5-24.816,-24.635,-24.487,-24.366,-24.266,-24.181,-24.109,-24.047,
     6        -28.910,-27.341,-26.327,-25.628,-25.123,-24.746,-24.457,
     6-24.230,-24.049,-23.902,-23.782,-23.681,-23.597,-23.525,-23.464,
     7        -28.517,-26.961,-25.955,-25.261,-24.760,-24.385,-24.098,
     7-23.873,-23.694,-23.548,-23.429,-23.329,-23.245,-23.174,-23.113,
     8        -28.213,-26.675,-25.680,-24.993,-24.497,-24.127,-23.843,
     8-23.620,-23.443,-23.299,-23.181,-23.082,-22.999,-22.929,-22.869,
     9        -27.942,-26.427,-25.446,-24.769,-24.280,-23.915,-23.635,
     9-23.416,-23.241,-23.100,-22.983,-22.887,-22.805,-22.736,-22.677,
     A        -27.706,-26.210,-25.241,-24.572,-24.088,-23.728,-23.451,
     A-23.235,-23.063,-22.923,-22.808,-22.713,-22.633,-22.565,-22.507/
      DATA C2/-27.475,-26.000,-25.043,-24.382,-23.905,-23.548,-23.275,
     1-23.062,-22.891,-22.753,-22.640,-22.546,-22.467,-22.400,-22.343,
     2        -27.221,-25.783,-24.844,-24.193,-23.723,-23.372,-23.102,
     2-22.892,-22.724,-22.588,-22.476,-22.384,-22.306,-22.240,-22.184,
     3        -26.863,-25.506,-24.607,-23.979,-23.523,-23.182,-22.919,
     3-22.714,-22.550,-22.417,-22.309,-22.218,-22.142,-22.078,-22.023,
     4        -26.685,-25.347,-24.457,-23.835,-23.382,-23.044,-22.784,
     4-22.580,-22.418,-22.286,-22.178,-22.089,-22.014,-21.950,-21.896,
     5        -26.085,-24.903,-24.105,-23.538,-23.120,-22.805,-22.561,
     5-22.370,-22.217,-22.093,-21.991,-21.906,-21.835,-21.775,-21.723,
     6        -25.902,-24.727,-23.936,-23.376,-22.964,-22.654,-22.415,
     6-22.227,-22.076,-21.955,-21.855,-21.772,-21.702,-21.644,-21.593,
     7        -25.215,-24.196,-23.510,-23.019,-22.655,-22.378,-22.163,
     7-21.992,-21.855,-21.744,-21.653,-21.577,-21.513,-21.459,-21.412,
     8        -24.914,-23.937,-23.284,-22.820,-22.475,-22.212,-22.007,
     8-21.845,-21.715,-21.609,-21.522,-21.449,-21.388,-21.336,-21.292,
     9        -24.519,-23.637,-23.039,-22.606,-22.281,-22.030,-21.834,
     9-21.678,-21.552,-21.450,-21.365,-21.295,-21.236,-21.185,-21.142,
     A        -24.086,-23.222,-22.650,-22.246,-21.948,-21.722,-21.546,
     A-21.407,-21.296,-21.205,-21.131,-21.070,-21.018,-20.974,-20.937/
      DATA C3/-23.850,-23.018,-22.472,-22.088,-21.805,-21.590,-21.422,
     1-21.289,-21.182,-21.095,-21.024,-20.964,-20.914,-20.872,-20.835,
     2        -23.136,-22.445,-21.994,-21.676,-21.440,-21.259,-21.117,
     2-21.004,-20.912,-20.837,-20.775,-20.723,-20.679,-20.642,-20.611,
     3        -23.199,-22.433,-21.927,-21.573,-21.314,-21.119,-20.969,
     3-20.851,-20.758,-20.682,-20.621,-20.571,-20.529,-20.493,-20.463,
     4        -22.696,-22.020,-21.585,-21.286,-21.071,-20.912,-20.791,
     4-20.697,-20.622,-20.563,-20.514,-20.475,-20.442,-20.414,-20.391,
     5        -22.119,-21.557,-21.194,-20.943,-20.761,-20.624,-20.518,
     5-20.434,-20.367,-20.313,-20.268,-20.231,-20.201,-20.175,-20.153,
     6        -21.855,-21.300,-20.931,-20.673,-20.485,-20.344,-20.235,
     6-20.151,-20.084,-20.031,-19.988,-19.953,-19.924,-19.900,-19.880,
     7        -21.126,-20.673,-20.382,-20.184,-20.044,-19.943,-19.868,
     7-19.811,-19.769,-19.736,-19.710,-19.690,-19.674,-19.662,-19.652,
     8        -20.502,-20.150,-19.922,-19.766,-19.657,-19.578,-19.520,
     8-19.478,-19.446,-19.422,-19.404,-19.390,-19.379,-19.371,-19.365,
     9        -20.030,-19.724,-19.530,-19.399,-19.309,-19.245,-19.199,
     9-19.166,-19.142,-19.125,-19.112,-19.103,-19.096,-19.091,-19.088,
     A        -19.640,-19.364,-19.189,-19.074,-18.996,-18.943,-18.906,
     A-18.881,-18.863,-18.852,-18.844,-18.839,-18.837,-18.836,-18.836/
      DATA C4/-19.333,-19.092,-18.939,-18.838,-18.770,-18.725,-18.695,
     1-18.675,-18.662,-18.655,-18.651,-18.649,-18.649,-18.651,-18.653,
     2        -19.070,-18.880,-18.756,-18.674,-18.621,-18.585,-18.562,
     2-18.548,-18.540,-18.536,-18.536,-18.537,-18.539,-18.542,-18.546,
     3        -18.851,-18.708,-18.617,-18.558,-18.521,-18.498,-18.484,
     3-18.477,-18.475,-18.476,-18.478,-18.482,-18.487,-18.493,-18.498,
     4        -18.709,-18.599,-18.533,-18.494,-18.471,-18.459,-18.454,
     4-18.454,-18.457,-18.462,-18.469,-18.476,-18.483,-18.490,-18.498,
     5        -18.656,-18.572,-18.524,-18.497,-18.485,-18.480,-18.482,
     5-18.486,-18.493,-18.501,-18.510,-18.519,-18.527,-18.536,-18.544,
     6        -18.670,-18.613,-18.582,-18.566,-18.561,-18.562,-18.568,
     6-18.575,-18.583,-18.592,-18.601,-18.610,-18.619,-18.627,-18.635,
     7        -18.728,-18.700,-18.687,-18.683,-18.685,-18.691,-18.698,
     7-18.706,-18.715,-18.723,-18.731,-18.739,-18.745,-18.752,-18.758,
     8        -18.839,-18.835,-18.836,-18.842,-18.849,-18.857,-18.865,
     8-18.872,-18.878,-18.883,-18.888,-18.892,-18.895,-18.898,-18.900,
     9        -19.034,-19.041,-19.049,-19.057,-19.064,-19.069,-19.071,
     9-19.071,-19.070,-19.068,-19.065,-19.061,-19.058,-19.054,-19.051,
     A        -19.372,-19.378,-19.382,-19.380,-19.372,-19.359,-19.341,
     A-19.321,-19.300,-19.280,-19.261,-19.243,-19.227,-19.212,-19.199/
      DATA C5/-19.780,-19.777,-19.763,-19.732,-19.686,-19.631,-19.573,
     1-19.517,-19.465,-19.419,-19.379,-19.344,-19.314,-19.288,-19.265,
     2        -20.151,-20.133,-20.087,-20.009,-19.911,-19.810,-19.715,
     2-19.631,-19.559,-19.497,-19.446,-19.402,-19.365,-19.333,-19.306,
     3        -20.525,-20.454,-20.312,-20.138,-19.970,-19.825,-19.705,
     3-19.607,-19.528,-19.464,-19.411,-19.367,-19.330,-19.300,-19.274,
     4        -20.869,-20.655,-20.366,-20.104,-19.894,-19.731,-19.604,
     4-19.505,-19.426,-19.363,-19.312,-19.271,-19.236,-19.208,-19.184,
     5        -21.179,-20.768,-20.380,-20.081,-19.856,-19.686,-19.556,
     5-19.454,-19.375,-19.311,-19.260,-19.218,-19.184,-19.155,-19.131,
     6        -21.167,-20.601,-20.206,-19.925,-19.719,-19.565,-19.447,
     6-19.355,-19.283,-19.226,-19.180,-19.143,-19.112,-19.087,-19.066,
     7        -20.918,-20.348,-19.976,-19.720,-19.536,-19.401,-19.299,
     7-19.220,-19.159,-19.112,-19.073,-19.043,-19.018,-18.998,-18.981,
     8        -20.753,-20.204,-19.847,-19.602,-19.427,-19.299,-19.203,
     8-19.129,-19.072,-19.028,-18.993,-18.965,-18.942,-18.924,-18.909,
     9        -20.456,-19.987,-19.677,-19.460,-19.302,-19.186,-19.098,
     9-19.030,-18.978,-18.937,-18.904,-18.878,-18.857,-18.841,-18.827,
     A        -20.154,-19.734,-19.461,-19.272,-19.136,-19.035,-18.960,
     A-18.902,-18.858,-18.824,-18.797,-18.775,-18.759,-18.745,-18.735/
      DATA C6/-19.941,-19.544,-19.288,-19.114,-18.992,-18.903,-18.837,
     1-18.788,-18.751,-18.723,-18.701,-18.684,-18.671,-18.661,-18.654,
     2        -19.657,-19.321,-19.104,-18.956,-18.853,-18.779,-18.724,
     2-18.684,-18.655,-18.632,-18.615,-18.602,-18.592,-18.585,-18.579,
     3        -19.388,-19.109,-18.930,-18.810,-18.725,-18.664,-18.620,
     3-18.586,-18.562,-18.543,-18.529,-18.518,-18.510,-18.503,-18.498,
     4        -19.201,-18.953,-18.794,-18.686,-18.611,-18.556,-18.515,
     4-18.485,-18.462,-18.446,-18.433,-18.423,-18.416,-18.410,-18.406,
     5        -18.923,-18.719,-18.588,-18.500,-18.439,-18.396,-18.365,
     5-18.344,-18.328,-18.318,-18.311,-18.307,-18.304,-18.303,-18.302,
     6        -18.614,-18.458,-18.361,-18.298,-18.258,-18.232,-18.216,
     6-18.206,-18.202,-18.201,-18.202,-18.205,-18.208,-18.213,-18.218,
     7        -18.419,-18.295,-18.222,-18.178,-18.153,-18.139,-18.132,
     7-18.131,-18.133,-18.138,-18.143,-18.150,-18.157,-18.164,-18.172,
     8        -18.296,-18.201,-18.148,-18.118,-18.101,-18.094,-18.091,
     8-18.093,-18.096,-18.101,-18.107,-18.113,-18.120,-18.126,-18.132,
     9        -18.021,-17.992,-17.977,-17.970,-17.967,-17.968,-17.970,
     9-17.974,-17.978,-17.983,-17.989,-17.994,-18.000,-18.005,-18.011,
     A        -17.694,-17.686,-17.686,-17.691,-17.698,-17.708,-17.718,
     A-17.729,-17.740,-17.750,-17.761,-17.771,-17.781,-17.790,-17.798/
      DATA C7/-17.374,-17.384,-17.400,-17.420,-17.440,-17.462,-17.483,
     1-17.503,-17.523,-17.541,-17.558,-17.575,-17.590,-17.603,-17.616,
     2        -17.169,-17.199,-17.230,-17.262,-17.293,-17.323,-17.351,
     2-17.378,-17.404,-17.427,-17.449,-17.469,-17.488,-17.505,-17.520,
     3        -17.151,-17.184,-17.217,-17.250,-17.282,-17.313,-17.342,
     3-17.369,-17.395,-17.418,-17.440,-17.461,-17.480,-17.497,-17.513,
     4        -17.230,-17.260,-17.290,-17.320,-17.348,-17.375,-17.401,
     4-17.425,-17.448,-17.469,-17.489,-17.508,-17.525,-17.541,-17.556,
     5        -17.379,-17.403,-17.425,-17.446,-17.467,-17.486,-17.505,
     5-17.524,-17.541,-17.558,-17.574,-17.588,-17.602,-17.615,-17.627,
     6        -17.596,-17.604,-17.609,-17.612,-17.616,-17.622,-17.628,
     6-17.636,-17.644,-17.652,-17.661,-17.670,-17.679,-17.687,-17.695,
     7        -17.846,-17.823,-17.795,-17.770,-17.750,-17.735,-17.725,
     7-17.719,-17.716,-17.715,-17.716,-17.719,-17.722,-17.726,-17.730,
     8        -18.089,-18.015,-17.942,-17.882,-17.836,-17.802,-17.777,
     8-17.760,-17.748,-17.740,-17.736,-17.734,-17.733,-17.734,-17.736,
     9        -18.299,-18.156,-18.038,-17.947,-17.881,-17.833,-17.798,
     9-17.774,-17.757,-17.745,-17.738,-17.733,-17.730,-17.729,-17.729,
     A        -18.441,-18.243,-18.096,-17.991,-17.915,-17.860,-17.821,
     A-17.792,-17.772,-17.757,-17.746,-17.738,-17.733,-17.730,-17.728/
      DATA C8/-18.474,-18.262,-18.111,-18.004,-17.926,-17.869,-17.826,
     1-17.795,-17.771,-17.753,-17.740,-17.730,-17.722,-17.717,-17.713,
     2        -18.387,-18.191,-18.053,-17.952,-17.878,-17.823,-17.782,
     2-17.752,-17.729,-17.711,-17.698,-17.689,-17.681,-17.676,-17.672,
     3        -18.161,-17.990,-17.874,-17.793,-17.736,-17.696,-17.668,
     3-17.648,-17.634,-17.625,-17.619,-17.616,-17.614,-17.614,-17.615,
     4        -17.908,-17.774,-17.690,-17.637,-17.604,-17.583,-17.572,
     4-17.567,-17.566,-17.568,-17.571,-17.576,-17.581,-17.587,-17.593,
     5        -17.681,-17.589,-17.540,-17.515,-17.506,-17.505,-17.511,
     5-17.520,-17.530,-17.542,-17.554,-17.566,-17.578,-17.589,-17.600,
     6        -17.647,-17.606,-17.584,-17.575,-17.573,-17.576,-17.582,
     6-17.589,-17.597,-17.605,-17.614,-17.623,-17.631,-17.639,-17.646,
     7        -17.300,-17.291,-17.291,-17.297,-17.307,-17.319,-17.333,
     7-17.347,-17.361,-17.375,-17.389,-17.402,-17.415,-17.427,-17.438,
     8        -16.786,-16.802,-16.825,-16.853,-16.883,-16.914,-16.944,
     8-16.974,-17.003,-17.030,-17.055,-17.079,-17.101,-17.122,-17.141,
     9        -16.489,-16.533,-16.579,-16.625,-16.670,-16.713,-16.754,
     9-16.793,-16.830,-16.864,-16.896,-16.925,-16.952,-16.977,-17.000,
     A        -16.694,-16.724,-16.756,-16.789,-16.823,-16.856,-16.888,
     A-16.919,-16.949,-16.976,-17.002,-17.026,-17.048,-17.069,-17.088/
      DATA C9/-16.935,-16.951,-16.971,-16.993,-17.016,-17.040,-17.064,
     1-17.088,-17.111,-17.132,-17.153,-17.172,-17.190,-17.206,-17.222,
     2        -17.200,-17.208,-17.220,-17.235,-17.251,-17.269,-17.286,
     2-17.304,-17.322,-17.338,-17.354,-17.369,-17.384,-17.397,-17.409,
     3        -17.597,-17.591,-17.589,-17.590,-17.594,-17.600,-17.608,
     3-17.617,-17.626,-17.635,-17.645,-17.654,-17.662,-17.671,-17.679,
     4        -18.166,-18.134,-18.107,-18.085,-18.068,-18.056,-18.047,
     4-18.041,-18.038,-18.036,-18.035,-18.035,-18.036,-18.038,-18.039,
     5        -19.000,-18.917,-18.838,-18.770,-18.714,-18.669,-18.632,
     5-18.603,-18.579,-18.560,-18.545,-18.532,-18.522,-18.514,-18.507,
     6        -20.313,-19.982,-19.754,-19.592,-19.472,-19.380,-19.309,
     6-19.253,-19.208,-19.172,-19.143,-19.119,-19.099,-19.083,-19.069,
     7        -19.751,-19.611,-19.520,-19.461,-19.423,-19.398,-19.382,
     7-19.372,-19.366,-19.364,-19.363,-19.364,-19.366,-19.368,-19.371,
     8        -19.581,-19.431,-19.337,-19.277,-19.240,-19.218,-19.207,
     8-19.202,-19.203,-19.207,-19.212,-19.220,-19.228,-19.236,-19.245,
     9        -19.685,-19.506,-19.389,-19.311,-19.258,-19.222,-19.199,
     9-19.184,-19.175,-19.170,-19.168,-19.169,-19.171,-19.174,-19.177,
     A        -19.977,-19.756,-19.606,-19.501,-19.425,-19.370,-19.330,
     A-19.300,-19.278,-19.262,-19.250,-19.241,-19.235,-19.230,-19.227/
      DATAC10/-20.445,-20.158,-19.958,-19.815,-19.711,-19.633,-19.574,
     1-19.528,-19.493,-19.465,-19.442,-19.425,-19.410,-19.398,-19.389,
     2        -20.980,-20.625,-20.391,-20.229,-20.110,-20.020,-19.949,
     2-19.892,-19.846,-19.807,-19.775,-19.748,-19.724,-19.704,-19.687,
     3        -21.404,-21.023,-20.771,-20.594,-20.461,-20.358,-20.274,
     3-20.205,-20.148,-20.099,-20.058,-20.022,-19.991,-19.965,-19.942,
     4        -21.309,-20.970,-20.753,-20.603,-20.495,-20.412,-20.348,
     4-20.295,-20.252,-20.215,-20.185,-20.158,-20.135,-20.115,-20.098,
     5        -21.221,-20.906,-20.707,-20.574,-20.480,-20.412,-20.361,
     5-20.322,-20.292,-20.268,-20.249,-20.233,-20.221,-20.210,-20.201,
     6        -21.441,-21.097,-20.878,-20.728,-20.623,-20.546,-20.489,
     6-20.446,-20.413,-20.387,-20.368,-20.352,-20.340,-20.330,-20.322,
     7        -21.668,-21.305,-21.071,-20.911,-20.797,-20.713,-20.650,
     7-20.602,-20.565,-20.536,-20.514,-20.496,-20.481,-20.470,-20.460,
     8        -21.926,-21.556,-21.316,-21.150,-21.031,-20.942,-20.874,
     8-20.822,-20.782,-20.750,-20.724,-20.704,-20.687,-20.674,-20.663,
     9        -22.319,-21.937,-21.686,-21.510,-21.380,-21.282,-21.206,
     9-21.147,-21.099,-21.061,-21.031,-21.006,-20.985,-20.968,-20.954,
     A        -22.969,-22.561,-22.288,-22.092,-21.945,-21.832,-21.743,
     A-21.672,-21.616,-21.570,-21.533,-21.503,-21.477,-21.457,-21.439/
      DATAC11/-24.001,-23.527,-23.199,-22.957,-22.772,-22.629,-22.516,
     1-22.427,-22.355,-22.297,-22.250,-22.212,-22.180,-22.153,-22.131,
     2        -24.233,-23.774,-23.477,-23.273,-23.128,-23.022,-22.943,
     2-22.883,-22.837,-22.802,-22.774,-22.752,-22.735,-22.721,-22.710,
     3        -24.550,-23.913,-23.521,-23.266,-23.094,-22.976,-22.893,
     3-22.836,-22.796,-22.768,-22.750,-22.737,-22.730,-22.726,-22.725,
     4        -24.301,-23.665,-23.274,-23.019,-22.848,-22.730,-22.648,
     4-22.591,-22.552,-22.525,-22.507,-22.495,-22.489,-22.485,-22.485,
     5        -24.519,-23.883,-23.491,-23.237,-23.065,-22.948,-22.866,
     5-22.809,-22.770,-22.743,-22.724,-22.713,-22.706,-22.703,-22.702/
      DATA PARTCH/
     1     203.741,  249.643,  299.341,  353.477,  412.607,  477.237,
     2     547.817,  624.786,  708.543,  799.463,  897.912, 1004.227,
     3    1118.738, 1241.761, 1373.588, 1514.481, 1664.677, 1824.394,
     4    1993.801, 2173.050, 2362.234, 2561.424, 2770.674, 2989.930,
     5    3219.204, 3458.378, 3707.355, 3966.005, 4234.155, 4511.604,
     6    4798.135, 5093.554, 5397.593, 5709.948, 6030.401, 6358.646,
     7    6694.379, 7037.313, 7387.147, 7743.579, 8106.313/
      DATA FREQ1/0./
C
      SBFCH=0.
      IF(FR.EQ.FREQ1) GO TO 30
      FREQ1=FR
      WAVENO=FR/2.99792458E10
      EVOLT=WAVENO/8065.479
      N=int(EVOLT*10.)
      EN=FLOAT(N)*.1
      IF(N.LT.20) RETURN
      IF(N.GE.105) RETURN
c
      DO IT=1,15
         CROSSCHT(IT)=(CROSSCH(IT,N)+(CROSSCH(IT,N+1)-CROSSCH(IT,N))*
     *                (EVOLT-EN)*10.)
      END DO
c
c     interpolate to obtain partition function
c
   30 IF(T.GE.9000.) RETURN
      IF(N.LT.20) RETURN
      IF(N.GE.105) RETURN
      IT=int((T-1000.)*twhui+1.)
      IF(IT.LT.1) IT=1
      TN=FLOAT(IT)*twhu+800.
      PART=PARTCH(IT)+(PARTCH(IT+1)-PARTCH(IT))*(T-TN)*twhui
c
c     interpolate to obtain cross-section
c
      IT=int((T-2000.)*fihui+1.)
      IF(IT.LT.1) IT=1
      TN=FLOAT(IT)*fihu+1500.
      SBFCH=EXP((CROSSCHT(IT)+(CROSSCHT(IT+1)-CROSSCHT(IT))*
c    *     (T-TN)*fihui)*tenl)*PART
     *     (T-TN)*fihui)*tenl)
      RETURN
      END
C
C
C ***********************************************************************
C
C

      FUNCTION SBFOH(FR,T)
C     ====================
C
C     cross-section times partition function for OH
C
C     from Kurucz ATLAS9
C
      INCLUDE 'PARAMS.FOR'
      parameter (fihu=500.,fihui=1./fihu,
     *           twhu=200.,twhui=1./twhu,
     *           tenl=2.30258509299405E0)
      DIMENSION CROSSOH(15,130),PARTOH(41),CROSSOHT(15)
      DIMENSION C1(150),C2(150),C3(150),C4(150),C5(150)
      DIMENSION C6(150),C7(150),C8(150),C9(150),C10(150)
      DIMENSION C11(150),C12(150),C13(150)
      EQUIVALENCE (CROSSOH(1, 1),C1(1)),(CROSSOH(1,11),C2(1))
      EQUIVALENCE (CROSSOH(1,21),C3(1)),(CROSSOH(1,31),C4(1))
      EQUIVALENCE (CROSSOH(1,41),C5(1)),(CROSSOH(1,51),C6(1))
      EQUIVALENCE (CROSSOH(1,61),C7(1)),(CROSSOH(1,71),C8(1))
      EQUIVALENCE (CROSSOH(1,81),C9(1)),(CROSSOH(1,91),C10(1))
      EQUIVALENCE (CROSSOH(1,101),C11(1))
      EQUIVALENCE (CROSSOH(1,111),C12(1))
      EQUIVALENCE (CROSSOH(1,121),C13(1))
C
      DATA C1/-30.855,-29.121,-27.976,-27.166,-26.566,-26.106,-25.742,
     1-25.448,-25.207,-25.006,-24.836,-24.691,-24.566,-24.457,-24.363,
     2        -30.494,-28.760,-27.615,-26.806,-26.206,-25.745,-25.381,
     2-25.088,-24.846,-24.645,-24.475,-24.330,-24.205,-24.097,-24.002,
     3        -30.157,-28.425,-27.280,-26.472,-25.872,-25.411,-25.048,
     3-24.754,-24.513,-24.312,-24.142,-23.997,-23.872,-23.764,-23.669,
     4        -29.848,-28.117,-26.974,-26.165,-25.566,-25.105,-24.742,
     4-24.448,-24.207,-24.006,-23.836,-23.692,-23.567,-23.458,-23.364,
     5        -29.567,-27.837,-26.693,-25.885,-25.286,-24.826,-24.462,
     5-24.169,-23.928,-23.727,-23.557,-23.412,-23.287,-23.179,-23.084,
     6        -29.307,-27.578,-26.436,-25.628,-25.029,-24.569,-24.205,
     6-23.912,-23.671,-23.470,-23.300,-23.155,-23.031,-22.922,-22.828,
     7        -29.068,-27.341,-26.199,-25.391,-24.792,-24.332,-23.969,
     7-23.676,-23.435,-23.234,-23.064,-22.920,-22.795,-22.687,-22.592,
     8        -28.820,-27.115,-25.978,-25.172,-24.574,-24.115,-23.752,
     8-23.459,-23.218,-23.017,-22.848,-22.703,-22.579,-22.470,-22.376,
     9        -28.540,-26.891,-25.768,-24.968,-24.372,-23.914,-23.552,
     9-23.259,-23.019,-22.818,-22.649,-22.504,-22.380,-22.272,-22.177,
     A        -28.275,-26.681,-25.574,-24.779,-24.186,-23.729,-23.368,
     A-23.076,-22.836,-22.636,-22.467,-22.322,-22.198,-22.090,-21.996/
      DATA C2/-27.993,-26.470,-25.388,-24.602,-24.014,-23.560,-23.200,
     1-22.909,-22.669,-22.470,-22.301,-22.157,-22.033,-21.925,-21.831,
     2        -27.698,-26.252,-25.204,-24.433,-23.851,-23.401,-23.043,
     2-22.754,-22.515,-22.316,-22.148,-22.005,-21.881,-21.773,-21.679,
     3        -27.398,-26.026,-25.019,-24.267,-23.696,-23.251,-22.896,
     3-22.609,-22.372,-22.174,-22.007,-21.864,-21.741,-21.634,-21.540,
     4        -27.100,-25.791,-24.828,-24.102,-23.543,-23.106,-22.756,
     4-22.472,-22.238,-22.041,-21.875,-21.733,-21.611,-21.504,-21.411,
     5        -26.807,-25.549,-24.631,-23.933,-23.391,-22.964,-22.621,
     5-22.341,-22.109,-21.915,-21.751,-21.610,-21.488,-21.383,-21.290,
     6        -26.531,-25.310,-24.431,-23.761,-23.238,-22.823,-22.488,
     6-22.214,-21.986,-21.795,-21.633,-21.494,-21.374,-21.269,-21.178,
     7        -26.239,-25.066,-24.225,-23.585,-23.082,-22.681,-22.356,
     7-22.089,-21.866,-21.679,-21.520,-21.383,-21.265,-21.162,-21.072,
     8        -25.945,-24.824,-24.017,-23.405,-22.923,-22.538,-22.223,
     8-21.964,-21.748,-21.565,-21.410,-21.276,-21.160,-21.059,-20.970,
     9        -25.663,-24.587,-23.810,-23.222,-22.761,-22.391,-22.088,
     9-21.838,-21.629,-21.452,-21.300,-21.170,-21.057,-20.958,-20.872,
     A        -25.372,-24.350,-23.603,-23.038,-22.596,-22.241,-21.950,
     A-21.710,-21.508,-21.337,-21.190,-21.064,-20.954,-20.858,-20.774/
      DATA C3/-25.076,-24.111,-23.396,-22.853,-22.429,-22.088,-21.809,
     1-21.578,-21.384,-21.220,-21.078,-20.957,-20.851,-20.758,-20.676,
     2        -24.779,-23.870,-23.189,-22.669,-22.261,-21.934,-21.667,
     2-21.445,-21.259,-21.101,-20.965,-20.848,-20.746,-20.656,-20.578,
     3        -24.486,-23.629,-22.983,-22.486,-22.095,-21.781,-21.524,
     3-21.311,-21.132,-20.980,-20.850,-20.737,-20.639,-20.553,-20.478,
     4        -24.183,-23.382,-22.774,-22.302,-21.928,-21.627,-21.381,
     4-21.177,-21.005,-20.859,-20.734,-20.625,-20.531,-20.449,-20.376,
     5        -23.867,-23.127,-22.561,-22.116,-21.761,-21.474,-21.238,
     5-21.043,-20.878,-20.738,-20.617,-20.513,-20.423,-20.344,-20.274,
     6        -23.538,-22.862,-22.340,-21.926,-21.592,-21.320,-21.096,
     6-20.909,-20.751,-20.617,-20.502,-20.402,-20.315,-20.239,-20.172,
     7        -23.234,-22.604,-22.120,-21.734,-21.422,-21.166,-20.953,
     7-20.776,-20.625,-20.497,-20.387,-20.291,-20.208,-20.135,-20.071,
     8        -22.934,-22.347,-21.898,-21.541,-21.250,-21.010,-20.811,
     8-20.643,-20.500,-20.378,-20.273,-20.182,-20.102,-20.033,-19.971,
     9        -22.637,-22.092,-21.676,-21.345,-21.075,-20.853,-20.666,
     9-20.508,-20.374,-20.259,-20.159,-20.073,-19.997,-19.931,-19.872,
     A        -22.337,-21.835,-21.452,-21.147,-20.899,-20.693,-20.520,
     A-20.373,-20.247,-20.139,-20.046,-19.964,-19.892,-19.830,-19.774/
      DATA C4/-22.049,-21.584,-21.230,-20.950,-20.721,-20.531,-20.372,
     1-20.236,-20.119,-20.019,-19.931,-19.855,-19.788,-19.729,-19.676,
     2        -21.768,-21.337,-21.011,-20.754,-20.544,-20.370,-20.223,
     2-20.098,-19.991,-19.898,-19.817,-19.746,-19.683,-19.628,-19.579,
     3        -21.494,-21.096,-20.796,-20.559,-20.367,-20.208,-20.074,
     3-19.960,-19.861,-19.776,-19.701,-19.636,-19.578,-19.527,-19.482,
     4        -21.233,-20.861,-20.585,-20.368,-20.193,-20.048,-19.926,
     4-19.821,-19.732,-19.654,-19.586,-19.526,-19.473,-19.426,-19.384,
     5        -20.983,-20.635,-20.380,-20.181,-20.021,-19.889,-19.778,
     5-19.683,-19.602,-19.531,-19.469,-19.415,-19.367,-19.324,-19.286,
     6        -20.743,-20.418,-20.182,-19.999,-19.853,-19.733,-19.633,
     6-19.547,-19.474,-19.410,-19.354,-19.305,-19.261,-19.223,-19.189,
     7        -20.515,-20.210,-19.991,-19.824,-19.690,-19.581,-19.490,
     7-19.413,-19.347,-19.290,-19.240,-19.196,-19.157,-19.122,-19.092,
     8        -20.297,-20.011,-19.808,-19.654,-19.532,-19.434,-19.352,
     8-19.282,-19.223,-19.172,-19.127,-19.088,-19.054,-19.023,-18.996,
     9        -20.090,-19.822,-19.633,-19.491,-19.381,-19.291,-19.218,
     9-19.156,-19.103,-19.057,-19.018,-18.983,-18.952,-18.925,-18.901,
     A        -19.893,-19.642,-19.467,-19.337,-19.236,-19.155,-19.089,
     A-19.034,-18.987,-18.946,-18.912,-18.881,-18.854,-18.831,-18.810/
      DATA C5/-19.705,-19.472,-19.309,-19.190,-19.098,-19.025,-18.966,
     1-18.917,-18.876,-18.840,-18.810,-18.783,-18.760,-18.739,-18.721,
     2        -19.527,-19.310,-19.161,-19.051,-18.968,-18.903,-18.851,
     2-18.807,-18.771,-18.740,-18.713,-18.690,-18.670,-18.653,-18.637,
     3        -19.357,-19.159,-19.022,-18.922,-18.847,-18.789,-18.743,
     3-18.704,-18.673,-18.646,-18.623,-18.603,-18.586,-18.571,-18.558,
     4        -19.195,-19.016,-18.892,-18.803,-18.736,-18.684,-18.643,
     4-18.610,-18.583,-18.560,-18.540,-18.523,-18.509,-18.496,-18.485,
     5        -19.042,-18.883,-18.772,-18.693,-18.634,-18.589,-18.553,
     5-18.525,-18.501,-18.481,-18.465,-18.451,-18.438,-18.428,-18.419,
     6        -18.894,-18.758,-18.662,-18.593,-18.542,-18.503,-18.473,
     6-18.448,-18.428,-18.412,-18.398,-18.386,-18.376,-18.367,-18.359,
     7        -18.752,-18.639,-18.559,-18.501,-18.458,-18.426,-18.400,
     7-18.380,-18.363,-18.350,-18.338,-18.328,-18.320,-18.313,-18.306,
     8        -18.611,-18.523,-18.460,-18.415,-18.381,-18.355,-18.334,
     8-18.318,-18.304,-18.293,-18.284,-18.276,-18.269,-18.263,-18.258,
     9        -18.471,-18.408,-18.362,-18.329,-18.304,-18.285,-18.269,
     9-18.257,-18.247,-18.238,-18.231,-18.224,-18.219,-18.214,-18.210,
     A        -18.330,-18.290,-18.261,-18.239,-18.223,-18.211,-18.201,
     A-18.192,-18.185,-18.179,-18.174,-18.169,-18.165,-18.162,-18.159/
      DATA C6/-18.190,-18.168,-18.154,-18.143,-18.135,-18.129,-18.124,
     1-18.120,-18.116,-18.112,-18.109,-18.106,-18.104,-18.102,-18.100,
     2        -18.055,-18.047,-18.043,-18.042,-18.040,-18.039,-18.039,
     2-18.038,-18.037,-18.036,-18.035,-18.034,-18.033,-18.033,-18.032,
     3        -17.929,-17.931,-17.935,-17.939,-17.943,-17.946,-17.948,
     3-17.950,-17.952,-17.953,-17.955,-17.956,-17.957,-17.958,-17.959,
     4        -17.818,-17.826,-17.834,-17.842,-17.849,-17.855,-17.860,
     4-17.865,-17.869,-17.872,-17.875,-17.878,-17.881,-17.883,-17.886,
     5        -17.724,-17.736,-17.747,-17.758,-17.767,-17.775,-17.782,
     5-17.788,-17.793,-17.798,-17.803,-17.807,-17.811,-17.815,-17.819,
     6        -17.651,-17.665,-17.678,-17.690,-17.701,-17.710,-17.718,
     6-17.725,-17.732,-17.738,-17.744,-17.749,-17.755,-17.760,-17.765,
     7        -17.601,-17.615,-17.629,-17.642,-17.653,-17.663,-17.672,
     7-17.680,-17.688,-17.695,-17.701,-17.708,-17.714,-17.720,-17.726,
     8        -17.572,-17.587,-17.602,-17.614,-17.626,-17.636,-17.645,
     8-17.654,-17.662,-17.670,-17.677,-17.684,-17.691,-17.698,-17.704,
     9        -17.565,-17.581,-17.595,-17.607,-17.619,-17.629,-17.638,
     9-17.647,-17.656,-17.664,-17.671,-17.679,-17.686,-17.693,-17.700,
     A        -17.580,-17.594,-17.608,-17.620,-17.630,-17.640,-17.650,
     A-17.658,-17.667,-17.675,-17.682,-17.690,-17.697,-17.704,-17.711/
      DATA C7/-17.613,-17.626,-17.639,-17.649,-17.659,-17.669,-17.677,
     1-17.686,-17.694,-17.701,-17.709,-17.716,-17.723,-17.730,-17.737,
     2        -17.663,-17.675,-17.685,-17.695,-17.703,-17.711,-17.719,
     2-17.727,-17.734,-17.741,-17.748,-17.755,-17.761,-17.768,-17.774,
     3        -17.728,-17.737,-17.745,-17.752,-17.759,-17.766,-17.772,
     3-17.778,-17.785,-17.791,-17.797,-17.803,-17.808,-17.814,-17.820,
     4        -17.803,-17.809,-17.814,-17.818,-17.823,-17.828,-17.832,
     4-17.837,-17.842,-17.847,-17.852,-17.856,-17.861,-17.866,-17.871,
     5        -17.884,-17.886,-17.888,-17.889,-17.891,-17.893,-17.896,
     5-17.899,-17.902,-17.905,-17.908,-17.912,-17.915,-17.919,-17.922,
     6        -17.966,-17.964,-17.961,-17.959,-17.958,-17.958,-17.958,
     6-17.959,-17.960,-17.961,-17.963,-17.964,-17.966,-17.968,-17.970,
     7        -18.040,-18.034,-18.028,-18.023,-18.019,-18.016,-18.013,
     7-18.012,-18.010,-18.010,-18.009,-18.009,-18.009,-18.009,-18.010,
     8        -18.096,-18.087,-18.078,-18.071,-18.065,-18.059,-18.055,
     8-18.051,-18.047,-18.045,-18.042,-18.040,-18.039,-18.037,-18.036,
     9        -18.125,-18.115,-18.105,-18.097,-18.089,-18.082,-18.076,
     9-18.070,-18.065,-18.061,-18.057,-18.053,-18.051,-18.048,-18.046,
     A        -18.120,-18.112,-18.103,-18.095,-18.087,-18.079,-18.072,
     A-18.066,-18.060,-18.055,-18.050,-18.046,-18.042,-18.039,-18.036/
      DATA C8/-18.083,-18.078,-18.071,-18.064,-18.057,-18.050,-18.044,
     1-18.037,-18.032,-18.026,-18.022,-18.017,-18.014,-18.010,-18.007,
     2        -18.025,-18.022,-18.017,-18.012,-18.006,-18.000,-17.994,
     2-17.989,-17.984,-17.979,-17.975,-17.971,-17.968,-17.965,-17.963,
     3        -17.957,-17.955,-17.952,-17.948,-17.943,-17.938,-17.934,
     3-17.929,-17.925,-17.922,-17.918,-17.916,-17.913,-17.911,-17.910,
     4        -17.890,-17.889,-17.886,-17.882,-17.879,-17.875,-17.871,
     4-17.867,-17.864,-17.862,-17.860,-17.858,-17.857,-17.856,-17.855,
     5        -17.831,-17.829,-17.826,-17.822,-17.819,-17.815,-17.812,
     5-17.810,-17.807,-17.806,-17.804,-17.803,-17.803,-17.803,-17.803,
     6        -17.786,-17.782,-17.777,-17.773,-17.769,-17.766,-17.763,
     6-17.761,-17.759,-17.758,-17.757,-17.757,-17.757,-17.758,-17.759,
     7        -17.753,-17.747,-17.741,-17.735,-17.731,-17.727,-17.724,
     7-17.722,-17.721,-17.720,-17.720,-17.720,-17.721,-17.722,-17.724,
     8        -17.733,-17.724,-17.716,-17.709,-17.703,-17.699,-17.696,
     8-17.694,-17.693,-17.692,-17.692,-17.693,-17.694,-17.695,-17.697,
     9        -17.723,-17.711,-17.700,-17.691,-17.685,-17.680,-17.676,
     9-17.674,-17.673,-17.672,-17.673,-17.673,-17.675,-17.676,-17.678,
     A        -17.718,-17.702,-17.689,-17.679,-17.672,-17.667,-17.663,
     A-17.660,-17.659,-17.659,-17.659,-17.660,-17.661,-17.663,-17.665/
      DATA C9/-17.713,-17.695,-17.681,-17.670,-17.662,-17.656,-17.653,
     1-17.650,-17.649,-17.649,-17.649,-17.650,-17.651,-17.653,-17.655,
     2        -17.705,-17.686,-17.671,-17.660,-17.652,-17.647,-17.643,
     2-17.641,-17.640,-17.640,-17.640,-17.641,-17.643,-17.645,-17.647,
     3        -17.690,-17.671,-17.657,-17.647,-17.640,-17.635,-17.632,
     3-17.630,-17.630,-17.630,-17.631,-17.632,-17.634,-17.636,-17.639,
     4        -17.667,-17.649,-17.637,-17.629,-17.623,-17.619,-17.618,
     4-17.617,-17.617,-17.618,-17.619,-17.621,-17.623,-17.626,-17.628,
     5        -17.635,-17.621,-17.611,-17.605,-17.601,-17.600,-17.599,
     5-17.599,-17.601,-17.602,-17.604,-17.607,-17.609,-17.612,-17.615,
     6        -17.596,-17.585,-17.579,-17.576,-17.575,-17.575,-17.576,
     6-17.578,-17.580,-17.582,-17.585,-17.588,-17.591,-17.595,-17.598,
     7        -17.550,-17.544,-17.542,-17.542,-17.544,-17.546,-17.548,
     7-17.552,-17.555,-17.558,-17.562,-17.566,-17.570,-17.573,-17.577,
     8        -17.501,-17.500,-17.501,-17.504,-17.508,-17.513,-17.517,
     8-17.521,-17.526,-17.530,-17.535,-17.539,-17.544,-17.548,-17.553,
     9        -17.449,-17.452,-17.457,-17.463,-17.470,-17.476,-17.482,
     9-17.488,-17.493,-17.499,-17.504,-17.509,-17.514,-17.519,-17.524,
     A        -17.396,-17.403,-17.412,-17.420,-17.429,-17.437,-17.444,
     A-17.451,-17.458,-17.464,-17.470,-17.476,-17.481,-17.487,-17.492/
      DATAC10/-17.344,-17.355,-17.366,-17.377,-17.387,-17.396,-17.405,
     1-17.413,-17.420,-17.427,-17.434,-17.440,-17.446,-17.452,-17.458,
     2        -17.295,-17.307,-17.321,-17.333,-17.345,-17.355,-17.365,
     2-17.373,-17.382,-17.389,-17.397,-17.404,-17.410,-17.417,-17.423,
     3        -17.249,-17.264,-17.278,-17.292,-17.304,-17.316,-17.326,
     3-17.335,-17.344,-17.352,-17.360,-17.368,-17.375,-17.382,-17.389,
     4        -17.209,-17.225,-17.241,-17.255,-17.268,-17.280,-17.291,
     4-17.301,-17.310,-17.319,-17.327,-17.335,-17.343,-17.350,-17.357,
     5        -17.177,-17.194,-17.210,-17.225,-17.239,-17.251,-17.262,
     5-17.272,-17.282,-17.291,-17.300,-17.308,-17.316,-17.324,-17.331,
     6        -17.154,-17.172,-17.189,-17.204,-17.218,-17.230,-17.242,
     6-17.252,-17.262,-17.272,-17.280,-17.289,-17.298,-17.306,-17.314,
     7        -17.144,-17.162,-17.179,-17.194,-17.208,-17.220,-17.232,
     7-17.242,-17.253,-17.262,-17.271,-17.280,-17.289,-17.297,-17.306,
     8        -17.146,-17.164,-17.181,-17.196,-17.210,-17.222,-17.234,
     8-17.245,-17.255,-17.265,-17.274,-17.283,-17.292,-17.301,-17.309,
     9        -17.163,-17.180,-17.197,-17.212,-17.225,-17.237,-17.249,
     9-17.260,-17.270,-17.280,-17.289,-17.298,-17.307,-17.316,-17.325,
     A        -17.193,-17.211,-17.227,-17.241,-17.254,-17.266,-17.277,
     A-17.288,-17.298,-17.308,-17.317,-17.327,-17.336,-17.345,-17.353/
      DATAC11/-17.239,-17.256,-17.271,-17.284,-17.297,-17.309,-17.320,
     1-17.330,-17.340,-17.350,-17.359,-17.369,-17.378,-17.387,-17.395,
     2        -17.299,-17.315,-17.329,-17.342,-17.354,-17.365,-17.376,
     2-17.386,-17.396,-17.405,-17.415,-17.424,-17.433,-17.442,-17.451,
     3        -17.373,-17.388,-17.402,-17.414,-17.425,-17.436,-17.446,
     3-17.456,-17.466,-17.475,-17.484,-17.493,-17.502,-17.511,-17.520,
     4        -17.462,-17.476,-17.489,-17.500,-17.511,-17.521,-17.531,
     4-17.541,-17.550,-17.559,-17.569,-17.578,-17.587,-17.595,-17.604,
     5        -17.567,-17.581,-17.592,-17.603,-17.613,-17.623,-17.632,
     5-17.641,-17.651,-17.660,-17.669,-17.678,-17.686,-17.695,-17.704,
     6        -17.689,-17.701,-17.712,-17.722,-17.732,-17.741,-17.750,
     6-17.759,-17.768,-17.777,-17.786,-17.795,-17.803,-17.812,-17.821,
     7        -17.829,-17.840,-17.851,-17.860,-17.869,-17.878,-17.887,
     7-17.896,-17.904,-17.913,-17.922,-17.930,-17.939,-17.948,-17.956,
     8        -17.988,-18.000,-18.010,-18.019,-18.028,-18.036,-18.045,
     8-18.053,-18.062,-18.070,-18.079,-18.087,-18.096,-18.104,-18.112,
     9        -18.171,-18.183,-18.192,-18.201,-18.210,-18.218,-18.227,
     9-18.235,-18.243,-18.252,-18.260,-18.268,-18.277,-18.285,-18.293,
     A        -18.381,-18.393,-18.403,-18.413,-18.422,-18.430,-18.438,
     A-18.447,-18.455,-18.463,-18.471,-18.479,-18.487,-18.495,-18.503/
      DATAC12/-18.625,-18.638,-18.650,-18.660,-18.669,-18.678,-18.687,
     1-18.695,-18.703,-18.711,-18.719,-18.726,-18.734,-18.742,-18.750,
     2        -18.912,-18.929,-18.943,-18.955,-18.966,-18.975,-18.984,
     2-18.993,-19.001,-19.008,-19.016,-19.023,-19.031,-19.038,-19.045,
     3        -19.260,-19.283,-19.303,-19.320,-19.333,-19.345,-19.355,
     3-19.364,-19.372,-19.380,-19.387,-19.394,-19.400,-19.407,-19.413,
     4        -19.704,-19.740,-19.771,-19.796,-19.816,-19.832,-19.845,
     4-19.855,-19.863,-19.870,-19.876,-19.882,-19.887,-19.892,-19.897,
     5        -20.339,-20.386,-20.424,-20.454,-20.476,-20.492,-20.502,
     5-20.509,-20.513,-20.516,-20.518,-20.520,-20.521,-20.523,-20.524,
     6        -21.052,-21.075,-21.093,-21.105,-21.114,-21.120,-21.123,
     6-21.125,-21.126,-21.127,-21.128,-21.130,-21.131,-21.133,-21.135,
     7        -21.174,-21.203,-21.230,-21.255,-21.278,-21.299,-21.320,
     7-21.339,-21.357,-21.375,-21.392,-21.408,-21.424,-21.439,-21.454,
     8        -21.285,-21.317,-21.346,-21.372,-21.395,-21.416,-21.435,
     8-21.452,-21.468,-21.483,-21.497,-21.511,-21.524,-21.536,-21.548,
     9        -21.396,-21.429,-21.459,-21.486,-21.511,-21.532,-21.551,
     9-21.569,-21.585,-21.600,-21.614,-21.627,-21.640,-21.652,-21.663,
     A        -21.516,-21.549,-21.580,-21.609,-21.635,-21.658,-21.678,
     A-21.696,-21.713,-21.728,-21.742,-21.755,-21.767,-21.779,-21.790/
      DATAC13/-21.651,-21.681,-21.711,-21.738,-21.763,-21.785,-21.804,
     1-21.821,-21.837,-21.851,-21.864,-21.876,-21.887,-21.898,-21.908,
     2        -21.810,-21.831,-21.853,-21.874,-21.893,-21.910,-21.925,
     2-21.938,-21.950,-21.961,-21.971,-21.980,-21.989,-21.998,-22.006,
     3        -22.009,-22.016,-22.026,-22.037,-22.048,-22.058,-22.066,
     3-22.074,-22.081,-22.088,-22.094,-22.099,-22.105,-22.111,-22.117,
     4        -22.353,-22.317,-22.296,-22.284,-22.276,-22.270,-22.266,
     4-22.262,-22.260,-22.258,-22.257,-22.257,-22.257,-22.258,-22.259,
     5        -22.705,-22.609,-22.552,-22.515,-22.488,-22.468,-22.451,
     5-22.438,-22.427,-22.418,-22.410,-22.405,-22.400,-22.397,-22.395,
     6        -22.889,-22.791,-22.731,-22.690,-22.659,-22.634,-22.612,
     6-22.594,-22.579,-22.566,-22.555,-22.546,-22.539,-22.533,-22.528,
     7        -23.211,-23.109,-23.041,-22.989,-22.945,-22.906,-22.872,
     7-22.842,-22.816,-22.793,-22.774,-22.757,-22.743,-22.732,-22.722,
     8        -25.312,-24.669,-24.250,-23.959,-23.746,-23.587,-23.463,
     8-23.366,-23.288,-23.225,-23.173,-23.131,-23.095,-23.066,-23.041,
     9        -25.394,-24.752,-24.333,-24.041,-23.829,-23.669,-23.546,
     9-23.449,-23.371,-23.308,-23.256,-23.214,-23.178,-23.149,-23.124,
     A        -25.430,-24.787,-24.369,-24.077,-23.865,-23.705,-23.582,
     A-23.484,-23.407,-23.344,-23.292,-23.249,-23.214,-23.185,-23.160/
      DATA PARTOH/
     1   145.979,  178.033,  211.618,  247.053,  284.584,  324.398,
     2   366.639,  411.425,  458.854,  509.012,  561.976,  617.823,
     3   676.626,  738.448,  803.363,  871.437,  942.735, 1017.330,
     4  1095.284, 1176.654, 1261.510, 1349.898, 1441.875, 1537.483,
     5  1636.753, 1739.733, 1846.434, 1956.883, 2071.080, 2189.029,
     6  2310.724, 2436.155, 2565.283, 2698.103, 2834.571, 2974.627,
     7  3118.242, 3265.366, 3415.912, 3569.837, 3727.077/
      DATA FREQ1/0./
C
      SBFOH=0.
      IF(FR.EQ.FREQ1) GO TO 30
      FREQ1=FR
      WAVENO=FR/2.99792458E10
      EVOLT=WAVENO/8065.479
      N=int(EVOLT*10.-20.)
      EN=FLOAT(N)*.1+2.
      IF(N.LE.0) RETURN
      IF(N.GE.130) RETURN
      DO IT=1,15
         CROSSOHT(IT)=(CROSSOH(IT,N)+(CROSSOH(IT,N+1)-CROSSOH(IT,N))*
     *                (EVOLT-EN)*10.)
      END DO
c
c     interpolate to obtain partition function
c
   30 IF(T.GE.9000.) RETURN
      IF(N.LE.0) RETURN
      IF(N.GE.130) RETURN
      IT=int((T-1000.)*twhui+1.)
      IF(IT.LT.1) IT=1
      TN=FLOAT(IT)*twhu+800.
      PART=PARTOH(IT)+(PARTOH(IT+1)-PARTOH(IT))*(T-TN)*twhui
c
c     interpolate to obtain cross-section
c
      IT=int((T-2000.)*fihui+1.)
      IF(IT.LT.1) IT=1
      TN=FLOAT(IT)*fihu+1500.
      SBFOH=EXP((CROSSOHT(IT)+(CROSSOHT(IT+1)-CROSSOHT(IT))*
c    *     (T-TN)*fihui)*tenl)*PART
     *     (T-TN)*fihui)*tenl)
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE XENINI
C     =================
C
C     Initializes necessary arrays for evaluating hydrogen line profiles
C     from the XENOMORPH tables
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
C
      DO I=1,4
         DO J=1,22
            ILXEN(I,J)=0
         END DO
      END DO
      if(ihxenb.gt.0) then
         ihxenb=23
         ihxenr=ihxenb+1
         open(unit=ihxenb,file='xenomorph.blue.dat',status='old')
         open(unit=ihxenr,file='xenomorph.red.dat',status='old')
         write(6,641) ihxenb,ihxenr
       else
         return
      end if
c
  641 format(' -----------'/
     *       ' reading XENOMORPH tables; ihxen =',2i3,/
     *       ' -----------')
C
C ---------------------------------
C     read  tables - blue wing
C ---------------------------------
C
      ILINE=0
      READ(IHXENB,*) NTAB
      DO ITAB=1,NTAB
      ILINEB=ILINE
      READ(IHXENB,*) NLXEN
      DO ILI=1,NLXEN
         ILINE=ILINE+1
         READ(IHXENB,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,
     *                  NWL,NE,NT
         XNEMIN=ANEMIN
         ILXEN(I,J)=ILINE
         NWLXEN(ILINE)=NWL
         NTHXEN(ILINE)=NT
         NEHXEN(ILINE)=NE
         DO IWL=1,NWL
            ALXEN(ILINE,IWL)=ALMIN+(IWL-1)*DLA
         END DO
         DO INE=1,NE
            XNEXEN(INE,ILINE)=ANEMIN+(INE-1)*DLE
         END DO
         DO IT=1,NT
            XTXEN(IT,ILINE)=TMIN+(IT-1)*DLT
         END DO
      END DO
c
      DO ILI=1,NLXEN
         ILNE=ILINEB+ILI
         NWL=NWLXEN(ILNE)
         READ(IHXENB,500)
         DO INE=1,NEHXEN(ILNE)
            DO IT=1,NTHXEN(ILNE)
               READ(IHXENB,*) QLT,(PRFXB(ILNE,IWL,IT,INE),IWL=1,NWL)
            END DO
         END DO
C
      END DO
      END DO
  500 FORMAT(1X)
      CLOSE(IHXENB)
C
C ---------------------------------
C     read  tables - red wing
C ---------------------------------
C
      ILINE=0
      READ(IHXENR,*) NTAB
      DO ITAB=1,NTAB
      ILINEB=ILINE
      READ(IHXENR,*) NLXEN
      DO ILI=1,NLXEN
         ILINE=ILINE+1
         READ(IHXENR,*) I,J,ALMIN,ANEMIN,TMIN,DLA,DLE,DLT,
     *                  NWL,NE,NT
      END DO
c
      DO ILI=1,NLXEN
         ILNE=ILINEB+ILI
         NWL=NWLXEN(ILNE)
         READ(IHXENR,500)
         DO INE=1,NEHXEN(ILNE)
            DO IT=1,NTHXEN(ILNE)
               READ(IHXENR,*) QLT,(PRFXR(ILNE,IWL,IT,INE),IWL=1,NWL)
            END DO
         END DO
C
      END DO
      END DO
C
C     interpolation to the actual values of temperature and electron
C     density
C
      do id =1,nd
         tl=log10(temp(id))
         anel=log10(elec(id))
         do ili=1,nlxen
            iline=ilineb+ili
            nwl=nwlxen(iline)
            do iwl=1,nwl
               call intxen(prfb0,prfr0,tl,anel,iwl,iline)
               prfb(iline,id,iwl)=prfb0
               prfr(iline,id,iwl)=prfb0
            end do
         end do
      end do
      CLOSE(IHXENR)
c
      RETURN
      END
C
C
C ********************************************************************
C
C
      SUBROUTINE INTXEN(W0B,W0R,X0,Z0,IWL,ILINE)
C     ==========================================
C
C     Interpolation in temperature and electron density from the
C     Xenomorph tables for hydrogen lines to the actual valus of
C     temperature and electron density
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      DIMENSION ZZ(3),XX(3),WXB(3),WZB(3),WXR(3),WZR(3)
C
      NX=2
      NZ=2
      NT=NTHXEN(ILINE)
      NE=NEHXEN(ILINE)
C
      DO 10 IZZ=1,NE-1
         IPZ=IZZ
         IF(Z0.LE.XNEXEN(IZZ+1,ILINE)) GO TO 20
   10 CONTINUE
   20 N0Z=IPZ-NZ/2+1
      IF(N0Z.LT.1) N0Z=1
      IF(N0Z.GT.NE-NZ+1) N0Z=NE-NZ+1
      N1Z=N0Z+NZ-1
C
      DO 300 IZZ=N0Z,N1Z
         I0Z=IZZ-N0Z+1
         ZZ(I0Z)=XNEXEN(IZZ,ILINE)
         DO 30 IX=1,NT-1
            IPX=IX
            IF(X0.LE.XTXEN(IX+1,ILINE)) GO TO 40
   30    CONTINUE
   40    N0X=IPX-NX/2+1
         IF(N0X.LT.1) N0X=1
         IF(N0X.GT.NT-NX+1) N0X=NT-NX+1
         N1X=N0X+NX-1
         DO 200 IX=N0X,N1X
            I0=IX-N0X+1
            XX(I0)=XTXEN(IX,ILINE)
            WXB(I0)=PRFXB(ILINE,IWL,IX,IZZ)
            WXR(I0)=PRFXR(ILINE,IWL,IX,IZZ)
  200    CONTINUE
         WZB(I0Z)=YINT(XX,WXB,X0)
         WZR(I0Z)=YINT(XX,WXR,X0)
  300 CONTINUE
      W0B=YINT(ZZ,WZB,Z0)
      W0R=YINT(ZZ,WZR,Z0)
      RETURN
      END
C
C
C     ******************************************************************
C
C
      SUBROUTINE GOMINI
C     =================
C
C     Initialization and reading of the opacity table for thermal processe
C     and Rayleigh scattering
c     raytab: scattering opacities in cm^2/gm at 5.0872638d14 Hz (sodium D)
c     (NOTE: Quantities in rayleigh.tab are in log_e)
C
c     tempvec: array of temperatures
c     rhovec: array of densities (gm/cm^3)
c     nu:     array of frequencies
c     table:  absorptive opacities in cm^2/gm
c     (NOTE:  Quantities in absorption.tab are in log_e)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      COMMON/GOMOPA/frgtab(mfhtab),wlgtab(mfhtab),hydopg(mfhtab,mdepth),
     *              nugfreq
      common/gompar/hglim,ihgom
      dimension temvec(mtabth),elevec(mtabeh),
     *          hydcrs(mtabth,mtabeh,mfhtab)
c
      if(ihgom.eq.0) return
C
      open(53,file='gomhyd.dat',status='old')
c
      read(53,*) nugfreq,nugtemp,nugele
      read(53,*)
      read(53,*) (temvec(i),i=1,nugtemp)
      read(53,*)
      read(53,*) (elevec(j),j=1,nugele)
      do it=1,nugtemp
         temvec(it)=log(temvec(it)*1.161e4)
      end do
c     write(6,600) ihgom,nugfreq,nugtemp,nugele
c 600 format(' ihgom,nugfr,nugt,nuge ',4i4)
c
      EGTAB1 = elevec(1)
      EGTAB2 = elevec(nugele)
      TGTAB1 = temvec(1)
      TGTAB2 = temvec(nugtemp)
c
      do k = 1, nugfreq
         read(53,501) eneev
         frgtab(k)=3.28805e15/13.595*eneev
         wlgtab(k)=2.997925e18/frgtab(k)
         do i = 1, nugtemp
            read(53,*) (hydcrs(i,j,k),j=1,nugele)
         end do
      end do
      frg1=frgtab(1)
      frg2=frgtab(nugfreq)
c
  501 format(40x,f17.14)
      close(53)
C
c     Interpolate to the actual temperature and electron density
c     at the individual depth points
C
      do 10 id=1,nd
         if(elec(id).lt.HGLIM) go to 10
         rl=log(elec(id))
         tl=log(temp(id))
c
         DELTAR=(RL-EGTAB1)/(EGTAB2-EGTAB1)*FLOAT(nugele-1)
         JR = 1 + IDINT(DELTAR)
         IF(JR.LT.1) JR = 1
         IF(JR.GT.(nugele-1)) JR = nugele-1
         r1i=elevec(jr)
         r2i=elevec(jr+1)
         dri=(RL-R1i)/(R2i-R1i)
         if(JR .eq. 1) dri = 0.d0
C
         DELTAT=(TL-TGTAB1)/(TGTAB2-TGTAB1)*FLOAT(nugtemp-1)
         JP = 1 + IDINT(DELTAT)
         IF(JP.LT.1) JP = 1
         IF(JP.GT.nugtemp-1) JP = nugtemp-1
         t1i=temvec(jp)
         t2i=temvec(jp+1)
         dti=(TL-T1i)/(T2i-T1i)
         if(JP .eq. 1) dti = 0.d0
C
c        loop over tabular frequencies
c
         do jf=1,nugfreq
            opr1=hydcrs(jp,jr,jf)+dti*
     *           (hydcrs(jp+1,jr,jf)-hydcrs(jp,jr,jf))
            opr2=hydcrs(jp,jr+1,jf)+dti*
     *           (hydcrs(jp+1,jr+1,jf)-hydcrs(jp,jr+1,jf))
            opac=opr1+dri*(opr2-opr1)
            hydopg(jf,id)=opac+log(0.02654*4.1347e-15)
         end do
   10 continue
      return
      end
C
C     ****************************************************\
C
C
      subroutine ghydop(id,i0,i1,pj,absoh,emish)
c     ==========================================
c
c     hydrogen opacity -- lines + pseudocontinuum from Gomez tables
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      COMMON/GOMOPA/frgtab(mfhtab),wlgtab(mfhtab),hydopg(mfhtab,mdepth),
     *              nugfreq
      dimension absoh(mfreq),emish(mfreq),pj(40)
c
      frg1=frgtab(1)
      frg2=frgtab(nugfreq)
      do 20 ij=i0,i1
         fr=freq(ij)
         if(fr.lt.frg1.or.fr.gt.frg2) go to 20
         wla=2.997925e18/fr
         frl=log10(fr)
c
         if(ij.eq.i0) igf=nugfreq
   10    continue
         if(wla.gt.wlgtab(igf)) then
            igf=igf-1
            go to 10
         end if
         ig0=igf
         if(ig0.le.2) ig0=2
         ig1=igf-1
         abl=(hydopg(ig1,id)-hydopg(ig0,id))*(wla-wlgtab(ig0))/
     *       (wlgtab(ig1)-wlgtab(ig0))+hydopg(ig0,id)
c
         ii=1
         if(freq(ij).gt.8.22013e14) then
            pp=pj(1)*2.
          else
            pp=pj(2)*8.
         end if
c
         F15=FR*1.E-15
         XKF=EXP(-4.79928e-11*FR/TEMP(ID))
         XKFB=XKF*1.4743E-2*F15*F15*F15

         oph=exp(abl)*pp
         absoh(ij)=absoh(ij)+oph
         emish(ij)=emish(ij)+oph*xkfb/(1.-xkf)
   20 continue
c
      return
      end

C
C ********************************************************************
C
      subroutine ingrid(mode,inext,igrd)
C     ==================================
C
c     setting state parameters for the opacity grid calculations
c
c     input:
c           temp1 - lowest value of T
c           temp2 - largest value of T
c           ntemp - number of temperature values
c           dens1 - lowest value of the density parameter
c           dens2 - largest value of the density parameter
c           ndens - number of the density parameter values
c
c           isdens = 0 - density parameter is electron density
c                  > 0 - density parameter is mass density
c                  < 0 - density parameter is gas pressure
c
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      parameter (un=1.,ten15=1.e-15,c18=2.997925e18)
      real*4 absgrd(mttab,mrtab,mfgrid),dtim
      common/alsave/ALAM0s,ALASTs,CUTOF0s,CUTOFSs,RELOPs,SPACEs
      common/gridp0/tempg(mttab),densg(mrtab),elecgr(mttab,mrtab),
     *              densg0(mttab),temp1,ntemp,ndens
      common/gridf0/wlgrid(mfgrid),nfgrid
      common/fintab/absgrd
      common/prfrgr/ipfreq,indext,indexn
      common/igrddd/igrdd,irelin
      common/initab/absop(msftab),wltab(msftab),
     *       nfrtab(mttab,mrtab),inttab 
      common/elecm0/elecm(mdepth)
      common/timeta/dtim
      common/relabu/relabn(matom),popul0(mlevel,1)
      dimension abgrd(mfgrid),xli(3)
      character*(80) tabname
      common/tabout/tabname,ibingr
c
c     --------------
c     initialization
c     --------------
c
      igrdd=igrd
      if(mode.eq.0) then
c
      read(2,*) ntemp,temp1,temp2
      read(2,*) idens
      read(2,*) ndens,dens1,dens2
      if(ifeos.eq.0) read(2,*) nfgrid,inttab,wlam1,wlam2
      read(2,*) tabname,ibingr
c
      irsct=0
      irsche=0
      irsch2=0
c
      wl1=log(wlam1)
      wl2=log(wlam2)
      dwl=(wl2-wl1)/(nfgrid-1)
      do i=1,nfgrid
         wlgrid(i)=exp(wl1+(i-1)*dwl)
      end do
c
      if(temp1.gt.0.) then
      at1=log(temp1)
      at2=log(temp2)
      dt=0.
      if(ntemp.gt.1) dt=(at2-at1)/(ntemp-1)
      do i=1,ntemp
         tempg(i)=exp(at1+(i-1)*dt)
      end do
c
      at1=log(dens1)
      at2=log(dens2)
      dt=0.
      if(ndens.gt.1) dt=(at2-at1)/(ndens-1)
      do i=1,ndens
         densg(i)=exp(at1+(i-1)*dt)
      end do
c
      write(6,621) ntemp
      write(6,622) (tempg(i),i=1,ntemp)
      write(6,623) ndens
      write(6,622) (densg(i),i=1,ndens)
  621 format(/' COMPUTING AN OPACITY TABLE WITH GRID PARAMETERS:'/
     *' ===== ntemp, temp ',i4)
  622 format(1p8e10.2)
  623 format(/' ===== ndens, dens',i4)
      else
      call inpmod
      ntemp=nd
      ndens=1
      do it=1,ntemp
         tempg(it)=temp(it)
         densg0(it)=dens(it)
         elecm(it)=elec(it)
      end do
      densg(1)=densg0(1)
      if(ifeos.eq.0) then
      write(6,621) ntemp
      write(6,622) (tempg(i),i=1,ntemp)
      write(6,623) ndens
      write(6,622) (densg0(i),i=1,ntemp)
      end if
      ndens=1
      idens=2
      end if
c
      nd=1
      idstd=1
      inext=1
      frmx=0.
      frmn=1.e20
c
      indext=1
      indexn=1
      ipfreq=0
      irelin=1
      temp(1)=tempg(indext)
c
      write(6,646) indext,temp(1),
     *             indexn,densg(indexn)
  646 format(/' ************************************',
     *       /' GRID POINT OF THE OPACITY TABLE WITH:'/
     *        ' INDEX TEMP, T   ',i4,f10.1/
     *        ' INDEX DENS, DENS',I4,1PE10.1,
     *       /' ************************************'/)
c
      if(temp1.le.0.) elec(1)=elecm(indext)
      call densit(densg(indexn),idens)
      if(ntemp.eq.1.and.ndens.eq.1) inext=0
   
      elecgr(indext,indexn)=elec(1)
  
      call abnchn(0)

      return
c
c     ---------------------------------------------
c     after computing the table for one T-rho pair:
c     ---------------------------------------------
c
      else if(mode.eq.1) then
      if(ifeos.eq.0) then
c
      call timing(1,igrd+1)
c
      do i=1,3
         xli(i)=0.
      end do
      do i=1,nmlist
         xli(i)=float(nlinmt(i))*1.e-3
      end do
c        
      if(imode.ge.-5) then
         if(indext.eq.1.and.indexn.eq.1)
     *   write(29,625)       
         write(29,626) indext,indexn,temp(1),dens(1),elec(1),
     *              float(nlin0)*1.e-3,
     *              (xli(i),i=1,3),dtim
  625 format('  it   ir     t       rho      elec',6x,
     *   ' atomic   molec1  molec2  molec3      time'/) 
  626 format(2i4,f9.2,1p2e10.2,2x,0pf8.1,2x,3f8.1,2x,f8.2)
       else
         alam0=alam0s
         if(alam0s.eq.0.) alam0=5.e7/temp(1)/10.
         if(alam0s.lt.0.) alam0=-5.e7/temp(1)/alam0s
         alast=alasts
         if(alasts.eq.0.) alast=5.e7/temp(1)*20.
         if(alasts.lt.0.) alast=-5.e7/temp(1)*alasts
         if(alast.gt.1.e5) alast=1.e5
         write(29,629) temp(1),elec(1),dens(1),
     *                 alam0,alast   
      end if
  629 format(1p3e11.3,0pf9.3,0pf12.3)
c
c     ------------------------------------------------
c     interpolate and store previously computed table
c     ------------------------------------------------
c
      nfr=ipfreq
      nfrtab(indext,indexn)=ipfreq
      write(*,*) 'indext,indexn,nfreq',indext,indexn,ipfreq 
      write(*,*) 'nfr,nfgrid',nfr,nfgrid
c
      if(inttab.eq.1) then
c        call interp(wltab,absop,wlgrid,abgrd,nfr,nfgrid,2,0,0)
         call intrp(wltab,absop,wlgrid,abgrd,nfr,nfgrid)
       else
         ij=0
         ijgrd=0
   30    continue
         ijgrd=ijgrd+1
         wlgr=0.5*(wlgrid(ijgrd)+wlgrid(ijgrd+1))
         isum=0
         sum=0.
   40    continue
         ij=ij+1
         if(ij.gt.nfr) go to 50
         wlt=wltab(ij)
         abl=absop(ij)
         if(wlt.le.wlgr) then
            sum=sum+exp(abl)
            isum=isum+1
            go to 40
         end if
         if(isum.gt.0) then
            abgrd(ijgrd)=log(sum/float(isum))
          else
            abg=abl+(absop(ij+1)-abl)/(wltab(ij+1)-wlt)*(wlgr-wlt)
            abgrd(ijgrd)=abg
c           write(*,*) 'grd',ij,absop(ij+1),abl,wltab(ij+1),
c    *                 wlt,wlgr,abg,abgrd(ijgrd),ijgrd
         end if
         if(ijgrd.lt.nfgrid) then
            ij=ij-1
            go to 30
          else if(ijgrd.eq.nfgrid) then
            wlgr=wlgrid(nfgrid)
            sum=0.
            isum=0
            if(ij.lt.nfr) ij=ij-1
            go to 40
         end if
      end if
   50 continue
c
      do ij=1,nfgrid
         absgrd(indext,indexn,ij)=real(abgrd(ij))
c        write(28,628) wlgrid(ij),abgrd(ij)
c 628 format(2f15.5)
      end do
      absgrd(indext,indexn,nfgrid)=absgrd(indext,indexn,nfgrid-1)
      end if
c
c     ------------------------------
c     prepare values for a new table
c     ------------------------------
c
      ipfreq=0
      if(indexn.lt.ndens) then
         indexn=indexn+1
         rho=densg(indexn)
         write(6,646) indext,tempg(indext),
     *                indexn,densg(indexn)
         call densit(rho,idens)
         inext=1
       else 
         indexn=1
         irelin=1
         if(indext.lt.ntemp) then
            indext=indext+1
            temp(1)=tempg(indext)
            if(temp1.le.0.) then
               densg(indexn)=densg0(indext)
               elec(1)=elecm(indext)
            end if
            rho=densg(indexn)
            write(6,646) indext,tempg(indext),
     *                   indexn,densg(indexn)
            call densit(rho,idens)
            inext=1
          else
            inext=0
         end if
      end if
      if(inext.eq.1) then
         rewind(19)
         if(inlist.lt.0) rewind(19)
      end if
c
      elecgr(indext,indexn)=elec(1)
c
      call abnchn(0)
      id=1
         do i=1,4
            do j=i+1,22
               call hydtab(i,j,id)
            end do
         end do
      end if
c
      return
      end
C 
C
C ********************************************************************
C
C
      subroutine ougrid(abso,scat)
C     ============================
C
C     output of grid opacities
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      common/prfrgr/ipfreq,indext,indexn
      common/gridf0/wlgrid(mfgrid),nfgrid
      common/initab/absop(msftab),wltab(msftab),
     *              nfrtab(mttab,mrtab),inttab 
      parameter (un=1.,ten15=1.e-15,c18=2.997925e18)
      DIMENSION ABSO(MFREQ),SCAT(MFREQ)
c
      d1=un/dens(1)
      if (nfreq.le.3) return 
c
      if(iprin.lt.4) then
         do ij=3,nfreq-1
            abl=log(abso(ij)*d1)
            ipfreq=ipfreq+1
            absop(ipfreq)=abl
            wltab(ipfreq)=2.997925e18/freq(ij)
         end do
       else
         do ij=3,nfreq-1
            abl=log(abso(ij)*d1)
            ipfreq=ipfreq+1
            write(27,637) ipfreq,c18/freq(ij),abl
            absop(ipfreq)=abl
            wltab(ipfreq)=2.997925e18/freq(ij)
         end do
      end if
  637 format(i10,f14.5,0pf12.5)    
c
      return
      end
C 
C
C ********************************************************************
C
C
      subroutine fingrd
c     =================
c
c     storing the complete, interpolated, opacity table 
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      real*4 absgrd(mttab,mrtab,mfgrid)
      common/gridp0/tempg(mttab),densg(mrtab),elecgr(mttab,mrtab),
     *              densg0(mttab),temp1,ntemp,ndens
      common/gridf0/wlgrid(mfgrid),nfgrid
      common/fintab/absgrd
      common/relabu/relabn(matom),popul0(mlevel,1)
      character*(80) tabname
      common/tabout/tabname,ibingr
c
      if(ifeos.ne.0) return
c     read(2,*) tabname,ibingr
c
      iophmp=iophmi
      if(ielhm.gt.0.and.relabn(1).gt.0.) iophmp=1 
      if(ibingr.eq.0) then
         open(53,file=tabname,status='unknown')
         write(53,600) 
         do iat=1,30
            write(53,601) typat(iat),abnd(iat),abnd(iat)*relabn(iat)
         end do
         write(53,602) ifmol,tmolim
         write(53,603) iophmp,ioph2p,iophem,iopch,iopoh,ioph2m,
     *                 ioh2h2,ioh2he,ioh2h1,iohhe
         write(53,611) nfgrid,ntemp,ndens
         write(53,612) (log(tempg(i)),i=1,ntemp)
         write(53,613) (log(densg(j)),j=1,ndens)
         write(53,614) ((log(elecgr(i,j)),j=1,ndens),i=1,ntemp)
         do k = 1, nfgrid
            write(53,615) k,wlgrid(k),2.997925e18/wlgrid(k)
            do j = 1,ndens
               write(53,616) (absgrd(i,j,k),i=1,ntemp)
            end do
         end do
  600    format('opacity table with element abundances:'/
     *          'element   for EOS   for opacities')
  601    format('  ',a4,1p2e12.3)
  602    format(/'molecules - ifmol,tmolim:'/,i4,f10.1)
  603    format('additional opacities'/
     *   '  H-  H2+ He- CH  OH  H2- CIA: H2H2 H2He H2H  HHe'/
     *          6i4,4x,4i4)
  611    format(/'number of frequencies, temperatures, densities:'
     *          /10x,3i10)
  612    format('log temperatures'/(6F11.6))
  613    format('log densities'/(6F11.6))
  614    format('log electron densities from EOS'/(6f11.6))
  615    format(/' *** frequency # : ',i8,f15.5/1pe20.8)
  616    format((1p6e14.6))
       end if
         do iat=1,30
            write(63) typat(iat),abnd(iat),abnd(iat)*relabn(iat)
         end do
         write(63) ifmol,tmolim
         write(63) iophmp,ioph2p,iophem,iopch,iopoh,ioph2m,
     *             ioh2h2,ioh2he,ioh2h1,iohhe
         write(63) nfgrid,ntemp,ndens
         write(63) (log(tempg(i)),i=1,ntemp)
         write(63) (log(densg(j)),j=1,ndens)
         write(63) ((log(elecgr(i,j)),j=1,ndens),i=1,ntemp)
         do k = 1, nfgrid
            write(63) wlgrid(k)
            do j = 1, ndens
               write(63) (absgrd(i,j,k),i=1,ntemp)
            end do
         end do
c     end if
c
      return
      end
c
c
c     *************************************************************
c
c
      subroutine abnchn(mode)
c     =======================
c
c     changing abundances (eliminating) species for an
c     evaluating an opacity table
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      common/relabu/relabn(matom),popul0(mlevel,1)
      data iread/1/
c
      if(iread.eq.1) then
         do ia=1,matom
            relabn(ia)=1.
         end do
   10    continue
         read(2,*,err=20,end=20) iatom,rela
         relabn(iatom)=rela
         write(*,*) 'ABUNDANCES CHANGED (AT.NUMBER, ABUND):',iatom,rela
         go to 10
   20    continue
         if(relabn(1).eq.0.) then
            iophmi=0
            ioph2p=0
         end if
         iread=0
      end if
c
      if(mode.eq.0) then
         do iat=1,natom
            do ii=n0a(iat),nka(iat)
               popul0(ii,1)=popul(ii,1)
            end do
         end do
         return
      end if
c
      do iat=1,natom
         ia=numat(iat)
         do ii=n0a(iat),nka(iat)
            popul(ii,1)=popul0(ii,1)*relabn(ia)
         end do
      end do
c
      do ia=1,matom
         do io=1,mion0
            rrr(1,io,ia)=rrr(1,io,ia)*relabn(ia)
         end do
      end do
c
      return
      end
c
c
c     *************************************************************
c
c

      subroutine densit(rho,idens)    
C     ============================
C
C     determining the state parameters for the opacity grid
C     calculations
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'      
      DIMENSION ES(MLEVEL,MLEVEL),BS(MLEVEL),POPLTE(MLEVEL)
c
      id=1
      dm(id)=0.
      IF(IFMOL.EQ.0.OR.TEMP(ID).GT.TMOLIM) 
     *   WMM(ID)=WMY(ID)*HMASS/YTOT(ID)
         if(idens.eq.0) then
            ELEC(ID)=rho
            ane=elec(id)
            call todens(id,temp(id),an,ane)
            DENS(ID)=(an-ane)*wmm(id)
            p=an*bolk*temp(id)
c           WRITE(6,602) ID,TEMP(ID),DENS(ID),ELEC(ID)
          else if(idens.lt.0) then
            AN=rho/TEMP(ID)/BOLK
            CALL ELDENS(ID,TEMP(ID),AN,ANE)
            ELEC(ID)=ANE
            DENS(ID)=WMM(ID)*(AN-ELEC(ID)) 
c           WRITE(6,601) ID,TEMP(ID),DENS(ID),ELEC(ID),ane0,an
          else if(idens.eq.1) then
            DENS(ID)=RHO
            CALL RHONEN(ID,TEMP(ID),RHO,AN,ANE)
            ELEC(ID)=ANE
            DENS(ID)=RHO
            rho0=WMM(ID)*(AN-ANE) 
c           WRITE(6,601) IDens,TEMP(ID),DENS(ID),ane,rho0,an
          else if(idens.eq.2) then
            CALL RHONEN(ID,TEMP(ID),RHO,AN,ANE)
            DENS(ID)=RHO
            ANE=ELEC(ID)
            rho0=WMM(ID)*(AN-ANE) 
c           WRITE(6,601) idens,TEMP(ID),DENS(ID),ane,rho0,an
         end if 
c 601 FORMAT(' **densit** t,rho,ne,rho0,an',I3,0PF10.1,1P5D11.3)
c 602 FORMAT(' **densit** t,rho,ne',I3,0PF10.1,1P5D11.3)


      CALL INIMOD
c
      CALL WNSTOR(ID)
      CALL SABOLF(ID)
      CALL RATMAT(ID,ES,BS)
      CALL LEVSOL(ES,BS,POPLTE,NLEVEL) 
      DO J=1,NLEVEL
         POPUL(J,ID)=POPLTE(J)
      END DO
c      
      return
      end


C
C ********************************************************************
C
 
      SUBROUTINE TODENS(ID,T,AN,ANE)
C     ==============================
C
C     determines AN (and ANP, AHTOT, and AHMOL) from T and ANE
C
C     Input parameters:
C     T    - temperature
C     ANE  - electron number density
C
C     Output:
C     AN    - total particle density
C     ANP   - proton number density
C     AHTOT - total hydrogen number density
C     AHMOL - relative number of hydrogen molecules with respect to the
C             total number of hydrogens
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      common/hydmol/anhmi,ahmol,ih2,ih2p,ihm
      parameter (un=1.d0,two=2.d0,half=0.5d0)
C
      QM=0.
      Q2=0.
      QP=0.
      Q=0.
      DQN=0.
      TK=BOLK*T
      THET=5.0404D3/T
C
C     Coefficients entering ionization (dissociation) balance of:
C     atomic hydrogen          - QH;
C     negative hydrogen ion    - QM   (considered only if IHM>0);
C     hydrogen molecule        - QP   (considered only if IH2>0);
C     ion of hydrogen molecule - Q2   (considered only if IH2P>0).
C
      IF(T.LE.8000.)  QM=1.0353D-16/T/SQRT(T)*EXP(8762.9/T)
      IF(T.le.8000.) QP=TK*EXP((-11.206998+THET*(2.7942767+THET*
     *             (0.079196803-0.024790744*THET)))*2.30258509299405)
      IF(T.LE.8000.) Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*
     *            (-0.056191273+0.0032687661*THET)))*2.30258509299405)
      QH=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405)
C
C     procedure STATE determines Q (and DQN) - the total charge (and its
C     derivative wrt temperature) due to ionization of all atoms which
C     are considered (both explicit and non-explicit), by solving the set
C     of Saha equations for the current values of T and ANE
C
      CALL STATE(ID,T,ANE,Q)
C
C     Auxiliary parameters for evaluating the elements of matrix of
C     linearized equations.
C     Note that complexity of the matrix depends on whether the hydrogen
C     molecule is taken into account
C     Treatment of hydrogen ionization-dissociation is based on
C     Mihalas, in Methods in Comput. Phys. 7, p.10 (1967)
C
      G2=QH/ANE
      G3=0.
      G4=0.
      G5=0.
      D=0.
      E=0.
      G3=QM*ANE
      A=UN+G2+G3
      D=G2-G3
      IF(IT.GT.1) GO TO 60
      IF(IH2.EQ.0.AND.IH2P.EQ.0) GO TO 40
      IF(IH2.EQ.0) GO TO 20
      E=G2*QP/Q2
      B=TWO*(UN+E)
      GG=ANE*Q2
      GO TO 30
   20 B=TWO
      E=UN
      GG=G2*ANE*QP
   30 C1=B*(GG*B+A*D)-E*A*A
      C2=A*(TWO*E+B*Q)-D*B
      C3=-E-B*Q
      F1=(SQRT(C2*C2-4.*C1*C3)-C2)*HALF/C1
      FE=F1*D+E*(UN-A*F1)/B+Q
      GO TO 50
   40 F1=UN/A
      FE=D/A+Q
   50 AH=ANE/FE
      ANH=AH*F1
   60 AE=ANH/ANE
      GG=AE*QP
      E=ANH*Q2
      B=ANH*QM
C
c      S(1)=AN-ANE-YTOT(ID)*AH
c      S(2)=ANH*(D+GG)+Q*AH-ANE
c      S(3)=AH-ANH*(A+TWO*(E+GG))
c
      hhn=A+TWO*(E+GG)
      anh=ane/(d+gg+q*hhn)
      ah=anh*hhn
      an=ane+ytot(id)*ah
C
      AHTOT=AH
      AHMOL=TWO*ANH*(ANH*Q2+ANH/ANE*QP)/AH
      ANP=ANH/ANE*QH
      RETURN
      END
C
C
C ***********************************************************************
C
 
      subroutine rhonen(id,t,rho,an,ane)
c     ==================================
c
c     iterative determination of N and Ne from given T and RHO
c
C
C     Input:  T   - temperature
C             RHO - mass density
C     Output: AN  - total particle density
C             ANE - elctron density
C
      INCLUDE 'PARAMS.FOR'
      common/nerela/anerel
c
      it=0
      if(id.eq.1.and.anerel.eq.0.) then
         anerel=0.5
         if(t.lt.9000.) anerel=0.4
         if(t.lt.8000.) anerel=0.1
         if(t.lt.7000.) anerel=0.01
         if(t.lt.6000.) anerel=0.001
         if(t.lt.5500.) anerel=0.0001
c         if(t.lt.5000.) anerel=1.e-5
c         if(t.lt.4000.) anerel=1.e-6
      end if
   10 continue
      it=it+1
      an=rho/wmm(id)/(1.d0-anerel)
      ane0=anerel*an
      call eldens(id,t,an,ane)
      anerel=ane/an
      write(6,602) it,id,t,rho,an,ane,wmm(id),anerel
  602 format(/' **** rhonen it,id,t,r,N,Ne,wmm,ner',2i4,f7.0,1p5e11.4)
      if(abs((ane-ane0)/ane0).lt.1.e-5) go to 20
      if(it.lt.50) go to 10
c      write(6,601) an,ane,ane0
c 601 format(/' slow convergence of RHONEN - N,Ne,Nep=',1p3e11.3)
   20 continue
c
      return
      end   
C
C ********************************************************************
C

      SUBROUTINE ELDENS(ID,T,AN,ANE)
C     ==============================
C
C     Evaluation of the electron density and the total hydrogen
C     number density for a given total particle number density
C     and temperature;
C     by solving the set of Saha equations, charge conservation and
C     particle conservation equations (by a Newton-Raphson method)
C
C     Input parameters:
C     T    - temperature
C     AN   - total particle number density
C
C     Output:
C     ANE   - electron density
C     ANP   - proton number density
C     AHTOT - total hydrogen number density
C     AHMOL - relativer number of hydrogen molecules with respect to the
C             total number of hydrogens
C     ENERG - part of the internal energy: excitation and ionization
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      common/hydmol/anhmi,ahmol,ih2,ih2p,ihm
      common/hydato/ah,anh,anp
      common/nerela/anerel
      parameter (un=1.d0,two=2.d0,half=0.5d0)
      DIMENSION R(3,3),S(3),P(3)
C
      TK=BOLK*T
      if(ifmol.gt.0.and.t.lt.tmolim) then
         aein=an*anerel
         call moleq(id,t,an,aein,ane,0)
         anerel=ane/an
         return
      end if
c
      QM=0.
      Q2=0.
      QP=0.
      Q=0.
      DQN=0.
      TK=BOLK*T
      THET=5.0404D3/T
C
C     Coefficients entering ionization (dissociation) balance of:
C     atomic hydrogen          - QH;
C     negative hydrogen ion    - QM   (considered only if IHM>0);
C     hydrogen molecule        - QP   (considered only if IH2>0);
C     ion of hydrogen molecule - Q2   (considered only if IH2P>0).
C
      IF(IATREF.EQ.IATH) THEN
c      IF(T.LE.9000.)  QM=1.0353D-16/T/SQRT(T)*EXP(8762.9/T)
c      IF(T.le.9000.) QP=TK*EXP((-11.206998+THET*(2.7942767+THET*
c    *   (0.079196803-0.024790744*THET)))*2.30258509299405)
c      IF(T.LE.9000.) Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*
c    *   (-0.056191273+0.0032687661*THET)))*2.30258509299405)
       QM=1.0353D-16/T/SQRT(T)*EXP(8762.9/T)
       QP=TK*EXP((-11.206998+THET*(2.7942767+THET*
     *   (0.079196803-0.024790744*THET)))*2.30258509299405)
       Q2=TK*EXP((-12.533505+THET*(4.9251644+THET*
     *   (-0.056191273+0.0032687661*THET)))*2.30258509299405)
       QH0=EXP((15.38287+1.5*LOG10(T)-13.595*THET)*2.30258509299405)
      END IF
C
C     Initial estimate of the electron density
C
      if(anerel.le.0.) then
         if(t.gt.1.e4) then
            anerel=0.5
          else
            if(elec(id).gt.0..and.dens(id).gt.0.) then
               anerel=elec(id)/(elec(id)+dens(id)/wmm(id))
             else
               anerel=0.1
            end if
         end if
      end if
c
      ANE=AN*ANEREL
      IT=0
C
C     Basic Newton-Raphson loop - solution of the non-linear set
C     for the unknown vector P, consistiong of AH, ANH (neutral
C     hydrogen number density) and ANE.
C
   10 IT=IT+1
C
C     procedure STATE determines Q (and DQN) - the total charge (and its
C     derivative wrt temperature) due to ionization of all atoms which
C     are considered (both explicit and non-explicit), by solving the set
C     of Saha equations for the current values of T and ANE
C
      CALL STATE(ID,T,ANE,Q)
      QH=QH0*2./PFSTD(1,1)
C
C     Auxiliary parameters for evaluating the elements of matrix of
C     linearized equations.
C     Note that complexity of the matrix depends on whether the hydrogen
C     molecule is taken into account
C     Treatment of hydrogen ionization-dissociation is based on
C     Mihalas, in Methods in Comput. Phys. 7, p.10 (1967)
C
      IF(IATREF.EQ.IATH) THEN
      G2=QH/ANE
      G3=0.
      G4=0.
      G5=0.
      D=0.
      E=0.
      G3=QM*ANE
      A=UN+G2+G3
      D=G2-G3
      IF(IT.GT.1) GO TO 60
      IF(IH2.EQ.0.AND.IH2P.EQ.0) GO TO 40
      IF(IH2.EQ.0) GO TO 20
      E=G2*QP/Q2
      B=TWO*(UN+E)
      GG=ANE*Q2
      GO TO 30
   20 B=TWO
      E=UN
      GG=G2*ANE*QP
   30 C1=B*(GG*B+A*D)-E*A*A
      C2=A*(TWO*E+B*Q)-D*B
      C3=-E-B*Q
      F1=(SQRT(C2*C2-4.*C1*C3)-C2)*HALF/C1
      FE=F1*D+E*(UN-A*F1)/B+Q
      GO TO 50
   40 F1=UN/A
      FE=D/A+Q
   50 AH=ANE/FE
      if(ah.gt.3.*an) ah=dens(id)/hmass/wmy(id)
      ANH=AH*F1
   60 AE=ANH/ANE
      GG=AE*QP
      E=ANH*Q2
      B=ANH*QM
C
C     Matrix of the linearized system R, and the rhs vector S
C
      R(1,1)=YTOT(ID)
c     R(1,2)=0.
      r(1,2)=-two*(anh*q2+gg)
      R(1,3)=UN
      R(2,1)=-Q
      R(2,2)=-D-TWO*GG
      R(2,3)=UN+B+AE*(G2+GG)-DQN*AH
      R(3,1)=-UN
      R(3,2)=A+4.*(anh*q2+GG)
      R(3,3)=B-AE*(G2+TWO*GG)
      S(1)=AN-ANE-YTOT(ID)*AH+anh*(anh*q2+gg)
      S(2)=ANH*(D+GG)+Q*AH-ANE
      S(3)=AH-ANH*(A+TWO*(anh*q2+GG))
C
C     Solution of the linearized equations for the correction vector P
C
      CALL LINEQS(R,S,P,3,3)
C
C     New values of AH, ANH, and ANE
C
      AH=AH+P(1)
      ANH=ANH+P(2)
      DELNE=P(3)
      ANE=ANE+DELNE
C
C     hydrogen is not the reference atom
C
      ELSE
C
C     Matrix of the linearized system R, and the rhs vector S
C
      IF(IT.EQ.1) THEN
         ANE=AN*HALF
         AH=ANE/YTOT(ID)
      END IF
      R(1,1)=YTOT(ID)
      R(1,2)=UN
      R(2,1)=-Q-QREF
      R(2,2)=UN-(DQN+DQNR)*AH
      S(1)=AN-ANE-YTOT(ID)*AH
      S(2)=(Q+QREF)*AH-ANE
C
C     Solution of the linearized equations for the correction vector P
C
      CALL LINEQS(R,S,P,2,3)
      AH=AH+P(1)
      DELNE=P(2)
      ANE=ANE+DELNE
      END IF
C
C     Convergence criterion
C
      IF(ANE.LE.0.) ANE=1.D-7*AN
      IF(ABS(DELNE/ANE).GT.1.D-6.AND.IT.LE.20) GO TO 10
C
C     ANEREL is the exact ratio betwen electron density and total
C     particle density, which is going to be used in the subseguent
C     call of ELDENS
C
      ANEREL=ANE/AN
      AHTOT=AH
      IF(IATREF.EQ.IATH) THEN
c     AHMOL=TWO*ANH*(ANH*Q2+ANH/ANE*QP)/AH
      AHMOL=ANH*ANH*Q2
      ANP=ANH/ANE*QH
      ANHMI=ANH*ANE*QM
      anhn=anh+anp+anhmi+2.*ahmol
      wmm(id)=wmy(id)/(ytot(id)-ahmol/anhn)*hmass
      END IF
C
      RETURN
      END
C

 
C
C ********************************************************************
C
C
      SUBROUTINE TIMING(MOD,ITER)
C     ===========================
C
C     Timing procedure (call machine dependent routine!!)
C
      CHARACTER ROUT*6
      dimension dummy(2)
      common/timeta/dtim
      DATA T0/0./
      SAVE T0
C
      TIME=etime(dummy)
      DT=TIME-T0
      T0=TIME
      IP=ITER
      IF(MOD.EQ.1) THEN
         ROUT=' TABLE'
      ELSE IF(MOD.EQ.2) THEN
         ROUT=' FINAL'
      ENDIF
      WRITE(69,600) IP,TIME,DT,ROUT
      dtim=dt
  600 FORMAT(I6,2F11.2,2X,A6)
      RETURN
      END
C
C
C ********************************************************************
C
C
      subroutine eospri
c     =================
c
c     Outprint of Equation of State parameters
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      common/moltst/pfmol(500,mdepth),anmol(500,mdepth),
     *              pfato(100,mdepth),anato(100,mdepth),
     *              pfion(100,mdepth),anion(100,mdepth)
      common/hydmol/anhmi,ahmol,ih2,ih2p,ihm
      common/hydato/ah,anh,anp
      dimension hpo(mdepth),nelemx(38)
      dimension amh2(5)
      data nelemx/ 1, 2, 3, 4, 5, 6, 7, 8, 9,
     *            11,12,13,14,15,16,17,19,20,
     *            21,22,23,24,25,26,28,29,32,
     *            35,37,38,39,40,41,53,56,57,58,60/
      data amh2/1.13390E+01,-2.97499E+00,4.10842E-02,-3.58550E-03,
     *          1.31844E-04/
      data init/1/
c
c     id=idstd

      do id=1,nd

      t=temp(id)
      ane=elec(id)
      rho=dens(id)
      ann = dens(id)/wmm(id)+elec(id)
c
      if(ifmol.eq.0.or.t.gt.tmolim) then
         it=0
   10    continue
         ann0=ann
         it=it+1
         call eldens(id,t,ann,ane)
         anmol(1,id)=anhmi
         anmol(2,id)=ahmol
         anato(1,id)=anh
         anion(1,id)=anp
         hpop=dens(id)/wmy(id)/hmass
         do i=1,nmetal
            j=nelemx(i)
            anato(j,id)=anato(j,id)*hpop
            anion(j,id)=anion(j,id)*hpop
         end do
         anato(1,id)=anh
         anion(1,id)=anp
c        wmm(id)=(wmy(id)+2.*anmol(2,id)/hpop)/ytot(id)*hmass
         wmm(id)=wmy(id)/(ytot(id)-anmol(2,id)/hpop)*hmass
         ann=dens(id)/wmm(id)+ane
         if((ann-ann0)/ann0.gt.1.e-5) go to 10
      end if
c
         nmetal=38
         write(*,*) ''
         write(*,*) 'atomic number densities and partition functions'
         write(*,*) ''
         atot=0.
         do i=1,nmetal
            j=nelemx(i)
            if(j.le.28)
     *         write(6,621) j,typat(j),anato(j,id),pfato(j,id)
            atot=atot+anato(j,id)
         end do
         write(*,*) ''
         write(*,*) 'ionic number densities and partition functions'
         write(*,*) ''
         ctot=0.
         do i=1,nmetal
            j=nelemx(i)
            if(j.le.28)
     *         write(6,622) j,typat(j),anion(j,id),pfion(j,id)
            atot=atot+anion(j,id)
            ctot=ctot+anion(j,id)
         end do
  621 format(i4,a3,3x,1p2e12.4)
  622 format(i4,a3,'+',2x,1p2e12.4)
c
         if(ifmol.gt.0.and.t.le.tmolim) then
             write(6,600)
             do i=1,nmolec
                if(anmol(i,id).gt.ann*1.e-10)
     *             write(6,601) i, cmol(i), anmol(i,id), pfmol(i,id)
                atot=atot+anmol(i,id)
             end do
         end if
 600     format(/ 'Molecular number densities and partition functions'/)
 601     format(i4,1x,A8,1x,1pe12.4,1x,e12.4)
c
            ahmi=1.0353e-16/t/sqrt(t)*exp(8762.9/t)*
     *           anato(1,id)*ane
c
c     original B&C H2+
c
           APLOGJ=amh2(5)
           te=5040./t
           DO K=1,4
              KM5=5-K
              APLOGJ=APLOGJ*TE + amh2(KM5)
           END DO
           tk=1.38054e-16*t
           ph2=-aplogj+log10(anato(1,id)*anion(1,id))+2.*log10(tk)
           anh2b=(10.**ph2)/tk

            htot=anato(1,id)+anion(1,id)+anmol(1,id)+
     *           2.*(anmol(2,id)+anmol(3,id))+anmol(4,id)+anmol(5,id)+
     *           anmol(12,id)+2.*anmol(13,id)+anmol(14,id)+
     *           anmol(15,id)+
     *           anmol(16,id)+anmol(17,id)+anmol(32,id)+anmol(34,id)+
     *           4.*anmol(37,id)+2.*anmol(38,id)+3.*anmol(39,id)+
     *           2.*anmol(40,id)+3.*anmol(41,id)+2.*anmol(57,id)+
     *           anmol(118,id)+anmol(133,id)+
     *           2.*anmol(140,id)+3.*anmol(141,id)+4.*anmol(142,id)+
     *           anmol(148,id)+2.*anmol(149,id)+anmol(222,id)
            ahe=(anato(2,id)+anion(2,id))/htot
            ac= (anato(6,id)+anion(6,id)+anmol(5,id)+anmol(6,id)+
     *           anmol(7,id)+2.*(anmol(8,id)+2.*anmol(13,id))+
     *           anmol(14,id)+2.*anmol(15,id)+anmol(20,id)+
     *           anmol(37,id)+anmol(38,id)+anmol(39,id)+
     *           anmol(44,id)+anmol(118,id)
     *           )/htot
c    *           6.*anmol(310,id))/htot
            an= (anato(7,id)+anion(7,id)+anmol(7,id)+2.*anmol(9,id)+
     *           anmol(11,id)+anmol(12,id)+anmol(14,id)+
     *           anmol(24,id)+anmol(40,id)+anmol(41,id))/htot
            ao= (anato(8,id)+anion(8,id)+anmol(3,id)+anmol(4,id)+
     *           anmol(6,id)+2.*anmol(10,id)+anmol(11,id)+anmol(25,id)+
     *           anmol(26,id)+anmol(29,id)+anmol(30,id)+
     *           anmol(31,id)+anmol(132,id)+anmol(221,id))/htot
            write(6,623) t,dens(id),ann,atot+ane,ane,ctot-anmol(1,id),
     *           anato(1,id),anion(1,id),
     *           anmol(1,id),anmol(2,id),
     *           anmol(312,id),anmol(426,id),anh2b,
     *           htot,
     *           anmol(1,id),ahmi,anmol(1,id)/ahmi,
     *           anato(6,id),anion(6,id),anmol(6,id),anmol(37,id),
     *           anato(7,id),anion(7,id),anmol(9,id),anmol(41,id),
     *           anato(8,id),anion(8,id),anmol(3,id),anmol(6,id),
     *           ahe,ahe/abndd(2,id),
     *           ac,ac/abndd(6,id),
     *           an,an/abndd(7,id),
     *           ao,ao/abndd(8,id)
            act=ac*htot
            ant=an*htot
            aot=ao*htot
  623 format(/'EOS useful quantities - summary'//
     *       'T,rho       ',f13.2,1pe13.5/
     *       'N           ',1p2e13.5/
     *       'n_e         ',1p2e13.5/
     *       'H,H+,H-,H2  ',1p4e13.5/
     *       'H2-,H2+,H2+b',1p3e13.5/
     *       'Htot        ',1pe13.5/
     *       'H-          ',1p3e13.5/
     *       'C,C+,CO,CH4 ',1p4e13.5/
     *       'N,N+,N2,NH3 ',1p4e13.5/
     *       'O,O+,H2O,CO ',1p4e13.5/
     *       'He/H        ',1p2e13.5/
     *       'C/H         ',1p2e13.5/
     *       'N/H         ',1p2e13.5/
     *       'O/H         ',1p2e13.5/)
c
      if(init.eq.1) then
      write(51,625)
      write(52,626)
  625 format('    T      rho     w_mol    Ne/Ntot  N(Htot)    '
     * 'n(H)   n(H2)',6x,
     * 'a(He)   a(C)    a(N)    a(O)     n(C)     n(CO)   n(CH4)',5x,
     * 'n(N)     n(N2)    n(NH3)    n(O)     n(H2O)   n(CO)'/)
       init=0
       end if
c
c      write(51,624) t,dens(id),wmm(id)/hmass,ane/ann,
c    * htot,anato(1,id)/htot,2.*anmol(2,id)/htot,
c    * ahe/abndd(2,id),ac/abndd(6,id),an/abndd(7,id),ao/abndd(8,id),
c    * anato(6,id)/act,anmol(6,id)/act,anmol(37,id)/act,
c    * anato(7,id)/ant,2.*anmol(9,id)/ant,anmol(41,id)/ant,
c    * anato(8,id)/aot,anmol(3,id)/aot,anmol(6,id)/aot
       write(51,624) t,dens(id),wmm(id)/hmass,ane/ann,
     * htot,anato(1,id),2.*anmol(2,id),
     * ahe/abndd(2,id),ac/abndd(6,id),an/abndd(7,id),ao/abndd(8,id),
     * anato(6,id),anmol(6,id),anmol(37,id),
     * anato(7,id),anmol(9,id),anmol(41,id),
     * anato(8,id),anmol(3,id),anmol(6,id)
  624 format(f8.1,1pe9.2,0pf8.5,1x,1p4e9.2,1x,0p4f8.5,1x,1p3e9.2,1x,
     *       3e9.2,1x,3e9.2)
c
       write(52,627) t,dens(id),wmm(id)/hmass,ann,ane,htot,
     * anato(1,id),anion(1,id),anmol(1,id),anmol(2,id),anmol(312,id),
     * anmol(426,id),anh2b
  626 format('    T      rho     w_mol      N        Ne     N(Htot)   ',
     * 'N(H)    N(H+)    N(H-)   N(H2)    N(H2-)    N(H2+)   N(H2+b)'/)
  627 format(f8.1,1pe9.2,0pf8.5,1x,1p10e9.2)
c
      end do

      return
      end
C
C
C *******************************************************************
C
C

      subroutine cia_h2h2(t,ah2,ff,opac)
c     ===================--=============
c
c     CIA H2-H2 opacity
c     data from Borysow A., Jorgensen U.G., Fu Y. 2001, JQSRT 68, 235
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (nlines=1000)
      dimension freq(nlines),temp(7),alpha(nlines,7)
      parameter (amagat=2.6867774d+19,fac=1./amagat**2)
      data temp / 1000. , 2000. , 3000. , 4000. , 5000. , 6000. ,
     *            7000. /
      data ntemp /7/
      data ifirst /0/
      PARAMETER (CAS=2.997925D10)
c     input frequency in Hz but needed wave numbers in cm^-1
      f=ff/cas
c     read in CIA tables if this is the first call
      if (ifirst.eq.0) then
         write(*,'(a)') 'Reading in H2-H2 CIA opacity tables...'
         open(10,file="./data/CIA_H2H2.dat",status='old')
         do i=1,3
            read (10,*)
         enddo
         do i=1,nlines
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp)
         enddo
         close(10)

c     take logarithm of tables prior to doing linear interpolations

         do i=1,nlines
            do j=1,ntemp
               alpha(i,j)=log(alpha(i,j))
            enddo
         enddo

         ifirst=1
      endif

c     locate position in temperature array
      call locate(temp,ntemp,t,j,ntemp)

      if (j.eq.0) then
         write(*,*)
         write(*,'(a,f6.0,a)')
     *   'Warning: requested temperature is below',temp(1),' K'
         write(*,'(a)') 'CIA H2-H2 opacity set to 0'
         write(*,*)
         opac=0.
         return
      endif

c     locate position in frequency array
      call locate(freq,nlines,f,i,nlines)

c     linearly interpolate in frequency and temperature

      if (j.eq.ntemp) then
c     hold values constant if off high temperature end of table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         tt=(f-freq(i))/(freq(i+1)-freq(i))
         alp=(1.-tt)*y1 + tt*y2
      else if (i.eq.0 .or. i.eq.nlines) then
c     set values to a very small number if off frequency table
         alp=-50.
      else
c     interpolate linearly within table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         y3=alpha(i+1,j+1)
         y4=alpha(i,j+1)

         tt=(f-freq(i))/(freq(i+1)-freq(i))
         uu=(t-temp(j))/(temp(j+1)-temp(j))

         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +
     *       (1.-tt)*uu*y4
      endif

      alp=exp(alp)

c     final opacity

      opac=fac*ah2*ah2*alp
c
      return
      end

C
C
C
C ********************************************************************
C
C

      SUBROUTINE locate(xx,n,x,j,nxdim)
c     =================================
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension xx(nxdim)
c
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      if(x.eq.xx(1)) then
        j=1
      else if(x.eq.xx(n)) then
        j=n-1
      else
        j=jl
      endif
      return
      END

C
C
C ********************************************************************
C
C

      subroutine cia_h2he(t,ah2,ahe,ff,opac)
c     ======================================
c
c     CIA H2-He opacity
c     data from Jorgensen U.G., Hammer D., Borysow A., Falkesgaard J., 2000,
c     Astronomy & Astrophysics 361, 283
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (nlines=242)
      dimension freq(nlines),temp(7),alpha(nlines,7)
      parameter (amagat=2.6867774d+19,fac=1./amagat**2)
      data temp / 1000. , 2000. , 3000. , 4000. , 5000. , 6000. ,
     *            7000. /
      data ntemp /7/
      data ifirst /0/
      PARAMETER (CAS=2.997925D10)
c     input frequency in Hz but needed wave numbers in cm^-1
      f=ff/cas
c     read in CIA tables if this is the first call
      if (ifirst.eq.0) then
         write(*,'(a)') 'Reading in H2-He CIA opacity tables...'
         open(10,file="./data/CIA_H2He.dat",status='old')
         do i=1,3
            read (10,*)
         enddo
         do i=1,nlines
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp)
         enddo
         close(10)

c     take logarithm of tables prior to doing linear interpolations

         do i=1,nlines
            do j=1,ntemp
               alpha(i,j)=log(alpha(i,j))
            enddo
         enddo

         ifirst=1
      endif

c     locate position in temperature array
      call locate(temp,ntemp,t,j,ntemp)

      if (j.eq.0) then
         write(*,*)
         write(*,'(a,f6.0,a)')
     *   'Warning: requested temperature is below',temp(1),' K'
         write(*,'(a)') 'CIA H2-He opacity set to 0'
         write(*,*)
         opac=0.
         return
      endif

c     locate position in frequency array
      call locate(freq,nlines,f,i,nlines)

c     linearly interpolate in frequency and temperature

      if (j.eq.ntemp) then
c     hold values constant if off high temperature end of table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         tt=(f-freq(i))/(freq(i+1)-freq(i))
         alp=(1.-tt)*y1 + tt*y2
      else if (i.eq.0 .or. i.eq.nlines) then
c     set values to a very small number if off frequency table
         alp=-50.
      else
c     interpolate linearly within table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         y3=alpha(i+1,j+1)
         y4=alpha(i,j+1)

         tt=(f-freq(i))/(freq(i+1)-freq(i))
         uu=(t-temp(j))/(temp(j+1)-temp(j))

         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +
     *       (1.-tt)*uu*y4
      endif

      alp=exp(alp)

c     final opacity

      opac=fac*ah2*ahe*alp
c
      return
      end
C
C
C *******************************************************************
C
C

      subroutine cia_h2h(t,ah2,ah,ff,opac)
c     ====================================
c
c     CIA H2-H opacity - data taken from TURBOSPEC
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (nlines=67)
      dimension freq(nlines),temp(4),alpha(nlines,4)
      parameter (amagat=2.6867774d+19,fac=1./amagat**2)
      data temp / 1000. , 1500., 2000. , 2500. /
      data ntemp /4/
      data ifirst /0/
      PARAMETER (CAS=2.997925D10)
c     input frequency in Hz but needed wave numbers in cm^-1
      f=ff/cas
c     read in CIA tables if this is the first call
      if (ifirst.eq.0) then
         write(*,'(a)') 'Reading in H2-H CIA opacity tables...'
         open(10,file="./data/CIA_H2H.dat",status='old')
         do i=1,3
            read (10,*)
         enddo
         do i=1,nlines
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp)
         enddo
         close(10)

c     take logarithm of tables prior to doing linear interpolations

         do i=1,nlines
            do j=1,ntemp
               alpha(i,j)=log(alpha(i,j))
            enddo
         enddo

         ifirst=1
      endif

c     locate position in temperature array
      call locate(temp,ntemp,t,j,ntemp)

      if (j.eq.0) then
         write(*,*)
         write(*,'(a,f6.0,a)')
     *   'Warning: requested temperature is below',temp(1),' K'
         write(*,'(a)') 'CIA H2-H opacity set to 0'
         write(*,*)
         opac=0.
         return
      endif

c     locate position in frequency array
      call locate(freq,nlines,f,i,nlines)

c     linearly interpolate in frequency and temperature

      if (j.eq.ntemp) then
c     hold values constant if off high temperature end of table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         tt=(f-freq(i))/(freq(i+1)-freq(i))
         alp=(1.-tt)*y1 + tt*y2
      else if (i.eq.0 .or. i.eq.nlines) then
c     set values to a very small number if off frequency table
         alp=-50.
      else
c     interpolate linearly within table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         y3=alpha(i+1,j+1)
         y4=alpha(i,j+1)

         tt=(f-freq(i))/(freq(i+1)-freq(i))
         uu=(t-temp(j))/(temp(j+1)-temp(j))

         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +
     *       (1.-tt)*uu*y4
      endif

      alp=exp(alp)

c     final opacity

      opac=fac*ah2*ah*alp
c
      return
      end
C
C
C *******************************************************************
C
C

      subroutine cia_hhe(t,ah,ahe,ff,opac)
c     ====================================
c
c     CIA H-He opacity
c     data from Gustafsson M., Frommhold, L. 2001, ApJ 546, 1168
c
      IMPLICIT REAL*8(A-H,O-Z)
      parameter (nlines=43)
      dimension freq(nlines),temp(11),alpha(nlines,11)
      parameter (amagat=2.6867774d+19,fac=1./amagat**2)
      data temp / 1000.,  1500.,  2250., 3000.,  4000.,  5000.,
     *            6000.,  7000., 8000.,  9000., 10000./
      data ntemp /11/
      data ifirst /0/
      PARAMETER (CAS=2.997925D10)
c     input frequency in Hz but needed wave numbers in cm^-1
      f=ff/cas
c     read in CIA tables if this is the first call
      if (ifirst.eq.0) then
         write(*,'(a)') 'Reading in H-He CIA opacity tables...'
         open(10,file="./data/CIA_HHe.dat",status='old')
         do i=1,3
            read (10,*)
         enddo
         do i=1,nlines
            read (10,*) freq(i),(alpha(i,j),j=1,ntemp)
         enddo
         close(10)

c     take logarithm of tables prior to doing linear interpolations

         do i=1,nlines
            do j=1,ntemp
               alpha(i,j)=log(alpha(i,j))
            enddo
         enddo

         ifirst=1
      endif

c     locate position in temperature array
      call locate(temp,ntemp,t,j,ntemp)

      if (j.eq.0) then
         write(*,*)
         write(*,'(a,f6.0,a)')
     *   'Warning: requested temperature is below',temp(1),' K'
         write(*,'(a)') 'CIA H-He opacity set to 0'
         write(*,*)
         opac=0.
         return
      endif

c     locate position in frequency array
      call locate(freq,nlines,f,i,nlines)

c     linearly interpolate in frequency and temperature

      if (j.eq.ntemp) then
c     hold values constant if off high temperature end of table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         tt=(f-freq(i))/(freq(i+1)-freq(i))
         alp=(1.-tt)*y1 + tt*y2
      else if (i.eq.0 .or. i.eq.nlines) then
c     set values to a very small number if off frequency table
         alp=-50.
      else
c     interpolate linearly within table
         y1=alpha(i,j)
         y2=alpha(i+1,j)
         y3=alpha(i+1,j+1)
         y4=alpha(i,j+1)

         tt=(f-freq(i))/(freq(i+1)-freq(i))
         uu=(t-temp(j))/(temp(j+1)-temp(j))

         alp=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +
     *       (1.-tt)*uu*y4
      endif

      alp=exp(alp)

c     final opacity

      opac=fac*ah*ahe*alp
c
      return
      end
C
C
C *******************************************************************
C
C
      subroutine h2minus(t,anh2,ane,fr,oph2m)
C     =======================================
C
C     H- free-free opacity
C
C     data from K L Bell 1980 J. Phys. B: At. Mol. Phys. 13 1859, Table 1
C     The first column is theta=5040/T(K)
C     The first row are names for each row corresponding to lambda (angstroms)
C     The last row for 10.0 is linearly extrapolated
C     The units of everything else is 10^26 cm4/dyn-1
C 
      INCLUDE 'PARAMS.FOR'
      dimension FFthet(9),FFlamb(18),FFkapp(18,9)
      data FFthet / 0.5, 0.8, 1.0, 1.2, 1.6, 2.0,
     *     2.8, 3.6, 10.0 /
      data nthet /9/
      data FFlamb /151883., 113913.,  91130.,  60753.,
     *     45565.,  36452.,  30377.,  22783.,
     *     18226.,  15188.,  11391.,  9113.,  7594.,
     *     6509.,  5696.,  5063.,  4142.,  3505./
      data nlamb /18/
      data FFkapp /
     *     7.16e+01,4.03e+01,2.58e+01,1.15e+01,6.47e+00,
     *     4.15e+00,2.89e+00,1.63e+00,1.05e+00,7.36e-01,
     *     4.20e-01,2.73e-01,1.92e-01,1.43e-01,1.10e-01,
     *     8.70e-02,5.84e-02,4.17e-02,9.23e+01,5.20e+01,
     *     3.33e+01,1.48e+01,8.37e+00,5.38e+00,3.76e+00,
     *     2.14e+00,1.39e+00,9.75e-01,5.64e-01,3.71e-01,
     *     2.64e-01,1.98e-01,1.54e-01,1.24e-01,8.43e-02,
     *     6.10e-02,1.01e+02,5.70e+01,3.65e+01,1.63e+01,
     *     9.20e+00,5.92e+00,4.14e+00,2.36e+00,1.54e+00,
     *     1.09e+00,6.35e-01,4.22e-01,3.03e-01,2.30e-01,
     *     1.80e-01,1.46e-01,1.01e-01,7.34e-02,1.08e+02,
     *     6.08e+01,3.90e+01,1.74e+01,9.84e+00,6.35e+00,
     *     4.44e+00,2.55e+00,1.66e+00,1.18e+00,6.97e-01,
     *     4.67e-01,3.39e-01,2.59e-01,2.06e-01,1.67e-01,
     *     1.17e-01,8.59e-02,1.18e+02,6.65e+01,4.27e+01,
     *     1.91e+01,1.08e+01,6.99e+00,4.91e+00,2.84e+00,
     *     1.87e+00,1.34e+00,8.06e-01,5.52e-01,4.08e-01,
     *     3.17e-01,2.55e-01,2.10e-01,1.49e-01,1.11e-01,
     *     1.26e+02,7.08e+01,4.54e+01,2.04e+01,1.16e+01,
     *     7.50e+00,5.28e+00,3.07e+00,2.04e+00,1.48e+00,
     *     9.09e-01,6.33e-01,4.76e-01,3.75e-01,3.05e-01,
     *     2.53e-01,1.82e-01,1.37e-01,1.38e+02,7.76e+01,
     *     4.98e+01,2.24e+01,1.28e+01,8.32e+00,5.90e+00,
     *     3.49e+00,2.36e+00,1.74e+00,1.11e+00,7.97e-01,
     *     6.13e-01,4.92e-01,4.06e-01,3.39e-01,2.49e-01,
     *     1.87e-01,1.47e+02,8.30e+01,5.33e+01,2.40e+01,
     *     1.38e+01,9.02e+00,6.44e+00,3.90e+00,2.68e+00,
     *     2.01e+00,1.32e+00,9.63e-01,7.51e-01,6.09e-01,
     *     5.07e-01,4.27e-01,3.16e-01,2.40e-01,2.19e+02,
     *     1.26e+02,8.13e+01,3.68e+01,2.18e+01,1.46e+01,
     *     1.08e+01,7.18e+00,5.24e+00,4.17e+00,3.00e+00,
     *     2.29e+00,1.86e+00,1.55e+00,1.32e+00,1.13e+00,
     *     8.52e-01,6.64e-01/

c     locate position in temperature array
      theta=5040./t

      call locate(FFthet,nthet,theta,j,nthet)
      if (j.eq.0) then
 1       write(*,*)
         write(*,'(a,f6.0,a)')
     *   'Error: requested temperature is outside the ranges'
         write(*,'(a)') 'h2minus:Stop'
         write(*,*)
         stop
      endif
      flamb=CL*1.D8/fr
c     locate position in wavelength array
      call locate(FFlamb,nlamb,flamb,i,nlamb)

c     linearly interpolate in frequency and temperature
      if (j.eq.nthet) then
c     hold values constant if off high temperature end of table
         y1=FFkapp(i,j)
         y2=FFkapp(i+1,j)
         tt=(flamb-FFlamb(i))/(FFlamb(i+1)-FFlamb(i))
         Fkappa=(1.-tt)*y1 + tt*y2
      else if (i.eq.0 .or. i.eq.nlines) then
c     set values to 0 if off frequency table
         Fkappa=0.0
      else
c     interpolate linearly within table
         y1=FFkapp(i,j)
         y2=FFkapp(i+1,j)
         y3=FFkapp(i+1,j+1)
         y4=FFkapp(i,j+1)

         tt=(flamb-FFlamb(i))/(FFlamb(i+1)-FFlamb(i))
         uu=(theta-FFthet(j))/(FFthet(j+1)-FFthet(j))

         Fkappa=(1.-tt)*(1.-uu)*y1 + tt*(1.-uu)*y2 + tt*uu*y3 +
     *       (1.-tt)*uu*y4
      endif
      pe=ane*BOLK*t
      oph2m= anh2 * 1.0E-26 *pe * Fkappa
      return
      end
c
c
c   **********************************************************************
c
c
      subroutine h2opf(t,pf)
c   
c     partition function for H2Ofrom EXOMOILA data
c
      INCLUDE 'PARAMS.FOR'
      dimension ttab(10000),pftab(10000)
c   
      data init /1/
c
      if(init.eq.1) then
         open(67,file='./data/h2o_exomol.pf',status='old')
         do i=1,10000
            read(67,*) ttab(i),pftab(i)
         end do 
         close(67)
         init=0
      end if
c   
      itab=ifix(real(t))
      pf=pftab(itab)+(t-ttab(itab))*(pftab(itab+1)-pftab(itab))
      return
      end

C
C
C *******************************************************************
C
C
      function gvdw(il,ilist,id)
c     ==========================
c
c     evaluation of the Van der Waals broadening parameter
c
c     currently, two possibilities, determined by the value of the parameter
c     ivdwli(ilist) - the mode of evaluation is the same for the whole line list
c       = 0 - standard expression
c       > 0 - evaluation using EXOMOL data, assuming breadening by H2 and He
c
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'LINDAT.FOR'
      COMMON/PRFQUA/DOPA1(MATOM,MDEPTH),VDWC(MDEPTH)
c
c     clasical, original expression
c
      if(ivdwli(ilist).eq.0) then
         gvdw=gwm(il,ilist)*vdwc(id)
         return
      end if
c
c     EXOMOL form - broadening by H2 and He
c
c     con= 1.e-6*c*k
      con=4.1388e-12
      t=temp(id)
      anhe=rrr(id,1,2)
      gvdw=con*t*((296./t)**gexph2(il,ilist)*gvdwh2(il,ilist)*anh2(id)+
     *            (296./t)**gexphe(il,ilist)*gvdwhe(il,ilist)*anhe)
      return
      end
C
C
C *******************************************************************
C
C


