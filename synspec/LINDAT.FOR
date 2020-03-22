c     PARAMETER (MLIN0   = 200000, 
      PARAMETER (MLIN0   =1000000, 
c     PARAMETER (MLIN0   = 500000, 
     *           MGRIEM  =     10,
     *           MNLT    =   2000,
     *           MSPHE2  =     20,
c    *           MLIN    =  10000,
     *           MLIN    = 190000,
     *           MPRF    =  MLIN0)
C
c     PARAMETER (MLINM0  = 300000,
c      PARAMETER (MLINM0  =1000000,
      PARAMETER (MLINM0  =25000000,
     *           MLINM   = 500000,
     *           MMLIST  =      4)
c     PARAMETER (MLINM0  = 20,
c    *           MLINM   =  2,
c    *           MMLIST  =      2)
C
       REAL*4       EXCL0(MLIN0),
     *              EXCU0(MLIN0),
     *              GF0(MLIN0),
     *              EXTIN(MLIN0),
     *              BNUL(MLIN0),
     *              GAMR0(MPRF),
     *              GS0(MPRF),
     *              GW0(MPRF),
     *              WGR0(4,MGRIEM),
     *              EXCLM(MLINM0,MMLIST),
     *              GFM(MLINM0,MMLIST),
     *              EXTINM(MLINM0,MMLIST),
     *              GRM(MLINM0,MMLIST),
     *              GSM(MLINM0,MMLIST),
     *              GWM(MLINM0,MMLIST)
C     
      COMMON/LINTOT/FREQ0(MLIN0),
     *              EXCL0,
     *              EXCU0,
     *              GF0,
     *              EXTIN,
     *              BNUL,
     *              INDAT(MLIN0),
     *              INDNLT(MLIN0),
     *              ILOWN(MLIN0),
     *              IUPN(MLIN0),
     *              IJCONT(MLIN0),
     *              INDLIN(MLIN),
     *              INDLIP(MLIN),
     *              NLIN0,NLIN,IRLIST,
     *              NNLT,NGRIEM
C
      COMMON/MOLTOT/FREQM(MLINM0,MMLIST),
     *              EXCLM,
     *              GFM,
     *              EXTINM,
     *              GRM,
     *              GSM,
     *              GWM,
     *              INDATM(MLINM0,MMLIST),
     *              INMLIN(MLINM,MMLIST),
     *              INMLIP(MLINM,MMLIST),
     *              NLINM0(MMLIST),
     *              NLINML(MMLIST),
     *              IUNITM(MMLIST),
     *              NMLIST
C
      COMMON/LINPRF/GAMR0,
     *              GS0,
     *              GW0,
     *              WGR0,
     *              IPRF0(MPRF),
     *              ISPRF(MPRF),
     *              IGRIEM(MPRF),
     *              ISP0(MSPHE2),NSP
C
      COMMON/LINNLT/ABCENT(MNLT,MDEPTH),
     *              SLIN(MNLT,MDEPTH)
C
      COMMON/LINDEP/PLAN(MDEPTH),
     *              STIM(MDEPTH),
     *              EXHK(MDEPTH)
C
      COMMON/LINCTR/DFRCON,IJCNTR(MLIN),IJCMTR(MLINM,MMLIST)
      COMMON/MLINRE/FRLASM(MMLIST),ALASTM(MMLIST),
     *              NXTSEM(MMLIST),IPRSEM(MMLIST),IREADM(MMLIST)

