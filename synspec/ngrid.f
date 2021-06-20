      program ngrid
C     =============
c
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER (MFGRID  =  300000,
     *           MTTAB   =      25,
     *           MRTAB   =      15,
     *           MSFTAB  = 6000000)
c
      real*4 absgrd(mttab,mrtab,mfgrid)
      common/gridp0/tempg(mttab),densg(mrtab),elecgr(mttab,mrtab),
     *              densg0(mttab),temp1,ntemp,ndens
      common/gridf0/wlgrid(mfgrid),nfgrid
      common/fintab/absgrd
      common/initab/absop(msftab),wltab(msftab)
      common/elecm0/elecm(mttab)
      dimension tempvec(mttab),rhovec(mrtab)
      dimension abgrd(mfgrid),abunt(30),abuno(30)
      dimension typa(30)
      character*(80) tabname,optable
      character*4 typa
c
      istop=0
      read(5,*) nlamb,inttab,wlam1,wlam2
      read(5,*) tabname,ibingr
      if(nlamb.gt.mfgrid) then
         write(*,*) 'nlambda.gt.mfgrid - recompile with larger mfgrid'
         istop=1
      end if
      ibinop=0
      read(5,*,err=5,end=5) optable,ibinop
    5 continue
      write(*,*) 'original optable  ',optable
c
      nfgrid=nlamb
      wl1=log(wlam1)
      wl2=log(wlam2)
      dwl=(wl2-wl1)/(nfgrid-1)
      do i=1,nfgrid
         wlgrid(i)=exp(wl1+(i-1)*dwl)
      end do
c
c     read the header of the old opacity table
c
      if(ibinop.eq.0) then
         open(52,file=optable,status='old')
       else
         open(52,file=optable,form='unformatted',status='old')
      end if
c
      if(ibinop.eq.0) then
         read(52,*)
         read(52,*)
         do iat=1,30
            read(52,*) typa(iat),abunt(iat),abuno(iat)
         end do
         read(52,*)
         read(52,*)
         read(52,*) ifmolt,tmolit
         read(52,*)
         read(52,*)
         read(52,*) iophmt,ioph2t,iophet,iopcht,iopoht
         read(52,*)
         read(52,*)
         read(52,*) numfre0,numtem0,numrh0
         read(52,*)
         read(52,*) (tempvec(i),i=1,numtem0)
         read(52,*)
         read(52,*) (rhovec(j),j=1,numrh0)
         read(52,*)
         read(52,*) ((elecgr(i,j),j=1,numrh0),i=1,numtem0)
       else
         do iat=1,30
            read(52) typa(iat),abunt(iat),abuno(iat)
         end do
         read(52) ifmolt,tmolit
         read(52) iophmt,ioph2t,iophet,iopcht,iopoht
         read(52) numfre0,numtem0,numrh0
         read(52) (tempvec(i),i=1,numtem0)
         read(52) (rhovec(j),j=1,numrh0)
         read(52) ((elecgr(i,j),j=1,numrh0),i=1,numtem0)
      end if
      ntemp=numtem0
      ndens=numrh0
      write(*,*) 'ntemp,ndens',ntemp,ndens
      if(ntemp.gt.mttab) then
         write(*,*) 'ntemp.gt.mttab - recompile with larger mttab'
         istop=1
      end if
      if(ndens.gt.mrtab) then
         write(*,*) 'ndens.gt.mrtab - recompile with larger mrtab'
         istop=1
      end if
c
      if(istop.eq.1) stop
c
c     read file fort.27 (internal opacitis) and set the
c     arrays for the new opacity tables
c
      if(inttab.ne.1) then
c
c     1st possibility - approximate preserving an integral of opacity
c
         indext=0
         indexn=0
         nfr=0
   10    continue
         indext=indext+1
   20    continue
         indexn=indexn+1
         nfr=0
         ijgrd=0
   30    continue
         ijgrd=ijgrd+1
         wlgr=0.5*(wlgrid(ijgrd)+wlgrid(ijgrd+1))
         isum=0
         sum=0.
   40    continue
         read(27,*,err=50,end=50) ip,wl,abl
         if(wl.le.wlgr) then
            sum=sum+exp(abl)
            isum=isum+1
            if(ip.le.nfr) go to 50
            nfr=ip
            go to 40
         end if
         abgrd(ijgrd)=log(sum/float(isum))
         if(ijgrd.lt.nfgrid) go to 30
   50    continue
         write(*,*) 'it,ir,nf',indext,indexn,nfr,nfgrid
         do ij=1,nfgrid
            absgrd(indext,indexn,ij)=real(abgrd(ij))
         end do
      absgrd(indext,indexn,nfgrid)=absgrd(indext,indexn,nfgrid-1)

c
         if(indexn.lt.ndens) then
            go to 20
          else
            indexn=0
         end if
         if(indext.lt.ntemp) then
            go to 10
          else
            indext=0
         end if
c
       else
c
c        2nd possibility - an interpolation in wavelengths
c
         indext=0
         indexn=0
         nfr=0
  110    continue
         indext=indext+1
  120    continue
         indexn=indexn+1
         nfr=0
  130    continue
         read(27,*,err=150,end=150) ip,wl,abl
         if(ip.le.nfr) go to 150
         wltab(ip)=wl
         absop(ip)=abl
         nfr=ip
         go to 130
  150    continue
         write(*,*) 'it,ir,nf',indext,indexn,nfr,nfgrid
         call intrp(wltab,absop,wlgrid,abgrd,nfr,nfgrid)
         do ij=1,nfgrid
            absgrd(indext,indexn,ij)=real(abgrd(ij))
         end do
c
         if(indexn.lt.ndens) then
            go to 120
          else
            indexn=0
         end if
         if(indext.lt.ntemp) then
            go to 110
          else
            indext=0
         end if
      end if
c
c     store the opacities in the new table
c
      if(ibingr.eq.0) then
         open(53,file=tabname,status='unknown')
         write(53,600)
         do iat=1,30
            write(53,601) typa(iat),abunt(iat),abuno(iat)
         end do
         write(53,602) ifmolt,tmolit
         write(53,603) iophmt,ioph2t,iophet,iopcht,iopoht
         write(53,611) nfgrid,ntemp,ndens
         write(53,612) (tempvec(i),i=1,ntemp)
         write(53,613) (rhovec(j),j=1,ndens)
         write(53,614) ((elecgr(i,j),j=1,ndens),i=1,ntemp)
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
  603    format('additional opacities'/'  H-  H2+ He-  CH  OH'/
     *          5i4)
  611    format(/'number of frequencies, temperatures, densities:'
     *          /10x,3i10)
  612    format('log temperatures'/(6F11.6))
  613    format('log densities'/(6F11.6))
  614    format('log electron densities from EOS'/(6f11.6))
  615    format(/' *** frequency # : ',i8,f15.5/1pe20.8)
  616    format((1p6e14.6))
       else
         do iat=1,30
            write(63) typa(iat),abunt(iat),abuno(iat)
         end do
         write(63) ifmolt,tmolit
         write(63) iophmt,ioph2t,iophet,iopcht,iopoht
         write(63) nfgrid,ntemp,ndens
         write(63) (tempvec(i),i=1,ntemp)
         write(63) (rhovec(j),j=1,ndens)
         write(63) ((elecgr(i,j),j=1,ndens),i=1,ntemp)
         do k = 1, nfgrid
            write(63) wlgrid(k)
            do j = 1, ndens
               write(63) (absgrd(i,j,k),i=1,ntemp)
            end do
         end do
      end if
c
      end
C
C
C     ****************************************************************
C
C
      subroutine intrp(wltab,absop,wlgrid,abgrd,nfr,nfgrid)
c     =====================================================
c
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER (MFGRID  =  400000)
c
      dimension wltab(1),absop(1),wlgrid(1),abgrd(1)
      dimension yint(mfgrid),jint(mfgrid)
c
c      set up interpolation coefficients for frequency interpolation
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
c
      return
      end

