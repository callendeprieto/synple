      program list2bin
c     ================
c
c     transforms an ASCII line list to a binary
c
c     INPUT:  std input - ASCII list
c
c     OUTPUT: fort.12 - binary list
c             std. output gives an info on number of lines
c
      implicit real*8 (a-h,o-z)
      character*102 rec
      dimension x(10),j(7)
c
c     first, find the format of the line list by examining the 1st record
c
      n=1
      npr=10
      npi=7
      read(5,'(a102)') rec
      write(*,*) rec
      read(rec,*,iostat=kst) (x(i),i=1,10),(j(i),i=1,7)
      if(kst.ne.0) then
         read(rec,*,iostat=kst1) (x(i),i=1,10),j(1)
         npr=10
         npi=1
         if(kst1.ne.0) then
            read(rec,*,iostat=kst2) (x(i),i=1,9)
            npr=9
            npi=0
            if(kst2.ne.0) then
               read(rec,*,iostat=kst3) (x(i),i=1,7)
               npr=7
               npi=0
               if(kst3.ne.0) then
                  read(rec,*,iostat=kst4) (x(i),i=1,4)
                  npr=4
                  npi=0
                  if(kst4.ne.0) then
                     write(*,*) 
     *               'no applicable format of the line list detected'
                     stop
                  end if
               end if
            end if
         end if
      end if
c
      write(6,600) npr,npi
  600 format('number of parameters in the list; real,int:',2i4)
      write(12) (x(i),i=1,npr),(j(i),i=1,npi)
c
c     read the list and transport to binary
c
   10 continue
      read(5,*,err=10,end=20) (x(i),i=1,npr),(j(i),i=1,npi)
      n=n+1
      write(12) (x(i),i=1,npr),(j(i),i=1,npi)
      go to 10
   20 continue
c
      write(6,601) n
  601 format(i10,'  lines included in the binary file fort.12')
      end

