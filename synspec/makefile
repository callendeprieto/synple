#
#	makefile for synspec/rotin
#	Carlos 2019
#

F90 = gfortran 
OPT= -fno-automatic -mcmodel=medium

#F90 = g95
#OPT= -fstatic -mcmodel=medium -Wl,--no-relax

#F90 = sunf95
#OPT = -Bstatic -xmodel=medium 
 
#F90 =  ifort 
#OPT = -shared-intel -mcmodel medium 

OS1 = ${shell perl -e 'my $$osVersion = `uname `; chomp($$osVersion); print "$$osVersion"'}
OS2 = ${shell perl -e 'my $$osVersion = `uname -r`; chomp($$osVersion); print "$$osVersion"'}
OS3 = ${shell perl -e 'my $$osVersion = `machine`; chomp($$osVersion); print "$$osVersion"'}
OS4 = ${shell perl -e 'my $$osVersion = `hostname`; chomp($$osVersion); print "$$osVersion"'}


all:  synspec rotin install 

synspec: s54d.f  LINDAT.FOR  MODELP.FOR  PARAMS.FOR  SYNTHP.FOR  WINCOM.FOR
	${F90}  ${OPT} s54d.f -o s54d

rotin: rotin3.f
	${F90}  ${OPT} rotin3.f -o rotin3

install:
	mv s54d   ../bin/
	mv rotin3 ../bin
                                                        
clean: 
	rm -f s54d rotin3 ../bin/s54d ../bin/rotin3 


