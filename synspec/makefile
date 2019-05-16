#
#	makefile for synspec/rotin
#	Carlos 2019
#

F90 = gfortran 
OPT= -fno-automatic -mcmodel=medium  

#F90 = g95
#OPT= -fno-automatic -mcmodel=medium

#F90 = sunf95
#OPT = -Bstatic -xmodel=medium 
 
#F90 =  ifort 
#OPT = -shared-intel -mcmodel medium 

OS1 = ${shell perl -e 'my $$osVersion = `uname `; chomp($$osVersion); print "$$osVersion"'}
OS2 = ${shell perl -e 'my $$osVersion = `uname -r`; chomp($$osVersion); print "$$osVersion"'}
OS3 = ${shell perl -e 'my $$osVersion = `machine`; chomp($$osVersion); print "$$osVersion"'}
OS4 = ${shell perl -e 'my $$osVersion = `hostname`; chomp($$osVersion); print "$$osVersion"'}


all:  synspec rotin install 

synspec: synspec53p.f  LINDAT.FOR  MODELP.FOR  PARAMS.FOR  SYNTHP.FOR  WINCOM.FOR
	${F90}  ${OPT} synspec53p.f -o synspec53p

rotin: rotin3.f
	${F90}  ${OPT} rotin3.f -o rotin3

install:
	mv synspec53p   ../bin/
	mv rotin3 ../bin
                                                        
clean: 
	rm -f synspec53p rotin3


