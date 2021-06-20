#
#	makefile for synspec/rotin
#	Carlos 2019
#

F90 = gfortran 
OPT= -fno-automatic 

#F90 = g95
#OPT= -fstatic -Wl,--no-relax

#F90 = sunf95
#OPT = -Bstatic -xmodel=medium 
 
#F90 =  ifort 
#OPT = -shared-intel -mcmodel medium 

#F90 = nvfortran
#OPT = -Bstatic -mcmodel=medium

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    OPT += -mcmodel=medium 
endif


all:  compile install 

compile: synspec54 rotin list2bin

synspec54: synspec54.f  LINDAT.FOR  MODELP.FOR  PARAMS.FOR  SYNTHP.FOR  WINCOM.FOR
	${F90}  ${OPT} synspec54.f -o synspec54

rotin: rotin.f
	${F90}  ${OPT} rotin.f -o rotin

list2bin: list2bin.f
	${F90} ${OPT} list2bin.f -o list2bin

install: ../bin ../bin/synspec54 ../bin/rotin ../bin/list2bin

../bin/synspec54 ../bin/rotin ../bin/list2bin: synspec54 rotin list2bin
	cp synspec54   ../bin/
	cp rotin ../bin/
	cp list2bin ../bin/

../bin:
	mkdir ../bin
                                                        
clean: 
	rm -f synspec54 rotin list2bin ../bin/synspec54 ../bin/rotin ../bin/list2bin


