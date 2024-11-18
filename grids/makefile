#
#	makefile for downloading synple/bas/FERRE grids 
#	Carlos 2024
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/grids/
WGET = `which wget`

#general line lists
STANDARD =  n_sc2-BRZ-DESIrbf.pickle n_sc2-LAMOSTrbf.pickle n_sc2-STISrbf.pickle n_sc2-XPrbf.pickle  n_sc2-STIS-NICMOSrbf.pickle  n_sc2-OSIRIS-R2500Urbf.pickle  n_sc2-IDS-R900Vrbf.pickle n_DA-STISrbf.pickle



all: 
	$(call sequence, $(STANDARD))


clean: 
	rm -f *.pickle

define sequence
	for file in $(1) ; do \
		${WGET} ${FTP}$$file ; \
	done
endef



