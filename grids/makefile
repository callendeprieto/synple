#
#	makefile for downloading bas/ferre grids 
#	Carlos 2024
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/grids/
WGET = `which wget`

STANDARD = n_sc2-BRZ-DESIrbf.pickle n_sc2-BRZ-DESIrbf-EBV.pickle n_sc2-IDS-R900Vrbf.pickle n_sc2-LAMOSTrbf.pickle n_sc2-LAMOSTrbf-EBV.pickle n_sc2-MILESrbf.pickle n_sc2-OSIRIS-R2500Urbf.pickle n_sc2-STISrbf.pickle n_sc2-STISrbf-EBV.pickle n_sc2-STIS-NICMOSrbf.pickle n_sc2-STIS-NICMOSrbf-EBV.pickle n_sc2-XPrbf.pickle n_sc2-XPrbf-EBV.pickle


all: 
	$(call sequence, $(STANDARD))

clean: 
	rm -f *.pickle


define sequence
	for file in $(1) ; do \
		${WGET} ${FTP}$$file ; \
	done
endef

