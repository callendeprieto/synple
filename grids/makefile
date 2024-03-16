#
#	makefile for downloading synple/bas/FERRE grids 
#	Carlos 2024
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/grids/
WGET = `which wget`

#general line lists
STANDARD =  n_sc1-BRZdesirbf.pickle n_sc1-lamostrbf.pickle n_sc1-stisrbf.pickle n_sc1-xprbf.pickle 


all: 
	$(call sequence, $(STANDARD))


clean: 
	rm -f *.pickle

define sequence
	for file in $(1) ; do \
		${WGET} ${FTP}$$file ; \
	done
endef

