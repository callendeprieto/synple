#
#	makefile for downloading synple/synspec line lists
#	Carlos 2019
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/atlases/
WGET = `which wget`

#general line lists
STANDARD =  ardata.fits iag-nir.fits  iag-vis.fits  solarfluxintwl.fits  ardata-ReadMe.txt  iag-ReadMe.txt   solarfluxintwl-ReadMe.txt
#STANDARD = solarfluxintwl.fits solarfluxintwl-ReadMe.txt


all: 
	$(call sequence, $(STANDARD))

clean: 
	rm -f *.fits *.txt


define sequence
	for file in $(1) ; do \
		${WGET} ${FTP}$$file ; \
	done
endef

