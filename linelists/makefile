#
#	makefile for downloading synple/synspec line lists
#	Carlos 2019
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/linelists/
WGET = `which wget`
GUNZIP = `which gunzip`

#general line lists
#STANDARD =  gfall08oct17hix.19 kmol2019.20  #original
#STANDARD =  gfATO.19 gfMOL.20		     #gfMOL.20 reduced with redulis/gfelim>-8)
STANDARD =  gfATO.19 gfMOLsun.20 gfTiO.20 H2O-8.20 #TiO and H2O from EXOMOL(gfelim>-8)

#isotopic ratios more typical of a red giant (after 1st dredge up)
GIANT = gfMOLgiant.20 # to replace gfMOLsun (also reduced with redulis/gfelim>-8)

#Allende Prieto et al. 2018
OLD = gfallx3_bpo.19 kmol3_0.01_30.20

#Add TiO and H2O (older Schwenke, Partridge & Schwenke line lists)
SCHWENKE = tio.20 h2o.20

#apogee
HBAND = synspec.20200921nlte.19 synspec.20200921nltesun.20 synspec.20200921nltegiant.20 synspec.h2o-BC8.5V.20

#
#Yeisson-synple                     SDSS-SAS
#
#apogee.20180901.giants.20  synspec.20180901giant.molec
#apogee.20180901.solar.20   synspec.20180901sun.molec
#
#synspec20200709.lte.19     synspec.20200709.atoms
#synspec20200705.nlte.19    synspec.20200705nlte.atoms 


all: 
	$(call sequence, $(STANDARD))


old:
	$(call sequence, $(OLD))

schwenke:
	$(call sequence, $(SCHWENKE))

hband: 
	$(call sequence, $(HBAND))


clean: 
	rm -f *.19 *.20 *.11 *.00 *.gz


define sequence
	for file in $(1) ; do \
		${WGET} ${FTP}$$file.gz ; \
	       	${GUNZIP} -c $$file.gz > /tmp/$$file ; \
		wc /tmp/$$file       > $$file.00  ; \
		head -1 /tmp/$$file >> $$file.00  ; \
		tail -1 /tmp/$$file >> $$file.00  ; \
		rm -f $$file.gz ; \
		../bin/list2bin < /tmp/$$file ; \
		mv fort.12 $$file.11 ; \
		rm -f /tmp/$$file  ; \
	done
endef

