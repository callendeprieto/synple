#
#	makefile for downloading synple/synspec line lists
#	Carlos 2019
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/linelists/
WGET = `which wget`
GUNZIP = `which gunzip`

#general line lists
STANDARD = gfallx3_bpo.19 kmol3_0.01_30.20

#apogee
HBAND = synspec20200709.lte.19 synspec20200705.nlte.19 apogee.20180901.solar.20 apogee.20180901.giants.20 H2O_apogee_TSl3500-8.5_nobroad.20

#output binary files
TMP1 = $(STANDARD:.19=.11)
BINARIES = $(TMP1:.20=.11)
TMP2 = $(HBAND:.19=.11)
HBANDBINARIES = $(TMP2:.20=.11)



all:    $(STANDARD) $(BINARIES)

$(STANDARD):  
	for file in $(STANDARD) ; do \
		${WGET} ${FTP}$$file.gz ; \
                ${GUNZIP} $$file.gz ; \
	done

$(BINARIES):
	for file in $(STANDARD) ; do \
		echo $$file | ../bin/list2bin ; \
	done

apogee: $(HBAND) $(HBANDBINARIES)

$(HBAND):
	for file in ${HBAND} ; do \
		${WGET} ${FTP}$$file.gz ; \
                ${GUNZIP} $$file.gz ; \
	done

$(HBANDBINARIES): $(HBAND)
	for file in ${HBAND} ; do \
		echo $$file | ../bin/list2bin ; \
	done


