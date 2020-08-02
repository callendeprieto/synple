#
#	makefile for downloading synple/synspec line lists
#	Carlos 2019
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/linelists/
WGET = `which wget`

#general line lists
STANDARD = gfallx3_bpo.19 kmol3_0.01_30.20 

#apogee
HBAND = apogeeDR16.20180901.19 apogee_molec_aboo.20 apogee_molec_sun.20

#output binary files
TMP1 = $(STANDARD:.19=.11)
BINARIES = $(TMP1:.20=.11)
TMP2 = $(HBAND:.19=.11)
HBANDBINARIES = $(TMP2:.20=.11)



all:    $(STANDARD) $(BINARIES)

$(STANDARD):  
	for file in $(STANDARD) ; do \
		${WGET} ${FTP}$$file ; \
	done

$(BINARIES):
	for file in $(STANDARD) ; do \
		echo $$file | ../bin/list2bin ; \
	done

apogee: $(HBAND) $(HBANDBINARIES)

$(HBAND):
	for file in ${HBAND} ; do \
		echo $$file | ../bin/list2bin ; \
	done

$(HBANDBINARIES): $(HBAND)
	for file in ${HBAND} ; do \
		echo $$file | ../bin/list2bin ; \
	done


