#
#	makefile for downloading synple/synspec line lists
#	Carlos 2019
#

FTP = ftp://carlos:allende@ftp.ll.iac.es/linelists/
WGET = /bin/wget

STANDARD = gfallx3_bpo.19 kmol3_0.01_30.20
HBAND = apogeeDR16.20180901.19 apogeeDR16_arc.20 apogeeDR16_sun.20

all:
	for file in ${STANDARD} ; do \
		${WGET} ${FTP}$$file ; \
	done

apogee:

	for file in ${HBAND} ; do \
		${WGET} ${FTP}$$file ; \
	done