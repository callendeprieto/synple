
all: synple.pdf

clean:
	rm synple.aux synple.dvi synple.log synple.toc synple.pdf

synple.pdf: synple.tex 
	latex synple.tex
	latex synple.tex
	dvipdf synple
	rm synple.aux synple.dvi synple.log synple.toc

