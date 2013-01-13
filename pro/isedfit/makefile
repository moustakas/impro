ALL= \
	iSEDfit_cookbook.pdf

%.dvi: %.tex
	latex $*; latex $*

%.pdf : %.dvi
	dvips -Ppdf -N0 -o $*.ps $*.dvi
	ps2pdf $*.ps
	rm -f $*.ps $*.aux $*.log $*.dvi $*.tex~

all: $(ALL)

clean:
	rm *.log *.dvi *.aux *.bbl *.blg

fullclean:
	rm *.pdf *.log *.dvi *.aux *.bbl *.blg

dummy:
