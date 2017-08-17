LL=latex
PDF=dvipdfm
LFLAGS=-interaction=nonstopmode -shell-escape
FILE = contraction


all: latex latex dvi
#all: latex latex bib pdf

latex: 
	$(LL) $(FILE).tex $(LFLAGS) 

bib:
	bibtex $(FILE).aux

pdf:
	$(PDF) $(FILE).dvi

dvi:
	dvips -P pdf -q $(FILE).dvi
	ps2pdf $(FILE).ps

open:
	xdg-open $(FILE).pdf &

clean:
	rm -f *.aux *.dvi *.log *.out *.ps
