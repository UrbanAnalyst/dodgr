LFILE = dists-categorical

all: knith
#all: knith open 

knith: $(LFILE).Rmd
	echo "rmarkdown::render('$(LFILE).Rmd',output_file='$(LFILE).html')" | R --no-save -q

knitr: $(LFILE).Rmd
	echo "rmarkdown::render('$(LFILE).Rmd',rmarkdown::md_document(variant='gfm'))" | R --no-save -q

open: $(LFILE).html
	xdg-open $(LFILE).html &

clean:
	rm -rf *.html *.png 



fig1: pdflone convone rmone

pdflone: fig1.tex
	pdflatex fig1.tex

convone: fig1.pdf
	convert -density 300 fig1.pdf -quality 100 fig1.png

rmone:
	rm fig1.aux fig1.log fig1.pdf

fig2: pdfltwo convtwo rmtwo

pdfltwo: fig2.tex
	pdflatex fig2.tex

convtwo: fig2.pdf
	convert -density 300 fig2.pdf -quality 100 fig2.png

rmtwo:
	rm fig2.aux fig2.log fig2.pdf

fig3: pdflthree convthree rmthree

pdflthree: fig3.tex
	pdflatex fig3.tex

convthree: fig3.pdf
	convert -density 300 fig3.pdf -quality 100 fig3.png

rmthree:
	rm fig3.aux fig3.log fig3.pdf

figt1: pdflone_t convone_t rmone_t

pdflone_t: times-fig1.tex
	pdflatex times-fig1.tex

convone_t: times-fig1.pdf
	convert -density 300 times-fig1.pdf -quality 100 times-fig1.png

rmone_t:
	rm times-fig1.aux times-fig1.log times-fig1.pdf
