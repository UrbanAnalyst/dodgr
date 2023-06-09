RFILE = README
VIGNETTE = parallel

all: init vignette

init:
	echo "pkgdown::init_site()" | R --no-save -q

site:
	echo "pkgdown::build_site()" | R --no-save q

vignette:
	echo "pkgdown::build_article('$(VIGNETTE)',quiet=FALSE)" | R --no-save -q

knith: $(RFILE).Rmd
	echo "rmarkdown::render('$(RFILE).Rmd',output_file='$(RFILE).html')" | R --no-save -q

doc: clean
	Rscript -e 'devtools::document()'
	Rscript -e 'rmarkdown::render("$(RFILE).Rmd",rmarkdown::md_document(variant="gfm"))'

open:
	xdg-open docs/articles/$(VIGNETTE).html &

check:
	Rscript -e 'library(pkgcheck); checks <- pkgcheck(); print(checks); summary (checks)'

clean:
	rm -rf *.html *.png README_cache docs/
