LFILE = README
VIGNETTE = parallel

all: init vignette

init:
	echo "pkgdown::init_site()" | R --no-save -q

site:
	echo "pkgdown::build_site()" | R --no-save q

vignette:
	echo "pkgdown::build_article('$(VIGNETTE)',quiet=FALSE)" | R --no-save -q

knith: $(LFILE).Rmd
	echo "rmarkdown::render('$(LFILE).Rmd',output_file='$(LFILE).html')" | R --no-save -q

knitr: $(LFILE).Rmd
	echo "rmarkdown::render('$(LFILE).Rmd', output_format = rmarkdown::md_document (variant = 'gfm'))" | R --no-save -q

open:
	xdg-open docs/articles/$(VIGNETTE).html &

check:
	Rscript -e 'library(pkgcheck); checks <- pkgcheck(); print(checks); summary (checks)'

clean:
	rm -rf *.html *.png README_cache docs/
