#################
# documentation #
#################

RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)

doc: $(HTML)

%.html: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

%.docx: %.rmd $(COMMON)
	Rscript -e "library(rmarkdown); render('$<', 'word_document');"

