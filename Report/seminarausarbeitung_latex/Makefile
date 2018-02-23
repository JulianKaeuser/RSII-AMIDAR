TEX = pdflatex

TEXFILES := $(wildcard *.tex)

all: report.pdf

report.pdf: $(TEXFILES)
	pdflatex --halt-on-error "report.tex"
	pdflatex --halt-on-error "report.tex"