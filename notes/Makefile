PDF_LATEX=pdflatex
BIB_TEX=bibtex
REFS=$(wildcard *.bib)

# Figure out the right RM for win/unix
ifdef ComSpec
	RM=del /F /Q
else
	RM=rm -f
endif

all: notes.pdf midterm_report.pdf

midterm_report.pdf: midterm_report.tex
	$(PDF_LATEX) midterm_report

%.pdf: %.tex $(REFS)
	$(PDF_LATEX) $*
	$(BIB_TEX) $*  
	$(PDF_LATEX) $*
	$(PDF_LATEX) $*


.PHONY: clean
clean:
	$(RM) *.aux; $(RM) *.dvi; $(RM) *.log; $(RM) *.pdf; $(RM) *.gz; \
		$(RM) *.bbl; $(RM) *.blg; $(RM) *.out; $(RM) *.lbl; $(RM) *.brf


