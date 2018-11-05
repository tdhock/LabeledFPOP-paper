HOCKING-labeled-FPOP.pdf: HOCKING-labeled-FPOP.tex refs.bib
	rm -f *.aux *.bbl
	pdflatex HOCKING-labeled-FPOP
	bibtex HOCKING-labeled-FPOP
	pdflatex HOCKING-labeled-FPOP
	pdflatex HOCKING-labeled-FPOP
