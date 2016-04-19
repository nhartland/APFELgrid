default:
	pdflatex paper
	bibtex   paper
	pdflatex paper
	pdflatex paper

clean:
	rm -rf *.bbl *~ *.out *.toc *.log *.pdf *.aux *.blg *.spl
