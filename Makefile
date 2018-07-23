all: project.pdf

project.pdf: project.tex sources.bib fig.pdf rdf.py
	latexmk -pdflatex=lualatex -pdf
fig.pdf: fig.asy
	asy -maxtile "(400, 400)" fig.asy
clean:
	latexmk -c
	rm -rf __pycache__ pythontex-files-project *.pytxcode *.auxlock *.run.xml data *.bbl tmp/project-figure* *.figlist *.makefile latex.out *.mod *.xyz
