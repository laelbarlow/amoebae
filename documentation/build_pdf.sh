#!/usr/bin/env bash
DOCTEX=AMOEBAE_documentation
pdflatex --shell-escape $DOCTEX.tex &&\
bibtex $DOCTEX &&\
makeglossaries $DOCTEX &&\
pdflatex --shell-escape $DOCTEX.tex &&\
pdflatex --shell-escape $DOCTEX.tex &&\
open AMOEBAE_documentation.pdf
