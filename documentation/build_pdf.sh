#!/usr/bin/env bash
DOCTEX=AMOEBAE_documentation.tex
pdflatex --shell-escape $DOCTEX &&\
bibtex $DOCTEX &&\
makeglossaries $DOCTEX &&\
pdflatex --shell-escape $DOCTEX &&\
pdflatex --shell-escape $DOCTEX &&\
open AMOEBAE_documentation.pdf
