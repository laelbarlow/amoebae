#!/usr/bin/env bash
DOCTEX=amoebae_commands
pdflatex --shell-escape $DOCTEX.tex &&\
bibtex $DOCTEX &&\
makeglossaries $DOCTEX &&\
pdflatex --shell-escape $DOCTEX.tex &&\
pdflatex --shell-escape $DOCTEX.tex
