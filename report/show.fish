#!/usr/bin/fish

latexmk -pdf -pvc report.tex &
zathura report.pdf

