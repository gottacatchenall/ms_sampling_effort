#! /usr/bin/bash


mkdir -p tmp
cd tmp

pandoc ../suppl.md \
    -s -o report.tex \
    -F \
    ~/apps/pandoc-crossref  \
    -H ../header.tex

pdflatex report.tex

mv report.pdf ../supplement.pdf 

cd .. 
rm -rf tmp
