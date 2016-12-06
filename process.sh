#!/bin/sh

for f in `ls ./doc/documentation/*.ipynb`; do
	jupyter nbconvert --to html $f 
done

mv ./doc/documentation/*html ./website/doc/documentation_html/

for f in `ls ./doc/tutorials/*.ipynb`; do
	jupyter nbconvert --to html $f 
done

mv ./doc/tutorials/*html ./website/doc/tutorial_html/
