#!/bin/bash

source activate mlbd

python ~/workspace/tesseract/bin/visualize.py \
	trace.txt \
	runtime.png \
	--xaxis np \
	--yaxis runtime \
	--col method \
	--row bsize \
	--sharey
