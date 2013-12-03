#!/bin/bash
make
./test
dot -Tpdf output.dot -o output.pdf
