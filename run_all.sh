#!/bin/bash

for input_file in $(ls *.dat); do
	echo $input_file
	./vortex_solver.exe $input_file &
done
