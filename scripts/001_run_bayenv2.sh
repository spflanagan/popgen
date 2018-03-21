#!/bin/bash

for i in {1..10}
do
	echo "Running rep $i"
	./bayenv2 -i SNPSFILE -p 12 -k 100000 -r 628398 > matrix.$i.out
done
