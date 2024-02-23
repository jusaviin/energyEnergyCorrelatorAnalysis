#!/bin/bash

for i in $(seq 1 2); do
  for j in $(seq 0 0); do
    root -l -b -q 'plotting/determineNumberOfUnfoldingIterations.C('${i}','${j}')'
  done
done
