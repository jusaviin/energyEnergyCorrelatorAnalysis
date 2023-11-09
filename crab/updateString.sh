#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "Usage of the script:"
  echo "$0 inputString elementsToReplace replaceValues"
  echo "inputString = String that could benefit for some elements to be altered."
  echo "elementsToReplace = List of element indices in the string to be replaces."
  echo "replaceValues = Values used to replace the text in defined indices."
  exit
fi

INPUT=($1)
ELEMENTS=($2)
UPDATES=($3)

COUNTER=0
for INDEX in ${ELEMENTS[*]}; do
  INPUT[$INDEX]=${UPDATES[${COUNTER}]}
  ((COUNTER++))
done

OUTPUT="${INPUT[*]}"
echo $OUTPUT
