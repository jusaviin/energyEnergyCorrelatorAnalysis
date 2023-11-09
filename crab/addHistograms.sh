#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Script for adding all root files from a folder in EOS"
  echo "Usage of the script:"
  echo "$0 [folderName] [baseName]"
  echo "folderName = Name of the folder where the root files are"
  echo "baseName = Name given for the output file"
  exit
fi

FOLDERNAME=${1%/}
BASENAME=$2

hadd -ff ${BASENAME}.root `xrdfs root://cmseos.fnal.gov ls -u $FOLDERNAME | grep '\.root'`
