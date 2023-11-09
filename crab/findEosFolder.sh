#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Find the path in which the root files are in EOS"
  echo "Usage of the script:"
  echo "$0 folder"
  echo "folder = Lowest level folder at EOS"
  exit
fi

FOLDER=$1

FIRSTLEVEL=(`eos root://cmseos.fnal.gov ls /store/user/jviinika/${FOLDER}`)
SECONDLEVEL=(`eos root://cmseos.fnal.gov ls /store/user/jviinika/${FOLDER}/${FIRSTLEVEL}`)
THIRDLEVEL=(`eos root://cmseos.fnal.gov ls /store/user/jviinika/${FOLDER}/${FIRSTLEVEL}/${SECONDLEVEL}`)
FOURTHLEVEL=(`eos root://cmseos.fnal.gov ls /store/user/jviinika/${FOLDER}/${FIRSTLEVEL}/${SECONDLEVEL}/${THIRDLEVEL}`)


echo "/store/user/jviinika/${FOLDER}/${FIRSTLEVEL}/${SECONDLEVEL}/${THIRDLEVEL}/${FOURTHLEVEL}"
