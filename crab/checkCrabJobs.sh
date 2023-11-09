#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Script checking if an input list of job names contain failed jobs"
  echo "Usage of the script:"
  echo "$0 inputList [-r] [-t]"
  echo "inputList = List of all the jobs that were submitted"
  echo "-r = Try to resubmit all failed jobs"
  echo "-t = New time in minutes allowed for jobs"
  exit
fi

INPUTLIST=$1
TEMPFILENAME="veryTempFile.txt"
shift

# Read the optional arguments. (Semicolon after letter: expects argument)
while getopts ":rt:" opt; do
case $opt in
r) RESUBMIT=true
;;
t) TIME="$OPTARG"
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values for optional arguments
RESUBMIT=${RESUBMIT:-false}
TIME=${TIME:-0}

# Loop over the submitted jobs
while read -r JOBNAME
do
  crab status -d "${JOBNAME}/crab_${JOBNAME}" > $TEMPFILENAME

  # Print the job status summary table
  echo ${JOBNAME}
  awk '/Jobs status/ {p=1}; /^$/ {p=0}; p' "$TEMPFILENAME"

  # Option to resubmit failed jobs
  if $RESUBMIT; then

    # Only resubmit jobs if the status table contains failed jobs
    if grep -q "Error Summary" $TEMPFILENAME; then

      # Check if we want to give the jobs more time
      OPTIONSTRING=""
      if [ $TIME -gt 0 ]; then
        OPTIONSTRING+=" --maxjobruntime ${TIME}"
      fi

      echo "Resubmitting failed jobs!"
      if [ OPTIONSTRING != "" ]; then
        echo "Options for resubmit: ${OPTIONSTRING}"
      fi
      crab resubmit -d "${JOBNAME}/crab_${JOBNAME}" $OPTIONSTRING
    fi
  fi

  echo ""

done < <(cat $INPUTLIST)

# Remove the temporary helper file
rm $TEMPFILENAME
