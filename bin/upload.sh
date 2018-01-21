#!/bin/bash
set -e

FOLDER=/media/4TB1/blastocrithidia/annotation/p57
DATE=$(date +"%B_%d")
FOLDER=$FOLDER/$DATE

i=1
while ssh trypka "[ -d $FOLDER ]"
do
  i=$((i+1))
  if [ $i = 2 ]; then
    FOLDER=${FOLDER}_${i}
  else
    FOLDER=$(echo $FOLDER | sed "s#[0-9]\+\$#$i#")
  fi
done

ssh trypka "mkdir $FOLDER"
scp results/annotator/*.{txt,gff} trypka:$FOLDER
