#!/bin/bash -e

# this will log the environment, useful for debugging!
(set -o posix; set; ulimit -a; uname -a; lsb_release -a; hostname -A) 1>&2

module load phuluu/bedtools/2.29.2 elypar/ChromHMM/1.23 elypar/java/jdk1.8.0_25

CHROMHMM="/share/ScratchGeneral/elypar/ChromHMM/"
PWD=$(pwd)

unset DISPLAY

# Learn the model with 5-15 states
for i in `seq 2 6`
do
  CMD="java -d64 -Xmx16G -jar $CHROMHMM/ChromHMM.jar LearnModel -p 8 ${PWD}/Binarized_data/ "$i"_states "$i" hg38"
  echo $CMD & eval $CMD
done

# Compare models
for i in `seq 2 6`
do
  mkdir "$i"_comparison
  for j in `seq 2 $i`
  do
    cp "$j"_states/emissions_"$j".txt "$i"_comparison/
  done
  cd "$i"_comparison
  java -d64 -Xmx16G -jar $CHROMHMM/ChromHMM.jar CompareModels emissions_"$i".txt ./ "$i"_comparison
  cd ..
done  

