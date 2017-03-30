#!/bin/bash
#Replacing the Picard version due to error:Exception in thread "main" java.lang.SecurityException: Invalid signature file digest for Manifest main attributes, when using current version
loc=$1
id=`basename $loc |sed 's/.bam//' ` 
mkdir $id
java -jar -Xmx16g /seq/software/picard/current/bin/picard.jar SamToFastq \
INPUT=${loc} \
OUTPUT_PER_RG=TRUE \
OUTPUT_DIR=./$id \
NON_PF=True

echo "Started gzipping"
gzip ./$id/*
echo "Finished gzipping"
