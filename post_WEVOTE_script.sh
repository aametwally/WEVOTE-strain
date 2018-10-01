#!/bin/bash



#First, extract strain reads from the previous FASTA file
#TODO LIST: 
# Change 562 for species of interest.
# Skip the oneline part




#1. This line turns the fasta into a single line
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' DATAFILES/reads_converted.fasta > DATAFILES/oneline_reads.fa
echo 'converted fasta into a single line fasta!'

#2.This awk line will search the WEVOTE details for a given taxid (classified by WEVOTE) and print out the READ_ID. It will then search for that READ_ID in the original fasta file and feed it into a new fasta
awk '{ if ($NF == 562) {print $1} }' DATAFILES/WEVOTE_Details.txt  | grep -A 1 -f - DATAFILES/oneline_reads.fa   > DATAFILES/ecoli_reads.fa
echo 'done!'

#3 Execute the kmer code

cd $PWD/my_kraken/src/
./mladen_main ../
