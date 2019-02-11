#!/bin/bash



#First, extract strain reads from the previous FASTA file
#TODO LIST: 
# Change 562 for species of interest.
# Skip the oneline part
# Obviously need to generalize it and change names so they're not ecoli



rm DATAFILES/ecoli_reads.fa
taxid=562

#1. This line turns the fasta into a single line
echo "Converting fasta file to single-line fasta"
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' DATAFILES/simreads.fa > DATAFILES/oneline_reads.fa

echo "Extracting reads from taxid $taxid"
#2.This awk line will search the WEVOTE details for a given taxid (classified by WEVOTE) and print out the READ_ID. It will then search for that READ_ID in the original fasta file and feed it into a new fasta
awk -v taxid="$taxid" '{ if ($NF == taxid) {print $1} }' DATAFILES/WEVOTE_Details.txt  | grep -A 1 -f - DATAFILES/oneline_reads.fa   >> DATAFILES/ecoli_reads.fa
echo 'done!'

#3 Execute the kmer code

#echo "Executing main script"
#cd $PWD/src/
#./mladen_main ../DATAFILES/ecoli_reads.fa
