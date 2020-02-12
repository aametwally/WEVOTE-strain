#!/bin/bash

# Part 2:
# This script is part 2 in database creation: It's goal is to create jellyfish k-mer counted files of each of the genomes downloaded in part 1



dbdir=$1

if [ -z "$2" ]; then
    kmerlen=31
else
    kmerlen=$2
fi

echo $dbdir
cd $dbdir
ls



# Cleanup
if [ -e seq2tax.txt ];
then
    rm seq2tax.txt
fi

mkdir jelly

# Get taxID from seqID
for fasta in library/added/*.fna;
do
    seqID=$(head -n1 $fasta | awk '{print $1}' | sed 's/^.\(.*\)..$/\1/')
    echo -e $seqID'\t'>> seq2tax.tmp
done

if [ -f "seq2tax.txt" ]; then
    echo "seq2tax.txt already exists, delete to recreate mapping"
else
    touch seq2tax.txt
    echo -e "seqID\ttaxID" >> seq2tax.txt
    #  Faster to fgrep from a file than to place into loop
    LC_ALL=C fgrep -f seq2tax.tmp taxonomy/nucl_gb.accession2taxid | awk '{print $1,$3}'>>seq2tax.txt
fi
rm seq2tax.tmp


# Jellyfish Counting for each FASTA file
for fasta in library/added/*.fna;
do
    seqID=$(head -n1 $fasta | awk '{print $1}' | sed 's/^.\(.*\)..$/\1/')
    taxID=$(grep $seqID seq2tax.txt| awk '{print $2}' )
    echo $seqID $taxID
    #  Generate Kmers
    echo  "jellyfish dump for ${fasta} located in jelly/${taxID}.jdb"
    if [ ! -f jelly/$taxID.jdb ]; then
        jellyfish count -m $kmerlen -s 100M $fasta #make kmerdb
        jellyfish dump mer_counts_0 > jelly/$taxID.jdb #  dumb db to a file
        sed -i '/^>/d' jelly/$taxID.jdb #  remove the counter lines
    else
        echo "jdb file for ${taxid} exists, skipping"
    fi
done

