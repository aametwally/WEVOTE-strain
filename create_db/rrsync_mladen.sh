#!/bin/bash


#File to rsync all of the fna.gz files for a given taxid



TAXID=562
echo $TAXID
awk -v taxid="$TAXID" '{if ($7==taxid) print $NF;}' assembly_summary.txt > urls.tmp

#remove first part so rsync from file works
sed -i 's,ftp://ftp.ncbi.nlm.nih.gov/genomes/,,g' urls.tmp
#append last part to make a file
awk -F/ '{print $0"/"$6"_genomic.fna.gz"}' urls.tmp > urls.tmp.tmp && mv urls.tmp.tmp urls.tmp

#for debugging, only work on 10 files at a time
#head -n 10 urls.tmp > urls.tmp.tmp && mv urls.tmp.tmp urls.tmp

rsync --dry-run --no-motd --files-from=urls.tmp rsync://ftp.ncbi.nlm.nih.gov/genomes/ . 2> rsync.err

rsync --no-motd --files-from=urls.tmp rsync://ftp.ncbi.nlm.nih.gov/genomes/ .




