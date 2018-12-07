#!/bin/bash


# Previous file - download-taxonomy - same as kraken's no edits

#Current goal - Nov 28th - attempt to make a functional script that takes a SPECIES taxid as input, and returns the downloaded genomes of all children that can then be built into a database
#alternate options --implement genus download, --add option for representative genomes

#Next file - build database - kraken file as well


set -u
set -e


#LIBRARY_DIR = "$KRAKEN_DB_NAME/library"
#LIBRARY_DIR = "mladen/library"
NCBI_SERVER = "ftp.ncbi.nlm.nih.gov"
FTP_SERVER = "ftp://$NCBI_SERVER"
RSYNC_SERVER = "rsync://$NCBI_SERVER"
THIS_DIR = $PWD


library_name = $TAXID
library_file="library.fna"

#if [ -e "$LIBRARY_DIR/$library_name/.completed" ]; then
#    echo "Skipping $library_name, already completed library download"
#    exit 0
#fi

echo "$library_name"
mkdir -p $LIBRARY_DIR/$library_name #-p means to make parentdir if needed
cd $LIBRARY_DIR/$library_name
rm -f assembly_summary.txt
remote_dir_name = $library_name

if ! wget -q $FTP_SERVER/genomes/refseq/$library_name/assembly_summary.txt; then
    echo "can't get assembly summary file for $library_name, exiting" > /dev/fd/2
    exit 1
fi

rm -rf all/ library.f* manifest.txt rsync.err
#rsync_from_ncbi.pl assembly_summary.txt
#scan_fasta_file.pl $library_file > prelim_map.txt

#fetch assembly_summary.txt to match taxid and pull it into something.
TAXID= 562

awk -v taxid="$TAXID" '{if ($7==taxid) print $NF;}' assembly_summary.txt > urls.tmp






touch .completed







