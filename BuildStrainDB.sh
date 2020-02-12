#!/bin/bash


# Need Arguments
if [ $# -eq 0 ]
then

    echo "Missing Options!"
    echo "(run $0 -h for help)"
    echo ""
    exit 0
fi

# Input variables
taxid=$1
dbname=$2

if [ -z "$3" ]; then
    kmerlen=31
else
    kmerlen=$3
fi

dbfolder="${dbname}_dir"
dbfolder_jelly="${dbfolder}/jelly"


echo $taxid $dbname $dbfolder $dbfolder_jelly

# Go to folder containing script
scripts_dir="$(dirname "$0")"
scripts_dir="$(readlink -f $scripts_dir)" #make absolute
echo $scripts_dir



# Run
echo "Step 1: Download Taxonomic and Genomic References"
python $scripts_dir"/bsdb_dlGenomes.py"  $taxid $dbname
echo "Step 1 Complete: Files stored in ${dbname}_dir"

pwd

echo "Step 2: "Generate K-mers for each strain""
bash $scripts_dir"/bsdb_countKmers.sh" $dbfolder $kmerlen
echo "Step 2 Complete: K-mer dumps stored in $dbfolder_jelly"

echo "after step 2 location"
pwd

cd $dbfolder_jelly

echo "Step 3:"K-mer filtering and converting to python dictionary""
python $scripts_dir"/bsdb_hashKmers.py"
echo "Step 3 Complete: "
