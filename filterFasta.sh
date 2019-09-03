#!/bin/bash

#USAGE: ./filterFasta.sh $TAXID $WEVOTE_DETAILS_FILE $ORIGINAL_FASTA $OUTFILE


# Input and Exception Handling
while getopts "t:l:i:o:" opt; do
  case $opt in
    t)
        if ! [[ "$OPTARG" =~ ^[0-9]+$ ]]
            then
                echo "TaxID must be an integer. Exiting." >&2
                exit 0
        fi
        taxid=$OPTARG
      ;;
    l)
        wevote_file=$OPTARG
        ;;
    i) 
        og_fasta=$OPTARG
        ;;
    o) 
        fa=".fa"
        out_fasta_prefix=$OPTARG
        out_fasta="$out_fasta_prefix$fa"

        ;;
    \?)
        echo "Error: Unknown argument."
        echo "Usage: ./post_WEVOTE_script.sh -t <TAXID> -l <WEVOTE_DETAILS_FILE> -i <ORIGINAL_FASTA> -o <OUTFILE>" >&2
        exit 0
      ;;
  esac
done

if [ $OPTIND -ne 9 ]
then
    echo $OPTIND
    echo "Error: Must provide exactly four arguments, exiting."
    echo "Usage: ./post_WEVOTE_script.sh -t <TAXID> -l <WEVOTE_DETAILS_FILE> -i <ORIGINAL_FASTA> -o <OUTFILE>"

    exit 0
fi

# Clear old files
rm -f $out_fasta

# This line turns the fasta into a single line
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' \
    $og_fasta > linearfasta.tmp



# This awk line will search the WEVOTE details for a given taxid (classified by WEVOTE) and print out the READ_ID. It will then search for that READ_ID in the original fasta file and feed it into a new fasta
echo ""
echo "Extracting reads from taxid $taxid ...."



# Store into array to print out # of matches
myarr=($(awk -v taxid="$taxid" '{ if ($NF == taxid) {print $1} }' $wevote_file ))
len=${#myarr[@]}
echo " $len matches found.."
echo " Creating new fasta from $len hits"


# Take output and put into tmp file for grep
awk -v taxid="$taxid" '{ if ($NF == taxid) {print $1} }' $wevote_file > match_reads.tmp

# Use grep -F "Extremely fast"
grep -A 1 -F -w -f  match_reads.tmp linearfasta.tmp > $out_fasta

echo "Done! New Fasta file: $out_fasta"

rm linearfasta.tmp match_reads.tmp

