#!/bin/bash

ISSFILE='issr1.fastq'
readID=$(grep '^@' $ISSFILE) #for the actual file
searchID=$(grep '^@' issr1.fastq | sed 's/\..*//' | sed 's/^@//') #gets you the searchable readID from iss file
if [ -e mapping.txt ]; then
    rm mapping.txt
fi

for i in $searchID
do
    tmp=$(LC_ALL=C fgrep $i taxtest_hflu_dir/seq2tax.txt| awk '{print $2}')
    taxID+="$tmp\n"
done

#LC_ALL=C fgrep -f mapping.txt taxtest_hflu_dir/seq2tax.txt| awk '{print $2}'
paste <(printf %s "$readID") <(printf %s "$taxID")
#echo -e "$readID\t$taxID" > mapping.txt
