import os
import sys
import argparse
from collections import Counter
from pprint import pprint

try:
    import numpy as np
except ImportError:
    sys.stderr.write("Error! numpy python library not found!!\n")
    sys.exit(1)

try:
    import pickle
except ImportError:
    sys.stderr.write("Error! pickle package not found!!\n")
    sys.exit(1)
try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("Error! SeqIO not found!!\n")
    sys.exit(1)

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


kmerLen = 31


def getArgs():
    parser = argparse.ArgumentParser(description="", epilog="")
    parser.add_argument("input_fasta", help="This is the input fasta file", type=str)
    parser.add_argument(
        "input_fasta2",
        help="This is the second input fasta file, requires --paired to be designated.",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--db",
        help="This is the location of the .pickle file created in the database creation process",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--input_type",
        help="Designate between fasta and fastq",
        choices=["fasta", "fastq"],
        type=str,
        required=True,
    )
    parser.add_argument(
        "--paired",
        action="store_true",
        help="This is the paired file flag, requires 2 file inputs.",
    )

    return vars(parser.parse_args())


def paired(params):
    """function to concatenate two paired
    files into one and return it back into the script"""
    # Important to know that SeqIO parse is a generator function so once you loop through f1 and f2 you change the type of them.

    # Load files and zip
    fasta1 = SeqIO.parse(open(params["input_fasta"]), params["input_type"])
    fasta2 = SeqIO.parse(open(params["input_fasta2"]), params["input_type"])
    pairedFasta = zip(fasta1, fasta2)

    # Regenerate new fasta file with new string
    records = []
    for f1, f2 in pairedFasta:
        f1.seq = Seq(str(f1.seq) + "NNNNNNNN" + str(f2.seq))
        rec = SeqRecord(f1.seq, id=f1.id)
        records.append(rec)
    SeqIO.write(records, "pairedfasta.tmp", "fasta")


#### Classification steps begin here ####


def strainClassify():
    # Input File, DB directory
    params = getArgs()

    if params["paired"] and not params["input_fasta2"]:
        sys.stderr.write("Require two files for paired option\n")
        sys.exit(0)
    if not params["db"]:
        sys.stderr.write("Database .pickle file is required!\n")
        sys.exit(0)
    if not os.path.isfile(params["input_fasta"]):
        sys.stderr.write(
            " Cannot find the input file at: "
            + params["input_fasta"]
            + "\n Please located the correct file location"
        )
        sys.exit(1)


#  Load the Database
    print("Loading Database..")
    pickle_in = open(params["db"], "rb")
    kdb = pickle.load(pickle_in)

#  Classify to the Database
        print("Classify")
        annot = {}
    numreads = 0

    if params["input_fasta2"]:  # Paired Reads
        paired(params)
        fastaObject = SeqIO.parse(open("pairedfasta.tmp"), params["input_type"])
    else:  # Unpaired Reads
        fastaObject = SeqIO.parse(open(params["input_fasta"]), params["input_type"])

# Read Loop
    for fastaRead in fastaObject:
        numreads += 1
        read_id, read_seq = fastaRead.id, str(fastaRead.seq)
        readTaxList = []
        # K-mers Loop
        for i in range(0, (len(read_seq) - kmerLen + 1)):
            kmer = read_seq[i : i + kmerLen]
            currTaxList = kdb.get(kmer)
            if currTaxList is None:  # Aka k-mer not found
                break
            else:  # Add the taxid of the k-mer
                readTaxList.append(currTaxList)

        """Determine taxa designation based on k-mer hits"""

        # mostCommonTuples[Index][0 = taxid, 1 = count ]
        flattenedTaxa = [y for x in readTaxList for y in x]
        cc = Counter(flattenedTaxa)

        # Case 1: No k-mers found for the entire read
        if not len(flattenedTaxa):
            annot[read_id] = None

        # Case 2: Exactly 1 item for most common
        elif len(cc.most_common()) == 1:
            mostCommonTuples = cc.most_common()
            annot[read_id] = [mostCommonTuples[0][0]]

        # Case 3: Updated Case 2 attempt
        elif len(cc.most_common()) > 1:
            mostCommonTuples = cc.most_common()
            # print(cc.most_common)

            # The first one will have the highest count, set that as the base value
            topCount = mostCommonTuples[0][1]

            # Add this taxid to the dictionary with the read_id as the key
            annot[read_id] = [mostCommonTuples[0][0]]

            # Look at each tuple, check and see if the count is equal
            # if so add it, otherwise stop

            # print("Length of most common tuples", len(mostCommonTuples))
            for i in range(1, len(mostCommonTuples)):
                # if the count is equal to topcount, add, otherwise break
                if mostCommonTuples[i][1] == topCount:
                    annot[read_id].append(mostCommonTuples[i][0])
                else:
                    break

        # Case 4: Neither NA or 1 or >1, aka shouldn't happen
        else:
            print("Error:This shouldn't exist")
            exit(1)

            # annot[read_id] = cc.most_common(1)[0][0]

    fa = open('cc.txt', 'wb')
    pickle.dump(annot, fa)

    print('ANNOT_VALUES\n')
    print(annot.values())
    # Done with read-by-read loop now.

    # grab classified results and convert the values (taxid) to ints
    #priorProb = [int(x) for x in annot.values() if not isinstance(x, list) and isinstance(x, str)]
    # Flatten annot.values() from list of lists to single list
    # priorProb = annot.values()

    # convert flat list to ints
    priorProb = [] #[int(x) for x in priorProb if isinstance(x,str)]
    for x in annot.values():
        if isinstance(x, list) and len(x) is 1:
            for ii in x:
                if not isinstance(ii,type(None)):
                    priorProb.append(ii)


    # Probability dist: First remove the none's and identify unambiguous reads
    #priorProb = [int(x for x in annot.values() if isinstance(x, list) and len(x) == 1]

    # priorProb2 = [x for x in annot.values() if isinstance(x, list) and len(x) == 1]
    # Step 2: Flatten out the list of lists into just 1 list
    # priorProb = [int(indx) for row in priorProb for indx in row]

    """
    for sublist in annot.values():
        if isinstance(sublist,list) and len(sublist) is 1:
            for item in sublist:
                priorProb.append(item)
    """

    #priorProb = [int(x) for x in priorProb]

    print("PRIORPROB\n")
    print(priorProb)
        # if not isinstance(x, list) and isinstance(x, str)]
    # print("After priorProb", priorProb)

    ppCount = Counter(priorProb)
    print("PPCOUNT\n")
    print(ppCount)
    # Output

    # Inference model
    i=0
    for read in annot:
        if isinstance(annot[read], list):
            # annot[read] = [t1,t2,t3]
            probList = [ppCount[x] for x in annot[read]]  # [p1, p2, p3]
            if i%10000 == 0:
                print(type(ppCount[0]), type(annot[read][0]))
                print("problist",i)
                print(probList)
            ind = np.argmax(probList)
            annot[read] = annot[read][ind]  # now its equal to the max prob one
            i+=1
        print("read,annot[read]")
        print(read, annot[read])
    print("PROBLIST\n")
    print(probList)

# Calculate abundance
# convert values to list
    freq = {}
    taxList = list(annot.values())
    for item in taxList:
        if item in freq:
            freq[item] += 1
        else:
            freq[item] = 1

# Filter out low hits
    belowThreshKeys = list()
    THRESH = 0.001
# Generate list of taxa to delete
    for key, value in freq.items():
        if value < THRESH * numreads:
            belowThreshKeys.append(key)

# Go through dict and delete these
    for key in belowThreshKeys:
        if key in freq:
            del freq[key]

# Move Unclassified out of dict
    uReads = freq[None]
    del freq[None]


""" Output stuff """

"""
    TODO:
    This also needs to be organized and cleaned up
    Also add the list of k-mers for each read like kraken does
"""

# Calculate filtered abundances
print("Abundance Info")
for key, value in freq.items():
    abund = round(value / sum(freq.values()), 3)
    print("% s : % s" % (key, abund))

print("Total unfiltered counts")
print(sum(freq.values()))
# currrently numreads is const, but if im filtering i need to add up the values of all the filters and then do numreads - sum(values)

print("Unclassified Reads")
print(uReads)

print("Total # of Reads")
print(numreads)
