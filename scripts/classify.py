import os
import sys
import argparse
from collections import Counter

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



kmerLen = 31



def getArgs():
    parser = argparse.ArgumentParser(description='',epilog='')
    parser.add_argument("input_fasta", help = "This is the input fasta file",type=str)
    parser.add_argument("--db", help = "This is the location of the .pickle file created in the database creation process",type=str,required=True)
    args = parser.parse_args()
    return args.input_fasta, args.db

def paired():
    """function to concatenate two paired files into one and return it back into the script"""
    fasta1 = SeqIO.parse(open(file1),'fasta')
    fasta2 = SeqIO.parse(open(file2),'fasta')

    for i,j in zz:
        nn = str(i.seq) +"NNNNNNNN"+ str(j.seq)
        print(nn)
        # This works, but need to know what to return back






# Input File, DB directory
input_file,pickle_dir = getArgs()


#  Load the Database
print('Loading Database..')
pickle_in = open(pickle_dir, "rb")
kdb = pickle.load(pickle_in)

#  Classify to the Database
print('Classify')
annot = {}
numreads=0
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    numreads += 1
    read_id, read_seq = fasta.id, str(fasta.seq)
    readTaxList=[]
    for i in range(0,(len(read_seq)-kmerLen+1)):
        kmer = read_seq[i:i+kmerLen]
        currTaxList = kdb.get(kmer)
        if currTaxList is None:
            break
        else:
            readTaxList.append(currTaxList)

    flattenedTaxa = [y for x in readTaxList for y in x]
    cc = Counter(flattenedTaxa)
    if not len(flattenedTaxa): # no matching kmers
        annot[read_id] = None
    elif len(cc.most_common())>1 and (cc.most_common()[0][1] == cc.most_common()[1][1]): # Tie
        annot[read_id] = [cc.most_common()[0][0],cc.most_common()[1][0]]
    elif len(cc.most_common())>1 and (cc.most_common()[0][1] == cc.most_common()[1][1] == cc.most_common()[2][1]): # Tie
        annot[read_id] = [cc.most_common()[0][0],cc.most_common()[1][0],cc.most_common()[2][0]]
    else:
        annot[read_id] = cc.most_common(1)[0][0]

# grab classified results and convert to ints
priorProb = [ int(x) for x in annot.values() if not isinstance(x, list) and isinstance(x, str)]
ppCount = Counter(priorProb)
# Output

# Inference model
for read in annot:
    if isinstance(annot[read], list):
        #annot[read] = [t1,t2,t3]
        probList = [ppCount[int(x)] for x in annot[read]] #[p1, p2, p3]
        ind = np.argmax(probList)
        annot[read] = annot[read][ind] #now its equal to the max prob one
    print(read, annot[read])

# Calculate abundance
# convert values to list
freq = {}
taxList = list(annot.values())
for item in taxList:
    if (item in freq):
        freq[item] += 1
    else:
        freq[item] = 1

#Filter out low hits
belowThreshKeys = list()
THRESH = 0.01
# Generate list of taxa to delete
for key,value in freq.items():
    if value < THRESH*numreads:
        belowThreshKeys.append(key)

# Go through dict and delete these
for key in belowThreshKeys:
    if key in freq:
        del freq[key]

# Move Unclassified out of dict
uReads = freq[None]
del freq[None]



# Calculate filtered abundances
print("Abundance Info")
for key,value in freq.items():
    abund = round(value/sum(freq.values()),3)
    print ( "% s : % s" %(key,abund) )

print("Total unfiltered counts")
print(sum(freq.values()))
# currrently numreads is const, but if im filtering i need to add up the values of all the filters and then do numreads - sum(values)

print("Unclassified Reads")
print(uReads)

print("Total # of Reads")
print(numreads)






