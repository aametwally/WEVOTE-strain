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

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



kmerLen = 31



def getArgs():
    parser = argparse.ArgumentParser(description='',epilog='')
    parser.add_argument("input_fasta", help = "This is the input fasta file",type=str)
    parser.add_argument("input_fasta2", help = "This is the second input fasta file, requires --paired to be designated.",type=str, nargs='?')
    parser.add_argument("--db", help = "This is the location of the .pickle file created in the database creation process",type=str,required=True)
    parser.add_argument("--input_type", help ="Designate between fasta and fastq",choices=['fasta','fastq'],type=str,required=True)
    parser.add_argument("--paired",action='store_true', help = "This is the paired file flag, requires 2 file inputs.")

    return vars(parser.parse_args())

def paired(params):
    """function to concatenate two paired files into one and return it back into the script"""
    #Important to know that SeqIO parse is a generator function so once you loop through f1 and f2 you change the type of them.

    # Load files and zip
    fasta1 = SeqIO.parse(open(params['input_fasta']),params['input_type'])
    fasta2 = SeqIO.parse(open(params['input_fasta2']),params['input_type'])
    pairedFasta = zip(fasta1,fasta2)

    # Regenerate new fasta file with new string
    records = []
    for f1,f2 in pairedFasta:
        f1.seq = Seq(str(f1.seq) +"NNNNNNNN"+ str(f2.seq))
        rec = SeqRecord(f1.seq,id=f1.id)
        records.append(rec)
    SeqIO.write(records, 'pairedfasta.tmp', 'fasta')



#### Classification steps begin here ####


# Input File, DB directory
params = getArgs()

if params['paired'] and not params['input_fasta2']:
    sys.stderr.write('Require two files for paired option\n')
    sys.exit(0)
if not params['db']:
    sys.stderr.write('Database .pickle file is required!\n')
    sys.exit(0)
if not os.path.isfile(params['input_fasta']):
    sys.stderr.write(" Cannot find the input file at: " + params['input_fasta'] +
            "\n Please located the correct file location")
    sys.exit(1)






#  Load the Database
print('Loading Database..')
pickle_in = open(params['db'], "rb")
kdb = pickle.load(pickle_in)

#  Classify to the Database
print('Classify')
annot = {}
numreads=0
if params['input_fasta2']:
    paired(params)
    fasta_sequences = SeqIO.parse(open('pairedfasta.tmp'),params['input_type'])
else:
    fasta_sequences = SeqIO.parse(open(params['input_fasta']),params['input_type'])

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

    """Determine taxa designation based on k-mer hits"""
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






