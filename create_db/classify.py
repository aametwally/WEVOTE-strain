from Bio import SeqIO
import pickle
from collections import Counter

input_file = "hflu_filter.fa"
kmerLen = 31



#  Load the kmer dictionary
print('loading pickle')
pickle_in = open("hflu_genomes_dir/jelly/thin_dict.pickle", "rb")
kdb = pickle.load(pickle_in)

#  Classify
annot = {}
numreads=0
print('classify')
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    numreads += 1
    read_id, read_seq = fasta.id, str(fasta.seq)
    # print(read_id)
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
    # print(cc)  # Returns the highest occurring item
    if not len(flattenedTaxa): # no matching kmers
        annot[read_id] = None
    elif len(cc.most_common())>1 and (cc.most_common()[0][1] == cc.most_common()[1][1]): # Tie
        annot[read_id] = [cc.most_common()[0][0],cc.most_common()[1][0]]
    elif len(cc.most_common())>1 and (cc.most_common()[0][1] == cc.most_common()[1][1] == cc.most_common()[2][1]): # Tie
        annot[read_id] = [cc.most_common()[0][0],cc.most_common()[1][0],cc.most_common()[2][0]]
    else:
        annot[read_id] = cc.most_common(1)[0][0]

# Output
for read in annot:
    print(read, annot[read])
print(numreads)


