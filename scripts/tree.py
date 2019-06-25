from __future__ import division
import time

def generate_node_parent(dmp_file):
    """ Creates a dictionary with node::parent as key:value pairing"""
    d={}
    with open(dmp_file) as f:
        for line in f:
            row = line.split("|")
            tax_id = int(row[0])
            parent_id = int(row[1])
            d[(tax_id)] = parent_id
    return d

def generate_node_rank(dmp_file):
    """ Creates a dictionary with node::rank as key:value pairing"""
    d_rank={}
    with open(dmp_file) as f:
        for line in f:
            row = line.split("|")
            tax_id = int(row[0])
            rank_id = row[2]
            d_rank[(tax_id)] = rank_id
    return d_rank

def traverse_tree(taxid):
    """Given a taxid and a graph find each parent in the graph"""
    phylo_list=[]
    phylo_list.append(taxid)
    while taxid != 1:
        taxid = parent_tree[taxid]
        phylo_list.append(taxid)
    phylo_list.reverse() #needed to line up columns
    return phylo_list






#Future function to generate the true reads values
def write_trueread_file(reads_file,outfile):
    """ Given a file of reads and their IDs, write a file that lists
    their read,ID, parentID, parentID, etc. up to the root"""
    tax_id = []
    read_id = []
    outF = open(outfile, "w+")
    outF.write("Read\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n")

    with open(reads_file) as f:
        for line in f:
            colvals = line.split(" ")
            #Generate a list of read_id and tax_id
            read_id = colvals[0]
            tax_id = int(colvals[5])
            if tax_id == 0:
                output_list = [0]*6
            else:
                phylogeny_list = traverse_tree(tax_id)
                outF.write(read_id + "\t")

                output_list = [None]*6

                #Cleanup unranked items, Only output definitive ranks
                for item in phylogeny_list:
                    if rank_tree[item] == "\tphylum\t":
                        output_list[0] = item
                    elif rank_tree[item] == "\tclass\t":
                        output_list[1] = item
                    elif rank_tree[item] == "\torder\t":
                        output_list[2] = item
                    elif rank_tree[item] == "\tfamily\t":
                            output_list[3] = item
                    elif rank_tree[item] == "\tgenus\t":
                        output_list[4] = item
                    elif rank_tree[item] == "\tspecies\t":
                        output_list[5] = item
                    else:
                        phylogeny_list.remove(item)
            for i in output_list:
                outF.write(str(i))
                outF.write("\t")
            outF.write("\n")

def precisionSensitivity(outfile,truefile,rankTree,parentTree):
    """Takes in wevotestrain output and compare to true values"""
    TP = 0
    FP = 0
    FN = 0
    trueTaxID = 0

    #Generate dictionary of mapping with d[read_ID] = taxid
    trueDict = {}
    with open(truefile) as tf:
        for line in tf:
            items = line.strip().split("\t")
            key, values = items[0], items[1]
            trueDict[key] = values


    with open(outfile) as f:
        for line in f:
            strainFlag = 0 #1 if strain
            colvals = line.strip().split("\t")
            readID = colvals[1]
            taxID = colvals[2]
            if colvals[0] == 'None':
                FN += 1
                continue
            rank = rankTree[int(taxID)].strip()
            parentRank = rank
            parentID = taxID
            if rank == 'no rank': #determine if strain
                while parentID != 1 and strainFlag ==0:
                    parentRank = rankTree[int(parentID)].strip()
                    parentID = parentTree[int(parentID)]
                    #while you traverse up, if you hit species at any point then you started with a strain
                    if parentRank == 'species':
                        strainFlag = 1

                if strainFlag == 1:
                    #compare for match
                    if taxID == trueDict[readID]:
                        TP += 1
                    else:
                        FP += 1
            else:
                FN += 1

    #false negatives that were lost in wevote
    FN2 = 600 - (TP + FP + FN)
    FN = FN + FN2
    Prec = TP/(TP+FP)
    Sens = TP/(TP+FN)
    print TP,FP,FN,Prec,Sens






"""Run stuff here"""

s1 = time.time()
parentTree = generate_node_parent('nodes.dmp')
s2 = time.time()

s3 = time.time()
rankTree = generate_node_rank('nodes.dmp')
s4 = time.time()

precisionSensitivity('reads_hflu.out','mapping.uh',rankTree,parentTree)
'''
s5=time.time()
write_trueread_file('r1_only_cat_WEVOTE_DETAILS.txt','CLARKphylo.txt')
s6=time.time()
print(s6-s5)
'''

