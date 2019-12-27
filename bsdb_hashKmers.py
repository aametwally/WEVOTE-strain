#!/usr/bin/env python3
import os
import pickle


# Part 3 of database creation, turns k-mer files into searchable hash

#TODO:
# dynamic path
# better pickle name
# more organization

""" File navigation logistics """

#script_dir = os.path.dirname(os.path.realpath('__file__')) #<-- absolute dir the script is in
run_dir = os.getcwd() #Where we currently are, should be strain_dir
print(os.getcwd())
#rel_path= "bench_hflu_dir/jelly" #FIX
#abs_file_path = os.path.join(script_dir, rel_path)
jelly_path = os.getcwd()
jelly_files = [f for f in os.listdir(jelly_path) if f.endswith('.jdb')]
kdict = {}
#print(jelly_files)
print("this many jdb files",len(jelly_files))


""" Add kmer to dictionary """

for jf in jelly_files:
    with open(jf, 'r') as f:
        taxid = os.path.basename(f.name)
        # Note, file_extension is never used, delete it?
        taxid, file_extension = os.path.splitext(taxid)
        line = f.readline()
        # Go through every line in the file
        for line in f:
            line = line.strip()
            # If this is a unique line, add to dict
            if line not in kdict:
                kdict[line] = [taxid]
            # If this is already here, append the taxid to it
            else:
                kdict[line].append(taxid)


""" Filter out very common k-mers """

#key is kmer, value is taxid
# If the k-mer is found to exist in every jellyfish file, it is useless
# Note, this could have a threshold, you could set threshold to 90% and say if the kmer is found in .9 * len(jellyfiles) then remove it.

rmlist =[]
keeplist=[]
for key, value in kdict.items():
    if len(value) >= len(jelly_files):
        rmlist.append(key)
    else:
        keeplist.append(key)


for kmer in rmlist:
    del kdict[kmer]

print("K-mers that are not present in 100% of the files")
print(len(keeplist))

print("Removed K-mers (found in every genome)")
print(len(rmlist))


""" Export the new k-mer database as a python pickle """

pickle_out = open("strainDB.pickle","wb")
pickle.dump(kdict, pickle_out)
pickle_out.close()

