#!/usr/bin/env python3
import os
import pickle



# This script takes jellyfish files
# Needs to have dynamic path added so that the user can select the folder in command line
# At the end you will have a pickled folder that will be linked to the classify script




#TODO:
# dynamic path
# better pickle name
# more organization






script_dir = os.path.dirname(os.path.realpath('__file__')) #<-- absolute dir the script is in
rel_path= "bench_hflu_dir/jelly" #FIX
abs_file_path = os.path.join(script_dir, rel_path)
os.chdir(abs_file_path)
jelly_files = [f for f in os.listdir(abs_file_path) if f.endswith('.jdb')]
kdict = {}

print(jelly_files)
print("this many jdb files",len(jelly_files))
for jf in jelly_files:
    with open(jf, 'r') as f:
        taxid = os.path.basename(f.name)
        taxid, file_extension = os.path.splitext(taxid)
        line = f.readline()
        for line in f:
            line = line.strip()
            if line not in kdict:
                kdict[line] = [taxid]
            else:
                kdict[line].append(taxid)

rmlist =[]
keeplist=[]
for key, value in kdict.items():
    if len(value) >= len(jelly_files):
        rmlist.append(key)
    else:
        keeplist.append(key)


for kmer in rmlist:
    del kdict[kmer]

print('kept,removed')
print(len(keeplist),len(rmlist))

pickle_out = open("strainDB.pickle","wb")
pickle.dump(kdict, pickle_out)
pickle_out.close()

