#!/usr/bin/env python3
import os
import pickle



script_dir = os.path.dirname(os.path.realpath('__file__')) #<-- absolute dir the script is in
rel_path= "hflu_genomes_dir/jelly"
abs_file_path = os.path.join(script_dir, rel_path)
print(abs_file_path)
os.chdir(abs_file_path)
jelly_files = [f for f in os.listdir(abs_file_path) if f.endswith('.jdb')]
kdict = {}

print(jelly_files)
print("this many jdb files",len(jelly_files))
for jf in jelly_files:
    print(jf)
    with open(jf, 'r') as f:
        taxid = os.path.basename(f.name)
        taxid, file_extension = os.path.splitext(taxid)
        line = f.readline()
        for line in f:
            line = line.strip()
            if line not in kdict:
                kdict[line] = [taxid]
                # print('new kmer')
            else:
                kdict[line].append(taxid)
                # print(len(kdict[line]))

rmlist =[]
keeplist=[]
for key, value in kdict.items():
    if len(value) >= len(jelly_files):
        print("this kmer isn't unique enough" )
        rmlist.append(key)
    else:
        keeplist.append(key)


for kmer in rmlist:
    del kdict[kmer]

print('kept,removed')
print(len(keeplist),len(rmlist))

pickle_out = open("thin_dict.pickle","wb")
pickle.dump(kdict, pickle_out)
pickle_out.close()

