"""hello"""
import os
import subprocess as sp
import fileinput



#taxid = 562
#taxid = 1392 #bacillus anthracus
#taxid = 1263 #ruminococcus(genus)


#Get a list of species to look for
#Currently this part depends on a piece of awk code I ran (in my boostnotes) , needs to 
#incorporated into the code as to automatically pull that species. But right now just works for a list
#of taxids which are all ruminococcus species


#I lost the boostnotes to get this species list, but it was some awk stuff 
with open('rumino_species.txt') as rf:
    lines = rf.read().splitlines()
results = list(map(int,lines))
print(results)

# 1. Look at assembly summary and if you see the parent species == taxid, grab it.
tmpfile = open('urls.tmp','w')
for taxid in results:
    #grab
    tax_string = '{if ($7==%d) print $NF;}' % taxid
    cmd = ["awk",tax_string,"assembly_summary.txt"]
    sp.run(cmd,check=True, stdout=tmpfile)

#Delete header
findme="ftp://ftp.ncbi.nlm.nih.gov/genomes/"
with fileinput.input('urls.tmp',inplace=True) as f:
    for line in f:
        print(line.replace(findme,'').rstrip())

#Append file name
f2 = open('urls.txt','w')
cmd2=["awk","-F/",'{print $0"/" $6 "_genomic.fna.gz"}',"urls.tmp"]
sp.run(cmd2,check=True,stdout=f2)
#sp.run(['rm','urls.tmp'])

#2. Call rsync 
#sp.run(["rsync","--dry-run","--no-motd","--files-from=urls.txt","rsync://ftp.ncbi.nlm.nih.gov/genomes/",".","2>","rsync.err"],check=True)

sp.run(["rsync","--no-motd","--files-from=urls.txt","rsync://ftp.ncbi.nlm.nih.gov/genomes/","."])
