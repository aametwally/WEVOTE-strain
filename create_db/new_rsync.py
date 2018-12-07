"""hello"""
import os
import subprocess as sp
import fileinput



taxid = 1392

""" 1. Grab matching urls from the assembly summary"""
#grab
tmpfile = open('urls.tmp','w')
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
sp.run(['rm','urls.tmp'])

"""2. Call rsync """
#sp.run(["rsync","--dry-run","--no-motd","--files-from=urls.txt","rsync://ftp.ncbi.nlm.nih.gov/genomes/",".","2>","rsync.err"],check=True)

sp.run(["rsync","--no-motd","--files-from=urls.txt","rsync://ftp.ncbi.nlm.nih.gov/genomes/","."])
