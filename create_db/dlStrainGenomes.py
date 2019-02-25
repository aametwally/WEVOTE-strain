import os
import subprocess as sp
import fileinput
import sys
import argparse



def getArgs():
    """ Function to get input arguments, currently just <taxid> and <dbname"""
    parser = argparse.ArgumentParser(description='',epilog='')
    parser.add_argument("taxid", help="This is the taxonomic ID",type=int)
    parser.add_argument("dbname", help="This is the database name")
    args = parser.parse_args()
    print(args.taxid)
    taxid = args.taxid
    dbname = args.dbname
    return taxid,dbname

    

taxid, dbname = getArgs()



#Get a list of species to look for
#Currently this part depends on a piece of awk code I ran (in my boostnotes) , needs to 
#incorporated into the code as to automatically pull that species. But right now just works for a list
#of taxids which are all ruminococcus species

DB_directory = sys.argv[2] + "_dir"
#os.mkdir(DB_directory)
os.chdir(DB_directory)

#Run kraken-build, currently writing a lite version but this works best for now
krakencall = "kraken-build --download-taxonomy --db " + sys.argv[2]
os.system(krakencall)

"""
###TAXONOMY STUFF
#1. Get accession to taxon map
print("Accession to Taxon")
os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz")
os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz")
os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz")

print("Taxonomy Tree")
os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
print("Extracting... (This takes a while)")
os.system("gunzip *.gz")
os.system("tar -xvf taxdump.tar")
"""

#1.3 Get assembly_summary.txt
print("Download assembly summary")
os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")


#This is to generate a list of species, not currently used

# 2. Look at assembly summary and if you see the parent species == taxid, grab it.
tmpfile = open('urls.tmp','w')
#grab
tax_string = '{if ($7==%d) print $NF;}' % taxid
cmd = ["awk",tax_string,"assembly_summary.txt"]
sp.run(cmd,check=True, stdout=tmpfile) #only runs in python 3.5+

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
sp.run(["rsync","--no-motd","--files-from=urls.txt","rsync://ftp.ncbi.nlm.nih.gov/genomes/","."])

#3. Extract and build
print("Extracting Reads")
os.system('gunzip -r all/*')

os.system("find . -name '*.fna' -print0 |     xargs -0 -I{} -n1 kraken-build --add-to-library {} --db %s"%dbname)

#Final step: BUILD DATABASE
#This needs to be done on qsubmission
print("building...")
#os.system("kraken-build --build --db %s"%dbname)


