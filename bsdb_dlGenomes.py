import os
import subprocess as sp
import fileinput
import sys
import argparse


"""
HOWTO:

This script downloads the required files from the NCBI database, including the taxonomy files (same for all species), followed by the complete genomes in NCBI for the given TAXID
Folder named $dbname_dir and located in cwd

USAGE:

python dlStrainGenomes.py $taxid $dbname

"""



def getArgs():
    """ Function to get input arguments, currently just <taxid> and <dbname>"""
    parser = argparse.ArgumentParser(description='',epilog='')
    parser.add_argument("taxid", help="This is the taxonomic ID",type=int)
    parser.add_argument("dbname", help="This is the database name")
    args = parser.parse_args()
    taxid = args.taxid
    dbname = args.dbname
    return taxid,dbname



taxid, dbname = getArgs()

# Generate folder and taxonomy folder
DB_directory = sys.argv[2] + "_dir"
lib_DB_directory = "../" + DB_directory

if not os.path.exists(DB_directory):
    os.mkdir(DB_directory)
os.chdir(DB_directory)

if not os.path.exists("taxonomy"):
    os.mkdir("taxonomy")
os.chdir("taxonomy")


# Download accession files
if not os.path.isfile("nucl_gb.accession2taxid.gz"):
    print("Accession to Taxon")
    # os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_est.accession2taxid.gz")
    os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz")
    # os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gss.accession2taxid.gz")

# Download and extract taxonomy dump
if not os.path.isfile("taxdump.tar.gz"):
    print("Taxonomy Tree")
    os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")

if not os.path.isfile("citations.dmp"):
    print("Extracting... (This takes a while)")
    os.system("gunzip *.gz")
    os.system("tar -xvf taxdump.tar")

# Download assembly summary
if not os.path.isfile("assembly_summary.txt"):
    print("Download assembly summary")
    os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")


# 2. Look at assembly summary and if you see the parent species == taxid, grab it.
tmpfile = open('urls.tmp','w')
#grab all complete genomes belonging to that species id #$7 is species aka parent taxid
#explanation: $6 is taxid (if at strain level then shouldn't be equal to species taxid
#$7 is parent taxid (aka should be equal to species)
#$12 is assembly_level which should be a complete genome or chromosome
#20 is the ftp file location


# Run shell script pulling genome names
tax_string = 'BEGIN { FS="\t"} { if ($6!=%d && $7==%d && ($12=="Complete Genome" || $12 =="Chromosome" )) print $20;}' % (taxid,taxid)
cmd = ["awk",tax_string,"assembly_summary.txt"]
sp.run(cmd,check=True, stdout=tmpfile) #only runs in python 3.5+
print("grabbed genomes")

# Delete header
findme="ftp://ftp.ncbi.nlm.nih.gov/genomes/"
with fileinput.input('urls.tmp',inplace=True) as f:
    for line in f:
        print(line.replace(findme,'').rstrip())

# Append file name
f2 = open('urls.txt','w')
cmd2=["awk","-F/",'{print $0"/" $6 "_genomic.fna.gz"}',"urls.tmp"]
sp.run(cmd2,check=True,stdout=f2)
#sp.run(['rm','urls.tmp'])

# Call rsync and download all the genomes
sp.run(["rsync","--no-motd","--files-from=urls.txt","rsync://ftp.ncbi.nlm.nih.gov/genomes/","."])

# Extract fasta files
print("Extracting fasta files")
os.system('gunzip -r all/*')
os.chdir('../')

os.system("mkdir -p library/added")
os.system("find . -name '*.fna' | xargs mv -t library/added/")


