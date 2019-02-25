import argparse

def getArgs():
    """ Function to get input arguments"""
    parser = argparse.ArgumentParser(description='',epilog='')
    parser.add_argument("taxid", help="This is the taxonomic ID",type=int)
    parser.add_argument("wevote_file", help="The WEVOTE Log file containing the annotation of each read",type=int)
    parser.add_argument("original_fasta", help="This is the original fasta file that was an input to wevote",type=int)
    parser.add_argument("species_fasta", help="User-defined output name (fasta file of species)")
    args = parser.parse_args()
    taxid = args.taxid
    wevote_file= args.dbname
    og_fasta = args.original_fasta
    outfile = args.species_fasta
    return taxid,dbname,og_fasta,outfile
    return 0

    

taxid, wevote_file, og_fasta, outfile = getArgs()
getArgs()
