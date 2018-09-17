#include "kraken_headers.hpp"

using namespace std;

//Current goal:
//1. open fasta file
//2. parse fasta contents and loop through file
//3. take dnasequence, and run through kmer stuff (get this from kraken)
//
int main(int argc, char **argv){
    string current_line;
    string seq;
    string dnaseq;
    //open the file
    ifstream readsFile( argv[1] );
    // do stuff
    if (readsFile.is_open())
    {
        while ( getline (readsFile,current_line) )
        {
            cout<< current_line[0] << '\n';
            seq = current_line;
        }

    }
    dnaseq = seq.str();





























    return 0;
}


    























































































































































/*
#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "quickfile.hpp"
#include "seqreader.hpp"


const size_t DEF_WORK_UNIT_SIZE = 500000;

using namespace std;
using namespace kraken;

//function defs
void parse_command_line(int argc, char **argv);
void usage(int exit_code=EX_USAGE);
void process_file(char *filename);
void classify_sequence(DNASequence &dna, ostringstream &koss,
        ostringstream &coss, ostringstream &uoss,
        ostringstream &coss2, ostringstream &uoss2);
string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);
set<uint32_t> get_ancestry(uint32_t taxon);
void report_stats(struct timeval time1, struct timeval time2);

//variable defs - mostly switches
int Num_threads = 1;
string DB_filename, Index_filename, Nodes_filename;
//bool Quick_mode = false;
//bool Fastq_input = false;
//bool Fastq_output = false;
//bool Paired_input = false;
//bool Print_classified = false;
//bool Print_unclassified = false;
bool Print_kraken = true; //idk what this does really
bool Populate_memory = false;
bool Only_classified_kraken_output = false;

uint32_t Minimum_hit_count = 1;
map<uint32_t, uint32_t> Parent_map;
KrakenDB Database;
string Classified_output_file, Unclassified_output_file, Kraken_output_file;
string Output_format;
ostream *Classified_output;
ostream *Classified_output2;
ostream *Unclassified_output;
ostream *Unclassified_output2;
ostream *Kraken_output;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;
*/


